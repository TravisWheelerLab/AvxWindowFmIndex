#include "AwFmIndexStruct.h"
#include <stdlib.h>
#include <string.h>
#include "AwFmIndex.h"
#include "FastaVector.h"

#define AW_FM_BWT_BYTE_ALIGNMENT 32

struct AwFmIndex *
awFmIndexAlloc(const struct AwFmIndexConfiguration *_RESTRICT_ const config,
               const size_t bwtLength) {

  // allocate the index
  struct AwFmIndex *index = malloc(sizeof(struct AwFmIndex));
  if (index == NULL) {
    return NULL;
  }
  // initialize all bytes in the index to 0.
  memset(index, 0, sizeof(struct AwFmIndex));
  memcpy(&index->config, config, sizeof(struct AwFmIndexConfiguration));
  index->bwtLength = bwtLength;

  // allocate the prefixSums
  size_t prefixSumsLength = awFmGetPrefixSumsLength(config->alphabetType);
  index->prefixSums = malloc(prefixSumsLength * sizeof(uint64_t));
  if (index->prefixSums == NULL) {
    awFmDeallocIndex(index);
    return NULL;
  }

  // allocate the blockLists
  size_t numBlocksInBwt = awFmNumBlocksFromBwtLength(bwtLength);
  size_t sizeOfBwtBlock = config->alphabetType == AwFmAlphabetAmino
                              ? sizeof(struct AwFmAminoBlock)
                              : sizeof(struct AwFmNucleotideBlock);

  // alloc the backward bwt
  index->bwtBlockList.asNucleotide =
      aligned_alloc(AW_FM_BWT_BYTE_ALIGNMENT, numBlocksInBwt * sizeOfBwtBlock);
  if (index->bwtBlockList.asNucleotide == NULL) {
    awFmDeallocIndex(index);
    return NULL;
  }

  const size_t kmerSeedTableSize = awFmGetKmerTableLength(index);
  // allocate the kmerSeedTable
  index->kmerSeedTable =
      malloc(kmerSeedTableSize * sizeof(struct AwFmSearchRange));
  if (index->kmerSeedTable == NULL) {
    awFmDeallocIndex(index);
    return NULL;
  }

  return index;
}

void awFmDeallocIndex(struct AwFmIndex *index) {
  if (index != NULL) {
    fclose(index->fileHandle);
    free(index->bwtBlockList.asNucleotide);
    free(index->prefixSums);
    free(index->kmerSeedTable);
    free(index->suffixArray.values);
    if (index->fastaVector != NULL) {
      fastaVectorDealloc(index->fastaVector);
      free(index->fastaVector);
    }
    free(index);
  }
}

uint_fast8_t awFmGetAlphabetCardinality(const enum AwFmAlphabetType alphabet) {
  return (alphabet == AwFmAlphabetAmino) ? AW_FM_AMINO_CARDINALITY
                                         : AW_FM_NUCLEOTIDE_CARDINALITY;
}

size_t awFmGetKmerTableLength(const struct AwFmIndex *_RESTRICT_ index) {
  const size_t multiplier =
      awFmGetAlphabetCardinality(index->config.alphabetType);
  size_t length = 1;
  for (size_t i = 0; i < index->config.kmerLengthInSeedTable; i++) {
    length *= multiplier;
  }

  return length;
}

bool awFmBwtPositionIsSampled(const struct AwFmIndex *_RESTRICT_ const index,
                              const uint64_t position) {
  return (position % index->config.suffixArrayCompressionRatio) == 0;
}

uint64_t awFmGetCompressedSuffixArrayLength(
    const struct AwFmIndex *_RESTRICT_ const index) {
  return 1 +
         ((index->bwtLength - 1) / index->config.suffixArrayCompressionRatio);
}

bool awFmSearchRangeIsValid(
    const struct AwFmSearchRange *_RESTRICT_ const searchRange) {
  return searchRange->startPtr <= searchRange->endPtr;
}

size_t awFmNumBlocksFromBwtLength(const size_t suffixArrayLength) {
  return 1 + ((suffixArrayLength - 1) / AW_FM_POSITIONS_PER_FM_BLOCK);
}

uint8_t awFmGetPrefixSumsLength(const enum AwFmAlphabetType alphabet) {
  return awFmGetAlphabetCardinality(alphabet) +
         2; // 1 for sentinel count, 1 for bwt length
}

bool awFmReturnCodeSuccess(const enum AwFmReturnCode returnCode) {
  return returnCode >= 0;
}

size_t awFmGetBlockIndexFromGlobalPosition(const size_t globalQueryPosition) {
  return globalQueryPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
}

uint_fast8_t
awFmGetBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition) {
  return globalQueryPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
}

size_t
awFmSearchRangeLength(const struct AwFmSearchRange *_RESTRICT_ const range) {
  uint64_t length = range->endPtr - range->startPtr;
  return (range->startPtr <= range->endPtr) ? length + 1 : 0;
}

bool awFmIndexIsVersionValid(const uint16_t versionNumber) {
  return versionNumber == AW_FM_CURRENT_VERSION_NUMBER;
}

bool awFmIndexContainsFastaVector(
    const struct AwFmIndex *_RESTRICT_ const index) {
  return index->featureFlags & (1 << AW_FM_FEATURE_FLAG_BIT_FASTA_VECTOR);
}

inline bool awFmReturnCodeIsFailure(const enum AwFmReturnCode rc) {
  return rc < 0;
}

inline bool awFmReturnCodeIsSuccess(const enum AwFmReturnCode rc) {
  return !awFmReturnCodeIsFailure(rc);
}

uint32_t awFmGetNumSequences(const struct AwFmIndex *_RESTRICT_ const index) {
  if (index->fastaVector) {
    return index->fastaVector->metadata.count;
  } else {
    return 1;
  }
}