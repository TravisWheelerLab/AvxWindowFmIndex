#include "AwFmIndex.h"
#include <string.h>
#include <stdlib.h>


struct AwFmIndex *awFmIndexAlloc(const struct AwFmIndexMetadata *restrict const metadata,
  const size_t bwtLength){

  //allocate the index
  struct AwFmIndex *index = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
    sizeof(struct AwFmIndex));
  if(index == NULL){
    return NULL;
  }
  //initialize all bytes in the index to 0.
  memset(index, 0, sizeof(struct AwFmIndex));
  memcpy(&index->metadata, metadata, sizeof(struct AwFmIndexMetadata));
  index->bwtLength = bwtLength;

  //allocate the prefixSums
  size_t prefixSumsLength = awFmGetPrefixSumsLength(metadata->alphabetType);
  index->prefixSums = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
    prefixSumsLength * sizeof(uint64_t));
  if(index->prefixSums == NULL){
    awFmDeallocIndex(index);
    return NULL;
  }

  //allocate the blockLists
  size_t numBlocksInBwt = awFmNumBlocksFromBwtLength(bwtLength);
  size_t sizeOfBwtBlock = metadata->alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

  //alloc the backward bwt
  index->bwtBlockList.asNucleotide = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
    numBlocksInBwt * sizeOfBwtBlock);
  if(index->bwtBlockList.asNucleotide == NULL){
    awFmDeallocIndex(index);
    return NULL;
  }

  const size_t kmerSeedTableSize = awFmGetKmerTableLength(index);
  //allocate the kmerSeedTable
  printf("allocating kmer table of length %zu\n", kmerSeedTableSize);
  index->kmerSeedTable.table = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
    kmerSeedTableSize * sizeof(struct AwFmSearchRange));
  if(index->kmerSeedTable.table == NULL){
    awFmDeallocIndex(index);
    return NULL;
  }

  //allocate the sequence ending table
  index->kmerSeedTable.sequenceEndingKmerEncodings = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
    (metadata->kmerLengthInSeedTable - 1) * sizeof(uint64_t));
  if(index->kmerSeedTable.sequenceEndingKmerEncodings == NULL){
    awFmDeallocIndex(index);
    return NULL;
  }

  return index;
}


void awFmDeallocIndex(struct AwFmIndex *index){
  if(index != NULL){
    fclose(index->fileHandle);
    free(index->bwtBlockList.asNucleotide);
    free(index->prefixSums);
    free(index->kmerSeedTable.table);
    free(index->kmerSeedTable.sequenceEndingKmerEncodings);
    free(index);
  }
}


uint_fast8_t awFmGetAlphabetCardinality(const enum AwFmAlphabetType alphabet){
  return (alphabet == AwFmAlphabetNucleotide)?
    AW_FM_NUCLEOTIDE_CARDINALITY: AW_FM_AMINO_CARDINALITY;
}


size_t awFmGetKmerTableLength(const struct AwFmIndex *restrict index){
  const size_t multiplier = awFmGetAlphabetCardinality(index->metadata.alphabetType);
  size_t length = 1;
  for(size_t i = 0; i < index->metadata.kmerLengthInSeedTable; i++){
    length *= multiplier;
  }

  return length;
}


 bool awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index, const uint64_t position){
  return (position % index->metadata.suffixArrayCompressionRatio )== 0;
}


 uint64_t awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index){
  return index->bwtLength / index->metadata.suffixArrayCompressionRatio;
}


 bool awFmSearchRangeIsValid(const struct AwFmSearchRange *restrict const searchRange){
  return searchRange->startPtr <= searchRange->endPtr;
}


size_t awFmNumBlocksFromBwtLength(const size_t suffixArrayLength){
  return  1 + ((suffixArrayLength -1) / AW_FM_POSITIONS_PER_FM_BLOCK);
}

uint8_t awFmGetPrefixSumsLength(const enum AwFmAlphabetType alphabet){
  return  awFmGetAlphabetCardinality(alphabet) +1;
}

bool awFmReturnCodeSuccess(const enum AwFmReturnCode returnCode){
  return returnCode >= 0;
}


size_t awFmGetBlockIndexFromGlobalPosition(const size_t globalQueryPosition){
 return globalQueryPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
}


uint_fast8_t awFmGetBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition){
  return globalQueryPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
}


 size_t awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range){
  uint64_t length = range->endPtr - range->startPtr;
  return (range->startPtr <= range->endPtr)? length + 1: 0;
}
