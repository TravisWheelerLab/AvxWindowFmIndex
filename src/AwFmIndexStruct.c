#include "AwFmIndexStruct.h"
#include "AwFmIndex.h"
#include <stdlib.h>
#include <string.h>

#define AW_FM_BWT_BYTE_ALIGNMENT 32

struct AwFmIndex *awFmIndexAlloc(const struct AwFmIndexMetadata *restrict const metadata,
                                 const size_t bwtLength) {

    // allocate the index
    struct AwFmIndex *index = malloc(sizeof(struct AwFmIndex));
    if (index == NULL) {
        return NULL;
    }
    // initialize all bytes in the index to 0.
    memset(index, 0, sizeof(struct AwFmIndex));
    memcpy(&index->metadata, metadata, sizeof(struct AwFmIndexMetadata));
    index->bwtLength = bwtLength;

    // allocate the prefixSums
    size_t prefixSumsLength = awFmGetPrefixSumsLength(metadata->alphabetType);
    index->prefixSums = malloc(prefixSumsLength * sizeof(uint64_t));
    if (index->prefixSums == NULL) {
        awFmDeallocIndex(index);
        return NULL;
    }

    // allocate the blockLists
    size_t numBlocksInBwt = awFmNumBlocksFromBwtLength(bwtLength);
    size_t sizeOfBwtBlock = metadata->alphabetType == AwFmAlphabetNucleotide
                                ? sizeof(struct AwFmNucleotideBlock)
                                : sizeof(struct AwFmAminoBlock);

    // alloc the backward bwt
    index->bwtBlockList.asNucleotide =
        aligned_alloc(AW_FM_BWT_BYTE_ALIGNMENT, numBlocksInBwt * sizeOfBwtBlock);
    if (index->bwtBlockList.asNucleotide == NULL) {
        awFmDeallocIndex(index);
        return NULL;
    }

    const size_t kmerSeedTableSize = awFmGetKmerTableLength(index);
    // allocate the kmerSeedTable
    index->kmerSeedTable = malloc(kmerSeedTableSize * sizeof(struct AwFmSearchRange));
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
        free(index->inMemorySuffixArray);
        free(index);
    }
}

uint_fast8_t awFmGetAlphabetCardinality(const enum AwFmAlphabetType alphabet) {
    return (alphabet == AwFmAlphabetNucleotide) ? AW_FM_NUCLEOTIDE_CARDINALITY
                                                : AW_FM_AMINO_CARDINALITY;
}

size_t awFmGetKmerTableLength(const struct AwFmIndex *restrict index) {
    const size_t multiplier = awFmGetAlphabetCardinality(index->metadata.alphabetType);
    size_t length = 1;
    for (size_t i = 0; i < index->metadata.kmerLengthInSeedTable; i++) {
        length *= multiplier;
    }

    return length;
}

bool awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index,
                              const uint64_t position) {
    return (position % index->metadata.suffixArrayCompressionRatio) == 0;
}

uint64_t awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index) {
    return 1 + ((index->bwtLength - 1) / index->metadata.suffixArrayCompressionRatio);
}

bool awFmSearchRangeIsValid(const struct AwFmSearchRange *restrict const searchRange) {
    return searchRange->startPtr <= searchRange->endPtr;
}

size_t awFmNumBlocksFromBwtLength(const size_t suffixArrayLength) {
    return 1 + ((suffixArrayLength - 1) / AW_FM_POSITIONS_PER_FM_BLOCK);
}

uint8_t awFmGetPrefixSumsLength(const enum AwFmAlphabetType alphabet) {
    return awFmGetAlphabetCardinality(alphabet) + 2; // 1 for sentinel count, 1 for bwt length
}

bool awFmReturnCodeSuccess(const enum AwFmReturnCode returnCode) { return returnCode >= 0; }

size_t awFmGetBlockIndexFromGlobalPosition(const size_t globalQueryPosition) {
    return globalQueryPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
}

uint_fast8_t awFmGetBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition) {
    return globalQueryPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
}

size_t awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range) {
    uint64_t length = range->endPtr - range->startPtr;
    return (range->startPtr <= range->endPtr) ? length + 1 : 0;
}
