#ifndef AW_FM_INDEX_STRUCTS_H
#define AW_FM_INDEX_STRUCTS_H

#include <stdint.h>
#include <stdbool.h>
#include <immintrin.h>
#include "AwFmGlobals.h"


#define AW_FM_INDEX_METADATA_BYTE_SIZE 64

/*Structs*/
struct AwFmBlock{
  __m256i   letterBitVectors[5];
  uint64_t  baseOccupancies[20];
};

struct AwFmIndexMetadata{
  uint32_t  versionNumber;
};

union AwFmIndexPaddedMetadata{
  struct AwFmIndexMetadata data;
  uint8_t padding[AW_FM_INDEX_METADATA_BYTE_SIZE];
};

struct AwFmIndex{
  struct AwFmBlock *blockList;
  uint64_t  numBlocks;
  uint16_t  suffixArrayCompressionRatio;
  uint64_t  rankPrefixSums[AMINO_CARDINALITY + 1]; //last position acts as BWT length
  char      *fileSrc;
  uint8_t   *databaseSequence;  //usually NULL, used in construction and saving to file
  uint64_t  *fullSuffixArray;   //usually NULL, used in construction and saving to file
  union AwFmIndexPaddedMetadata metadata;
};

struct AwFmSearchRange{
  uint64_t startPtr;
  uint64_t endPtr;
};


/*Public Functions*/

//API Function

struct  AwFmIndex *alignedAllocAwFmIndex(const char *restrict const fileSrc);
struct  AwFmBlock *alignedAllocBlockList(const size_t numBlocks);
        void      deallocateFmIndex(struct AwFmIndex *restrict index);
inline  bool      awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index, const uint64_t position);
inline  size_t    awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range);
inline  uint64_t  awFmGetBwtLength(const struct AwFmIndex *restrict const index);
inline  uint64_t  awFmGetDbSequenceLength(const struct AwFmIndex *restrict const index);
inline  uint64_t  awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index);
inline bool       awFmSearchRangeIsValid(const struct AwFmSearchRange *restrict const searchRange);

#endif /* end of include guard: AW_FM_INDEX_STRUCTS_H */
