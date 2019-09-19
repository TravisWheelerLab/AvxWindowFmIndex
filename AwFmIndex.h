#ifndef AW_FM_INDEX_STRUCTS_H
#define AW_FM_INDEX_STRUCTS_H

#include <stdint.h>

/*Structs*/
struct AwFmBlock{
  __m256i   letterBitVectors[5];
  uint64_t  baseOccupancies[20];
};

struct AwFmIndexMetadata{
  uint32_t  versionNumber;
  uint8_t   IArrayPositionWidth;
};

union AwFmIndexPaddedMetadata{
  struct AwFmIndexMetadata data;
  uint8_t padding[FM_INDEX_METADATA_BYTE_SIZE];
};

struct AwFmIndex{
  struct AwFmBlock *blockList;
  uint64_t numBlocks;
  uint64_t rankPrefixSums[AMINO_CARDINALITY + 1]; //last position acts as BWT length
  uint64_t nullTerminatorPosition;
  size_t *fullTextPositionLookupTable;
  union AwFmIndexPaddedMetadata metadata;
};


/*Public Functions*/
//API Function
struct  AwFmIndex *alignedAllocAwFmIndex(void);
//API Function
struct  AwFmBlock *alignedAllocBlockList(const size_t numBlocks);
        size_t    *allocAndInitFullTextPositionLookupTable(const size_t numPositions);
        void      deallocateFmIndex(const struct FmIndex *restrict const index);

#endif /* end of include guard: AW_FM_INDEX_STRUCTS_H */
