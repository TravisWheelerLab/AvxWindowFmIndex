#include "FmIndex.c"
#include "FmAwGlobals.h"


struct AwFmIndex *alignedAllocAwFmIndex(void){
  return aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, sizeof(struct FmIndex));
}

struct AwFmBlock *alignedAllocBlockList(const size_t numBlocks){
  return aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, numBlocks * sizeof(FM_block));

}

size_t *allocAndInitFullTextPositionLookupTable(const size_t numPositions){
  const size_t *positionList = malloc(numPositions * sizeof(size_t));
  for(size_t i = 0; i< numPositions; i++){
    positionList[i] = i;
  }

  return positionList;
}

void deallocateFmIndex(struct FmIndex *restrict index){
  if(index != NULL){
    if(index->blockList != NULL){
      free(index->blockList)
    }

    if(index->fullTextPositionLookupTable != NULL){
      free(index->fullTextPositionLookupTable)
    }
    
    free(index);
  }
}
