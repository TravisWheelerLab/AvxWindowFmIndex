#include "FmIndex.c"
#include "FmAwGlobals.h"


//TODO: when creating rankPrefixSums, make sure that the prefix sum for A ia 1 (for the $ terminator)

/*
 * Function:  alignedAllocAwFmIndex
 * --------------------
 * Dynamically allocates memory for the AwFmIndex struct, aligned to a cache line (64 bytes).
 *
 *  Returns:
 *    Pointer to the allocated aligned memory. As per C11 specs, will return NULL on allocation failure.
 */
struct AwFmIndex *alignedAllocAwFmIndex(void){
  return aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, sizeof(struct FmIndex));
}


/*
 * Function:  alignedAllocBlockList
 * --------------------
 * Dynamically allocates memory for the AwFmBlock array, aligned to a cache line (64 bytes).
 *
 *  Inputs:
 *    numBlocks: size of the allocated AwFmBlock array.
 *
 *  Returns:
 *    Pointer to the allocated aligned memory. As per C11 specs, will return NULL on allocation failure.
 */
struct AwFmBlock *alignedAllocBlockList(const size_t numBlocks){
  return aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, numBlocks * sizeof(FM_block));

}


/*
 * Function:  deallocateFmIndex
 * --------------------
 * Frees the memory associated with the given AwFmIndex. This also deallocates the block list.
 *
 *  Inputs:
 *    index: Index struct to be deallocated.
 */
void deallocateFmIndex(struct FmIndex *restrict index){
  if(index != NULL){
    if(index->blockList != NULL){
      free(index->blockList)
    }

    free(index);
  }
}

/*
 * Function:  awFmGetBwtLength
 * --------------------
 * Gets the length of the BWT from the AwFmIndex data structure.
 *   The AwFmIndex stores the BWT length at the end of the rankPrefixSums array.
 *
 *  Inputs:
 *    index: AwFmIndex struct to query.
 *
 *  Outputs:
*     length of the BWT array for the given AwFmIndex.
 */
uint64_t awFmGetBwtLength(const struct AwFmIndex *restrict const index){
  //the bwt length is stored in the prefix sums after all the actual amino acid prefix sums.
  //conceptually, the bwt length == prefix sum for a character after every alphabet letter.
  return index->rankPrefixSums[AMINO_CARDINALITY];
}

/*
 * Function:  awFmGetDbSequenceLength
 * --------------------
 * Gets the length of the database sequence from the AwFmIndex data structure.
 *   Since the database sequence doesn't contain a null terminator ($), it is
 *   one less than the BWT length.
 *
 *  Inputs:
 *    index: AwFmIndex struct to query.
 *
 *  Outputs:
 *    length of the database sequence for the given AwFmIndex.
 */
//the database sequence contains every character from the bwt, minus the null terminator $.
uint64_t awFmGetDbSequenceLength(const struct AwFmIndex *restrict const index){
  return awFmGetBwtLength(index) - 1;
}
