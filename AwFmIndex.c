#include "AwFmIndex.h"
#include "AwFmGlobals.h"
#include <string.h>
#include <stdlib.h>


/*
 * Function:  awFmIndexSetFileSrc
 * --------------------
 * Sets the fileSrc member variable inside the index to the given fileSrc.
 *  On success, the fileSrc argument will be copied into a newly allocated string in the index.
 *  As such, the given fileSrc argument may be deallocated freely after calling this function.   
 *
 *  Inputs:
 *    index:      Pointer to an allocated AwFmIndex struct that will have its fileSrc set.
 *    fileSrc:    Null-terminated string representing the absolute location of the file
 *      referenced by this FM-index.
 *
 *  Returns:
 *    AwFmReturnCode showing the result of this action. Possible returns are:
 *      AwFmSuccess on success,
 *      AwFmNoFileSrcGiven if the fileSrc was null,
 *      AwFmGeneralFailure if a null terminator is not found in the fileSrc in a reasonable
 *        number of characters. This function will look for as many characters as the value
 *        in the MAXIMUM_FILE_PATH_LENGTH macro.
 *      AwFmAllocationFailure if the dynamic allocation of the copy of the fileSrc.
 */
enum AwFmReturnCode awFmIndexSetFileSrc(struct AwFmIndex *restrict const index, const char *restrict const fileSrc){
  //check to make sure the fileSrc isn't null
  if(__builtin_expect(fileSrc == NULL, 0)){
    return AwFmNoFileSrcGiven;
  }

  const void *nullTerminatorPosition = memchr(fileSrc, 0, MAXIMUM_FILE_PATH_LENGTH);
  if(__builtin_expect(nullTerminatorPosition == NULL, 0)){
    return AwFmGeneralFailure;
  }

  //allocate a copy of the fileSrc string, and copy it over (null terminator included)
  const size_t fileSrcLength = (char*)nullTerminatorPosition - fileSrc;
  index->fileSrc = malloc((fileSrcLength + 1) * sizeof(char));
  if(__builtin_expect(index->fileSrc == NULL, 0)){
    return AwFmAllocationFailure;
  }

  memcpy(index->fileSrc, fileSrc, (fileSrcLength + 1) * sizeof(char));

  return AwFmSuccess;
}


/*
 * Function:  alignedAllocAwFmIndex
 * --------------------
 * Dynamically allocates memory for the AwFmIndex struct, aligned to a cache line (64 bytes).
 *  Inputs:
 *    fileSrc:    Null-terminated string representing the absolute location
 *      of the file to read index from, or write index to.
 *      If for some reason file reading/writing is not needed, fileSrc may be left NULL.
 *      This string will be copied into a dynamically allocated string, and the AwFmIndex struct
 *      will be responsible for the lifetime of this copy.
 *
 *  Returns:
 *    Pointer to the allocated aligned memory. Returns NULL on allocation failure.
 */
struct AwFmIndex *alignedAllocAwFmIndex(void){
  struct AwFmIndex *index = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, sizeof(struct AwFmIndex));
  if(index != NULL){
    memset(index, 0, sizeof(struct AwFmIndex));
  }
  return index;
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
  return aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, numBlocks * sizeof(struct AwFmBlock));
}


/*
 * Function:  deallocateFmIndex
 * --------------------
 * Frees the memory associated with the given AwFmIndex. This also deallocates the block list.
 *
 *  Inputs:
 *    index: Index struct to be deallocated.
 */
 void deallocateFmIndex(struct AwFmIndex *restrict index){

   if(index != NULL){
     if(index->fileSrc != NULL){
       free(index->fileSrc);
     }
     if(index->blockList != NULL){
       free(index->blockList);
     }

     free(index);
   }
 }


/*
 * Function:  awFmSearchRangeLength
 * --------------------
 * Gets the number of positions included in the given AwFmSearchRange
 *
 *  Inputs:
 *    range: Range of positions in the BWT that corresponds to some number of
 *      instances of a given kmer.
 *
 *  Outputs:
 *    Number of positions in the given range if the range is valid (startPtr < endPtr),
 *      or 0 otherwise, as that would imply that no instances of that kmer were found.
 */
 size_t awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range){
  uint64_t length = range->endPtr - range->startPtr;
  return (range->startPtr < range->endPtr)? length: 0;
}

/*
 * Function:  awFmBwtPositionIsSampled
 * --------------------
 * Determines if the given BWT position is sampled under the given AwFmIndex's
 *  suffix array compression ratio.
 *
 *  Inputs:
 *    index:      Pointer to the BWT's AwFmIndex.
 *    position:   Position in the BWT.
 *
 *  Outputs:
 *    True if the position is sampled in the compressedSuffixArray, or false otherwise.
 */
 bool awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index, const uint64_t position){
  return position % index->suffixArrayCompressionRatio;
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


 uint64_t awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index){
  return awFmGetBwtLength(index) / index->suffixArrayCompressionRatio;
}


/*
 * Function:  awFmSearchRangeIsValid
 * --------------------
 * Compares the start pointer and end pointer to determine if this range
 *  contains a valid range of the relevant kmer suffix.
 *
 *  Inputs:
 *    searchRange: AwFmSearchRange struct to compare.
 *
 *  Outputs:
 *    True if this range represents a range of positions for the given kmer suffix in the database,
 *      or false if the database does not contain the given kmer suffix.
 */
 bool awFmSearchRangeIsValid(const struct AwFmSearchRange *restrict const searchRange){
  return searchRange->startPtr < searchRange->endPtr;
}
