#include "AwFmIndex.h"
#include "AwFmOccurrence.h"
#include "AwFmSearch.h"

/*
 * Function:  AwFmIterativeRangeBackwardSearch
 * --------------------
 *  Queries the FM-Index to iteratively prepend an additional letter to the
 *    implicit kmer suffix, and return a new range of positions in the BWT where
 *    the extended kmer suffix may be found.
 *
 *  If the given kmer is not found, the AwFmSearch Range will result in a false value
 *   when given to the awFmSearchRangeIsValid function.
 *
 *  Inputs:
 *    index:          Pointer to the valid AwFmIndex struct.
 *    currentRange:   Range of positions in the BWT where the current kmer suffix may be found.
 *    prefixLetter:   Letter that will be prepended to the implicit kmer to extend the kmer suffix.
 *
 *  Returns:
 *    AwFmSearch range representing the range of BWT positions where the implicit
 *    kmer suffixmay be found, as long as startPtr < endPtr. Otherwise (startPtr >= endPtr),
 *    the given kmer suffix does not exist in the database sequence.
 */
struct AwFmSearchRange AwFmIterativeRangeBackwardSearch(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const currentRange, const uint8_t prefixLetter){

  struct AwFmSearchRange range;
  const uint64_t letterRankPrefixSum = index->rankPrefixSums[prefixLetter];

  range.startPtr = awFmGetOccurances(index, currentRange->startPtr, prefixLetter) + letterRankPrefixSum;
  awFmOccurancesDataPrefetch(index, range.startPtr);

  range.endPtr = awFmGetOccurances(index, currentRange->endPtr, prefixLetter) + letterRankPrefixSum;
  awFmOccurancesDataPrefetch(index, range.endPtr);

  return range;
}

//TODO: check this functionality against the literature.
/*
 * Function:  AwFmIterativeRangeForwardSearch
 * --------------------
 *  Queries the FM-Index to iteratively append an additional letter to the
 *    implicit kmer suffix, and return a new range of positions in the BWT where
 *    the extended kmer suffix may be found.
 *
 *  If the given kmer is not found, the AwFmSearch Range will result in a false value
 *   when given to the awFmSearchRangeIsValid function.
 *
 *  Inputs:
 *    index:          Pointer to the valid AwFmIndex struct.
 *    currentRange:   Range of positions in the BWT where the current kmer suffix may be found.
 *    suffixLetter:   Letter that will be appended to the implicit kmer to extend the kmer suffix.
 *
 *  Returns:
 *    AwFmSearch range representing the range of BWT positions where the implicit
 *    kmer suffixmay be found, as long as startPtr < endPtr. Otherwise (startPtr >= endPtr),
 *    the given kmer suffix does not exist in the database sequence.
 */
struct AwFmSearchRange AwFmIterativeRangeForwardSearch(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const currentRange, const uint8_t suffixLetter){

      struct AwFmSearchRange range;
      const uint64_t letterRankPrefixSum = index->rankPrefixSums[suffixLetter];

      range.startPtr = awFmGetOccurances(index, currentRange->startPtr, suffixLetter) + letterRankPrefixSum;
      awFmOccurancesDataPrefetch(index, range.startPtr);

      range.endPtr = awFmGetOccurances(index, currentRange->endPtr, suffixLetter) + letterRankPrefixSum;
      awFmOccurancesDataPrefetch(index, range.endPtr);

      return range;
}



/*
 * Function:  awFmFindDatabaseHitPositionsFromSearchRange
 * --------------------
 *  Takes a range of BWT positions, backtraces each position to find the nearest sample in the
 *    compressed suffix array, and looks up those suffix array positions on disk to
 *    determine the corresponding database sequence position for each BWT position
 *    between the searchRange's pointers (inclusive startPtr, exclusive endPtr).
 *
 *  It is the caller's responsibility to free() the returned sequence position array.
 *
 *  This function will overwrite the data in the positionArray, returning the
 *    database sequence positions in the corresponding elements of the positionArray.
 *
 *  Note: positionArray and offsetArray should be of equal length. Having either shorter
 *    than positionArrayLength will result in undefined behavior.
 *
 *  Inputs:
 *    index:              Pointer to the valid AwFmIndex struct.
 *    positionArray:          Array of positions in the implied full suffix array to load.
 *      This function will convert these positions to indices in the compressed suffix array,
 *      as long as the positions are multiples of the suffix array compression ratio.
 *    offsetArray:  array of offsets to be added to the database sequence positions.
 *      This is needed because we can only query the suffix array on positions that are sampled.
 *    positionArrayLength:    Length of the positionArray and offsetArray.
 *
 *  Returns:
 *    AwFmReturnCode detailing the result of the read attempt. Possible return values:
 *      AwFmFileReadOkay on success,
 *      AwFmFileOpenFail on failure to open the AwFm Index file
 *      AwFmFileReadFail on failure to read as many characters as was expected by the sequence.
 *      AwFmIllegalPositionError on a suffix array position being out of bounds of
 *        the file's compressed suffix array.
 */
uint64_t *awFmFindDatabaseHitPositionsFromSearchRange(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const searchRange,
  enum AwFmReturnCode *restrict fileAccessResult){

  const uint64_t numPositionsInRange            = awFmSearchRangeLength(searchRange);

  //allocate the position array and offsetArray (uses one malloc call, only needs 1 free call)
  //restrict keywords should be okay, since they shouldn't read/write to overlapping elements
  uint64_t *const restrict positionArray  = malloc(numPositionsInRange * sizeof(uint64_t) * 2);
  uint64_t *const restrict offsetArray    = positionArray + numPositionsInRange;
  //check for allocation failures
  if(__builtin_expect(positionArray == NULL, 0)){
    *fileAccessResult = AwFmAllocationFailure;
    return NULL;
  }

  //call a prefetch for each block that contains the positions that we need to start querying
  for(uint64_t i = searchRange->startPtr; i < searchRange->endPtr; i += POSITIONS_PER_FM_BLOCK){
    awFmOccurancesDataPrefetch(index, i);
  }

  //backtrace each position until we have a list of the positions in the database sequence.
  for(uint64_t i=0; i < numPositionsInRange; i++){
    uint64_t databaseSequenceOffset = 0;
    uint64_t backtracePosition      = searchRange->startPtr + i;

    while(!awFmBwtPositionIsSampled(index, backtracePosition)){
      backtracePosition = awFmBackstepBwtPosition(index, backtracePosition);
      databaseSequenceOffset++;
    }

    positionArray[i]  = backtracePosition;
    offsetArray[i]    = databaseSequenceOffset;
  }

  *fileAccessResult = awFmDbSequencePositionsFromSuffixArrayFile(index, positionArray,
    offsetArray, numPositionsInRange);

  return positionArray;
}



/*
 * Function:  awFmDatabaseSingleKmerExactMatch
 * --------------------
 *  Queries the FM-Index for the range of BWT positions that represent instances
 *    of the given Kmer found in the database.
 *
 *  If the given kmer is not found, the AwFmSearch Range will result in a false value
 *   when given to the awFmSearchRangeIsValid function.
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur, including
 *      creating potential segfauts.
 *    kmerLength:   Length of the kmer to be queried. Undefined behavior may occur if
 *      the function is given a kmerLength of 0.
 *
 *  Returns:
 *    AwFmSearch range representing the range of BWT positions where the given
 *    kmer may be found, as long as startPtr < endPtr. Otherwise (startPtr >= endPtr),
 *    the given kmer does not exist in the database sequence.
 */
struct AwFmSearchRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
const char *restrict const kmer, const uint16_t kmerLength){
  uint_fast16_t kmerQueryPosition = kmerLength - 1;

  //set the initial range from the rankPrefixSums
  const uint8_t firstSuffixLetter = kmer[kmerQueryPosition];
  struct AwFmSearchRange searchRange = {index->rankPrefixSums[firstSuffixLetter], index->rankPrefixSums[firstSuffixLetter + 1]};

  //request the data be prefetched for each occurance call.
  awFmOccurancesDataPrefetch(index, searchRange.startPtr);
  awFmOccurancesDataPrefetch(index, searchRange.endPtr);
  do{
    kmerQueryPosition--;
    searchRange = AwFmIterativeRangeBackwardSearch(index, &searchRange, kmer[kmerQueryPosition]);
  }while(__builtin_expect(awFmSearchRangeIsValid(&searchRange) && kmerQueryPosition, 1));

  return searchRange;
}


/*
 * Function:  awFmSingleKmerExists
 * --------------------
 *  Queries the FM-Index to determine if the database sequence contains any instances of the given kmer
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur, including
 *      creating potential segfauts.
 *    kmerLength:   Length of the kmer to be queried. Undefined behavior may occur if
 *      the function is given a kmerLength of 0.
 *
 *  Returns:
 *    True if the kmer exists in the database sequence, or false if it
 *      cannot be found in the database sequence.
 */
bool awFmSingleKmerExists(const struct AwFmIndex *restrict const index, const char *restrict const kmer,
  const uint16_t kmerLength){
    struct AwFmSearchRange kmerRange = awFmDatabaseSingleKmerExactMatch(index, kmer, kmerLength);
    return kmerRange.startPtr < kmerRange.endPtr;
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
