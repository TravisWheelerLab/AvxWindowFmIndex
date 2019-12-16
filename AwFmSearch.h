#ifndef AW_FM_INDEX_SEARCH_H
#define AW_FM_INDEX_SEARCH_H

#include "AwFmFile.h"
#include "AwFmIndex.h"
#include <stdbool.h>
#include <stdint.h>


/*
 * Function:  awFmIterativeStepBiDirectionalSearch
 * --------------------
 * Performs a single bi-directional search step on the given index. This can either search
 *  backward or forward, depending on the searchDirection argument.
 *  In lieu of returning an additional value, this function updates the data pointed to by
 *  the range ptr.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    searchDirection: determines whether to search backwards (prepending a prefix character), or
 *    forward (appending a suffix character)
 *    range: range in the BWT that corresponds to the implicit kmer that is about to be extended.
 *      this acts as an out-parameter, and will update to the newly extended range once finished.
 *    letter: letter of the prefix or suffix character.
 *
 */
void awFmIterativeStepBiDirectionalSearch(const struct AwFmIndex *restrict const index,
  const enum AwFmSearchDirection searchDirection, struct AwFmBiDirectionalRange *restrict const range,
  const uint8_t letter);


/*
 * Function:  awFmIterativeStepBackwardSearch
 * --------------------
 * Performs a single backwared search step on the given index.
 *  In lieu of returning an additional value, this function updates the data pointed to by
 *  the range ptr.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    range: range in the BWT that corresponds to the implicit kmer that is about to be extended.
 *      this acts as an out-parameter, and will update to the newly extended range once finished.
 *    letter: letter of the prefix or suffix character.
 */
void awFmIterativeStepBackwardSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBackwardRange *restrict const range, const uint8_t letter);


/*
 * Function:  awFmSeedKmerRangeFromTable
 * --------------------
 * Given a kmer, queries the kmerSeedTable  for the partially completed range in the backward suffix array.
 * The length of the kmer's stored in the seed table can be found in the index metadata's kmerLengthInSeedTable attribute.
 * The given kmer must be at least this long.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: character string to search for in the kmerSeedTable.
 *
 *  Returns:
 *    pointer to the backwardsRange struct in the kmerSeedTable.
 */
struct AwFmBackwardRange *awFmSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer);


/*
 * Function:  awFmFindDatabaseHitPositions
 * --------------------
 * Returns an array of positions in the database sequence that are represented by the given searchRange.
 * This function returns a dynamically allocated array, and deallocation is the caller's responsiblity.
 *
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    searchRange: range of positions in the backward BWT. Each valid element in this range
 *      will be represented in the returned array.
 *    fileAccessResult: out parameter returning the success or failure of the function.
 *      Possible returns are:
 *        AwFmGeneralFailure on search range with no valid positions (start > end)
 *        AwFmAllocationFailure on failure to allocate additional memory.
 *        AwFmFileReadOkay on success
 *
 *  Returns:
 *    Dynamically allocated array of sequence positions. The length of this array
 *      is the number of valid positions in the search range ((startPtr - endPtr) + 1).
 *      This function will return NULL if the search range is invalid (if startPtr > endPtr).
 *      It is the caller's responsiblity to deallocate this array.
 */
uint64_t *awFmFindDatabaseHitPositions(const struct AwFmIndex *restrict const index,
  const struct AwFmBackwardRange *restrict const searchRange,
  enum AwFmReturnCode *restrict fileAccessResult);


/*
 * Function:  awFmFindDatabaseHitPositions
 * --------------------
 *  Takes a range of BWT positions, backtraces each position to find the nearest sample in the
 *    compressed suffix array, and looks up those suffix array positions on disk to
 *    determine the corresponding database sequence position for each BWT position
 *    between the searchRange's pointers (inclusive startPtr, exclusive endPtr).
 *
 *  Note: When using a bi-directional FM-index, the range given should correspond the the
 *    traditional backward BWT, not the forward BWT.
 *
 *  It is the caller's responsibility to free() the returned sequence position array and offset array.
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
 uint64_t *awFmFindDatabaseHitPositions(const struct AwFmIndex *restrict const index,
   const struct AwFmBackwardRange *restrict const searchRange, enum AwFmReturnCode *restrict fileAccessResult);


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
struct AwFmBackwardRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint16_t kmerLength);


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
  const uint16_t kmerLength);


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
size_t awFmSearchRangeLength(const struct AwFmBackwardRange *restrict const range);


/*
 * Function:  awFmSwapBiDirectionalRangePointerDirection
 * --------------------
 * Swaps the direction of the bidirectional range struct in order to search in the opposite direction.
 *  This function works in place on the given range, swapping the startPtr and startPrimePtr, and
 *  updating the endPtr to reflect the correct end for the new startPtr.
 *
 *  Inputs:
 *    range: bidirectional range in the Bwt. This act as an out-argument, and will be updated by this function.
 */
void awFmSwapBiDirectionalRangePointerDirection(struct AwFmBiDirectionalRange *restrict const range);

#endif /* end of include guard: AW_FM_INDEX_SEARCH_H */
