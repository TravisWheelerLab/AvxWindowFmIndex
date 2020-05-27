#ifndef AW_FM_INDEX_SEARCH_H
#define AW_FM_INDEX_SEARCH_H

#include "AwFmFile.h"
#include "AwFmIndexStruct.h"
#include <stdbool.h>
#include <stdint.h>


/*
 * Function:  awFmIterativeStepBackwardSearch
 * --------------------
 * Performs a single backward search step on the given index.
 *  In lieu of returning an additional value, this function updates the data pointed to by
 *  the range ptr.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    range: range in the BWT that corresponds to the implicit kmer that is about to be extended.
 *      this acts as an out-parameter, and will update to the newly extended range once finished.
 *    letter: letter of the prefix or suffix character.
 */
void awFmNucleotideIterativeStepBackwardSearch(const struct AwFmIndex *restrict const index,
  struct AwFmSearchRange *restrict const range, const uint8_t letter);

/*
 * Function:  awFmIterativeStepBackwardSearch
 * --------------------
 * Performs a single backward search step on the given index.
 *  In lieu of returning an additional value, this function updates the data pointed to by
 *  the range ptr.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    range: range in the BWT that corresponds to the implicit kmer that is about to be extended.
 *      this acts as an out-parameter, and will update to the newly extended range once finished.
 *    letter: letter of the prefix or suffix character.
 */
void awFmAminoIterativeStepBackwardSearch(const struct AwFmIndex *restrict const index,
  struct AwFmSearchRange *restrict const range, const uint8_t letter);


uint64_t awFmNucleotideBackwardSearchSingle(const struct AwFmIndex *restrict const index,
  uint64_t queryPosition, const uint8_t letter);


uint64_t awFmAminoBackwardSearchSingle(const struct AwFmIndex *restrict const index,
  uint64_t queryPosition, const uint8_t letter);


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
   const struct AwFmSearchRange *restrict const searchRange, enum AwFmReturnCode *restrict fileAccessResult);


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
  *    AwFmSearchRange representing the range of BWT positions where the given
  *    kmer may be found, as long as startPtr < endPtr. Otherwise (startPtr >= endPtr),
  *    the given kmer does not exist in the database sequence.
  */
struct AwFmSearchRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint16_t kmerLength);


/*
 * Function:  awFmBacktraceBwtPosition
 * --------------------
 * Given a specified Bwt position, backsteps to find the position one before in
 *  original sequence.
 *
 *  Inputs:
 *    index: Index to backstep
 *    alphabet: alphabet of the index, either Nucleotide or Amino
 *    bwtPosition: Position of the character to be returned.
 *
 *  Returns:
 *    Position in the suffix array of the character in the sequence immediately preceeding the one
 *      found at the given bwtPosition.
 */
size_t awFmNucleotideBacktraceBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition);

/*
 * Function:  awFmBacktraceBwtPosition
 * --------------------
 * Given a specified Bwt position, backsteps to find the position one before in
 *  original sequence.
 *
 *  Inputs:
 *    index: Index to backstep
 *    alphabet: alphabet of the index, either Nucleotide or Amino
 *    bwtPosition: Position of the character to be returned.
 *
 *  Returns:
 *    Position in the suffix array of the character in the sequence immediately preceeding the one
 *      found at the given bwtPosition.
 */
size_t awFmAminoBacktraceBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition);

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
 * Function:  awFmNucleotideNonSeededSearch
 * --------------------
 *  Finds the range in the bwt corresponding to the entire given nucleotide kmer. While this
 *    function will function on any input kmer, this function explicitly exists
 *    to query for kmers that are too short to be memoized in the kmerTable.
 *    If the kmer you are searching for is at least as long as those in the kmerTable,
 *    do not use this function. Instead, start with the memozied range for the kmer suffix
 *    and extend from there.
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur, including
 *      creating potential segfauts.
 *    kmerLength:   Length of the kmer to be queried. Undefined behavior may occur if
 *      the function is given a kmerLength of 0.
 *    range:        Pointer to a search range where the output will be written.
 *      Since gcc doesn't seem to perform NRVO correctly, this performs slightly better.
 *
 */
void awFmNucleotideNonSeededSearch(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength, struct AwFmSearchRange *range);


/*
 * Function:  awFmAminoNonSeededSearch
 * --------------------
 *  Finds the range in the bwt corresponding to the entire given amino kmer. While this
 *    function will function on any input kmer, this function explicitly exists
 *    to query for kmers that are too short to be memoized in the kmerTable.
 *    If the kmer you are searching for is at least as long as those in the kmerTable,
 *    do not use this function. Instead, start with the memozied range for the kmer suffix
 *    and extend from there.
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur, including
 *      creating potential segfauts.
 *    kmerLength:   Length of the kmer to be queried. Undefined behavior may occur if
 *      the function is given a kmerLength of 0.
 *    range:        Pointer to a search range where the output will be written.
 *      Since gcc doesn't seem to perform NRVO correctly, this performs slightly better.
 *
 */
void awFmAminoNonSeededSearch(const struct AwFmIndex *restrict const index,
 const char *restrict const kmer, const uint8_t kmerLength, struct AwFmSearchRange *range);


#endif /* end of include guard: AW_FM_INDEX_SEARCH_H */
