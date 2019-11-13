#ifndef AW_FM_INDEX_SEARCH_H
#define AW_FM_INDEX_SEARCH_H

#include "AwFmFile.h"
#include "AwFmIndex.h"
#include <stdbool.h>
#include <stdint.h>


//


/*Given the range of BWT positions where current kmer suffix is found, returns the
    range where the suffix can be found after prepending the next prefix letter.*/
struct AwFmSearchRange AwFmIterativeAminoAcidRangeBackwardSearch(const struct AwFmAminoBlock *restrict const blockList,
  const uint64_t *restrict const rankPrefixSums, const struct AwFmSearchRange *restrict const currentRange,
  const uint8_t prefixLetter);

/*Given the range of BWT positions where current kmer suffix is found, returns the
    range where the suffix can be found after appending the next suffix letter.*/
struct AwFmSearchRange AwFmIterativeRangeForwardSearch(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const currentRange, const uint8_t queryLetter);

/*Returns an array of positions in the database sequence that are represented by the
    given searchRange. The length of the returned array is equal to the difference
    between the searchRange pointers. It is the caller's responsibility to free the
    array that is returned.*/
uint64_t *awFmFindAminoAcidDatabaseHitPositionsForward(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const searchRange,
  enum AwFmReturnCode *restrict fileAccessResult);

/*Returns an array of positions in the database sequence that are represented by the
    given searchRange. The length of the returned array is equal to the difference
    between the searchRange pointers. It is the caller's responsibility to free the
    array that is returned.*/
uint64_t *awFmFindAminoAcidDatabaseHitPositionsBackward(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const searchRange,
  enum AwFmReturnCode *restrict fileAccessResult);

/*queries the AwFmIndex, returning a range of positions in the BWT where the given kmer is found.*/
struct AwFmSearchRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint16_t kmerLength);

/*queries the AwFmIndex, returning true if the given kmer is found in the database sequence.*/
bool awFmSingleKmerExists(const struct AwFmIndex *restrict const index, const char *restrict const kmer,
  const uint16_t kmerLength);

size_t awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range);

#endif /* end of include guard: AW_FM_INDEX_SEARCH_H */
