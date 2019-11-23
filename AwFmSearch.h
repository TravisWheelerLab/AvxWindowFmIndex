#ifndef AW_FM_INDEX_SEARCH_H
#define AW_FM_INDEX_SEARCH_H

#include "AwFmFile.h"
#include "AwFmIndex.h"
#include <stdbool.h>
#include <stdint.h>



struct AwFmBackwardRange{
  uint64_t startPtr;
  uint64_t endPtr;
};

struct AwFmBiDirectionalRange{
  uint64_t startPtr;
  uint64_t endPtr;
  uint64_t startPrimePtr;
};



void AwFmIterativeStepBidirectionalNucleotideSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBiDirectionalRange *restrict const range, const uint8_t letter);

void AwFmIterativeStepBidirectionalAminoAcidSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBiDirectionalRange *restrict const range, const uint8_t letter);

void AwFmIterativeStepBackwardNucleotideSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBackwardRange *restrict const range, const uint8_t letter);

void AwFmIterativeStepBackwardAminoAcidSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBackwardRange *restrict const range, const uint8_t letter);


/*Returns an array of positions in the database sequence that are represented by the
    given searchRange. The length of the returned array is equal to the difference
    between the searchRange pointers. It is the caller's responsibility to free the
    array that is returned.*/
uint64_t *awFmFindDatabaseHitPositions(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const searchRange,
  enum AwFmReturnCode *restrict fileAccessResult);

/*queries the AwFmIndex, returning a range of positions in the BWT where the given kmer is found.*/
struct AwFmBackwardRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint16_t kmerLength);

/*queries the AwFmIndex, returning true if the given kmer is found in the database sequence.*/
bool awFmSingleKmerExists(const struct AwFmIndex *restrict const index, const char *restrict const kmer,
  const uint16_t kmerLength);

size_t awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range);

#endif /* end of include guard: AW_FM_INDEX_SEARCH_H */
