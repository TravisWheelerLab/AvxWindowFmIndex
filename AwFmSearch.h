#ifndef AW_FM_INDEX_SEARCH_H
#define AW_FM_INDEX_SEARCH_H

#include "AwFmIndex.h"
#include "AwFmFile.h"
#include <stdint.h>
#include <stdbool.h>

struct AwFmBacktracePosition{
  uint64_t position;
  uint64_t offset;
};

//todo: what is this?
struct AwFmDbPositionArray *awFmDatabaseSearchForKmers(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint16_t kmerLength);

struct AwFmSearchRange awFmSearchKmerSuffix(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const currentRange, const uint8_t queryLetter);

uint64_t *awFmFindDatabaseHitPositionsFromSearchRange(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const searchRange,
  enum AwFmFileAccessCode *restrict fileAccessResult);

struct AwFmSearchRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
const char *restrict const kmer, const uint16_t kmerLength);

bool awFmSingleKmerExists(const struct AwFmIndex *restrict const index, const char *restrict const kmer,
  const uint16_t kmerLength);




#endif /* end of include guard: AW_FM_INDEX_SEARCH_H */
