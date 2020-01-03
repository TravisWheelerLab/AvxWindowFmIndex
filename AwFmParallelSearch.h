#ifndef AW_FM_PARALLEL_SEARCH_H
#define AW_FM_PARALLEL_SEARCH_H

#include "AwFmIndex.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>


struct AwFmKmer{
  uint16_t length;
  char *string;
};

struct AwFmParallelSearchData{
  struct  AwFmKmer            *kmerList;
  struct  AwFmBacktraceVector *sequencePositionLists;
          size_t              capacity;
          size_t              count;
          uint_fast16_t       numThreads;
};


struct AwFmParallelSearchData *awFmCreateParallelSearchData(const size_t capacity,
  const uint_fast8_t numThreads);

void awFmDeallocParallelSearchData(struct AwFmParallelSearchData *restrict const searchData);

void awFmParallelSearch(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData);

#endif /* end of include guard: AW_FM_PARALLEL_SEARCH_H */
