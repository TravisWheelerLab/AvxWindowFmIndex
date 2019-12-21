#ifndef AW_FM_PARALLEL_SEARCH_H
#define AW_FM_PARALLEL_SEARCH_H

#include "AwFmIndex.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>


struct AwFmParallelSearchMetadata{
bool  useOmpMultithreading;
uint16_t  numThreads;
uint16_t  numConcurrentQueries;
};

struct AwFmKmer{
  uint16_t length;
  char *string;
};

struct AwFmParallelSearchData{
  struct  AwFmParallelSearchMetadata  metadata;
  struct  AwFmKmer                    *kmerList;
  struct  AwFmVector                  *sequencePositionLists;
          size_t                      capacity;
          size_t                      count;
};

struct AwFmBacktraceData{
  size_t position;
  size_t backtraceOffset;
};


struct AwFmParallelSearchData *awFmCreateParallelSearchData(const size_t capacity,
  const struct AwFmParallelSearchMetadata *restrict const metadata);

void awFmDeallocParallelSearchData(struct AwFmParallelSearchData *restrict const searchData);

void awFmSearchDataAddKmer(struct AwFmParallelSearchData *restrict const searchData,
  char *restrict const kmer, const uint16_t kmerLength);

enum AwFmReturnCode awFmParallelSearch(const struct AwFmIndex *restrict const index,
  const struct AwFmParallelSearchData *restrict searchData);

#endif /* end of include guard: AW_FM_PARALLEL_SEARCH_H */
