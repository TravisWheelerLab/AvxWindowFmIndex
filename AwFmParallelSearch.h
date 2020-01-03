#ifndef AW_FM_PARALLEL_SEARCH_H
#define AW_FM_PARALLEL_SEARCH_H

#include "AwFmIndex.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>


struct AwFmKmer{
  uint16_t  length;
  char      *string;
};

struct AwFmParallelSearchData{
  struct  AwFmKmer            *kmerList;
  struct  AwFmBacktraceVector *sequencePositionLists;
          size_t              capacity;
          size_t              count;
          uint_fast16_t       numThreads;
};


/*
 * Function:  awFmCreateParallelSearchData
 * --------------------
 *  Allocates and initializes an AwFmParallelSearchData struct to be used to search
 *  for groups of kmers in a thread-parallel manner that also hides memory read latency
 *  through multiple concurrent queries per thread.
 *
 *  Note that the kmers inside the searchData are not allocated, and only contain
 *  char pointers that can be set to the kmers you want to query for.
 *
 *  Inputs:
 *    capacity:     How many kmers the searchData struct can hold.
 *    numThreads:   number of threads to use to query the searchData.
 *
 *  Returns:
 *    Pointer to the allocated searchData struct, or null on failure.
 */
struct AwFmParallelSearchData *awFmCreateParallelSearchData(const size_t capacity,
  const uint_fast8_t numThreads);


/*
 * Function:  awFmDeallocParallelSearchData
 * --------------------
 *  Deallocates the given search data, and all internally stored dynamically allocated memory.
 *
 *  Note that, since the searchData doesn't own the kmer char strings, it will not try to
 *  deallocate them. If kmers were dynamically allocated externally, it is the caller's responsibility
 *  to deallocate them. This means that if these pointers are the only pointer to the data,
 *  and if they were dynamically allocated, forgetting to deallocate them before calling this function
 *  will leak the data.
 *
 *  Inputs:
 *    searchData:   pointer to the searchData struct to deallocate
 */
void awFmDeallocParallelSearchData(struct AwFmParallelSearchData *restrict const searchData);


/*
 * Function:  awFmParallelSearch
 * --------------------
 *  Using the given index and a searchData struct preloaded with kmers, query the kmers
 *  in a concurrent, thread-parallel manner. The suggested use case for this function is as follows:
 *
 *    1. allocate a searchData struct with awFmCreateParallelSearchData().
 *    2. for searching for n kmers, set count to n (must be smaller than capacity!), and
 *      set the first n kmers with the correct char pointer and length.
 *    3. call this function using the index to search.
 *    4. use the database sequence positions returned in the sequencePositionLists member.
 *    5. to query for additional kmers, reuse the searchData struct, starting with step (2).
 *    6. deallocate the searchData with awFmDeallocParallelSearchData when finished.
 *
 *
 *  Inputs:
 *    index:        pointer to the index to search.
 *    searchData:   pointer to the searchData loaded with kmers to search for.
 */
void awFmParallelSearch(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData);

#endif /* end of include guard: AW_FM_PARALLEL_SEARCH_H */
