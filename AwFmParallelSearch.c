#include "AwFmParallelSearch.h"
#include "AwFmVector.h"

#include <string.h>

//TODO: use pread for thread safe read access to suffix array

#define DEFAULT_POSITION_LIST_CAPACITY  256


struct AwFmParallelSearchData *awFmCreateParallelSearchData(const size_t capacity,
  const struct AwFmParallelSearchMetadata *restrict const metadata){

  struct AwFmParallelSearchData *searchData = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
                                                      sizeof(struct AwFmParallelSearchData));
  if(searchData == NULL){
    return NULL;
  }

  searchData->capacity  = capacity;
  searchData->count     = 0;
  memcpy(&searchData->metadata, metadata, sizeof(struct AwFmParallelSearchMetadata));

  searchData->kmerList = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES, capacity * sizeof(struct AwFmKmer));
  if(searchData->kmerList == NULL){
    free(searchData);
    return NULL;
  }
  searchData->sequencePositionLists = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
                                                capacity * sizeof(struct AwFmVector));

  if(searchData->sequencePositionLists == NULL){
    free(searchData->kmerList);
    free(searchData);
    return NULL;
  }

  //initialize all the elements in the positionLists to obviously null values
  for(size_t i = 0; i < searchData->capacity; i++){
    awFmVectorCreate(DEFAULT_POSITION_LIST_CAPACITY, sizeof(uint64_t), searchData->sequencePositionLists);
  }

  return searchData;
}


void awFmDeallocParallelSearchData(struct AwFmParallelSearchData *restrict const searchData){
  free(searchData->kmerList);

  for(size_t i = 0; i < searchData->capacity; i++){
    awFmVectorDealloc(&searchData->sequencePositionLists[i]);
    free(searchData->sequencePositionLists);
  }

  free(searchData);
}


void awFmSearchDataAddKmer(struct AwFmParallelSearchData *restrict const searchData,
  char *restrict const kmer, const uint16_t kmerLength){
  searchData->kmerList[searchData->count].length = kmerLength;
  searchData->kmerList[searchData->count].string = kmer;
  searchData->count++;
}

enum AwFmReturnCode awFmParallelSearch(const struct AwFmIndex *restrict const index,
  const struct AwFmParallelSearchData *restrict searchData){

  //make local copies of the control data to encourage them to be used in a thread-shared manner.
  const size_t    searchDataCount       = searchData->count;
  const bool      useOmpMultithreading  = searchData->metadata.useOmpMultithreading;
  const uint16_t  numThreads            = searchData->metadata.numThreads;
  const uint16_t  numConcurrentQueries  = searchData->metadata.numConcurrentQueries;

  #pragma omp parallel for num_threads(numThreads) schedule(dynamic) if(useOmpMultithreading)
  for(size_t threadBlockStartIndex = 0;
    threadBlockStartIndex < searchDataCount;
    threadBlockStartIndex += numConcurrentQueries){
      size_t threadBlockEndIndex = threadBlockStartIndex + numConcurrentQueries;

      struct AwFmBackwardRange ranges[numConcurrentQueries];
      struct AwFmBacktraceData backtraces[numConcurrentQueries];
      //make a struct AwFmBacktraceData?

      //make struct for searchData, ranges, threadBlockStartIndex, numConcurrentQueries?
      parallelSearchFindKmerSeedsForBlock(searchData, ranges, threadBlockStartIndex, numConcurrentQueries);

      parallelSearchExtendKmersInBlock(searchData, ranges, threadBlockStartIndex, numConcurrentQueries);

      //these position lists can serve as the temporary home for the SA positions
      parallelSearchCreatePositionLists(searchData, ranges, threadBlockStartIndex, numConcurrentQueries);


      parallelSearchTracebackPositionLists(searchData, ranges, threadBlockEndIndex, numConcurrentQueries);

      //during backtrace, we're backtracing 'numConcurrentQueries' positions
      //when one finishes, replace the 'ptr to the one we're backtracing' with the next one that needs to be
      // backtraced




      //for each, dynamically allocate the integer arrays, set in searchData
      //perform all the backtraces. here, this can be concurrent across the block again.
      //we might need 'offsets' here again
      //it might be a good idea to use {} to scope the ranges to make them fall off stack
      //so we can reuse the stack memory for offsets

      //once all the SA positions are found,
      //lock the threadHandle,
      //for each hit, seek into SA file, read each position
      //unlock threadHandle

    }


  }
