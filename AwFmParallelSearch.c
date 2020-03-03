#include "AwFmParallelSearch.h"
#include "AwFmSearch.h"
#include "AwFmKmerTable.h"
#include "AwFmBacktraceVector.h"
#include "AwFmFile.h"
#include "AwFmLetter.h"

#include <string.h>


#define NUM_CONCURRENT_QUERIES          32


void parallelSearchFindKmerSeedsForBlock(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

void parallelSearchExtendKmersInBlock(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

void parallelSearchTracebackPositionLists(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

void parallelSearchSuffixArrayLookup(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);



struct AwFmParallelSearchData *awFmCreateParallelSearchData(const size_t capacity,
  const uint_fast8_t numThreads){

  struct AwFmParallelSearchData *searchData = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
                                                      sizeof(struct AwFmParallelSearchData));
  if(searchData == NULL){
    return NULL;
  }

  searchData->capacity    = capacity;
  searchData->count       = 0;
  searchData->numThreads  = numThreads;

  searchData->kmerList = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES, capacity * sizeof(struct AwFmKmer));
  if(searchData->kmerList == NULL){
    free(searchData);
    return NULL;
  }

  searchData->sequencePositionLists = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
                                                capacity * sizeof(struct AwFmBacktraceVector));

  if(searchData->sequencePositionLists == NULL){
    free(searchData->kmerList);
    free(searchData);
    return NULL;
  }

  //initialize all the elements in the positionLists to obviously null values
  for(size_t i = 0; i < searchData->capacity; i++){
    awFmBacktraceVectorCreate(&searchData->sequencePositionLists[i]);
  }

  return searchData;
}


void awFmDeallocParallelSearchData(struct AwFmParallelSearchData *restrict const searchData){
  free(searchData->kmerList);

  for(size_t i = 0; i < searchData->capacity; i++){
    awFmBacktraceVectorDealloc(&searchData->sequencePositionLists[i]);
  }

  free(searchData->sequencePositionLists);
  free(searchData);
}



void awFmParallelSearch(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData){
  //make local copies of the control data to encourage them to be used in a thread-shared manner.
  const size_t searchDataCount = searchData->count;

  #pragma omp parallel for num_threads(8) schedule(dynamic)
  for(size_t threadBlockStartIndex = 0; threadBlockStartIndex < searchDataCount; threadBlockStartIndex += NUM_CONCURRENT_QUERIES){
    const size_t threadBlockEndIndex = threadBlockStartIndex + NUM_CONCURRENT_QUERIES > searchData->count?
      searchData->count: threadBlockStartIndex + NUM_CONCURRENT_QUERIES;
    // printf("from tbsi %zu\n", threadBlockStartIndex);
    struct AwFmSearchRange ranges[NUM_CONCURRENT_QUERIES];

    parallelSearchFindKmerSeedsForBlock(index, searchData, ranges, threadBlockStartIndex, threadBlockEndIndex);


    parallelSearchExtendKmersInBlock(index, searchData, ranges, threadBlockStartIndex, threadBlockEndIndex);

    parallelSearchTracebackPositionLists(index, searchData, ranges, threadBlockStartIndex, threadBlockEndIndex);

    parallelSearchSuffixArrayLookup(index, searchData, threadBlockStartIndex, threadBlockEndIndex);
  }
}



/*private function prototypes*/
void parallelSearchFindKmerSeedsForBlock(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){

  for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
    const uint8_t kmerLength  = searchData->kmerList[kmerIndex].length;
    const char    *kmerString = searchData->kmerList[kmerIndex].string;

    const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;
    if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
      ranges[rangesIndex] = awFmNucleotideKmerSeedRangeFromTable(index, kmerString, kmerLength);
    }
    else{
      ranges[rangesIndex] = awFmAminoKmerSeedRangeFromTable(index, kmerString, kmerLength);
    }
  }
}


void parallelSearchExtendKmersInBlock(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){

  bool hasActiveQueries = true;
  uint8_t currentKmerLetterIndex = index->metadata.kmerLengthInSeedTable;

  while(hasActiveQueries){
    currentKmerLetterIndex++;
    hasActiveQueries = false;

    for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
      const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;
      const struct AwFmKmer *kmerPtr = &searchData->kmerList[kmerIndex];

      if((kmerPtr->length >= currentKmerLetterIndex) && awFmSearchRangeIsValid(&ranges[rangesIndex])){
        hasActiveQueries = true;
        const uint8_t currentQueryLetterIndex = kmerPtr->length - currentKmerLetterIndex;
        // const char    currentQueryLetter      = kmerPtr->string[currentQueryLetterIndex];
        const uint8_t queryLetterIndex = awFmAsciiNucleotideToLetterIndex(kmerPtr->string[currentQueryLetterIndex]);

        if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
          const uint8_t queryLetterIndex = awFmAsciiNucleotideToLetterIndex(kmerPtr->string[currentQueryLetterIndex]);
          awFmNucleotideIterativeStepBackwardSearch(index, &ranges[rangesIndex], queryLetterIndex);
        }
        else{
          const uint8_t queryLetterIndex = awFmAsciiAminoAcidToLetterIndex(kmerPtr->string[currentQueryLetterIndex]);
          awFmAminoIterativeStepBackwardSearch(index, &ranges[rangesIndex], queryLetterIndex);
        }
      }
    }
  }
}



void parallelSearchTracebackPositionLists(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){

  for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
  const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;
    const size_t rangeLength  = awFmSearchRangeLength(&ranges[rangesIndex]);
    struct AwFmBacktrace *restrict const backtraceArray = searchData->sequencePositionLists[kmerIndex].backtraceArray;
    awFmBacktraceVectorSetCount(&searchData->sequencePositionLists[kmerIndex], rangeLength);

    size_t positionInRangeToBacktrace = 0;
    while(positionInRangeToBacktrace < rangeLength){
      uint64_t position = ranges[rangesIndex].startPtr + positionInRangeToBacktrace;
      uint64_t offset   = 0;

      if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
        while(!awFmBwtPositionIsSampled(index, position)){
          if(__builtin_expect(position == index->sentinelCharacterPosition, 0)){
            position = 0;
          }else{
            position = awFmNucleotideBacktraceBwtPosition(index, position);
          }
          offset++;
        }
      }
      else{
        while(!awFmBwtPositionIsSampled(index, position)){
          if(__builtin_expect(position == index->sentinelCharacterPosition, 0)){
            position = 0;
          }else{
            position = awFmAminoBacktraceBwtPosition(index, position);
          }
          offset++;
        }
      }
      backtraceArray[positionInRangeToBacktrace].position = position;
      backtraceArray[positionInRangeToBacktrace].offset   = offset;
      positionInRangeToBacktrace++;
    }
  }
}


void parallelSearchSuffixArrayLookup(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){
    for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){

    const struct AwFmBacktraceVector *backtraceVectorPtr = &searchData->sequencePositionLists[kmerIndex];
    const size_t numSuffixArrayPositions = backtraceVectorPtr->count;
    for(size_t backtraceIndex = 0; backtraceIndex < numSuffixArrayPositions; backtraceIndex++){
      struct AwFmBacktrace *restrict const backtracePtr = &backtraceVectorPtr->backtraceArray[backtraceIndex];
      awFmSuffixArrayReadPositionParallel(index, backtracePtr);
    }
  }
}
