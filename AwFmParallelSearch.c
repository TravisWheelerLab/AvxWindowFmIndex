#include "AwFmParallelSearch.h"
#include "AwFmSearch.h"
#include "AwFmKmerTable.h"
#include "AwFmBacktraceVector.h"
#include "AwFmFile.h"

#include <string.h>


#define NUM_CONCURRENT_QUERIES          8


void parallelSearchFindKmerSeedsForBlock(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex);

void parallelSearchExtendKmersInBlock(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex);

void parallelSearchTracebackPositionLists(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex);

void parallelSearchSuffixArrayLookup(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, const size_t threadBlockStartIndex);



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
    awFmBacktraceVectorCreate(searchData->sequencePositionLists);
  }

  return searchData;
}


void awFmDeallocParallelSearchData(struct AwFmParallelSearchData *restrict const searchData){
  free(searchData->kmerList);

  for(size_t i = 0; i < searchData->capacity; i++){
    awFmBacktraceVectorDealloc(&searchData->sequencePositionLists[i]);
    free(searchData->sequencePositionLists);
  }

  free(searchData);
}



void awFmParallelSearch(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData){
  //make local copies of the control data to encourage them to be used in a thread-shared manner.
  const size_t searchDataCount = searchData->count;

  #pragma omp parallel for num_threads(searchData->numThreads) schedule(dynamic)
  for(size_t threadBlockStartIndex = 0;
    threadBlockStartIndex < searchDataCount;
    threadBlockStartIndex += NUM_CONCURRENT_QUERIES){
      struct AwFmSearchRange ranges[NUM_CONCURRENT_QUERIES];

      parallelSearchFindKmerSeedsForBlock(index, searchData, ranges, threadBlockStartIndex);

      parallelSearchExtendKmersInBlock(index, searchData, ranges, threadBlockStartIndex);

      parallelSearchTracebackPositionLists(index, searchData, ranges, threadBlockStartIndex);

      parallelSearchSuffixArrayLookup(index, searchData, threadBlockStartIndex);
    }
  }



/*private function prototypes*/
void parallelSearchFindKmerSeedsForBlock(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex){

  for(size_t i = 0; i < NUM_CONCURRENT_QUERIES; i++){
    const size_t  kmerIndex   = i + threadBlockStartIndex;
    const uint8_t kmerLength  = searchData->kmerList[kmerIndex].length;
    const char    *kmerString = searchData->kmerList[kmerIndex].string;
    
    const struct AwFmSearchRange range = index->metadata.alphabetType == AwFmAlphabetNucleotide?
    awFmNucleotideSeedKmerRangeFromTable(index, kmerString, kmerLength):
    awFmAminoSeedKmerRangeFromTable(index, kmerString, kmerLength);
    memcpy(&ranges[i], &range, sizeof(struct AwFmSearchRange));
  }
}


void parallelSearchExtendKmersInBlock(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex){

  bool hasActiveQueries = true;
  uint8_t currentKmerLetterIndex = index->metadata.kmerLengthInSeedTable;

  while(hasActiveQueries){
    currentKmerLetterIndex++;
    hasActiveQueries = false;

    for(uint8_t i = 0; i < NUM_CONCURRENT_QUERIES; i++){
      const size_t kmerIndex = i + threadBlockStartIndex;
      const struct AwFmKmer *kmerPtr = &searchData->kmerList[kmerIndex];

      if((kmerPtr->length >= currentKmerLetterIndex) && awFmSearchRangeIsValid(&ranges[i])){
        hasActiveQueries = true;
        const uint8_t currentQueryLetterIndex = kmerPtr->length - currentKmerLetterIndex;
        const char    currentQueryLetter      = kmerPtr->string[currentQueryLetterIndex];

        if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
          awFmNucleotideIterativeStepBackwardSearch(index, &ranges[i], currentQueryLetter);
        }
        else{
          awFmAminoIterativeStepBackwardSearch(index, &ranges[i], currentQueryLetter);
        }
      }
    }
  }
}



void parallelSearchTracebackPositionLists(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex){

  for(uint8_t i = 0; i < NUM_CONCURRENT_QUERIES; i++){
    const size_t kmerIndex    = i + threadBlockStartIndex;
    const size_t rangeLength  = awFmSearchRangeLength(&ranges[i]);
    struct AwFmBacktrace *restrict const backtraceArray = searchData->sequencePositionLists[kmerIndex].backtraceArray;
    awFmBacktraceVectorSetCount(&searchData->sequencePositionLists[kmerIndex], rangeLength);

    size_t positionInRangeToBacktrace = 0;
    while(positionInRangeToBacktrace < rangeLength){
      uint64_t position = ranges[i].startPtr + positionInRangeToBacktrace;
      uint64_t offset   = 0;

      if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
        while(!awFmBwtPositionIsSampled(index, position)){
          position = awFmNucleotideBacktraceBwtPosition(index, position);
          offset++;
        }
      }
      else{
        while(!awFmBwtPositionIsSampled(index, position)){
          position = awFmAminoBacktraceBwtPosition(index, position);
          offset++;
        }
      }
      backtraceArray[positionInRangeToBacktrace].position = position;
      backtraceArray[positionInRangeToBacktrace].offset   = offset;
    }
  }
}


void parallelSearchSuffixArrayLookup(const struct AwFmIndex *restrict const index,
  struct AwFmParallelSearchData *restrict const searchData, const size_t threadBlockStartIndex){
  for(uint8_t i = threadBlockStartIndex; i < (threadBlockStartIndex + NUM_CONCURRENT_QUERIES); i++){

    const struct AwFmBacktraceVector *backtraceVectorPtr = &searchData->sequencePositionLists[i];
    const size_t numSuffixArrayPositions = backtraceVectorPtr->count;
    for(size_t backtraceIndex = 0; backtraceIndex < numSuffixArrayPositions; backtraceIndex++){

      struct AwFmBacktrace *restrict const backtracePtr = &backtraceVectorPtr->backtraceArray[backtraceIndex];
      awFmSuffixArrayReadPositionParallel(index, backtracePtr);
    }
  }
}
