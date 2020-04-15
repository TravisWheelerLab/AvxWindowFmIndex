#include "AwFmParallelSearch.h"
#include "AwFmSearch.h"
#include "AwFmKmerTable.h"
#include "AwFmBacktraceVector.h"
#include "AwFmFile.h"
#include "AwFmLetter.h"

#include <string.h>


#define NUM_CONCURRENT_QUERIES          32
#define DEFAULT_POSITION_LIST_CAPACITY  256


void parallelSearchFindKmerSeedsForBlock(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

void parallelSearchExtendKmersInBlock(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

void parallelSearchTracebackPositionLists(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

void parallelSearchSuffixArrayLookup(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

bool setPositionListCount(struct AwFmKmerSearchData *restrict const searchData, uint32_t count);


struct AwFmKmerSearchList *awFmCreateKmerSearchList(const size_t capacity){
  struct AwFmKmerSearchList *searchList = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
    sizeof(struct AwFmKmerSearchList));
    if(searchList == NULL){
      return NULL;
    }

    searchList->capacity = capacity;
    searchList->count = 0;

    searchList->kmerSearchData = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
      capacity * sizeof(struct AwFmKmerSearchData));

    if(searchList->kmerSearchData == NULL){
      free(searchList);
      return NULL;
    }

    bool positionListAllocationFailed = false;
    for(size_t i = 0; i < capacity; i++){
      searchList->kmerSearchData[i].kmerString            = NULL;
      searchList->kmerSearchData[i].kmerLength            = 0;
      searchList->kmerSearchData[i].capacity              = DEFAULT_POSITION_LIST_CAPACITY;
      searchList->kmerSearchData[i].count                 = 0;
      searchList->kmerSearchData[i].positionBacktraceList = aligned_alloc(
        AW_FM_CACHE_LINE_SIZE_IN_BYTES,
        DEFAULT_POSITION_LIST_CAPACITY * sizeof(struct AwFmBacktrace));

      //check for an allocation failure
      positionListAllocationFailed |= (searchList->kmerSearchData[i].positionBacktraceList == NULL);
    }

    //if any of the allocations failed, dealloc everything and return NULL.
    if(positionListAllocationFailed){
      for(size_t i = 0; i < capacity; i++){
        free(searchList->kmerSearchData[i].positionBacktraceList);
      }
      free(searchList->kmerSearchData);
      free(searchList);
      return NULL;
    }

  return searchList;
}


void awFmDeallocKmerSearchList(struct AwFmKmerSearchList *restrict const searchList){
  for(size_t i = 0; i < searchList->capacity; i++){
    free(searchList->kmerSearchData[i].positionBacktraceList);
  }
  free(searchList->kmerSearchData);
  free(searchList);
}


void awFmParallelSearchLocate(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, uint8_t numThreads){

  #pragma omp parallel for num_threads(numThreads)
  for(size_t threadBlockStartIndex = 0; threadBlockStartIndex < searchList->count; threadBlockStartIndex += AW_FM_NUM_CONCURRENT_QUERIES){

    const size_t threadBlockEndIndex = threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES > searchList->count?
      searchList->count: threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES;
    struct AwFmSearchRange ranges[AW_FM_NUM_CONCURRENT_QUERIES];

    parallelSearchFindKmerSeedsForBlock(  index, searchList, ranges,  threadBlockStartIndex, threadBlockEndIndex);
    parallelSearchExtendKmersInBlock(     index, searchList, ranges,  threadBlockStartIndex, threadBlockEndIndex);
    parallelSearchTracebackPositionLists( index, searchList, ranges,  threadBlockStartIndex, threadBlockEndIndex);
    parallelSearchSuffixArrayLookup(      index, searchList,          threadBlockStartIndex, threadBlockEndIndex);
  }
}

void awFmParallelSearchCount(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, uint8_t numThreads){

  #pragma omp parallel for num_threads(numThreads)
  for(size_t threadBlockStartIndex = 0; threadBlockStartIndex < searchList->count; threadBlockStartIndex += AW_FM_NUM_CONCURRENT_QUERIES){

    const size_t threadBlockEndIndex = threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES > searchList->count?
      searchList->count: threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES;
    struct AwFmSearchRange ranges[AW_FM_NUM_CONCURRENT_QUERIES];

    parallelSearchFindKmerSeedsForBlock(  index, searchList, ranges,  threadBlockStartIndex, threadBlockEndIndex);
    parallelSearchExtendKmersInBlock(     index, searchList, ranges,  threadBlockStartIndex, threadBlockEndIndex);

    //load the range lengths into the count member variables.
    for(size_t i = threadBlockStartIndex; i < threadBlockEndIndex; i++){
      searchList->kmerSearchData[i].count = awFmSearchRangeLength(&ranges[i-threadBlockStartIndex]);
    }
  }
}


void parallelSearchFindKmerSeedsForBlock(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){

  for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
    const uint8_t kmerLength  = searchList->kmerSearchData[kmerIndex].kmerLength;
    const char    *kmerString = searchList->kmerSearchData[kmerIndex].kmerString;

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
  struct AwFmKmerSearchList *restrict const searchList, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){
  bool hasActiveQueries = true;
  uint8_t currentKmerLetterIndex = index->metadata.kmerLengthInSeedTable;

  while(hasActiveQueries){
    currentKmerLetterIndex++;
    hasActiveQueries = false;

    for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
      const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;
      const uint8_t kmerLength  = searchList->kmerSearchData[kmerIndex].kmerLength;
      const char    *kmerString = searchList->kmerSearchData[kmerIndex].kmerString;

      if((kmerLength >= currentKmerLetterIndex) && awFmSearchRangeIsValid(&ranges[rangesIndex])){
        hasActiveQueries = true;
        const uint8_t currentQueryLetterIndex = kmerLength - currentKmerLetterIndex;

        if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
          const uint8_t queryLetterIndex = awFmAsciiNucleotideToLetterIndex(kmerString[currentQueryLetterIndex]);
          awFmNucleotideIterativeStepBackwardSearch(index, &ranges[rangesIndex], queryLetterIndex);
        }
        else{
          const uint8_t queryLetterIndex = awFmAsciiAminoAcidToLetterIndex(kmerString[currentQueryLetterIndex]);
          awFmAminoIterativeStepBackwardSearch(index, &ranges[rangesIndex], queryLetterIndex);
        }
      }
    }
  }
}


void parallelSearchTracebackPositionLists(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, struct AwFmSearchRange *restrict const ranges,
  const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){

  for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
    const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;

    const size_t rangeLength  = awFmSearchRangeLength(&ranges[rangesIndex]);
    setPositionListCount(&searchList->kmerSearchData[kmerIndex], rangeLength);
    struct AwFmBacktrace *restrict const backtracePositionList = searchList->kmerSearchData[kmerIndex].positionBacktraceList;

    size_t positionInRangeToBacktrace = 0;
    while(positionInRangeToBacktrace < rangeLength){
      //initialize the offset.
      uint64_t offset = 0;
      uint64_t position = ranges[rangesIndex].startPtr + positionInRangeToBacktrace;
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

      backtracePositionList[positionInRangeToBacktrace].position = position;
      backtracePositionList[positionInRangeToBacktrace]._offset = offset;

      if(kmerIndex == 819){
        printf("\n\n\n420GAY420 CRITICALFAILURE: position %zu greater than bwt length %zu pos in range %zu, backtrace list ptr: %p\n",
        position, index->bwtLength, positionInRangeToBacktrace, backtracePositionList);
      }
      positionInRangeToBacktrace++;
    }
  }
}


void parallelSearchSuffixArrayLookup(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){
    for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){

    const size_t numSuffixArrayPositions = searchList->kmerSearchData[kmerIndex].count;
    for(size_t backtraceIndex = 0; backtraceIndex < numSuffixArrayPositions; backtraceIndex++){
      struct AwFmBacktrace *restrict const backtracePtr = &searchList->kmerSearchData[kmerIndex].positionBacktraceList[backtraceIndex];

        if(kmerIndex == 819){
          printf("\n\n\nin sal: position %zu greater than bwt length %zu (backtrace index %zu), kmer index %zu\n",
            backtracePtr->position, index->bwtLength, backtraceIndex, kmerIndex);
        }
      awFmSuffixArrayReadPositionParallel(index, backtracePtr);
    }
  }
}


bool setPositionListCount(struct AwFmKmerSearchData *restrict const searchData, uint32_t newCount){
    if(__builtin_expect(searchData->capacity >= newCount, 1)){
      searchData->count = newCount;
    }
    else{
      const size_t oldLengthInBytes = searchData->capacity * sizeof(uint64_t);
      const size_t newCapacity = newCount * 2;
      const size_t newLengthInBytes = newCapacity * sizeof(uint64_t);
      void *tmpPtr = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES, newLengthInBytes);
      if(__builtin_expect(tmpPtr == 0, 0)){
        return false;
      }

      memcpy(tmpPtr, searchData->positionBacktraceList, oldLengthInBytes);
      free(searchData->positionBacktraceList);
      searchData->capacity = newCapacity;
      searchData->count = newCount;
      searchData->positionBacktraceList = tmpPtr;
    }

  return true;
}


//
// //todo: deprecated
// struct AwFmParallelSearchData *awFmCreateParallelSearchData(const size_t capacity,
// const uint_fast8_t numThreads){
//
//   struct AwFmParallelSearchData *searchData = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
//                                                       sizeof(struct AwFmParallelSearchData));
//   if(searchData == NULL){
//     return NULL;
//   }
//
//   searchData->capacity    = capacity;
//   searchData->count       = 0;
//   searchData->numThreads  = numThreads;
//
//   searchData->kmerList = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES, capacity * sizeof(struct AwFmKmer));
//   if(searchData->kmerList == NULL){
//     free(searchData);
//     return NULL;
//   }
//
//   searchData->sequencePositionLists = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
//                                                 capacity * sizeof(struct AwFmBacktraceVector));
//
//   if(searchData->sequencePositionLists == NULL){
//     free(searchData->kmerList);
//     free(searchData);
//     return NULL;
//   }
//
//   //initialize all the elements in the positionLists to obviously null values
//   for(size_t i = 0; i < searchData->capacity; i++){
//     awFmBacktraceVectorCreate(&searchData->sequencePositionLists[i]);
//   }
//
//   return searchData;
// }
// void awFmDeallocParallelSearchData(struct AwFmParallelSearchData *restrict const searchData){
//   free(searchData->kmerList);
//
//   for(size_t i = 0; i < searchData->capacity; i++){
//     awFmBacktraceVectorDealloc(&searchData->sequencePositionLists[i]);
//   }
//
//   free(searchData->sequencePositionLists);
//   free(searchData);
// }
// void awFmParallelSearch(const struct AwFmIndex *restrict const index,
//   struct AwFmParallelSearchData *restrict const searchData){
//   //make local copies of the control data to encourage them to be used in a thread-shared manner.
//   const size_t searchDataCount = searchData->count;
//
//   #pragma omp parallel for num_threads(searchData->numThreads)
//   for(size_t threadBlockStartIndex = 0; threadBlockStartIndex < searchDataCount; threadBlockStartIndex += AW_FM_NUM_CONCURRENT_QUERIES){
//     const size_t threadBlockEndIndex = threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES > searchData->count?
//       searchData->count: threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES;
//     struct AwFmSearchRange ranges[AW_FM_NUM_CONCURRENT_QUERIES];
//     parallelSearchFindKmerSeedsForBlock(  index, searchData, ranges,  threadBlockStartIndex, threadBlockEndIndex);
//     parallelSearchExtendKmersInBlock(     index, searchData, ranges,  threadBlockStartIndex, threadBlockEndIndex);
//     parallelSearchTracebackPositionLists( index, searchData, ranges,  threadBlockStartIndex, threadBlockEndIndex);
//     parallelSearchSuffixArrayLookup(      index, searchData,          threadBlockStartIndex, threadBlockEndIndex);
//   }
// }
//


/*private function prototypes*/

/*private function prototypes*/


//TODO: deprecated
// void parallelSearchFindKmerSeedsForBlock(const struct AwFmIndex *restrict const index,
//   struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
//   const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){
//
//   for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
//     const uint8_t kmerLength  = searchData->kmerList[kmerIndex].length;
//     const char    *kmerString = searchData->kmerList[kmerIndex].string;
//
//     const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;
//     if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
//       ranges[rangesIndex] = awFmNucleotideKmerSeedRangeFromTable(index, kmerString, kmerLength);
//     }
//     else{
//       ranges[rangesIndex] = awFmAminoKmerSeedRangeFromTable(index, kmerString, kmerLength);
//     }
//   }
// }


// //todo: deprecated
// void parallelSearchExtendKmersInBlock(const struct AwFmIndex *restrict const index,
//   struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
//   const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){
//   bool hasActiveQueries = true;
//   uint8_t currentKmerLetterIndex = index->metadata.kmerLengthInSeedTable;
//
//   while(hasActiveQueries){
//     currentKmerLetterIndex++;
//     hasActiveQueries = false;
//
//     for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
//       const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;
//       const struct AwFmKmer *kmerPtr = &searchData->kmerList[kmerIndex];
//
//       if((kmerPtr->length >= currentKmerLetterIndex) && awFmSearchRangeIsValid(&ranges[rangesIndex])){
//         hasActiveQueries = true;
//         const uint8_t currentQueryLetterIndex = kmerPtr->length - currentKmerLetterIndex;
//
//         if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
//           const uint8_t queryLetterIndex = awFmAsciiNucleotideToLetterIndex(kmerPtr->string[currentQueryLetterIndex]);
//           awFmNucleotideIterativeStepBackwardSearch(index, &ranges[rangesIndex], queryLetterIndex);
//         }
//         else{
//           const uint8_t queryLetterIndex = awFmAsciiAminoAcidToLetterIndex(kmerPtr->string[currentQueryLetterIndex]);
//           awFmAminoIterativeStepBackwardSearch(index, &ranges[rangesIndex], queryLetterIndex);
//         }
//       }
//     }
//   }
// }


//todo: deprecated
// void parallelSearchTracebackPositionLists(const struct AwFmIndex *restrict const index,
//   struct AwFmParallelSearchData *restrict const searchData, struct AwFmSearchRange *restrict const ranges,
//   const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){
//
//   for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
//     const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;
//
//     const size_t rangeLength  = awFmSearchRangeLength(&ranges[rangesIndex]);
//     awFmBacktraceVectorSetCount(&searchData->sequencePositionLists[kmerIndex], rangeLength);
//     struct AwFmBacktrace *restrict const backtraceArray = searchData->sequencePositionLists[kmerIndex].backtraceArray;
//
//     size_t positionInRangeToBacktrace = 0;
//     while(positionInRangeToBacktrace < rangeLength){
//       uint64_t position = ranges[rangesIndex].startPtr + positionInRangeToBacktrace;
//       uint64_t offset   = 0;
//
//       if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
//         while(!awFmBwtPositionIsSampled(index, position)){
//           if(__builtin_expect(position == index->sentinelCharacterPosition, 0)){
//             position = 0;
//           }else{
//             position = awFmNucleotideBacktraceBwtPosition(index, position);
//           }
//           offset++;
//         }
//       }
//       else{
//         while(!awFmBwtPositionIsSampled(index, position)){
//           if(__builtin_expect(position == index->sentinelCharacterPosition, 0)){
//             position = 0;
//           }else{
//             position = awFmAminoBacktraceBwtPosition(index, position);
//           }
//           offset++;
//         }
//       }
//       backtraceArray[positionInRangeToBacktrace].position = position;
//       backtraceArray[positionInRangeToBacktrace].offset   = offset;
//       positionInRangeToBacktrace++;
//     }
//   }
// }

// //todo: deprecated
// void parallelSearchSuffixArrayLookup(const struct AwFmIndex *restrict const index,
//   struct AwFmParallelSearchData *restrict const searchData, const size_t threadBlockStartIndex, const size_t threadBlockEndIndex){
//     for(size_t kmerIndex = threadBlockStartIndex; kmerIndex < threadBlockEndIndex; kmerIndex++){
//
//     const struct AwFmBacktraceVector *backtraceVectorPtr = &searchData->sequencePositionLists[kmerIndex];
//     const size_t numSuffixArrayPositions = backtraceVectorPtr->count;
//     for(size_t backtraceIndex = 0; backtraceIndex < numSuffixArrayPositions; backtraceIndex++){
//       struct AwFmBacktrace *restrict const backtracePtr = &backtraceVectorPtr->backtraceArray[backtraceIndex];
//       awFmSuffixArrayReadPositionParallel(index, backtracePtr);
//     }
//   }
// }
