#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "AwFmIndex.h"
#include "AwFmIndexStruct.h"
#include "AwFmKmerTable.h"
#include "AwFmLetter.h"
#include "AwFmSearch.h"
#include "AwFmSuffixArray.h"

#define NUM_CONCURRENT_QUERIES 32
#define DEFAULT_POSITION_LIST_CAPACITY 4

void parallelSearchFindKmerSeedsForBlock(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmKmerSearchList *_RESTRICT_ const searchList,
    struct AwFmSearchRange *_RESTRICT_ const ranges,
    const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

void parallelSearchExtendKmersInBlock(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmKmerSearchList *_RESTRICT_ const searchList,
    struct AwFmSearchRange *_RESTRICT_ const ranges,
    const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

enum AwFmReturnCode parallelSearchTracebackPositionLists(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmKmerSearchList *_RESTRICT_ const searchList,
    struct AwFmSearchRange *_RESTRICT_ const ranges,
    const size_t threadBlockStartIndex, const size_t threadBlockEndIndex);

bool setPositionListCount(
    struct AwFmKmerSearchData *_RESTRICT_ const searchData, uint32_t count);

struct AwFmKmerSearchList *awFmCreateKmerSearchList(const size_t capacity) {
  // struct AwFmKmerSearchList *searchList =
  // aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
  //   sizeof(struct AwFmKmerSearchList));
  struct AwFmKmerSearchList *searchList =
      malloc(sizeof(struct AwFmKmerSearchList));
  if (searchList == NULL) {
    return NULL;
  }

  searchList->capacity = capacity;
  searchList->count = 0;

  // searchList->kmerSearchData = aligned_alloc(AW_FM_CACHE_LINE_SIZE_IN_BYTES,
  //   capacity * sizeof(struct AwFmKmerSearchData));
  searchList->kmerSearchData =
      malloc(capacity * sizeof(struct AwFmKmerSearchData));

  if (searchList->kmerSearchData == NULL) {
    free(searchList);
    return NULL;
  }

  bool positionListAllocationFailed = false;
  for (size_t i = 0; i < capacity; i++) {
    searchList->kmerSearchData[i].kmerString = NULL;
    searchList->kmerSearchData[i].kmerLength = 0;
    searchList->kmerSearchData[i].capacity = DEFAULT_POSITION_LIST_CAPACITY;
    searchList->kmerSearchData[i].count = 0;
    searchList->kmerSearchData[i].positionList =
        malloc(DEFAULT_POSITION_LIST_CAPACITY * sizeof(uint64_t));

    // check for an allocation failure
    positionListAllocationFailed |=
        (searchList->kmerSearchData[i].positionList == NULL);
  }

  // if any of the allocations failed, dealloc everything and return NULL.
  if (positionListAllocationFailed) {
    for (size_t i = 0; i < capacity; i++) {
      free(searchList->kmerSearchData[i].positionList);
    }
    free(searchList->kmerSearchData);
    free(searchList);
    return NULL;
  }

  return searchList;
}

void awFmDeallocKmerSearchList(
    struct AwFmKmerSearchList *_RESTRICT_ const searchList) {
  for (size_t i = 0; i < searchList->capacity; i++) {
    free(searchList->kmerSearchData[i].positionList);
  }
  free(searchList->kmerSearchData);
  free(searchList);
}

enum AwFmReturnCode
awFmParallelSearchLocate(const struct AwFmIndex *_RESTRICT_ const index,
                         struct AwFmKmerSearchList *_RESTRICT_ const searchList,
                         uint32_t numThreads) {

  const uint32_t searchListCount = searchList->count;
  enum AwFmReturnCode atomicReturnCode = AwFmSuccess;
  if (numThreads > 1) {
#pragma omp parallel for num_threads(numThreads)
    for (size_t threadBlockStartIndex = 0;
         threadBlockStartIndex < searchListCount;
         threadBlockStartIndex += AW_FM_NUM_CONCURRENT_QUERIES) {

      const size_t threadBlockEndIndex =
          threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES >
                  searchList->count
              ? searchList->count
              : threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES;

      struct AwFmSearchRange ranges[AW_FM_NUM_CONCURRENT_QUERIES];

      parallelSearchFindKmerSeedsForBlock(index, searchList, ranges,
                                          threadBlockStartIndex,
                                          threadBlockEndIndex);
      parallelSearchExtendKmersInBlock(index, searchList, ranges,
                                       threadBlockStartIndex,
                                       threadBlockEndIndex);
      enum AwFmReturnCode rc = parallelSearchTracebackPositionLists(
          index, searchList, ranges, threadBlockStartIndex,
          threadBlockEndIndex);
      if (__builtin_expect(awFmReturnCodeIsFailure(rc), 0)) {
#pragma omp atomic write
        atomicReturnCode = AwFmFileReadFail;
      }
    }
    return atomicReturnCode;
  } else {
    for (size_t threadBlockStartIndex = 0;
         threadBlockStartIndex < searchListCount;
         threadBlockStartIndex += AW_FM_NUM_CONCURRENT_QUERIES) {
      const size_t threadBlockEndIndex =
          threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES >
                  searchList->count
              ? searchList->count
              : threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES;
      struct AwFmSearchRange ranges[AW_FM_NUM_CONCURRENT_QUERIES];

      parallelSearchFindKmerSeedsForBlock(index, searchList, ranges,
                                          threadBlockStartIndex,
                                          threadBlockEndIndex);
      parallelSearchExtendKmersInBlock(index, searchList, ranges,
                                       threadBlockStartIndex,
                                       threadBlockEndIndex);
      enum AwFmReturnCode rc = parallelSearchTracebackPositionLists(
          index, searchList, ranges, threadBlockStartIndex,
          threadBlockEndIndex);
      if (__builtin_expect(awFmReturnCodeIsFailure(rc), 0)) {
        return rc;
      }
    }
    return AwFmSuccess;
  }
}

void awFmParallelSearchCount(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmKmerSearchList *_RESTRICT_ const searchList,
    uint32_t numThreads) {

  const uint32_t searchListCount = searchList->count;

  if (numThreads > 1) {
#pragma omp parallel for num_threads(numThreads)
    for (size_t threadBlockStartIndex = 0;
         threadBlockStartIndex < searchListCount;
         threadBlockStartIndex += AW_FM_NUM_CONCURRENT_QUERIES) {

      const size_t threadBlockEndIndex =
          threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES >
                  searchList->count
              ? searchList->count
              : threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES;
      struct AwFmSearchRange ranges[AW_FM_NUM_CONCURRENT_QUERIES];

      parallelSearchFindKmerSeedsForBlock(index, searchList, ranges,
                                          threadBlockStartIndex,
                                          threadBlockEndIndex);
      parallelSearchExtendKmersInBlock(index, searchList, ranges,
                                       threadBlockStartIndex,
                                       threadBlockEndIndex);

      // load the range lengths into the count member variables.
      for (size_t i = threadBlockStartIndex; i < threadBlockEndIndex; i++) {
        searchList->kmerSearchData[i].count =
            awFmSearchRangeLength(&ranges[i - threadBlockStartIndex]);
      }
    }
  } else {
    // exact duplicate of above code, without the omp pragma, so it doesn't kill
    // performance with only 1 thread.
    for (size_t threadBlockStartIndex = 0;
         threadBlockStartIndex < searchListCount;
         threadBlockStartIndex += AW_FM_NUM_CONCURRENT_QUERIES) {

      const size_t threadBlockEndIndex =
          threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES >
                  searchList->count
              ? searchList->count
              : threadBlockStartIndex + AW_FM_NUM_CONCURRENT_QUERIES;
      struct AwFmSearchRange ranges[AW_FM_NUM_CONCURRENT_QUERIES];

      parallelSearchFindKmerSeedsForBlock(index, searchList, ranges,
                                          threadBlockStartIndex,
                                          threadBlockEndIndex);
      parallelSearchExtendKmersInBlock(index, searchList, ranges,
                                       threadBlockStartIndex,
                                       threadBlockEndIndex);

      // load the range lengths into the count member variables.
      for (size_t i = threadBlockStartIndex; i < threadBlockEndIndex; i++) {
        searchList->kmerSearchData[i].count =
            awFmSearchRangeLength(&ranges[i - threadBlockStartIndex]);
      }
    }
  }
}

void parallelSearchFindKmerSeedsForBlock(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmKmerSearchList *_RESTRICT_ const searchList,
    struct AwFmSearchRange *_RESTRICT_ const ranges,
    const size_t threadBlockStartIndex, const size_t threadBlockEndIndex) {

  for (size_t kmerIndex = threadBlockStartIndex;
       kmerIndex < threadBlockEndIndex; kmerIndex++) {
    const struct AwFmKmerSearchData *searchData =
        &searchList->kmerSearchData[kmerIndex];
    const uint64_t kmerLength = searchData->kmerLength;
    const char *kmerString = searchData->kmerString;

    const bool queryCanUseKmerTable =
        awFmQueryCanUseKmerTable(index, kmerString, kmerLength);

    const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;

    // these are used if the kmer is ineligible for using the kmerSeedTable
    const uint64_t kmerStringNonSeededStart =
        kmerLength < index->config.kmerLengthInSeedTable
            ? 0
            : kmerLength - index->config.kmerLengthInSeedTable;
    const uint64_t kmerStringNonSeededLength =
        kmerLength < index->config.kmerLengthInSeedTable
            ? kmerLength
            : index->config.kmerLengthInSeedTable;

    if (index->config.alphabetType != AwFmAlphabetAmino) {
      // TODO: reimplement partial seeded search when it's implementable
      if (queryCanUseKmerTable) {
        ranges[rangesIndex] =
            awFmNucleotideKmerSeedRangeFromTable(index, kmerString, kmerLength);
      } else {
        awFmNucleotideNonSeededSearch(
            index, kmerString + kmerStringNonSeededStart,
            kmerStringNonSeededLength, &ranges[rangesIndex]);
      }
    } else {
      if (queryCanUseKmerTable) {
        ranges[rangesIndex] =
            awFmAminoKmerSeedRangeFromTable(index, kmerString, kmerLength);
      } else {
        awFmAminoNonSeededSearch(index, kmerString + kmerStringNonSeededStart,
                                 kmerStringNonSeededLength,
                                 &ranges[rangesIndex]);
      }
    }
  }
}

void parallelSearchExtendKmersInBlock(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmKmerSearchList *_RESTRICT_ const searchList,
    struct AwFmSearchRange *_RESTRICT_ const ranges,
    const size_t threadBlockStartIndex, const size_t threadBlockEndIndex) {
  bool hasActiveQueries = true;
  uint64_t currentKmerLetterIndex = index->config.kmerLengthInSeedTable;

  while (hasActiveQueries) {
    currentKmerLetterIndex++;
    hasActiveQueries = false;

    for (size_t kmerIndex = threadBlockStartIndex;
         kmerIndex < threadBlockEndIndex; kmerIndex++) {
      const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;
      const struct AwFmKmerSearchData *_RESTRICT_ const searchData =
          &searchList->kmerSearchData[kmerIndex];
      const uint64_t kmerLength = searchData->kmerLength;
      const char *kmerString = searchData->kmerString;

      if ((kmerLength >= currentKmerLetterIndex) &&
          awFmSearchRangeIsValid(&ranges[rangesIndex])) {
        hasActiveQueries = true;
        const uint64_t currentQueryLetterIndex =
            kmerLength - currentKmerLetterIndex;

        if (index->config.alphabetType != AwFmAlphabetAmino) {
          const uint8_t queryLetterIndex = awFmAsciiNucleotideToLetterIndex(
              kmerString[currentQueryLetterIndex]);
          awFmNucleotideIterativeStepBackwardSearch(index, &ranges[rangesIndex],
                                                    queryLetterIndex);
        } else {
          const uint8_t queryLetterIndex = awFmAsciiAminoAcidToLetterIndex(
              kmerString[currentQueryLetterIndex]);
          awFmAminoIterativeStepBackwardSearch(index, &ranges[rangesIndex],
                                               queryLetterIndex);
        }
      }
    }
  }
}

enum AwFmReturnCode parallelSearchTracebackPositionLists(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmKmerSearchList *_RESTRICT_ const searchList,
    struct AwFmSearchRange *_RESTRICT_ const ranges,
    const size_t threadBlockStartIndex, const size_t threadBlockEndIndex) {

  for (size_t kmerIndex = threadBlockStartIndex;
       kmerIndex < threadBlockEndIndex; kmerIndex++) {
    const uint64_t rangesIndex = kmerIndex - threadBlockStartIndex;

    struct AwFmKmerSearchData *searchData =
        &searchList->kmerSearchData[kmerIndex];
    const size_t rangeLength = awFmSearchRangeLength(&ranges[rangesIndex]);
    setPositionListCount(searchData, rangeLength);
    // struct AwFmBacktrace *_RESTRICT_ const backtracePositionList = {
    //   .position = searchList->kmerSearchData[kmerIndex].positionList;
    // };searchList->kmerSearchData[kmerIndex].positionBacktraceList;

    for (size_t indexOfPositionToBacktrace = 0;
         indexOfPositionToBacktrace < rangeLength;
         indexOfPositionToBacktrace++) {

      // initialize the offset.
      struct AwFmBacktrace backtrace = {
          .position = ranges[rangesIndex].startPtr + indexOfPositionToBacktrace,
          .offset = 0};

      if (index->config.alphabetType != AwFmAlphabetAmino) {
        while (!awFmBwtPositionIsSampled(index, backtrace.position)) {
          backtrace.position =
              awFmNucleotideBacktraceBwtPosition(index, backtrace.position);
          backtrace.offset++;
        }
      } else {
        while (!awFmBwtPositionIsSampled(index, backtrace.position)) {
          backtrace.position =
              awFmAminoBacktraceBwtPosition(index, backtrace.position);
          backtrace.offset++;
        }
      }

      if (__builtin_expect(
              awFmSuffixArrayReadPositionParallel(index, &backtrace), 0) !=
          AwFmSuccess) {
        return AwFmFileReadFail;
      }
      searchData->positionList[indexOfPositionToBacktrace] = backtrace.position;
    }
  }
  return AwFmSuccess;
}

bool setPositionListCount(
    struct AwFmKmerSearchData *_RESTRICT_ const searchData, uint32_t newCount) {
  if (__builtin_expect(searchData->capacity >= newCount, 1)) {
    searchData->count = newCount;
  } else {
    const size_t newCapacity = newCount;
    const size_t newLengthInBytes = newCapacity * sizeof(uint64_t);
    void *tmpPtr = realloc(searchData->positionList, newLengthInBytes);
    if (__builtin_expect(tmpPtr == 0, 0)) {
      fprintf(stderr, "Critical memory failure: could not allocate memory for "
                      "position list.\n");
      return false;
    }

    searchData->capacity = newCapacity;
    searchData->count = newCount;
    searchData->positionList = tmpPtr;
  }

  return true;
}
