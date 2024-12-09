#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../build/divsufsort64.h"
#include "../../src/AwFmIndex.h"
#include "../test.h"

char buffer[2048];
uint8_t aminoLookup[20] = {'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l',
                           'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'};
uint8_t nucleotideLookup[4] = {'a', 'g', 'c', 't'};

void inMemorySaUncompressedTest(void);
void inMemorySaCompressedTest(void);
void inMemoryFromFileTest(void);

int main(int argc, char **argv) {
  srand(time(NULL));
  inMemorySaUncompressedTest();
  inMemorySaCompressedTest();
  inMemoryFromFileTest();
}

void inMemorySaUncompressedTest(void) {
  for (size_t testNum = 0; testNum < 500; testNum++) {
    size_t sequenceLength = 200 + rand() % 1000;
    uint8_t *sequence = malloc((sequenceLength + 100) * sizeof(uint8_t));
    printf("sequence length %zu\n", sequenceLength);
    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = nucleotideLookup[rand() % 4];
    }
    // add zeros to the end to make sure data is value, but won't compare true
    // to any strings.
    memset(sequence + sequenceLength, 0, 100);

    struct AwFmIndex *index;
    struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 1,
                                            .kmerLengthInSeedTable = 8,
                                            .alphabetType = AwFmAlphabetDna,
                                            .keepSuffixArrayInMemory = true,
                                            .storeOriginalSequence = true};

    enum AwFmReturnCode returnCode = awFmCreateIndex(
        &index, &config, sequence, sequenceLength, "testIndex.awfmi");

    if (returnCode < 0) {
      printf("error: create returned code %i!!!\n", returnCode);
    }
    size_t numKmers = 100;
    uint8_t kmerLength = 8;
    struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numKmers);
    searchList->count = numKmers;
    for (size_t i = 0; i < numKmers; i++) {
      searchList->kmerSearchData[i].kmerLength = kmerLength;
      searchList->kmerSearchData[i].kmerString = (char *)&sequence[i];
    }
    awFmParallelSearchLocate(index, searchList, 4);

    for (size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++) {
      uint64_t *positionList =
          searchList->kmerSearchData[kmerIndex].positionList;

      for (size_t sequencePosition = 0; sequencePosition < sequenceLength;
           sequencePosition++) {
        bool kmerFoundAtPosition =
            strncmp(searchList->kmerSearchData[kmerIndex].kmerString,
                    (char *)&sequence[sequencePosition], kmerLength) == 0;
        bool positionFoundInList = false;
        for (size_t backtraceIndex = 0;
             backtraceIndex < searchList->kmerSearchData[kmerIndex].count;
             backtraceIndex++) {
          if (positionList[backtraceIndex] == sequencePosition) {
            positionFoundInList = true;
          }
        }

        sprintf(buffer,
                "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in "
                "backtraceVector? %i.",
                kmerIndex, kmerLength, &sequence[sequencePosition],
                sequencePosition, kmerLength,
                searchList->kmerSearchData[kmerIndex].kmerString,
                kmerFoundAtPosition, positionFoundInList);
        if (kmerFoundAtPosition != positionFoundInList) {
          printf("not matched: position %zu,  backtrace count %u, backtrace "
                 "vector:\n",
                 sequencePosition, searchList->kmerSearchData[kmerIndex].count);
          for (size_t i = 0; i < searchList->kmerSearchData[kmerIndex].count;
               i++) {
            printf("%zu, ", positionList[i]);
          }
          printf("\n");
        }
        testAssertString(kmerFoundAtPosition == positionFoundInList, buffer);
      }
    }
    free(sequence);
    awFmDeallocKmerSearchList(searchList);
    awFmDeallocIndex(index);
  }
}

void inMemorySaCompressedTest(void) {
  for (size_t testNum = 0; testNum < 500; testNum++) {
    size_t sequenceLength = 120 + rand() % 40;
    printf("sequence length %zu\n", sequenceLength);
    uint8_t *sequence = malloc((sequenceLength + 100) * sizeof(uint8_t));
    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = nucleotideLookup[rand() % 4];
    }
    // add zeros to the end to make sure data is value, but won't compare true
    // to any strings.
    memset(sequence + sequenceLength, 0, 100);

    struct AwFmIndex *index;
    // as an extreme edge case, compress the SA to only one value.
    struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio =
                                                sequenceLength - 1,
                                            .kmerLengthInSeedTable = 8,
                                            .alphabetType = AwFmAlphabetDna,
                                            .keepSuffixArrayInMemory = true,
                                            .storeOriginalSequence = true};

    awFmCreateIndex(&index, &config, sequence, sequenceLength,
                    "testIndex.awfmi");
    size_t numKmers = 100;
    uint8_t kmerLength = 8;
    struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numKmers);
    searchList->count = numKmers;
    for (size_t i = 0; i < numKmers; i++) {
      searchList->kmerSearchData[i].kmerLength = kmerLength;
      searchList->kmerSearchData[i].kmerString = (char *)&sequence[i];
    }

    awFmParallelSearchLocate(index, searchList, 4);

    for (size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++) {
      const struct AwFmKmerSearchData *searchData =
          &searchList->kmerSearchData[kmerIndex];
      uint64_t *positionList = searchData->positionList;

      for (size_t sequencePosition = 0; sequencePosition < sequenceLength;
           sequencePosition++) {
        bool kmerFoundAtPosition =
            strncmp(searchData->kmerString, (char *)&sequence[sequencePosition],
                    kmerLength) == 0;
        bool positionFoundInList = false;
        for (size_t backtraceIndex = 0; backtraceIndex < searchData->count;
             backtraceIndex++) {
          if (positionList[backtraceIndex] == sequencePosition) {
            positionFoundInList = true;
          }
        }

        sprintf(buffer,
                "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in "
                "backtraceVector? %i.",
                kmerIndex, kmerLength, &sequence[sequencePosition],
                sequencePosition, kmerLength,
                searchList->kmerSearchData[kmerIndex].kmerString,
                kmerFoundAtPosition, positionFoundInList);

        if (kmerFoundAtPosition != positionFoundInList) {
          printf("not matched: position %zu,  backtrace count %u, backtrace "
                 "vector:\n",
                 sequencePosition, searchList->kmerSearchData[kmerIndex].count);
          for (size_t i = 0; i < searchList->kmerSearchData[kmerIndex].count;
               i++) {
            printf("%zu, ", positionList[i]);
          }
          printf("\n");
        }

        testAssertString(kmerFoundAtPosition == positionFoundInList, buffer);
      }
    }
    free(sequence);
    awFmDeallocKmerSearchList(searchList);
    awFmDeallocIndex(index);
  }
}

void inMemoryFromFileTest(void) {
  for (size_t testNum = 0; testNum < 500; testNum++) {
    size_t sequenceLength = 160;
    printf("sequence length %zu\n", sequenceLength);
    uint8_t *sequence = malloc((sequenceLength + 100) * sizeof(uint8_t));
    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = nucleotideLookup[rand() % 4];
    }
    // add zeros to the end to make sure data is value, but won't compare true
    // to any strings.
    memset(sequence + sequenceLength, 0, 100);

    struct AwFmIndex *index;
    struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 16,
                                            .kmerLengthInSeedTable = 8,
                                            .alphabetType = AwFmAlphabetDna,
                                            .keepSuffixArrayInMemory = true,
                                            .storeOriginalSequence = true};

    awFmCreateIndex(&index, &config, sequence, sequenceLength,
                    "testIndex.awfmi");
    awFmDeallocIndex(index);

    awFmReadIndexFromFile(&index, "testIndex.awfmi", true);

    size_t numKmers = 100;
    uint8_t kmerLength = 8;
    struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numKmers);
    searchList->count = numKmers;
    for (size_t i = 0; i < numKmers; i++) {
      searchList->kmerSearchData[i].kmerLength = kmerLength;
      searchList->kmerSearchData[i].kmerString = (char *)&sequence[i];
    }

    awFmParallelSearchLocate(index, searchList, 4);

    for (size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++) {
      uint64_t *positionList =
          searchList->kmerSearchData[kmerIndex].positionList;

      for (size_t sequencePosition = 0; sequencePosition < sequenceLength;
           sequencePosition++) {
        bool kmerFoundAtPosition =
            strncmp(searchList->kmerSearchData[kmerIndex].kmerString,
                    (char *)&sequence[sequencePosition], kmerLength) == 0;
        bool positionFoundInList = false;
        for (size_t backtraceIndex = 0;
             backtraceIndex < searchList->kmerSearchData[kmerIndex].count;
             backtraceIndex++) {
          if (positionList[backtraceIndex] == sequencePosition) {
            positionFoundInList = true;
          }
        }

        sprintf(buffer,
                "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in "
                "backtraceVector? %i.",
                kmerIndex, kmerLength, &sequence[sequencePosition],
                sequencePosition, kmerLength,
                searchList->kmerSearchData[kmerIndex].kmerString,
                kmerFoundAtPosition, positionFoundInList);
        if (kmerFoundAtPosition != positionFoundInList) {
          printf("not matched: position %zu,  backtrace count %u, backtrace "
                 "vector:\n",
                 sequencePosition, searchList->kmerSearchData[kmerIndex].count);
          for (size_t i = 0; i < searchList->kmerSearchData[kmerIndex].count;
               i++) {
            printf("%zu, ", positionList[i]);
          }
          printf("\n");
        }
        testAssertString(kmerFoundAtPosition == positionFoundInList, buffer);
      }
    }
    free(sequence);
    awFmDeallocKmerSearchList(searchList);
    awFmDeallocIndex(index);
  }
}
