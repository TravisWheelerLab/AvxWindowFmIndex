#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../build/divsufsort64.h"
#include "../../src/AwFmCreate.h"
#include "../../src/AwFmIndex.h"
#include "../../src/AwFmIndexStruct.h"
#include "../../src/AwFmKmerTable.h"
#include "../../src/AwFmLetter.h"
#include "../../src/AwFmOccurrence.h"
#include "../../src/AwFmParallelSearch.h"
#include "../../src/AwFmSearch.h"
#include "../../src/AwFmSuffixArray.h"
#include "../test.h"

char buffer[2048];
uint8_t aminoLookup[21] = {'a', 'c', 'd', 'e', 'f', 'g', 'h',
                           'i', 'k', 'l', 'm', 'n', 'p', 'q',
                           'r', 's', 't', 'v', 'w', 'y', 'z'};
uint8_t nucleotideLookup[5] = {'a', 'g', 'c', 't', 'x'};

void testParallelSearchNucleotide();
void testParallelSearchAmino();
void testParallelCount();

const uint8_t saCompressionRatio = 8;

int main(int argc, char **argv) {
  size_t seedTime = time(NULL);
  // printf("***seed time %zu\n", seedTime);
  // size_t seedTime = 1628294143;
  srand(seedTime);
  testParallelCount();
  testParallelSearchNucleotide();
  testParallelSearchAmino();

  printf("parallel search testing finished.\n");
}

void testParallelSearchAmino(void) {
  struct AwFmIndex *index;
  struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio =
                                              saCompressionRatio,
                                          .kmerLengthInSeedTable = 5,
                                          .alphabetType = AwFmAlphabetAmino,
                                          .keepSuffixArrayInMemory = true,
                                          .storeOriginalSequence = false};

  const uint64_t sequenceLength = 10000 + rand() % 50000;
  printf("creating amino sequence of length %zu.\n", sequenceLength);
  uint8_t *sequence = malloc((sequenceLength + 11) * sizeof(uint8_t));
  if (sequence == NULL) {
    printf("critical error: could not allocate sequence\n");
    exit(-1);
  }

  for (uint64_t i = 0; i < sequenceLength; i++) {
    sequence[i] = aminoLookup[rand() % 21];
  }
  // null terminate the sequence for easy printing
  sequence[sequenceLength] = 0;

  // create a reference Suffix Array
  uint64_t *referenceSuffixArray =
      malloc((sequenceLength + 1) * sizeof(uint64_t));
  int64_t divSufSortReturnCode = divsufsort64(
      sequence, (int64_t *)(referenceSuffixArray), sequenceLength + 1);
  if (divSufSortReturnCode < 0) {
    free(referenceSuffixArray);
    sprintf(buffer,
            "reference suffix array could not be generated, error code %zi.\n",
            divSufSortReturnCode);
    testAssertString(false, buffer);
  }

  awFmCreateIndex(&index, &config, sequence, sequenceLength, "testIndex.awfmi");
  printf("index generated\n");

  for (size_t saIndex = 0; saIndex < sequenceLength + 1;
       saIndex += saCompressionRatio) {
    uint64_t refValue = referenceSuffixArray[saIndex];
    uint64_t position = saIndex;
    enum AwFmReturnCode rc =
        awFmReadPositionsFromSuffixArray(index, &position, 1);
    testAssertString(rc == AwFmSuccess || rc == AwFmFileReadOkay,
                     "reading value from index SA failed.");
    sprintf(buffer,
            "index compressed SA value @ index %zu (%zu) did not match "
            "reference val %zu.",
            saIndex, position, refValue);
    testAssertString(refValue == position, buffer);

    sprintf(buffer,
            "suffix array value %zu was greater than sequence length %zu\n.",
            position, sequenceLength);
    testAssertString(position <= sequenceLength, buffer);
  }

  const size_t kmerSearchListCapacity = 1000 + (rand() % 6000);
  for (uint8_t numThreads = 1; numThreads < 16; numThreads++) {
    struct AwFmKmerSearchList *searchList =
        awFmCreateKmerSearchList(kmerSearchListCapacity);

    if (searchList == NULL) {
      printf("critical error: parallel search data could not be allocated\n");
      exit(-2);
    }

    // fill the search data struct
    searchList->count = kmerSearchListCapacity - (rand() % 200);
    printf("beginning test with  %zu kmers, with %d threads.\n",
           searchList->count, numThreads);
    for (size_t i = 0; i < searchList->count; i++) {
      uint16_t kmerLength = 4 + rand() % 20;
      searchList->kmerSearchData[i].kmerLength = kmerLength;

      char *const kmer = malloc(kmerLength * sizeof(char));
      searchList->kmerSearchData[i].kmerString = kmer;

      for (uint16_t kmerLetterIndex = 0; kmerLetterIndex < kmerLength;
           kmerLetterIndex++) {
        kmer[kmerLetterIndex] = aminoLookup[rand() % 20];
      }
    }

    // query for the sequence positions
    clock_t startTime = clock();
    awFmParallelSearchLocate(index, searchList, numThreads);
    clock_t endTime = clock();
    clock_t totalTime = ((endTime - startTime));
    printf("total time: %zu ticks\n", totalTime);

    printf("par search data count %zu\n", searchList->count);
    for (size_t kmerIndex = 0; kmerIndex < searchList->count; kmerIndex++) {
      const struct AwFmKmerSearchData *searchData =
          &searchList->kmerSearchData[kmerIndex];
      for (size_t hit = 0; hit < searchData->count; hit++) {
        if (searchData->positionList[hit] > sequenceLength) {
          sprintf(buffer,
                  "position %zu was greater than the sequenceLength %zu!\n",
                  searchData->positionList[hit], sequenceLength);
          testAssertString(searchData->positionList[hit] < sequenceLength,
                           buffer);
        }
      }

      sprintf(buffer,
              "backtrace vector's capacity was lower than it's count!\n");
      testAssertString(searchData->capacity >= searchData->count, buffer);

      // search the sequence, and at each position, ensure:
      // if the kmer is at this location, make sure it's in the kmer list.
      // if it's not, make sure it isn't in the list
      for (size_t sequencePosition = 0; sequencePosition < sequenceLength;
           sequencePosition++) {
        bool kmerFoundAtThisPosition =
            (strncmp((char *)&sequence[sequencePosition],
                     searchData->kmerString, searchData->kmerLength)) == 0;

        const uint64_t *positionList = searchData->positionList;
        bool thisPositionInList = false;
        for (size_t positionIndexInList = 0;
             positionIndexInList < searchData->count; positionIndexInList++) {
          const size_t realPosition = positionList[positionIndexInList];
          thisPositionInList |= (realPosition == sequencePosition);
        }

        if (kmerFoundAtThisPosition && !thisPositionInList) {

          sprintf(buffer,
                  "kmer %.*s (%.*s) was found at position %zu, but was not "
                  "represented in the position list.",
                  (uint32_t)searchData->kmerLength, searchData->kmerString,
                  (uint32_t)searchData->kmerLength, sequence + sequencePosition,
                  sequencePosition);
          testAssertString(false, buffer);
          exit(-1);
        } else if (!kmerFoundAtThisPosition && thisPositionInList) {
          printf("kmer %zu has count %u\n", kmerIndex, searchData->count);
          printf("position list: ");
          for (size_t i = 0; i < searchData->count; i++) {
            printf("pos: %zu, ", positionList[i]);
          }
          printf("\n");
          sprintf(buffer,
                  "kmer %.*s was in the position list, but was not found at "
                  "the expected sequence position %zu (%.*s found "
                  "there instead).",
                  (uint32_t)searchData->kmerLength, searchData->kmerString,
                  sequencePosition, (uint32_t)searchData->kmerLength,
                  &sequence[sequencePosition]);
          testAssertString(false, buffer);
          exit(-4);
          printf("sequence length: %zu\n", sequenceLength);
        }
      }
    }

    // deallocate the kmer strings
    for (size_t i = 0; i < searchList->count; i++) {
      free(searchList->kmerSearchData[i].kmerString);
    }
    // check the search data to make sure it was accurate.
    awFmDeallocKmerSearchList(searchList);
  }
  free(referenceSuffixArray);
  free(sequence);
  awFmDeallocIndex(index);
}

void testParallelSearchNucleotide() {
  struct AwFmIndex *index;
  struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio =
                                              saCompressionRatio,
                                          .kmerLengthInSeedTable = 9,
                                          .alphabetType = AwFmAlphabetDna,
                                          .keepSuffixArrayInMemory = true,
                                          .storeOriginalSequence = false};

  const uint64_t sequenceLength = 5000 + rand() % 5000;
  printf("creating nucleotide sequence of length %zu.\n", sequenceLength);
  uint8_t *sequence = malloc((sequenceLength + 11) * sizeof(uint8_t));
  if (sequence == NULL) {
    printf("critical error: could not allocate sequence\n");
    exit(-1);
  }

  for (uint64_t i = 0; i < sequenceLength; i++) {
    sequence[i] = nucleotideLookup[rand() % 5];
  }
  // null terminate the sequence for easy printing
  sequence[sequenceLength] = 0;

  awFmCreateIndex(&index, &config, sequence, sequenceLength, "testIndex.awfmi");

  const size_t kmerSearchListCapacity = 100 + (rand() % 10000);
  for (uint8_t numThreads = 1; numThreads < 16; numThreads++) {

    struct AwFmKmerSearchList *searchList =
        awFmCreateKmerSearchList(kmerSearchListCapacity);

    if (searchList == NULL) {
      printf("critical error: parallel search data could not be allocated\n");
      exit(-2);
    }

    // fill the search data struct
    searchList->count = kmerSearchListCapacity - (rand() % 200);
    printf("beginning test with  %zu kmers, with %d threads.\n",
           searchList->count, numThreads);
    for (size_t i = 0; i < searchList->count; i++) {
      uint16_t kmerLength = 7 + rand() % 30;
      searchList->kmerSearchData[i].kmerLength = kmerLength;

      char *const kmer = malloc(kmerLength * sizeof(char));
      searchList->kmerSearchData[i].kmerString = kmer;

      for (uint16_t kmerLetterIndex = 0; kmerLetterIndex < kmerLength;
           kmerLetterIndex++) {
        kmer[kmerLetterIndex] = nucleotideLookup[rand() % 4];
      }
    }

    // query for the sequence positions
    clock_t startTime = clock();
    awFmParallelSearchLocate(index, searchList, numThreads);
    clock_t endTime = clock();
    clock_t totalTime = ((endTime - startTime));
    printf("total time: %zu ticks\n", totalTime);
    // printf("parallel search completed, checking answers...\n");

    for (size_t kmerIndex = 0; kmerIndex < searchList->count; kmerIndex++) {
      const struct AwFmKmerSearchData *searchData =
          &searchList->kmerSearchData[kmerIndex];

      sprintf(buffer,
              "position vector's capacity was lower than it's count!\n");
      testAssertString(searchData->capacity >= searchData->count, buffer);

      // search the sequence, and at each position, ensure:
      // if the kmer is at this location, make sure it's in the kmer list.
      // if it's not, make sure it isn't in the list
      for (size_t sequencePosition = 0; sequencePosition < sequenceLength;
           sequencePosition++) {
        bool kmerFoundAtThisPosition =
            (strncmp((char *)&sequence[sequencePosition],
                     searchData->kmerString, searchData->kmerLength)) == 0;

        const uint64_t *positionList = searchData->positionList;
        bool thisPositionInList = false;
        for (size_t positionIndexInList = 0;
             positionIndexInList < searchData->count; positionIndexInList++) {
          const size_t realPosition = positionList[positionIndexInList];
          thisPositionInList |= (realPosition == sequencePosition);
        }

        if (kmerFoundAtThisPosition && !thisPositionInList) {
          printf("kmer %zu has count %u\n", kmerIndex, searchData->count);
          printf("position list: ");
          for (size_t i = 0; i < searchData->count; i++) {
            printf("pos: %zu, ", positionList[i]);
          }
          printf("from seq: %.*s\n", (uint32_t)searchData->kmerLength,
                 sequence + sequencePosition);
          sprintf(buffer,
                  "kmer %.*s (%.*s) was found at position %zu, but was not "
                  "represented in the position list.",
                  (uint32_t)searchData->kmerLength, searchData->kmerString,
                  (uint32_t)searchData->kmerLength, sequence + sequencePosition,
                  sequencePosition);
          printf("kmer @ position %zu in range: %.*s. in sequence: %.*s\n",
                 sequencePosition, (uint32_t)searchData->kmerLength,
                 searchData->kmerString, (uint32_t)searchData->kmerLength,
                 sequence + sequencePosition);
          printf("(seq pos = %zu), pos list count %u\n", sequencePosition,
                 searchData->count);
          for (size_t i = 0; i < searchData->count; i++) {
            printf("pos: %zu, ", positionList[i]);
          }
          printf("\n");

          testAssertString(false, buffer);
          exit(-2);
        } else if (!kmerFoundAtThisPosition && thisPositionInList) {
          printf("kmer %zu has count %u\n", kmerIndex, searchData->count);
          printf("position list: ");
          for (size_t i = 0; i < searchData->count; i++) {
            printf("pos: %zu, ", positionList[i]);
          }
          printf("\n");
          sprintf(buffer,
                  "kmer %.*s was in the position list, but was not found at "
                  "the expected sequence position %zu (%.*s found "
                  "there instead).",
                  (uint32_t)searchData->kmerLength, searchData->kmerString,
                  sequencePosition, (uint32_t)searchData->kmerLength,
                  &sequence[sequencePosition]);
          testAssertString(false, buffer);
          exit(-3);
          printf("sequence length: %zu\n", sequenceLength);
        }
      }
    }

    // deallocate the kmer strings
    for (size_t i = 0; i < searchList->count; i++) {
      free(searchList->kmerSearchData[i].kmerString);
    }
    // check the search data to make sure it was accurate.
    awFmDeallocKmerSearchList(searchList);
  }
  free(sequence);
  awFmDeallocIndex(index);
}

void testParallelCount(void) {
  struct AwFmIndex *index;
  struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio =
                                              saCompressionRatio,
                                          .kmerLengthInSeedTable = 5,
                                          .alphabetType = AwFmAlphabetAmino,
                                          .keepSuffixArrayInMemory = true,
                                          .storeOriginalSequence = false};

  const uint64_t sequenceLength = 10000 + rand() % 50000;
  printf("creating amino sequence of length %zu.\n", sequenceLength);
  uint8_t *sequence = malloc((sequenceLength + 11) * sizeof(uint8_t));
  if (sequence == NULL) {
    printf("critical error: could not allocate sequence\n");
    exit(-1);
  }

  for (uint64_t i = 0; i < sequenceLength; i++) {
    sequence[i] = aminoLookup[rand() % 4];
  }
  // null terminate the sequence for easy printing
  sequence[sequenceLength] = 0;

  awFmCreateIndex(&index, &config, sequence, sequenceLength, "testIndex.awfmi");
  printf("index generated\n");

  const size_t kmerSearchListCapacity = 1000 + (rand() % 6000);
  for (uint8_t numThreads = 1; numThreads < 16; numThreads++) {

    struct AwFmKmerSearchList *searchList =
        awFmCreateKmerSearchList(kmerSearchListCapacity);

    if (searchList == NULL) {
      printf("critical error: parallel search data could not be allocated\n");
      exit(-2);
    }

    // fill the search data struct
    searchList->count = kmerSearchListCapacity - (rand() % 200);
    printf("beginning test with  %zu kmers, with %d threads.\n",
           searchList->count, numThreads);
    for (size_t i = 0; i < searchList->count; i++) {
      uint16_t kmerLength = 4 + rand() % 20;
      searchList->kmerSearchData[i].kmerLength = kmerLength;

      char *const kmer = malloc(kmerLength * sizeof(char));
      searchList->kmerSearchData[i].kmerString = kmer;

      for (uint16_t kmerLetterIndex = 0; kmerLetterIndex < kmerLength;
           kmerLetterIndex++) {
        kmer[kmerLetterIndex] = aminoLookup[rand() % 20];
      }
    }

    // query for the sequence positions
    clock_t startTime = clock();
    awFmParallelSearchCount(index, searchList, numThreads);
    clock_t endTime = clock();
    clock_t totalTime = ((endTime - startTime));
    printf("total time: %zu ticks\n", totalTime);
    // printf("parallel search completed, checking answers...\n");

    for (size_t kmerIndex = 0; kmerIndex < searchList->count; kmerIndex++) {

      const struct AwFmKmerSearchData *searchData =
          &searchList->kmerSearchData[kmerIndex];

      // search the sequence, and at each position, ensure:
      // if the kmer is at this location, make sure it's in the kmer list.
      // if it's not, make sure it isn't in the list
      size_t officialKmerCount = 0;
      for (size_t sequencePosition = 0;
           sequencePosition < sequenceLength - searchData->kmerLength;
           sequencePosition++) {
        bool kmerFoundAtThisPosition =
            (strncmp((char *)&sequence[sequencePosition],
                     searchData->kmerString, searchData->kmerLength)) == 0;
        officialKmerCount += kmerFoundAtThisPosition ? 1 : 0;
      }
      sprintf(buffer,
              "for kmer index %zu, returned count %u did not match actual "
              "count %zu.",
              kmerIndex, searchData->count, officialKmerCount);
      testAssertString(officialKmerCount == searchData->count, buffer);
    }

    // deallocate the kmer strings
    for (size_t i = 0; i < searchList->count; i++) {
      free(searchList->kmerSearchData[i].kmerString);
    }
    // check the search data to make sure it was accurate.
    awFmDeallocKmerSearchList(searchList);
  }
  free(sequence);
  awFmDeallocIndex(index);
}
