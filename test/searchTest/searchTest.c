#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../build/divsufsort64.h"
#include "../../src/AwFmIndex.h"
#include "../../src/AwFmIndexStruct.h"
#include "../test.h"

char buffer[2048];
// randomize sequence.
// make random kmer, search for it.
// check to make sure it's in all the positions reported, and not in any not
// reported. ezpz

char buffer[2048];
uint8_t aminoLookup[21] = {'a', 'c', 'd', 'e', 'f', 'g', 'h',
                           'i', 'k', 'l', 'm', 'n', 'p', 'q',
                           'r', 's', 't', 'v', 'w', 'y', 'z'};
uint8_t nucleotideLookup[5] = {'a', 'c', 'g', 't', 'x'};

void generateRandomIndex(struct AwFmIndex **index, uint8_t **sequence,
                         size_t sequenceLength, uint64_t **suffixArray);
void testSearchForRandomKmers(const struct AwFmIndex *index, uint64_t numKmers,
                              uint8_t *sequence, size_t sequenceLength,
                              uint64_t *suffixArray);
bool sequencePositionInRange(const struct AwFmSearchRange *searchRange,
                             const uint64_t position,
                             const uint64_t *suffixArray);
void testRangeForCorrectness(const struct AwFmSearchRange *range,
                             const uint8_t *sequence,
                             const size_t sequenceLength,
                             const uint64_t *suffixArray, const char *kmer,
                             const uint8_t kmerLength);
struct AwFmSearchRange
findRangeForKmer(const struct AwFmIndex *_RESTRICT_ const index,
                 const char *kmer, const uint64_t kmerLength);

#define GENERATE_INDEX_FROM_SCRATCH

int main(int argc, char **argv) {
  srand(time(NULL));

  struct AwFmIndex *index = NULL;
  uint64_t *suffixArray = NULL;
  uint8_t *sequence = NULL;

  for (uint64_t numIndicesToTest = 0; numIndicesToTest < 25;
       numIndicesToTest++) {

    uint64_t sequenceLength = 2000 + rand() % 4000;
    printf("testing index %zu, sequence length %zu.\n", numIndicesToTest,
           sequenceLength);
    generateRandomIndex(&index, &sequence, sequenceLength, &suffixArray);
    testSearchForRandomKmers(index, 1000, sequence, sequenceLength,
                             suffixArray);

    awFmDeallocIndex(index);
  }

  printf("tests finished\n");
  free(suffixArray);
  free(sequence);
}

void generateRandomIndex(struct AwFmIndex **index, uint8_t **sequence,
                         size_t sequenceLength, uint64_t **suffixArray) {
  printf("generating random index.\n");
  // decide if we're building a nucleotide or amino index.
  enum AwFmAlphabetType alphabetType =
      (rand() & 1) ? AwFmAlphabetDna : AwFmAlphabetAmino;
  struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 200,
                                          .kmerLengthInSeedTable = 4,
                                          .alphabetType = alphabetType,
                                          .keepSuffixArrayInMemory = rand() & 1,
                                          .storeOriginalSequence = false};

  // allocate sequence and suffix array
  *sequence =
      realloc((*sequence),
              (sequenceLength + 100) *
                  sizeof(uint8_t)); //+11 is added for padding when using strcmp
  *suffixArray =
      realloc((*suffixArray), (sequenceLength + 1) * sizeof(uint64_t));

  if (sequence == NULL) {
    printf("critical failure: could not allocate sequence in unit test.");
    exit(-1);
  } else if (suffixArray == NULL) {
    printf("critical failure: could not allocate suffix array in unit test.");
    exit(-2);
  }

  uint8_t *characterLookupTable =
      alphabetType == AwFmAlphabetDna ? nucleotideLookup : aminoLookup;
  uint8_t alphabetCardinalty = awFmGetAlphabetCardinality(alphabetType);
  for (size_t i = 0; i < sequenceLength; i++) {
    uint8_t characterIndex =
        rand() % (alphabetCardinalty + 1); //+1 is for ambiguity character
    (*sequence)[i] = characterLookupTable[characterIndex];
  }
  (*sequence)[sequenceLength] = '$';

  int64_t divSufSortReturnCode =
      divsufsort64((*sequence), (int64_t *)(*suffixArray), sequenceLength + 1);
  if (divSufSortReturnCode < 0) {
    printf("critical failure: divsufsort returned error code %li\n",
           divSufSortReturnCode);
    exit(-3);
  }
  enum AwFmReturnCode awFmReturnCode = awFmCreateIndex(
      index, &config, *sequence, sequenceLength, "testIndex.awfmi");
  if (awFmReturnCode < 0) {
    printf("critical failure: create index returned error code %i\n",
           awFmReturnCode);
    exit(-4);
  }
}

void testSearchForRandomKmers(const struct AwFmIndex *index, uint64_t numKmers,
                              uint8_t *sequence, size_t sequenceLength,
                              uint64_t *suffixArray) {
  for (size_t kmerNum = 0; kmerNum < numKmers; kmerNum++) {
    char kmer[13];
    uint8_t kmerLength = 1 + rand() % 10;

    // build the kmer.
    // printf("nuclookup ptr %p, amino %p, index p: %p\n", nucleotideLookup,
    // aminoLookup, index);
    uint8_t *characterLookupTable =
        (index->config.alphabetType == AwFmAlphabetDna) ? &nucleotideLookup[0]
                                                        : &aminoLookup[0];
    uint8_t alphabetCardinalty =
        awFmGetAlphabetCardinality(index->config.alphabetType);
    for (size_t i = 0; i < kmerLength; i++) {
      kmer[i] = characterLookupTable[rand() % alphabetCardinalty];
    }

    struct AwFmSearchRange range = findRangeForKmer(index, kmer, kmerLength);
    testRangeForCorrectness(&range, sequence, sequenceLength, suffixArray, kmer,
                            kmerLength);
  }
}

bool sequencePositionInRange(const struct AwFmSearchRange *searchRange,
                             const uint64_t position,
                             const uint64_t *suffixArray) {

  for (size_t i = searchRange->startPtr; i <= searchRange->endPtr; i++) {
    if (suffixArray[i] == position) {
      return true;
    }
  }
  return false;
}

void testRangeForCorrectness(const struct AwFmSearchRange *range,
                             const uint8_t *sequence,
                             const size_t sequenceLength,
                             const uint64_t *suffixArray, const char *kmer,
                             const uint8_t kmerLength) {
  for (size_t position = 0; position < sequenceLength; position++) {
    const bool kmerFoundAtPosition =
        strncmp((char *)(sequence + position), kmer, kmerLength) == 0;

    bool kmerExpectedAtPosition =
        sequencePositionInRange(range, position, suffixArray);

    if (kmerFoundAtPosition && !kmerExpectedAtPosition) {
      sprintf(buffer,
              "kmer  \"%.*s\" was found at position %zu, (%.*s) but was not "
              "expected in range [%zu, %zu]:",
              kmerLength, kmer, position, kmerLength, sequence + position,
              range->startPtr, range->endPtr);

      testAssertString(false, buffer);

      printf("[");
      for (size_t i = range->startPtr; i <= range->endPtr; i++) {
        printf("%zu, ", suffixArray[i]);
      }
      printf("]\n");
    } else if (!kmerFoundAtPosition && kmerExpectedAtPosition) {
      sprintf(buffer,
              "kmer  \"%.*s\" returned range [%zu, %zu], which contains "
              "position %zu, but that position was not in range:",
              kmerLength, kmer, range->startPtr, range->endPtr, position);

      printf("[");
      for (size_t i = range->startPtr; i <= range->endPtr; i++) {
        printf("%zu, ", suffixArray[i]);
      }
      printf("]\n");
    }
  }
}

struct AwFmSearchRange
findRangeForKmer(const struct AwFmIndex *_RESTRICT_ const index,
                 const char *kmer, const uint64_t kmerLength) {
  return awFmFindSearchRangeForString(index, kmer, kmerLength);
}
