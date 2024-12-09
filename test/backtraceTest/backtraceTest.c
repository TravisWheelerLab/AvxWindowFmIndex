#include "../../build/FastaVector.h"
#include "../../build/divsufsort64.h"
#include "../../src/AwFmIndex.h"
#include "../../src/AwFmSearch.h"
#include "../test.h"
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

char buffer[2048];
uint8_t aminoLookup[21] = {'a', 'c', 'd', 'e', 'f', 'g', 'h',
                           'i', 'k', 'l', 'm', 'n', 'p', 'q',
                           'r', 's', 't', 'v', 'w', 'y', 'z'};
uint8_t nucleotideLookup[5] = {'a', 'g', 'c', 't', 'x'};

void testNucleotideBacktrace();
void testAminoBacktrace();

int main(int argc, char **argv) {
  srand(time(NULL));
  // srand(0);
  testNucleotideBacktrace();
  testAminoBacktrace();

  printf("backtrace function testing finished.\n");
}

void testNucleotideBacktrace() {
  for (uint8_t testNum = 0; testNum < 100; testNum++) {
    printf("nucleotide testnum: %u\n", testNum);
    const uint64_t sequenceLength = 15 + rand() % 100;
    printf("sequence length: %zu\n", sequenceLength);
    struct AwFmIndex *index;
    uint64_t *suffixArray;
    uint8_t *sequence;

    sequence = malloc((sequenceLength + 2) * sizeof(uint8_t));
    if (sequence == NULL) {
      printf("critical failure: sequence could not be allocated\n");
      exit(-1);
    }

    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = nucleotideLookup[rand() % 5];
    }
    sequence[sequenceLength] = '$';
    sequence[sequenceLength + 1] = 0;
    printf("sequence: %s\n", sequence);

    suffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
    if (suffixArray == NULL) {
      free(sequence);
      printf("critical failure: suffix array could not be allocated\n");
      exit(-2);
    }

    divsufsort64(sequence, (int64_t *)suffixArray, sequenceLength + 1);

    struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 1,
                                            .kmerLengthInSeedTable = 8,
                                            .alphabetType = AwFmAlphabetDna,
                                            .keepSuffixArrayInMemory = false,
                                            .storeOriginalSequence = false};
    awFmCreateIndex(&index, &config, sequence, sequenceLength,
                    "testIndex.awfmi");

    if (index == NULL) {
      free(sequence);
      free(suffixArray);
      printf("critical failure: index could not be created\n");
      exit(-3);
    }

    uint64_t sequencePosition = sequenceLength;
    // find the suffix array position of the sequence position
    uint64_t suffixArrayPosition = 0;

    while (suffixArray[suffixArrayPosition] != 0) {
      uint64_t newSuffixArrayPosition =
          awFmNucleotideBacktraceBwtPosition(index, suffixArrayPosition);
      uint64_t newSequencePosition = suffixArray[newSuffixArrayPosition];

      testAssertString(newSequencePosition == sequencePosition - 1, buffer);

      if (newSequencePosition != sequencePosition - 1) {
        exit(-123);
      }
      sequencePosition = newSequencePosition;
      suffixArrayPosition = newSuffixArrayPosition;
    }

    awFmDeallocIndex(index);
    free(sequence);
    free(suffixArray);
  }
}

void testAminoBacktrace(void) {
  printf("beginning Amino backtrace test.\n");
  for (uint8_t testNum = 0; testNum < 100; testNum++) {
    printf("amino testnum: %u\n", testNum);
    const uint64_t sequenceLength = 10 + rand() % 100;
    printf("amino sequence length %zu.\n", sequenceLength);
    struct AwFmIndex *index;
    uint64_t *suffixArray;
    uint8_t *sequence;

    sequence = malloc((sequenceLength + 2) * sizeof(uint8_t));
    if (sequence == NULL) {
      printf("critical failure: sequence could not be allocated\n");
      exit(-1);
    }

    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = aminoLookup[rand() % 21];
    }
    sequence[sequenceLength] = '$';
    sequence[sequenceLength + 1] = 0;
    printf("sequence: %s\n", sequence);

    suffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
    if (suffixArray == NULL) {
      free(sequence);
      printf("critical failure: suffix array could not be allocated\n");
      exit(-2);
    }

    divsufsort64(sequence, (int64_t *)(suffixArray), sequenceLength + 1);
    suffixArray[0] = sequenceLength;

    struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 1,
                                            .kmerLengthInSeedTable = 4,
                                            .alphabetType = AwFmAlphabetAmino,
                                            .keepSuffixArrayInMemory = true,
                                            .storeOriginalSequence = false};
    awFmCreateIndex(&index, &config, sequence, sequenceLength,
                    "testIndex.awfmi");

    if (index == NULL) {
      free(sequence);
      free(suffixArray);
      printf("critical failure: index could not be created\n");
      exit(-3);
    }

    // find the last character that's not a sentinel. that'll be our starting
    // backtrace position. since we can't query for sentinels, and backtrace
    // will never wrap around to the end of the sequence, backtracing through
    // these last sentinels will never happen in practice.
    uint64_t sequencePosition = sequenceLength;

    uint64_t suffixArrayPosition = 0;
    while (suffixArray[suffixArrayPosition] != 0) {

      uint64_t newSuffixArrayPosition =
          awFmAminoBacktraceBwtPosition(index, suffixArrayPosition);
      uint64_t newSequencePosition = suffixArray[newSuffixArrayPosition];

      sprintf(buffer,
              "sequence index did not decrease by 1. SA ptr %zu -> %zu, seq "
              "pos %zu -> %zu\n",
              suffixArrayPosition, newSuffixArrayPosition, sequencePosition,
              newSequencePosition);
      testAssertString(newSequencePosition == sequencePosition - 1, buffer);

      if (newSequencePosition != sequencePosition - 1) {
        exit(-123);
      }
      sequencePosition = newSequencePosition;
      suffixArrayPosition = newSuffixArrayPosition;
    }

    awFmDeallocIndex(index);
    free(sequence);
    free(suffixArray);
  }
}
