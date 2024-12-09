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
#include "../../src/AwFmKmerTable.h"
#include "../../src/AwFmLetter.h"
#include "../test.h"

char buffer[4096];
const char nucleotideLookup[5] = {"XACGT"};
const uint8_t aminoLookup[21] = {'a', 'c', 'd', 'e', 'f', 'g', 'h',
                                 'i', 'k', 'l', 'm', 'n', 'p', 'q',
                                 'r', 's', 't', 'v', 'w', 'y', 'z'};

void testSufSort() {
  printf("testing suf sort\n");
  // intput data
  char *text = "abracadabra";
  int n = 11;
  int i, j;

  // allocate
  int64_t *SA = malloc(n * sizeof(int64_t));
  printf("before sorting\n");
  // sort
  divsufsort64((uint8_t *)text, SA, n);
  printf("after sortingn");
  // output
  for (i = 0; i < n; ++i) {
    printf("SA[%2d] = %2zu: ", i, SA[i]);
    for (j = SA[i]; j < n; ++j) {
      printf("%c", text[j]);
    }
    printf("$\n");
  }

  // deallocate
  free(SA);
}

void testDisallowOverwrite();

void testMetadataCheck();

void testKmerTableLengths(const struct AwFmIndex *_RESTRICT_ const index,
                          const uint8_t *sequence, const size_t sequenceLength);

struct AwFmIndex *testCreateNucleotideIndex(const uint8_t *sequence,
                                            const size_t sequenceLength);

void testCreateAminoIndex();
void testPrefixSums(const struct AwFmIndex *_RESTRICT_ const index,
                    const uint8_t *sequence, const size_t sequenceLength);

int main(int argc, char **argv) {
  srand(time(NULL));

  testCreateAminoIndex();
  // testSufSort();

  // create the sequence to test
  const size_t sequenceLength = 600;
  uint8_t *sequence = malloc(sequenceLength * sizeof(char));
  if (sequence == NULL) {
    printf("critical error: sequence null!\n");
    exit(-1);
  }
  for (uint64_t i = 0; i < sequenceLength; i++) {
    sequence[i] = nucleotideLookup[rand() % 4];
  }
  // printf("sequence: ");
  // for(size_t i = 0; i < sequenceLength; i++){
  //   printf("%c, ", sequence[i]);
  // }
  printf("\n");

  struct AwFmIndex *_RESTRICT_ const index =
      testCreateNucleotideIndex(sequence, sequenceLength);
  testPrefixSums(index, sequence, sequenceLength);
  testKmerTableLengths(index, sequence, sequenceLength);

  printf("index creation tests finished, deallocating and exiting.\n");
  awFmDeallocIndex(index);

  free(sequence);
}

struct AwFmIndex *testCreateNucleotideIndex(const uint8_t *sequence,
                                            const size_t sequenceLength) {
  printf("beginning nucleotide index creation test\n");
  const char *fileSrc = "testIndex.awfmi";

  const bool keepSuffixArrayInMemory = false;
  // not using initilize because valgrind complains about the padding bytes when
  // I do it.
  struct AwFmIndexConfiguration config = {0};
  config.suffixArrayCompressionRatio = 1;
  config.kmerLengthInSeedTable = 4;
  config.alphabetType = AwFmAlphabetDna;
  config.keepSuffixArrayInMemory = false;
  config.storeOriginalSequence = false;

  const uint32_t expectedVersionNumber = AW_FM_CURRENT_VERSION_NUMBER;
  struct AwFmIndex *_RESTRICT_ index;
  enum AwFmReturnCode returnCode =
      awFmCreateIndex(&index, &config, sequence, sequenceLength, fileSrc);
  sprintf(buffer, "return code was not successful, returned %d", returnCode);
  testAssertString(returnCode >= 0, buffer);

  sprintf(buffer, "version number expected %d, got %d", index->versionNumber,
          expectedVersionNumber);
  testAssertString(index->versionNumber == expectedVersionNumber, buffer);

  sprintf(buffer, "suffix array compression ratio expected %d, got %d",
          config.suffixArrayCompressionRatio,
          index->config.suffixArrayCompressionRatio);
  testAssertString(config.suffixArrayCompressionRatio ==
                       index->config.suffixArrayCompressionRatio,
                   buffer);

  sprintf(buffer, "kmerLengthInSeedTable expected %d, got %d",
          config.kmerLengthInSeedTable, index->config.kmerLengthInSeedTable);
  testAssertString(config.kmerLengthInSeedTable ==
                       index->config.kmerLengthInSeedTable,
                   buffer);

  sprintf(buffer, "alphabetType expected %d, got %d", config.alphabetType,
          index->config.alphabetType);
  testAssertString(config.alphabetType == index->config.alphabetType, buffer);

  sprintf(buffer, "in memory suffix array flag should be (%u), but was %u",
          keepSuffixArrayInMemory, index->config.keepSuffixArrayInMemory);
  testAssertString(
      index->config.keepSuffixArrayInMemory == keepSuffixArrayInMemory, buffer);

  sprintf(buffer,
          "bwt length should be equal to sequence length + 1 (%zu), got %zu",
          sequenceLength + 1, index->bwtLength);
  testAssertString(index->bwtLength == sequenceLength + 1, buffer);

  return index;
}

void testPrefixSums(const struct AwFmIndex *_RESTRICT_ const index,
                    const uint8_t *sequence, const size_t sequenceLength) {
  const uint8_t alphabetSize =
      awFmGetAlphabetCardinality(index->config.alphabetType);
  size_t *letterCounts = malloc((alphabetSize + 2) * sizeof(size_t));
  memset(letterCounts, 0, alphabetSize * sizeof(size_t));

  printf("prefix sums: ");
  for (uint8_t i = 0; i < alphabetSize + 2; i++) {
    printf("%zu, ", index->prefixSums[i]);
  }
  printf("\n");

  for (size_t seqPos = 0; seqPos < sequenceLength; seqPos++) {
    uint8_t letterAsIndex =
        index->config.alphabetType == AwFmAlphabetAmino
            ? awFmAsciiAminoAcidToLetterIndex(sequence[seqPos])
            : awFmAsciiNucleotideToLetterIndex(sequence[seqPos]);
    letterCounts[letterAsIndex]++;
  }

  for (uint8_t i = 0; i < alphabetSize + 1; i++) {
    size_t prefixSum = 1;
    for (uint8_t j = 0; j < i; j++) {
      prefixSum += letterCounts[j];
    }
    sprintf(buffer, "prefix sum for letter index %u expected %zu, got %zu\n", i,
            prefixSum, index->prefixSums[i]);
    testAssertString(prefixSum == index->prefixSums[i], buffer);
  }

  free(letterCounts);
}

void testKmerTableLengths(const struct AwFmIndex *_RESTRICT_ const index,
                          const uint8_t *sequence,
                          const size_t sequenceLength) {
  printf("beginning kmer table lengths test\n");
  const uint8_t kmerLength = index->config.kmerLengthInSeedTable;
  for (size_t sequencePosition = 0;
       sequencePosition <= sequenceLength - kmerLength; sequencePosition++) {
    const uint8_t *kmerPtr1 = &sequence[sequencePosition];
    size_t kmerCount = 0;

    // if this kmer contains an ambiguity character, skip it.
    bool kmerContainsSentinel = false;
    for (size_t i = 0; i < kmerLength; i++) {
      kmerContainsSentinel |= kmerPtr1[i] == 'z' || kmerPtr1[i] == 'X';
    }
    if (kmerContainsSentinel) {
      continue;
    }

    for (size_t seqPos2 = 0; seqPos2 <= sequenceLength - kmerLength;
         seqPos2++) {
      const uint8_t *kmerPtr2 = &sequence[seqPos2];
      if (strncmp((char *)kmerPtr1, (char *)kmerPtr2, kmerLength) == 0) {
        kmerCount++;
      }
    }

    // printf("for kmer %.*s, found %zu.\n", kmerLength, kmerPtr1, kmerCount);
    const struct AwFmSearchRange rangeInTable =
        index->config.alphabetType == AwFmAlphabetAmino
            ? awFmAminoKmerSeedRangeFromTable(index, (char *)kmerPtr1,
                                              kmerLength)
            : awFmNucleotideKmerSeedRangeFromTable(index, (char *)kmerPtr1,
                                                   kmerLength);

    const size_t numInRange = rangeInTable.endPtr - rangeInTable.startPtr + 1;

    // printf("range start: %zu, range end: %zu\n", rangeInTable.startPtr,
    // rangeInTable.endPtr);
    sprintf(buffer,
            "for kmer %.*s, found %zu occurrences, but table stored %zu. "
            "(range %zu-%zu)\n",
            kmerLength, kmerPtr1, kmerCount, numInRange, rangeInTable.startPtr,
            rangeInTable.endPtr);
    testAssertString(kmerCount == numInRange, buffer);
  }
}

void testCreateAminoIndex(void) {
  printf("beginning amino index creation test\n");
  uint64_t sequenceLength = 40000;
  uint8_t *sequence = malloc((sequenceLength + 1) * sizeof(uint8_t));
  for (size_t i = 0; i < sequenceLength; i++) {
    sequence[i] = aminoLookup[rand() % 21];
  }
  // null terminate the sequence
  sequence[sequenceLength] = 0;
  // printf("sequence: %s\n", sequence);

  uint64_t *suffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
  divsufsort64(sequence, (int64_t *)(suffixArray + 1), sequenceLength);
  suffixArray[0] = sequenceLength;

  struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 100,
                                          .kmerLengthInSeedTable = 3,
                                          .alphabetType = AwFmAlphabetAmino,
                                          .keepSuffixArrayInMemory = false,
                                          .storeOriginalSequence = false};
  struct AwFmIndex *index;
  awFmCreateIndex(&index, &config, sequence, sequenceLength, "testIndex.awfmi");

  size_t numSentinelsEncountered = 0;
  for (size_t i = 0; i < sequenceLength + 1; i++) {
    uint8_t letterFromSequence =
        (suffixArray[i] == 0) ? '$' : sequence[suffixArray[i] - 1];

    uint8_t encodedLetterFromSequence =
        (suffixArray[i] == 0)
            ? 0
            : awFmAminoAcidLetterIndexToCompressedVector(
                  awFmAsciiAminoAcidToLetterIndex(letterFromSequence));

    size_t blockIndex = i / AW_FM_POSITIONS_PER_FM_BLOCK;
    uint8_t positionInBlock = i % AW_FM_POSITIONS_PER_FM_BLOCK;
    uint8_t byteInBlock = positionInBlock / 8;
    uint8_t bitInByte = positionInBlock % 8;

    uint8_t *letterBitVectorsAsBytes =
        (uint8_t *)index->bwtBlockList.asAmino[blockIndex].letterBitVectors;
    uint8_t letterFromBwt =
        ((letterBitVectorsAsBytes[byteInBlock] >> bitInByte) & 1) |
        (((letterBitVectorsAsBytes[32 + byteInBlock] >> bitInByte) & 1) << 1) |
        (((letterBitVectorsAsBytes[(32 * 2) + byteInBlock] >> bitInByte) & 1)
         << 2) |
        (((letterBitVectorsAsBytes[(32 * 3) + byteInBlock] >> bitInByte) & 1)
         << 3) |
        (((letterBitVectorsAsBytes[(32 * 4) + byteInBlock] >> bitInByte) & 1)
         << 4);

    if (letterFromBwt == 0x10) {
      numSentinelsEncountered++;
    }
    sprintf(buffer,
            "sequence @ i=%zu holds character encoding %#04x (letter %c), but "
            "bwt held encoding %#04x",
            i, encodedLetterFromSequence, letterFromSequence, letterFromBwt);
    testAssertString(letterFromBwt == encodedLetterFromSequence, buffer);
    if (letterFromBwt != encodedLetterFromSequence) {
      printf("num prev sentinels: %zu\n", numSentinelsEncountered);
    }
  }

  awFmDeallocIndex(index);
  free(suffixArray);
  free(sequence);
  printf("amino creation test completed\n");
}
