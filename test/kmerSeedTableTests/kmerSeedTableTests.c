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
#include "../../src/AwFmSearch.h"
#include "../test.h"

char buffer[2048];
uint8_t aminoLookup[21] = {'a', 'c', 'd', 'e', 'f', 'g', 'h',
                           'i', 'k', 'l', 'm', 'n', 'p', 'q',
                           'r', 's', 't', 'v', 'w', 'y', 'z'};
uint8_t nucleotideLookup[5] = {'a', 'g', 'c', 't', 'x'};

void checkRangeForCorrectness(
    const struct AwFmSearchRange *_RESTRICT_ const range,
    const uint64_t *_RESTRICT_ const suffixArray, const char *kmer,
    const uint8_t kmerLength, const uint8_t *_RESTRICT_ const sequence,
    const uint64_t sequenceLength);
bool rangeCompare(struct AwFmSearchRange range1, struct AwFmSearchRange range2);
void testAllKmerRanges(struct AwFmIndexConfiguration *config,
                       uint64_t sequenceLength);

int main(int argc, char **argv) {
  srand(time(NULL));
  printf("main\n");

  struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 240,
                                          .kmerLengthInSeedTable = 5,
                                          .alphabetType = AwFmAlphabetDna,
                                          .keepSuffixArrayInMemory = false};
  testAllKmerRanges(&config, 2000);

  printf("testing amino ranges\n");
  config.alphabetType = AwFmAlphabetAmino;
  config.kmerLengthInSeedTable = 2;
  testAllKmerRanges(&config, 2000);

  printf("end\n");
}

void checkRangeForCorrectness(
    const struct AwFmSearchRange *_RESTRICT_ const range,
    const uint64_t *_RESTRICT_ const suffixArray, const char *kmer,
    const uint8_t kmerLength, const uint8_t *_RESTRICT_ const sequence,
    const uint64_t sequenceLength) {

  // make a null terminated copy of the kmer
  char kmerBuffer[kmerLength + 1];
  memcpy(kmerBuffer, kmer, kmerLength);
  kmerBuffer[kmerLength] = 0;

  for (uint64_t position = 0; position <= sequenceLength - kmerLength;
       position++) {
    bool kmerFoundAtPosition =
        strncmp(kmer, (char *)sequence + position, kmerLength) == 0;

    // printf("kmer %.*s at pos %zu, matches kmer %.*s\n", kmerLength,
    // sequence+position, position, kmerLength, kmer);
    bool positionFoundInRange = false;

    for (uint64_t suffixArrayPosition = range->startPtr;
         suffixArrayPosition <= range->endPtr; suffixArrayPosition++) {
      // printf("SAP %zu, POS %zu\n", suffixArray[suffixArrayPosition],
      // position);
      positionFoundInRange |= (suffixArray[suffixArrayPosition] == position);
    }

    if (!positionFoundInRange) {
      char sequenceBuffer[kmerLength + 1];
      memcpy(sequenceBuffer, sequence + position, kmerLength);
      sequenceBuffer[kmerLength] = 0;

      sprintf(buffer,
              "range had position %zu for kmer %s, but that position holds "
              "kmer %s.\n",
              position, kmerBuffer, sequenceBuffer);
    } else {
      sprintf(buffer,
              "position %zu matches kmer %s, but was not represented in the "
              "range.\n",
              position, kmerBuffer);
    }
    testAssertString((kmerFoundAtPosition == positionFoundInRange), buffer);
    if (kmerFoundAtPosition != positionFoundInRange) {
      printf("failure: position: %zu. kmer found at position %u, found in "
             "range %u\n",
             position, kmerFoundAtPosition, positionFoundInRange);
      printf("kmer @pos: %.*s\n", 2, &sequence[position]);
      printf("kmer: %s, \n", kmerBuffer);
      printf("seq len: %zu, range [%zu, %zu]\n", sequenceLength,
             range->startPtr, range->endPtr);
      // printf("range values: ");
      printf("sa positions: ");
      for (uint64_t i = range->startPtr; i <= range->endPtr; i++) {
        printf("%zu, ", suffixArray[i]);
      }
      printf("\n");
      for (uint64_t ptr = range->startPtr; ptr <= range->endPtr; ptr++) {
        printf(" %zu, ", suffixArray[ptr]);
      }
      printf("\n");

      exit(-3);
    }
  }
}

bool rangeCompare(struct AwFmSearchRange range1,
                  struct AwFmSearchRange range2) {
  return (range1.startPtr == range2.startPtr) &&
         (range1.endPtr == range2.endPtr);
}

void testAllKmerRanges(struct AwFmIndexConfiguration *config,
                       uint64_t sequenceLength) {
  const uint8_t kmerLength = config->kmerLengthInSeedTable;
  struct AwFmIndex *index;
  for (uint8_t testNum = 0; testNum < 4; testNum++) {
    uint8_t *sequence = malloc((sequenceLength + 1) * sizeof(uint8_t));
    printf("generating index\n");
    if (sequence == NULL) {
      printf("CRITICAL FAILURE: could not allocate sequence\n");
      exit(-1);
    }
    for (uint64_t i = 0; i < sequenceLength; i++) {
      sequence[i] = (config->alphabetType == AwFmAlphabetAmino)
                        ? aminoLookup[rand() % 21]
                        : nucleotideLookup[rand() % 5];
    }
    // null terminate the sequence
    sequence[sequenceLength] = 0;

    printf("sequence: %s\n", sequence);
    // create the reference suffix array
    uint64_t *suffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
    if (suffixArray == NULL) {
      printf("CRITICAL FAILURE: could not allocate suffix array\n");
      exit(-1);
    }

    // set the suffix array. remember, the first position isn't included by
    // libdivsufsort, and needs to be manually set to the sentinel character
    // position.
    suffixArray[0] = sequenceLength;
    divsufsort64(sequence, (int64_t *)(suffixArray + 1), sequenceLength);

    printf("suffix array generated.\n");
    // printf("suffix array: ");
    // for(uint64_t i = 0; i < sequenceLength+1; i++){
    //   printf("%zu,",  suffixArray[i]);
    //   if(((i+1) % 10) == 0){
    //     printf("  ");
    //   }
    // }
    // printf("\n");
    //
    // printf("bwt: ");
    // for(uint64_t i = 0; i < sequenceLength+1; i++){
    //   printf("%c, ", (suffixArray[i] == 0)? '$' :
    //   sequence[suffixArray[i]-1]); if(((i+1) % 10) == 0){
    //     printf("  ");
    //   }
    // }
    // printf("\n");

    awFmCreateIndex(&index, config, sequence, sequenceLength,
                    "testIndex.awfmi");

    struct AwFmSearchRange range = {0, 0};
    char *kmer2 = "agaaa";
    if (config->alphabetType != AwFmAlphabetAmino) {
      awFmNucleotideNonSeededSearch(index, kmer2, 5, &range);
    } else {
      awFmAminoNonSeededSearch(index, kmer2, 5, &range);
    }
    printf("range for agaaa: %zu, %zu\n", range.startPtr, range.endPtr);
    printf("index generated.\n");
    // printf("prefix sums: ");
    // for(uint8_t i = 0; i <20; i++){
    //   printf("[%d]: %zu', ", i, index->prefixSums[i]);
    // }
    // printf("\n");

    // printf("kmer table:\n");
    // printf("kmer table len: %zu\n", awFmGetKmerTableLength(index));
    // for(uint64_t i = 0; i < awFmGetKmerTableLength(index); i++){
    // for(uint64_t i = 5600; i < 5620; i++){
    //   printf("%zu: [%zu, %zu]\n", i, index->kmerSeedTable.table[i].startPtr,
    //   index->kmerSeedTable.table[i].endPtr);
    // }

    printf("created index\n");
    clock_t startTime = clock();

    size_t tableLength = awFmGetKmerTableLength(index);
    char kmer[kmerLength];
    printf("checking ranges...\n");
    for (size_t kmerIndex = 0; kmerIndex < tableLength; kmerIndex++) {
      if (kmerIndex % 100 == 0) {
        printf("checking kmer index %zu\n", kmerIndex);
      }
      size_t kmerIndexCopy = kmerIndex;
      for (uint8_t letterInKmer = 0; letterInKmer < kmerLength;
           letterInKmer++) {
        if (index->config.alphabetType != AwFmAlphabetAmino) {
          kmer[letterInKmer] = nucleotideLookup[(kmerIndexCopy % 4)];
          kmerIndexCopy /= 4;
        } else {
          kmer[letterInKmer] = aminoLookup[(kmerIndexCopy % 20)];
          kmerIndexCopy /= 20;
        }
      }
      struct AwFmSearchRange kmerRange =
          config->alphabetType == AwFmAlphabetAmino
              ? awFmAminoKmerSeedRangeFromTable(index, kmer, kmerLength)
              : awFmNucleotideKmerSeedRangeFromTable(index, kmer, kmerLength);

      checkRangeForCorrectness(&kmerRange, suffixArray, kmer, kmerLength,
                               sequence, sequenceLength);
    }

    clock_t endTime = clock();
    clock_t elapsedTime = (endTime - startTime) / (CLOCKS_PER_SEC / 1000);
    printf("elapsed time for 8,000 kmers: %zu\n", elapsedTime);
    printf("test  %u complete\n", testNum);

    awFmDeallocIndex(index);
    free(sequence);
    free(suffixArray);
  }
}
