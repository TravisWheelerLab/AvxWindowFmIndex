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
#include "../../src/AwFmSuffixArray.h"
#include "../test.h"

char buffer[2048];
uint8_t aminoLookup[21] = {'a', 'c', 'd', 'e', 'f', 'g', 'h',
                           'i', 'k', 'l', 'm', 'n', 'p', 'q',
                           'r', 's', 't', 'v', 'w', 'y', '$'};
uint8_t nucleotideLookup[5] = {'a', 'g', 'c', 't', '$'};

void suffixArrayTest(void);
void sequenceRecallTest(void);
void indexReadTest(void);

int main(int argc, char **argv) {
  srand(time(NULL));
  sequenceRecallTest();
  suffixArrayTest();
  indexReadTest();
}

void sequenceRecallTest(void) {
  printf("beginning sequence recall test.");
  for (size_t testNum = 0; testNum < 100; testNum++) {
    printf("    test num: %zu\n", testNum);
    size_t sequenceLength = 4000 + rand() % 400;
    uint8_t *sequence = malloc(sequenceLength * sizeof(uint8_t));
    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = nucleotideLookup[rand() % 4];
    }

    struct AwFmIndex *index;
    struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 240,
                                            .kmerLengthInSeedTable = 2,
                                            .alphabetType = AwFmAlphabetDna,
                                            .keepSuffixArrayInMemory = false,
                                            .storeOriginalSequence = true};
    awFmCreateIndex(&index, &config, sequence, sequenceLength,
                    "testIndex.awfmi");

    char sequenceBuffer[2048];
    for (size_t sequencePosition = 0; sequencePosition < sequenceLength;
         sequencePosition++) {
      size_t sequenceSegmentLength = rand() % 10 + 1;
      if (sequencePosition + sequenceSegmentLength >= index->bwtLength) {
        sequenceSegmentLength = (index->bwtLength - 1) - sequencePosition;
      }
      // size_t priorFlankLength = rand()%10+1;
      // size_t postFlankLength = rand()%10+1;
      enum AwFmReturnCode returnCode = awFmReadSequenceFromFile(
          index, sequencePosition, sequenceSegmentLength, sequenceBuffer);

      sprintf(buffer,
              "awFmReadSequenceFromFile returned failure code %i for position "
              "%zu (segment length %zu, total len %zu)",
              returnCode, sequencePosition, sequenceSegmentLength,
              sequenceLength);
      testAssertString(returnCode > 0, buffer);

      int compareResult =
          strncmp(sequenceBuffer, (char *)(sequence + sequencePosition),
                  sequenceSegmentLength);

      sprintf(buffer,
              "sequence segment %.*s  length %zu at position %zu did not match "
              "what was returned in the buffer %.*s",
              (int)sequenceSegmentLength, (char *)(sequence + sequencePosition),
              sequenceSegmentLength, sequencePosition,
              (int)sequenceSegmentLength, sequenceBuffer);

      testAssertString(compareResult == 0, buffer);
    }
    awFmDeallocIndex(index);

    // test amino recall
    config.alphabetType = AwFmAlphabetAmino;
    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = aminoLookup[rand() % 20];
    }
    awFmCreateIndex(&index, &config, sequence, sequenceLength,
                    "testIndex.awfmi");
    for (size_t sequencePosition = 0; sequencePosition < sequenceLength;
         sequencePosition++) {
      size_t sequenceSegmentLength = rand() % 10 + 1;
      if (sequencePosition + sequenceSegmentLength >= index->bwtLength) {
        sequenceSegmentLength = (index->bwtLength - 1) - sequencePosition;
      }
      enum AwFmReturnCode returnCode = awFmReadSequenceFromFile(
          index, sequencePosition, sequenceSegmentLength, sequenceBuffer);

      sprintf(buffer, "awFmReadSequenceFromFile returned failure code %i",
              returnCode);
      testAssertString(returnCode > 0, buffer);

      int compareResult =
          strncmp(sequenceBuffer, (char *)(sequence + sequencePosition),
                  sequenceSegmentLength);
      if (compareResult != 0) {
        printf("seq seg len: %zu, res %i, seq: %.*s, buf %.*s, first l %u, %u, "
               "return code %u\n",
               sequenceSegmentLength, compareResult, (int)sequenceSegmentLength,
               (char *)(sequence + sequencePosition),
               (int)sequenceSegmentLength, sequenceBuffer,
               sequence[sequencePosition], sequenceBuffer[0], returnCode);
      }
      sprintf(buffer,
              "sequence segment %.*s did not match what was returned in the "
              "buffer %.*s",
              (int)sequenceSegmentLength, (char *)(sequence + sequencePosition),
              (int)sequenceSegmentLength, sequenceBuffer);

      testAssertString(compareResult == 0, buffer);
    }
    awFmDeallocIndex(index);
    free(sequence);
  }
}

void suffixArrayTest(void) {
  printf("beginning suffix array test\n");

  for (size_t compressionRatio = 1; compressionRatio < 32; compressionRatio++) {
    printf("    compression ratio %zu, nucleotide\n", compressionRatio);
    struct AwFmIndex *index;
    struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio =
                                                compressionRatio,
                                            .kmerLengthInSeedTable = 2,
                                            .alphabetType = AwFmAlphabetDna,
                                            .keepSuffixArrayInMemory = false,
                                            .storeOriginalSequence = true};

    const size_t sequenceLength = 4000 + rand() % 2000;
    uint8_t *sequence = malloc((sequenceLength + 1) * sizeof(uint8_t));
    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = nucleotideLookup[rand() % 5];
    }
    sequence[sequenceLength] = 0;

    int64_t *referenceSuffixArray =
        malloc((sequenceLength + 1) * sizeof(uint64_t));
    divsufsort64(sequence, (int64_t *)(referenceSuffixArray + 1),
                 sequenceLength);
    referenceSuffixArray[0] = sequenceLength;

    awFmCreateIndex(&index, &config, sequence, sequenceLength,
                    "testIndex.awfmi");
    size_t numElementsInSuffixArray = sequenceLength / compressionRatio;

    for (size_t i = 0; i < numElementsInSuffixArray; i++) {
      size_t bwtPosition = i * compressionRatio;
      uint64_t positionFromSuffixArray = referenceSuffixArray[bwtPosition];

      struct AwFmBacktrace backtrace = {.position = bwtPosition, .offset = 0};
      enum AwFmReturnCode returnCode =
          awFmSuffixArrayReadPositionParallel(index, &backtrace);
      sprintf(buffer,
              "on position %zu, suffix array read position parallel returned "
              "error code %i.",
              i, returnCode);
      testAssertString(returnCode > 0, buffer);

      sprintf(buffer,
              "on position %zu, backtrace position %zu did not match value in "
              "suffix array (%zu).",
              i, backtrace.position, positionFromSuffixArray);

      testAssertString(backtrace.position == positionFromSuffixArray, buffer);
    }

    free(referenceSuffixArray);
    awFmDeallocIndex(index);
    printf("    compression ratio %zu, amino\n", compressionRatio);
    // rereandomize the sequence with amino acids.
    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = aminoLookup[rand() % 21];
    }

    referenceSuffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
    divsufsort64(sequence, (referenceSuffixArray + 1), sequenceLength);
    referenceSuffixArray[0] = sequenceLength;

    config.alphabetType = AwFmAlphabetAmino;

    awFmCreateIndex(&index, &config, sequence, sequenceLength,
                    "testIndex.awfmi");

    numElementsInSuffixArray = sequenceLength / compressionRatio;

    for (size_t i = 0; i < numElementsInSuffixArray; i++) {
      size_t bwtPosition = i * compressionRatio;
      uint64_t positionFromSuffixArray = referenceSuffixArray[bwtPosition];

      struct AwFmBacktrace backtrace = {.position = bwtPosition, .offset = 0};
      enum AwFmReturnCode returnCode =
          awFmSuffixArrayReadPositionParallel(index, &backtrace);
      sprintf(buffer,
              "on position %zu, suffix array read position parallel returned "
              "error code %i.",
              i, returnCode);
      testAssertString(returnCode > 0, buffer);

      sprintf(buffer,
              "on position %zu, backtrace position %zu did not match value in "
              "suffix array (%zu).",
              i, backtrace.position, positionFromSuffixArray);

      testAssertString(backtrace.position == positionFromSuffixArray, buffer);
    }

    free(sequence);
    free(referenceSuffixArray);
    awFmDeallocIndex(index);
  }
}

// tests to make sure that the index read with awFmReadIndexFromFile creates
// an index identical one to the one it's built from.
void indexReadTest(void) {

  for (size_t testNum = 0; testNum < 100; testNum++) {
    printf("test %zu\n", testNum);
    const uint8_t kmerLengthInSeedTable = rand() % 3 + 2;
    const enum AwFmAlphabetType alphabetType =
        rand() & 1 ? AwFmAlphabetDna : AwFmAlphabetAmino;

    const uint64_t sequenceLength = 10000 + (rand() % 1000);
    uint8_t *sequence = malloc(sequenceLength * sizeof(uint8_t));
    if (sequence == NULL) {
      printf("ERROR: sequence was null on index read test, num %zu\n", testNum);
      exit(-444);
    }

    // randomize the sequence
    for (size_t i = 0; i < sequenceLength; i++) {
      sequence[i] = alphabetType == AwFmAlphabetDna
                        ? nucleotideLookup[rand() % 5]
                        : aminoLookup[rand() % 21];
    }

    struct AwFmIndex *index;
    struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 255,
                                            .kmerLengthInSeedTable =
                                                kmerLengthInSeedTable,
                                            .alphabetType = alphabetType,
                                            .keepSuffixArrayInMemory = false,
                                            .storeOriginalSequence = false};

    enum AwFmReturnCode returnCode = awFmCreateIndex(
        &index, &config, sequence, sequenceLength, "testindex.awfmi");
    if (returnCode < 0) {
      printf("ERROR: creating initial index returned error code %i\n",
             returnCode);
    }
    free(sequence);
    fclose(index->fileHandle);

    // load another copy of the index
    struct AwFmIndex *indexFromFile;
    returnCode =
        awFmReadIndexFromFile(&indexFromFile, "testindex.awfmi", false);
    if (returnCode < 0) {
      printf("ERROR: reading index copy returned error code %i\n", returnCode);
    }

    // compare the bwt lengths
    sprintf(buffer,
            "the bwt lengths of the original (%zu) and the one from file (%zu) "
            "did not match.",
            index->bwtLength, indexFromFile->bwtLength);
    testAssertString(index->bwtLength == indexFromFile->bwtLength, buffer);

    // version number
    sprintf(buffer,
            "the version numbers of the original (%i) and the one from file "
            "(%i) did not match.",
            index->versionNumber, indexFromFile->versionNumber);
    testAssertString(index->versionNumber == indexFromFile->versionNumber,
                     buffer);

    // suffix array compression ratio
    sprintf(buffer,
            "the suffix array compression ratios of the original (%i) and the "
            "one from file (%i) did not match.",
            index->config.suffixArrayCompressionRatio,
            indexFromFile->config.suffixArrayCompressionRatio);
    testAssertString(index->config.suffixArrayCompressionRatio ==
                         indexFromFile->config.suffixArrayCompressionRatio,
                     buffer);

    // kmer length in seed table
    sprintf(buffer,
            "the kmer length in seed table values of the original (%i) and the "
            "one from file (%i) did not match.",
            index->config.kmerLengthInSeedTable,
            indexFromFile->config.kmerLengthInSeedTable);
    testAssertString(index->config.kmerLengthInSeedTable ==
                         indexFromFile->config.kmerLengthInSeedTable,
                     buffer);

    // alphabet type
    sprintf(buffer,
            "the alphabet type of the original (%i) and the one from file (%i) "
            "did not match.",
            index->config.alphabetType, indexFromFile->config.alphabetType);
    testAssertString(index->config.alphabetType ==
                         indexFromFile->config.alphabetType,
                     buffer);

    // suffix array file offset
    sprintf(buffer,
            "the suffixArrayFileOffset of the original (%zu) and the one from "
            "file (%zu) did not match.",
            index->suffixArrayFileOffset, indexFromFile->suffixArrayFileOffset);
    testAssertString(index->suffixArrayFileOffset ==
                         indexFromFile->suffixArrayFileOffset,
                     buffer);

    // sequenceFileOffset
    sprintf(buffer,
            "the sequenceFileOffset of the original (%zu) and the one from "
            "file (%zu) did not match.",
            index->sequenceFileOffset, indexFromFile->sequenceFileOffset);
    testAssertString(
        index->sequenceFileOffset == indexFromFile->sequenceFileOffset, buffer);

    // bwtBlockList
    size_t blockListLengthInBytes =
        alphabetType == AwFmAlphabetAmino
            ? awFmNumBlocksFromBwtLength(index->bwtLength) *
                  sizeof(struct AwFmAminoBlock)
            : awFmNumBlocksFromBwtLength(index->bwtLength) *
                  sizeof(struct AwFmNucleotideBlock);

    sprintf(buffer, "the bwt block lists of the original and the one from file "
                    "did not match.");
    testAssertString(memcmp(index->bwtBlockList.asNucleotide,
                            indexFromFile->bwtBlockList.asNucleotide,
                            blockListLengthInBytes) == 0,
                     buffer);

    // prefix sums
    sprintf(
        buffer,
        "the prefix sums of the original and the one from file did not match.");
    testAssertString(
        memcmp(index->prefixSums, indexFromFile->prefixSums,
               awFmGetPrefixSumsLength(alphabetType) * sizeof(uint64_t)) == 0,
        buffer);

    // kmerSeedTable
    size_t kmerTableLengthInBytes =
        awFmGetKmerTableLength(index) * sizeof(struct AwFmSearchRange);

    sprintf(
        buffer,
        "the kmer table of the original and the one from file did not match.");
    testAssertString(memcmp(index->kmerSeedTable, indexFromFile->kmerSeedTable,
                            kmerTableLengthInBytes) == 0,
                     buffer);

    // dealloc both indices
    // here, we're cheating by doing what dealloc does, except not trying to
    // close the already closed file
    free(index->bwtBlockList.asNucleotide);
    free(index->prefixSums);
    free(index->kmerSeedTable);
    free(index->suffixArray.values);
    free(index);
    awFmDeallocIndex(indexFromFile);
  }
}
