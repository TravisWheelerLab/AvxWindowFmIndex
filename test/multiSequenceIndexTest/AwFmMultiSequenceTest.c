#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../lib/FastaVector/src/FastaVector.h"
#include "../../src/AwFmIndex.h"
#include "../../src/AwFmIndexStruct.h"
#include "../../src/AwFmSearch.h"
#include "../../src/AwFmSuffixArray.h"
#include "../test.h"
// #include "../../lib/libdivsufsort/include/divsufsort64.h"

char buffer[2048];
char aminoLookup[21] = {'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm',
                        'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y', 'x'};
char nucleotideLookup[5] = {'a', 'g', 'c', 't', 'x'};

struct FastaVector *
generateMultiSequenceFastaVector(const size_t numSequences,
                                 const size_t maxSequenceLength,
                                 const bool isAmino);

void randomizeSequenceBuffer(char *sequenceBuffer, size_t sequenceLength,
                             bool isAmino);
void getRawSequenceFromFastaVector(const struct FastaVector *const fastaVector,
                                   char **compositeSequenceBufferPtr,
                                   size_t *compositeSequenceLength);

void testFromFasta();

struct AwFmIndexConfiguration generateReasonableRandomMetadata();

void testAwFmIndexIdenticalForFastaVector(void);
void testAwFmIndexFastaVectorDataMatchesExpected(void);
void testAwFmIndexGivesCorrectLocalPositions(void);
void testAwFmIndexGivesCorrectHeaders(void);
void compareIndicesForEqualityIgnoreVersion(const struct AwFmIndex *index1,
                                            const struct AwFmIndex *index2);
void checkAllGlobalPositionsForCorrectLocalPositions(
    const struct FastaVector *_RESTRICT_ const fastaVector,
    const struct AwFmIndex *_RESTRICT_ const fastaVectorIndex);

int main(int argc, char **argv) {
  srand(time(NULL));
  testAwFmIndexIdenticalForFastaVector();
  testAwFmIndexFastaVectorDataMatchesExpected();
  testAwFmIndexGivesCorrectLocalPositions();
  testFromFasta();
}

struct FastaVector *
generateMultiSequenceFastaVector(const size_t numSequences,
                                 const size_t maxSequenceLength,
                                 const bool isAmino) {

  char headerBuffer[128] = {0};
  struct FastaVector *fastaVector = malloc(sizeof(struct FastaVector));
  testAssertString(fastaVector != NULL,
                   "in generateMultiSequenceFastaVector, fasta vector malloc "
                   "failed (returned NULL).");

  enum FastaVectorReturnCode rc = fastaVectorInit(fastaVector);
  testAssertString(rc == FASTA_VECTOR_OK,
                   "in generateMultiSequenceFastaVector, fastaVector init did "
                   "not return OK.");

  for (size_t sequenceIndex = 0; sequenceIndex < numSequences;
       sequenceIndex++) {
    size_t thisSequenceLength = (rand() % maxSequenceLength) + 1;
    char *sequenceBuffer = malloc(thisSequenceLength * sizeof(char));
    randomizeSequenceBuffer(sequenceBuffer, thisSequenceLength, isAmino);
    sprintf(headerBuffer, "header%zu| config", sequenceIndex);
    rc = fastaVectorAddSequenceToList(fastaVector, headerBuffer,
                                      strlen(headerBuffer), sequenceBuffer,
                                      thisSequenceLength);
    testAssertString(rc == FASTA_VECTOR_OK,
                     "in generateMultiSequenceFastaVector, fastaVector adding "
                     "sequence to list did not return OK.");
    free(sequenceBuffer);
  }

  return fastaVector;
}

void randomizeSequenceBuffer(char *sequenceBuffer, size_t sequenceLength,
                             bool isAmino) {
  const char *const letterLookupTable =
      isAmino ? aminoLookup : nucleotideLookup;
  const size_t letterLookupTableLength = isAmino ? 21 : 5;
  for (size_t letterIndex = 0; letterIndex < sequenceLength; letterIndex++) {
    const uint8_t chosenLetter = rand() % letterLookupTableLength;
    sequenceBuffer[letterIndex] = letterLookupTable[chosenLetter];
  }
}

void getRawSequenceFromFastaVector(const struct FastaVector *const fastaVector,
                                   char **compositeSequenceBufferPtr,
                                   size_t *compositeSequenceLength) {

  *compositeSequenceLength = fastaVector->sequence.count;
  *compositeSequenceBufferPtr = fastaVector->sequence.charData;
}

struct AwFmIndexConfiguration generateReasonableRandomMetadata() {
  struct AwFmIndexConfiguration config;
  config.suffixArrayCompressionRatio = (rand() % 20) + 1;
  config.alphabetType = (rand() % 2) == 0 ? AwFmAlphabetAmino : AwFmAlphabetDna;
  config.kmerLengthInSeedTable = config.alphabetType == AwFmAlphabetDna
                                     ? (rand() % 10) + 2
                                     : (rand() % 4) + 1;
  config.keepSuffixArrayInMemory = true;
  config.storeOriginalSequence = true;

  return config;
}

void testAwFmIndexIdenticalForFastaVector() {
  for (size_t testNum = 0; testNum < 10; testNum++) {
    struct FastaVector *fastaVector = NULL;
    // generating the config also sets the test's alphabet type.
    struct AwFmIndexConfiguration config = generateReasonableRandomMetadata();
    printf("test identical # %zu, config generated: sacr %u, klist %u, "
           "alphabet %u, ksim %u\n",
           testNum, config.suffixArrayCompressionRatio,
           config.kmerLengthInSeedTable, config.alphabetType,
           config.keepSuffixArrayInMemory);

    const bool isAmino = config.alphabetType == AwFmAlphabetAmino;
    const size_t numSequences = (rand() % 10) + 1;
    const size_t maxSequenceLength = 10000;
    const char *fastaFileSrc = "sequences.fasta";
    const uint32_t fileLineLength = (rand() % 522) + 1;
    fastaVector = generateMultiSequenceFastaVector(numSequences,
                                                   maxSequenceLength, isAmino);

    enum FastaVectorReturnCode rc =
        fastaVectorWriteFasta(fastaFileSrc, fastaVector, fileLineLength);
    testAssertString(
        rc == FASTA_VECTOR_OK,
        "writing fasta vector to file did not return FASTA_VECTOR_OK.");

    char *sequenceCollectionPtr = NULL;
    size_t sequenceCollectionLength = 0;
    getRawSequenceFromFastaVector(fastaVector, &sequenceCollectionPtr,
                                  &sequenceCollectionLength);
    testAssertString(sequenceCollectionPtr != NULL,
                     "sequence collection ptr was null!");

    char *defaultIndexFileSrc = "defaultIndex.awfmi";
    char *fastaVectorIndexFileSrc = "fastaVectorIndex.awfmi";
    struct AwFmIndex *defaultIndex = NULL;
    struct AwFmIndex *fastaVectorIndex = NULL;

    enum AwFmReturnCode awFmRc = awFmCreateIndex(
        &defaultIndex, &config, (uint8_t *)sequenceCollectionPtr,
        sequenceCollectionLength, defaultIndexFileSrc);
    testAssertString(awFmRc == AwFmFileWriteOkay,
                     "creating default index did not return AwFmSuccess.");
    awFmRc = awFmCreateIndexFromFasta(&fastaVectorIndex, &config, fastaFileSrc,
                                      fastaVectorIndexFileSrc);
    testAssertString(awFmRc == AwFmFileWriteOkay,
                     "creating fastaVector index did not return AwFmSuccess");

    compareIndicesForEqualityIgnoreVersion(defaultIndex, fastaVectorIndex);
    // TODO! why is fastVectorIndex having zeros in most of it's data?

    fastaVectorDealloc(fastaVector);
    free(fastaVector);
    awFmDeallocIndex(fastaVectorIndex);
    awFmDeallocIndex(defaultIndex);
  }
}

void testAwFmIndexFastaVectorDataMatchesExpected(void) {
  for (size_t testNum = 0; testNum < 10; testNum++) {
    struct FastaVector *fastaVector = NULL;
    // generating the config also sets the test's alphabet type.
    struct AwFmIndexConfiguration config = generateReasonableRandomMetadata();
    const bool isAmino = config.alphabetType == AwFmAlphabetAmino;
    const size_t numSequences = (rand() % 10) + 1;
    const size_t maxSequenceLength = 10000;
    const char *fastaFileSrc = "sequences.fasta";
    const uint32_t fileLineLength = (rand() % 522) + 1;
    fastaVector = generateMultiSequenceFastaVector(numSequences,
                                                   maxSequenceLength, isAmino);
    printf("test fastaVector data # %zu, config generated: sacr %u, klist %u, "
           "alphabet %u, ksim %u\n",
           testNum, config.suffixArrayCompressionRatio,
           config.kmerLengthInSeedTable, config.alphabetType,
           config.keepSuffixArrayInMemory);

    enum FastaVectorReturnCode rc =
        fastaVectorWriteFasta(fastaFileSrc, fastaVector, fileLineLength);
    testAssertString(
        rc == FASTA_VECTOR_OK,
        "writing fasta vector to file did not return FASTA_VECTOR_OK.");

    char *fastaVectorIndexFileSrc = "fastaVectorIndex.awfmi";
    struct AwFmIndex *fastaVectorIndex = NULL;

    enum AwFmReturnCode awFmRc = awFmCreateIndexFromFasta(
        &fastaVectorIndex, &config, fastaFileSrc, fastaVectorIndexFileSrc);
    testAssertString(awFmRc == AwFmFileWriteOkay,
                     "creating fastaVector index did not return AwFmSuccess.");

    // compare the fasta sequence
    sprintf(buffer,
            "sequence count %zu did not match original fastaVector count %zu.",
            fastaVectorIndex->fastaVector->sequence.count,
            fastaVector->sequence.count);
    testAssertString(fastaVectorIndex->fastaVector->sequence.count ==
                         fastaVector->sequence.count,
                     buffer);

    sprintf(buffer, "sequence count %zu greater than capacity %zu",
            fastaVectorIndex->fastaVector->sequence.count,
            fastaVectorIndex->fastaVector->sequence.capacity);
    testAssertString(fastaVectorIndex->fastaVector->sequence.capacity >=
                         fastaVectorIndex->fastaVector->sequence.count,
                     buffer);

    for (size_t letterIndex = 0; letterIndex < fastaVector->sequence.count;
         letterIndex++) {
      sprintf(buffer,
              "letter in sequence at index %zu from index %u did not match "
              "from fastaVector %u",
              letterIndex,
              fastaVectorIndex->fastaVector->sequence.charData[letterIndex],
              fastaVector->sequence.charData[letterIndex]);
      testAssertString(
          fastaVectorIndex->fastaVector->sequence.charData[letterIndex] ==
              fastaVector->sequence.charData[letterIndex],
          buffer);
    }

    // compare the header
    sprintf(buffer,
            "header count %zu did not match original fastaVector count %zu.",
            fastaVectorIndex->fastaVector->header.count,
            fastaVector->header.count);
    testAssertString(fastaVectorIndex->fastaVector->header.count ==
                         fastaVector->header.count,
                     buffer);

    sprintf(buffer, "header count %zu greater than capacity %zu",
            fastaVectorIndex->fastaVector->header.count,
            fastaVectorIndex->fastaVector->header.capacity);
    testAssertString(fastaVectorIndex->fastaVector->header.capacity >=
                         fastaVectorIndex->fastaVector->header.count,
                     buffer);

    for (size_t letterIndex = 0; letterIndex < fastaVector->header.count;
         letterIndex++) {
      sprintf(buffer,
              "letter in header at index %zu from index %u did not match from "
              "fastaVector %u",
              letterIndex,
              fastaVectorIndex->fastaVector->header.charData[letterIndex],
              fastaVector->header.charData[letterIndex]);
      testAssertString(
          fastaVectorIndex->fastaVector->header.charData[letterIndex] ==
              fastaVector->header.charData[letterIndex],
          buffer);
    }

    // compare the config
    sprintf(buffer,
            "config count %zu did not match original fastaVector count %zu.",
            fastaVectorIndex->fastaVector->metadata.count,
            fastaVector->metadata.count);
    testAssertString(fastaVectorIndex->fastaVector->metadata.count ==
                         fastaVector->metadata.count,
                     buffer);

    sprintf(buffer, "config count %zu greater than capacity %zu",
            fastaVectorIndex->fastaVector->metadata.count,
            fastaVectorIndex->fastaVector->metadata.capacity);
    testAssertString(fastaVectorIndex->fastaVector->metadata.capacity >=
                         fastaVectorIndex->fastaVector->metadata.count,
                     buffer);

    for (size_t configIndex = 0; configIndex < fastaVector->metadata.count;
         configIndex++) {
      size_t indexMetadataHeaderEnd =
          fastaVectorIndex->fastaVector->metadata.data[configIndex]
              .headerEndPosition;
      size_t indexMetadataSequenceEnd =
          fastaVectorIndex->fastaVector->metadata.data[configIndex]
              .sequenceEndPosition;
      size_t originalMetadataHeaderEnd =
          fastaVector->metadata.data[configIndex].headerEndPosition;
      size_t originalMetadataSequenceEnd =
          fastaVector->metadata.data[configIndex].sequenceEndPosition;
      sprintf(buffer,
              "header end position in config at index %zu from index %zu did "
              "not match from fastaVector %zu",
              configIndex, indexMetadataHeaderEnd, originalMetadataHeaderEnd);
      testAssertString(indexMetadataHeaderEnd == originalMetadataHeaderEnd,
                       buffer);
      sprintf(buffer,
              "sequence end position in config at index %zu from index %zu did "
              "not match from fastaVector %zu",
              configIndex, indexMetadataSequenceEnd,
              originalMetadataSequenceEnd);
      testAssertString(indexMetadataSequenceEnd == originalMetadataSequenceEnd,
                       buffer);
    }

    fastaVectorDealloc(fastaVector);
    free(fastaVector);
    awFmDeallocIndex(fastaVectorIndex);
  }
}

void testAwFmIndexGivesCorrectLocalPositions(void) {
  // generate fastavector, save to file.
  // build fm index from fastavector
  // for each position in the index bwt, ensure that the global position
  // gives the correct local position + sequence.
  for (size_t testNum = 0; testNum < 10; testNum++) {
    struct FastaVector *fastaVector = NULL;
    // generating the config also sets the test's alphabet type.
    struct AwFmIndexConfiguration config = generateReasonableRandomMetadata();
    const bool isAmino = config.alphabetType == AwFmAlphabetAmino;
    const size_t numSequences = (rand() % 10) + 1;
    const size_t maxSequenceLength = 10000;
    const char *fastaFileSrc = "sequences.fasta";
    const uint32_t fileLineLength = (rand() % 522) + 1;
    fastaVector = generateMultiSequenceFastaVector(numSequences,
                                                   maxSequenceLength, isAmino);
    printf("test gives correct locations # %zu, config generated: sacr %u, "
           "klist %u, alphabet %u, ksim %u\n",
           testNum, config.suffixArrayCompressionRatio,
           config.kmerLengthInSeedTable, config.alphabetType,
           config.keepSuffixArrayInMemory);

    enum FastaVectorReturnCode rc =
        fastaVectorWriteFasta(fastaFileSrc, fastaVector, fileLineLength);
    testAssertString(
        rc == FASTA_VECTOR_OK,
        "writing fasta vector to file did not return FASTA_VECTOR_OK.");

    char *fastaVectorIndexFileSrc = "fastaVectorIndex.awfmi";
    struct AwFmIndex *fastaVectorIndex = NULL;

    enum AwFmReturnCode awFmRc = awFmCreateIndexFromFasta(
        &fastaVectorIndex, &config, fastaFileSrc, fastaVectorIndexFileSrc);
    testAssertString(awFmRc == AwFmFileWriteOkay,
                     "creating fastaVector index did not return AwFmSuccess.");

    checkAllGlobalPositionsForCorrectLocalPositions(fastaVector,
                                                    fastaVectorIndex);

    fastaVectorDealloc(fastaVector);
    free(fastaVector);
    awFmDeallocIndex(fastaVectorIndex);
  }
}

void testAwFmIndexGivesCorrectHeaders(void) {
  for (size_t testNum = 0; testNum < 10; testNum++) {
    struct FastaVector *fastaVector = NULL;
    // generating the config also sets the test's alphabet type.
    struct AwFmIndexConfiguration config = generateReasonableRandomMetadata();
    const bool isAmino = config.alphabetType == AwFmAlphabetAmino;
    const size_t numSequences = (rand() % 10) + 1;
    const size_t maxSequenceLength = 10000;
    const char *fastaFileSrc = "sequences.fasta";
    const uint32_t fileLineLength = (rand() % 522) + 1;
    fastaVector = generateMultiSequenceFastaVector(numSequences,
                                                   maxSequenceLength, isAmino);
    printf("test gives correct headers # %zu, config generated: sacr %u, klist "
           "%u, alphabet %u, ksim %u\n",
           testNum, config.suffixArrayCompressionRatio,
           config.kmerLengthInSeedTable, config.alphabetType,
           config.keepSuffixArrayInMemory);

    enum FastaVectorReturnCode rc =
        fastaVectorWriteFasta(fastaFileSrc, fastaVector, fileLineLength);
    testAssertString(
        rc == FASTA_VECTOR_OK,
        "writing fasta vector to file did not return FASTA_VECTOR_OK.");

    char *fastaVectorIndexFileSrc = "fastaVectorIndex.awfmi";
    struct AwFmIndex *fastaVectorIndex = malloc(sizeof(struct AwFmIndex));

    enum AwFmReturnCode awFmRc = awFmCreateIndexFromFasta(
        &fastaVectorIndex, &config, fastaFileSrc, fastaVectorIndexFileSrc);
    testAssertString(awFmRc == AwFmSuccess,
                     "creating fastaVector index did not return AwFmSuccess.");

    for (size_t sequenceIndex = 0; sequenceIndex < fastaVector->metadata.count;
         sequenceIndex++) {
      char *headerBufferPtrFromFastaVector;
      if (sequenceIndex == 0) {
        headerBufferPtrFromFastaVector = fastaVector->header.charData;
      } else {
        size_t headerOffset =
            fastaVector->metadata.data[sequenceIndex - 1].headerEndPosition;
        headerBufferPtrFromFastaVector =
            fastaVector->header.charData + headerOffset;
      }
      size_t headerLengthFromFastaVector;
      if (sequenceIndex == 0) {
        headerLengthFromFastaVector =
            fastaVector->metadata.data[0].headerEndPosition;
      } else {
        size_t startPosition =
            fastaVector->metadata.data[sequenceIndex - 1].headerEndPosition;
        size_t endPosition =
            fastaVector->metadata.data[sequenceIndex].headerEndPosition;
        headerLengthFromFastaVector = endPosition - startPosition;
      }

      char *headerBuffer = NULL;
      size_t headerLength = 0;
      enum AwFmReturnCode rc = awFmGetHeaderStringFromSequenceNumber(
          fastaVectorIndex, sequenceIndex, &headerBuffer, &headerLength);

      sprintf(
          buffer,
          "get header string from sequence number did not return AwFmSuccess.");
      testAssertString(rc == AwFmSuccess, buffer);

      sprintf(buffer,
              "pointer from header buffer for sequence %zu (%p) did not match "
              "expected (%p).",
              sequenceIndex, headerBuffer, headerBufferPtrFromFastaVector);
      testAssertString(headerBuffer == headerBufferPtrFromFastaVector, buffer);

      sprintf(buffer,
              "header length for seuqence %zu (%zu) did not match expected %zu",
              sequenceIndex, headerLength, headerLengthFromFastaVector);
    }

    fastaVectorDealloc(fastaVector);
    awFmDeallocIndex(fastaVectorIndex);
  }
}

void compareIndicesForEqualityIgnoreVersion(const struct AwFmIndex *index1,
                                            const struct AwFmIndex *index2) {
  // test the easy stuff member data (the scalar values)
  sprintf(buffer, "index 1 bwt length %zu did not match index 2 length %zu.",
          index1->bwtLength, index2->bwtLength);
  testAssertString(index1->bwtLength == index2->bwtLength, buffer);

  sprintf(buffer,
          "index 1 SA compression ratio %i did not match index 2 ratio %i.",
          index1->config.suffixArrayCompressionRatio,
          index2->config.suffixArrayCompressionRatio);
  testAssertString(index1->config.suffixArrayCompressionRatio ==
                       index2->config.suffixArrayCompressionRatio,
                   buffer);

  sprintf(buffer, "index 1 SA kmer table len %i did not match index 2 len %i.",
          index1->config.kmerLengthInSeedTable,
          index2->config.kmerLengthInSeedTable);
  testAssertString(index1->config.kmerLengthInSeedTable ==
                       index2->config.kmerLengthInSeedTable,
                   buffer);

  sprintf(buffer, "index 1 alphabet type %i did not match index 2 alphabet %i.",
          index1->config.alphabetType, index2->config.alphabetType);
  testAssertString(index1->config.alphabetType == index2->config.alphabetType,
                   buffer);

  sprintf(buffer, "index 1 SA in mem bool %i did not match index 2 bool %i.",
          index1->config.keepSuffixArrayInMemory,
          index2->config.keepSuffixArrayInMemory);
  testAssertString(index1->config.keepSuffixArrayInMemory ==
                       index2->config.keepSuffixArrayInMemory,
                   buffer);

  sprintf(buffer,
          "index 1 SA file offset %zu did not match index 2 file offset %zu.",
          index1->suffixArrayFileOffset, index2->suffixArrayFileOffset);
  testAssertString(
      index1->suffixArrayFileOffset == index2->suffixArrayFileOffset, buffer);

  sprintf(
      buffer,
      "index 1 sequence file offset %zu did not match index 2 file offset %zu.",
      index1->sequenceFileOffset, index2->sequenceFileOffset);
  testAssertString(index1->sequenceFileOffset == index2->sequenceFileOffset,
                   buffer);

  // compare bwts
  const size_t numBwtBlocks = index1->bwtLength / 256;
  for (size_t blockIndex = 0; blockIndex < numBwtBlocks; blockIndex++) {
    if (index1->config.alphabetType == AwFmAlphabetDna) {
      bool blocksAreEqual =
          memcmp(&index1->bwtBlockList.asNucleotide[blockIndex],
                 &index2->bwtBlockList.asNucleotide[blockIndex],
                 sizeof(struct AwFmNucleotideBlock)) == 0;
      sprintf(buffer,
              "index 1 bwt block %zu did not compare equal to matching index 2 "
              "block.",
              blockIndex);
      testAssertString(blocksAreEqual, buffer);
      if (!blocksAreEqual) {
        struct AwFmNucleotideBlock *index1Block =
            &index1->bwtBlockList.asNucleotide[blockIndex];
        struct AwFmNucleotideBlock *index2Block =
            &index2->bwtBlockList.asNucleotide[blockIndex];
        printf("printing data for  block 1\n");
        for (size_t letterIndex = 0; letterIndex < 4; letterIndex++) {
          printf("\tcount for letter index %zu: \t %zu - %zu\n", letterIndex,
                 index1Block->baseOccurrences[letterIndex],
                 index2Block->baseOccurrences[letterIndex]);
        }
        for (size_t blockByte = 0; blockByte < 32 * 5; blockByte++) {
          if (((uint8_t *)index1Block->letterBitVectors)[blockByte] !=
              ((uint8_t *)index2Block->letterBitVectors)[blockByte]) {
            printf("byte %zu mismatch, %02X - %02X\n", blockByte,
                   ((uint8_t *)index1Block->letterBitVectors)[blockByte],
                   ((uint8_t *)index2Block->letterBitVectors)[blockByte]);
          }
        }
      }
    } else {
      bool blocksAreEqual = memcmp(&index1->bwtBlockList.asAmino[blockIndex],
                                   &index2->bwtBlockList.asAmino[blockIndex],
                                   sizeof(struct AwFmAminoBlock)) == 0;
      sprintf(buffer,
              "index 1 bwt block %zu did not compare equal to matching index 2 "
              "block.",
              blockIndex);
      testAssertString(blocksAreEqual, buffer);
    }
  }

  // compare prefix sums
  const size_t numPrefixSums =
      index1->config.alphabetType == AwFmAlphabetDna ? 6 : 22;
  for (size_t prefixSumIndex = 0; prefixSumIndex < numPrefixSums;
       prefixSumIndex++) {
    sprintf(buffer,
            "index 1 prefix sum #%zu (%zu) did not compare equal to matching "
            "index 2 prefix sum (%zu).",
            prefixSumIndex, index1->prefixSums[prefixSumIndex],
            index2->prefixSums[prefixSumIndex]);
    testAssertString(index1->prefixSums[prefixSumIndex] ==
                         index2->prefixSums[prefixSumIndex],
                     buffer);
  }

  // //compare the kmer tables
  // size_t kmerTableLength = awFmGetKmerTableLength(index1);
  // for(size_t kmerIndex = 0; kmerIndex < kmerTableLength; kmerIndex++){
  //   sprintf(buffer, "kmer table at index %zu did not match (index 1 [%zu,
  //   %zu], index 2 [%zu, %zu])",
  //     kmerIndex, index1->kmerSeedTable[kmerIndex].startPtr,
  //     index1->kmerSeedTable[kmerIndex].endPtr,
  //     index2->kmerSeedTable[kmerIndex].startPtr,
  //     index2->kmerSeedTable[kmerIndex].endPtr);
  //
  //   bool kmerTableItemMatches =  index1->kmerSeedTable[kmerIndex].startPtr ==
  //   index2->kmerSeedTable[kmerIndex].startPtr &&
  //     index1->kmerSeedTable[kmerIndex].endPtr ==
  //     index2->kmerSeedTable[kmerIndex].endPtr;
  //   testAssertString(kmerTableItemMatches, buffer);
  // }

  // compare the suffix arrays
  size_t compressedSaLength =
      index1->bwtLength / index1->config.suffixArrayCompressionRatio;
  for (size_t saIndex = 0; saIndex < compressedSaLength; saIndex++) {
    size_t index1Position =
        awFmGetValueFromCompressedSuffixArray(&index1->suffixArray, saIndex);
    size_t index2Position =
        awFmGetValueFromCompressedSuffixArray(&index2->suffixArray, saIndex);
    sprintf(
        buffer,
        "suffix array at position %zu did not match (index 1 %zu, index 2 %zu)",
        saIndex, index1Position, index2Position);
    bool saElementMatches = index1Position == index2Position;
    testAssertString(saElementMatches, buffer);
  }
}

void checkAllGlobalPositionsForCorrectLocalPositions(
    const struct FastaVector *_RESTRICT_ const fastaVector,
    const struct AwFmIndex *_RESTRICT_ const fastaVectorIndex) {

  for (size_t sequenceIndex = 0; sequenceIndex < fastaVector->metadata.count;
       sequenceIndex++) {
    size_t sequenceBeginPosition =
        sequenceIndex == 0
            ? 0
            : fastaVector->metadata.data[sequenceIndex - 1].sequenceEndPosition;
    size_t sequenceEndPosition =
        fastaVector->metadata.data[sequenceIndex].sequenceEndPosition;

    size_t sequenceLength = sequenceEndPosition - sequenceBeginPosition;
    for (size_t localPosition = 0; localPosition < sequenceLength;
         localPosition++) {
      const size_t globalPosition = sequenceBeginPosition + localPosition;
      size_t sequenceNumberOut = 0;
      size_t localSequencePositionOut = 0;
      enum AwFmReturnCode rc = awFmGetLocalSequencePositionFromIndexPosition(
          fastaVectorIndex, globalPosition, &sequenceNumberOut,
          &localSequencePositionOut);

      sprintf(buffer, "getting local position did not return code OK.");
      testAssertString(rc == AwFmSuccess, buffer);

      sprintf(buffer,
              "localSequencePosition expected %zu, but got %zu from global "
              "position %zu.",
              localSequencePositionOut, localPosition, globalPosition);
      testAssertString(localSequencePositionOut == localPosition, buffer);

      sprintf(buffer,
              "sequenceNumber %zu did not match expected %zu for global "
              "position %zu.",
              sequenceNumberOut, sequenceIndex, globalPosition);
      testAssertString(sequenceNumberOut == sequenceIndex, buffer);
    }
  }
}

void testFromFasta() {
  size_t localPosition, sequenceNumber;
  uint64_t *positions;
  struct AwFmIndex *index;
  struct AwFmIndexConfiguration config;
  config.suffixArrayCompressionRatio = 2;
  config.kmerLengthInSeedTable = 2;
  config.alphabetType = AwFmAlphabetAmino;
  config.keepSuffixArrayInMemory = true;
  config.storeOriginalSequence = false;
  enum AwFmReturnCode awfmrc =
      awFmCreateIndexFromFasta(&index, &config, "test2.fa", "test2.awfmi");

  testAssertString(awfmrc == AwFmFileWriteOkay,
                   "Error: create from fasta did not return file write okay");
  struct AwFmSearchRange searchRange;

  bool foundHeader = awFmSingleKmerExists(index, "t", 1);
  testAssertString((!foundHeader),
                   "Error: awfm found header 1 while looking for sequence");

  foundHeader = awFmSingleKmerExists(index, "v", 1);
  testAssertString((!foundHeader),
                   "Error: awfm found header 2 while looking for sequence");

  foundHeader = awFmSingleKmerExists(index, "w", 1);
  testAssertString((!foundHeader),
                   "Error: awfm found header 2 while looking for sequence");

  foundHeader = awFmSingleKmerExists(index, "y", 1);
  testAssertString((!foundHeader),
                   "Error: awfm found header 2 while looking for sequence");

  bool foundSequence = awFmSingleKmerExists(index, "acdef", 5);
  testAssertString(foundSequence, "Error: did not find sequence 1");

  foundSequence = awFmSingleKmerExists(index, "g", 1);
  testAssertString(foundSequence, "Error: did not find sequence 2");

  foundSequence = awFmSingleKmerExists(index, "hikl", 4);
  testAssertString(foundSequence, "Error: did not find sequence 3");

  foundSequence = awFmSingleKmerExists(index, "m", 1);
  testAssertString(foundSequence, "Error: did not find sequence 1");

  // check to make sure they're in the right position;
  awFmAminoNonSeededSearch(index, "acdef", 5, &searchRange);
  testAssertString(awFmSearchRangeLength(&searchRange) == 1,
                   "Error: seq 1 did not return exactly 1 hit");
  awFmAminoNonSeededSearch(index, "g", 1, &searchRange);
  testAssertString(awFmSearchRangeLength(&searchRange) == 1,
                   "Error: seq 2 did not return exactly 1 hit");
  awFmAminoNonSeededSearch(index, "hikl", 4, &searchRange);
  testAssertString(awFmSearchRangeLength(&searchRange) == 1,
                   "Error: seq 3 did not return exactly 1 hit");
  awFmAminoNonSeededSearch(index, "m", 1, &searchRange);
  testAssertString(awFmSearchRangeLength(&searchRange) == 1,
                   "Error: seq 4 did not return exactly 1 hit");

  // check to make sure the local positions make sense
  awFmAminoNonSeededSearch(index, "acdef", 5, &searchRange);
  positions = awFmFindDatabaseHitPositions(index, &searchRange, &awfmrc);
  testAssertString(awfmrc == AwFmFileReadOkay,
                   "getting db hit position did not return read okay.");
  awfmrc = awFmGetLocalSequencePositionFromIndexPosition(
      index, positions[0], &sequenceNumber, &localPosition);
  testAssertString(awfmrc == AwFmSuccess,
                   "get local position did not return success");
  testAssertString(sequenceNumber == 0,
                   "first seq did not return seq number 0");
  testAssertString(localPosition == 0,
                   "first seq did not return local position 0");
  free(positions);

  awFmAminoNonSeededSearch(index, "g", 1, &searchRange);
  positions = awFmFindDatabaseHitPositions(index, &searchRange, &awfmrc);
  testAssertString(awfmrc == AwFmFileReadOkay,
                   "getting db hit position did not return read okay.");
  awfmrc = awFmGetLocalSequencePositionFromIndexPosition(
      index, positions[0], &sequenceNumber, &localPosition);
  testAssertString(awfmrc == AwFmSuccess,
                   "get local position did not return success");
  testAssertString(sequenceNumber == 1,
                   "second seq did not return seq number 1");
  testAssertString(localPosition == 0,
                   "second seq did not return local position 0");
  free(positions);

  awFmAminoNonSeededSearch(index, "hikl", 4, &searchRange);
  positions = awFmFindDatabaseHitPositions(index, &searchRange, &awfmrc);
  testAssertString(awfmrc == AwFmFileReadOkay,
                   "getting db hit position did not return read okay.");
  awfmrc = awFmGetLocalSequencePositionFromIndexPosition(
      index, positions[0], &sequenceNumber, &localPosition);
  testAssertString(awfmrc == AwFmSuccess,
                   "get local position did not return success");
  testAssertString(sequenceNumber == 2,
                   "third seq did not return seq number 2");
  testAssertString(localPosition == 0,
                   "third seq did not return local position 0");
  free(positions);

  awFmAminoNonSeededSearch(index, "m", 1, &searchRange);
  positions = awFmFindDatabaseHitPositions(index, &searchRange, &awfmrc);
  testAssertString(awfmrc == AwFmFileReadOkay,
                   "getting db hit position did not return read okay.");
  awfmrc = awFmGetLocalSequencePositionFromIndexPosition(
      index, positions[0], &sequenceNumber, &localPosition);
  testAssertString(awfmrc == AwFmSuccess,
                   "get local position did not return success");
  testAssertString(sequenceNumber == 3,
                   "fourth seq did not return seq number 3");
  testAssertString(localPosition == 0,
                   "fourth seq did not return local position 0");
  free(positions);

  // check to make sure no hits overlap the sequences

  awFmAminoNonSeededSearch(index, "fg", 2, &searchRange);
  testAssertString(awFmSearchRangeLength(&searchRange) == 0,
                   "searching for overlap between seqs 1 and 2 found a hit");
  awFmAminoNonSeededSearch(index, "gh", 2, &searchRange);
  testAssertString(awFmSearchRangeLength(&searchRange) == 0,
                   "searching for overlap between seqs 2 and 3 found a hit");
  awFmAminoNonSeededSearch(index, "lm", 2, &searchRange);
  testAssertString(awFmSearchRangeLength(&searchRange) == 0,
                   "searching for overlap between seqs 3 and 4 found a hit");
}