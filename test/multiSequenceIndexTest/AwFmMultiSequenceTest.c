#include <string.h>
#include "../../src/AwFmIndex.h"
#include "../../src/AwFmIndexStruct.h"
#include "../../src/AwFmOccurrence.h"
#include "../../src/AwFmCreate.h"
#include "../../src/AwFmSearch.h"
#include "../../src/AwFmKmerTable.h"
#include "../../src/AwFmLetter.h"
#include "../../src/AwFmParallelSearch.h"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include "FastaVector.h"
// #include "../../libdivsufsort/build/include/divsufsort64.h"



char buffer[2048];
char aminoLookup[21]     = {'a','c','d','e','f',
                              'g','h','i','k','l',
                              'm','n','p','q','r',
                              's','t','v','w','y','x'};
char nucleotideLookup[5] = {'a','g','c','t','x'};

struct FastaVector* generateMultiSequenceFastaVector(const size_t numSequences,
  const size_t maxSequenceLength,const bool isAmino);

void randomizeSequenceBuffer(char *sequenceBuffer, size_t sequenceLength, bool isAmino);
void getRawSequenceFromFastaVector(const struct FastaVector *const fastaVector,
  char **compositeSequenceBufferPtr, size_t *compositeSequenceLength);

struct AwFmIndexMetadata generateReasonableRandomMetadata();

void testAwFmIndexIdenticalForFastaVector(void);
void testAwFmIndexFastaVectorDataMatchesExpected(void);
void testAwFmIndexGivesCorrectLocalPositions(void);
void testAwFmIndexGivesCorrectHeaders(void);
void compareIndicesForEqualityIgnoreVersion(const struct AwFmIndex *index1, const struct AwFmIndex *index2);
void checkAllGlobalPositionsForCorrectLocalPositions(const struct FastaVector *restrict const fastaVector,
  const struct AwFmIndex *restrict const fastaVectorIndex);


int main (int argc, char **argv){
  srand(time(NULL));
  testAwFmIndexIdenticalForFastaVector();
  testAwFmIndexFastaVectorDataMatchesExpected();
  testAwFmIndexGivesCorrectLocalPositions();
}

struct FastaVector* generateMultiSequenceFastaVector(const size_t numSequences,
  const size_t maxSequenceLength,const bool isAmino){

  char headerBuffer[128] = {0};
  struct FastaVector *fastaVector = malloc(sizeof(struct FastaVector));
  testAssertString(fastaVector != NULL, "in generateMultiSequenceFastaVector, fasta vector malloc failed (returned NULL).");

  enum FastaVectorReturnCode rc = fastaVectorInit(fastaVector);
  testAssertString(rc == FASTA_VECTOR_OK, "in generateMultiSequenceFastaVector, fastaVector init did not return OK.");

  for(size_t sequenceIndex = 0; sequenceIndex < numSequences; sequenceIndex++){
    size_t thisSequenceLength = (rand()% maxSequenceLength) + 1;
    char *sequenceBuffer = malloc(thisSequenceLength * sizeof(char));
    randomizeSequenceBuffer(sequenceBuffer, thisSequenceLength, isAmino);
    sprintf(headerBuffer, "header%zu| metadata", sequenceIndex);
    rc = fastaVectorAddSequenceToList(fastaVector, headerBuffer,
      strlen(headerBuffer), sequenceBuffer, thisSequenceLength);
    testAssertString(rc == FASTA_VECTOR_OK, "in generateMultiSequenceFastaVector, fastaVector adding sequence to list did not return OK.");
    free(sequenceBuffer);
  }


  return fastaVector;
}


void randomizeSequenceBuffer(char *sequenceBuffer, size_t sequenceLength, bool isAmino){
  const char *const letterLookupTable = isAmino? aminoLookup: nucleotideLookup;
  const size_t letterLookupTableLength = isAmino? 21: 5;
  for(size_t letterIndex = 0; letterIndex < sequenceLength; letterIndex++){
    const uint8_t chosenLetter = rand() % letterLookupTableLength;
    sequenceBuffer[letterIndex] = letterLookupTable[chosenLetter];
  }
}

void getRawSequenceFromFastaVector(const struct FastaVector *const fastaVector,
  char **compositeSequenceBufferPtr, size_t *compositeSequenceLength){

  *compositeSequenceLength =  fastaVector->sequence.count;
  *compositeSequenceBufferPtr = fastaVector->sequence.charData;
}


struct AwFmIndexMetadata generateReasonableRandomMetadata(){
  struct AwFmIndexMetadata metadata;
  metadata.suffixArrayCompressionRatio = (rand() % 20) + 1;
  metadata.alphabetType = (rand() % 2) == 0? AwFmAlphabetAmino : AwFmAlphabetNucleotide;
  metadata.kmerLengthInSeedTable = metadata.alphabetType == AwFmAlphabetNucleotide?
    (rand()%10) + 2:
    (rand()%4) + 1;
  metadata.keepSuffixArrayInMemory = true;

  return metadata;
}

void testAwFmIndexIdenticalForFastaVector(){
  for(size_t testNum = 0; testNum < 10; testNum++){
    struct FastaVector *fastaVector = NULL;
    //generating the metadata also sets the test's alphabet type.
    struct AwFmIndexMetadata metadata = generateReasonableRandomMetadata();
    printf("test identical # %zu, metadata generated: sacr %u, klist %u, alphabet %u, ksim %u\n", testNum, metadata.suffixArrayCompressionRatio, metadata.kmerLengthInSeedTable, metadata.alphabetType, metadata.keepSuffixArrayInMemory);


    const bool isAmino              = metadata.alphabetType == AwFmAlphabetAmino;
    const size_t numSequences       = (rand() % 10) + 1;
    const size_t maxSequenceLength  = 10000;
    const char *fastaFileSrc        = "sequences.fasta";
    const uint32_t fileLineLength   = (rand() % 522) + 1;
    fastaVector = generateMultiSequenceFastaVector(numSequences, maxSequenceLength, isAmino);


  enum FastaVectorReturnCode rc = fastaVectorWriteFasta(fastaFileSrc, fastaVector, fileLineLength);
  testAssertString(rc == FASTA_VECTOR_OK, "writing fasta vector to file did not return FASTA_VECTOR_OK.");


  char *sequenceCollectionPtr = NULL;
  size_t sequenceCollectionLength = 0;
  getRawSequenceFromFastaVector(fastaVector, &sequenceCollectionPtr, &sequenceCollectionLength);
  testAssertString(sequenceCollectionPtr != NULL, "sequence collection ptr was null!");

  char *defaultIndexFileSrc           = "defaultIndex.awfmi";
  char *fastaVectorIndexFileSrc       = "fastaVectorIndex.awfmi";
  struct AwFmIndex* defaultIndex      = NULL;
  struct AwFmIndex* fastaVectorIndex  = NULL;

  enum AwFmReturnCode awFmRc = awFmCreateIndex(&defaultIndex, &metadata, (uint8_t*)sequenceCollectionPtr,
    sequenceCollectionLength, defaultIndexFileSrc, true);
  testAssertString(awFmRc == AwFmFileWriteOkay, "creating default index did not return AwFmSuccess.");
  awFmRc = awFmCreateIndexFromFasta(&fastaVectorIndex, &metadata, fastaFileSrc, fastaVectorIndexFileSrc, true);
  testAssertString(awFmRc == AwFmFileWriteOkay, "creating fastaVector index did not return AwFmSuccess");

  compareIndicesForEqualityIgnoreVersion(defaultIndex, fastaVectorIndex);
  //TODO! why is fastVectorIndex having zeros in most of it's data?

  fastaVectorDealloc(fastaVector);
  free(fastaVector);
  awFmDeallocIndex(fastaVectorIndex);
  awFmDeallocIndex(defaultIndex);
  }
}


void testAwFmIndexFastaVectorDataMatchesExpected(void){
  for(size_t testNum = 0; testNum < 10; testNum++){
    struct FastaVector *fastaVector = NULL;
    //generating the metadata also sets the test's alphabet type.
    struct AwFmIndexMetadata metadata = generateReasonableRandomMetadata();
    const bool isAmino              = metadata.alphabetType == AwFmAlphabetAmino;
    const size_t numSequences       = (rand() % 10) + 1;
    const size_t maxSequenceLength  = 10000;
    const char *fastaFileSrc        = "sequences.fasta";
    const uint32_t fileLineLength   = (rand() % 522) + 1;
    fastaVector = generateMultiSequenceFastaVector(numSequences, maxSequenceLength, isAmino);
    printf("test fastaVector data # %zu, metadata generated: sacr %u, klist %u, alphabet %u, ksim %u\n", testNum, metadata.suffixArrayCompressionRatio, metadata.kmerLengthInSeedTable, metadata.alphabetType, metadata.keepSuffixArrayInMemory);



  enum FastaVectorReturnCode rc = fastaVectorWriteFasta(fastaFileSrc, fastaVector, fileLineLength);
  testAssertString(rc == FASTA_VECTOR_OK, "writing fasta vector to file did not return FASTA_VECTOR_OK.");

  char *fastaVectorIndexFileSrc       = "fastaVectorIndex.awfmi";
  struct AwFmIndex* fastaVectorIndex  = NULL;

  enum AwFmReturnCode awFmRc = awFmCreateIndexFromFasta(&fastaVectorIndex, &metadata, fastaFileSrc, fastaVectorIndexFileSrc, true);
  testAssertString(awFmRc == AwFmFileWriteOkay, "creating fastaVector index did not return AwFmSuccess.");

  //compare the fasta sequence
  sprintf(buffer, "sequence count %zu did not match original fastaVector count %zu.", fastaVectorIndex->fastaVector->sequence.count, fastaVector->sequence.count);
  testAssertString(fastaVectorIndex->fastaVector->sequence.count == fastaVector->sequence.count, buffer);

  sprintf(buffer, "sequence count %zu greater than capacity %zu", fastaVectorIndex->fastaVector->sequence.count, fastaVectorIndex->fastaVector->sequence.capacity);
  testAssertString(fastaVectorIndex->fastaVector->sequence.capacity >= fastaVectorIndex->fastaVector->sequence.count, buffer);

  for(size_t letterIndex = 0; letterIndex < fastaVector->sequence.count; letterIndex++){
    sprintf(buffer, "letter in sequence at index %zu from index %u did not match from fastaVector %u", letterIndex, fastaVectorIndex->fastaVector->sequence.charData[letterIndex], fastaVector->sequence.charData[letterIndex]);
    testAssertString( fastaVectorIndex->fastaVector->sequence.charData[letterIndex] == fastaVector->sequence.charData[letterIndex], buffer);
  }

  //compare the header
  sprintf(buffer, "header count %zu did not match original fastaVector count %zu.", fastaVectorIndex->fastaVector->header.count, fastaVector->header.count);
  testAssertString(fastaVectorIndex->fastaVector->header.count == fastaVector->header.count, buffer);

  sprintf(buffer, "header count %zu greater than capacity %zu", fastaVectorIndex->fastaVector->header.count, fastaVectorIndex->fastaVector->header.capacity);
  testAssertString(fastaVectorIndex->fastaVector->header.capacity >= fastaVectorIndex->fastaVector->header.count, buffer);

  for(size_t letterIndex = 0; letterIndex < fastaVector->header.count; letterIndex++){
    sprintf(buffer, "letter in header at index %zu from index %u did not match from fastaVector %u", letterIndex, fastaVectorIndex->fastaVector->header.charData[letterIndex], fastaVector->header.charData[letterIndex]);
    testAssertString( fastaVectorIndex->fastaVector->header.charData[letterIndex] == fastaVector->header.charData[letterIndex], buffer);
  }


  //compare the metadata
  sprintf(buffer, "metadata count %zu did not match original fastaVector count %zu.", fastaVectorIndex->fastaVector->metadata.count, fastaVector->metadata.count);
  testAssertString(fastaVectorIndex->fastaVector->metadata.count == fastaVector->metadata.count, buffer);

  sprintf(buffer, "metadata count %zu greater than capacity %zu", fastaVectorIndex->fastaVector->metadata.count, fastaVectorIndex->fastaVector->metadata.capacity);
  testAssertString(fastaVectorIndex->fastaVector->metadata.capacity >= fastaVectorIndex->fastaVector->metadata.count, buffer);

  for(size_t metadataIndex = 0; metadataIndex < fastaVector->metadata.count; metadataIndex++){
    size_t indexMetadataHeaderEnd = fastaVectorIndex->fastaVector->metadata.data[metadataIndex].headerEndPosition;
    size_t indexMetadataSequenceEnd = fastaVectorIndex->fastaVector->metadata.data[metadataIndex].sequenceEndPosition;
    size_t originalMetadataHeaderEnd = fastaVector->metadata.data[metadataIndex].headerEndPosition;
    size_t originalMetadataSequenceEnd = fastaVector->metadata.data[metadataIndex].sequenceEndPosition;
    sprintf(buffer, "header end position in metadata at index %zu from index %zu did not match from fastaVector %zu", metadataIndex, indexMetadataHeaderEnd, originalMetadataHeaderEnd);
    testAssertString(indexMetadataHeaderEnd == originalMetadataHeaderEnd, buffer);
    sprintf(buffer, "sequence end position in metadata at index %zu from index %zu did not match from fastaVector %zu", metadataIndex, indexMetadataSequenceEnd, originalMetadataSequenceEnd);
    testAssertString(indexMetadataSequenceEnd == originalMetadataSequenceEnd, buffer);
  }

  fastaVectorDealloc(fastaVector);
  free(fastaVector);
  awFmDeallocIndex(fastaVectorIndex);
  }
}


void testAwFmIndexGivesCorrectLocalPositions(void){
  //generate fastavector, save to file.
  //build fm index from fastavector
  //for each position in the index bwt, ensure that the global position
  // gives the correct local position + sequence.
  for(size_t testNum = 0; testNum < 10; testNum++){
    struct FastaVector *fastaVector = NULL;
    //generating the metadata also sets the test's alphabet type.
    struct AwFmIndexMetadata metadata = generateReasonableRandomMetadata();
    const bool isAmino              = metadata.alphabetType == AwFmAlphabetAmino;
    const size_t numSequences       = (rand() % 10) + 1;
    const size_t maxSequenceLength  = 10000;
    const char *fastaFileSrc        = "sequences.fasta";
    const uint32_t fileLineLength   = (rand() % 522) + 1;
    fastaVector = generateMultiSequenceFastaVector(numSequences, maxSequenceLength, isAmino);
    printf("test gives correct locations # %zu, metadata generated: sacr %u, klist %u, alphabet %u, ksim %u\n", testNum, metadata.suffixArrayCompressionRatio, metadata.kmerLengthInSeedTable, metadata.alphabetType, metadata.keepSuffixArrayInMemory);


    enum FastaVectorReturnCode rc = fastaVectorWriteFasta(fastaFileSrc, fastaVector, fileLineLength);
    testAssertString(rc == FASTA_VECTOR_OK, "writing fasta vector to file did not return FASTA_VECTOR_OK.");


    char *fastaVectorIndexFileSrc       = "fastaVectorIndex.awfmi";
    struct AwFmIndex* fastaVectorIndex  = NULL;

    enum AwFmReturnCode awFmRc = awFmCreateIndexFromFasta(&fastaVectorIndex, &metadata, fastaFileSrc, fastaVectorIndexFileSrc, true);
    testAssertString(awFmRc == AwFmFileWriteOkay, "creating fastaVector index did not return AwFmSuccess.");

    checkAllGlobalPositionsForCorrectLocalPositions(fastaVector, fastaVectorIndex);


    fastaVectorDealloc(fastaVector);
    free(fastaVector);
    awFmDeallocIndex(fastaVectorIndex);
  }
}


void testAwFmIndexGivesCorrectHeaders(void){
  for(size_t testNum = 0; testNum < 10; testNum++){
    struct FastaVector *fastaVector = NULL;
    //generating the metadata also sets the test's alphabet type.
    struct AwFmIndexMetadata metadata = generateReasonableRandomMetadata();
    const bool isAmino              = metadata.alphabetType == AwFmAlphabetAmino;
    const size_t numSequences       = (rand() % 10) + 1;
    const size_t maxSequenceLength  = 10000;
    const char *fastaFileSrc        = "sequences.fasta";
    const uint32_t fileLineLength   = (rand() % 522) + 1;
    fastaVector = generateMultiSequenceFastaVector(numSequences, maxSequenceLength, isAmino);
    printf("test gives correct headers # %zu, metadata generated: sacr %u, klist %u, alphabet %u, ksim %u\n", testNum, metadata.suffixArrayCompressionRatio, metadata.kmerLengthInSeedTable, metadata.alphabetType, metadata.keepSuffixArrayInMemory);

    enum FastaVectorReturnCode rc = fastaVectorWriteFasta(fastaFileSrc, fastaVector, fileLineLength);
    testAssertString(rc == FASTA_VECTOR_OK, "writing fasta vector to file did not return FASTA_VECTOR_OK.");


    char *fastaVectorIndexFileSrc       = "fastaVectorIndex.awfmi";
    struct AwFmIndex* fastaVectorIndex  = malloc(sizeof(struct AwFmIndex));

    enum AwFmReturnCode awFmRc = awFmCreateIndexFromFasta(&fastaVectorIndex, &metadata, fastaFileSrc, fastaVectorIndexFileSrc, true);
    testAssertString(awFmRc == AwFmSuccess, "creating fastaVector index did not return AwFmSuccess.");

    for(size_t sequenceIndex = 0; sequenceIndex < fastaVector->metadata.count; sequenceIndex++){
      char *headerBufferPtrFromFastaVector;
      if(sequenceIndex == 0){
        headerBufferPtrFromFastaVector = fastaVector->header.charData;
      }
      else{
        size_t headerOffset = fastaVector->metadata.data[sequenceIndex-1].headerEndPosition;
        headerBufferPtrFromFastaVector = fastaVector->header.charData + headerOffset;
      }
      size_t headerLengthFromFastaVector;
      if(sequenceIndex == 0){
        headerLengthFromFastaVector = fastaVector->metadata.data[0].headerEndPosition;
      }
      else{
        size_t startPosition = fastaVector->metadata.data[sequenceIndex-1].headerEndPosition;
        size_t endPosition = fastaVector->metadata.data[sequenceIndex].headerEndPosition;
        headerLengthFromFastaVector = endPosition - startPosition;
      }

      char *headerBuffer = NULL;
      size_t headerLength = 0 ;
      enum AwFmReturnCode rc = awFmGetHeaderStringFromSequenceNumber(fastaVectorIndex,
        sequenceIndex, &headerBuffer, &headerLength);

      sprintf(buffer, "get header string from sequence number did not return AwFmSuccess.");
      testAssertString(rc == AwFmSuccess, buffer);

      sprintf(buffer, "pointer from header buffer for sequence %zu (%p) did not match expected (%p).",
        sequenceIndex, headerBuffer, headerBufferPtrFromFastaVector);
      testAssertString(headerBuffer == headerBufferPtrFromFastaVector, buffer);

      sprintf(buffer, "header length for seuqence %zu (%zu) did not match expected %zu",
        sequenceIndex, headerLength, headerLengthFromFastaVector);
    }



    fastaVectorDealloc(fastaVector);
    awFmDeallocIndex(fastaVectorIndex);
  }
}

void compareIndicesForEqualityIgnoreVersion(const struct AwFmIndex *index1, const struct AwFmIndex *index2){
  //test the easy stuff member data (the scalar values)
  sprintf(buffer, "index 1 bwt length %zu did not match index 2 length %zu.", index1->bwtLength, index2->bwtLength);
  testAssertString(index1->bwtLength == index2->bwtLength, buffer);

  sprintf(buffer,  "index 1 SA compression ratio %i did not match index 2 ratio %i.", index1->metadata.suffixArrayCompressionRatio, index2->metadata.suffixArrayCompressionRatio);
  testAssertString(index1->metadata.suffixArrayCompressionRatio == index2->metadata.suffixArrayCompressionRatio, buffer);


  sprintf(buffer,  "index 1 SA kmer table len %i did not match index 2 len %i.", index1->metadata.kmerLengthInSeedTable, index2->metadata.kmerLengthInSeedTable);
  testAssertString(index1->metadata.kmerLengthInSeedTable == index2->metadata.kmerLengthInSeedTable, buffer);


  sprintf(buffer,  "index 1 alphabet type %i did not match index 2 alphabet %i.", index1->metadata.alphabetType, index2->metadata.alphabetType);
  testAssertString(index1->metadata.alphabetType == index2->metadata.alphabetType, buffer);


  sprintf(buffer,  "index 1 SA in mem bool %i did not match index 2 bool %i.", index1->metadata.keepSuffixArrayInMemory, index2->metadata.keepSuffixArrayInMemory);
  testAssertString(index1->metadata.keepSuffixArrayInMemory == index2->metadata.keepSuffixArrayInMemory, buffer);

  sprintf(buffer,  "index 1 SA file offset %zu did not match index 2 file offset %zu.", index1->suffixArrayFileOffset, index2->suffixArrayFileOffset);
  testAssertString(index1->suffixArrayFileOffset == index2->suffixArrayFileOffset, buffer);

  sprintf(buffer,  "index 1 sequence file offset %zu did not match index 2 file offset %zu.", index1->sequenceFileOffset, index2->sequenceFileOffset);
  testAssertString(index1->sequenceFileOffset == index2->sequenceFileOffset, buffer);


  //compare bwts
  const size_t numBwtBlocks = index1->bwtLength / 256;
  for(size_t blockIndex = 0; blockIndex < numBwtBlocks; blockIndex++){
    if(index1->metadata.alphabetType == AwFmAlphabetNucleotide){
      bool blocksAreEqual = memcmp(&index1->bwtBlockList.asNucleotide[blockIndex], &index2->bwtBlockList.asNucleotide[blockIndex], sizeof(struct AwFmNucleotideBlock)) == 0;
      sprintf(buffer, "index 1 bwt block %zu did not compare equal to matching index 2 block.", blockIndex);
      testAssertString(blocksAreEqual, buffer);
      if(!blocksAreEqual){
        struct AwFmNucleotideBlock *index1Block = &index1->bwtBlockList.asNucleotide[blockIndex];
        struct AwFmNucleotideBlock *index2Block = &index2->bwtBlockList.asNucleotide[blockIndex];
        printf("printing data for  block 1\n");
        for(size_t letterIndex = 0; letterIndex < 4; letterIndex++){
          printf("\tcount for letter index %zu: \t %zu - %zu\n", letterIndex,
            index1Block->baseOccurrences[letterIndex], index2Block->baseOccurrences[letterIndex]);
        }
        for(size_t blockByte = 0; blockByte < 32 * 5; blockByte++){
          if(((uint8_t*)index1Block->letterBitVectors)[blockByte] != ((uint8_t*)index2Block->letterBitVectors)[blockByte]){
            printf("byte %zu mismatch, %02X - %02X\n", blockByte, ((uint8_t*)index1Block->letterBitVectors)[blockByte], ((uint8_t*)index2Block->letterBitVectors)[blockByte]);
          }
        }
      }
    }
    else{
      bool blocksAreEqual = memcmp(&index1->bwtBlockList.asAmino[blockIndex], &index2->bwtBlockList.asAmino[blockIndex], sizeof(struct AwFmAminoBlock)) == 0;
      sprintf(buffer, "index 1 bwt block %zu did not compare equal to matching index 2 block.", blockIndex);
      testAssertString(blocksAreEqual, buffer);
    }
  }

  //compare prefix sums
  const size_t numPrefixSums = index1->metadata.alphabetType == AwFmAlphabetNucleotide? 6 :22;
  for(size_t prefixSumIndex = 0; prefixSumIndex < numPrefixSums; prefixSumIndex++){
    sprintf(buffer, "index 1 prefix sum #%zu (%zu) did not compare equal to matching index 2 prefix sum (%zu).", prefixSumIndex, index1->prefixSums[prefixSumIndex], index2->prefixSums[prefixSumIndex]);
    testAssertString(index1->prefixSums[prefixSumIndex] == index2->prefixSums[prefixSumIndex], buffer);
  }

  // //compare the kmer tables
  // size_t kmerTableLength = awFmGetKmerTableLength(index1);
  // for(size_t kmerIndex = 0; kmerIndex < kmerTableLength; kmerIndex++){
  //   sprintf(buffer, "kmer table at index %zu did not match (index 1 [%zu, %zu], index 2 [%zu, %zu])",
  //     kmerIndex, index1->kmerSeedTable[kmerIndex].startPtr, index1->kmerSeedTable[kmerIndex].endPtr,
  //     index2->kmerSeedTable[kmerIndex].startPtr, index2->kmerSeedTable[kmerIndex].endPtr);
  //
  //   bool kmerTableItemMatches =  index1->kmerSeedTable[kmerIndex].startPtr == index2->kmerSeedTable[kmerIndex].startPtr &&
  //     index1->kmerSeedTable[kmerIndex].endPtr == index2->kmerSeedTable[kmerIndex].endPtr;
  //   testAssertString(kmerTableItemMatches, buffer);
  // }

  //compare the suffix arrays
  size_t compressedSaLength = index1->bwtLength / index1->metadata.suffixArrayCompressionRatio;
  for(size_t saIndex = 0; saIndex < compressedSaLength; saIndex++){
    sprintf(buffer, "suffix array at position %zu did not match (index 1 %zu, index 2 %zu)", saIndex, index1->inMemorySuffixArray[saIndex], index2->inMemorySuffixArray[saIndex]);
    bool saElementMatches = index1->inMemorySuffixArray[saIndex] == index2->inMemorySuffixArray[saIndex];
    testAssertString(saElementMatches, buffer);
  }
}

void checkAllGlobalPositionsForCorrectLocalPositions(const struct FastaVector *restrict const fastaVector,
  const struct AwFmIndex *restrict const fastaVectorIndex){

  for(size_t sequenceIndex = 0; sequenceIndex < fastaVector->metadata.count; sequenceIndex++){
    size_t sequenceBeginPosition = sequenceIndex == 0? 0: fastaVector->metadata.data[sequenceIndex-1].sequenceEndPosition;
    size_t sequenceEndPosition = fastaVector->metadata.data[sequenceIndex].sequenceEndPosition;

    size_t sequenceLength = sequenceEndPosition- sequenceBeginPosition;
    for(size_t localPosition = 0; localPosition < sequenceLength; localPosition++){
      const size_t globalPosition = sequenceBeginPosition + localPosition;
      size_t sequenceNumberOut = 0;
      size_t localSequencePositionOut = 0;
      enum AwFmReturnCode rc = awFmGetLocalSequencePositionFromIndexPosition(fastaVectorIndex,
        globalPosition, &sequenceNumberOut, &localSequencePositionOut);

      sprintf(buffer, "getting local position did not return code OK.");
      testAssertString(rc == AwFmSuccess, buffer);

      sprintf(buffer, "localSequencePosition expected %zu, but got %zu from global position %zu.",
        localSequencePositionOut, localPosition, globalPosition);
      testAssertString(localSequencePositionOut == localPosition, buffer);

      sprintf(buffer, "sequenceNumber %zu did not match expected %zu for global position %zu.",
        sequenceNumberOut, sequenceIndex, globalPosition);
      testAssertString(sequenceNumberOut == sequenceIndex, buffer);
    }
  }
}
