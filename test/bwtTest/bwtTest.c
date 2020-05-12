#include "../../src/AwFmIndexStruct.h"
#include "../../src/AwFmCreate.h"
#include "../../src/AwFmIndex.h"
#include "../../src/AwFmKmerTable.h"
#include "../../src/AwFmLetter.h"
#include "../../src/AwFmOccurrence.h"
#include "../test.h"
#include "divsufsort64.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

char buffer[4096];
const char nucleotideLookup[5]  = {"$ACGT"};
const uint8_t aminoLookup[21]   = {'$','a','c','d','e','f',
                                  'g','h','i','k','l',
                                  'm','n','p','q','r',
                                  's','t','v','w','y'};

void testNucletotideBwtGeneration(void);
void testAminoBwtGeneration(void);

int main (int argc, char **argv){
  srand(8);

  testNucletotideBwtGeneration();
  testAminoBwtGeneration();

}

void testNucletotideBwtGeneration(void){
  printf("beginning nucleotide bwt test\n");
  for(size_t testNum = 0; testNum < 100; testNum++){
    printf("testnum %zu\n", testNum);
    const uint64_t sequenceLength = 4000 + rand()%4000;
    uint8_t *sequence = malloc((sequenceLength+1) * sizeof(uint8_t));
    uint64_t *suffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));

    for(size_t i = 0; i < sequenceLength; i++){
      sequence[i] = nucleotideLookup[rand()%5];
    }

    divsufsort64(sequence, (int64_t*)(suffixArray + 1), sequenceLength);
    suffixArray[0] = sequenceLength;
    sequence[sequenceLength] = 0;

    // printf("sequence:");
    // for(size_t i = 0; i < sequenceLength+1; i++){
    //   if((i) % 10 == 0){
    //     printf("\n");
    //   }
    //   printf("%c, ", sequence[i]);
    // }
    // printf("\n");
    //
    // printf("suffix array:");
    // for(size_t i = 0; i < sequenceLength+1; i++){
    //   if((i) % 10 == 0){
    //     printf("\n");
    //   }
    //   if(suffixArray[i] == 0){
    //     printf("\033[1;31m0\033[0m, ");
    //   }
    //   else{
    //     printf("%zu, ", suffixArray[i]);
    //   }
    // }
    // printf("\n");
    //
    // printf("bwt:");
    // for(size_t i = 0; i < sequenceLength; i++){
    //   if((i) % 10 == 0){
    //     printf("\n");
    //   }
    //
    //   if(suffixArray[i] == 0){
    //     printf("\033[1;31m$\033[0m");
    //   }
    //   else{
    //     printf("%c", sequence[suffixArray[i]-1]);
    //
    //   }
    // }
    // printf("\n");

    struct AwFmIndex *index;
    struct AwFmIndexMetadata metadata = {.versionNumber=1, .suffixArrayCompressionRatio=240, .kmerLengthInSeedTable=2, .alphabetType = AwFmAlphabetNucleotide, .keepSuffixArrayInMemory=false};
    awFmCreateIndex(&index, &metadata, sequence, sequenceLength, "testIndex.awfmi", true);

    for(size_t i = 0; i <= sequenceLength; i++){
      uint64_t suffixArrayValue = suffixArray[i];
      char letterAtBwt = (suffixArrayValue == 0)? '$': sequence[suffixArrayValue-1];
      uint8_t letterIndex = awFmAsciiNucleotideToLetterIndex(letterAtBwt);

      //grab the letter from the AwFmIndex bwt
      size_t blockIndex = i / 256;
      uint8_t byteInBlock = (i/8) % 32;
      uint8_t bitInByte = i % 8;

      uint8_t *letterBitVectorsAsBytes = (uint8_t*)index->bwtBlockList.asNucleotide[blockIndex].letterBitVectors;

      uint8_t letterIndexInAwFm = ((letterBitVectorsAsBytes[byteInBlock]>> bitInByte) & 1) |
                                  ((letterBitVectorsAsBytes[byteInBlock+32]>>bitInByte) & 1) << 1 |
                                  ((letterBitVectorsAsBytes[byteInBlock+64]>>bitInByte) & 1) << 2;
      // printf("index %zu, letter %u, expected %u\n", i, letterIndexInAwFm, letterIndex);
      sprintf(buffer, "at SA position %zu, bwt letter index %u did not match the index stored in AwFm (%d)", i, letterIndex, letterIndexInAwFm);
      testAssertString(letterIndex == letterIndexInAwFm, buffer);
      if(letterIndex != letterIndexInAwFm){
        exit(-04);
      }
    }
    free(sequence);
    free(suffixArray);
    awFmDeallocIndex(index);
  }
}


void testAminoBwtGeneration(void){
  printf("beginning amino bwt test\n");
  for(size_t testNum = 0; testNum < 20; testNum++){
    printf("testnum = %zu\n", testNum);
    const uint64_t sequenceLength = 3000 + rand()%5000;
    uint8_t *sequence = malloc((sequenceLength+1) * sizeof(uint8_t));
    uint64_t *suffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));

    for(size_t i = 0; i < sequenceLength; i++){
      sequence[i] = aminoLookup[rand()%21];
    }

    divsufsort64(sequence, (int64_t*)(suffixArray + 1), sequenceLength);
    suffixArray[0] = sequenceLength;
    sequence[sequenceLength] = 0;

    // printf("sequence:");
    // for(size_t i = 0; i < sequenceLength+1; i++){
    //   if((i) % 10 == 0){
    //     printf("\n");
    //   }
    //   printf("%c, ", sequence[i]);
    // }
    // printf("\n");
    //
    // printf("suffix array:");
    // for(size_t i = 0; i < sequenceLength+1; i++){
    //   if((i) % 10 == 0){
    //     printf("\n");
    //   }
    //   printf("%zu, ", suffixArray[i]);
    // }
    // printf("\n");
    //
    // printf("bwt:");
    // for(size_t i = 0; i < sequenceLength; i++){
    //   if((i) % 10 == 0){
    //     printf("\n");
    //   }
    //
    //   if(suffixArray[i] == 0){
    //     printf("\033[1;31m$\033[0m");
    //   }
    //   else{
    //     printf("%c", sequence[suffixArray[i]-1]);
    //
    //   }
    // }
    // printf("\n");

    struct AwFmIndex *index;
    struct AwFmIndexMetadata metadata = {.versionNumber=1, .suffixArrayCompressionRatio=240, .kmerLengthInSeedTable=4, .alphabetType = AwFmAlphabetAmino, .keepSuffixArrayInMemory=false};
    awFmCreateIndex(&index, &metadata, sequence, sequenceLength, "testIndex.awfmi", true);


    for(size_t i = 0; i <= sequenceLength; i++){
      uint64_t suffixArrayValue = suffixArray[i];
      char letterAtBwt = (suffixArrayValue == 0)? '$': sequence[suffixArrayValue-1];
      uint8_t letterIndex = awFmAsciiAminoAcidToLetterIndex(letterAtBwt);

      //grab the letter from the AwFmIndex bwt
      size_t blockIndex = i / 256;
      uint8_t positionInBlock = i%256;
      // uint8_t byteInBlock = (i/8) % 32;
      // uint8_t bitInByte = i % 8;
      //
      // uint8_t *letterBitVectorsAsBytes = (uint8_t*)index->bwtBlockList.asNucleotide[blockIndex].letterBitVectors;

      uint8_t encodedLetterInAwFm = awFmGetAminoLetterAtBwtPosition(&index->bwtBlockList.asAmino[blockIndex], positionInBlock);
      uint8_t letterIndexInAwFm   = awFmAminoAcidCompressedVectorToLetterIndex(encodedLetterInAwFm);

      sprintf(buffer, "at SA position %zu, bwt letter index %u did not match the index stored in AwFm (%d)", i, letterIndex, letterIndexInAwFm);

      bool lettersMatch = (letterIndex == letterIndexInAwFm) || (encodedLetterInAwFm == 0 && letterAtBwt == '$');
      testAssertString(lettersMatch, buffer);
    }
    free(sequence);
    free(suffixArray);
    awFmDeallocIndex(index);
  }
}
