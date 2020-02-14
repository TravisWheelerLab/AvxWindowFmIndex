#include "../../AwFmIndex.h"
#include "../../AwFmCreate.h"
#include "../../AwFmKmerTable.h"
#include "../../AwFmLetter.h"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <divsufsort64.h>

char buffer[2048];
uint8_t aminoLookup[20]     = {'a','c','d','e','f',
                              'g','h','i','k','l',
                              'm','n','p','q','r',
                              's','t','v','w','y'};
uint8_t nucleotideLookup[4] = {'a','g','c','t'};

void checkRangeForCorrectness(const struct AwFmSearchRange *restrict const range,
  const uint64_t *restrict const suffixArray, const char *kmer, const uint8_t kmerLength,
  const uint8_t *restrict const sequence, const uint64_t sequenceLength);

bool rangeCompare(struct AwFmSearchRange range1, struct AwFmSearchRange range2);

void testPartialKmerRanges(const struct AwFmIndexMetadata *metadata);


int main(int argc, char **argv){
  srand(3);
  printf("main\n");

  struct AwFmIndexMetadata metadata = {.versionNumber=1, .suffixArrayCompressionRatio = 240, .kmerLengthInSeedTable = 8, .alphabetType =AwFmAlphabetNucleotide};
  testPartialKmerRanges(&metadata);
  // metadata.alphabetType = AwFmAlphabetAmino;
  // testPartialKmerRanges(&metadata);

  printf("end\n");
}


void checkRangeForCorrectness(const struct AwFmSearchRange *restrict const range,
  const uint64_t *restrict const suffixArray, const char *kmer, const uint8_t kmerLength,
  const uint8_t *restrict const sequence, const uint64_t sequenceLength){

    //make a null terminated copy of the kmer
  char kmerBuffer[kmerLength+1];
  memcpy(kmerBuffer, kmer, kmerLength);
  kmerBuffer[kmerLength] = 0;

  uint64_t lengthOfRange = 1 + range->endPtr - range->startPtr;

  for(uint64_t position = 0; position <= sequenceLength - kmerLength; position++){
    bool kmerFoundAtPosition = strncmp(kmer, (char*)sequence + position, kmerLength) == 0;
    bool positionFoundInRange = false;

    for(uint64_t suffixArrayPosition = range->startPtr; suffixArrayPosition <= range->endPtr; suffixArrayPosition++){
      // printf("SAP %zu, POS %zu\n", suffixArray[suffixArrayPosition], position);
      positionFoundInRange |= (suffixArray[suffixArrayPosition] == position);
    }

    if(!positionFoundInRange){
      char sequenceBuffer[kmerLength+1];
      memcpy(sequenceBuffer, sequence+position, kmerLength);
      sequenceBuffer[kmerLength] = 0;

      sprintf(buffer, "range had position %zu for kmer %s, but that position holds kmer %s.\n", position, kmerBuffer, sequenceBuffer);
    }
    else{
      sprintf(buffer, "position %zu matches kmer %s, but was not represented in the range.\n", position, kmerBuffer);
    }
    testAssertString( (kmerFoundAtPosition == positionFoundInRange), buffer);
    if(kmerFoundAtPosition != positionFoundInRange){
      printf("failure: position: %zu. kmer found at position %u, found in range %u\n", position, kmerFoundAtPosition, positionFoundInRange);
      printf("kmer @pos: %.*s\n", 5, &sequence[position]);
      printf("kmer: %s, ", kmerBuffer);
      printf("seq len: %zu, range [%zu, %zu]\n", sequenceLength, range->startPtr, range->endPtr);
      // printf("range values: ");
      printf("sa positions: ");
      for(uint64_t i = range->startPtr; i <= range->endPtr; i++){
        printf("%zu, ",suffixArray[i]);
      }
      printf("\n");
      // for(uint64_t ptr = range->startPtr; ptr <= range->endPtr;ptr++){
      //   printf(" %zu, ", suffixArray[ptr]);
      //
      // }
      printf("\n");
      exit(-3);
    }
  }
}

bool rangeCompare(struct AwFmSearchRange range1, struct AwFmSearchRange range2){
  return (range1.startPtr == range2.startPtr) && (range1.endPtr == range2.endPtr);
}



void testPartialKmerRanges(const struct AwFmIndexMetadata *metadata){
  struct AwFmIndex *index;
  for(uint8_t testNum = 0; testNum < 10; testNum++){
    const uint64_t sequenceLength = (rand() % 6000) + 100;
    printf("sequence length: %zu\n", sequenceLength);
    uint8_t *sequence = malloc((sequenceLength + 1) * sizeof(uint8_t));
    if(sequence == NULL){
      printf("CRITICAL FAILURE: could not allocate sequence\n");
      exit(-1);
    }
    for(uint64_t i = 0; i < sequenceLength; i++){
      sequence[i] = (metadata->alphabetType == AwFmAlphabetNucleotide?
        nucleotideLookup[rand()%4]: aminoLookup[rand()%20]);
    }
    //null terminate the sequence
    sequence[sequenceLength] = 0;

    printf("sequence: %s\n", sequence);
    //create the reference suffix array
    uint64_t *suffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
    if(suffixArray == NULL){
      printf("CRITICAL FAILURE: could not allocate suffix array\n");
      exit(-1);
    }

    //set the suffix array. remember, the first position isn't included by libdivsufsort,
    //and needs to be manually set to the sentinel character position.
    suffixArray[0] = sequenceLength;
    divsufsort64(sequence, (int64_t*)(suffixArray + 1), sequenceLength);

    awFmCreateIndex(&index, metadata, sequence, sequenceLength, "testIndex.awfmi", true);


    for(uint8_t kmerTestNum = 0; kmerTestNum < 100; kmerTestNum++){
      for(uint8_t kmerLength = 4; kmerLength < index->metadata.kmerLengthInSeedTable; kmerLength++){
        char kmer[kmerLength+1];
        for(uint8_t letterIndex = 0; letterIndex < kmerLength; letterIndex++){
          kmer[letterIndex]  = metadata->alphabetType == AwFmAlphabetNucleotide?
            nucleotideLookup[rand()%4]: aminoLookup[rand()%20];
        }
        kmer[kmerLength] = 0;
        printf("checking kmer %s...\n", kmer);

        struct AwFmSearchRange kmerRange = metadata->alphabetType == AwFmAlphabetNucleotide?
          awFmNucleotideKmerSeedRangeFromTable(index, kmer, kmerLength):
          awFmAminoKmerSeedRangeFromTable(index, kmer, kmerLength);

        checkRangeForCorrectness(&kmerRange, suffixArray, kmer, kmerLength, sequence, sequenceLength);
      }
    }

    awFmDeallocIndex(index);
    free(sequence);
    free(suffixArray);
  }
}
