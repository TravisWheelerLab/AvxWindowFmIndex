#include "../../src/AwFmIndexStruct.h"
#include "../../src/AwFmIndex.h"
#include "../../src/AwFmOccurrence.h"
#include "../../src/AwFmCreate.h"
#include "../../src/AwFmSearch.h"
#include "../../src/AwFmKmerTable.h"
#include "../../src/AwFmLetter.h"
#include "../../src/AwFmParallelSearch.h"
#include "../../src/AwFmBacktraceVector.h"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include "divsufsort64.h"



char buffer[2048];
uint8_t aminoLookup[20]     = {'a','c','d','e','f',
                              'g','h','i','k','l',
                              'm','n','p','q','r',
                              's','t','v','w','y'};
uint8_t nucleotideLookup[4] = {'a','g','c','t'};

void inMemorySaUncompressedTest(void);
void inMemorySaCompressedTest(void);
void inMemoryFromFileTest(void);

int main (int argc, char **argv){
  srand(0);
    inMemorySaUncompressedTest();
    inMemorySaCompressedTest();
  inMemoryFromFileTest();

}

void inMemorySaUncompressedTest(void){
  for(size_t testNum = 0; testNum < 500; testNum++){
    size_t sequenceLength = 200 + rand()%1000;
    uint8_t *sequence = malloc((sequenceLength+100) * sizeof(uint8_t));
    printf("sequence length %zu\n", sequenceLength);
    for(size_t i = 0; i < sequenceLength; i++){
      sequence[i] = nucleotideLookup[rand()%4];
    }
    //add zeros to the end to make sure data is value, but won't compare true to any strings.
    memset(sequence + sequenceLength, 0, 100);

    struct AwFmIndex *index;
    struct AwFmIndexMetadata metadata = {.versionNumber=1, .suffixArrayCompressionRatio = 1,
      .kmerLengthInSeedTable=8, .alphabetType=AwFmAlphabetNucleotide, .keepSuffixArrayInMemory=true};


    enum AwFmReturnCode returnCode = awFmCreateIndex(&index, &metadata, sequence, sequenceLength, "testIndex.awfmi", true);

    if(returnCode < 0){
      printf("error: create returned code %i!!!\n", returnCode);
    }
    size_t numKmers = 100;
    uint8_t kmerLength = 8;
    struct AwFmParallelSearchData *searchData = awFmCreateParallelSearchData(numKmers, 8);
    searchData->count = numKmers;
    for(size_t i = 0; i < numKmers; i++){
      searchData->kmerList[i].length = kmerLength;
      searchData->kmerList[i].string = (char*)&sequence[i];
    }

    awFmParallelSearch(index, searchData);

    for(size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++){
      struct AwFmBacktraceVector *backtraceVector = &searchData->sequencePositionLists[kmerIndex];

      for(size_t sequencePosition = 0; sequencePosition < sequenceLength; sequencePosition++){
        bool kmerFoundAtPosition = strncmp(searchData->kmerList[kmerIndex].string, (char*)&sequence[sequencePosition], kmerLength) == 0;
        bool positionInBacktraceVector = false;
        for(size_t backtraceIndex = 0; backtraceIndex < backtraceVector->count; backtraceIndex++){
          if(backtraceVector->backtraceArray[backtraceIndex].position == sequencePosition){
            positionInBacktraceVector = true;
          }
        }
        sprintf(buffer, "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in backtraceVector? %i.",
          kmerIndex, kmerLength, &sequence[sequencePosition], sequencePosition,
          kmerLength, searchData->kmerList[kmerIndex].string, kmerFoundAtPosition, positionInBacktraceVector);
          if(kmerFoundAtPosition != positionInBacktraceVector){
            printf("not matched: position %zu,  backtrace count %zu, backtrace vector:\n", sequencePosition, backtraceVector->count);
            for(size_t i = 0; i < backtraceVector->count; i++){
              printf("[%zu, (%zu)], ", backtraceVector->backtraceArray[i].position, backtraceVector->backtraceArray[i].offset);

            }
            printf("\n");
          }
        testAssertString(kmerFoundAtPosition == positionInBacktraceVector, buffer);
      }
    }
    free(sequence);
    awFmDeallocParallelSearchData(searchData);
    awFmDeallocIndex(index);
  }
}



void inMemorySaCompressedTest(void){
  for(size_t testNum = 0; testNum < 500; testNum++){
    size_t sequenceLength = 120 + rand()%40;
    printf("sequence length %zu\n", sequenceLength);
    uint8_t *sequence = malloc((sequenceLength+100) * sizeof(uint8_t));
    for(size_t i = 0; i < sequenceLength; i++){
      sequence[i] = nucleotideLookup[rand()%4];
    }
    //add zeros to the end to make sure data is value, but won't compare true to any strings.
    memset(sequence + sequenceLength, 0, 100);

    struct AwFmIndex *index;
    //as an extreme edge case, compress the SA to only one value.
    struct AwFmIndexMetadata metadata = {.versionNumber=1, .suffixArrayCompressionRatio = sequenceLength-1,
      .kmerLengthInSeedTable=8, .alphabetType=AwFmAlphabetNucleotide, .keepSuffixArrayInMemory=true};


    awFmCreateIndex(&index, &metadata, sequence, sequenceLength, "testIndex.awfmi", true);
    size_t numKmers = 100;
    uint8_t kmerLength = 8;
    struct AwFmParallelSearchData *searchData = awFmCreateParallelSearchData(numKmers, 8);
    searchData->count = numKmers;
    for(size_t i = 0; i < numKmers; i++){
      searchData->kmerList[i].length = kmerLength;
      searchData->kmerList[i].string = (char*)&sequence[i];
    }

    awFmParallelSearch(index, searchData);

    for(size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++){
      struct AwFmBacktraceVector *backtraceVector = &searchData->sequencePositionLists[kmerIndex];

      for(size_t sequencePosition = 0; sequencePosition < sequenceLength; sequencePosition++){
        bool kmerFoundAtPosition = strncmp(searchData->kmerList[kmerIndex].string, (char*)&sequence[sequencePosition], kmerLength) == 0;
        bool positionInBacktraceVector = false;
        for(size_t backtraceIndex = 0; backtraceIndex < backtraceVector->count; backtraceIndex++){
          if(backtraceVector->backtraceArray[backtraceIndex].position == sequencePosition){
            positionInBacktraceVector = true;
          }
        }
        sprintf(buffer, "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in backtraceVector? %i.",
          kmerIndex, kmerLength, &sequence[sequencePosition], sequencePosition,
          kmerLength, searchData->kmerList[kmerIndex].string, kmerFoundAtPosition, positionInBacktraceVector);
          if(kmerFoundAtPosition != positionInBacktraceVector){
            printf("not matched: position %zu,  backtrace count %zu, backtrace vector:\n", sequencePosition, backtraceVector->count);
            for(size_t i = 0; i < backtraceVector->count; i++){
              printf("[%zu, (%zu)], ", backtraceVector->backtraceArray[i].position, backtraceVector->backtraceArray[i].offset);

            }
            printf("\n");
          }
        testAssertString(kmerFoundAtPosition == positionInBacktraceVector, buffer);
      }
    }
    free(sequence);
    awFmDeallocParallelSearchData(searchData);
    awFmDeallocIndex(index);
  }
}

void inMemoryFromFileTest(void){
  for(size_t testNum = 0; testNum < 500; testNum++){
    size_t sequenceLength = 160;
    printf("sequence length %zu\n", sequenceLength);
    uint8_t *sequence = malloc((sequenceLength + 100) * sizeof(uint8_t));
    for(size_t i = 0; i < sequenceLength; i++){
      sequence[i] = nucleotideLookup[rand()%4];
    }
    //add zeros to the end to make sure data is value, but won't compare true to any strings.
    memset(sequence + sequenceLength, 0, 100);

    struct AwFmIndex *index;
    struct AwFmIndexMetadata metadata = {.versionNumber=1, .suffixArrayCompressionRatio = 16,
      .kmerLengthInSeedTable=8, .alphabetType=AwFmAlphabetNucleotide, .keepSuffixArrayInMemory=true};


    awFmCreateIndex(&index, &metadata, sequence, sequenceLength, "testIndex.awfmi", true);
    awFmDeallocIndex(index);

    awFmReadIndexFromFile(&index, "testIndex.awfmi", true);

    size_t numKmers = 100;
    uint8_t kmerLength = 8;
    struct AwFmParallelSearchData *searchData = awFmCreateParallelSearchData(numKmers, 8);
    searchData->count = numKmers;
    for(size_t i = 0; i < numKmers; i++){
      searchData->kmerList[i].length = kmerLength;
      searchData->kmerList[i].string = (char*)&sequence[i];
    }

    awFmParallelSearch(index, searchData);

    for(size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++){
      struct AwFmBacktraceVector *backtraceVector = &searchData->sequencePositionLists[kmerIndex];

      for(size_t sequencePosition = 0; sequencePosition < sequenceLength; sequencePosition++){
        bool kmerFoundAtPosition = strncmp(searchData->kmerList[kmerIndex].string, (char*)&sequence[sequencePosition], kmerLength) == 0;
        bool positionInBacktraceVector = false;
        for(size_t backtraceIndex = 0; backtraceIndex < backtraceVector->count; backtraceIndex++){
          if(backtraceVector->backtraceArray[backtraceIndex].position == sequencePosition){
            positionInBacktraceVector = true;
          }
        }
        sprintf(buffer, "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in backtraceVector? %i.",
          kmerIndex, kmerLength, &sequence[sequencePosition], sequencePosition,
          kmerLength, searchData->kmerList[kmerIndex].string, kmerFoundAtPosition, positionInBacktraceVector);
          if(kmerFoundAtPosition != positionInBacktraceVector){
            printf("not matched: position %zu,  backtrace count %zu, backtrace vector:\n", sequencePosition, backtraceVector->count);
            for(size_t i = 0; i < backtraceVector->count; i++){
              printf("[%zu, (%zu)], ", backtraceVector->backtraceArray[i].position, backtraceVector->backtraceArray[i].offset);

            }
            printf("\n");
          }
        testAssertString(kmerFoundAtPosition == positionInBacktraceVector, buffer);
      }
    }
    free(sequence);
    awFmDeallocParallelSearchData(searchData);
    awFmDeallocIndex(index);
  }
}
