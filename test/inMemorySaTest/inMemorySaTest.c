#include "../../src/AwFmIndexStruct.h"
#include "../../src/AwFmIndex.h"
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
  srand(time(NULL));
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
    struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numKmers);
    searchList->count = numKmers;
    for(size_t i = 0; i < numKmers; i++){
      searchList->kmerSearchData[i].kmerLength = kmerLength;
      searchList->kmerSearchData[i].kmerString = (char*)&sequence[i];
    }
    awFmParallelSearchLocate(index, searchList, 4);

    for(size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++){
      struct AwFmBacktrace *backtraceList = searchList->kmerSearchData[kmerIndex].positionBacktraceList;

      for(size_t sequencePosition = 0; sequencePosition < sequenceLength; sequencePosition++){
        bool kmerFoundAtPosition = strncmp(searchList->kmerSearchData[kmerIndex].kmerString, (char*)&sequence[sequencePosition], kmerLength) == 0;
        bool positionInBacktraceList = false;
        for(size_t backtraceIndex = 0; backtraceIndex < searchList->kmerSearchData[kmerIndex].count; backtraceIndex++){
          if(backtraceList[backtraceIndex].position == sequencePosition){
            positionInBacktraceList = true;
          }
        }

        sprintf(buffer, "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in backtraceVector? %i.",
          kmerIndex, kmerLength, &sequence[sequencePosition], sequencePosition,
          kmerLength, searchList->kmerSearchData[kmerIndex].kmerString, kmerFoundAtPosition, positionInBacktraceList);
          if(kmerFoundAtPosition != positionInBacktraceList){
            printf("not matched: position %zu,  backtrace count %u, backtrace vector:\n",
              sequencePosition, searchList->kmerSearchData[kmerIndex].count);
            for(size_t i = 0; i < searchList->kmerSearchData[kmerIndex].count; i++){
              printf("[%zu, (%zu)], ", backtraceList[i].position, backtraceList[i]._offset);

            }
            printf("\n");
          }
        testAssertString(kmerFoundAtPosition == positionInBacktraceList, buffer);
      }
    }
    free(sequence);
    awFmDeallocKmerSearchList(searchList);
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
    struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numKmers);
    searchList->count = numKmers;
    for(size_t i = 0; i < numKmers; i++){
      searchList->kmerSearchData[i].kmerLength = kmerLength;
      searchList->kmerSearchData[i].kmerString = (char*)&sequence[i];
    }

    awFmParallelSearchLocate(index, searchList, 4);

    for(size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++){
      struct AwFmBacktrace *backtraceList = searchList->kmerSearchData[kmerIndex].positionBacktraceList;

      for(size_t sequencePosition = 0; sequencePosition < sequenceLength; sequencePosition++){
        bool kmerFoundAtPosition = strncmp(searchList->kmerSearchData[kmerIndex].kmerString, (char*)&sequence[sequencePosition], kmerLength) == 0;
        bool positionInBacktraceList = false;
        for(size_t backtraceIndex = 0; backtraceIndex < searchList->kmerSearchData[kmerIndex].count; backtraceIndex++){
          if(backtraceList[backtraceIndex].position == sequencePosition){
            positionInBacktraceList = true;
          }
        }

        sprintf(buffer, "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in backtraceVector? %i.",
          kmerIndex, kmerLength, &sequence[sequencePosition], sequencePosition,
          kmerLength, searchList->kmerSearchData[kmerIndex].kmerString, kmerFoundAtPosition, positionInBacktraceList);

        if(kmerFoundAtPosition != positionInBacktraceList){
          printf("not matched: position %zu,  backtrace count %u, backtrace vector:\n",
            sequencePosition, searchList->kmerSearchData[kmerIndex].count);
          for(size_t i = 0; i < searchList->kmerSearchData[kmerIndex].count; i++){
            printf("[%zu, (%zu)], ", backtraceList[i].position, backtraceList[i]._offset);
          }
          printf("\n");
        }

        testAssertString(kmerFoundAtPosition == positionInBacktraceList, buffer);
      }
    }
    free(sequence);
    awFmDeallocKmerSearchList(searchList);
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
    struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numKmers);
    searchList->count = numKmers;
    for(size_t i = 0; i < numKmers; i++){
      searchList->kmerSearchData[i].kmerLength = kmerLength;
      searchList->kmerSearchData[i].kmerString = (char*)&sequence[i];
    }

    awFmParallelSearchLocate(index, searchList, 4);

    for(size_t kmerIndex = 0; kmerIndex < numKmers; kmerIndex++){
      struct AwFmBacktrace *backtraceList = searchList->kmerSearchData[kmerIndex].positionBacktraceList;

      for(size_t sequencePosition = 0; sequencePosition < sequenceLength; sequencePosition++){
        bool kmerFoundAtPosition = strncmp(searchList->kmerSearchData[kmerIndex].kmerString, (char*)&sequence[sequencePosition], kmerLength) == 0;
        bool positionInBacktraceList = false;
        for(size_t backtraceIndex = 0; backtraceIndex < searchList->kmerSearchData[kmerIndex].count; backtraceIndex++){
          if(backtraceList[backtraceIndex].position == sequencePosition){
            positionInBacktraceList = true;
          }
        }

        sprintf(buffer, "kmer index %zu  (%.*s) at position %zu  (%.*s) found? %i. in backtraceVector? %i.",
          kmerIndex, kmerLength, &sequence[sequencePosition], sequencePosition,
          kmerLength, searchList->kmerSearchData[kmerIndex].kmerString, kmerFoundAtPosition, positionInBacktraceList);
          if(kmerFoundAtPosition != positionInBacktraceList){
            printf("not matched: position %zu,  backtrace count %u, backtrace vector:\n",
              sequencePosition, searchList->kmerSearchData[kmerIndex].count);
            for(size_t i = 0; i < searchList->kmerSearchData[kmerIndex].count; i++){
              printf("[%zu, (%zu)], ", backtraceList[i].position, backtraceList[i]._offset);

            }
            printf("\n");
          }
        testAssertString(kmerFoundAtPosition == positionInBacktraceList, buffer);
      }
    }
    free(sequence);
    awFmDeallocKmerSearchList(searchList);
    awFmDeallocIndex(index);
  }
}
