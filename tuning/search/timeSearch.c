#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <time.h>
#include "AwFmIndex.h"

uint8_t aminoLookup[20]     = {'a','c','d','e','f',
                              'g','h','i','k','l',
                              'm','n','p','q','r',
                              's','t','v','w','y'};
uint8_t nucleotideLookup[4] = {'a','g','c','t'};

char indexFilenameBuffer[256];
int numKmersToQuery;
int kmerLength;
int numThreadsInParallelQuery;
int randSeed;
bool randSeedSetInArgs;
bool keepSuffixArrayInMemory;
bool useCountFunction;

void parseArgs(int argc, char **argv);

int main(int argc, char **argv){
  //set default parameters
  strcpy(indexFilenameBuffer, "index.awfmi");
  numKmersToQuery = 10000;
  kmerLength = 12;
  randSeedSetInArgs = false;
  keepSuffixArrayInMemory = false;
  useCountFunction = false;

  parseArgs(argc, argv);

  //set the random seed
  if(!randSeedSetInArgs){
    randSeed = time(NULL);
  }
  srand(randSeed);

  struct AwFmIndex *index;
//  printf("reading AwFmIndex file...");
  enum AwFmReturnCode returnCode = awFmReadIndexFromFile(&index, indexFilenameBuffer, keepSuffixArrayInMemory);
//  printf(" index file read.\n");
  if(returnCode < 0){
    printf("Error during index read: awFmReadIndexFromFile returned error code %i\n", returnCode);
    exit(-1);
  }

  struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numKmersToQuery);
  searchList->count = numKmersToQuery;
 //   printf(" search data struct generated.\n");

  if(searchList == NULL){
    printf("Error: could not allocate memory for the search data struct.\n");
    exit(-2);
  }

  char *kmerQuerySequence = malloc(numKmersToQuery + kmerLength * sizeof(char));
  if(kmerQuerySequence == NULL){
    printf("error: could not allocate memory for query kmers.\n");
    exit(-3);
  }

  //randomize the kmer query sequence
  for(size_t i = 0 ; i < numKmersToQuery + kmerLength; i++){
    kmerQuerySequence[i] = index->metadata.alphabetType == AwFmAlphabetNucleotide?
      nucleotideLookup[rand()%4]:
      aminoLookup[rand()%20];
  }

  //set the kmers in the search data.
  for(size_t i = 0; i < numKmersToQuery;i++){
    searchList->kmerSearchData[i].kmerLength = kmerLength;
    searchList->kmerSearchData[i].kmerString = kmerQuerySequence+i;
  }

  //printf("beginning parallel search...");
  clock_t searchStartTime = clock();
  //search for the kmers
  if(useCountFunction){
    awFmParallelSearchCount(index, searchList, numThreadsInParallelQuery);
  }
  else{
    awFmParallelSearchLocate(index, searchList, numThreadsInParallelQuery);
  }
  clock_t searchEndTime = clock();
  clock_t elapsedSearchTime = searchEndTime - searchStartTime;

  //printf("parallel search complete.\n");
  //printf("elapsed search time for %i kmers: %zu ticks (%f) seconds\n",
//  numKmersToQuery, elapsedSearchTime, (float)elapsedSearchTime/ (float)CLOCKS_PER_SEC);
//	printf("start time %zu, end time %zu\n", searchStartTime, searchEndTime);
  //dealloc the searchData struct
  awFmDeallocKmerSearchList(searchList);
  awFmDeallocIndex(index);
  free(kmerQuerySequence);

}


void parseArgs(int argc, char **argv){
  int option = 0;
  while((option = getopt(argc, argv, "mcf:n:k:t:s:")) != -1){
    //printf("option: %c ", option);
    switch(option){
      case 'm':
        keepSuffixArrayInMemory = true;
        break;
      case 'c':
        useCountFunction = true;
        break;
      case 'f':
        strcpy(indexFilenameBuffer, optarg);
        break;
      case 'n':
        sscanf(optarg, "%i", &numKmersToQuery);
        break;
      case 'k':
        sscanf(optarg, "%i", &kmerLength);
        break;
      case 't':
        sscanf(optarg, "%i", &numThreadsInParallelQuery);
        break;
      case 's':
        sscanf(optarg, "%i", &randSeed);
        randSeedSetInArgs = true;
        break;
    }
  }
}
