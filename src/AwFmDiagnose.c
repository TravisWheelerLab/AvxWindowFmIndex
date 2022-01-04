#include <ctype.h>
#include <getopt.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "AwFmIndex.h"
#include "AwFmIndexStruct.h"
#include "AwFmSearch.h"
#include "AwFmParallelSearch.h"
#include "AwFmSuffixArray.h"


/*
Diagnostic tools, useful for debugging the project, or just getting a view into the repo
*/

void displayGeneralIndexData(struct AwFmIndex *index);
void displayDetailedIndexData(struct AwFmIndex *index);
void displayFastaVectorData(struct AwFmIndex *index);
void displayQuerySearchResults(struct AwFmIndex *index, char *query);
void displaySuffixArray(struct AwFmIndex *index);

void displayDiagnosticInfoForIndex(char *indexFileSrc, bool displayGeneral,
  bool dispalyDetailed, bool saInMemory, char *queryString){

	struct AwFmIndex *index;
	enum AwFmReturnCode returnCode = awFmReadIndexFromFile(&index, indexFileSrc, saInMemory);
  if(returnCode < 0){
    printf("error: reading index generated error code %i\n", returnCode);
    exit(3);
  }
  if(displayGeneral){
    displayGeneralIndexData(index);
  }
  if(dispalyDetailed){
    displayDetailedIndexData(index);
    displaySuffixArray(index);
    displayFastaVectorData(index);
  }
  if(queryString != 0){
    printf("Querying for single query %s\n", queryString);
    displayQuerySearchResults(index, queryString);
  }
}

void displayGeneralIndexData(struct AwFmIndex *index){
  printf("Displaying index Configuration: \n");
  printf("\t SA compressionRatio: %i\n", index->config.suffixArrayCompressionRatio);
  printf("\t SA kmer table seed length: %i\n", index->config.kmerLengthInSeedTable);
  printf("\t alphabet: %s\n", index->config.alphabetType == 1? "Amino": "Nucleotide");
  printf("\t SA in memory: %s\n", index->config.keepSuffixArrayInMemory? "True": "False");
  printf("\t storing orig sequence: %s\n", index->config.storeOriginalSequence? "True": "False");
  printf("\n");

  printf("\tIndex VersionNumber: %i\n", index->versionNumber);
  bool isVersionNumberValid = awFmIndexIsVersionValid(index->versionNumber);
  if(!isVersionNumberValid){
    printf("\033[31m ERROR: version number %i is invalid\033[0m\n", index->versionNumber);
  }
  if(index->versionNumber != AW_FM_CURRENT_VERSION_NUMBER){
    printf("\033[33mVersion number %i is not the most recent version (%i).\033[0m\n", index->versionNumber, AW_FM_CURRENT_VERSION_NUMBER);
  }
  printf("FeatureFlags:\n");
  printf("\t FastaVector flag: %s\n", (index->featureFlags & 1)? "True": "False");
  printf("\n");
}

void displayDetailedIndexData(struct AwFmIndex *index){
  printf("BWT length: %zu\n", index->bwtLength);
  //display the prefix sums
  uint8_t prefixSumsLength = awFmGetPrefixSumsLength(index->config.alphabetType);
  printf("prefix sums:\n");
  for(int i = 0; i < prefixSumsLength; i++){
    printf("\n symbol %i: %zu\n", i, index->prefixSums[i]);
  }
  printf("\n");

  printf("First 10 values in the kmer seedTable:\n");
  for(int i = 0; i < 10; i++){
    printf("\t kst index %i: [%zu, %zu]\n", i, index->kmerSeedTable[i].startPtr, index->kmerSeedTable[i].endPtr);
  }
  printf("\n");
}

void displayFastaVectorData(struct AwFmIndex *index){
  if(index->fastaVector == NULL){
    printf("No FastaVector used in this index\n");
  }
  else{
    printf("FastaVector header contents:\n");
    printf("\t count of sequence used to generate index: %zu\n", index->fastaVector->metadata.count);
    printf("\t fastaVectorMetadata capacity: %zu\n", index->fastaVector->metadata.capacity);
    printf("\t fastaVector sequence length: %zu\n", index->fastaVector->sequence.count);
  }
}


void displaySuffixArray(struct AwFmIndex *index){
  size_t saLength = index->bwtLength;
  size_t numValuesToPrint = saLength < 100? saLength: 100;
  printf("Suffix Array, first %zu values:\n", numValuesToPrint);
  for(size_t i = 0; i > numValuesToPrint; i++){
    size_t positionFromSa = awFmGetValueFromCompressedSuffixArray(&index->suffixArray, i);
    printf("\t SA pos %zu: %zu\n", i, positionFromSa);
  }
}

void displayQuerySearchResults(struct AwFmIndex *index, char *query){

  int queryLen = strlen(query);
  struct AwFmSearchRange nonSeededRange;
  if(index->config.alphabetType == AwFmAlphabetNucleotide){
    awFmNucleotideNonSeededSearch(index, query, queryLen, &nonSeededRange);
  }
  else{
    awFmAminoNonSeededSearch(index, query, queryLen, &nonSeededRange);
  }
  printf("Non seeded Kmer Search range: [%zu, %zu], length %zu\n", nonSeededRange.startPtr, nonSeededRange.endPtr, awFmSearchRangeLength(&nonSeededRange));

  struct AwFmSearchRange seededRange = awFmDatabaseSingleKmerExactMatch(index, query, queryLen);
  printf("Seeded Kmer Search range: [%zu, %zu], length %zu\n", seededRange.startPtr, seededRange.endPtr,
    awFmSearchRangeLength(&seededRange));

  struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(1);
    if(searchList == NULL){
      printf("Error: could not allocate memory for the search data struct.\n");
      exit(-2);
    }
  searchList->kmerSearchData[0].kmerString = query;
  searchList->kmerSearchData[0].kmerLength = queryLen;
  searchList->count = 1;
  int numThreadsInParallelQuery = 1;
  awFmParallelSearchCount(index, searchList, numThreadsInParallelQuery);
  printf("parallel search count returned %u hits\n", searchList->kmerSearchData[0].count);


  awFmParallelSearchLocate(index, searchList, numThreadsInParallelQuery);
  printf("parallel search locate returned %u hits\n", searchList->kmerSearchData[0].count);

  size_t numToPrint = 0;
  if(searchList->kmerSearchData[0].count > 100){
    printf("query returned lots of matches, printing to 100.\n");
    numToPrint = 100;
  }
  else{
    printf("printing all %u matches\n", searchList->kmerSearchData[0].count);
    numToPrint = searchList->kmerSearchData[0].count;
  }

  for(size_t i = 0; i < numToPrint; i++){
    printf("\t hit %zu: %zu\n", numToPrint, searchList->kmerSearchData[0].positionList[i]);
  }



}
