#include <ctype.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "AwFmIndex.h"

uint8_t aminoLookup[20] = {
		'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'};
uint8_t nucleotideLookup[4] = {'a', 'g', 'c', 't'};

char indexFilenameBuffer[1024];
char sequenceFilenameBuffer[1024];
int numKmersToQuery;
int kmerLength;
int numThreadsInParallelQuery;
bool randSeedSetInArgs;
bool keepSuffixArrayInMemory;
bool useCountFunction;

void parseArgs(int argc, char **argv);
void allocateQueryBuffers(struct AwFmKmerSearchList *searchList, size_t numKmers, uint8_t kmerLength);
void fillSearchListwithKmers(
		struct AwFmKmerSearchList *searchList, size_t numKmers, uint8_t kmerLength, struct AwFmIndex *index);
bool charIsAmbiguous(char c, enum AwFmAlphabetType alphabet);
bool getStringFromSequenceFile(
		FILE *openSequenceFile, size_t position, char *buffer, int kmerLength, enum AwFmAlphabetType alphabet);
void makeValidQueryFromSequenceFile(FILE *openSequenceFile, size_t sequenceLength, uint8_t kmerLength,
		enum AwFmAlphabetType alphabet, char *queryBuffer);


int main(int argc, char **argv) {
	// set default parameters
	strcpy(indexFilenameBuffer, "index.awfmi");
	numKmersToQuery					= 10000;
	kmerLength							= 12;
	randSeedSetInArgs				= false;
	keepSuffixArrayInMemory = false;
	useCountFunction				= false;

	parseArgs(argc, argv);


	struct AwFmIndex *index;
	enum AwFmReturnCode returnCode = awFmReadIndexFromFile(&index, indexFilenameBuffer, keepSuffixArrayInMemory);
	if(returnCode < 0) {
		printf("Error during index read: awFmReadIndexFromFile returned error code %i\n", returnCode);
		exit(-1);
	}

	struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numKmersToQuery);
	if(searchList == NULL) {
		printf("Error: could not allocate memory for the search data struct.\n");
		exit(-2);
	}
	allocateQueryBuffers(searchList, numKmersToQuery, kmerLength);


	clock_t totalTime = 0;
	// run the test 4 times, averaging the result.
	for(size_t i = 0; i < 4; i++) {
		fillSearchListwithKmers(searchList, numKmersToQuery, kmerLength, index);

		// printf("beginning parallel search...");
		clock_t searchStartTime = clock();
		// search for the kmers
		if(useCountFunction) {
			awFmParallelSearchCount(index, searchList, numThreadsInParallelQuery);
		}
		else {
			awFmParallelSearchLocate(index, searchList, numThreadsInParallelQuery);
		}
		clock_t searchEndTime			= clock();
		clock_t elapsedSearchTime = searchEndTime - searchStartTime;

		totalTime += elapsedSearchTime;
	}

	totalTime /= 4;
	float totalTimeAsFloat = (float)totalTime / (float)CLOCKS_PER_SEC;
	printf("average:%f\n", totalTimeAsFloat);

	// deallocate the query strings.
	for(size_t i = 0; i < numKmersToQuery; i++) {
		free(searchList->kmerSearchData[i].kmerString);
	}
	awFmDeallocKmerSearchList(searchList);
	awFmDeallocIndex(index);
}


void parseArgs(int argc, char **argv) {
	int option = 0;
	while((option = getopt(argc, argv, "mcf:j:n:k:t:")) != -1) {
		// printf("option: %c ", option);
		switch(option) {
			case 'm': keepSuffixArrayInMemory = true; break;
			case 'c': useCountFunction = true; break;
			case 'f': strcpy(indexFilenameBuffer, optarg); break;
			case 'j': strcpy(sequenceFilenameBuffer, optarg); break;
			case 'n': sscanf(optarg, "%i", &numKmersToQuery); break;
			case 'k': sscanf(optarg, "%i", &kmerLength); break;
			case 't': sscanf(optarg, "%i", &numThreadsInParallelQuery); break;
		}
	}
}

void fillSearchListwithKmers(
		struct AwFmKmerSearchList *searchList, size_t numKmers, uint8_t kmerLength, struct AwFmIndex *index) {
	FILE *sequenceFile = fopen(sequenceFilenameBuffer, "r");

	for(size_t i = 0; i < numKmers; i++) {
		searchList->kmerSearchData[i].kmerLength = kmerLength;

		makeValidQueryFromSequenceFile(sequenceFile, index->bwtLength - 1, kmerLength, index->metadata.alphabetType,
				searchList->kmerSearchData[i].kmerString);
	}
}


bool charIsAmbiguous(char c, enum AwFmAlphabetType alphabet) {
	if(alphabet == AwFmAlphabetNucleotide) {
		for(uint8_t i = 0; i < 4; i++) {
			if(tolower(c) == nucleotideLookup[i]) {
				return false;
			}
		}
	}
	else {
		for(uint8_t i = 0; i < 20; i++) {
			if(tolower(c) == aminoLookup[i]) {
				return false;
			}
		}
	}
	return true;
}

void allocateQueryBuffers(struct AwFmKmerSearchList *searchList, size_t numKmers, uint8_t kmerLength) {
	searchList->count = numKmers;
	for(size_t i = 0; i < numKmers; i++) {
		searchList->kmerSearchData[i].kmerString = calloc(kmerLength + 1, sizeof(char));
	}
}

void makeValidQueryFromSequenceFile(FILE *openSequenceFile, size_t sequenceLength, uint8_t kmerLength,
		enum AwFmAlphabetType alphabet, char *queryBuffer) {
	size_t position = rand() % (sequenceLength - 100) + 100;	// extra math to skip the first header.

	bool hasValidQuery = false;
	while(!hasValidQuery) {
		hasValidQuery = getStringFromSequenceFile(openSequenceFile, position, queryBuffer, kmerLength, alphabet);
	}
}


bool getStringFromSequenceFile(
		FILE *openSequenceFile, size_t position, char *buffer, int kmerLength, enum AwFmAlphabetType alphabet) {
	fseek(openSequenceFile, position, SEEK_SET);
	int charactersGrabbed = 0;
	while(charactersGrabbed < kmerLength && !feof(openSequenceFile)) {
		char charFromFile = fgetc(openSequenceFile);

		// if the character is at least a character...
		if(charFromFile >= 'A' && charFromFile < 'z') {
			if(charIsAmbiguous(charFromFile, alphabet)) {
				return false;
			}
			else {
				buffer[charactersGrabbed++] = charFromFile;
			}
		}
	}

	return true;
}
