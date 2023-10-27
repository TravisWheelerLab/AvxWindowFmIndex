//this test is deprecated
// #include <getopt.h>
// #include <stdbool.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <time.h>

// #include "../../src/AwFmIndex.h"
// #include "../../src/AwFmIndexStruct.h"
// #include "../../src/AwFmOccurrence.h"


// uint32_t numQueries							 = 0;
// uint32_t dbSizeInWindows				 = 0;
// bool requestOccupancyPrefetching = false;


// void parseArgs(int argc, char **argv) {
// 	int option = 0;
// 	while((option = getopt(argc, argv, "pq:s:")) != -1) {
// 		printf("option: %c ", option);
// 		switch(option) {
// 			case 'p':
// 				printf("prefetching requested\n");
// 				requestOccupancyPrefetching = true;
// 				break;
// 			case 'q':
// 				numQueries = atoi(optarg);
// 				printf("num queries: %d\n", numQueries);
// 				break;
// 			case 's':
// 				dbSizeInWindows = atoi(optarg);
// 				printf("db size in windows: %d\n", dbSizeInWindows);
// 				break;
// 			case '?': printf("unknown option flag given: %c is not a supported flag.\n", optopt); break;
// 		}
// 	}
// }

// void checkArgs() {
// 	bool allArgumentsParsed = true;
// 	if(numQueries == 0) {
// 		printf(
// 				"Num queries was not given (or was given 0). please provide the number of queries to perform with the -q "
// 				"flag\n");
// 		allArgumentsParsed = false;
// 	}
// 	if(dbSizeInWindows == 0) {
// 		printf(
// 				"database size was not given (or was given 0). please provide the size of the db to test with the -s flag\n");
// 		allArgumentsParsed = false;
// 	}
// 	if(!allArgumentsParsed) {
// 		exit(-1);
// 	}
// }

// void performDbQueries(const struct AwFmIndex *_RESTRICT_ const index, uint64_t positionsInDb) {
// 	size_t ptrs[2];
// 	ptrs[0] = rand() % positionsInDb;
// 	ptrs[1] = rand() % positionsInDb;


// 	if(requestOccupancyPrefetching) {
// 		for(uint32_t i = 0; i < numQueries; i++) {
// 			uint8_t letter = rand() % 20;
// 			ptrs[0]				 = awFmGetOccupancy(index, ptrs[0], letter) + index->rankPrefixSums[letter];
// 			ptrs[0]				 = ptrs[0] % positionsInDb;
// 			awFmOccupancyDataPrefetch(index, ptrs[0]);

// 			ptrs[1] = awFmGetOccupancy(index, ptrs[1], letter) + index->rankPrefixSums[letter];
// 			ptrs[1] = ptrs[1] % positionsInDb;
// 			awFmOccupancyDataPrefetch(index, ptrs[1]);
// 		}
// 	}
// 	else {
// 		for(uint32_t i = 0; i < numQueries; i++) {
// 			uint8_t letter = rand() % 20;
// 			ptrs[0]				 = awFmGetOccupancy(index, ptrs[0], letter) + index->rankPrefixSums[letter];
// 			ptrs[0]				 = ptrs[0] % positionsInDb;

// 			ptrs[1] = awFmGetOccupancy(index, ptrs[1], letter) + index->rankPrefixSums[letter];
// 			ptrs[1] = ptrs[1] % positionsInDb;
// 		}
// 	}
// }

// int main(int argc, char **argv) {
// 	srand(time(NULL));
// 	parseArgs(argc, argv);
// 	checkArgs();
// 	uint64_t positionsInDb = dbSizeInWindows * POSITIONS_PER_FM_BLOCK;
// 	// create the "database"
// 	struct AwFmBlock *blockList = awFmAlignedAllocBlockList(dbSizeInWindows);
// 	if(blockList == NULL) {
// 		printf("could not allocate block list... maybe window count of %d was too big?\n", dbSizeInWindows);
// 		exit(-2);
// 	}

// 	// randomize the data in the blocklist
// 	for(size_t i = 0; i < (sizeof(struct AwFmBlock) * dbSizeInWindows) / sizeof(uint32_t); i++) {
// 		((uint32_t *)blockList)[i] = rand();
// 	}

// 	struct AwFmIndex *index = awFmAlignedAllocAwFmIndex();
// 	if(index == NULL) {
// 		printf("could not allocate index... wtf?\n");
// 		exit(-2);
// 	}
// 	// set the block list in or fake index.
// 	index->blockList = blockList;
// 	for(size_t i = 0; i <= AMINO_CARDINALITY; i++) {
// 		size_t prefixSum														 = (positionsInDb / AMINO_CARDINALITY) * i;
// 		index->rankPrefixSums[AMINO_CARDINALITY - i] = prefixSum;
// 	}

// 	clock_t before = clock();
// 	performDbQueries(index, positionsInDb);
// 	clock_t timeElapsed = clock() - before;

// 	double seconds = (double)timeElapsed / CLOCKS_PER_SEC;
// 	printf("completed test in %f seconds. queries: %d, numWindows: %d\n\n\n", seconds, numQueries, dbSizeInWindows);
// 	free(index->blockList);
// 	free(index);
// }
