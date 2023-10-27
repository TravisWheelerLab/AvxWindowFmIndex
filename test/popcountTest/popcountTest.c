//this test is deprecated.
// #include <getopt.h>
// #include <immintrin.h>
// #include <stdbool.h>
// #include <stdint.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <time.h>

// #include "../../src/AwFmIndexStruct.h"


// // without cache stuff (using same vector), awFmVectorPopcountShift is ~8% faster than mula
// size_t length;
// volatile size_t offset = 0;

// uint8_t extractPopcnt(volatile __m256i vec) {
// 	return _mm_popcnt_u64(_mm256_extract_epi64(vec, 0)) + _mm_popcnt_u64(_mm256_extract_epi64(vec, 1)) +
// 				 _mm_popcnt_u64(_mm256_extract_epi64(vec, 2)) + _mm_popcnt_u64(_mm256_extract_epi64(vec, 3));
// }

// uint8_t storePopcnt(volatile __m256i vec) {
// 	uint8_t unalignedBuffer[32 + 31];
// 	uint8_t *alignedBuffer = (uint8_t *)((intptr_t)(unalignedBuffer + 31) & ~31);
// 	_mm256_store_si256((__m256i *)alignedBuffer, vec);

// 	uint64_t *bufferAs64_t = (uint64_t *)alignedBuffer;
// 	return _mm_popcnt_u64(bufferAs64_t[0]) + _mm_popcnt_u64(bufferAs64_t[1]) + _mm_popcnt_u64(bufferAs64_t[2]) +
// 				 _mm_popcnt_u64(bufferAs64_t[3]);
// }


// inline uint_fast8_t awFmVectorPopcountShift(const __m256i occurrenceVector) {
// 	const __m256i lowBitsLookupTable =
// 			_mm256_setr_epi8(4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0, 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0);
// 	const __m256i highBitsLookupTable = _mm256_setr_epi8(-4, -3, -3, -2, -3, -2, -2, -1, -3, -2, -2, -1, -2, -1, -1, 0,
// 			-4, -3, -3, -2, -3, -2, -2, -1, -3, -2, -2, -1, -2, -1, -1, 0);

// 	const __m256i lowNybbleBitmasked				 = _mm256_and_si256(occurrenceVector, _mm256_set1_epi8(0x0F));
// 	const __m256i lowNybbleBitCount					 = _mm256_shuffle_epi8(lowBitsLookupTable, lowNybbleBitmasked);
// 	const __m256i highNybbleBits						 = _mm256_srli_si256(occurrenceVector, 4);
// 	const __m256i highNybbleBitmasked				 = _mm256_and_si256(highNybbleBits, _mm256_set1_epi8(0x0F));
// 	const __m256i highNybbleNegativeBitCount = _mm256_shuffle_epi8(highBitsLookupTable, highNybbleBitmasked);
// 	const __m256i sadCountVector						 = _mm256_sad_epu8(lowNybbleBitCount, highNybbleNegativeBitCount);
// 	// todo: try keeping nybble counts seperate, shifting both to be in the lower 128 bit regs, and performing separate
// 	// hadds, and extracting out values at end from both nybble vectors.

// 	// TODO: try 4 extracts, see if it's faster
// 	// shift and add, placing the final two 16-bit sums in the least significant bits of each 128-bit lane.
// 	const __m256i laneVectorSums = _mm256_add_epi16(_mm256_slli_si256(sadCountVector, 4), sadCountVector);
// 	const uint16_t finalSum			 = _mm256_extract_epi16(laneVectorSums, 0) + _mm256_extract_epi16(laneVectorSums, 1);
// 	return finalSum;
// }


// inline uint_fast8_t awFmVectorPopcountQuadAdd(const __m256i occurrenceVector) {
// 	const __m256i lowBitsLookupTable =
// 			_mm256_setr_epi8(4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0, 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0);
// 	const __m256i highBitsLookupTable = _mm256_setr_epi8(-4, -3, -3, -2, -3, -2, -2, -1, -3, -2, -2, -1, -2, -1, -1, 0,
// 			-4, -3, -3, -2, -3, -2, -2, -1, -3, -2, -2, -1, -2, -1, -1, 0);

// 	const __m256i lowNybbleBitmasked				 = _mm256_and_si256(occurrenceVector, _mm256_set1_epi8(0x0F));
// 	const __m256i lowNybbleBitCount					 = _mm256_shuffle_epi8(lowBitsLookupTable, lowNybbleBitmasked);
// 	const __m256i highNybbleBits						 = _mm256_srli_si256(occurrenceVector, 4);
// 	const __m256i highNybbleBitmasked				 = _mm256_and_si256(highNybbleBits, _mm256_set1_epi8(0x0F));
// 	const __m256i highNybbleNegativeBitCount = _mm256_shuffle_epi8(highBitsLookupTable, highNybbleBitmasked);
// 	const __m256i sadCountVector						 = _mm256_sad_epu8(lowNybbleBitCount, highNybbleNegativeBitCount);
// 	// todo: try keeping nybble counts seperate, shifting both to be in the lower 128 bit regs, and performing separate
// 	// hadds, and extracting out values at end from both nybble vectors.

// 	// TODO: try 4 extracts, see if it's faster
// 	// shift and add, placing the final two 16-bit sums in the least significant bits of each 128-bit lane.
// 	const uint16_t finalSum = _mm256_extract_epi16(sadCountVector, 0) + _mm256_extract_epi16(sadCountVector, 1) +
// 														_mm256_extract_epi16(sadCountVector, 2) + _mm256_extract_epi16(sadCountVector, 3);
// 	return finalSum;
// }

// inline uint_fast8_t mulaCount(__m256i v) {
// 	__m256i lookup =
// 			_mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
// 	__m256i low_mask	= _mm256_set1_epi8(0x0f);
// 	__m256i lo				= _mm256_and_si256(v, low_mask);
// 	__m256i hi				= _mm256_and_si256(_mm256_srli_epi32(v, 4), low_mask);
// 	__m256i popcnt1		= _mm256_shuffle_epi8(lookup, lo);
// 	__m256i popcnt2		= _mm256_shuffle_epi8(lookup, hi);
// 	__m256i total			= _mm256_add_epi8(popcnt1, popcnt2);
// 	__m256i quadSum		= _mm256_sad_epu8(total, _mm256_setzero_si256());
// 	uint16_t finalSum = _mm256_extract_epi16(quadSum, 0) + _mm256_extract_epi16(quadSum, 1) +
// 											_mm256_extract_epi16(quadSum, 2) + _mm256_extract_epi16(quadSum, 3);
// 	return finalSum;
// }


// void parseArgs(int argc, char **argv) {
// 	int option = 0;
// 	while((option = getopt(argc, argv, "l:o:")) != -1) {
// 		printf("option: %c ", option);
// 		switch(option) {

// 			case 'l':
// 				sscanf(optarg, "%zu", &length);
// 				printf("sequence length: dbSequenceLength\n");
// 				break;
// 			case 'o': sscanf(optarg, "%zu", &offset); break;
// 		}
// 	}
// }


// int main(int argc, char **argv) {
// 	srand(time(NULL));


// 	parseArgs(argc, argv);

// 	if(length == 0) {
// 		printf("error: length is required, and must not be 0.\n");
// 		exit(EXIT_FAILURE);
// 	}

// 	uint8_t *buffer = malloc(1L << 30L);


// 	size_t accumulator = 0;
// 	size_t acc2				 = 0;


// 	clock_t beforeMula = clock();

// 	for(size_t i = 0; i < length; i++) {
// 		// printf("q\n");
// 		uint8_t *vecPtr = buffer + offset;
// 		__m256i vec			= _mm256_loadu_si256((__m256i *)vecPtr);
// 		accumulator += mulaCount(vec);
// 		acc2++;
// 	}
// 	clock_t afterMula		= clock();
// 	clock_t beforeShift = clock();


// 	for(size_t i = 0; i < length; i++) {
// 		// printf("q\n");
// 		uint8_t *vecPtr = buffer + offset;
// 		__m256i vec			= _mm256_loadu_si256((__m256i *)vecPtr);
// 		accumulator += awFmVectorPopcountShift(vec);
// 		acc2++;
// 	}
// 	clock_t afterShift = clock();


// 	clock_t beforeQuad = clock();

// 	printf("length: %zu\n", length);
// 	for(size_t i = 0; i < length; i++) {
// 		// printf("s\n");
// 		uint8_t *vecPtr = buffer + offset;
// 		__m256i vec			= _mm256_loadu_si256((__m256i *)vecPtr);
// 		accumulator += awFmVectorPopcountQuadAdd(vec);
// 		acc2++;
// 	}
// 	clock_t afterQuad = clock();


// 	// store
// 	clock_t beforeStore = clock();

// 	printf("length: %zu\n", length);
// 	uint8_t *vecPtr = buffer + offset;
// 	__m256i vec			= _mm256_loadu_si256((__m256i *)vecPtr);
// 	for(size_t i = 0; i < length; i++) {
// 		// printf("s\n");
// 		accumulator += storePopcnt(vec);
// 		acc2++;
// 	}
// 	clock_t afterStore = clock();

// 	// extract
// 	clock_t beforeExtract = clock();

// 	printf("length: %zu\n", length);
// 	for(size_t i = 0; i < length; i++) {
// 		// printf("s\n");
// 		accumulator += extractPopcnt(vec);
// 		acc2++;
// 	}
// 	clock_t afterExtract = clock();


// 	printf("acc2: %zu\n", acc2);
// 	acc2 = 0;

// 	size_t shiftTicks		= (afterShift - beforeShift);
// 	size_t quadTicks		= (afterQuad - beforeQuad);
// 	size_t mulaTicks		= (afterMula - beforeMula);
// 	size_t extractTicks = (afterExtract - beforeExtract);
// 	size_t storeTicks		= (afterStore - beforeStore);
// 	printf("results: \nShift: \t\t%zu\nQuad: \t\t%zu\nMula: \t\t%zu\nextract:\t%zu\nstore:\t\t%zu\naccu: %zu\n",
// 			shiftTicks, quadTicks, mulaTicks, extractTicks, storeTicks, accumulator);
// }
