#include "AwFmCreate.h"

#include <stdlib.h>
#include <string.h>

#include "AwFmFile.h"
#include "AwFmIndex.h"
#include "AwFmIndexStruct.h"
#include "AwFmLetter.h"
#include "AwFmSearch.h"
#include "AwFmSimdConfig.h"
#include "AwFmSuffixArray.h"
#include "FastaVector.h"
#include "divsufsort64.h"


/*private function prototypes*/
void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t sequenceLength,
		const uint8_t *restrict const sequence, const uint64_t *restrict const unsampledSuffixArray);

void populateKmerSeedTable(struct AwFmIndex *restrict const index);

void populateKmerSeedTableRecursive(struct AwFmIndex *restrict const index, struct AwFmSearchRange range,
		uint8_t currentKmerLength, uint64_t currentKmerIndex, uint64_t letterIndexMultiplier);


void fullSequenceSanitize(const uint8_t *const sequence, uint8_t *const sanitizedSequenceCopy,
		const size_t sequenceLength, const enum AwFmAlphabetType alphabetType);


/*function implementations*/
enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex *restrict *index,
		struct AwFmIndexConfiguration *restrict const config, const uint8_t *restrict const sequence,
		const size_t sequenceLength, const char *restrict const fileSrc, const bool allowFileOverwrite) {

	// first, do a sanity check on inputs
	if(config == NULL) {
		return AwFmNullPtrError;
	}
	if(sequence == NULL) {
		return AwFmNullPtrError;
	}
	if(fileSrc == NULL) {
		return AwFmNullPtrError;
	}

	// set the index out arg initally to NULL, if this function fully completes this will get overwritten
	*index = NULL;

	const size_t suffixArrayLength = sequenceLength + 1;
	// create a sanitized copy of the input sequence
	uint8_t *sanitizedSequenceCopy = malloc(suffixArrayLength);
	if(sanitizedSequenceCopy == NULL) {
		return AwFmAllocationFailure;
	}

	// sanitize the sequence, turning ambiguity characters into the singular ambiguity character
	// character (x for nucleotide, z for amino)
	fullSequenceSanitize(sequence, sanitizedSequenceCopy, sequenceLength, config->alphabetType);

	// append the final sentinel character as a terminator.
	sanitizedSequenceCopy[suffixArrayLength - 1] = '$';

	// allocate the index and all internal arrays.
	struct AwFmIndex *restrict indexData = awFmIndexAlloc(config, suffixArrayLength);
	if(indexData == NULL) {
		return AwFmAllocationFailure;
	}
	indexData->versionNumber = AW_FM_CURRENT_VERSION_NUMBER;
	indexData->featureFlags	 = 0;
	indexData->fastaVector	 = NULL;	// set the fastaVector struct to null, since we aren't using it for this version.
	memcpy(&indexData->config, config, sizeof(struct AwFmIndexConfiguration));

	// init the in memory suffix array to NULL, to be safe. this will get overwritten on success,
	// if the metadata demands in memory SA. If not, this will be left NULL.
	indexData->suffixArray.values = NULL;

	// set the bwtLength
	indexData->bwtLength = suffixArrayLength;

	uint64_t *suffixArray = malloc((suffixArrayLength) * sizeof(uint64_t));

	if(suffixArray == NULL) {
		free(sanitizedSequenceCopy);
		awFmDeallocIndex(indexData);
		return AwFmAllocationFailure;
	}

	// create the suffix array, storing it starting in the second element of the suffix array we allocated.
	// this doesn't clobber the sentinel we added earlier, and makes for easier bwt creation.
	int64_t divSufSortReturnCode = divsufsort64(sanitizedSequenceCopy, (int64_t *)(suffixArray), suffixArrayLength);
	if(divSufSortReturnCode < 0) {
		free(sanitizedSequenceCopy);
		free(suffixArray);
		awFmDeallocIndex(indexData);
		return AwFmSuffixArrayCreationFailure;
	}

	// set the bwt and prefix sums
	setBwtAndPrefixSums(indexData, indexData->bwtLength, sanitizedSequenceCopy, suffixArray);
	// after generating the bwt, the sequence copy is no longer needed.
	free(sanitizedSequenceCopy);

	populateKmerSeedTable(indexData);

	// initialize the compressed suffix array
	enum AwFmReturnCode returnCode = awFmInitCompressedSuffixArray(
			suffixArray, suffixArrayLength, &indexData->suffixArray, config->suffixArrayCompressionRatio);
	indexData->suffixArrayFileOffset = awFmGetSuffixArrayFileOffset(indexData);
	indexData->sequenceFileOffset		 = awFmGetSequenceFileOffset(indexData);
	// file descriptor will be set in awFmWriteIndexToFile

	// create the file
	returnCode = awFmWriteIndexToFile(indexData, sequence, sequenceLength, fileSrc, allowFileOverwrite);

	if(!config->keepSuffixArrayInMemory) {
		free(indexData->suffixArray.values);
		indexData->suffixArray.values = NULL;
	}

	// set the index as an out argument.
	*index = indexData;

	return returnCode;
}

/*function implementations*/
enum AwFmReturnCode awFmCreateIndexFromFasta(struct AwFmIndex *restrict *index,
		struct AwFmIndexConfiguration *restrict const config, const char *fastaSrc, const char *restrict const indexFileSrc,
		const bool allowFileOverwrite) {
			#ifdef _DEBUG
				printf("\x1b[36mDEBUG: awFmCreateIndexFromFasta\x1b[0m\n");
			#endif

	// first, do a sanity check on inputs
	if(config == NULL) {
		return AwFmNullPtrError;
	}
	if(fastaSrc == NULL) {
		return AwFmNullPtrError;
	}
	if(indexFileSrc == NULL) {
		return AwFmNullPtrError;
	}

	// set the index out arg initally to NULL, if this function fully completes this will get overwritten
	*index = NULL;

	// read in the sequence file with FastaVector
	struct FastaVector *fastaVector = malloc(sizeof(struct FastaVector));
	if(fastaVector == NULL) {
		return AwFmAllocationFailure;
	}
	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: fastaVector Allocated\x1b[0m\n");
	#endif
	enum FastaVectorReturnCode fastaVectorReturnCode = fastaVectorInit(fastaVector);
	if(fastaVectorReturnCode == FASTA_VECTOR_ALLOCATION_FAIL) {
		return AwFmAllocationFailure;
	}
	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: fastaVector Init\x1b[0m\n");
	#endif
	fastaVectorReturnCode = fastaVectorReadFasta(fastaSrc, fastaVector);
	if(fastaVectorReturnCode == FASTA_VECTOR_FILE_OPEN_FAIL) {
		return AwFmFileOpenFail;
	}

	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: fasta finished reading\x1b[0m\n");
	#endif
	else if(fastaVectorReturnCode == FASTA_VECTOR_ALLOCATION_FAIL) {
		return AwFmAllocationFailure;
	}

	const size_t fullSequenceLength = fastaVector->sequence.count;
	const char *fullSequencePtr			= fastaVector->sequence.charData;

	const size_t suffixArrayLength = fullSequenceLength + 1;
	// create a sanitized copy of the input sequence
	uint8_t *sanitizedSequenceCopy = malloc(suffixArrayLength);
	if(sanitizedSequenceCopy == NULL) {
		return AwFmAllocationFailure;
	}

	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: sequence copied\x1b[0m"\n);
	#endif

	// sanitize the sequence, turning ambiguity characters into the singular ambiguity character
	// character (x for nucleotide, z for amino)
	// giving the sequencePtr as the in and out arrays makes this an in-place algorithm
	fullSequenceSanitize((uint8_t *)fullSequencePtr, sanitizedSequenceCopy, fullSequenceLength, config->alphabetType);


	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: sequence sanitized\x1b[0m\n");
	#endif
	// FastaVector will always have one extra storage at the end of the string, so this shouldn't ever
	// cause a buffer overflow
	sanitizedSequenceCopy[suffixArrayLength - 1] = '$';

	// allocate the index and all internal arrays.
	struct AwFmIndex *restrict indexData = awFmIndexAlloc(config, suffixArrayLength);
	if(indexData == NULL) {
		return AwFmAllocationFailure;
	}

	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: index allocated\x1b[0m\n");
	#endif

	indexData->versionNumber = AW_FM_CURRENT_VERSION_NUMBER;

	indexData->featureFlags = 0 | (1 << AW_FM_FEATURE_FLAG_BIT_FASTA_VECTOR);
	indexData->fastaVector	= fastaVector;
	memcpy(&indexData->config, config, sizeof(struct AwFmIndexConfiguration));


	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: config allocated\x1b[0m\n");
	#endif
	// init the in memory suffix array to NULL, to be safe. this will get overwritten on success,
	// if the metadata demands in memory SA. If not, this will be left NULL.
	indexData->suffixArray.values = NULL;

	// set the bwtLength
	indexData->bwtLength = suffixArrayLength;

	uint64_t *suffixArray = malloc((suffixArrayLength) * sizeof(uint64_t));

	if(suffixArray == NULL) {
		free(sanitizedSequenceCopy);
		awFmDeallocIndex(indexData);
		return AwFmAllocationFailure;
	}


	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: SA allocated\x1b[0m\n");
	#endif
	// create the suffix array, storing it starting in the second element of the suffix array we allocated.
	// this doesn't clobber the sentinel we added earlier, and makes for easier bwt creation.
	int64_t divSufSortReturnCode = divsufsort64(sanitizedSequenceCopy, (int64_t *)(suffixArray), suffixArrayLength);
	if(divSufSortReturnCode < 0) {
		free(sanitizedSequenceCopy);
		free(suffixArray);
		awFmDeallocIndex(indexData);
		return AwFmSuffixArrayCreationFailure;
	}

	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: SA generated, printing first 10 values:\x1b[0m\n");
		for(size_t i = 0; i < 10; i++){
			printf("\x1b[32mSA index %zu: %zu\x1b[0m\n", i, suffixArray[i]);
		}
	#endif

	// set the bwt and prefix sums
	setBwtAndPrefixSums(indexData, indexData->bwtLength, sanitizedSequenceCopy, suffixArray);
	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: BWT and prefix sums set, printing prefix sums:\x1b[0m\n");
		for(size_t i = 0; i < 6; i++){
			printf("\x1b[32mchars less than symbol %zu index: %zu\x1b[0m\n", i, indexData->prefixSums[i]);
		}
	#endif
	// after generating the bwt, the sequence copy is no longer needed.
	free(sanitizedSequenceCopy);
	#ifdef _DEBUG
		printf("\x1b[36mfreed sequence copy\x1b[0m\n");
	#endif

	populateKmerSeedTable(indexData);
	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: kmer table populated, first 10 values:\x1b[0m\n");
		for(size_t i = 0; i < 10; i++){
			printf("\x1b[32mkmer @ index %zu: [%zu, %zu]\x1b[0m\n", i,
				indexData->kmerSeedTable[i].startPtr, indexData->kmerSeedTable[i].endPtr);
		}
	#endif
	// initialize the compressed suffix array
	enum AwFmReturnCode returnCode = awFmInitCompressedSuffixArray(
			suffixArray, suffixArrayLength, &indexData->suffixArray, config->suffixArrayCompressionRatio);
	if(returnCode != AwFmSuccess) {
		return returnCode;
	}

	#ifdef _DEBUG
		printf("\x1b[36mDEBUG: compressed SA generated, first 10 values:\x1b[0m\n");
		for(size_t i = 0; i < 10; i++){
			printf("\x1b[32m compressed SA index %i: %zu\x1b[0m\n", i, awFmGetValueFromCompressedSuffixArray(indexData->suffixArray, i));
		}
	#endif

	indexData->suffixArrayFileOffset = awFmGetSuffixArrayFileOffset(indexData);
	indexData->sequenceFileOffset		 = awFmGetSequenceFileOffset(indexData);
	// file descriptor will be set in awFmWriteIndexToFile

		#ifdef _DEBUG
			printf("\x1b[36mfile offsets generated\x1b[0m\n");
		#endif

	// create the file
	returnCode =
			awFmWriteIndexToFile(indexData, (uint8_t *)fullSequencePtr, fullSequenceLength, fastaSrc, allowFileOverwrite);

			#ifdef _DEBUG
				printf("\x1b[36mindex written to file\x1b[0m\n");
			#endif
	if(!config->keepSuffixArrayInMemory) {
		// if it's kept in memory, the suffixArray array is now used in the compressedSuffixArray.
		// so, don't free() it now, but it'll need to be free()d at dealloc time.
		free(suffixArray);
		indexData->suffixArray.values = NULL;
		#ifdef _DEBUG
			printf("\x1b[36msuffix array freed (not supposed to be kept in memory at the moment)\x1b[0m\n");
		#endif
	}


	// set the index as an out argument.
	*index = indexData;

	#ifdef _DEBUG
		printf("\x1b[36mgeneration from fasta finished, returning\x1b[0m\n");
	#endif
	return returnCode;
}

void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t bwtLength,
		const uint8_t *restrict const sequence, const uint64_t *restrict const unsampledSuffixArray) {
	if(index->config.alphabetType == AwFmAlphabetNucleotide) {
		uint64_t baseOccurrences[8] = {0};
		// baseOccurrences is length 8 because that's how long the signpost baseOccurrences in
		// each window need to be to keep alignment to 32B AVX2 boundries.

		for(uint64_t suffixArrayPosition = 0; suffixArrayPosition < bwtLength; suffixArrayPosition++) {
			const size_t blockIndex												 = suffixArrayPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
			const uint8_t positionInBlock									 = suffixArrayPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
			const uint8_t byteInVector										 = positionInBlock / 8;
			const uint8_t bitInVectorByte									 = positionInBlock % 8;
			struct AwFmNucleotideBlock *nucleotideBlockPtr = &index->bwtBlockList.asNucleotide[blockIndex];
			uint8_t *restrict const letterBitVectorBytes	 = (uint8_t *)nucleotideBlockPtr->letterBitVectors;

			if(__builtin_expect(positionInBlock == 0, 0)) {
				// when we start a new block, copy over the base occurrences, and initialize the bit vectors
				// while we only use 5 elements, copy over all 8 (to preserve padding and so valgrind
				// doesn't complain about invalid writes)
				memcpy(nucleotideBlockPtr->baseOccurrences, baseOccurrences, 8 * sizeof(uint64_t));
				memset(nucleotideBlockPtr->letterBitVectors, 0, sizeof(AwFmSimdVec256) * AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW);
			}

			uint64_t sequencePositionInSuffixArray = unsampledSuffixArray[suffixArrayPosition];
			uint8_t letterIndex;
			if(__builtin_expect(sequencePositionInSuffixArray == 0, 0)) {
				// set to the sentinel value (index 5) if we're looking at the letter before the first character.
				letterIndex = 5;
			}
			else {
				uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
				letterIndex						 = awFmAsciiNucleotideToLetterIndex(sequence[positionInBwt]);
			}

			uint8_t letterAsCompressedVector = awFmNucleotideLetterIndexToCompressedVector(letterIndex);
			baseOccurrences[letterIndex]++;
			letterBitVectorBytes[byteInVector] |= ((letterAsCompressedVector >> 0) & 0x1) << bitInVectorByte;
			letterBitVectorBytes[byteInVector + 32] |= ((letterAsCompressedVector >> 1) & 0x1) << bitInVectorByte;
			letterBitVectorBytes[byteInVector + 64] |= ((letterAsCompressedVector >> 2) & 0x1) << bitInVectorByte;
		}

		// set the prefix sums
		index->prefixSums[0] = 1;	 // 1 is for the sentinel
		baseOccurrences[0]++;			 // add the sentinel to the count of a's
		for(uint8_t i = 1; i < AW_FM_NUCLEOTIDE_CARDINALITY + 2; i++) {
			index->prefixSums[i] = baseOccurrences[i - 1];
			baseOccurrences[i] += baseOccurrences[i - 1];
		}
	}
	else {
		uint64_t baseOccurrences[24] = {0};

		for(uint64_t suffixArrayPosition = 0; suffixArrayPosition < bwtLength; suffixArrayPosition++) {
			const size_t blockIndex																	= suffixArrayPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
			const uint8_t positionInBlock														= suffixArrayPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
			const uint8_t byteInVector															= positionInBlock / 8;
			const uint8_t bitInVectorByte														= positionInBlock % 8;
			struct AwFmAminoBlock *restrict const aminoBlockPointer = &index->bwtBlockList.asAmino[blockIndex];
			uint8_t *restrict const letterBitVectorBytes						= (uint8_t *)aminoBlockPointer->letterBitVectors;

			if(__builtin_expect(positionInBlock == 0, 0)) {
				// when we start a new block, copy over the base occurrences, and initialize the bit vectors
				// while we only use 21 elements, copy over all 24 (to preserve padding and so valgrind
				// doesn't complain about invalid writes)
				memcpy(aminoBlockPointer->baseOccurrences, baseOccurrences, 24 * sizeof(uint64_t));
				memset(aminoBlockPointer->letterBitVectors, 0, sizeof(AwFmSimdVec256) * AW_FM_AMINO_VECTORS_PER_WINDOW);
			}

			uint64_t sequencePositionInSuffixArray = unsampledSuffixArray[suffixArrayPosition];
			uint8_t letterIndex;
			if(__builtin_expect(sequencePositionInSuffixArray != 0, 1)) {
				uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
				letterIndex						 = awFmAsciiAminoAcidToLetterIndex(sequence[positionInBwt]);
			}
			else {
				// set to the sentinel value (index 21) if we're looking at the letter before the first character.
				letterIndex = 21;
			}
			uint8_t letterAsVectorFormat = awFmAminoAcidLetterIndexToCompressedVector(letterIndex);
			baseOccurrences[letterIndex]++;
			letterBitVectorBytes[byteInVector] |= (((letterAsVectorFormat >> 0) & 0x1) << bitInVectorByte);
			letterBitVectorBytes[byteInVector + 32] |= (((letterAsVectorFormat >> 1) & 0x1) << bitInVectorByte);
			letterBitVectorBytes[byteInVector + 64] |= (((letterAsVectorFormat >> 2) & 0x1) << bitInVectorByte);
			letterBitVectorBytes[byteInVector + 96] |= (((letterAsVectorFormat >> 3) & 0x1) << bitInVectorByte);
			letterBitVectorBytes[byteInVector + 128] |= (((letterAsVectorFormat >> 4) & 0x1) << bitInVectorByte);
		}
		// set the prefix sums
		index->prefixSums[0] = 1;
		baseOccurrences[0]++;	 // add the sentinel into the count of a's
		for(uint8_t i = 1; i < AW_FM_AMINO_CARDINALITY + 2; i++) {
			index->prefixSums[i] = baseOccurrences[i - 1];
			baseOccurrences[i] += baseOccurrences[i - 1];
		}
	}
}

void populateKmerSeedTable(struct AwFmIndex *restrict const index) {
	const uint8_t alphabetCardinality = awFmGetAlphabetCardinality(index->config.alphabetType);
	for(uint8_t i = 0; i < alphabetCardinality; i++) {
		struct AwFmSearchRange range = {.startPtr = index->prefixSums[i], .endPtr = index->prefixSums[i + 1] - 1};
		populateKmerSeedTableRecursive(index, range, 1, i, alphabetCardinality);
	}
}


void populateKmerSeedTableRecursive(struct AwFmIndex *restrict const index, struct AwFmSearchRange range,
		uint8_t currentKmerLength, uint64_t currentKmerIndex, uint64_t letterIndexMultiplier) {
	const uint8_t alphabetSize = awFmGetAlphabetCardinality(index->config.alphabetType);

	const uint8_t kmerLength = index->config.kmerLengthInSeedTable;

	// base case
	if(kmerLength == currentKmerLength) {
		index->kmerSeedTable[currentKmerIndex] = range;
		return;
	}

	// recursive case
	for(uint8_t extendedLetter = 0; extendedLetter < alphabetSize; extendedLetter++) {
		struct AwFmSearchRange newRange = range;
		if(index->config.alphabetType == AwFmAlphabetNucleotide) {
			awFmNucleotideIterativeStepBackwardSearch(index, &newRange, extendedLetter);
		}
		else {
			awFmAminoIterativeStepBackwardSearch(index, &newRange, extendedLetter);
		}

		uint64_t newKmerIndex = currentKmerIndex + (extendedLetter * letterIndexMultiplier);
		populateKmerSeedTableRecursive(
				index, newRange, currentKmerLength + 1, newKmerIndex, letterIndexMultiplier * alphabetSize);
	}
}


inline void fullSequenceSanitize(const uint8_t *const sequence, uint8_t *const sanitizedSequenceCopy,
		const size_t sequenceLength, const enum AwFmAlphabetType alphabetType) {

	if(alphabetType == AwFmAlphabetNucleotide) {
		for(size_t i = 0; i < sequenceLength; i++) {
			sanitizedSequenceCopy[i] = awFmAsciiNucleotideLetterSanitize(sequence[i]);
		}
	}
	else {
		for(size_t i = 0; i < sequenceLength; i++) {
			sanitizedSequenceCopy[i] = awFmAsciiAminoLetterSanitize(sequence[i]);
		}
	}
}
