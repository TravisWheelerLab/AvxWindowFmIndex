#include "AwFmSearch.h"

#include "AwFmLetter.h"
#include "AwFmOccurrence.h"
#include "AwFmSuffixArray.h"


struct AwFmSearchRange awFmCreateInitialQueryRange(
		const struct AwFmIndex *_RESTRICT_ const index, const char *_RESTRICT_ const query, const uint8_t queryLength) {

	uint8_t finalLetterIndexInQuery;
	if(index->config.alphabetType == AwFmAlphabetNucleotide) {
		finalLetterIndexInQuery = awFmAsciiNucleotideToLetterIndex(query[queryLength - 1]);
	}
	else {
		finalLetterIndexInQuery = awFmAsciiAminoAcidToLetterIndex(query[queryLength - 1]);
	}

	struct AwFmSearchRange searchRange;
	searchRange.startPtr = index->prefixSums[finalLetterIndexInQuery],
	searchRange.endPtr	 = index->prefixSums[finalLetterIndexInQuery + 1] - 1;

	return searchRange;
}


void awFmNucleotideIterativeStepBackwardSearch(
		const struct AwFmIndex *_RESTRICT_ const index, struct AwFmSearchRange *_RESTRICT_ const range, const uint8_t letter) {

	// query for the start pointer
	uint64_t queryPosition				 = range->startPtr - 1;
	const uint64_t letterPrefixSum = index->prefixSums[letter];

	uint64_t blockIndex				 = awFmGetBlockIndexFromGlobalPosition(queryPosition);
	uint8_t localQueryPosition = awFmGetBlockQueryPositionFromGlobalPosition(queryPosition);

	// before needing the bit vectors, we can figure out if they sentinel character will be added.
	uint64_t newStartPointer = letterPrefixSum;
	uint64_t baseOccurrence	 = index->bwtBlockList.asNucleotide[blockIndex].baseOccurrences[letter];
	AwFmSimdVec256 occurrenceVector =
			awFmMakeNucleotideOccurrenceVector(&(index->bwtBlockList.asNucleotide[blockIndex]), letter);

	uint_fast16_t vectorPopcount = AwFmMaskedVectorPopcount(occurrenceVector, localQueryPosition);
	newStartPointer += vectorPopcount + baseOccurrence;

	range->startPtr = newStartPointer;

	// prefetch the next start ptr
	uint64_t newStartBlock = (newStartPointer - 1) / AW_FM_POSITIONS_PER_FM_BLOCK;
	uint8_t *newStartBlockPtr =
			((uint8_t *)index->bwtBlockList.asNucleotide) + (newStartBlock * sizeof(struct AwFmNucleotideBlock));
	AwFmSimdPrefetch(newStartBlockPtr);
	AwFmSimdPrefetch(newStartBlockPtr + 64);

	// query for the new end pointer
	queryPosition			 = range->endPtr;
	blockIndex				 = awFmGetBlockIndexFromGlobalPosition(queryPosition);
	localQueryPosition = awFmGetBlockQueryPositionFromGlobalPosition(queryPosition);

	// the -1 is because of the formula u=Cx[a] + Occ(a,u) -1.
	// we can subtract it here to kill time before needing the block from memory
	uint64_t newEndPointer = letterPrefixSum - 1;

	baseOccurrence	 = index->bwtBlockList.asNucleotide[blockIndex].baseOccurrences[letter];
	occurrenceVector = awFmMakeNucleotideOccurrenceVector(&(index->bwtBlockList.asNucleotide[blockIndex]), letter);
	vectorPopcount	 = AwFmMaskedVectorPopcount(occurrenceVector, localQueryPosition);

	newEndPointer += vectorPopcount + baseOccurrence;

	// prefetch the next start ptr
	uint64_t newEndBlock = (newEndPointer) / AW_FM_POSITIONS_PER_FM_BLOCK;
	uint8_t *newEndBlockPtr =
			((uint8_t *)index->bwtBlockList.asNucleotide) + (newEndBlock * sizeof(struct AwFmNucleotideBlock));
	AwFmSimdPrefetch(newEndBlockPtr);
	AwFmSimdPrefetch(newEndBlockPtr + 64);

	range->endPtr = newEndPointer;
}


void awFmAminoIterativeStepBackwardSearch(
		const struct AwFmIndex *_RESTRICT_ const index, struct AwFmSearchRange *_RESTRICT_ const range, const uint8_t letter) {

	// query for the start pointer
	uint64_t queryPosition				 = range->startPtr - 1;
	const uint64_t letterPrefixSum = index->prefixSums[letter];

	uint64_t blockIndex				 = awFmGetBlockIndexFromGlobalPosition(queryPosition);
	uint8_t localQueryPosition = awFmGetBlockQueryPositionFromGlobalPosition(queryPosition);
	uint64_t baseOccurrence		 = index->bwtBlockList.asAmino[blockIndex].baseOccurrences[letter];
	AwFmSimdVec256 occurrenceVector =
			awFmMakeAminoAcidOccurrenceVector(&(index->bwtBlockList.asAmino[blockIndex]), letter);
	uint16_t vectorPopcount	 = AwFmMaskedVectorPopcount(occurrenceVector, localQueryPosition);
	uint64_t newStartPointer = letterPrefixSum + vectorPopcount + baseOccurrence;

	// prefetch the next start ptr
	uint64_t newStartBlock = (newStartPointer - 1) / AW_FM_POSITIONS_PER_FM_BLOCK;
	uint8_t *newStartBlockPtr =
			((uint8_t *)index->bwtBlockList.asAmino) + (newStartBlock * sizeof(struct AwFmAminoBlock));
	for(size_t cacheLine = 0; cacheLine < 5; cacheLine++) {
		AwFmSimdPrefetch(newStartBlockPtr + (cacheLine * 64));
	}

	range->startPtr = newStartPointer;

	// query for the new end pointer
	queryPosition			 = range->endPtr;
	blockIndex				 = awFmGetBlockIndexFromGlobalPosition(queryPosition);
	localQueryPosition = awFmGetBlockQueryPositionFromGlobalPosition(queryPosition);

	baseOccurrence	 = index->bwtBlockList.asAmino[blockIndex].baseOccurrences[letter];
	occurrenceVector = awFmMakeAminoAcidOccurrenceVector(&(index->bwtBlockList.asAmino[blockIndex]), letter);
	vectorPopcount	 = AwFmMaskedVectorPopcount(occurrenceVector, localQueryPosition);

	const uint64_t newEndPointer = letterPrefixSum + vectorPopcount + baseOccurrence - 1;

	// prefetch the next start ptr
	uint64_t newEndBlock		= (newEndPointer - 1) / AW_FM_POSITIONS_PER_FM_BLOCK;
	uint8_t *newEndBlockPtr = ((uint8_t *)index->bwtBlockList.asAmino) + (newEndBlock * sizeof(struct AwFmAminoBlock));
	for(size_t cacheLine = 0; cacheLine < 5; cacheLine++) {
		AwFmSimdPrefetch(newEndBlockPtr + (cacheLine * 64));
	}

	range->endPtr = newEndPointer;
}


uint64_t *awFmFindDatabaseHitPositions(const struct AwFmIndex *_RESTRICT_ const index,
		const struct AwFmSearchRange *_RESTRICT_ const searchRange, enum AwFmReturnCode *_RESTRICT_ fileAccessResult) {

	const uint64_t numPositionsInRange = awFmSearchRangeLength(searchRange);

	// if there were no elements in the search range, abandon the query.
	if(__builtin_expect(numPositionsInRange == 0, 0)) {
		*fileAccessResult = AwFmGeneralFailure;
		return NULL;
	}

	uint64_t *const _RESTRICT_ positionArray = malloc(numPositionsInRange * sizeof(uint64_t));
	uint64_t *const _RESTRICT_ offsetArray	 = malloc(numPositionsInRange * sizeof(uint64_t));
	// check for allocation failures
	if(__builtin_expect(positionArray == NULL, 0)) {
		*fileAccessResult = AwFmAllocationFailure;
		return NULL;
	}
	if(__builtin_expect(offsetArray == NULL, 0)) {
		free(positionArray);
		*fileAccessResult = AwFmAllocationFailure;
		return NULL;
	}

	// call a prefetch for each block that contains the positions that we need to start querying
	const uint_fast16_t blockWidth = index->config.alphabetType == AwFmAlphabetNucleotide ?
																			 sizeof(struct AwFmNucleotideBlock) :
																			 sizeof(struct AwFmAminoBlock);

	for(uint64_t i = searchRange->startPtr; i < searchRange->endPtr; i += AW_FM_POSITIONS_PER_FM_BLOCK) {
		awFmBlockPrefetch((uint8_t *)index->bwtBlockList.asAmino, blockWidth, i);
	}

	// backtrace each position until we have a list of the positions in the database sequence.
	for(uint64_t i = 0; i < numPositionsInRange; i++) {
		uint64_t databaseSequenceOffset = 0;
		uint64_t backtracePosition			= searchRange->startPtr + i;

		if(index->config.alphabetType == AwFmAlphabetNucleotide) {
			while(!awFmBwtPositionIsSampled(index, backtracePosition)) {
				backtracePosition = awFmNucleotideBacktraceBwtPosition(index, backtracePosition);
				databaseSequenceOffset++;
			}
		}
		else {
			while(!awFmBwtPositionIsSampled(index, backtracePosition)) {
				backtracePosition = awFmAminoBacktraceBwtPosition(index, backtracePosition);
				databaseSequenceOffset++;
			}
		}

		// position is divided by compression ratio to get the index in the suffix array
		positionArray[i] = backtracePosition / index->config.suffixArrayCompressionRatio;
		offsetArray[i]	 = databaseSequenceOffset;
	}

	// get the positions from the suffix array.
	*fileAccessResult = awFmReadPositionsFromSuffixArray(index, positionArray, numPositionsInRange);

	// make sure that reading from the suffix array actually succeeded
	if(*fileAccessResult == AwFmFileReadFail) {
		free(offsetArray);
		return positionArray;
	}

	// add the offsets to the returned positions to get the actual positions of the hits
	for(size_t i = 0; i < numPositionsInRange; i++) {
		positionArray[i] += offsetArray[i];
		positionArray[i] %= index->bwtLength;	 // mod by the length so that the sentinel wraps to zero.
	}
	free(offsetArray);

	*fileAccessResult = AwFmFileReadOkay;
	return positionArray;
}


enum AwFmReturnCode awFmGetLocalSequencePositionFromIndexPosition(const struct AwFmIndex *_RESTRICT_ const index,
		size_t globalPosition, size_t *sequenceNumber, size_t *localSequencePosition) {
	if(!index->fastaVector) {
		return AwFmUnsupportedVersionError;
	}
	else {
		if(globalPosition >= index->bwtLength) {
			return AwFmIllegalPositionError;
		}

		const size_t numSequencesInIndex = index->fastaVector->metadata.count;
		for(size_t sequenceIndex = 0; sequenceIndex < numSequencesInIndex; sequenceIndex++) {
			size_t sequenceEndPosition = index->fastaVector->metadata.data[sequenceIndex].sequenceEndPosition;
			if(sequenceEndPosition > globalPosition) {
				// generate the sequenceNumber and localSequencePosition
				size_t sequenceStartPosition;
				if(sequenceIndex == 0) {
					sequenceStartPosition = 0;
				}
				else {
					sequenceStartPosition = index->fastaVector->metadata.data[sequenceIndex - 1].sequenceEndPosition;
				}
				*sequenceNumber				 = sequenceIndex;
				*localSequencePosition = globalPosition - sequenceStartPosition;
				return AwFmSuccess;
			}
		}

		return AwFmIllegalPositionError;
	}
}


enum AwFmReturnCode awFmGetHeaderStringFromSequenceNumber(
		const struct AwFmIndex *_RESTRICT_ const index, size_t sequenceNumber, char **headerBuffer, size_t *headerLength) {
	if(sequenceNumber > index->fastaVector->metadata.count) {
		return AwFmIllegalPositionError;
	}

	size_t headerStartOffset = 0;
	if(sequenceNumber != 0) {
		headerStartOffset = index->fastaVector->metadata.data[sequenceNumber - 1].headerEndPosition;
	}
	size_t headerEndPosition = index->fastaVector->metadata.data[sequenceNumber].headerEndPosition;

	*headerBuffer = index->fastaVector->header.charData + headerStartOffset;
	*headerLength = headerEndPosition - headerStartOffset;
	return AwFmSuccess;
}


struct AwFmSearchRange awFmDatabaseSingleKmerExactMatch(
		const struct AwFmIndex *_RESTRICT_ const index, const char *_RESTRICT_ const kmer, const uint16_t kmerLength) {
	int8_t kmerLetterPosition = kmerLength - 1;
	uint16_t bwtBlockWidth;
	uint8_t kmerLetterIndex;

	if(index->config.alphabetType == AwFmAlphabetNucleotide) {
		bwtBlockWidth		= sizeof(struct AwFmNucleotideBlock);
		kmerLetterIndex = awFmAsciiNucleotideToLetterIndex(kmer[kmerLetterPosition]);
	}
	else {
		bwtBlockWidth		= sizeof(struct AwFmAminoBlock);
		kmerLetterIndex = awFmAsciiAminoAcidToLetterIndex(kmer[kmerLetterPosition]);
	}

	// create the inital range from the first suffix letter.
	struct AwFmSearchRange range = {
			.startPtr = index->prefixSums[kmerLetterIndex], .endPtr = index->prefixSums[kmerLetterIndex + 1] - 1};

	// start by prefetching the endptr
	awFmBlockPrefetch(index->bwtBlockList.asNucleotide, bwtBlockWidth, range.endPtr);
	if(index->config.alphabetType == AwFmAlphabetNucleotide) {
		while(__builtin_expect(awFmSearchRangeIsValid(&range) && (kmerLetterPosition--), 1)) {
			kmerLetterIndex = awFmAsciiNucleotideToLetterIndex(kmer[kmerLetterPosition]);
			awFmNucleotideIterativeStepBackwardSearch(index, &range, kmerLetterIndex);
		}
	}
	else {
		while(__builtin_expect(awFmSearchRangeIsValid(&range) && (kmerLetterPosition--), 1)) {
			kmerLetterIndex = awFmAsciiAminoAcidToLetterIndex(kmer[kmerLetterPosition]);
			awFmAminoIterativeStepBackwardSearch(index, &range, kmerLetterIndex);
		}
	}
	return range;
}


bool awFmSingleKmerExists(
		const struct AwFmIndex *_RESTRICT_ const index, const char *_RESTRICT_ const kmer, const uint16_t kmerLength) {

	struct AwFmSearchRange kmerRange = awFmDatabaseSingleKmerExactMatch(index, kmer, kmerLength);
	return kmerRange.startPtr <= kmerRange.endPtr;
}


inline size_t awFmNucleotideBacktraceBwtPosition(
		const struct AwFmIndex *_RESTRICT_ const index, const uint64_t bwtPosition) {
	const uint64_t *prefixSums			 = index->prefixSums;
	const uint64_t blockIndex				 = awFmGetBlockIndexFromGlobalPosition(bwtPosition);
	const uint8_t localQueryPosition = awFmGetBlockQueryPositionFromGlobalPosition(bwtPosition);

	const struct AwFmNucleotideBlock *_RESTRICT_ const blockPtr = &index->bwtBlockList.asNucleotide[blockIndex];
	const uint8_t letterIndex = awFmGetNucleotideLetterAtBwtPosition(blockPtr, localQueryPosition);

	// if we encountered the sentinel, we know the position and can stop backtracing
	if(__builtin_expect(letterIndex == 5, 0)) {
		return 0;
	}

	const uint64_t baseOccurrence					= blockPtr->baseOccurrences[letterIndex];
	const AwFmSimdVec256 occurrenceVector = awFmMakeNucleotideOccurrenceVector(blockPtr, letterIndex);


	uint16_t vectorPopcount				= AwFmMaskedVectorPopcount(occurrenceVector, localQueryPosition);
	uint64_t backtraceBwtPosition = prefixSums[letterIndex] + baseOccurrence + vectorPopcount - 1;

	return backtraceBwtPosition;
}


inline size_t awFmAminoBacktraceBwtPosition(const struct AwFmIndex *_RESTRICT_ const index, const uint64_t bwtPosition) {
	const uint64_t *prefixSums			 = index->prefixSums;
	const uint64_t blockIndex				 = awFmGetBlockIndexFromGlobalPosition(bwtPosition);
	const uint8_t localQueryPosition = awFmGetBlockQueryPositionFromGlobalPosition(bwtPosition);

	const struct AwFmAminoBlock *_RESTRICT_ const blockPtr = &index->bwtBlockList.asAmino[blockIndex];
	uint8_t letterIndex																	 = awFmGetAminoLetterAtBwtPosition(blockPtr, localQueryPosition);

	// if we encountered the sentinel, we know the position and can stop backtracing
	if(__builtin_expect(letterIndex == 21, 0)) {
		return 0;
	}

	AwFmSimdVec256 occurrenceVector = awFmMakeAminoAcidOccurrenceVector(blockPtr, letterIndex);
	const uint64_t baseOccurrence		= blockPtr->baseOccurrences[letterIndex];

	const uint16_t vectorPopcount				= AwFmMaskedVectorPopcount(occurrenceVector, localQueryPosition);
	const uint64_t backtraceBwtPosition = prefixSums[letterIndex] + baseOccurrence + vectorPopcount - 1;

	return backtraceBwtPosition;
}


inline void awFmNucleotideNonSeededSearch(const struct AwFmIndex *_RESTRICT_ const index, const char *_RESTRICT_ const kmer,
		const uint8_t kmerLength, struct AwFmSearchRange *range) {

	uint8_t indexInKmerString = kmerLength - 1;
	uint8_t queryLetterIndex	= awFmAsciiNucleotideToLetterIndex(kmer[indexInKmerString]);
	range->startPtr						= index->prefixSums[queryLetterIndex];
	range->endPtr							= index->prefixSums[queryLetterIndex + 1] - 1;

	while(indexInKmerString-- != 0 && awFmSearchRangeIsValid(range)) {
		awFmNucleotideIterativeStepBackwardSearch(index, range, awFmAsciiNucleotideToLetterIndex(kmer[indexInKmerString]));
	}
}


inline void awFmAminoNonSeededSearch(const struct AwFmIndex *_RESTRICT_ const index, const char *_RESTRICT_ const kmer,
		const uint8_t kmerLength, struct AwFmSearchRange *range) {

	uint8_t indexInKmerString = kmerLength - 1;
	uint8_t queryLetterIndex	= awFmAsciiAminoAcidToLetterIndex(kmer[indexInKmerString]);
	range->startPtr = index->prefixSums[queryLetterIndex], range->endPtr = index->prefixSums[queryLetterIndex + 1] - 1;

	while(indexInKmerString-- != 0 && awFmSearchRangeIsValid(range)) {
		awFmAminoIterativeStepBackwardSearch(index, range, awFmAsciiAminoAcidToLetterIndex(kmer[indexInKmerString]));
	}
}
