#include "AwFmIndex.h"
#include "AwFmOccurrence.h"
#include "AwFmSearch.h"


void AwFmIterativeStepBidirectionalNucleotideSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBiDirectionalRange *restrict const range, const uint8_t letter){

  uint64_t letterRankPrefixSum        = AwFmGetRankPrefixSum(index, letter);
  uint64_t sentinelCharacterPosition  = index->sentinelCharacterPosition;

  //query for the start pointer
  //-1 relates to the "Occ(a,s-1) and OccLt(a,s-1) in the literature"
  uint64_t queryPosition = range->startPtr - 1;
  const struct AwFmNucleotideBlock *restrict blockPtr = AwFmGetNucleotideBlockPtr(index, queryPosition);
  struct AwFmOccurrenceVectorPair occurrenceVectors;
  AwFmMakeNucleotideOccurrenceVectorPair(blockPtr, occurrenceVectors, letter, sentinelCharacterPosition);

  uint64_t baseOccurrence     = AwFmGetNucleotideBaseOccurrence(blockPtr, letter);
  uint64_t vectorPopcount     = AwFmNucleotideVectorPopcount(occurrenceVectors.occurrenceVector, queryPosition, sentinelCharacterPosition);
  uint64_t newStartPointer    = vectorPopcount + baseOccurrence + letterRankPrefixSum;
  //like above, the next start query will get data from the position 1 before the startPtr
  AwFmNucleotideBlockPrefetch(newStartPointer - 1);

  uint64_t baseOccurrenceGte  = AwFmGetNucleotideBaseOccurrenceGte(blockPtr, letter);
  uint64_t occurrenceGte      = AwFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t startOccurrenceLt  = queryPosition - occurrenceGte;


  //query for the end pointer
   queryPosition = range->endPtr;
  blockPtr = AwFmGetNucleotideBlockPtr(index, queryPosition);
  AwFmMakeNucleotideOccurrenceVectorPair(blockPtr, occurrenceVectors, letter, sentinelCharacterPosition);

  baseOccurrence = AwFmGetNucleotideBaseOccurrence(blockPtr, letter);
  vectorPopcount = AwFmNucleotideVectorPopcount(occurrenceVectors.occurrenceVector, queryPosition, sentinelCharacterPosition);
  //In the literature, e = C[a] + Occ(a,e) - 1, so the -1 in the next line corresponds to the -1 here.
  uint64_t newEndPointer     = vectorPopcount + baseOccurrence + letterRankPrefixSum - 1;
  AwFmNucleotideBlockPrefetch(newEndPointer);

  baseOccurrenceGte  = AwFmGetNucleotideBaseOccurrenceGte(blockPtr, letter);
  occurrenceGte      = AwFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t endOccurrenceLt       = queryPosition - occurrenceGte;

  uint64_t newStartPrimePtr = range->startPrimePtr  + endOccurrenceLt - startOccurrenceLt;

  range->startPtr       = newStartPointer;
  range->endPtr         = newEndPointer;
  range->startPrimePtr  = newStartPrimePtr;
}


void AwFmIterativeStepBidirectionalAminoAcidSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBiDirectionalRange *restrict const range, const uint8_t letter){

  uint64_t letterRankPrefixSum        = AwFmGetRankPrefixSum(index, letter);

  //query for the start pointer
  //-1 relates to the "Occ(a,s-1) and OccLt(a,s-1) in the literature"
  uint64_t queryPosition = range->startPtr - 1;
  const struct AwFmNucleotideBlock *restrict blockPtr = AwFmGetAminoAcidBlockPtr(index, queryPosition);
  struct AwFmOccurrenceVectorPair occurrenceVectors;
  AwFmMakeAminoAcidOccurrenceVectorPair(blockPtr, occurrenceVectors, letter);

  uint64_t startBaseOccurrence = AwFmGetAminoAcidBaseOccurrence(blockPtr, letter);
  uint64_t startVectorPopcount = AwFmAminoAcidVectorPopcount(occurrenceVectors.occurrenceVector, queryPosition);
  uint64_t newStartPointer     = startVectorPopcount + startBaseOccurrence + letterRankPrefixSum;
  //like above, the next start query will get data from the position 1 before the startPtr
  AwFmAminoAcidBlockPrefetch(newStartPointer - 1);

  uint64_t baseOccurrenceGte  = AwFmGetAminoAcidBaseOccurrenceGte(blockPtr, letter);
  uint64_t occurrenceGte      = AwFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t startOccurrenceLt       = queryPosition - occurrenceGte;


  //query for the end pointer
   queryPosition = range->endPtr;
  blockPtr = AwFmGetAminoAcidBlockPtr(index, queryPosition);
  AwFmMakeAminoAcidOccurrenceVectorPair(blockPtr, occurrenceVectors, letter);

  startBaseOccurrence = AwFmGetAminoAcidBaseOccurrence(blockPtr, letter);
  startVectorPopcount = AwFmVectorPopcount(occurrenceVectors.occurrenceVector, queryPosition);
  //In the literature, e = C[a] + Occ(a,e) - 1, so the -1 in the next line corresponds to the -1 here.
  uint64_t newEndPointer     = startVectorPopcount + startBaseOccurrence + letterRankPrefixSum - 1;
  AwFmAminoAcidBlockPrefetch(newStartPointer);

  baseOccurrenceGte  = AwFmGetAminoAcidBaseOccurrenceGte(blockPtr, letter);
  occurrenceGte      = AwFmAminoAcidVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t endOccurrenceLt = queryPosition - occurrenceGte;

  uint64_t newStartPrimePtr = range->startPrimePtr  + endOccurrenceLt - startOccurrenceLt;

  range->startPtr       = newStartPointer;
  range->endPtr         = newEndPointer;
  range->startPrimePtr  = newStartPrimePtr;
}



void AwFmIterativeStepBackwardNucleotideSearch(const struct AwFmIndex *restrict const index,
struct AwFmBackwardRange *restrict const range, const uint8_t letter){

  uint64_t letterRankPrefixSum        = AwFmGetRankPrefixSum(index, letter);
  uint64_t sentinelCharacterPosition  = index->sentinelCharacterPosition;


  //query for the start pointer
  uint64_t queryPosition = range->startPtr - 1;
  const struct AwFmNucleotideBlock *restrict blockPtr = AwFmGetNucleotideBlockPtr(index, queryPosition);
  uint64_t baseOccurrence = AwFmGetNucleotideBaseOccurrence(blockPtr, letter);
  __m256i occurrenceVector = AwFmGetNucleotideOccurrenceVector(blockPtr, letter);

  uint64_t vectorPopcount = AwFmNucleotideVectorPopcount(occurrenceVector, queryPosition, sentinelCharacterPosition);
  range->startPtr = vectorPopcount + baseOccurrence + letterRankPrefixSum;
  AwFmNucleotideBlockPrefetch(range->startPtr - 1);


  //query for the end pointer
  queryPosition = range->endPtr;
  blockPtr = AwFmGetNucleotideBlockPtr(index, queryPosition);
  baseOccurrence = AwFmGetNucleotideBaseOccurrence(blockPtr, letter);
  occurrenceVector = AwFmGetNucleotideOccurrenceVector(blockPtr, letter);

  vectorPopcount = AwFmNucleotideVectorPopcount(occurrenceVector, queryPosition, sentinelCharacterPosition);
  range->startPtr = vectorPopcount + baseOccurrence + letterRankPrefixSum - 1;
  AwFmNucleotideBlockPrefetch(range->startPtr);
}


void AwFmIterativeStepBackwardAminoAcidSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBackwardRange *restrict const range, const uint8_t letter){

  uint64_t letterRankPrefixSum        = AwFmGetRankPrefixSum(index, letter);

  //query for the start pointer
  uint64_t queryPosition = range->startPtr - 1;
  const struct AwFmNucleotideBlock *restrict blockPtr = AwFmGetAminoAcidBlockPtr(index, queryPosition);
  uint64_t baseOccurrence = AwFmGetAminoAcidBaseOccurrence(blockPtr, letter);
  __m256i occurrenceVector = AwFmGetAminoAcidOccurrenceVector(blockPtr, letter);

  uint64_t vectorPopcount = AwFmAminoAcidVectorPopcount(occurrenceVector, queryPosition);
  range->startPtr = vectorPopcount + baseOccurrence + letterRankPrefixSum;
  AwFmAminoAcidBlockPrefetch(range->startPtr - 1);


  //query for the end pointer
  queryPosition = range->endPtr;
  blockPtr = AwFmGetAminoAcidBlockPtr(index, queryPosition);
  baseOccurrence = AwFmGetAminoAcidBaseOccurrence(blockPtr, letter);
  occurrenceVector = AwFmGetAminoAcidOccurrenceVector(blockPtr, letter);

  vectorPopcount = AwFmAminoAcidVectorPopcount(occurrenceVector, queryPosition);
  range->startPtr = vectorPopcount + baseOccurrence + letterRankPrefixSum - 1;
  AwFmAminoAcidBlockPrefetch(range->startPtr);
}


/*
 * Function:  awFmFindDatabaseHitPositions
 * --------------------
 *  Takes a range of BWT positions, backtraces each position to find the nearest sample in the
 *    compressed suffix array, and looks up those suffix array positions on disk to
 *    determine the corresponding database sequence position for each BWT position
 *    between the searchRange's pointers (inclusive startPtr, exclusive endPtr).
 *
 *  Note: When using a bi-directional FM-index, the range given should correspond the the
 *    traditional backward BWT, not the forward BWT.
 *
 *  It is the caller's responsibility to free() the returned sequence position array.
 *
 *  This function will overwrite the data in the positionArray, returning the
 *    database sequence positions in the corresponding elements of the positionArray.
 *
 *  Note: positionArray and offsetArray should be of equal length. Having either shorter
 *    than positionArrayLength will result in undefined behavior.
 *
 *  Inputs:
 *    index:              Pointer to the valid AwFmIndex struct.
 *    positionArray:          Array of positions in the implied full suffix array to load.
 *      This function will convert these positions to indices in the compressed suffix array,
 *      as long as the positions are multiples of the suffix array compression ratio.
 *    offsetArray:  array of offsets to be added to the database sequence positions.
 *      This is needed because we can only query the suffix array on positions that are sampled.
 *    positionArrayLength:    Length of the positionArray and offsetArray.
 *
 *  Returns:
 *    AwFmReturnCode detailing the result of the read attempt. Possible return values:
 *      AwFmFileReadOkay on success,
 *      AwFmFileOpenFail on failure to open the AwFm Index file
 *      AwFmFileReadFail on failure to read as many characters as was expected by the sequence.
 *      AwFmIllegalPositionError on a suffix array position being out of bounds of
 *        the file's compressed suffix array.
 */
uint64_t *awFmFindDatabaseHitPositions(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const searchRange,
  enum AwFmReturnCode *restrict fileAccessResult){

  const bool isNucleotideDatabase     = index->metadata.alphabetType == AwFmAlphabetTypeNucleotide;
  const uint64_t numPositionsInRange  = awFmSearchRangeLength(searchRange);

  //allocate the position array and offsetArray (uses one malloc call, only needs 1 free call)
  //restrict keywords should be okay, since they shouldn't read/write to overlapping elements
  uint64_t *const restrict positionArray  = malloc(numPositionsInRange * sizeof(uint64_t) * 2);
  uint64_t *const restrict offsetArray    = positionArray + numPositionsInRange;
  //check for allocation failures
  if(__builtin_expect(positionArray == NULL, 0)){
    *fileAccessResult = AwFmAllocationFailure;
    return NULL;
  }

  //call a prefetch for each block that contains the positions that we need to start querying
  if(isNucleotideDatabase){
    for(uint64_t i = searchRange->startPtr; i < searchRange->endPtr; i += POSITIONS_PER_FM_BLOCK){
      AwFmNucleotideBlockPrefetch(index->backwardBwtBlockList.asAmino);
    }
  }
  else{
    for(uint64_t i = searchRange->startPtr; i < searchRange->endPtr; i += POSITIONS_PER_FM_BLOCK){
      AwFmAminoAcidBlockPrefetch(index->backwardBwtBlockList.asAmino);
    }
  }

  //backtrace each position until we have a list of the positions in the database sequence.
  for(uint64_t i=0; i < numPositionsInRange; i++){
    uint64_t databaseSequenceOffset = 0;
    uint64_t backtracePosition      = searchRange->startPtr + i;

    while(!awFmBwtPositionIsSampled(index, backtracePosition)){
      backtracePosition = awFmBackstepBwtPosition(index, backtracePosition);
      databaseSequenceOffset++;
    }

    positionArray[i]  = backtracePosition;
    offsetArray[i]    = databaseSequenceOffset;
  }

  *fileAccessResult = awFmDbSequencePositionsFromSuffixArrayFile(index, positionArray,
    offsetArray, numPositionsInRange);

  return positionArray;
}


/*
 * Function:  awFmDatabaseSingleKmerExactMatch
 * --------------------
 *  Queries the FM-Index for the range of BWT positions that represent instances
 *    of the given Kmer found in the database.
 *
 *  If the given kmer is not found, the AwFmSearch Range will result in a false value
 *   when given to the awFmSearchRangeIsValid function.
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur, including
 *      creating potential segfauts.
 *    kmerLength:   Length of the kmer to be queried. Undefined behavior may occur if
 *      the function is given a kmerLength of 0.
 *
 *  Returns:
 *    AwFmSearch range representing the range of BWT positions where the given
 *    kmer may be found, as long as startPtr < endPtr. Otherwise (startPtr >= endPtr),
 *    the given kmer does not exist in the database sequence.
 */
struct AwFmBackwardRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
const char *restrict const kmer, const uint16_t kmerLength){

  bool isNucleotideDatabase = index->metadata.alphabetType == AwFmAlphabetNucleotide;

    uint_fast16_t kmerQueryPosition = kmerLength - 1;

    //set the initial range from the rankPrefixSums
    const uint8_t firstSuffixLetter = kmer[kmerQueryPosition];
    struct AwFmBackwardRange searchRange = {index->rankPrefixSums[firstSuffixLetter] + 1, index->rankPrefixSums[firstSuffixLetter + 1]};

    //request the data be prefetched for each occurance call.
    if(isNucleotideDatabase){
      awFmNucleotideDataPrefetch(index->backwardBwtBlockList.asAmino, searchRange.startPtr);
      awFmNucleotideDataPrefetch(index->backwardBwtBlockList.asAmino, searchRange.endPtr);

      while(__builtin_expect(awFmSearchRangeIsValid(&searchRange) && kmerQueryPosition, 1)){
        kmerQueryPosition--;
        const uint8_t letter = kmer[kmerQueryPosition];
        AwFmIterativeStepBackwardNucleotideSearch(index, &searchRange,letter);
      }
    }
    else{
      awFmAminoAcidDataPrefetch(index->backwardBwtBlockList.asAmino, searchRange.startPtr);
      awFmAminoAcidDataPrefetch(index->backwardBwtBlockList.asAmino, searchRange.endPtr);

      while(__builtin_expect(awFmSearchRangeIsValid(&searchRange) && kmerQueryPosition, 1)){
        kmerQueryPosition--;
        const uint8_t letter = kmer[kmerQueryPosition];
        AwFmIterativeStepBackwardAminoAcidSearch(index, &searchRange,letter);
      }
    }

  return searchRange;
}


/*
 * Function:  awFmSingleKmerExists
 * --------------------
 *  Queries the FM-Index to determine if the database sequence contains any instances of the given kmer
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur, including
 *      creating potential segfauts.
 *    kmerLength:   Length of the kmer to be queried. Undefined behavior may occur if
 *      the function is given a kmerLength of 0.
 *
 *  Returns:
 *    True if the kmer exists in the database sequence, or false if it
 *      cannot be found in the database sequence.
 */
bool awFmSingleKmerExists(const struct AwFmIndex *restrict const index, const char *restrict const kmer,
  const uint16_t kmerLength){
    struct AwFmBackwardRange kmerRange = awFmDatabaseSingleKmerExactMatch(index, kmer, kmerLength);
    return kmerRange.startPtr < kmerRange.endPtr;
}



/*
 * Function:  awFmSearchRangeLength
 * --------------------
 * Gets the number of positions included in the given AwFmSearchRange
 *
 *  Inputs:
 *    range: Range of positions in the BWT that corresponds to some number of
 *      instances of a given kmer.
 *
 *  Outputs:
 *    Number of positions in the given range if the range is valid (startPtr < endPtr),
 *      or 0 otherwise, as that would imply that no instances of that kmer were found.
 */
 size_t awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range){
  uint64_t length = range->endPtr - range->startPtr;
  return (range->startPtr < range->endPtr)? length: 0;
}
