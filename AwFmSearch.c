#include "AwFmIndex.h"
#include "AwFmOccurrence.h"
#include "AwFmSearch.h"


void awFmIterativeStepBiDirectionalSearch(const struct AwFmIndex *restrict const index,
  const enum AwFmSearchDirection searchDirection, struct AwFmBiDirectionalRange *restrict const range,
  const uint8_t letter){
  const bool alphabetIsNucleotide     = index->metadata.alphabetType == AwFmAlphabetNucleotide;
  const size_t blockByteWidth         = alphabetIsNucleotide? sizeof(struct AwFmNucleotideBlock): sizeof( struct AwFmAminoBlock);
  uint64_t letterPrefixSum            = index->prefixSums[letter];

  union AwFmBwtBlockList blockList;
  uint64_t sentinelCharacterPosition;
  if(searchDirection == AwFmSearchDirectionBackward){
    blockList                 = index->backwardBwtBlockList;
    sentinelCharacterPosition = index->backwardSentinelCharacterPosition;
  }
  else{
    blockList                 = index->forwardBwtBlockList;
    sentinelCharacterPosition = index->forwardSentinelCharacterPosition;
  }

  //query for the start pointer
  //-1 relates to the "Occ(a,s-1) and OccLt(a,s-1) in the literature"
  uint64_t queryPosition = range->startPtr - 1;
  uint64_t blockIndex = queryPosition % POSITIONS_PER_FM_BLOCK;

  uint64_t baseOccurrence;
  struct AwFmOccurrenceVectorPair occurrenceVectors;

  uint64_t baseOccurrenceGte = 0;

  if(alphabetIsNucleotide){
    baseOccurrence = blockList.asNucleotide[blockIndex].baseOccurrences[letter];
    awFmMakeNucleotideOccurrenceVectorPair(&(blockList.asNucleotide[blockIndex]),
      queryPosition, letter, sentinelCharacterPosition, &occurrenceVectors);

    //compute the base occurrence of any letter greater than or equal to our query letter
    for(uint_fast8_t i = 0; i < 4; i++){
      baseOccurrenceGte += blockList.asNucleotide[blockIndex].baseOccurrences[i];
    }

  }
  else{
    baseOccurrence = blockList.asAmino[blockIndex].baseOccurrences[letter];
    awFmMakeAminoAcidOccurrenceVectorPair(&(blockList.asAmino[blockIndex]),
      queryPosition, letter, &occurrenceVectors);

    //compute the base occurrence of any letter greater than or equal to our query letter
    for(uint_fast8_t i = 0; i < 20; i++){
      baseOccurrenceGte += blockList.asAmino[blockIndex].baseOccurrences[i];
    }
  }

  uint_fast8_t  vectorPopcount  = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  uint64_t      newStartPointer = letterPrefixSum + vectorPopcount + baseOccurrence;

  //prefetch the next start ptr
  uint64_t newStartBlock    = (newStartPointer - 1) % POSITIONS_PER_FM_BLOCK;
  uint8_t *newStartBlockPtr = ((uint8_t*)blockList.asNucleotide) + (newStartBlock * blockByteWidth);
  for(size_t cacheLine = 0; cacheLine < 5; cacheLine++){
    _mm_prefetch(newStartBlockPtr + cacheLine, _MM_HINT_T2);
  }

  uint64_t occurrenceGte      = awFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t startOccurrenceLt  = queryPosition - occurrenceGte;


  //get the new end pointer
  queryPosition = range->endPtr;
  blockIndex = queryPosition % POSITIONS_PER_FM_BLOCK;

  baseOccurrenceGte = 0;

  if(alphabetIsNucleotide){
    baseOccurrence = blockList.asNucleotide[blockIndex].baseOccurrences[letter];
    awFmMakeNucleotideOccurrenceVectorPair(&(blockList.asNucleotide[blockIndex]),
      queryPosition, letter, sentinelCharacterPosition, &occurrenceVectors);

    //compute the base occurrence of any letter greater than or equal to our query letter
    for(uint_fast8_t i = 0; i < 4; i++){
      baseOccurrenceGte += blockList.asNucleotide[blockIndex].baseOccurrences[i];
    }
  }
  else{
    baseOccurrence = blockList.asAmino[blockIndex].baseOccurrences[letter];
    awFmMakeAminoAcidOccurrenceVectorPair(&(blockList.asAmino[blockIndex]),
      queryPosition, letter, &occurrenceVectors);

    //compute the base occurrence of any letter greater than or equal to our query letter
    for(uint_fast8_t i = 0; i < 20; i++){
      baseOccurrenceGte += blockList.asAmino[blockIndex].baseOccurrences[i];
    }
  }
  vectorPopcount = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  uint64_t newEndPointer    = letterPrefixSum + vectorPopcount + baseOccurrence - 1;

  //prefetch the next end ptr
  uint64_t newEndBlock    = newEndPointer % POSITIONS_PER_FM_BLOCK;
  uint8_t *newEndBlockPtr = ((uint8_t*)blockList.asNucleotide) + (newEndBlock * blockByteWidth);
  for(size_t cacheLine = 0; cacheLine < 5; cacheLine++){
    _mm_prefetch(newEndBlockPtr + cacheLine, _MM_HINT_T2);
  }

  occurrenceGte             = awFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t endOccurrenceLt  = queryPosition - occurrenceGte;
  uint64_t newStartPrimePtr = range->startPrimePtr  + endOccurrenceLt - startOccurrenceLt;

  range->startPtr       = newStartPointer;
  range->endPtr         = newEndPointer;
  range->startPrimePtr  = newStartPrimePtr;
}



void awFmIterativeStepBackwardSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBackwardRange *restrict const range, const uint8_t letter){

  const bool alphabetIsNucleotide     = index->metadata.alphabetType == AwFmAlphabetNucleotide;
  const size_t blockByteWidth         = alphabetIsNucleotide? sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);
  uint64_t letterPrefixSum            = index->prefixSums[letter];
  uint64_t sentinelCharacterPosition  = index->backwardSentinelCharacterPosition;

  const union AwFmBwtBlockList blockList = index->backwardBwtBlockList;

  //query for the start pointer
  uint64_t queryPosition  = range->startPtr - 1;
  uint64_t blockIndex = queryPosition % POSITIONS_PER_FM_BLOCK;

  uint64_t baseOccurrence;
  struct AwFmOccurrenceVectorPair occurrenceVectors;

  if(alphabetIsNucleotide){
    baseOccurrence = blockList.asNucleotide[blockIndex].baseOccurrences[letter];
    awFmMakeNucleotideOccurrenceVectorPair(&(blockList.asNucleotide[blockIndex]),
    queryPosition, letter, sentinelCharacterPosition, &occurrenceVectors);
  }
  else{
    baseOccurrence = blockList.asAmino[blockIndex].baseOccurrences[letter];
    awFmMakeAminoAcidOccurrenceVectorPair(&(blockList.asAmino[blockIndex]),
      queryPosition, letter, &occurrenceVectors);
  }
  uint_fast8_t  vectorPopcount  = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  uint64_t      newStartPointer = letterPrefixSum + vectorPopcount + baseOccurrence;

  //prefetch the next start ptr
  uint64_t newStartBlock    = (newStartPointer - 1) % POSITIONS_PER_FM_BLOCK;
  uint8_t *newStartBlockPtr = ((uint8_t*)blockList.asNucleotide) + (newStartBlock * blockByteWidth);
  for(size_t cacheLine = 0; cacheLine < 5; cacheLine++){
    _mm_prefetch(newStartBlockPtr + cacheLine, _MM_HINT_T2);
  }

  //query for the new end pointer
  queryPosition   = range->endPtr;
  blockIndex = queryPosition % POSITIONS_PER_FM_BLOCK;

  if(alphabetIsNucleotide){
    baseOccurrence = blockList.asNucleotide[blockIndex].baseOccurrences[letter];
    awFmMakeNucleotideOccurrenceVectorPair(&(blockList.asNucleotide[blockIndex]),
    queryPosition, letter, sentinelCharacterPosition, &occurrenceVectors);
  }
  else{
    baseOccurrence = blockList.asAmino[blockIndex].baseOccurrences[letter];
    awFmMakeAminoAcidOccurrenceVectorPair(&(blockList.asAmino[blockIndex]),
      queryPosition, letter, &occurrenceVectors);
  }

  vectorPopcount  = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  const uint64_t newEndPointer = letterPrefixSum + vectorPopcount + baseOccurrence - 1;

  //prefetch the next start ptr
  uint64_t newEndBlock    = (newEndPointer - 1) % POSITIONS_PER_FM_BLOCK;
  uint8_t *newEndBlockPtr = ((uint8_t*)blockList.asNucleotide) + (newEndBlock * blockByteWidth);
  for(size_t cacheLine = 0; cacheLine < 5; cacheLine++){
    _mm_prefetch(newEndBlockPtr + cacheLine, _MM_HINT_T2);
  }

  range->startPtr = newStartPointer;
  range->endPtr   = newEndPointer;
}


uint64_t *awFmFindDatabaseHitPositions(const struct AwFmIndex *restrict const index,
  const struct AwFmBackwardRange *restrict const searchRange, enum AwFmReturnCode *restrict fileAccessResult){

  const uint64_t numPositionsInRange  = awFmSearchRangeLength(searchRange);

  uint64_t *const restrict positionArray  = malloc(numPositionsInRange * sizeof(uint64_t));
  uint64_t *const restrict offsetArray    = malloc(numPositionsInRange * sizeof(uint64_t));
  //check for allocation failures
  if(__builtin_expect(positionArray == NULL, 0)){
    *fileAccessResult = AwFmAllocationFailure;
    return NULL;
  }
  if(__builtin_expect(offsetArray == NULL, 0)){
    free(positionArray);
    *fileAccessResult = AwFmAllocationFailure;
    return NULL;
  }

  //call a prefetch for each block that contains the positions that we need to start querying
  const uint_fast16_t blockWidth = index->metadata.alphabetType == AwFmAlphabetTypeNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

  for(uint64_t i = searchRange->startPtr; i < searchRange->endPtr; i += POSITIONS_PER_FM_BLOCK){
    awFmBlockPrefetch((uint8_t*)index->backwardBwtBlockList.asAmino, blockWidth, i);
  }


  //backtrace each position until we have a list of the positions in the database sequence.
  for(uint64_t i=0; i < numPositionsInRange; i++){
    uint64_t databaseSequenceOffset = 0;
    uint64_t backtracePosition      = searchRange->startPtr + i;

    while(!awFmBwtPositionIsSampled(index, backtracePosition)){
      backtracePosition = awFmBackstepBwtPosition(index, backtracePosition);
      databaseSequenceOffset++;
    }

    //position is divided by compression ratio to get the index in the suffix array
    positionArray[i]  = backtracePosition / index->metadata.suffixArrayCompressionRatio;
    offsetArray[i]    = databaseSequenceOffset;
  }

  //get the positions from the suffix array.
  *fileAccessResult = awFmReadPositionsFromSuffixArray(index, positionArray,
    numPositionsInRange);


  //add the offsets to the returned positions to get the actual positions of the hits
  for(size_t i = 0; i < numPositionsInRange; i++){
    positionArray[i] += offsetArray[i];
  }
  free(offsetArray);

  return positionArray;
}




struct AwFmBackwardRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint16_t kmerLength){

  uint_fast16_t kmerQueryPosition = kmerLength - 1;

  //set the initial range from the prefixSums
  const uint8_t firstSuffixLetter = kmer[kmerQueryPosition];
  struct AwFmBackwardRange searchRange = {index->prefixSums[firstSuffixLetter], index->prefixSums[firstSuffixLetter + 1] - 1};

  const uint_fast16_t blockWidth = (index->metadata.alphabetType == AwFmAlphabetNucleotide)?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

  //request the data be prefetched for each occurance call.
  awFmBlockPrefetch((uint8_t*)index->backwardBwtBlockList.asAmino, blockWidth, searchRange.startPtr - 1);
  awFmBlockPrefetch((uint8_t*)index->backwardBwtBlockList.asAmino, blockWidth, searchRange.endPtr);



  while(__builtin_expect(awFmSearchRangeIsValid(&searchRange) && kmerQueryPosition, 1)){
    kmerQueryPosition--;
    const uint8_t letter = kmer[kmerQueryPosition];
    awFmIterativeStepBackwardSearch(index, &searchRange, letter);
  }

  return searchRange;
}



bool awFmSingleKmerExists(const struct AwFmIndex *restrict const index, const char *restrict const kmer,
  const uint16_t kmerLength){

  struct AwFmBackwardRange kmerRange = awFmDatabaseSingleKmerExactMatch(index, kmer, kmerLength);
  return kmerRange.startPtr <= kmerRange.endPtr;
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
 size_t awFmSearchRangeLength(const struct AwFmBackwardRange *restrict const range){
  uint64_t length = range->endPtr - range->startPtr;
  return (range->startPtr < range->endPtr)? length + 1: 0;
}


/*
 * Function:  awFmSwapBiDirectionalRangePointerDirection
 * --------------------
 * Swaps the direction of the bidirectional range struct in order to search in the opposite direction.
 *  This function works in place on the given range, swapping the startPtr and startPrimePtr, and
 *  updating the endPtr to reflect the correct end for the new startPtr.
 *
 *  Inputs:
 *    range: bidirectional range in the Bwt. This act as an out-argument, and will be updated by this function.
 */
void awFmSwapBiDirectionalRangePointerDirection(struct AwFmBiDirectionalRange *restrict const range){
  uint64_t rangeSize          = range->endPtr - range->startPtr;
  uint64_t tempStartPrimePtr  = range->startPtr;

  range->startPtr       = range->startPrimePtr;
  range->endPtr         = range->startPtr + rangeSize;
  range->startPrimePtr  = tempStartPrimePtr;
}
