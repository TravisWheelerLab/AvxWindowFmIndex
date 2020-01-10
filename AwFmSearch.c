#include "AwFmSearch.h"
#include "AwFmOccurrence.h"
#include "AwFmLetter.h"
#include <assert.h>


bool awFmNucleotidePopcountIncludesSentinelCharacter(const uint8_t queryLetter,
  const uint64_t globalSentinelPosition, const uint64_t globalQueryPosition);

void awFmNucleotideIterativeStepBackwardSearch(const struct AwFmIndex *restrict const index,
  struct AwFmSearchRange *restrict const range, const uint8_t letter){

  //query for the start pointer
  uint64_t queryPosition  = range->startPtr - 1;
  uint64_t blockIndex = queryPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
  const uint64_t letterPrefixSum = index->prefixSums[letter];

  uint64_t baseOccurrence = index->bwtBlockList.asNucleotide[blockIndex].baseOccurrences[letter];
  __m256i occurrenceVector = awFmMakeNucleotideOccurrenceVector(&(index->bwtBlockList.asNucleotide[blockIndex]),
    queryPosition, letter);
  uint_fast8_t  vectorPopcount  = awFmVectorPopcount(occurrenceVector);

  //TODO: check for sentinel here and subtract 1 if appropriate
  if(awFmNucleotidePopcountIncludesSentinelCharacter(letter, index->sentinelCharacterPosition, queryPosition)){
    vectorPopcount--;
  }

  uint64_t newStartPointer = letterPrefixSum + vectorPopcount + baseOccurrence;

  //prefetch the next start ptr
  uint64_t newStartBlock    = (newStartPointer - 1) % AW_FM_POSITIONS_PER_FM_BLOCK;
  uint8_t *newStartBlockPtr = ((uint8_t*)index->bwtBlockList.asNucleotide) + (newStartBlock * sizeof(struct AwFmNucleotideBlock));
  for(size_t cacheLine = 0; cacheLine < 5; cacheLine++){
    _mm_prefetch(newStartBlockPtr + cacheLine, _MM_HINT_T2);
  }

  //query for the new end pointer
  queryPosition   = range->endPtr;
  blockIndex = queryPosition % AW_FM_POSITIONS_PER_FM_BLOCK;

  baseOccurrence = index->bwtBlockList.asNucleotide[blockIndex].baseOccurrences[letter];
  occurrenceVector = awFmMakeNucleotideOccurrenceVector(&(index->bwtBlockList.asNucleotide[blockIndex]),
    queryPosition, letter);
  vectorPopcount  = awFmVectorPopcount(occurrenceVector);

  if(awFmNucleotidePopcountIncludesSentinelCharacter(letter, index->sentinelCharacterPosition, queryPosition)){
    vectorPopcount--;
  }
  const uint64_t newEndPointer = letterPrefixSum + vectorPopcount + baseOccurrence - 1;

  //prefetch the next start ptr
  uint64_t newEndBlock    = (newEndPointer - 1) % AW_FM_POSITIONS_PER_FM_BLOCK;
  uint8_t *newEndBlockPtr = ((uint8_t*)index->bwtBlockList.asNucleotide) + (newEndBlock * sizeof(struct AwFmNucleotideBlock));
  for(size_t cacheLine = 0; cacheLine < 5; cacheLine++){
    _mm_prefetch(newEndBlockPtr + cacheLine, _MM_HINT_T2);
  }

  range->startPtr = newStartPointer;
  range->endPtr   = newEndPointer;
}


void awFmAminoIterativeStepBackwardSearch(const struct AwFmIndex *restrict const index,
  struct AwFmSearchRange *restrict const range, const uint8_t letter){

  //query for the start pointer
  uint64_t queryPosition  = range->startPtr - 1;
  uint64_t blockIndex     = queryPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
  const uint64_t letterPrefixSum = index->prefixSums[letter];

  uint64_t baseOccurrence   = index->bwtBlockList.asAmino[blockIndex].baseOccurrences[letter];
  __m256i occurrenceVector  = awFmMakeAminoAcidOccurrenceVector(&(index->bwtBlockList.asAmino[blockIndex]),
    queryPosition, letter);
  uint_fast8_t  vectorPopcount  = awFmVectorPopcount(occurrenceVector);

  uint64_t      newStartPointer = letterPrefixSum + vectorPopcount + baseOccurrence;

  //prefetch the next start ptr
  uint64_t newStartBlock    = (newStartPointer - 1) % AW_FM_POSITIONS_PER_FM_BLOCK;
  uint8_t *newStartBlockPtr = ((uint8_t*)index->bwtBlockList.asNucleotide) + (newStartBlock * sizeof(struct AwFmAminoBlock));
  for(size_t cacheLine = 0; cacheLine < 5; cacheLine++){
    _mm_prefetch(newStartBlockPtr + cacheLine, _MM_HINT_T2);
  }

  //query for the new end pointer
  queryPosition = range->endPtr;
  blockIndex    = queryPosition % AW_FM_POSITIONS_PER_FM_BLOCK;

  baseOccurrence    = index->bwtBlockList.asAmino[blockIndex].baseOccurrences[letter];
  occurrenceVector  = awFmMakeAminoAcidOccurrenceVector(&(index->bwtBlockList.asAmino[blockIndex]),
    queryPosition, letter);
  vectorPopcount    = awFmVectorPopcount(occurrenceVector);

  const uint64_t newEndPointer = letterPrefixSum + vectorPopcount + baseOccurrence - 1;

  //prefetch the next start ptr
  uint64_t newEndBlock    = (newEndPointer - 1) % AW_FM_POSITIONS_PER_FM_BLOCK;
  uint8_t *newEndBlockPtr = ((uint8_t*)index->bwtBlockList.asNucleotide) + (newEndBlock * sizeof(struct AwFmAminoBlock));
  for(size_t cacheLine = 0; cacheLine < 5; cacheLine++){
    _mm_prefetch(newEndBlockPtr + cacheLine, _MM_HINT_T2);
  }

  range->startPtr = newStartPointer;
  range->endPtr   = newEndPointer;
}


uint64_t *awFmFindDatabaseHitPositions(const struct AwFmIndex *restrict const index,
  const struct AwFmSearchRange *restrict const searchRange, enum AwFmReturnCode *restrict fileAccessResult){

  const uint64_t numPositionsInRange  = awFmSearchRangeLength(searchRange);

  // if there were no elements in the search range, abandon the query.
  if(__builtin_expect(numPositionsInRange == 0, 0)){
    *fileAccessResult = AwFmGeneralFailure;
    return NULL;
  }

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
  const uint_fast16_t blockWidth = index->metadata.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

  for(uint64_t i = searchRange->startPtr; i < searchRange->endPtr; i += AW_FM_POSITIONS_PER_FM_BLOCK){
    awFmBlockPrefetch((uint8_t*)index->bwtBlockList.asAmino, blockWidth, i);
  }

  //backtrace each position until we have a list of the positions in the database sequence.
  for(uint64_t i=0; i < numPositionsInRange; i++){
    uint64_t databaseSequenceOffset = 0;
    uint64_t backtracePosition      = searchRange->startPtr + i;

    if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
      while(!awFmBwtPositionIsSampled(index, backtracePosition)){
        backtracePosition = awFmNucleotideBacktraceBwtPosition(index, backtracePosition);
        databaseSequenceOffset++;
      }
    }
    else{
      while(!awFmBwtPositionIsSampled(index, backtracePosition)){
        backtracePosition = awFmAminoBacktraceBwtPosition(index, backtracePosition);
        databaseSequenceOffset++;
      }
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

  *fileAccessResult = AwFmFileReadOkay;
  return positionArray;
}


struct AwFmSearchRange awFmDatabaseSingleKmerExactMatch(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint16_t kmerLength){
    //todo: use kmer lookup table
    assert(false);  //todo: implement

  uint_fast16_t kmerQueryPosition = kmerLength - 1;

  //set the initial range from the prefixSums
  const uint8_t firstSuffixLetter = kmer[kmerQueryPosition];
  struct AwFmSearchRange searchRange = {index->prefixSums[firstSuffixLetter], index->prefixSums[firstSuffixLetter + 1] - 1};

  const uint_fast16_t blockWidth = (index->metadata.alphabetType == AwFmAlphabetNucleotide)?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

  //request the data be prefetched for each occurance call.
  awFmBlockPrefetch((uint8_t*)index->bwtBlockList.asAmino, blockWidth, searchRange.startPtr - 1);
  awFmBlockPrefetch((uint8_t*)index->bwtBlockList.asAmino, blockWidth, searchRange.endPtr);

  if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
    while(__builtin_expect(awFmSearchRangeIsValid(&searchRange) && kmerQueryPosition, 1)){
      kmerQueryPosition--;
      const uint8_t letter = kmer[kmerQueryPosition];
      awFmNucleotideIterativeStepBackwardSearch(index, &searchRange, letter);
    }
  }
  else{
    while(__builtin_expect(awFmSearchRangeIsValid(&searchRange) && kmerQueryPosition, 1)){
      kmerQueryPosition--;
      const uint8_t letter = kmer[kmerQueryPosition];
      awFmAminoIterativeStepBackwardSearch(index, &searchRange, letter);
    }
  }

  return searchRange;
}


bool awFmSingleKmerExists(const struct AwFmIndex *restrict const index, const char *restrict const kmer,
  const uint16_t kmerLength){

  struct AwFmSearchRange kmerRange = awFmDatabaseSingleKmerExactMatch(index, kmer, kmerLength);
  return kmerRange.startPtr <= kmerRange.endPtr;
}


inline size_t awFmNucleotideBacktraceBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition){
  const uint64_t *prefixSums               = index->prefixSums;
  const uint64_t blockIndex                = awFmGetBlockIndexFromGlobalPosition(bwtPosition);
  const uint8_t localQueryPosition         = awFmGetBlockQueryPositionFromGlobalPosition(bwtPosition);

  const struct AwFmNucleotideBlock *restrict const blockPtr = &index->bwtBlockList.asNucleotide[blockIndex];

  const uint8_t frequencyIndexLetter  = awFmGetNucleotideLetterAtBwtPosition(blockPtr, localQueryPosition);
  const __m256i occurrenceVector      = awFmMakeNucleotideOccurrenceVector(blockPtr, localQueryPosition,
    frequencyIndexLetter);


  uint_fast8_t vectorPopcount = awFmVectorPopcount(occurrenceVector);
  if(awFmNucleotidePopcountIncludesSentinelCharacter(frequencyIndexLetter,
    index->sentinelCharacterPosition, bwtPosition)){
    vectorPopcount--;
  }
  uint64_t backtraceBwtPosition = prefixSums[frequencyIndexLetter] + vectorPopcount;

  return backtraceBwtPosition;
}


inline size_t awFmAminoBacktraceBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition){
  const uint64_t  *prefixSums         = index->prefixSums;
  const uint64_t  blockIndex          = awFmGetBlockIndexFromGlobalPosition(bwtPosition);
  const uint8_t   localQueryPosition  = awFmGetBlockQueryPositionFromGlobalPosition(bwtPosition);

  const struct AwFmAminoBlock *restrict const blockPtr = &index->bwtBlockList.asAmino[blockIndex];

  uint8_t frequencyIndexLetter  = awFmGetAminoLetterAtBwtPosition(blockPtr, bwtPosition);
  __m256i occurrenceVector      = awFmMakeAminoAcidOccurrenceVector(blockPtr,
    localQueryPosition, frequencyIndexLetter);

  const uint_fast8_t vectorPopcount   = awFmVectorPopcount(occurrenceVector);
  const uint64_t backtraceBwtPosition = prefixSums[frequencyIndexLetter] + vectorPopcount;

  return backtraceBwtPosition;
}


inline bool awFmNucleotidePopcountIncludesSentinelCharacter(const uint8_t queryLetter,
  const uint64_t globalSentinelPosition, const uint64_t globalQueryPosition){
    const uint64_t  sentinelBlockIndex    = awFmGetBlockIndexFromGlobalPosition(globalSentinelPosition);
    const uint64_t  queryBlockIndex       = awFmGetBlockIndexFromGlobalPosition(globalQueryPosition);
    const uint8_t   localSentinelPosition = awFmGetBlockQueryPositionFromGlobalPosition(globalSentinelPosition);
    const uint8_t   localQueryPosition    = awFmGetBlockQueryPositionFromGlobalPosition(globalQueryPosition);

  return  (queryLetter == 0) &&
          (sentinelBlockIndex == queryBlockIndex) &&
          (localSentinelPosition <= localQueryPosition);
}
