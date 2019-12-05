#include "AwFmIndex.h"
#include "AwFmOccurrence.h"
#include "AwFmSearch.h"


void awFmIterativeStepBidirectionalNucleotideSearch(const struct AwFmIndex *restrict const index,
  enum AwFmSearchDirection searchDirection, struct AwFmBiDirectionalRange *restrict const range, const uint8_t letter){

  uint64_t letterprefixSum            = index->prefixSums[letter];
  uint64_t sentinelCharacterPosition  = index->sentinelCharacterPosition;
  const struct AwFmNucleotideBlock *restrict const blockList = (searchDirection == AwFmSearchDirectionBackward)?
    index->backwardBwtBlockList.asNucleotide:
    index->forwardBwtBlockList.asNucleotide;

  //query for the start pointer
  //-1 relates to the "Occ(a,s-1) and OccLt(a,s-1) in the literature"
  uint64_t queryPosition = range->startPtr - 1;
  const struct AwFmNucleotideBlock *restrict blockPtr = awFmGetNucleotideBlockPtr(index, searchDirection, queryPosition);
  struct AwFmOccurrenceVectorPair occurrenceVectors;
  awFmMakeNucleotideOccurrenceVectorPair(blockPtr, queryPosition, letter, sentinelCharacterPosition, &occurrenceVectors);

  uint64_t baseOccurrence     = blockPtr->baseOccurrences[letter];
  uint_fast8_t vectorPopcount = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  uint64_t newStartPointer    = letterprefixSum + vectorPopcount + baseOccurrence;
  //like above, the next start query will get data from the position 1 before the startPtr
  awFmBlockPrefetch((uint8_t*)blockList, sizeof(struct AwFmNucleotideBlock), newStartPointer - 1);


  uint64_t baseOccurrenceGte  = awFmGetNucleotideBaseOccurrenceGte(blockPtr, letter);
  uint64_t occurrenceGte      = awFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t startOccurrenceLt  = queryPosition - occurrenceGte;


  //query for the end pointer
  queryPosition = range->endPtr;
  blockPtr = awFmGetNucleotideBlockPtr(index, searchDirection, queryPosition);
  awFmMakeNucleotideOccurrenceVectorPair(blockPtr, queryPosition, letter, sentinelCharacterPosition, &occurrenceVectors);


  baseOccurrence = blockPtr->baseOccurrences[letter];
  vectorPopcount = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  //In the literature, e = C[a] + Occ(a,e) - 1, so the -1 in the next line corresponds to the -1 here.
  uint64_t newEndPointer    = letterprefixSum + vectorPopcount + baseOccurrence - 1;
  awFmBlockPrefetch((uint8_t*) blockList, sizeof(struct AwFmNucleotideBlock), newEndPointer);

  baseOccurrenceGte         = awFmGetNucleotideBaseOccurrenceGte(blockPtr, letter);
  occurrenceGte             = awFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t endOccurrenceLt  = queryPosition - occurrenceGte;
  uint64_t newStartPrimePtr = range->startPrimePtr  + endOccurrenceLt - startOccurrenceLt;

  range->startPtr           = newStartPointer;
  range->endPtr             = newEndPointer;
  range->startPrimePtr      = newStartPrimePtr;
}


void awFmIterativeStepBidirectionalAminoAcidSearch(const struct AwFmIndex *restrict const index,
  enum AwFmSearchDirection searchDirection, struct AwFmBiDirectionalRange *restrict const range, const uint8_t letter){

  uint64_t letterprefixSum     = index->prefixSums[letter];

  const struct AwFmAminoBlock *restrict const blockList = (searchDirection == AwFmSearchDirectionBackward)?
    index->backwardBwtBlockList.asAmino:
    index->forwardBwtBlockList.asAmino;

  //query for the start pointer
  //-1 relates to the "Occ(a,s-1) and OccLt(a,s-1) in the literature"
  uint64_t queryPosition = range->startPtr - 1;
  const struct AwFmAminoBlock *restrict blockPtr = awFmGetAminoBlockPtr(index, searchDirection, queryPosition);
  struct AwFmOccurrenceVectorPair occurrenceVectors;
  awFmMakeAminoAcidOccurrenceVectorPair(blockPtr, queryPosition, letter, &occurrenceVectors);

  uint64_t startBaseOccurrence  = blockPtr->baseOccurrences[letter];
  uint64_t vectorPopcount       = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  uint64_t newStartPointer      = letterprefixSum + vectorPopcount + startBaseOccurrence;
  //like above, the next start query will get data from the position 1 before the startPtr
  awFmBlockPrefetch((uint8_t*) blockList, sizeof(struct AwFmAminoBlock), newStartPointer - 1);

  uint64_t baseOccurrenceGte    = awFmGetAminoAcidBaseOccurrenceGte(blockPtr, letter);
  uint64_t occurrenceGte        = awFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t startOccurrenceLt    = queryPosition - occurrenceGte;

  //query for the end pointer
  queryPosition = range->endPtr;
  blockPtr = awFmGetAminoBlockPtr(index, searchDirection, queryPosition);
  awFmMakeAminoAcidOccurrenceVectorPair(blockPtr, queryPosition, letter , &occurrenceVectors);

  startBaseOccurrence = blockPtr->baseOccurrences[letter];
  vectorPopcount = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  //In the literature, e = C[a] + Occ(a,e) - 1, so the -1 in the next line corresponds to the -1 here.
  uint64_t newEndPointer    = letterprefixSum + vectorPopcount + startBaseOccurrence - 1;
  awFmBlockPrefetch((uint8_t*) blockList, sizeof(struct AwFmAminoBlock), newEndPointer);

  baseOccurrenceGte         = awFmGetAminoAcidBaseOccurrenceGte(blockPtr, letter);
  occurrenceGte             = awFmVectorPopcount(occurrenceVectors.occurrenceGteVector) + baseOccurrenceGte;
  uint64_t endOccurrenceLt  = queryPosition - occurrenceGte;

  uint64_t newStartPrimePtr = range->startPrimePtr  + endOccurrenceLt - startOccurrenceLt;

  range->startPtr           = newStartPointer;
  range->endPtr             = newEndPointer;
  range->startPrimePtr      = newStartPrimePtr;
}


void awFmIterativeStepBackwardNucleotideSearch(const struct AwFmIndex *restrict const index,
struct AwFmBackwardRange *restrict const range, const uint8_t letter){

  uint64_t letterprefixSum        = index->prefixSums[letter];
  uint64_t sentinelCharacterPosition  = index->sentinelCharacterPosition;

  //query for the start pointer
  uint64_t queryPosition  = range->startPtr - 1;
  const struct AwFmNucleotideBlock *restrict blockPtr = awFmGetNucleotideBlockPtr(index, AwFmSearchDirectionBackward, queryPosition);
  uint64_t baseOccurrence = blockPtr->baseOccurrences[letter];
  struct AwFmOccurrenceVectorPair occurrenceVectors;
  awFmMakeNucleotideOccurrenceVectorPair(blockPtr, queryPosition, letter, sentinelCharacterPosition, &occurrenceVectors);

  uint64_t vectorPopcount = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  range->startPtr         = letterprefixSum + vectorPopcount + baseOccurrence;
  awFmBlockPrefetch((uint8_t*)index->backwardBwtBlockList.asNucleotide, sizeof(struct AwFmNucleotideBlock), range->startPtr - 1);

  //query for the end pointer
  queryPosition   = range->endPtr;
  blockPtr        = awFmGetNucleotideBlockPtr(index, AwFmSearchDirectionBackward, queryPosition);
  baseOccurrence  = blockPtr->baseOccurrences[letter];
  awFmMakeNucleotideOccurrenceVectorPair(blockPtr, queryPosition, letter, sentinelCharacterPosition, &occurrenceVectors);

  vectorPopcount  = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  range->endPtr   = letterprefixSum + vectorPopcount + baseOccurrence - 1;
  awFmBlockPrefetch((uint8_t*) index->backwardBwtBlockList.asNucleotide, sizeof(struct AwFmNucleotideBlock), range->endPtr);
}


void awFmIterativeStepBackwardAminoAcidSearch(const struct AwFmIndex *restrict const index,
  struct AwFmBackwardRange *restrict const range, const uint8_t letter){

  uint64_t letterprefixSum        = index->prefixSums[letter];

  //query for the start pointer
  uint64_t queryPosition = range->startPtr - 1;
  const struct AwFmAminoBlock *restrict blockPtr = awFmGetAminoBlockPtr(index, AwFmSearchDirectionBackward, queryPosition);
  uint64_t baseOccurrence = blockPtr->baseOccurrences[letter];
  struct AwFmOccurrenceVectorPair occurrenceVectors;
  awFmMakeAminoAcidOccurrenceVectorPair(blockPtr, queryPosition, letter, &occurrenceVectors);

  uint64_t vectorPopcount = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  range->startPtr = letterprefixSum + vectorPopcount + baseOccurrence;
  awFmBlockPrefetch((uint8_t*) index->backwardBwtBlockList.asAmino, sizeof(struct AwFmAminoBlock), range->startPtr - 1);


  //query for the end pointer
  queryPosition     = range->endPtr;
  blockPtr          = awFmGetAminoBlockPtr(index, AwFmSearchDirectionBackward, queryPosition);
  baseOccurrence    = blockPtr->baseOccurrences[letter];
  awFmMakeAminoAcidOccurrenceVectorPair(blockPtr, queryPosition, letter, &occurrenceVectors);

  vectorPopcount = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
  range->endPtr = letterprefixSum + vectorPopcount + baseOccurrence - 1;
  awFmBlockPrefetch((uint8_t*) index->backwardBwtBlockList.asAmino, sizeof(struct AwFmAminoBlock), range->endPtr);
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
  const struct AwFmBackwardRange *restrict const searchRange,
  enum AwFmReturnCode *restrict fileAccessResult){

  const uint64_t numPositionsInRange  = awFmSearchRangeLength(searchRange);

  uint64_t *const restrict positionArray  = malloc(numPositionsInRange * sizeof(uint64_t));
  uint64_t *const restrict offsetArray    = malloc(numPositionsInRange * sizeof(uint64_t));
  //check for allocation failures
  if(__builtin_expect(positionArray == NULL, 0)){
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

    //set the initial range from the prefixSums
    const uint8_t firstSuffixLetter = kmer[kmerQueryPosition];
    struct AwFmBackwardRange searchRange = {index->prefixSums[firstSuffixLetter] + 1, index->prefixSums[firstSuffixLetter + 1]};

    const uint_fast16_t blockWidth = (index->metadata.alphabetType == AwFmAlphabetNucleotide)?
      sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

    //request the data be prefetched for each occurance call.
    awFmBlockPrefetch((uint8_t*)index->backwardBwtBlockList.asAmino, blockWidth, searchRange.startPtr);
    awFmBlockPrefetch((uint8_t*)index->backwardBwtBlockList.asAmino, blockWidth, searchRange.endPtr);

    if(isNucleotideDatabase){
      while(__builtin_expect(awFmSearchRangeIsValid(&searchRange) && kmerQueryPosition, 1)){
        kmerQueryPosition--;
        const uint8_t letter = kmer[kmerQueryPosition];
        awFmIterativeStepBackwardNucleotideSearch(index, &searchRange,letter);
      }
    }
    else{
      while(__builtin_expect(awFmSearchRangeIsValid(&searchRange) && kmerQueryPosition, 1)){
        kmerQueryPosition--;
        const uint8_t letter = kmer[kmerQueryPosition];
        awFmIterativeStepBackwardAminoAcidSearch(index, &searchRange,letter);
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
 size_t awFmSearchRangeLength(const struct AwFmBackwardRange *restrict const range){
  uint64_t length = range->endPtr - range->startPtr;
  return (range->startPtr < range->endPtr)? length: 0;
}
