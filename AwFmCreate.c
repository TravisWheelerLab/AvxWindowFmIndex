#include "AwFmCreate.h"
#include "AwFmLetter.h"
#include "AwFmIndex.h"
#include <stdlib.h>
#include <string.h>
#include <divsufsort64.h>

/*Creates the fully-sampled suffix array from the given database sequence.*/
// enum AwFmReturnCode awFmCreateFullSuffixArray(const uint8_t *databaseSequence,
//   const uint64_t databaseSequenceLength, uint64_t **fullSuffixArrayOut);

/*Creates the AwFmIndex blockList from the given database sequence.*/
// enum AwFmReturnCode awFmCreateBlockList(struct AwFmIndex *restrict const index,
//   const uint8_t *restrict const databaseSequence, const uint64_t databaseSequenceLength,
//   uint64_t *totalOccupancies);

/*sets the rankPrefixSums array. This function requires that the blockList and suffix array already be set..*/
// void AwFmSetRankPrefixSums(struct AwFmIndex *restrict const index, const uint64_t *restrict const totalOccupancies);

// void awFmInitBlock(struct AwFmIndex *const restrict index, const uint64_t blockIndex, uint64_t *totalOccupanciesSoFar);

/*
 * Function:  awFmCreateIndex
 * --------------------
 * Create an AwFmIndex struct from a given database sequence. This function will deallocate
 *  all data it allocated if any part fails.
 *  Note: databaseSequence is not copied, instead, the AwFmIndex copies the pointer to the string.
 *  As such, the caller is responsible for the lifetime of this string, and it should not be
 *  deallocated until the AwFmIndex is written to file.
 *
 *  Inputs:
 *    indexPtr:                     Pointer to the AwFmIndex ptr to be allocated.
 *    databaseSequence:             ASCII-encoded string representing the
 *      amino-acid sequence to build the database from.
 *      NOTE: this is passes by reference (pointer copy). The caller is expected to maintain
 *      the lifetime of the databaseSequence string, and assuming that the AwFmIndex struct
 *      is going to be written to a file, the caller is expected to not deallocate this string until
 *      the AwFmIndex is written to file.
 *    databaseSequenceLength:       length of the databaseSequence string.
 *    suffixArrayCompressionRatio:  how often to sample the suffix array. Larger values
 *      reduce the size of the Fm-index data file, but increase the time it takes to reconstruct
 *      the database sequence positions when hits are found.
 *    fileSrc: null terminated string representing the file path to the file associated with this
 *      AwFmIndex. While this function does not write the AwFmIndex to a file, this file path
 *      will be what's used to write it.
 *
 *  Returns:
 *    AwFmReturnCode showing the result of this action. Possible returns are:
 *      AwFmSuccess on success,
 *      AwFmNullPtrError if the indexPtr, databaseSequence, or fileSrc were null.
 *      AwFmAllocationFailure on failure to allocate the AwFmIndex, suffix array copy,
 *        or block list.
 *      AwFmNoFileSrcGiven if the fileSrc is null.
 *      AwFmGeneralFailure if a null terminator is not found in the fileSrc in a
 *        reasonable window (defined by macro MAXIMUM_FILE_PATH_LENGTH).
 *      AwFmSuffixArrayCreationFailure on failure in building the suffix array (this error
 *        is caused by internal failure in the libdivsufsort library).
*/
enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex **indexPtr, const uint8_t *restrict const databaseSequence,
  const size_t databaseSequenceLength, const uint16_t suffixArrayCompressionRatio, const char *restrict const fileSrc){

  if(indexPtr == NULL || databaseSequence == NULL || fileSrc == NULL){
    return AwFmNullPtrError;
  }

  struct AwFmIndex *index = awFmAlignedAllocAwFmIndex();
  if(index == NULL){
    return AwFmAllocationFailure;
  }

  //set the suffix array compression ratio
  index->suffixArrayCompressionRatio = suffixArrayCompressionRatio;

  //set the fileSrc, or return the error code if it failed.
  if(fileSrc != NULL){
    const enum AwFmReturnCode setFileSrcReturnCode = awFmIndexSetFileSrc(index, fileSrc);
    if(__builtin_expect(setFileSrcReturnCode < 0, 0)){
      awFmDeallocateFmIndex(index);
      return setFileSrcReturnCode;
    }
  }

  //set the database sequence by reference
  if(databaseSequence == NULL){
    awFmDeallocateFmIndex(index);
    return AwFmNullPtrError;
  }
  index->databaseSequence = databaseSequence;

  //create the full suffix array
  enum AwFmReturnCode suffixArrayCreationReturnCode = awFmCreateFullSuffixArray(databaseSequence,
    databaseSequenceLength, &index->fullSuffixArray);
  if(__builtin_expect(suffixArrayCreationReturnCode < 0, 0)){
    awFmDeallocateFmIndex(index);
    return suffixArrayCreationReturnCode;
  }

  //the occupancy values will accumulate into this array.
  // the +1 leaves one slot at the end for unrecognized letters to accumulate into,
  // and will be treated as a garbage data slot so that illegal characters don't crash the application.
  uint64_t totalOccupancies[AMINO_CARDINALITY + 1];
  //create the block list
  enum AwFmReturnCode blockListCreationReturnCode = awFmCreateBlockList(index,
    databaseSequence, databaseSequenceLength, totalOccupancies);

  if(__builtin_expect(blockListCreationReturnCode < 0, 0)){
    awFmDeallocateFmIndex(index);
    return blockListCreationReturnCode;
  }

  AwFmSetRankPrefixSums(index, totalOccupancies);

  //set the out pointer for the index.
  *indexPtr = index;

  return AwFmSuccess;
}


void AwFmSetRankPrefixSums(struct AwFmIndex *restrict const index, const uint64_t *restrict const totalOccupancies){
  //set the rank prefix sums
  //accumulator starts at 1 because of the null terminator $ being before all other letters
  uint64_t rankPrefixAccumulator = 1;

  index->rankPrefixSums[0] = rankPrefixAccumulator;
  for(uint_fast8_t i = 0; i < AMINO_CARDINALITY; i++){
    rankPrefixAccumulator += totalOccupancies[i];
    index->rankPrefixSums[i + 1] = rankPrefixAccumulator;
  }
}



//if failure, the caller may need to deallocate index.
/*
 * Function:  awFmCreateBlockList
 * --------------------
 * Allocates and initializes the BlockList with the suffix array created by
 *  the given database sequence. If this function returns a failure value (any
 *  value < 0), the given index should be deallocated.
 *
 *  Inputs:
 *    index              AwFmIndex struct that's block list will be set.
 *    databaseSequence:             ASCII-encoded string representing the
 *      amino-acid sequence to build the blockList from.
 *    databaseSequenceLength:       length of the databaseSequence string to be copied
 *      into the AwFmIndex struct.
 *    totalOccupancies:   Array of at least length 21 (AMINO_CARDINALITY + 1) whose
 *      values will be assigned the full occupancies for each amino over the entire database.
 *      the array has one final element at the end so that the null terminator has
 *      somewhere to accumulate to, which will later be thrown away.
 *
 *  Returns:
 *    AwFmReturnCode showing the result of this action. Possible returns are:
 *      AwFmSuccess on success,
 *      AwFmNullPtrError if the indexPtr, databaseSequence, or totalOccupancies were null.
 *      AwFmAllocationFailure on failure to allocate the block list.
*/
enum AwFmReturnCode awFmCreateBlockList(struct AwFmIndex *const restrict index,
  const uint8_t *restrict const databaseSequence, const uint64_t databaseSequenceLength,
  uint64_t *totalOccupancies){

  if(index == NULL || databaseSequence == NULL || totalOccupancies == NULL){
    return AwFmNullPtrError;
  }

  const uint64_t suffixArrayLength = databaseSequenceLength + 1;

  //ceiling function to compute the number of blocks needed
  //I would normally write this in the form (x+(y-1)/y), but written like this
  // avoids the potential of overflow.
  index->numBlocks = awFmNumBlocksFromSuffixArrayLength(suffixArrayLength);
  index->blockList = awFmAlignedAllocBlockList(index->numBlocks);
  if(index->blockList == NULL){
    return AwFmAllocationFailure;
  }

  //keep an accumulator for the baseOccupancies and current position in suffix array
  //by adding 1 to the length, anything that parses incorrectly won't crash the program,
  //but will accumulate into an unused value
  memset(totalOccupancies, 0,  AMINO_CARDINALITY * sizeof(uint64_t));

  for(uint64_t blockIndex = 0; blockIndex < index->numBlocks; blockIndex++){
    awFmInitBlock(index, blockIndex, totalOccupancies, suffixArrayLength);
  }

  return AwFmSuccess;
}


void awFmInitBlock(struct AwFmIndex *const restrict index, const uint64_t blockIndex,
  uint64_t *totalOccupanciesSoFar, const size_t suffixArrayLength){

  //copy the base occupancy over, and set the bit vectors to all zeros.
  memcpy(index->blockList[blockIndex].baseOccupancies, totalOccupanciesSoFar, AMINO_CARDINALITY * sizeof(uint64_t));
  memset(index->blockList[blockIndex].letterBitVectors, 0, sizeof(__m256i) * 5);

  //create a uint8_t pointer to the letter bit vectors to more easily set individual bytes.
  void *letterBitVectorsAsRawPtr = &index->blockList[blockIndex].letterBitVectors[0];
  uint64_t *letterBitVectorsAsBytes = letterBitVectorsAsRawPtr;

  bool isLastBlock = blockIndex == index->numBlocks;
  //if we're in the last block, we need to stop at the actual end of the suffix array.
  uint_fast8_t lengthOfThisBlock = isLastBlock? (suffixArrayLength % POSITIONS_PER_FM_BLOCK) - 1: POSITIONS_PER_FM_BLOCK;

  //if we're on the first block, start at position 1, since 0 will be the null terminator.
  const uint_fast8_t blockStartingPosition = (blockIndex == 0)? 1: 0;
  for(uint_fast8_t i = blockStartingPosition; i < lengthOfThisBlock; i++){

    const uint8_t byteInBlock             = i / 8;
    const uint8_t bitInBlockByte          = i % 8;
    const uint64_t positionInSuffixArray  = (blockIndex * POSITIONS_PER_FM_BLOCK) + i;

    const uint64_t suffixArrayValue       = index->fullSuffixArray[positionInSuffixArray];
    const uint8_t letterAsAscii           = index->databaseSequence[suffixArrayValue];
    const uint8_t letterAsFrequencyIndex  = awFmAsciiLetterToLetterIndex(letterAsAscii);
    uint8_t letterAsVectorFormat          = awFmAsciiLetterToCompressedVectorFormat(letterAsAscii);

    //accumulate the letter's occupancy.
    totalOccupanciesSoFar[letterAsFrequencyIndex]++;

    //set the correct bits in the letterBitVectors
    for(uint_fast8_t bitInVectorLetter = 0; bitInVectorLetter < 5; bitInVectorLetter++){
      const uint8_t letterBit = (letterAsVectorFormat >> bitInVectorLetter) & 1;
      const uint8_t byteInBlockVectors = (letterBit * BYTES_PER_AVX2_REGISTER) + byteInBlock;
      //or in the shifted bit.
      letterBitVectorsAsBytes[byteInBlockVectors] |= (letterBit << bitInBlockByte);
    }
  }
}


/*
 * Function:  awFmCreateFullSuffixArray
 * --------------------
 * Creates a fully-sampled suffix array from the given database sequence.
 *
 *  Inputs:
 *    databaseSequence:       ASCII-encoded string representing the
 *      amino-acid sequence to build the blockList from.
 *    databaseSequenceLength: Length of the databaseSequence string to be copied
 *      into the AwFmIndex struct.
 *    fullSuffixArrayOut:     Pointer to the suffixArray ptr to be allocated and set
 *      to the computed suffix array.
 *
 *  Returns:
 *    AwFmReturnCode showing the result of this action. Possible returns are:
 *      AwFmSuccess on success,
 *      AwFmAllocationFailure on failure to allocate the suffix array.
 *      AwFmNullPtrError if the indexPtr or databaseSequence were null.
 *      AwFmSuffixArrayCreationFailure on failure in building the suffix array (this error
 *        is caused by internal failure in the libdivsufsort library).
*/
enum AwFmReturnCode awFmCreateFullSuffixArray(const uint8_t *databaseSequence,
  const uint64_t databaseSequenceLength, uint64_t **fullSuffixArrayOut){

  //suffix array is one longer than database sequence, since it needs to store the null terminating $.
  const uint64_t suffixArrayLength = databaseSequenceLength + 1;
  uint64_t *fullSuffixArray = malloc(suffixArrayLength * sizeof(uint64_t));
  if(fullSuffixArray == NULL){
    return AwFmAllocationFailure;
  }
  int64_t *suffixArrayAfterTerminator = (int64_t*)(fullSuffixArray + 1);
  uint64_t divSufSortReturnCode = divsufsort64(databaseSequence, suffixArrayAfterTerminator, databaseSequenceLength);

  //the first position will be the location of the null terminator.
  fullSuffixArray[0] = databaseSequenceLength;

  *fullSuffixArrayOut = fullSuffixArray;

  if(divSufSortReturnCode < 0){
    free(fullSuffixArray);
    return AwFmSuffixArrayCreationFailure;
  }
  else{
    return AwFmSuccess;
  }
}
