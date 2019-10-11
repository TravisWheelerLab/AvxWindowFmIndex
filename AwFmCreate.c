#include "AwFmCreate.h"
#include "AwFmLetter.h"
#include <stdlib.h>
#include <string.h>


//temp extern while libdivsufsort isn't installed
extern uint64_t divsufsort(const uint8_t *T, uint64_t *SA, uint64_t n);


enum AwFmReturnCode awFmCreateFullSuffixArray(const uint8_t * databaseSequence,
  const uint64_t databaseSequenceLength, uint64_t **fullSuffixArrayOut);

enum AwFmReturnCode awFmCreateBlockList(const uint8_t *restrict const databaseSequence,
  const uint64_t databaseSequenceLength, struct AwFmIndex *restrict const index);



enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex **indexPtr, const uint8_t *restrict const databaseSequence,
  const size_t databaseSequenceLength, const uint16_t suffixArrayCompressionRatio, const char *restrict const fileSrc){

  struct AwFmIndex *index = alignedAllocAwFmIndex();
  if(index == NULL){
    return AwFmAllocationFailure;
  }

  //set the suffix array compression ratio
  index->suffixArrayCompressionRatio = suffixArrayCompressionRatio;

  //set the fileSrc, or return the error code if it failed.
  if(fileSrc != NULL){
    const enum AwFmReturnCode setFileSrcReturnCode = awFmIndexSetFileSrc(index, fileSrc);
    if(__builtin_expect(setFileSrcReturnCode < 0, 0)){
      deallocateFmIndex(index);
      return setFileSrcReturnCode;
    }
  }

  //set the database sequence by reference
  if(databaseSequence == NULL){
    deallocateFmIndex(index);
    return AwFmAllocationFailure;
  }
  index->databaseSequence = databaseSequence;

  //create the full suffix array
  enum AwFmReturnCode suffixArrayCreationReturnCode = awFmCreateFullSuffixArray(databaseSequence,
    databaseSequenceLength, &index->fullSuffixArray);
  if(__builtin_expect(suffixArrayCreationReturnCode < 0, 0)){
    deallocateFmIndex(index);
    return suffixArrayCreationReturnCode;
  }



  enum AwFmReturnCode blockListCreationReturnCode = awFmCreateBlockList(databaseSequence, databaseSequenceLength, index);
  if(__builtin_expect(blockListCreationReturnCode < 0, 0)){
    deallocateFmIndex(index);
    return blockListCreationReturnCode;
  }

  //set the rank prefix sums
  //accumulator starts at 1 because of the null terminator $ being before all other letters
  uint64_t rankPrefixAccumulator = 1;
  const uint64_t lastBlockIndex = index->numBlocks - 1;

  index->rankPrefixSums[0] = rankPrefixAccumulator;
  for(uint_fast8_t i = 0; i < AMINO_CARDINALITY; i++){
    rankPrefixAccumulator += index->blockList[lastBlockIndex].baseOccupancies[i];
    index->rankPrefixSums[i + 1] = rankPrefixAccumulator;
  }

  //TODO; check that this is a valid thing to return right now, I need to error check the function.
  return AwFmSuccess;
}


//if failure, the caller may need to deallocate index.
enum AwFmReturnCode awFmCreateBlockList(const uint8_t *restrict const databaseSequence,
  const uint64_t databaseSequenceLength, struct AwFmIndex *const restrict index){

  if(index == NULL){
    return AwFmNullPtrError;
  }

  const uint64_t suffixArrayLength = databaseSequenceLength + 1;
  //ceiling function to compute the number of blocks needed

  index->numBlocks = (suffixArrayLength + (POSITIONS_PER_FM_BLOCK - 1)) / POSITIONS_PER_FM_BLOCK;
  index->blockList = malloc(index->numBlocks * sizeof(struct AwFmBlock));
  if(index->blockList == NULL){
    return AwFmAllocationFailure;
  }

  //keep an accumulator for the baseOccupancies and current position in suffix array
  //by adding 1 to the length, anything that parses incorrectly won't crash the program,
  //but will accumulate into an unused value
  uint64_t accumulatedOccupancy[AMINO_CARDINALITY + 1] = {0};

  for(uint64_t blockIndex = 0; blockIndex < index->numBlocks; blockIndex++){
    //copy the base occupancy over, and set the bit vectors to all zeros.
    memcpy(index->blockList[blockIndex].baseOccupancies, accumulatedOccupancy, AMINO_CARDINALITY * sizeof(uint64_t));
    memset(index->blockList[blockIndex].letterBitVectors, 0, sizeof(__m256i) * 5);

    //create a uint8_t pointer to the letter bit vectors to more easily set individual bytes.
    void *letterBitVectorsAsRawPtr = &index->blockList[blockIndex].letterBitVectors[0];
    uint64_t *letterBitVectorsAsBytes = letterBitVectorsAsRawPtr;


    //if we're on the first block, start at position 1, since 0 will be the null terminator.
    const uint_fast8_t blockStartingPosition = (blockIndex == 0)? 1: 0;
    for(uint_fast8_t i = blockStartingPosition; i < POSITIONS_PER_FM_BLOCK; i++){

      const uint8_t byteInBlock             = i / 8;
      const uint8_t bitInBlockByte          = i % 8;
      const uint64_t positionInSuffixArray  = (blockIndex * POSITIONS_PER_FM_BLOCK) + i;

      const uint64_t suffixArrayValue       = index->fullSuffixArray[positionInSuffixArray];
      const uint8_t letterAsAscii           = databaseSequence[suffixArrayValue];
      const uint8_t letterAsFrequencyIndex  = awFmAsciiLetterToLetterIndex(letterAsAscii);
      uint8_t letterAsVectorFormat          = awFmAsciiLetterToCompressedVectorFormat(letterAsAscii);

      //accumulate the letter's occupancy.
      accumulatedOccupancy[letterAsFrequencyIndex]++;

      //set the correct bits in the letterBitVectors
      for(uint_fast8_t bitInVectorLetter = 0; bitInVectorLetter < 5; bitInVectorLetter++){
        const uint8_t letterBit = (letterAsVectorFormat >> bitInVectorLetter) & 1;
        const uint8_t byteInBlockVectors = (letterBit * BYTES_PER_AVX2_REGISTER) + byteInBlock;
        //or in the shifted bit.
        letterBitVectorsAsBytes[byteInBlockVectors] |= (letterBit << bitInBlockByte);
      }
    }
  }

  return AwFmSuccess;
}


enum AwFmReturnCode awFmCreateFullSuffixArray(const uint8_t * databaseSequence,
  const uint64_t databaseSequenceLength, uint64_t **fullSuffixArrayOut){

    //suffix array is one longer than database sequence, since it needs to store the null terminating $.
    const uint64_t suffixArrayLength = databaseSequenceLength + 1;
    uint64_t *fullSuffixArray = malloc(suffixArrayLength * sizeof(uint64_t));
    if(fullSuffixArray == NULL){
      return AwFmAllocationFailure;
    }
    uint64_t *suffixArrayAfterTerminator = fullSuffixArray + 1;
    uint64_t divSufSortReturnCode = divsufsort(databaseSequence, suffixArrayAfterTerminator, databaseSequenceLength);

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
