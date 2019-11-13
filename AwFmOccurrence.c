#include "AwFmLetter.h"
#include "AwFmOccurrence.h"
#include "AwFmGlobals.h"
#include <immintrin.h>
#include <stdbool.h>

/*Private Function Prototypes*/
/*Computes the index in the blockList where the given BWT position is found.*/
static inline size_t        getBlockIndexFromGlobalPosition(const uint64_t globalQueryPosition);

/*computes the position where the given BWT position can be found in the relevant AwFmBlock.*/
static inline uint_fast8_t  getBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition);

/*constructs the AVX2 vector representing the positions in the AwFmBlock where the given letter can be found.*/
inline __m256i              createAminoAcidOccurrenceVector(const struct AwFmAminoBlock *restrict const blockList,
  const size_t blockIndex, const uint8_t letter);

/*Creates and applies a mask vector to the occurrence vector to zero out all positions after the query position.*/
static inline __m256i       applyQueryPositionBitmask(const __m256i OccurrenceVector,
  const uint8_t localQueryPosition);

/*performs a popcount operation on the given maskedOccurrenceVector to determine the Occurrences
  of this block before query position.*/
static inline uint_fast8_t  countBitsInMaskedOccurrenceVector(const __m256i maskedOccurrenceVector,
  const __m256i lowBitsLookupTable, const __m256i highBitsLookupTable);


/*
 * Function:  awFmGetAminoAcidOccurrence
 * --------------------
 * Calculates the BWT's occurrency at the given position for the given letter.
 *    Occurrency is calculated by creating an AVX2 vector representing the positions in the block that contain the specified letter.
 *    Then, a position mask is applied to make sure no positions after the given query position are included.
 *
 *  Inputs:
 *    blockList:      Pointer to the array of blocks that make up the BWT array.
 *    queryPosition:  Position in the BWT that the occurrency call will calculate.
 *    letter:         Frequency-indexed letter that's occurrency is being queried.
 */
uint64_t awFmGetAminoAcidOccurrence(const struct AwFmAminoBlock *restrict const blockList, const size_t queryPosition, const uint8_t letter){
  const uint_fast8_t queryLocalPosition = getBlockQueryPositionFromGlobalPosition(queryPosition);
  const size_t blockIndex               = getBlockIndexFromGlobalPosition(queryPosition);
  const __m256i occurrenceVector         = createAminoAcidOccurrenceVector(blockList, blockIndex, letter);
  const __m256i maskedOccVector         = applyQueryPositionBitmask(occurrenceVector, queryLocalPosition);
  //make lookup tables for computing the occurrence popcount.
  //high bitmask contains negative numbers because the SAD intrinsic we use to horizontal add actual subtracts. so subtracting negative == adding!
  const __m256i lowBitsLookupTable      = _mm256_setr_epi8( 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0, 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0);
  const __m256i highBitsLookupTable     = _mm256_setr_epi8(-4,-3,-3,-2,-3,-2,-2,-1,-3,-2,-2,-1,-2,-1,-1, 0,-4,-3,-3,-2,-3,-2,-2,-1,-3,-2,-2,-1,-2,-1,-1, 0);
  const uint_fast8_t blockOccurrency     = countBitsInMaskedOccurrenceVector(maskedOccVector, lowBitsLookupTable, highBitsLookupTable);
  const uint64_t blockBaseOccurrency     = blockList[blockIndex].baseOccurrences[letter];

  const uint64_t occurrency              = blockOccurrency + blockBaseOccurrency;
  return occurrency;
}


/*
 * Function:  awFmAminoAcidDataPrefetch
 * --------------------
 *  Manually request cache prefetching for the amino acid data to be read for the occurrence call.
 *  Prefetching is done using either the _mm_prefetch intel intrinsic, or the __builtin_prefetch GCC extension.
 *    The correct block to prefetch is computed from the position.
 *
 *  Inputs:
 *    blockList:      Pointer to the array of blocks that make up the BWT array.
 *    queryPosition:  Position in the BWT that the occurrence call will calculate.
 *
 */
void awFmAminoAcidDataPrefetch(const struct AwFmAminoBlock *restrict const blockList, const uint64_t queryPosition){
  const uint64_t blockIndex    = getBlockIndexFromGlobalPosition(queryPosition);
  //make the blockAddress pointer as a uint8_t* to make clean and easy pointer arithmetic when defining cache line boundries.
  const uint8_t *blockAddress  = (uint8_t*)(blockList + blockIndex);

  for(uint_fast16_t prefetchOffset = 0; prefetchOffset < sizeof(struct AwFmAminoBlock); prefetchOffset += CACHE_LINE_SIZE_IN_BYTES){
    _mm_prefetch(blockAddress + prefetchOffset, _MM_HINT_T2);
  }
}



/*
 * Function:  awFmBackstepBwtPosition
 * --------------------
 *  Backtraces the the given position to find the BWT position of the letter that
 *    comes immediately before it in the database sequence.
 *
 *   Inputs:
 *     index:          Pointer to the AwFmIndex to query.
 *     bwtPosition:    Position of the requsted letter in the BWT.
 *
 *   Returns:
 *     Position of the letter immediately before it in the database sequence.
 */
inline size_t awFmBackstepAminoAcidBwtPosition(const struct AwFmAminoBlock *restrict const blockList,
  const uint64_t *restrict const rankPrefixSums, const uint16_t suffixArrayCompressionRatio, const uint64_t bwtPosition){

  const uint8_t frequencyIndexLetter  = awFmGetAminoAcidLetterAtBwtPosition(blockList, bwtPosition);
  const uint64_t occurrence           = awFmGetAminoAcidOccurrence(blockList, bwtPosition, frequencyIndexLetter);
  const uint64_t backtraceBwtPosition = occurrence + rankPrefixSums[frequencyIndexLetter];

  //prefetch the next occurrency call if needed.
  const bool backtraceNotFinished = backtraceBwtPosition % suffixArrayCompressionRatio;

  //using __builtin_expect helps the branch prediction expect that it's likely to need to prefetch more.
  if(__builtin_expect(backtraceNotFinished, 1)){
    awFmAminoAcidDataPrefetch(blockList, backtraceBwtPosition);
  }

  return backtraceBwtPosition;
}



 uint8_t awFmGetAminoAcidLetterAtBwtPosition(const struct AwFmAminoBlock *restrict const blockList, const uint64_t bwtPosition){
   const size_t blockIndex      = getBlockIndexFromGlobalPosition(bwtPosition);
   const size_t positionInBlock = getBlockQueryPositionFromGlobalPosition(bwtPosition);
   const uint8_t byteInBlock    = positionInBlock / 8;
   const uint8_t bitInBlockByte = positionInBlock / 8;

   uint8_t letter = 0;
   for(uint8_t letterBit = 5; letterBit >= 0; letterBit--){
     const uint8_t *restrict const blockVectorAsByteArray = (uint8_t*)(blockList[blockIndex].letterBitVectors + letterBit);
     const uint8_t bit = (blockVectorAsByteArray[byteInBlock] >> bitInBlockByte) & 1;

     letter = (letter << 1) | bit;
   }

   return letter;
}

/*
 * Function:  getBlockIndexFromGlobalPosition
 * --------------------
 *  Computes the block index, given the full BWT query position.
 *  Inputs:
 *    globalQueryPosition: Position in the BWT that the occurrency function is requesting
 *   Returns:
 *     Index of the block where the given query position resides.
 */
inline size_t getBlockIndexFromGlobalPosition(const size_t globalQueryPosition){
  return globalQueryPosition / POSITIONS_PER_FM_BLOCK;
}


/*
 * Function:  getBlockQueryPositionFromGlobalPosition
 * --------------------
 *  Computes bit position inside the block that represents the given full BWT position.
 *  Inputs:
 *    globalQueryPosition: Position in the BWT that the occurrency function is requesting
 *   Returns:
 *     Bit position into the block's AVX2 vectors where the query position lands.
 */
inline uint_fast8_t getBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition){
  return globalQueryPosition % POSITIONS_PER_FM_BLOCK;
}


/*
 * Function:  createAminoAcidOccurrenceVector
 * --------------------
 *  Loads the block's letter bit vectors into YMM registers, and uses a switch table to compress them
 *  into a vector where each bit represents whether the letter was found at that position.
 *
 *  Inputs:
 *    blockList: ptr to the array of BWT blocks.
 *    blockIndex: Index of the block we will load.
 *    letter: letter to generate an occurrence vector of.
 *
 *   Returns:
 *    256-bit AVX2 vector with bits set at every position where the given letter was found.
 */
inline __m256i createAminoAcidOccurrenceVector(const struct AwFmAminoBlock *restrict const blockList, const size_t blockIndex, const uint8_t letter){
  //load the letter bit vectors
  const __m256i *restrict const blockVectorPtr = blockList[blockIndex].letterBitVectors;
  const __m256i bit0Vector = _mm256_load_si256(blockVectorPtr);
  const __m256i bit1Vector = _mm256_load_si256(blockVectorPtr + 1);
  const __m256i bit2Vector = _mm256_load_si256(blockVectorPtr + 2);
  const __m256i bit3Vector = _mm256_load_si256(blockVectorPtr + 3);
  const __m256i bit4Vector = _mm256_load_si256(blockVectorPtr + 4);

  //GCC and ICC-complient hint to the compiler that the first letter is the most common,
  //so if you're going to do branch target prediction, predict Leucine.
  //Also, let the compiler know that the last 8 letters are unlikely.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
  __builtin_expect(letter == 0, 1);     //0 is very likely, expect it
  __builtin_expect(letter >= 12, 0);    //12-19 are very unlikely, don't expect them.
#pragma GCC diagnostic pop

  switch(letter){
    case 0:   /*L (leucine). Group One, encoding 0b11100*/
      return  _mm256_andnot_si256(bit1Vector, _mm256_andnot_si256(bit0Vector, bit4Vector));
    case 1:   /*A (alanine). Group One, encoding 0b11010*/
      return _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit0Vector, bit4Vector));
    case 2:   /*G (glycine). Group One, encoding 0b10110*/
      return _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit0Vector, bit4Vector));
    case 3:   /*V (valine). Group One, encoding 0b11001*/
      return _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit1Vector, bit4Vector));
    case 4:   /*E (glutamic acid). Group One, encoding 0b10101*/
      return _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit1Vector, bit4Vector));
    case 5:   /*S (serine). Group One, encoding 0b10011*/
      return _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit2Vector, bit4Vector));
    case 6:   /*I (isoleucine). Group One, encoding 0b00011*/
      return _mm256_and_si256(bit1Vector, _mm256_andnot_si256(bit4Vector, bit0Vector));
    case 7:   /*K (lysine). Group One, encoding 0b00101*/
      return _mm256_and_si256(bit2Vector, _mm256_andnot_si256(bit4Vector, bit0Vector));
    case 8:   /*R (arginine). Group One, encoding 0b01001*/
      return _mm256_and_si256(bit3Vector, _mm256_andnot_si256(bit4Vector, bit0Vector));
    case 9:   /*D (aspartic acid). Group One, encoding 0b00110*/
      return _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit4Vector, bit1Vector));
    case 10:  /*T (threonine). Group One, encoding 0b01010*/
      return _mm256_and_si256(bit3Vector, _mm256_andnot_si256(bit4Vector, bit1Vector));
    case 11:  /*P (proline). Group One, encoding 0b01100*/
      return _mm256_and_si256(bit3Vector, _mm256_andnot_si256(bit4Vector, bit2Vector));
    case 12:  /*N (asparagine). Group Two, encoding 0b11110*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_andnot_si256(bit0Vector, bit1Vector))));
    case 13:  /*Q (glutamine). Group Two, encoding 0b11101*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_andnot_si256(bit1Vector, bit0Vector))));
    case 14:  /*F (phenylalanine). Group Two, encoding 0b11011*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit1Vector, _mm256_andnot_si256(bit2Vector, bit0Vector))));
    case 15:  /*Y (tyrosine), Group Two, encoding 0b10111*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit2Vector, _mm256_and_si256(bit1Vector, _mm256_andnot_si256(bit3Vector, bit0Vector))));
    case 16:  /*M (methionine). Group Two, encoding 0b00001*/
      return _mm256_andnot_si256(bit4Vector, _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit1Vector, bit0Vector))));
    case 17:  /*H (histidine). Group Two, encoding 0b00010*/
      return _mm256_andnot_si256(bit4Vector, _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit0Vector, bit1Vector))));
    case 18:  /*C (cysteine). Group Two, encoding 0b00100*/
      return _mm256_andnot_si256(bit4Vector, _mm256_andnot_si256(bit3Vector, _mm256_and_si256(bit1Vector, _mm256_andnot_si256(bit0Vector, bit2Vector))));
    case 19:  /*W (tryptophan), Group Two, encoding 0b01000*/
      return _mm256_andnot_si256(bit4Vector, _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit1Vector, _mm256_andnot_si256(bit3Vector, bit0Vector))));
    default: __builtin_unreachable();   //GCC respects this, doesn't check for letters that aren't valid
    }
}



/*
 * Function:  createAminoAcidOccurrenceGteVector
 * --------------------
 *  Loads the block's letter bit vectors into YMM registers, and uses a switch table to compress them
 *  into a vector where each bit represents positions where the character was
 *  greater than or equal to the given character.
 *
 *  Inputs:
 *    blockList: ptr to the array of BWT blocks.
 *    blockIndex: Index of the block we will load.
 *    letter: letter to generate an occurrence vector of.
 *
 *   Returns:
 *    256-bit AVX2 vector with bits set at every position that represents a character greater than or equal to letter.
 */
inline __m256i createAminoAcidOccurrenceLessThanVector(const struct AwFmAminoBlock *restrict const blockList, const size_t blockIndex, const uint8_t letter){
  //load the letter bit vectors
  const __m256i *restrict const blockVectorPtr = blockList[blockIndex].letterBitVectors;
  const __m256i bit4Vector = _mm256_load_si256(blockVectorPtr + 4);
  const __m256i bit3Vector = _mm256_load_si256(blockVectorPtr + 3);
  const __m256i bit2Vector = _mm256_load_si256(blockVectorPtr + 2);
  const __m256i bit1Vector = _mm256_load_si256(blockVectorPtr + 1);
  const __m256i bit0Vector = _mm256_load_si256(blockVectorPtr);

  //GCC and ICC-complient hint to the compiler that the first letter is the most common,
  //so if you're going to do branch target prediction, predict Leucine.
  //Also, let the compiler know that the last 8 letters are unlikely.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
  __builtin_expect(letter == 0, 1);     //0 is very likely, expect it
  __builtin_expect(letter >= 12, 0);    //12-19 are very unlikely, don't expect them.
#pragma GCC diagnostic pop

  //each bit, AND if bit is set, or if bit is clear, down to the last set bit.
  switch(letter){
    case 0:   /*L (leucine). Group One, encoding 0b11100*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, bit2Vector));
    case 1:   /*A (alanine). Group One, encoding 0b11010*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, bit1Vector)));
    case 2:   /*G (glycine). Group One, encoding 0b10110*/
      return _mm256_and_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, bit1Vector)));
    case 3:   /*V (valine). Group One, encoding 0b11001*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
    case 4:   /*E (glutamic acid). Group One, encoding 0b10101*/
      return _mm256_and_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
    case 5:   /*S (serine). Group One, encoding 0b10011*/
      return _mm256_and_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_and_si256(bit1Vector, bit0Vector))));
    case 6:   /*I (isoleucine). Group One, encoding 0b00011*/
      return _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_and_si256(bit1Vector, bit0Vector))));
    case 7:   /*K (lysine). Group One, encoding 0b00101*/
      return _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
    case 8:   /*R (arginine). Group One, encoding 0b01001*/
      return _mm256_or_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
    case 9:   /*D (aspartic acid). Group One, encoding 0b00110*/
      return _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, bit1Vector)));
    case 10:  /*T (threonine). Group One, encoding 0b01010*/
      return _mm256_or_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, bit1Vector)));
    case 11:  /*P (proline). Group One, encoding 0b01100*/
      return _mm256_or_si256(bit4Vector, _mm256_and_si256(bit3Vector, bit2Vector));
    case 12:  /*N (asparagine). Group Two, encoding 0b11110*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit2Vector, bit1Vector)));
    case 13:  /*Q (glutamine). Group Two, encoding 0b11101*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
    case 14:  /*F (phenylalanine). Group Two, encoding 0b11011*/
      return _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_and_si256(bit1Vector, bit0Vector))));
    case 15:  /*Y (tyrosine), Group Two, encoding 0b10111*/
      return _mm256_and_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_and_si256(bit1Vector, bit0Vector))));
    case 16:  /*M (methionine). Group Two, encoding 0b00001*/
      return _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
    case 17:  /*H (histidine). Group Two, encoding 0b00010*/
      return _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_or_si256(bit2Vector, bit1Vector)));
    case 18:  /*C (cysteine). Group Two, encoding 0b00100*/
      return _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, bit2Vector));
    case 19:  /*W (tryptophan), Group Two, encoding 0b01000*/
      return _mm256_or_si256(bit4Vector, bit3Vector);
    default: __builtin_unreachable();   //GCC respects this, doesn't check for letters that aren't valid
    }
}


/*
 * Function:  applyQueryPositionBitmask
 * --------------------
 *  creates a bitmask based on the block's local query position, and sets
 *  all positions after the query position to zero in the occurrence vector.
 *
 *  Inputs:
 *    occurrenceVector: AVX2 vector generated from the createAminoAcidOccurrenceVector function.
 *    localQueryPosition: position in the block queried in the occurrence function.
 *
 *   Returns:
 *    256-bit vector where only bits before this position are allowed to be set (all bits >= query position are cleared).
 */
inline __m256i applyQueryPositionBitmask(const __m256i occurrenceVector, const uint8_t localQueryPosition){
  const uint8_t numBytesToPreserveAllOnes = localQueryPosition / 8;
  const uint8_t lastQueryByteBitmask      = (1 << (localQueryPosition % 8)) - 1;

  uint8_t bitmaskArray[BYTES_PER_AVX2_REGISTER] = {0};
  for(uint_fast8_t i = 0; i < numBytesToPreserveAllOnes; i++){
    bitmaskArray[i] = 0xFF;
  }
  bitmaskArray[numBytesToPreserveAllOnes] = lastQueryByteBitmask;
  for( uint_fast8_t i = numBytesToPreserveAllOnes; i < BYTES_PER_AVX2_REGISTER; i++){
    bitmaskArray[i] = 0;
  }

  //load and apply the bitmask
  __m256i bitmask = _mm256_lddqu_si256((__m256i*)bitmaskArray);
  return _mm256_and_si256(occurrenceVector, bitmask);
}


/*
 * Function:  countBitsInMaskedOccurrenceVector
 * --------------------
 *  Performs a popcount on the given AVX2 masked occurrence vector.
 *  Loosely based on technique in Mula et.al, 2016
 *
 *  Significant performance gains are realized over Mula by keeping the MSb of each byte unpopulated,
 *    and by using two lookup table vectors to allow the SAD intrinsic to add the nybble counts while it horizontal adds.
 *
 *  Inputs:
 *    maskedOccurrenceVector: Position masked AVX2 popcount vector generated from the applyQueryPositionBitmask function.
 *    lowBitsLookupTable: AVX2 register where each byte contains the number of bits set in the integer that shares the lower nybble with it's index.
 *    highBitsLookupTable: same as lowBitsLookupTable, but negated, so the SAD intrinsic adds them (SAD finds the absolute difference,
 *      so subtracting a negative adds them together).
 *
 *   Returns:
 *    Count of how many bits were set in the maskedOccurrenceVector.
 */
 //NOTE: check high bit count to make sure 1 doesn't get shifted into bit7 of the bytes.
inline uint_fast8_t countBitsInMaskedOccurrenceVector(const __m256i maskedOccurrenceVector, const __m256i lowBitsLookupTable, const __m256i highBitsLookupTable){

  const __m256i lowNybbleBitmasked          = _mm256_and_si256(maskedOccurrenceVector, _mm256_set1_epi8(0x0F));
  const __m256i lowNybbleBitCount           = _mm256_shuffle_epi8(lowBitsLookupTable, lowNybbleBitmasked);
  const __m256i highNybbleBits              = _mm256_srli_si256(maskedOccurrenceVector, 4);
  const __m256i highNybbleBitmasked         = _mm256_and_si256(highNybbleBits, _mm256_set1_epi8(0x0F));
  const __m256i highNybbleNegativeBitCount  = _mm256_shuffle_epi8(highBitsLookupTable, highNybbleBitmasked);
  const __m256i sadCountVector              = _mm256_sad_epu8(lowNybbleBitCount, highNybbleNegativeBitCount);

  //shift and add, placing the final two 16-bit sums in the least significant bits of each 128-bit lane.
  const __m256i laneVectorSums  = _mm256_add_epi16(_mm256_slli_si256 (sadCountVector, 8), sadCountVector);
  const uint16_t finalSum       = _mm256_extract_epi16(laneVectorSums, 0) + _mm256_extract_epi16(laneVectorSums, 8);
  return finalSum;
}
