#include "AwFmLetter.h"
#include "AwFmOccupancy.h"
#include "AwFmGlobals.h"
#include <immintrin.h>
#include <stdbool.h>

/*Private Function Prototypes*/
/*Computes the index in the blockList where the given BWT position is found.*/
static inline size_t        getBlockIndexFromGlobalPosition(const uint64_t globalQueryPosition);

/*computes the position where the given BWT position can be found in the relevant AwFmBlock.*/
static inline uint_fast8_t  getBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition);

/*constructs the AVX2 vector representing the positions in the AwFmBlock where the given letter can be found.*/
static inline __m256i       createBlockOccupancyVector(const struct AwFmBlock *const restrict blockList,
  const size_t blockIndex, const uint8_t letter);

/*Creates and applies a mask vector to the occupancy vector to zero out all positions after the query position.*/
static inline __m256i       applyQueryPositionBitmask(const __m256i occupancyVector,
  const uint8_t localQueryPosition);

/*performs a popcount operation on the given maskedOccupancyVector to determine the occupancy
  of this block before query position.*/
static inline uint_fast8_t  countBitsInMaskedOccupancyVector(const __m256i maskedOccupancyVector,
  const __m256i lowBitsLookupTable, const __m256i highBitsLookupTable);


/*
 * Function:  awFmGetOccupancy
 * --------------------
 * Calculates the BWT's occupancy at the given position for the given letter.
 *    Occupancy is calculated by creating an AVX2 vector representing the positions in the block that contain the specified letter.
 *    Then, a position mask is applied to make sure no positions after the given query position are included.
 *
 *  Inputs:
 *    blockList:      Pointer to the array of blocks that make up the BWT array.
 *    queryPosition:  Position in the BWT that the occupancy call will calculate.
 *    letter:         Frequency-indexed letter that's occupancy is being queried.
 */
uint64_t awFmGetOccupancy(const struct AwFmIndex *const restrict index, const size_t queryPosition, const uint8_t letter){
  const uint_fast8_t queryLocalPosition = getBlockQueryPositionFromGlobalPosition(queryPosition);
  const size_t blockIndex               = getBlockIndexFromGlobalPosition(queryPosition);
  const __m256i occupancyVector         = createBlockOccupancyVector(index->blockList, blockIndex, letter);
  const __m256i maskedOccVector         = applyQueryPositionBitmask(occupancyVector, queryLocalPosition);
  //make lookup tables for computing the occupancy popcount.
  //high bits is negative because the SAD intrinsic we use to horizontal add actual subtracts. so subtracting negative == adding!
  const __m256i lowBitsLookupTable      = _mm256_setr_epi8( 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0, 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0);
  const __m256i highBitsLookupTable     = _mm256_setr_epi8(-4,-3,-3,-2,-3,-2,-2,-1,-3,-2,-2,-1,-2,-1,-1, 0,-4,-3,-3,-2,-3,-2,-2,-1,-3,-2,-2,-1,-2,-1,-1, 0);
  const uint_fast8_t blockOccupancy     = countBitsInMaskedOccupancyVector(maskedOccVector, lowBitsLookupTable, highBitsLookupTable);
  const uint64_t blockBaseOccupancy     = index->blockList[blockIndex].baseOccupancies[letter];

  const uint64_t occupancy              = blockOccupancy + blockBaseOccupancy;
  return occupancy;
}


/*
 * Function:  awFmOccupancyDataPrefetch
 * --------------------
 *  Manually request cache prefetching for the data to be read for the occupancy call.
 *  Prefetching is done using either the _mm_prefetch intel intrinsic, or the __builtin_prefetch GCC extension.
 *    The correct block to prefetch is computed from the position.
 *
 *  Inputs:
 *    blockList:      Pointer to the array of blocks that make up the BWT array.
 *    queryPosition:  Position in the BWT that the occupancy call will calculate.
 *
 */
 void awFmOccupancyDataPrefetch(const struct AwFmIndex *restrict const index, const uint64_t queryPosition){
   for(uint_fast16_t prefetchOffset = 0; prefetchOffset < sizeof(struct AwFmBlock); prefetchOffset += CACHE_LINE_SIZE_IN_BYTES){
     #ifdef AVX_VECTOR_PREFETCH
        //make the blockAddress pointer as a uint8_t* to make clean and easy pointer arithmetic when defining cache line boundries.
        const uint64_t blockIndex    = getBlockIndexFromGlobalPosition(queryPosition);
        const uint8_t *blockAddress  = (uint8_t*)(index->blockList + blockIndex);
       _mm_prefetch(blockAddress + prefetchOffset, _MM_HINT_NTA);
     #else
     #ifdef AVX_VECTOR_GCC_BUILTIN_PREFETCH
       //make the blockAddress pointer as a uint8_t* to make clean and easy pointer arithmetic when defining cache line boundries.
       const uint64_t blockIndex    = getBlockIndexFromGlobalPosition(queryPosition);
       const uint8_t *blockAddress  = (uint8_t*)(index->blockList + blockIndex);
       const int PREFETCH_HINT_READ_ONLY             = 0;
       const int PREFETCH_HINT_NO_TEMPORAL_LOCALITY  = 0;
       __builtin_prefetch(blockAddress + prefetchOffset, PREFETCH_HINT_READ_ONLY, PREFETCH_HINT_NO_TEMPORAL_LOCALITY);
     #endif
     #endif
   }
 }


 /*
  * Function:  awFmGetLetterAtBwtPosition
  * --------------------
  *  Reads the AwFmBlock containing this position, and returns the frequency-index
  *   format letter stored at this position.
  *
  *   Inputs:
  *     index:          Pointer to the AwFmIndex to query.
  *     bwtPosition:    Position of the requsted letter in the BWT.
  *
  *   Returns:
  *     Frequency-indexed format letter that can be used to backtrace from this position.
  */
uint_fast8_t awFmGetLetterAtBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition){
  const uint64_t      blockIndex          = getBlockIndexFromGlobalPosition(bwtPosition);
  const uint_fast8_t  blockQueryPosition  = getBlockQueryPositionFromGlobalPosition(bwtPosition);
  const uint8_t       byteInAvxVector     = blockQueryPosition / 7;
  const uint_fast8_t  bitInVectorByte     = blockQueryPosition % 7;

  const __m256i *restrict const blockVectorPtr = index->blockList[blockIndex].letterBitVectors;

  uint8_t assembledVectorLetter = 0;
  for(uint_fast8_t i = 0; i < 5; i++){
    const uint8_t *ptrToFirstByteInVector = (uint8_t*)(blockVectorPtr + i);
    const uint8_t byteContainingLetterBit = *(ptrToFirstByteInVector + byteInAvxVector);
    const uint_fast8_t letterBitInVector = ((byteContainingLetterBit >> bitInVectorByte) & 1);
    assembledVectorLetter = (assembledVectorLetter << 1) | letterBitInVector;
  }

  return awFmCompressedVectorLetterToLetterIndex(assembledVectorLetter);
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
inline size_t awFmBackstepBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition){
  const uint8_t frequencyIndexLetter = awFmGetLetterAtBwtPosition(index, bwtPosition);
  const uint64_t occupancy = awFmGetOccupancy(index, bwtPosition, frequencyIndexLetter);
  const uint64_t backtraceBwtPosition = occupancy + index->rankPrefixSums[frequencyIndexLetter];

  //prefetch the next occupancy call if needed.
  const bool backtraceNotFinished = backtraceBwtPosition % index->suffixArrayCompressionRatio;

  //using __builtin_expect helps the branch prediction expect that it's likely to need to prefetch more.
  if(__builtin_expect(backtraceNotFinished, 1)){
    awFmOccupancyDataPrefetch(index, backtraceBwtPosition);
  }

  return backtraceBwtPosition;
}


/*
 * Function:  getBlockIndexFromGlobalPosition
 * --------------------
 *  Computes the block index, given the full BWT query position.
 *
 *  Inputs:
 *    globalQueryPosition: Position in the BWT that the occupancy function is requesting
 *
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
 *
 *  Inputs:
 *    globalQueryPosition: Position in the BWT that the occupancy function is requesting
 *
 *   Returns:
 *     Bit position into the block's AVX2 vectors where the query position lands.
 */
inline uint_fast8_t getBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition){
  return globalQueryPosition % POSITIONS_PER_FM_BLOCK;
}


/*
 * Function:  createBlockOccupancyVector
 * --------------------
 *  Loads the block's letter bit vectors into YMM registers, and uses a switch table to compress them
 *  into a vector where each bit represents whether the letter was found at that position.
 *
 *  Inputs:
 *    blockList: ptr to the array of BWT blocks.
 *    blockIndex: Index of the block we will load.
 *    letter: letter to generate an occupancy vector of.
 *
 *   Returns:
 *    256-bit AVX2 vector with bits set at every position where the given letter was found.
 */
inline __m256i createBlockOccupancyVector(const struct AwFmBlock *restrict const blockList, const size_t blockIndex, const uint8_t letter){
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
 * Function:  applyQueryPositionBitmask
 * --------------------
 *  creates a bitmask based on the block's local query position, and sets
 *  all positions after the query position to zero in the occupancy vector.
 *
 *  Inputs:
 *    occupancyVector: AVX2 vector generated from the createBlockOccupancyVector function.
 *    localQueryPosition: position in the block queried in the occupancy function.
 *
 *   Returns:
 *    256-bit vector where only bits before this position are allowed to be set (all bits >= query position are cleared).
 */
inline __m256i applyQueryPositionBitmask(const __m256i occupancyVector, const uint8_t localQueryPosition){
  const uint8_t numBytesToPreserveAllOnes = localQueryPosition / POSITIONS_PER_BYTE_IN_FM_BLOCK;
  const uint8_t lastQueryByteBitmask      = (1 << (localQueryPosition % POSITIONS_PER_BYTE_IN_FM_BLOCK)) - 1;


  uint8_t bitmaskArray[BYTES_PER_AVX2_REGISTER* 2];
  for(size_t i = 0; i < BYTES_PER_AVX2_REGISTER; i+= sizeof(uint64_t)){
    ((uint64_t*) bitmaskArray)[i] = ~0ULL;
  }

  //set the bits in the bite where the query position lands.
  bitmaskArray[BYTES_PER_AVX2_REGISTER] = lastQueryByteBitmask;
  __m256i * ptrToBitmask = (__m256i*)(bitmaskArray + POSITIONS_PER_BYTE_IN_FM_BLOCK - numBytesToPreserveAllOnes);
  return _mm256_lddqu_si256(ptrToBitmask);

}


/*
 * Function:  countBitsInMaskedOccupancyVector
 * --------------------
 *  Performs a popcount on the given AVX2 masked occupancy vector.
 *  Loosely based on technique in Mula et.al, 2016
 *
 *  Significant performance gains are realized over Mula by keeping the MSb of each byte unpopulated,
 *    and by using two lookup table vectors to allow the SAD intrinsic to add the nybble counts while it horizontal adds.
 *
 *  Inputs:
 *    maskedOccupancyVector: Position masked AVX2 popcount vector generated from the applyQueryPositionBitmask function.
 *    lowBitsLookupTable: AVX2 register where each byte contains the number of bits set in the integer that shares the lower nybble with it's index.
 *    highBitsLookupTable: same as lowBitsLookupTable, but negated, so the SAD intrinsic adds them (SAD finds the absolute difference,
 *      so subtracting a negative adds them together).
 *
 *   Returns:
 *    Count of how many bits were set in the maskedOccupancyVector.
 */
inline uint_fast8_t countBitsInMaskedOccupancyVector(const __m256i maskedOccupancyVector, const __m256i lowBitsLookupTable, const __m256i highBitsLookupTable){
  const __m256i lowNybbleBitCount           = _mm256_shuffle_epi8(lowBitsLookupTable, maskedOccupancyVector);
  const __m256i highNybbleBits              = _mm256_srli_si256(maskedOccupancyVector, 4);
  const __m256i highNybbleNegativeBitCount  = _mm256_shuffle_epi8(highBitsLookupTable, highNybbleBits);
  const __m256i sadCountVector              = _mm256_sad_epu8(lowNybbleBitCount, highNybbleNegativeBitCount);

  //shift and add, placing the final two 16-bit sums in the least significant bits of each 128-bit lane.
  const __m256i laneVectorSums  = _mm256_add_epi16(_mm256_slli_si256 (sadCountVector, 8), sadCountVector);
  const uint16_t finalSum       = _mm256_extract_epi16(laneVectorSums, 0) + _mm256_extract_epi16(laneVectorSums, 8);
  return finalSum;
}
