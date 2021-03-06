#include "AwFmIndex.h"
#include "AwFmOccurrence.h"
#include "AwFmLetter.h"

#include <stdbool.h>


#define BYTES_PER_AVX2_REGISTER         32



AwFmSimdVec256 awFmMakeNucleotideOccurrenceVector(const struct AwFmNucleotideBlock *restrict const blockPtr,
  const uint8_t letter){
  //load the letter bit vectors
  const AwFmSimdVec256 *restrict const blockVectorPtr = blockPtr->letterBitVectors;
  const AwFmSimdVec256 bit0Vector          = AwFmSimdVecLoad(blockVectorPtr);
  const AwFmSimdVec256 bit1Vector          = AwFmSimdVecLoad(blockVectorPtr + 1);
  const AwFmSimdVec256 bit2Vector          = AwFmSimdVecLoad(blockVectorPtr + 2);

  switch(letter){
    case 0: //Nucleotide A 0b110
      return AwFmSimdVecAnd(bit2Vector, bit1Vector);
    case 1://Nucleotide C 0b101
      return AwFmSimdVecAnd(bit2Vector, bit0Vector);
    case 2://Nucleotide G 0b011
      return AwFmSimdVecAnd(bit1Vector, bit0Vector);
    case 3://Nucletoide T 0b001
      return  AwFmSimdVecAndNot(bit2Vector, AwFmSimdVecAndNot(bit1Vector, bit0Vector));
    case 4: //ambiguity character 'X' 0b010
      return  AwFmSimdVecAndNot(bit2Vector, AwFmSimdVecAndNot(bit0Vector, bit1Vector));
      //0b100 is sentinel, but since you can't search for sentinels, it is not included here.
    default:
    __builtin_unreachable();
  }
}



/*
 * Function:  awFmMakeAminoAcidOccurrenceVector
 * --------------------
 * Computes the vector of characters before the given position equal to the given letter.
 *
 *  Inputs:
 *    blockPtr: Pointer to the AwFmNucleotideBlock in which the query position resides.
 *    localQueryPosition: The local position of the query in the block
 *    letter: letter for which the occurrence request is for.
 *
 *  Returns:
 *   Vector with bits set at every position the given letter was found.
 */
AwFmSimdVec256 awFmMakeAminoAcidOccurrenceVector(const struct AwFmAminoBlock *restrict const blockPtr,
  const uint8_t letter){


  //load the letter bit vectors
  const AwFmSimdVec256 *restrict const blockVectorPtr = blockPtr->letterBitVectors;
  const AwFmSimdVec256 bit0Vector = AwFmSimdVecLoad(blockVectorPtr);
  const AwFmSimdVec256 bit1Vector = AwFmSimdVecLoad(blockVectorPtr + 1);
  const AwFmSimdVec256 bit2Vector = AwFmSimdVecLoad(blockVectorPtr + 2);
  const AwFmSimdVec256 bit3Vector = AwFmSimdVecLoad(blockVectorPtr + 3);
  const AwFmSimdVec256 bit4Vector = AwFmSimdVecLoad(blockVectorPtr + 4);

  switch(__builtin_expect(letter, 0)){
    case 0:   /*A (alanine) encoding 0b01100*/
      return AwFmSimdVecAnd(bit3Vector, AwFmSimdVecAndNot(bit4Vector, bit2Vector));
    case 1:   /*C (cysteine) encoding 0b10111*/
      return  AwFmSimdVecAnd(AwFmSimdVecAndNot(bit3Vector, bit2Vector), AwFmSimdVecAnd(bit1Vector, bit0Vector));
    case 2:   /*D (aspartic acid) encoding 0b00011*/
      return AwFmSimdVecAnd(bit1Vector, AwFmSimdVecAndNot(bit4Vector, bit0Vector));
    case 3:   /*E (Glutamic acid) encoding 0b00110*/
      return AwFmSimdVecAndNot(bit4Vector, AwFmSimdVecAnd(bit2Vector, bit1Vector));
    case 4:   /*F (Phenylalanine) encoding 0b11110*/
      return AwFmSimdVecAnd(AwFmSimdVecAndNot(bit0Vector, bit3Vector), AwFmSimdVecAnd(bit2Vector, bit1Vector));
    case 5:   /*G (Glycine) encoding 0b11010*/
      return AwFmSimdVecAndNot(bit2Vector, AwFmSimdVecAndNot(bit0Vector, bit4Vector));
    case 6:   /*H (Histidine) encoding 0b11011*/
      return AwFmSimdVecAnd(AwFmSimdVecAndNot(bit2Vector, bit3Vector), AwFmSimdVecAnd(bit1Vector, bit0Vector));
    case 7:   /*I (Isoleucine) encoding 0b11001*/
      return AwFmSimdVecAndNot(bit2Vector, AwFmSimdVecAndNot(bit1Vector, bit4Vector));
    case 8:   /*K (Lysine) encoding 0b10101*/
      return AwFmSimdVecAndNot(bit3Vector, AwFmSimdVecAndNot(bit1Vector, bit4Vector));
    case 9:   /*L (Leucine) encoding 0b11100*/
      return AwFmSimdVecAndNot(bit1Vector, AwFmSimdVecAndNot(bit0Vector, bit4Vector));
    case 10:  /*M (Methionine) encoding 0b11101*/
      return AwFmSimdVecAnd(AwFmSimdVecAndNot(bit1Vector, bit3Vector), AwFmSimdVecAnd(bit2Vector, bit0Vector));
    case 11:  /*N (Asparagine) encoding 0b01000*/
      return AwFmSimdVecAndNot(AwFmSimdVecOr(bit0Vector, bit1Vector), AwFmSimdVecAndNot(bit2Vector, bit3Vector));
    case 12:  /*P (Proline) encoding 0b01001*/
      return AwFmSimdVecAnd(bit3Vector, AwFmSimdVecAndNot(bit4Vector, bit0Vector));
    case 13:  /*Q (glutamine) encoding 0b00100*/
      return AwFmSimdVecAndNot(AwFmSimdVecOr(bit3Vector, bit1Vector), AwFmSimdVecAndNot(bit0Vector, bit2Vector));
    case 14:  /*R (Arginine) encoding 0b10011*/
      return AwFmSimdVecAndNot(bit3Vector, AwFmSimdVecAndNot(bit2Vector, bit4Vector));
    case 15:  /*S (Serine) encoding 0b01010*/
      return AwFmSimdVecAnd(bit3Vector, AwFmSimdVecAndNot(bit4Vector, bit1Vector));
    case 16:  /*T (Threonine) encoding 0b00101*/
      return AwFmSimdVecAnd(bit2Vector, AwFmSimdVecAndNot(bit4Vector, bit0Vector));
    case 17:  /*V (Valine) encoding 0b10110*/
      return AwFmSimdVecAndNot(bit3Vector, AwFmSimdVecAndNot(bit0Vector, bit4Vector));
    case 18:  /*W (Tryptophan) encoding 0b00001*/
      return AwFmSimdVecAndNot(AwFmSimdVecOr(bit3Vector, bit2Vector), AwFmSimdVecAndNot(bit1Vector, bit0Vector));
    case 19:  /*Y (Tyrosine) encoding 0b00010*/
      return AwFmSimdVecAndNot(AwFmSimdVecOr(bit0Vector, bit2Vector), AwFmSimdVecAndNot(bit3Vector, bit1Vector));
    case 20: /*ambiguity character Z 0b11111 */
      return AwFmSimdVecAnd(AwFmSimdVecAnd(bit3Vector, bit2Vector),AwFmSimdVecAnd(bit1Vector, bit0Vector));
    //0b00000 is sentinel, but since you can't search for sentinels, it is not included here.
    default: __builtin_unreachable();   //GCC respects this, doesn't check for letters that aren't valid
  }
}



// inline uint16_t awFmVectorPopcount(const AwFmSimdVec256 occurrenceVector, const uint8_t localQueryPosition){
//   return awFmVectorPopcountBuiltin(occurrenceVector, localQueryPosition);
  //deprecated
// const AwFmSimdVec256 lowBitsLookupTable = _mm256_set_epi8( 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0, 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0);
// const AwFmSimdVec256 highBitsLookupTable = _mm256_set_epi8(-5,-4,-4,-3,-4,-3,-3,-2,-4,-3,-3,-2,-3,-2,-2, -1,-5,-4,-4,-3,-4,-3,-3,-2,-4,-3,-3,-2,-3,-2,-2, -1);
// const AwFmSimdVec256 lowNybbleBitmasked = AwFmSimdVecAnd(occurrenceVector, _mm256_set1_epi8(0x0F));
// const AwFmSimdVec256 lowNybbleBitCount = _mm256_shuffle_epi8(lowBitsLookupTable, lowNybbleBitmasked);
// const AwFmSimdVec256 highNybbleBits = _mm256_srli_epi16(occurrenceVector, 4);
// const AwFmSimdVec256 highNybbleBitmasked = AwFmSimdVecAnd(highNybbleBits, _mm256_set1_epi8(0x0F));
// const AwFmSimdVec256 highNybbleNegativeBitCount = _mm256_shuffle_epi8(highBitsLookupTable, highNybbleBitmasked);
// const AwFmSimdVec256 sadCountVector = _mm256_sad_epu8(lowNybbleBitCount, highNybbleNegativeBitCount);
//
// const uint8_t finalSum = 224 -
//   ((_mm256_extract_epi8(sadCountVector, 0)) + (_mm256_extract_epi8(sadCountVector, 8)) +
//   (_mm256_extract_epi8(sadCountVector, 16)) + (_mm256_extract_epi8(sadCountVector, 24)));
// //224 = 256 - 32, where the 32 is from subtracting one extra from each in the highBitsLookupTable
// // keeping all the high bitmask one less than actual forces every add to overflow, making this consistent.
//
// return finalSum;
// }


// inline uint16_t awFmVectorPopcountBuiltin(const AwFmSimdVec256 occurrenceVector, const uint8_t localQueryPosition){
//   uint64_t bitmasks[4] = {0};
//   uint8_t  bitmaskedQuadWordIndex = localQueryPosition / 64;
//
//   for(int8_t i = 0; i < bitmaskedQuadWordIndex; i++){
//       bitmasks[i] = ~0UL;
//   }
//   bitmasks[bitmaskedQuadWordIndex] = ~0UL >> (63 - (localQueryPosition % 64));
//
//   uint16_t popcount = _mm_popcnt_u64(_mm256_extract_epi64(occurrenceVector, 0) & bitmasks[0]);
//    popcount +=        _mm_popcnt_u64(_mm256_extract_epi64(occurrenceVector, 1) & bitmasks[1]);
//    popcount +=        _mm_popcnt_u64(_mm256_extract_epi64(occurrenceVector, 2) & bitmasks[2]);
//    popcount +=        _mm_popcnt_u64(_mm256_extract_epi64(occurrenceVector, 3) & bitmasks[3]);
//
//   return popcount;
// }



  //alternative:
  // uint8_t unalignedBuffer[32 + 31];
  // uint8_t *alignedBuffer = (uint8_t *)((intptr_t)( unalignedBuffer + 31 ) & ~31);
  // _mm256_store_si256((AwFmSimdVec256*)alignedBuffer, vec);
  //
  // uint64_t *bufferAs64_t = (uint64_t*) alignedBuffer;
  // return _mm_popcnt_u64(bufferAs64_t[0]) +
  // _mm_popcnt_u64(bufferAs64_t[1]) +
  // _mm_popcnt_u64(bufferAs64_t[2]) +
  // _mm_popcnt_u64(bufferAs64_t[3]);
// }


inline void awFmBlockPrefetch(const void *restrict const baseBlockListPtr, const uint64_t blockByteWidth,
  const uint64_t nextQueryPosition){

  const uint64_t blockIndex    = awFmGetBlockIndexFromGlobalPosition(nextQueryPosition);
  //make the blockAddress pointer as a uint8_t* to make clean and easy pointer arithmetic when defining cache line boundries.
  const uint8_t *blockAddress  = ((uint8_t*)baseBlockListPtr + (blockIndex * blockByteWidth));

  for(uint_fast16_t prefetchOffset = 0; prefetchOffset < blockByteWidth; prefetchOffset += AW_FM_CACHE_LINE_SIZE_IN_BYTES){
    AwFmSimdPrefetch(blockAddress + prefetchOffset);
  }
}


// inline uint8_t awFmGetLetterAtBwtPosition(const union AwFmBwtBlockList blockList,
//   const enum AwFmAlphabetType alphabet, const uint64_t bwtPosition){
//   const size_t blockIndex      = awFmGetBlockIndexFromGlobalPosition(bwtPosition);
//   const size_t positionInBlock = awFmGetBlockQueryPositionFromGlobalPosition(bwtPosition);
//   const uint8_t byteInBlock    = positionInBlock / 8;
//   const uint8_t bitInBlockByte = positionInBlock / 8;
//   uint8_t letter = 0;
//   uint8_t letterBitWidth = alphabet == AwFmAlphabetNucleotide? 2: 5;
//   const AwFmSimdVec256 *restrict const letterBitVectorPtr = alphabet == AwFmAlphabetNucleotide?
//     blockList.asNucleotide[blockIndex].letterBitVectors:
//     blockList.asAmino[blockIndex].letterBitVectors;
//
//   for(uint8_t letterBit = letterBitWidth; letterBit >= 0; letterBit--){
//     const uint8_t *restrict const blockVectorAsByteArray = (uint8_t*)(letterBitVectorPtr +letterBit);
//     const uint8_t bit = (blockVectorAsByteArray[byteInBlock] >> bitInBlockByte) & 1;
//
//     letter = (letter << 1) | bit;
//   }
//
//   return letter;
// }


/*
 * Function:  awFmGetNucleotideLetterAtBwtPosition
 * --------------------
 * Given a specific position in the nucleotide BWT, returns the letter in the BWT at this position.
 *
 *  Inputs:
 *    blockList: Block to query into. Note, this is usually not the pointer to the bwt,
 *      but a pointer to the block the query will be found in.
 *    bwtPosition: Position of the character to be returned, relative to the start of the block to query.
 *
 *  Returns:
 *    Letter index of the character at the bwtPosition in the specified blockList.
 */
uint8_t awFmGetNucleotideLetterAtBwtPosition(const struct AwFmNucleotideBlock *blockPtr, const uint8_t localPosition){
  const uint8_t byteInBlock    = localPosition / 8;
  const uint8_t bitInBlockByte = localPosition % 8;

  const uint8_t *restrict const letterBytePointer = &((uint8_t*) &blockPtr->letterBitVectors)[byteInBlock];
  const uint8_t letterAsCompressedVector = ((letterBytePointer[0]  >> bitInBlockByte) & 1) |
          ((letterBytePointer[32] >> bitInBlockByte) & 1) << 1 |
          ((letterBytePointer[64] >> bitInBlockByte) & 1) << 2;

  return awFmNucleotideCompressedVectorToLetterIndex(letterAsCompressedVector);
  }


/*
 * Function:  awFmGetLetterAtBwtPosition
 * --------------------
 * Given a specific position in the amino BWT, returns the letter in the BWT at this position.
 *
 *  Inputs:
 *    blockList: Block to query into. Note, this is usually not the pointer to the bwt,
 *      but a pointer to the block the query will be found in.
 *    bwtPosition: Position of the character to be returned, relative to the start of the block to query.
 *
 *  Returns:
 *    Letter index of the character at the bwtPosition in the specified blockList.
 */
uint8_t awFmGetAminoLetterAtBwtPosition(const struct AwFmAminoBlock *blockPtr, const uint8_t localPosition){
  const uint8_t byteInBlock    = localPosition / 8;
  const uint8_t bitInBlockByte = localPosition % 8;

  const uint8_t *restrict const letterBytePointer = &((uint8_t*) &blockPtr->letterBitVectors)[byteInBlock];
  const uint8_t letterAsCompressedVector = ((letterBytePointer[0] >> bitInBlockByte) & 1)       |
    ((letterBytePointer[32]     >> bitInBlockByte) & 1) << 1  |
    ((letterBytePointer[64]     >> bitInBlockByte) & 1) << 2  |
    ((letterBytePointer[96]     >> bitInBlockByte) & 1) << 3  |
    ((letterBytePointer[128]    >> bitInBlockByte) & 1) << 4;

  return awFmAminoAcidCompressedVectorToLetterIndex(letterAsCompressedVector);
}


// //I might only need to mask one quad-word, then I can just skip popcounting the ones after!
// //okay, new strategy: compress to bit vector with AVX2. popcount and accumulate up to mask word,
// //then shift the mask word to zero out as many bits as needed., popcount and add that, then done.
// //if given position 0, one bit should be set., if given 256, all bits should be set.
// /*
//  * Function:  createQueryPositionBitmask
//  * --------------------
//  *  Creates an AVX2 vector that acts as a bitmask to remove positions after the query position,
//  *    and to remove the sentinel character if it exists in this vector.
//  *
//  *  Inputs:
//  *    localQueryPosition: position in this AVX2 vector to query. All bits after
//  *      this position are cleared in the returned vector.
//  *
//  *   Returns:
//  *     Bitmask Vector for the occurrence vector.
//  */
// inline AwFmSimdVec256 createQueryPositionBitmask(const uint8_t localQueryPosition){
//   //make the mask for the quad word where the localQueryPositions shows up in.
//   uint64_t quadWordMask = (1UL << ((uint64_t)(localQueryPosition) +1) % 64) - 1;
//   const uint8_t quadWordIndexForQueryPosition = localQueryPosition / 64;
//   AwFmSimdVec256 maskVector = _mm256_setzero_si256();
//
//   //starting in the quad word where the query position shows up, set the quad word mask,
//   // and then change the quad word mask to all 1s so everything below gets set to 1s as well.
//   // this switch-case intentionally falls through, and is supposed to not have 'break' statements.
//   switch(quadWordIndexForQueryPosition){
//     case 3:
//     maskVector = _mm256_insert_epi64(maskVector, quadWordMask, 3);
//     quadWordMask = -1UL;
//     case 2:
//     maskVector = _mm256_insert_epi64(maskVector, quadWordMask, 2);
//     quadWordMask = -1UL;
//     case 1:
//     maskVector = _mm256_insert_epi64(maskVector, quadWordMask, 1);
//     quadWordMask = -1UL;
//     case 0:
//     maskVector = _mm256_insert_epi64(maskVector, quadWordMask, 0);
//   }
//
//   return maskVector;
// }
