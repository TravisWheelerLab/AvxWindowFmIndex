#include "AwFmOccurrence.h"
#include "AwFmIndex.h"
#include "AwFmLetter.h"
#include "AwFmGlobals.h"

#include <immintrin.h>
#include <stdbool.h>


#define CACHE_LINE_SIZE_IN_BYTES        64
#define BYTES_PER_AVX2_REGISTER         32

/*Private Function Prototypes*/
/*Computes the index in the blockList where the given BWT position is found.*/
inline size_t getBlockIndexFromGlobalPosition(const uint64_t globalQueryPosition);
/*computes the position where the given BWT position can be found in the relevant AwFmBlock.*/
inline uint_fast8_t getBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition);
/*creates a bitmask vector for masking out characters that are either after the query location, or containing the sentinel character.*/
inline __m256i createQueryPositionBitmask(const uint8_t localQueryPosition, bool containsSentinelCharacter, uint8_t sentinelCharacterPosition);


inline void awFmMakeNucleotideOccurrenceVectorPair(const struct AwFmNucleotideBlock *restrict const blockPtr,
  const uint64_t queryPosition, const uint8_t letter, const uint64_t sentinelCharacterPosition,
  struct AwFmOccurrenceVectorPair *vectorPair){

  //create the mask vector to mask away bits after the query position and the sentinel character (if present)
  uint64_t blockIndex                   = getBlockIndexFromGlobalPosition(queryPosition);
  uint64_t sentinelCharacterBlockIndex  = getBlockIndexFromGlobalPosition(sentinelCharacterPosition);
  uint8_t localQueryPosition            = getBlockQueryPositionFromGlobalPosition(queryPosition);
  uint8_t sentinelLocalPosition         = getBlockQueryPositionFromGlobalPosition(sentinelCharacterPosition);
  bool blockContainsSentinel            = sentinelCharacterBlockIndex == blockIndex;
  __m256i bitmask = createQueryPositionBitmask(localQueryPosition, blockContainsSentinel, sentinelLocalPosition);

  //load the letter bit vectors
  const __m256i *restrict const blockVectorPtr = blockPtr->letterBitVectors;
  const __m256i bit1Vector = _mm256_load_si256(blockVectorPtr + 1);
  const __m256i bit0Vector = _mm256_load_si256(blockVectorPtr);

  //create the occurrence vectors
  switch(letter){
    case 0: //Nucleotide A
      vectorPair->occurrenceVector     = _mm256_xor_si256(_mm256_or_si256(bit1Vector, bit0Vector), _mm256_set1_epi8(0xFF));
      vectorPair->occurrenceGteVector  = _mm256_set1_epi8(0xFF);
    case 1://Nucleotide C
      vectorPair->occurrenceVector     = _mm256_andnot_si256(bit1Vector, bit0Vector);
      vectorPair->occurrenceGteVector  = _mm256_or_si256(bit1Vector, bit0Vector);
    case 2://Nucleotide G
      vectorPair->occurrenceVector     = _mm256_andnot_si256(bit0Vector, bit1Vector);
      vectorPair->occurrenceGteVector  = bit1Vector;
    case 3://Nucletoide T
      vectorPair->occurrenceVector     = _mm256_and_si256(bit1Vector, bit0Vector);
      vectorPair->occurrenceGteVector  = vectorPair->occurrenceVector;
    default:
    __builtin_unreachable();
  }

  vectorPair->occurrenceVector    = _mm256_and_si256(bitmask, vectorPair->occurrenceVector);
  vectorPair->occurrenceGteVector = _mm256_and_si256(bitmask, vectorPair->occurrenceGteVector);
}


/*
 * Function:  awFmMakeAminoAcidOccurrenceVectorPair
 * --------------------
 * Computes the vector of characters before the given position equal to the given letter and
 *  the number of characters before the position greater than or equal to the character, and returns
 *  them in the vectorPair out-argument.
 *
 *  Inputs:
 *    blockPtr: address of the AwFmNucleotideBlock in which the query position resides.
 *    queryPosition: The global position of the query
 *    letter: letter for which the occurrence request is for.
 *    vectorPair: out-argument that will be used to return the occurrence vectors.
 */
inline void awFmMakeAminoAcidOccurrenceVectorPair(const struct AwFmAminoBlock *restrict const blockPtr,
  const uint64_t queryPosition, const uint8_t letter, struct AwFmOccurrenceVectorPair *vectorPair){

  //create the mask vector to mask away bits after the query position
  uint8_t localQueryPosition = getBlockQueryPositionFromGlobalPosition(queryPosition);
  __m256i bitmask = createQueryPositionBitmask(localQueryPosition, false, 0);

  //load the letter bit vectors
  const __m256i *restrict const blockVectorPtr = blockPtr->letterBitVectors;
  const __m256i bit0Vector = _mm256_load_si256(blockVectorPtr);
  const __m256i bit1Vector = _mm256_load_si256(blockVectorPtr + 1);
  const __m256i bit2Vector = _mm256_load_si256(blockVectorPtr + 2);
  const __m256i bit3Vector = _mm256_load_si256(blockVectorPtr + 3);
  const __m256i bit4Vector = _mm256_load_si256(blockVectorPtr + 4);

  //create the occurrence vectors
  switch(__builtin_expect(letter, 0)){
    case 0:   /*A (alanine) encoding 0b01100*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit3Vector, _mm256_andnot_si256(bit4Vector, bit2Vector));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_and_si256(bit3Vector, bit2Vector));
      break;
    case 1:   /*C (cysteine) encoding 0b10111*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit2Vector, _mm256_and_si256(bit1Vector, _mm256_andnot_si256(bit3Vector, bit0Vector))));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_and_si256(bit1Vector, bit0Vector))));
      break;
    case 2:   /*D (aspartic acid) encoding 0b00011*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit1Vector, _mm256_andnot_si256(bit4Vector, bit0Vector));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_and_si256(bit1Vector, bit0Vector))));
      break;
    case 3:   /*E (Glutamic acid) encoding 0b00110*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit4Vector, bit1Vector));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, bit1Vector)));
      break;
    case 4:   /*F (Phenylalanine) encoding 0b11110*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_andnot_si256(bit0Vector, bit1Vector))));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit2Vector, bit1Vector)));
      break;
    case 5:   /*G (Glycine) encoding 0b11010*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit0Vector, bit4Vector));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, bit1Vector)));
      break;
    case 6:   /*H (Histidine) encoding 0b11011*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit1Vector, _mm256_andnot_si256(bit2Vector, bit0Vector))));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_and_si256(bit1Vector, bit0Vector))));
      break;
    case 7:   /*I (Isoleucine) encoding 0b11001*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit1Vector, bit4Vector));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
      break;
    case 8:   /*K (Lysine) encoding 0b10101*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit1Vector, bit4Vector));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
      break;
    case 9:   /*L (aspartic acid) encoding 0b11100*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit1Vector, _mm256_andnot_si256(bit0Vector, bit4Vector));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, bit2Vector));
      break;
    case 10:  /*M (Methionine) encoding 0b11101*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_andnot_si256(bit1Vector, bit0Vector))));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
      break;
    case 11:  /*N (Asparagine) encoding 0b01000*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit4Vector, _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit1Vector, _mm256_andnot_si256(bit3Vector, bit0Vector))));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, bit3Vector);
      break;
    case 12:  /*P (Proline) encoding 0b01001*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit3Vector, _mm256_andnot_si256(bit4Vector, bit0Vector));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
      break;
    case 13:  /*Q (glutamine) encoding 0b00100*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit1Vector, _mm256_andnot_si256(bit0Vector, _mm256_andnot_si256(bit4Vector, bit2Vector))));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, bit2Vector));
      break;
    case 14:  /*R (Arginine) encoding 0b10011*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit2Vector, bit4Vector));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_and_si256(bit1Vector, bit0Vector))));
      break;
    case 15:  /*S (Serine) encoding 0b01010*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit3Vector, _mm256_andnot_si256(bit4Vector, bit1Vector));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_and_si256(bit3Vector, _mm256_or_si256(bit2Vector, bit1Vector)));
      break;
    case 16:  /*T (Threonine) encoding 0b00101*/
      vectorPair->occurrenceVector    = _mm256_and_si256(bit2Vector, _mm256_andnot_si256(bit4Vector, bit0Vector));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
      break;
    case 17:  /*V (Valine) encoding 0b10110*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit0Vector, bit4Vector));
      vectorPair->occurrenceGteVector = _mm256_and_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_and_si256(bit2Vector, bit1Vector)));
      break;
    case 18:  /*W (Tryptophan) encoding 0b00001*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit1Vector, _mm256_andnot_si256(bit4Vector, bit0Vector))));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_or_si256(bit2Vector, _mm256_or_si256(bit1Vector, bit0Vector))));
      break;
    case 19:  /*Y (Tyrosine) encoding 0b00010*/
      vectorPair->occurrenceVector    = _mm256_andnot_si256(bit3Vector, _mm256_andnot_si256(bit2Vector, _mm256_andnot_si256(bit0Vector, _mm256_andnot_si256(bit4Vector, bit1Vector))));
      vectorPair->occurrenceGteVector = _mm256_or_si256(bit4Vector, _mm256_or_si256(bit3Vector, _mm256_or_si256(bit2Vector, bit1Vector)));
      break;
    default: __builtin_unreachable();   //GCC respects this, doesn't check for letters that aren't valid
  }

  vectorPair->occurrenceVector    = _mm256_and_si256(bitmask, vectorPair->occurrenceVector);
  vectorPair->occurrenceGteVector = _mm256_and_si256(bitmask, vectorPair->occurrenceGteVector);
}

//TODO: try to have these be static and file scope, see if it can keep these in reg.
inline uint_fast8_t awFmVectorPopcount(const __m256i occurrenceVector){
  const __m256i lowBitsLookupTable  = _mm256_setr_epi8( 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0, 4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0);
  const __m256i highBitsLookupTable = _mm256_setr_epi8(-4,-3,-3,-2,-3,-2,-2,-1,-3,-2,-2,-1,-2,-1,-1, 0,-4,-3,-3,-2,-3,-2,-2,-1,-3,-2,-2,-1,-2,-1,-1, 0);

  const __m256i lowNybbleBitmasked          = _mm256_and_si256(occurrenceVector, _mm256_set1_epi8(0x0F));
  const __m256i lowNybbleBitCount           = _mm256_shuffle_epi8(lowBitsLookupTable, lowNybbleBitmasked);
  const __m256i highNybbleBits              = _mm256_srli_si256(occurrenceVector, 4);
  const __m256i highNybbleBitmasked         = _mm256_and_si256(highNybbleBits, _mm256_set1_epi8(0x0F));
  const __m256i highNybbleNegativeBitCount  = _mm256_shuffle_epi8(highBitsLookupTable, highNybbleBitmasked);
  const __m256i sadCountVector              = _mm256_sad_epu8(lowNybbleBitCount, highNybbleNegativeBitCount);
  //todo: try keeping nybble counts seperate, shifting both to be in the lower 128 bit regs, and performing separate hadds,
  //and extracting out values at end from both nybble vectors.

  //shift and add, placing the final two 16-bit sums in the least significant bits of each 128-bit lane.
  const uint16_t finalSum = _mm256_extract_epi16(sadCountVector, 0) + _mm256_extract_epi16(sadCountVector, 4) +
    _mm256_extract_epi16(sadCountVector, 8) + _mm256_extract_epi16(sadCountVector, 12);
  return finalSum;
}


inline void awFmBlockPrefetch(const uint8_t *restrict const baseBlockListPtr, const uint64_t blockByteWidth,
  const uint64_t nextQueryPosition){

  const uint64_t blockIndex    = getBlockIndexFromGlobalPosition(nextQueryPosition);
  //make the blockAddress pointer as a uint8_t* to make clean and easy pointer arithmetic when defining cache line boundries.
  const uint8_t *blockAddress  = (baseBlockListPtr + (blockIndex * blockByteWidth));

  for(uint_fast16_t prefetchOffset = 0; prefetchOffset < blockByteWidth; prefetchOffset += CACHE_LINE_SIZE_IN_BYTES){
    _mm_prefetch(blockAddress + prefetchOffset, _MM_HINT_T2);
  }
}


inline uint8_t awFmGetLetterAtBwtPosition(const union AwFmBwtBlockList blockList,
  const enum AwFmAlphabetType alphabet, const uint64_t bwtPosition){
  const size_t blockIndex      = getBlockIndexFromGlobalPosition(bwtPosition);
  const size_t positionInBlock = getBlockQueryPositionFromGlobalPosition(bwtPosition);
  const uint8_t byteInBlock    = positionInBlock / 8;
  const uint8_t bitInBlockByte = positionInBlock / 8;
  uint8_t letter = 0;
  uint8_t letterBitWidth = alphabet == AwFmAlphabetNucleotide? 2: 5;
  const __m256i *restrict const letterBitVectorPtr = alphabet == AwFmAlphabetNucleotide?
    blockList.asNucleotide[blockIndex].letterBitVectors:
    blockList.asAmino[blockIndex].letterBitVectors;

  for(uint8_t letterBit = letterBitWidth; letterBit >= 0; letterBit--){
    const uint8_t *restrict const blockVectorAsByteArray = (uint8_t*)(letterBitVectorPtr +letterBit);
    const uint8_t bit = (blockVectorAsByteArray[byteInBlock] >> bitInBlockByte) & 1;

    letter = (letter << 1) | bit;
  }

  return letter;
}


inline size_t awFmBacksetpBwtPosition(const struct AwFmIndex *restrict const index,
  const uint64_t bwtPosition){
    const enum AwFmAlphabetType alphabet        = index->metadata.alphabetType;
    const uint64_t *prefixSums                  = index->prefixSums;
    const uint64_t sentinelCharacterPosition    = index->backwardSentinelCharacterPosition;
    const uint64_t  blockIndex                  = getBlockIndexFromGlobalPosition(bwtPosition);

    uint64_t backtraceBwtPosition;
    uint8_t frequencyIndexLetter;
    const union AwFmBwtBlockList blockList = index->backwardBwtBlockList;
    frequencyIndexLetter = awFmGetLetterAtBwtPosition(blockList, alphabet, bwtPosition);

    struct AwFmOccurrenceVectorPair occurrenceVectors;
    if(alphabet == AwFmAlphabetNucleotide){
      awFmMakeNucleotideOccurrenceVectorPair(&blockList.asNucleotide[blockIndex], bwtPosition,
        frequencyIndexLetter, sentinelCharacterPosition, &occurrenceVectors);
    }
    else{
      awFmMakeAminoAcidOccurrenceVectorPair(&blockList.asAmino[blockIndex], bwtPosition,
        frequencyIndexLetter, &occurrenceVectors);
    }

    const uint_fast8_t vectorPopcount = awFmVectorPopcount(occurrenceVectors.occurrenceVector);
    backtraceBwtPosition = prefixSums[frequencyIndexLetter] + vectorPopcount;

    return backtraceBwtPosition;
  }



/*Private functions*/
/*
* Function:  getBlockIndexFromGlobalPosition
* --------------------
*  Computes the block index, given the full BWT query position.
*  Inputs:
*    globalQueryPosition: Position in the BWT that the occurrence function is requesting
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
 *    globalQueryPosition: Position in the BWT that the occurrence function is requesting
 *   Returns:
 *     Bit position into the block's AVX2 vectors where the query position lands.
 */
inline uint_fast8_t getBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition){
  return globalQueryPosition % POSITIONS_PER_FM_BLOCK;
}


/*
 * Function:  createQueryPositionBitmask
 * --------------------
 *  Creates an AVX2 vector that acts as a bitmask to remove positions after the query position,
 *    and to remove the sentinel character if it exists in this vector.
 *
 *  Inputs:
 *    localQueryPosition: position in this AVX2 vector to query. All bits after
 *      this position are cleared in the returned vector.
 *    containsSentinelCharacter: determines if the sentinel character is represented in this vector.
 *    sentinelCharacterPosition: if the containsSentinelCharacter argument is true,
 *      the bit at this position will be cleared.
 *
 *   Returns:
 *     Bitmask Vector for the occurrence vector.
 */
inline __m256i createQueryPositionBitmask(const uint8_t localQueryPosition,
  bool containsSentinelCharacter, uint8_t sentinelCharacterPosition){
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

  //if the block contains the sentinel character, mask it away now.
  if(__builtin_expect(containsSentinelCharacter, 0)){
    uint8_t sentinelBytePosition = localQueryPosition / 8;
    uint8_t sentinelBitPosition = localQueryPosition % 8;
    bitmaskArray[sentinelBytePosition] &= ~(1 << sentinelBitPosition);
  }

  //load and apply the bitmask
  return _mm256_lddqu_si256((__m256i*)bitmaskArray);
}
