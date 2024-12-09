#include <stdbool.h>
#include "AwFmOccurrence.h"
#include "AwFmIndex.h"
#include "AwFmLetter.h"

#define BYTES_PER_AVX2_REGISTER 32

AwFmSimdVec256 awFmMakeNucleotideOccurrenceVector(
    const struct AwFmNucleotideBlock *_RESTRICT_ const blockPtr,
    const uint8_t letter) {
  // load the letter bit vectors
  const AwFmSimdVec256 *_RESTRICT_ const blockVectorPtr =
      blockPtr->letterBitVectors;
  const AwFmSimdVec256 bit0Vector = AwFmSimdVecLoad(blockVectorPtr);
  const AwFmSimdVec256 bit1Vector = AwFmSimdVecLoad(blockVectorPtr + 1);
  const AwFmSimdVec256 bit2Vector = AwFmSimdVecLoad(blockVectorPtr + 2);

  switch (letter) {
  case 0: // Nucleotide A 0b110
    return AwFmSimdVecAnd(bit2Vector, bit1Vector);
  case 1: // Nucleotide C 0b101
    return AwFmSimdVecAnd(bit2Vector, bit0Vector);
  case 2: // Nucleotide G 0b011
    return AwFmSimdVecAnd(bit1Vector, bit0Vector);
  case 3: // Nucletoide T (or U) 0b001
    return AwFmSimdVecAndNot(bit2Vector,
                             AwFmSimdVecAndNot(bit1Vector, bit0Vector));
  case 4: // ambiguity character 'X' 0b010
    return AwFmSimdVecAndNot(bit2Vector,
                             AwFmSimdVecAndNot(bit0Vector, bit1Vector));
    // 0b100 is sentinel, but since you can't search for sentinels, it is not
    // included here.
  default:
    __builtin_unreachable();
  }
}

/*
 * Function:  awFmMakeAminoAcidOccurrenceVector
 * --------------------
 * Computes the vector of characters before the given position equal to the
 * given letter.
 *
 *  Inputs:
 *    blockPtr: Pointer to the AwFmNucleotideBlock in which the query position
 * resides. localQueryPosition: The local position of the query in the block
 *    letter: letter for which the occurrence request is for.
 *
 *  Returns:
 *   Vector with bits set at every position the given letter was found.
 */
AwFmSimdVec256 awFmMakeAminoAcidOccurrenceVector(
    const struct AwFmAminoBlock *_RESTRICT_ const blockPtr,
    const uint8_t letter) {

  // load the letter bit vectors
  const AwFmSimdVec256 *_RESTRICT_ const blockVectorPtr =
      blockPtr->letterBitVectors;
  const AwFmSimdVec256 bit0Vector = AwFmSimdVecLoad(blockVectorPtr);
  const AwFmSimdVec256 bit1Vector = AwFmSimdVecLoad(blockVectorPtr + 1);
  const AwFmSimdVec256 bit2Vector = AwFmSimdVecLoad(blockVectorPtr + 2);
  const AwFmSimdVec256 bit3Vector = AwFmSimdVecLoad(blockVectorPtr + 3);
  const AwFmSimdVec256 bit4Vector = AwFmSimdVecLoad(blockVectorPtr + 4);

  switch (__builtin_expect(letter, 0)) {
  case 0: /*A (alanine) encoding 0b01100*/
    return AwFmSimdVecAnd(bit3Vector,
                          AwFmSimdVecAndNot(bit4Vector, bit2Vector));
  case 1: /*C (cysteine) encoding 0b10111*/
    return AwFmSimdVecAnd(AwFmSimdVecAndNot(bit3Vector, bit2Vector),
                          AwFmSimdVecAnd(bit1Vector, bit0Vector));
  case 2: /*D (aspartic acid) encoding 0b00011*/
    return AwFmSimdVecAnd(bit1Vector,
                          AwFmSimdVecAndNot(bit4Vector, bit0Vector));
  case 3: /*E (Glutamic acid) encoding 0b00110*/
    return AwFmSimdVecAndNot(bit4Vector,
                             AwFmSimdVecAnd(bit2Vector, bit1Vector));
  case 4: /*F (Phenylalanine) encoding 0b11110*/
    return AwFmSimdVecAnd(AwFmSimdVecAndNot(bit0Vector, bit3Vector),
                          AwFmSimdVecAnd(bit2Vector, bit1Vector));
  case 5: /*G (Glycine) encoding 0b11010*/
    return AwFmSimdVecAndNot(bit2Vector,
                             AwFmSimdVecAndNot(bit0Vector, bit4Vector));
  case 6: /*H (Histidine) encoding 0b11011*/
    return AwFmSimdVecAnd(AwFmSimdVecAndNot(bit2Vector, bit3Vector),
                          AwFmSimdVecAnd(bit1Vector, bit0Vector));
  case 7: /*I (Isoleucine) encoding 0b11001*/
    return AwFmSimdVecAndNot(bit2Vector,
                             AwFmSimdVecAndNot(bit1Vector, bit4Vector));
  case 8: /*K (Lysine) encoding 0b10101*/
    return AwFmSimdVecAndNot(bit3Vector,
                             AwFmSimdVecAndNot(bit1Vector, bit4Vector));
  case 9: /*L (Leucine) encoding 0b11100*/
    return AwFmSimdVecAndNot(bit1Vector,
                             AwFmSimdVecAndNot(bit0Vector, bit4Vector));
  case 10: /*M (Methionine) encoding 0b11101*/
    return AwFmSimdVecAnd(AwFmSimdVecAndNot(bit1Vector, bit3Vector),
                          AwFmSimdVecAnd(bit2Vector, bit0Vector));
  case 11: /*N (Asparagine) encoding 0b01000*/
    return AwFmSimdVecAndNot(AwFmSimdVecOr(bit0Vector, bit1Vector),
                             AwFmSimdVecAndNot(bit2Vector, bit3Vector));
  case 12: /*P (Proline) encoding 0b01001*/
    return AwFmSimdVecAnd(bit3Vector,
                          AwFmSimdVecAndNot(bit4Vector, bit0Vector));
  case 13: /*Q (glutamine) encoding 0b00100*/
    return AwFmSimdVecAndNot(AwFmSimdVecOr(bit3Vector, bit1Vector),
                             AwFmSimdVecAndNot(bit0Vector, bit2Vector));
  case 14: /*R (Arginine) encoding 0b10011*/
    return AwFmSimdVecAndNot(bit3Vector,
                             AwFmSimdVecAndNot(bit2Vector, bit4Vector));
  case 15: /*S (Serine) encoding 0b01010*/
    return AwFmSimdVecAnd(bit3Vector,
                          AwFmSimdVecAndNot(bit4Vector, bit1Vector));
  case 16: /*T (Threonine) encoding 0b00101*/
    return AwFmSimdVecAnd(bit2Vector,
                          AwFmSimdVecAndNot(bit4Vector, bit0Vector));
  case 17: /*V (Valine) encoding 0b10110*/
    return AwFmSimdVecAndNot(bit3Vector,
                             AwFmSimdVecAndNot(bit0Vector, bit4Vector));
  case 18: /*W (Tryptophan) encoding 0b00001*/
    return AwFmSimdVecAndNot(AwFmSimdVecOr(bit3Vector, bit2Vector),
                             AwFmSimdVecAndNot(bit1Vector, bit0Vector));
  case 19: /*Y (Tyrosine) encoding 0b00010*/
    return AwFmSimdVecAndNot(AwFmSimdVecOr(bit0Vector, bit2Vector),
                             AwFmSimdVecAndNot(bit3Vector, bit1Vector));
  case 20: /*ambiguity character Z 0b11111 */
    return AwFmSimdVecAnd(AwFmSimdVecAnd(bit3Vector, bit2Vector),
                          AwFmSimdVecAnd(bit1Vector, bit0Vector));
  // 0b00000 is sentinel, but since you can't search for sentinels, it is not
  // included here.
  default:
    __builtin_unreachable(); // GCC respects this, doesn't check for letters
                             // that aren't valid
  }
}

inline void awFmBlockPrefetch(const void *_RESTRICT_ const baseBlockListPtr,
                              const uint64_t blockByteWidth,
                              const uint64_t nextQueryPosition) {

  const uint64_t blockIndex =
      awFmGetBlockIndexFromGlobalPosition(nextQueryPosition);
  // make the blockAddress pointer as a uint8_t* to make clean and easy pointer
  // arithmetic when defining cache line boundries.
  const uint8_t *blockAddress =
      ((uint8_t *)baseBlockListPtr + (blockIndex * blockByteWidth));

  for (uint_fast16_t prefetchOffset = 0; prefetchOffset < blockByteWidth;
       prefetchOffset += AW_FM_CACHE_LINE_SIZE_IN_BYTES) {
    AwFmSimdPrefetch(blockAddress + prefetchOffset);
  }
}

/*
 * Function:  awFmGetNucleotideLetterAtBwtPosition
 * --------------------
 * Given a specific position in the nucleotide BWT, returns the letter in the
 * BWT at this position.
 *
 *  Inputs:
 *    blockList: Block to query into. Note, this is usually not the pointer to
 * the bwt, but a pointer to the block the query will be found in. bwtPosition:
 * Position of the character to be returned, relative to the start of the block
 * to query.
 *
 *  Returns:
 *    Letter index of the character at the bwtPosition in the specified
 * blockList.
 */
uint8_t
awFmGetNucleotideLetterAtBwtPosition(const struct AwFmNucleotideBlock *blockPtr,
                                     const uint8_t localPosition) {
  const uint8_t byteInBlock = localPosition / 8;
  const uint8_t bitInBlockByte = localPosition % 8;

  const uint8_t *_RESTRICT_ const letterBytePointer =
      &((uint8_t *)&blockPtr->letterBitVectors)[byteInBlock];
  const uint8_t letterAsCompressedVector =
      ((letterBytePointer[0] >> bitInBlockByte) & 1) |
      ((letterBytePointer[32] >> bitInBlockByte) & 1) << 1 |
      ((letterBytePointer[64] >> bitInBlockByte) & 1) << 2;

  return awFmNucleotideCompressedVectorToLetterIndex(letterAsCompressedVector);
}

/*
 * Function:  awFmGetLetterAtBwtPosition
 * --------------------
 * Given a specific position in the amino BWT, returns the letter in the BWT at
 * this position.
 *
 *  Inputs:
 *    blockList: Block to query into. Note, this is usually not the pointer to
 * the bwt, but a pointer to the block the query will be found in. bwtPosition:
 * Position of the character to be returned, relative to the start of the block
 * to query.
 *
 *  Returns:
 *    Letter index of the character at the bwtPosition in the specified
 * blockList.
 */
uint8_t awFmGetAminoLetterAtBwtPosition(const struct AwFmAminoBlock *blockPtr,
                                        const uint8_t localPosition) {
  const uint8_t byteInBlock = localPosition / 8;
  const uint8_t bitInBlockByte = localPosition % 8;

  const uint8_t *_RESTRICT_ const letterBytePointer =
      &((uint8_t *)&blockPtr->letterBitVectors)[byteInBlock];
  const uint8_t letterAsCompressedVector =
      ((letterBytePointer[0] >> bitInBlockByte) & 1) |
      ((letterBytePointer[32] >> bitInBlockByte) & 1) << 1 |
      ((letterBytePointer[64] >> bitInBlockByte) & 1) << 2 |
      ((letterBytePointer[96] >> bitInBlockByte) & 1) << 3 |
      ((letterBytePointer[128] >> bitInBlockByte) & 1) << 4;

  return awFmAminoAcidCompressedVectorToLetterIndex(letterAsCompressedVector);
}
