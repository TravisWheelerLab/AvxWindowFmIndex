#ifndef AW_FM_OCCURANCE_H
#define AW_FM_OCCURANCE_H

#include "AwFmIndex.h"
// #include "AwFmSearch.h"
#include <stdint.h>


struct AwFmOccurrenceVectorPair{
  __m256i occurrenceVector;
  __m256i occurrenceGteVector;
};


/*
 * Function:  awFmMakeNucleotideOccurrenceVectorPair
 * --------------------
 * Computes the vector of characters before the given position equal to the given letter and
 *  the number of characters before the position greater than or equal to the character, and returns
 *  them in the vectorPair out-argument.
 *
 *  Inputs:
 *    blockPtr: Pointer to the AwFmNucleotideBlock in which the query position resides.
 *    localQueryPosition: The local position of the query in the block
 *    letter: letter for which the occurrence request is for.
 *    sentinelCharacterPosition: Global position of the sentinel character.
 *
 *  Returns:
 *    Pair of the occurrence and occurrenceGte vectors
 */
struct AwFmOccurrenceVectorPair awFmMakeNucleotideOccurrenceVectorPair(const struct AwFmNucleotideBlock *restrict const blockPtr,
  const uint8_t localQueryPosition, const uint8_t letter);

/*
 * Function:  awFmMakeNucleotideOccurrenceVector
 * --------------------
 * Computes the vector of characters before the given position equal to the given letter.
 *  Note that this function does not check for the sentinel character being in this vector.
 *
 *  Inputs:
 *    blockPtr: Pointer to the AwFmNucleotideBlock in which the query position resides.
 *    localQueryPosition: The local position of the query in the block
 *    letter: letter for which the occurrence request is for.
 *
 *  Returns:
 *   Vector with bits set at every position the given letter was found.
 */
__m256i awFmMakeNucleotideOccurrenceVector(const struct AwFmNucleotideBlock *restrict const blockPtr,
  const uint8_t localQueryPosition, const uint8_t letter);



/*
 * Function:  awFmMakeAminoAcidOccurrenceVectorPair
 * --------------------
 * Computes the vector of characters before the given position equal to the given letter and
 *  the number of characters before the position greater than or equal to the character, and returns
 *  them in the vectorPair out-argument.
 *
 *  Inputs:
 *    blockPtr: Pointer to the AwFmNucleotideBlock in which the query position resides.
 *    localQueryPosition: The local position of the query in the block
 *    letter: letter for which the occurrence request is for.
 *    vectorPair: out-argument that will be used to return the occurrence vectors.
 *
 *  Returns:
 *    Pair of the occurrence and occurrenceGte vectors
 */
struct AwFmOccurrenceVectorPair awFmMakeAminoAcidOccurrenceVectorPair(const struct AwFmAminoBlock *restrict const blockPtr,
  const uint8_t localQueryPosition, const uint8_t letter);


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
__m256i awFmMakeAminoAcidOccurrenceVector(const struct AwFmAminoBlock *restrict const blockPtr,
  const uint8_t localQueryPosition, const uint8_t letter);


/*
 * Function:  awFmVectorPopcount
 * --------------------
 * Computes the number of bits set in the given AVX2 vector.
 *  This function is similar to the one proposed in Mula(2018), with slight improvements.
 *
 *  Inputs:
 *    occurrenceVector: vector to perform the popcount on.
 *
 *  Returns:
 *    Count of the bits set in the occurrenceVector.
 */
uint16_t awFmVectorPopcount(const __m256i occurrenceVector);

/*
 * Function:  awFmVectorPopcountBuiltin
 * --------------------
 * Computes the number of bits set in the given AVX2 vector.
 *  This just uses the builtin _mm_popcnt_u64 instruction.
 *   on arrays of this size, this is like 4 times faster than Mula
 *
 *  Inputs:
 *    occurrenceVector: vector to perform the popcount on.
 *
 *  Returns:
 *    Count of the bits set in the occurrenceVector.
 */
uint16_t awFmVectorPopcountBuiltin(const __m256i occurrenceVector);

/*
 * Function:  awFmBlockPrefetch
 * --------------------
 * Requests that the CPU prefetch the block at the given query position into cache.
 *
 *  Inputs:
 *    baseBlockListPtr: pointer to the blockList to prefetch into.
 *    blockByteWidth: width of the block, being either sizeof(struct AwFmNucleotideBlock) or sizeof(struct AwFmAminoBlock).
 *      This is left as an argument to help unnecessary branching.
 *    nextQueryPosition: position in the blockList that contains the block that should be prefetched.
 */
void awFmBlockPrefetch(const void *restrict const baseBlockListPtr, const uint64_t blockByteWidth,
  const uint64_t nextQueryPosition);


/*
 * Function:  awFmGetLetterAtBwtPosition
 * --------------------
 * Given a specific position in the BWT, returns the letter in the BWT at this position.
 *
 *  Inputs:
 *    blockList: blockList to be queried.
 *    alphabet: alphabet of the index, either Nucleotide or Amino
 *    bwtPosition: Position of the character to be returned.
 *
 *  Returns:
 *    letter at the bwtPosition in the specified blockList.
 */
uint8_t awFmGetLetterAtBwtPosition(const union AwFmBwtBlockList blockList, const enum AwFmAlphabetType alphabet, const uint64_t bwtPosition);



#endif /* end of include guard: AW_FM_OCCURANCE_H */
