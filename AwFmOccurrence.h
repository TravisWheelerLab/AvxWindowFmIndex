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
 *    blockPtr: address of the AwFmNucleotideBlock in which the query position resides.
 *    queryPosition: The global position of the query
 *    letter: letter for which the occurrence request is for.
 *    sentinelCharacterPosition: Global position of the sentinel character.
 *    vectorPair: out-argument that will be used to return the occurrence vectors.
 */
void awFmMakeNucleotideOccurrenceVectorPair(const struct AwFmNucleotideBlock *restrict const blockPtr,
  const uint64_t queryPosition, const uint8_t letter, const uint64_t sentinelCharacterPosition,
  struct AwFmOccurrenceVectorPair *occurrenceVectors);


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
void awFmMakeAminoAcidOccurrenceVectorPair(const struct AwFmAminoBlock *restrict const blockPtr,
  const uint64_t queryPosition, const uint8_t letter, struct AwFmOccurrenceVectorPair *vectorPair);


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
uint_fast8_t awFmVectorPopcount(const __m256i occurrenceVector);


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
