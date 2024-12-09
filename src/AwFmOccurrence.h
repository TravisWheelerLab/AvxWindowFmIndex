#ifndef AW_FM_OCCURANCE_H
#define AW_FM_OCCURANCE_H

#include <stdint.h>
#include "AwFmIndexStruct.h"
#include "AwFmSimdConfig.h"

/*
 * Function:  awFmMakeNucleotideOccurrenceVector
 * --------------------
 * Computes the vector of characters before the given position equal to the
 * given letter. Note that this function does not check for the sentinel
 * character being in this vector.
 *
 *  Inputs:
 *    blockPtr: Pointer to the AwFmNucleotideBlock in which the query position
 * resides. localQueryPosition: The local position of the query in the block
 *    letter: letter for which the occurrence request is for.
 *
 *  Returns:
 *   Vector with bits set at every position the given letter was found.
 */
AwFmSimdVec256 awFmMakeNucleotideOccurrenceVector(
    const struct AwFmNucleotideBlock *_RESTRICT_ const blockPtr,
    const uint8_t letter);

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
    const uint8_t letter);

/*
 * Function:  awFmVectorPopcount
 * --------------------
 * Computes the number of bits set in the given AVX2 vector.
 *  This function is similar to the one proposed in Mula(2018), with slight
 * improvements.
 *
 *  Inputs:
 *    occurrenceVector: vector to perform the popcount on.
 *
 *  Returns:
 *    Count of the bits set in the occurrenceVector.
 */
uint16_t awFmVectorPopcount(const AwFmSimdVec256 occurrenceVector,
                            const uint8_t localQueryPosition);

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
uint16_t awFmVectorPopcountBuiltin(const AwFmSimdVec256 occurrenceVector,
                                   const uint8_t localQueryPosition);

/*
 * Function:  awFmBlockPrefetch
 * --------------------
 * Requests that the CPU prefetch the block at the given query position into
 * cache.
 *
 *  Inputs:
 *    baseBlockListPtr: pointer to the blockList to prefetch into.
 *    blockByteWidth: width of the block, being either sizeof(struct
 * AwFmNucleotideBlock) or sizeof(struct AwFmAminoBlock). This is left as an
 * argument to help unnecessary branching. nextQueryPosition: position in the
 * blockList that contains the block that should be prefetched.
 */
void awFmBlockPrefetch(const void *_RESTRICT_ const baseBlockListPtr,
                       const uint64_t blockByteWidth,
                       const uint64_t nextQueryPosition);

/*
 * Function:  awFmGetNucleotideLetterAtBwtPosition
 * --------------------
 * Given a specific position in the nucleotide BWT, returns the letter in the
 * BWT at this position.
 *
 *  Inputs:
 *    blockList: blockList to be queried.
 *    localPosition: Position of the character to be returned.
 *
 *  Returns:
 *    letter at the bwtPosition in the specified blockList.
 */
uint8_t
awFmGetNucleotideLetterAtBwtPosition(const struct AwFmNucleotideBlock *blockPtr,
                                     const uint8_t localPosition);

/*
 * Function:  awFmGetLetterAtBwtPosition
 * --------------------
 * Given a specific position in the BWT, returns the letter in the BWT at this
 * position.
 *
 *  Inputs:
 *    blockList: blockList to be queried.
 *    localPosition: Position of the character to be returned.
 *
 *  Returns:
 *    letter at the bwtPosition in the specified blockList.
 */
uint8_t awFmGetAminoLetterAtBwtPosition(const struct AwFmAminoBlock *blockPtr,
                                        const uint8_t localPosition);

#endif /* end of include guard: AW_FM_OCCURANCE_H */
