#ifndef AWFM_SIMD_CONFIG_H
#define AWFM_SIMD_CONFIG_H

#include <stdint.h>
#include "AwFmIndex.h"

/*
 * Function:  AwFmSimdVecLoad
 * --------------------
 *  Loads a 256-bit SIMD vector in an architecture-depedent way.
 *    This function either loads a single 256-bit AVX2 register on x86,
 *    or 2 128-bit vectors on ARM neon.
 *
 *  Inputs:
 *    memAddr: Memory address to load, this should be aligned to a 32B boundary
 *
 *  Returns:
 *    AwFmSimdVec256 containing the loaded data.
 */
AwFmSimdVec256 AwFmSimdVecLoad(const AwFmSimdVec256 *memAddr);

/*
 * Function:  AwFmSimdVecAnd
 * --------------------
 *  Performs a bitwise AND between the data in the two given vectors.

 *
 *  Inputs:
 *    v1 and v2 are the vectors to bitwise AND together
 *
 *  Returns:
 *    AwFmSimdVec256 result of the bitwise AND
 */
AwFmSimdVec256 AwFmSimdVecAnd(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2);

/*
 * Function:  AwFmSimdVecOr
 * --------------------
 *  Performs a bitwise OR between the data in the two given vectors.

 *
 *  Inputs:
 *    v1 and v2 are the vectors to bitwise OR together
 *
 *  Returns:
 *    AwFmSimdVec256 result of the bitwise OR
 */
AwFmSimdVec256 AwFmSimdVecOr(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2);

/*
 * Function:  AwFmSimdVecAndNot
 * --------------------
 *  This function bitwise NOT's the vector v1, and then bitwise ANDs the result
 * with v2.
 *
 *  Inputs:
 *    v1 and v2 are the vectors to to perform the operation on.
 *
 *  Returns:
 *    AwFmSimdVec256 result of NOTing v1 and ANDing it with v2.
 */
AwFmSimdVec256 AwFmSimdVecAndNot(const AwFmSimdVec256 v1,
                                 const AwFmSimdVec256 v2);

/*
 * Function:  AwFmMaskedVectorPopcount
 * --------------------
 *  Creates a 256-bit bitmask that masks away all bits after localQueryPosition,
 * and ANDs the result with vec. A population count is then taken of the result.
 *
 *  Inputs:
 *    vec: populated vector to take the masked popcount of.
 *    localQueryPosition: position in the vector for which all 1s up to this
 * position are to be counted.
 *
 *  Returns:
 *    a population count of the vector up to the localQueryPosition.
 */
uint32_t AwFmMaskedVectorPopcount(const AwFmSimdVec256 vec,
                                  const uint8_t localQueryPosition);

/*
 * Function:  AwFmSimdPrefetch
 * --------------------
 *  Performs a prefetch of the data at the given address, if able.
 *    In the current implementation, a prefetch is only performed for x86,
 *    the ARM implementation of this function does not perform a prefetch.
 *  Inputs:
 *    memAddr: address of memory to prefetch.
 */
void AwFmSimdPrefetch(const void *memAddr);

#endif
