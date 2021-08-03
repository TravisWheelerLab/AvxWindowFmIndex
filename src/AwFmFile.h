#ifndef AW_FM_INDEX_FILE_H
#define AW_FM_INDEX_FILE_H

#include "AwFmIndexStruct.h"

#include <stdbool.h>
#include <stdint.h>


/*
 * Function:  awFmReadPositionsFromSuffixArray
 * --------------------
 * Given an array of positions in the compressed suffix array, replaces each element in the positionArray
 *  with the sequence position found at that location in the compressed suffix array.
 *  Note that the positions in the position array are positions into the 'unsampled' suffix array,
 *  and will be divided by the compression ratio to find the index into the compressed
 *  suffix array stored in the file.
 *
 *  Inputs:
 *    index:          Pointer to the AwFmIndex struct that corresponds to the
 *      index file that stores the compressed suffix array.
 *    positionArray:  Array of indices in the compressed suffix array. The data
 *      at each index in this array will be replaced with the value stored at that
 *      position in the suffix array.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the read. Possible returns are:
 *      AwFmFileReadOkay on success.
 *      AwFmFileReadFail if the file could not be read sucessfully.
 */
enum AwFmReturnCode awFmReadPositionsFromSuffixArray(const struct AwFmIndex *restrict const index,
  uint64_t *restrict const positionArray, const size_t positionArrayLength);


/*
 * Function:  awFmSuffixArrayReadPositionParallel
 * --------------------
 * Reads the database sequence position from the suffix array using a backtrace.
 *  The backtrace contains the position in the full suffix array and the offset from backtracing.
 *
 *
 *  Inputs:
 *    index:   AwFmIndex struct containing the suffix array.
 *    backtracePtr: struct containing the Position and offset to look up in the suffix array.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the read. Possible returns are:
 *      AwFmFileReadOkay on success.
 *      AwFmFileReadFail if the file could not be read sucessfully.
 */
enum AwFmReturnCode awFmSuffixArrayReadPositionParallel(const struct AwFmIndex *restrict const index,
 struct AwFmBacktrace *restrict const backtracePtr);


 /*
  * Function:  awFmGetSuffixArrayValueFromFile
  * --------------------
  * Reads the database sequence position from the suffix array inside the index file.
  *   This function is used when the suffix array isn't being stored in memory.  If it
  *   is being stored in memory, instead use the awFmGetValueFromCompressedSuffixArray function
  *   in AwFmSuffixArray.h
  *
  *
  *  Inputs:
  *    index:   AwFmIndex struct containing the suffix array.
  *    positionInArray: position in the downsampled suffix array to retrieve from file.
  *     Note that since this is the downsampled position, you should get this value from
  *     backtracing until a position that is evenly divisible by the SA sampling ratio,
  *     and then dividing by said ratio. Remember to add the offset back to the value returned!
  *   valueOut: out-parameter for the value returned from the file
  *
  *  Returns:
  *    AwFmReturnCode represnting the result of the read. Possible returns are:
  *      AwFmFileReadOkay on success.
  *      AwFmFileReadFail if the file could not be read sucessfully.
  */
enum AwFmReturnCode awFmGetSuffixArrayValueFromFile(const struct AwFmIndex *restrict const index,
  const size_t positionInArray, size_t *valueOut);


/*
 * Function:  awFmGetSequenceFileOffset
 * --------------------
 * Computes the file offset for the start of the sequence in the AwFmIndex File.
 *
 *  Inputs:
 *    index: Pointer to the index struct to create the file offset for.
 *
 *  Returns:
 *    Offset into the file, in bytes, where the sequence starts.
 */
size_t awFmGetSequenceFileOffset(const struct AwFmIndex *restrict const index);


/*
 * Function:  awFmGetSuffixArrayFileOffset
 * --------------------
 * Computes the file offset for the start of the compressed suffix array in the AwFmIndex File.
 *
 *  Inputs:
 *    index: Pointer to the index struct to create the file offset for.
 *
 *  Returns:
 *    Offset into the file, in bytes, where the compressed suffix array starts.
 */
size_t awFmGetSuffixArrayFileOffset(const struct AwFmIndex *restrict const index);


/*
 * Function:  awFmGetFastaVectorFileOffset
 * --------------------
 * Computes the file offset for the start of the compressed suffix array in the AwFmIndex File.
 *
 *  Inputs:
 *    index: Pointer to the index struct to create the file offset for.
 *
 *  Returns:
 *    Offset into the file, in bytes, where the FastaVector data starts, if it exists.
 */
size_t awFmGetFastaVectorFileOffset(const struct AwFmIndex *restrict const index);

#endif /* end of include guard: AW_FM_INDEX_FILE_H */
