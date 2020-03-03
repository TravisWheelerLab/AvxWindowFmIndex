#ifndef AW_FM_INDEX_FILE_H
#define AW_FM_INDEX_FILE_H

#include "AwFmIndex.h"
#include "AwFmBacktraceVector.h"

#include <stdbool.h>
#include <stdint.h>

#define AW_FM_INDEX_FILE_EXTENSION      "awfmi"

/*
 * Function:  awFmWriteIndexToFile
 * --------------------
 * With a given AwFmIndex and associated data, writes the index to the file.
 *  If you, the user, want to create a new index file, use the awFmCreateIndex
 *  function in AwFmCreate.h
 *
 *  Inputs:
 *    index:          AwFmIndex struct to be written.
 *    suffixArray:    Full (uncompressed) suffix array that the AwFmIndex is built from.
 *    sequence:       Database sequence that the AwFmIndex is built from.
 *    sequenceLength: Length of the sequence.
 *    fileSrc:        File path to write the Index file to.
 *    allowOverwrite: If set, will allow overwriting the file at the given fileSrc.
 *      If allowOverwite is false, will return error code AwFmFileAlreadyExists.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the write. Possible returns are:
 *      AwFmFileWriteOkay on success.
 *      AwFmFileAlreadyExists if a file exists at the given fileSrc, but allowOverwite was false.
 *      AwFmFileWriteFail if a file write failed.
 */
enum AwFmReturnCode awFmWriteIndexToFile(struct AwFmIndex *restrict const index,
  const uint64_t *restrict const suffixArray, const uint8_t *restrict const sequence,
  const uint64_t sequenceLength, const char *restrict const fileSrc, const bool allowOverwrite);


/*
 * Function:  awFmReadIndexFromFile
 * --------------------
 * Reads the AwFmIndex file from the given fileSrc
 *
 *  Inputs:
 *    index:          double pointer to an unallocated AwFmIndex to be allocated and populated
 *      by this function.
 *    fileSrc:        Path to the file containing the AwFmIndex.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the write. Possible returns are:
 *      AwFmFileReadOkay on success.
 *      AwFmFileAlreadyExists if no file could be opened at the given fileSrc.
 *      AwFmFileFormatError if the header was not correct. This suggests that the file at this location
 *        is not the correct format.
 *      AwFmAllocationFailure on failure to allocated the necessary memory for the index.
 */
enum AwFmReturnCode awFmReadIndexFromFile(struct AwFmIndex *restrict *restrict index, const char *fileSrc);


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
 * Function:  awFmReadSequenceFromFile
 * --------------------
 * Given a sequence position, reads a section of sequence surrounding that position from the
 *  corresponding index file.
 *
 *  Inputs:
 *    index:          Pointer to the AwFmIndex struct that contains the file handle to read.
 *      index file that stores the compressed suffix array.
 *    sequencePosition:        Position in the sequence to read.
 *    priorFlankLength: How many character before the given sequence position to include
 *      in the sequence segment.
 *    postFlankLength:  How many characters after the given sequence position to include
 *      in the sequence segement.
 *    sequenceBuffer: Pointer to the buffer to read the sequence segment into. This
 *      buffer  must be large enough to hold (postFlankLength + priorFlankLength +1) characters.
 *      the + 1 on this length is for the null terminator added to the end of the string.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the read. Possible returns are:
 *      AwFmFileReadOkay on success.
 *      AwFmFileReadFail if the file could not be read sucessfully.
 */
enum AwFmReturnCode awFmReadSequenceFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankLength, const size_t postFlankLength,
  char *const sequenceBuffer);


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


#endif /* end of include guard: AW_FM_INDEX_FILE_H */
