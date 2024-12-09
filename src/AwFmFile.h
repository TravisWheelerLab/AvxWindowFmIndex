#ifndef AW_FM_INDEX_FILE_H
#define AW_FM_INDEX_FILE_H

#include <stdbool.h>
#include <stdint.h>
#include "AwFmIndexStruct.h"

/*
 * Function:  awFmGetSuffixArrayValueFromFile
 * --------------------
 * Reads the database sequence position from the suffix array inside the index
 * file. This function is used when the suffix array isn't being stored in
 * memory.  If it is being stored in memory, instead use the
 * awFmGetValueFromCompressedSuffixArray function in AwFmSuffixArray.h
 *
 *
 *  Inputs:
 *    index:   AwFmIndex struct containing the suffix array.
 *    positionInArray: position in the downsampled suffix array to retrieve from
 * file. Note that since this is the downsampled position, you should get this
 * value from backtracing until a position that is evenly divisible by the SA
 * sampling ratio, and then dividing by said ratio. Remember to add the offset
 * back to the value returned! valueOut: out-parameter for the value returned
 * from the file
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the read. Possible returns are:
 *      AwFmFileReadOkay on success.
 *      AwFmFileReadFail if the file could not be read sucessfully.
 */
enum AwFmReturnCode
awFmGetSuffixArrayValueFromFile(const struct AwFmIndex *_RESTRICT_ const index,
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
size_t
awFmGetSequenceFileOffset(const struct AwFmIndex *_RESTRICT_ const index);

/*
 * Function:  awFmGetSuffixArrayFileOffset
 * --------------------
 * Computes the file offset for the start of the compressed suffix array in the
 * AwFmIndex File.
 *
 *  Inputs:
 *    index: Pointer to the index struct to create the file offset for.
 *
 *  Returns:
 *    Offset into the file, in bytes, where the compressed suffix array starts.
 */
size_t
awFmGetSuffixArrayFileOffset(const struct AwFmIndex *_RESTRICT_ const index);

/*
 * Function:  awFmGetFastaVectorFileOffset
 * --------------------
 * Computes the file offset for the start of the compressed suffix array in the
 * AwFmIndex File.
 *
 *  Inputs:
 *    index: Pointer to the index struct to create the file offset for.
 *
 *  Returns:
 *    Offset into the file, in bytes, where the FastaVector data starts, if it
 * exists.
 */
size_t
awFmGetFastaVectorFileOffset(const struct AwFmIndex *_RESTRICT_ const index);

#endif /* end of include guard: AW_FM_INDEX_FILE_H */
