#ifndef AW_FM_SUFFIX_ARRAY_H
#define AW_FM_SUFFIX_ARRAY_H

#include "AwFmIndex.h"

struct AwFmSuffixArrayOffset {
	size_t byteOffset;
	uint8_t bitOffset;
};

/*
 * Function:  awFmInitCompressedSuffixArray
 * --------------------
 * Initializes the bit-compressed suffix array from a fully-sampled, uncompressed uint64_t* suffix array,
 *  Note that the fullSa must be dynamically allocated, and this function will realloc that array
 *  into the final array that will be kept to store the SA. This is an in-place algorithm, and will
 *  clobber anything inside fullSa. While these caveats may be annoying, they can save substantial amounts
 *  of memory in practice while building an fm-index.
 *
 *  Inputs:
 *  fullSa: An uncompressed, fully-sampled suffix array for the entire sequence.
 *  saLength: length of the suffix array, in elements.
 *  compressedSuffixArray: suffixArray struct that will manage the new compressed suffix array.
 *  samplingRatio: how often to sample the suffix array. High values result in more sparse suffix arrays,
 *    and therefore less memory usage but slower locate() functionality.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the construction. Possible results are
 *      AwFmSuccess on success.
 *      AwFmAllocationFailure if memory could not be properly allocated. If this happens,
 *        make sure that fullSa was dynamically allocated, and therefore can be realloced by this function.
 */
enum AwFmReturnCode awFmInitCompressedSuffixArray(
		uint64_t *fullSa, size_t saLength, struct AwFmCompressedSuffixArray *compressedSuffixArray, uint8_t samplingRatio);


/*
 * Function:  awFmGetValueFromCompressedSuffixArray
 * --------------------
 * Gets the value stored at the given (sampled) position in the suffix array.
 *  This function requires that the suffix array be stored in memory for the index.'
 *
 *    To use SA value from file (if SA isn't being stored in memory), use awFmGetSuffixArrayValueFromFile() in
 * AwFmFile.h
 *
 *  Inputs:
 *  suffixArray: struct containing the compressed, sampled suffix array.
 *  positionInArray: downsampled position in the suffix array to retrieve. As it is a downsampled position,
 *    the position * the sampling ratio should give you the original BWT position generated by backtracing.
 *
 *  Returns:
 *    Value at the given position in the compressed suffix array.
 */
size_t awFmGetValueFromCompressedSuffixArray(
		const struct AwFmCompressedSuffixArray *suffixArray, size_t positionInArray);
// to use SA value from file (if SA isn't being stored in memory), use awFmGetSuffixArrayValueFromFile in AwFmFile.h


/*
 * Function:  awFmComputeCompressedSaSizeInBytes
 * --------------------
 * computes the number of bytes required to store the bit-compressed, downsampled suffix array.
 *  Note that the suffix array adds 8 bytes to the end to act as padding to prevent buffer overflows,
 *  so this value will be 8 bytes larger than is directly needed to store all values.
 *
 *  Inputs:
 *  saLength: length of the unsampled suffix array, aka, the bwtLength.
 *  samplingRatio: how often to take samples from the suffix array. From index->config.suffixArrayCompressionRatio
 *
 *  Returns:
 *    Number of bytes required to store the compressed, downsampled suffix array.
 */
size_t awFmComputeCompressedSaSizeInBytes(size_t saLength, uint8_t samplingRatio);


/*
 * Function:  awFmComputeCompressedSaSizeInBytes
 * --------------------
 * computes the number of bytes required to store the bit-compressed, downsampled suffix array.
 *
 *  Inputs:
 *  saLength: length of the unsampled suffix array, aka, the bwtLength.
 *  samplingRatio: how often to take samples from the suffix array. From index->config.suffixArrayCompressionRatio
 *
 *  Returns:
 *    Number of bytes required to store the compressed, downsampled suffix array.
 */
uint8_t awFmComputeSuffixArrayValueMinWidth(const size_t saLength);


/*
 * Function:  awFmGetSampledSuffixArrayLength
 * --------------------
 * Returns the number of samples in a downsampled suffix array.
 *
 *  Inputs:
 *  bwtLength: number of elements in the bwt. This value is equal to 1 + the sequence length.
 *    This value can also be considered the length of the unsampled suffix array.
 *  compressionRatio: how often the suffix array is sampled.
 *
 *  Returns:
 *    Number of elements in the downsampled suffix array.
 */
size_t awFmGetSampledSuffixArrayLength(uint64_t bwtLength, uint64_t compressionRatio);


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
 *    positionArrayLength: number of values in the positionArray
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the read. Possible returns are:
 *      AwFmSuccess as long as the suffix array is stored in memory, not on disk.
 *      AwFmFileReadOkay on success.
 *      AwFmFileReadFail if the file could not be read sucessfully.
 */
enum AwFmReturnCode awFmReadPositionsFromSuffixArray(
		const struct AwFmIndex *const index, uint64_t *restrict const positionArray, const size_t positionArrayLength);


struct AwFmSuffixArrayOffset awFmGetOffsetIntoSuffixArrayByteArray(
		const uint8_t compressedValueBitWidth, const size_t indexOfValueInCompressedSa);


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
enum AwFmReturnCode awFmSuffixArrayReadPositionParallel(
		const struct AwFmIndex *restrict const index, struct AwFmBacktrace *restrict const backtracePtr);

#endif
