#ifndef AW_FM_SUFFIX_ARRAY_H
#define AW_FM_SUFFIX_ARRAY_H

#include "AwFmIndex.h"

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
enum AwFmReturnCode awFmInitCompressedSuffixArray(uint64_t *fullSa, size_t saLength,
  struct AwFmCompressedSuffixArray *compressedSuffixArray, uint8_t samplingRatio);


  /*
   * Function:  awFmGetValueFromCompressedSuffixArray
   * --------------------
   * Gets the value stored at the given (sampled) position in the suffix array.
   *  This function requires that the suffix array be stored in memory for the index.'
   *
   *    To use SA value from file (if SA isn't being stored in memory), use awFmGetSuffixArrayValueFromFile() in AwFmFile.h
   *
   *  Inputs:
   *  suffixArray: struct containing the compressed, sampled suffix array.
   *  positionInArray: downsampled position in the suffix array to retrieve. As it is a downsampled position,
   *    the position * the sampling ratio should give you the original BWT position generated by backtracing.
   *
   *  Returns:
   *    Value at the given position in the compressed suffix array.
   */
size_t awFmGetValueFromCompressedSuffixArray(struct AwFmCompressedSuffixArray suffixArray, size_t positionInArray);
// to use SA value from file (if SA isn't being stored in memory), use awFmGetSuffixArrayValueFromFile in AwFmFile.h


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


#endif