#ifndef AW_FM_INDEX_STRUCT_H
#define AW_FM_INDEX_STRUCT_H

#include "AwFmIndex.h"
#include <immintrin.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

/*
 * Function:  awFmIndexAlloc
 * --------------------
 * Dynamically allocates memory for the AwFmIndex struct and all internally stored arrays.
 *
 *  Inputs:
 *    metadata:         metadata struct that describes the format and parameters of the index.
 *      The metadata struct will be memcpy'd directly into the index.
 *    bwtLength:   Length of the BWT, in positions, that the index will hold
 *
 *  Returns:
 *    Allocated AwFmIndex struct, or NULL on an allocation failure.
 *      If any dynamic allocation fails, all data used in the AwFmIndex will be deallocated, too.
 */
struct AwFmIndex *awFmIndexAlloc(const struct AwFmIndexMetadata *restrict const metadata,
                                 const size_t bwtLength);

/*
 * Function:  awFmGetAlphabetCardinality
 * --------------------
 * Returns the number of letters in the given alphabet.
 *  If AwFmAlphabetNucleotide is given, returns 4.
 *  If AwFmAlphabetAmino  is given, returns 20.
 *  Inputs:
 *    alphabet: Alphabet to query.
 *
 *  Returns:
 *    Cardinality of the alphabet.
 */
uint_fast8_t awFmGetAlphabetCardinality(const enum AwFmAlphabetType alphabet);

/*
 * Function:  awFmGetKmerTableLength
 * --------------------
 * Computes the number of AwFmSearchRange structs in the kmerSeedTable.
 *  This value is equal to |A|^k, where |A| is the cardinalty of the alphabet as
 *    set in the given metadata, and k is the length of the kmers in the lookup table,
 *    also as set in the metadata.
 *  Inputs:
 *    index: AwFmIndex struct that contains the kmerSeedTable.
 *
 *  Returns:
 *    Number of AwFmSearchRange structs in the table.
 */
size_t awFmGetKmerTableLength(const struct AwFmIndex *restrict index);

/*
 * Function:  awFmNumBlocksFromBwtLength
 * --------------------
 * Computes the number of BWT blocks that represent the BWT or suffix array of the given length.
 *  This function essentially acts as the ceiling of (suffixArrayLength/ positions per block)
 *
 *  Inputs:
 *    suffixArrayLength: Length of the suffix array, and implicitly, the BWT, in blocks..
 *
 *  Returns:
 *    Number of blocks required to store the BWT.
 */
size_t awFmNumBlocksFromBwtLength(const size_t suffixArrayLength);

/*
 * Function:  awFmGetPrefixSumsLength
 * --------------------
 * Computes the length of the prefix sums array for the AwFmIndex struct.
 *  The length is based on the alphabet type. Nucleotide indices have a length
 *  of 5, and amino indices have length 21.
 *
 *  Inputs:
 *    alphabet: The alphabet of the index.
 *
 *  Returns:
 *    Number of uint64_t elements in the prefixSums array.
 */
uint8_t awFmGetPrefixSumsLength(const enum AwFmAlphabetType alphabet);

/*
 * Function:  awFmBwtPositionIsSampled
 * --------------------
 * Determines if the given position in the BWT is sampled in the compressed suffix array.
 *
 *  Inputs:
 *    index: AwFmIndex struct representing the BWT and suffix array.
 *    position: position in the BWT to query if it is sampled in the suffix array
 *
 *  Returns:
 *    True if the given position is sampled in the suffix array, false otherwise.
 */
bool awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index,
                              const uint64_t position);

/*
 * Function:  awFmGetCompressedSuffixArrayLength
 * --------------------
 * Calculates the length of the compressed suffix array.
 *
 *  Inputs:
 *    index: AwFmIndex struct representing the BWT and suffix array.
 *
 *  Returns:
 *    Number of positions in the compressed suffix array.
 */
uint64_t awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index);

/*
 * Function:  awFmSearchRangeIsValid
 * --------------------
 * Determines if the given AwFmSearchRange represents represents a range of positions, or
 *  if no positions are represented by the search range.
 *  A search range is valid if and only if the start ptr is less than or equal to the end pointer.
 *
 *  Inputs:
 *    searchRange: Pointer to the search range to query.
 *
 *  Returns:
 *    True if the search range represents a valid range of positions, or false if it represents no
 * elements.
 */
bool awFmSearchRangeIsValid(const struct AwFmSearchRange *restrict const searchRange);

/*
 * Function:  awFmReturnCodeSuccess
 * --------------------
 * With a given return code, returns true if the return code represented a successful state,
 *  or false if it represeted a failure state.
 *
 *  Inputs:
 *    returnCode: Return code to query for success.
 *
 *  Returns:
 *    True if the return code represents a successful action, otherwise returns false.
 */
bool awFmReturnCodeSuccess(const enum AwFmReturnCode returnCode);

/*
 * Function:  getBlockIndexFromGlobalPosition
 * --------------------
 *  Computes the block index, given the full BWT query position.
 *  Inputs:
 *    globalQueryPosition: Position in the BWT that the occurrence function is requesting
 *   Returns:
 *     Index of the block where the given query position resides.
 */
size_t awFmGetBlockIndexFromGlobalPosition(const size_t globalQueryPosition);

/*
 * Function:  getBlockQueryPositionFromGlobalPosition
 * --------------------
 *  Computes bit position inside the block that represents the given full BWT position.
 *  Inputs:
 *    globalQueryPosition: Position in the BWT that the occurrence function is requesting
 *   Returns:
 *     Bit position into the block's AVX2 vectors where the query position lands.
 */
uint_fast8_t awFmGetBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition);

/*
 * Function:  awFmSearchRangeLength
 * --------------------
 * Gets the number of positions included in the given AwFmSearchRange
 *
 *  Inputs:
 *    range: Range of positions in the BWT that corresponds to some number of
 *      instances of a given kmer.
 *
 *  Returns:
 *    Number of positions in the given range if the range is valid (startPtr < endPtr),
 *      or 0 otherwise, as that would imply that no instances of that kmer were found.
 */
size_t awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range);

#endif /* end of include guard: AW_FM_INDEX_STRUCT_H */
