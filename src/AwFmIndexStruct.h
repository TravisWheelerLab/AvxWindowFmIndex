#ifndef AW_FM_INDEX_STRUCT_H
#define AW_FM_INDEX_STRUCT_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include "AwFmIndex.h"

#define AW_FM_CURRENT_VERSION_NUMBER 8
#define AW_FM_FEATURE_FLAG_BIT_FASTA_VECTOR 0

/*
 * Function:  awFmIndexAlloc
 * --------------------
 * Dynamically allocates memory for the AwFmIndex struct and all internally
 * stored arrays.
 *
 *  Inputs:
 *    config:         configuration struct that describes the format and
 * parameters of the index. The config struct will be memcpy'd directly into the
 * index. bwtLength:   Length of the BWT, in positions, that the index will hold
 *
 *  Returns:
 *    Allocated AwFmIndex struct, or NULL on an allocation failure.
 *      If any dynamic allocation fails, all data used in the AwFmIndex will be
 * deallocated, too.
 */
struct AwFmIndex *
awFmIndexAlloc(const struct AwFmIndexConfiguration *_RESTRICT_ const config,
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
 *    set in the given configuration, and k is the length of the kmers in the
 * lookup table, also as set in the config. Inputs: index: AwFmIndex struct that
 * contains the kmerSeedTable.
 *
 *  Returns:
 *    Number of AwFmSearchRange structs in the table.
 */
size_t awFmGetKmerTableLength(const struct AwFmIndex *_RESTRICT_ index);

/*
 * Function:  awFmNumBlocksFromBwtLength
 * --------------------
 * Computes the number of BWT blocks that represent the BWT or suffix array of
 * the given length. This function essentially acts as the ceiling of
 * (suffixArrayLength/ positions per block)
 *
 *  Inputs:
 *    suffixArrayLength: Length of the suffix array, and implicitly, the BWT, in
 * blocks..
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
 * Determines if the given position in the BWT is sampled in the compressed
 * suffix array.
 *
 *  Inputs:
 *    index: AwFmIndex struct representing the BWT and suffix array.
 *    position: position in the BWT to query if it is sampled in the suffix
 * array
 *
 *  Returns:
 *    True if the given position is sampled in the suffix array, false
 * otherwise.
 */
bool awFmBwtPositionIsSampled(const struct AwFmIndex *_RESTRICT_ const index,
                              const uint64_t position);

/*
 * Function:  awFmGetCompressedSuffixArrayLength
 * --------------------
 * Calculates the number of elements in the the sampled suffix array.
 *
 *  Inputs:
 *    index: AwFmIndex struct representing the BWT and suffix array.
 *
 *  Returns:
 *    Number of positions in the compressed suffix array.
 */
uint64_t awFmGetCompressedSuffixArrayLength(
    const struct AwFmIndex *_RESTRICT_ const index);

/*
 * Function:  awFmSearchRangeIsValid
 * --------------------
 * Determines if the given AwFmSearchRange represents represents a range of
 * positions, or if no positions are represented by the search range. A search
 * range is valid if and only if the start ptr is less than or equal to the end
 * pointer.
 *
 *  Inputs:
 *    searchRange: Pointer to the search range to query.
 *
 *  Returns:
 *    True if the search range represents a valid range of positions, or false
 * if it represents no elements.
 */
bool awFmSearchRangeIsValid(
    const struct AwFmSearchRange *_RESTRICT_ const searchRange);

/*
 * Function:  awFmReturnCodeSuccess
 * --------------------
 * With a given return code, returns true if the return code represented a
 * successful state, or false if it represeted a failure state.
 *
 *  Inputs:
 *    returnCode: Return code to query for success.
 *
 *  Returns:
 *    True if the return code represents a successful action, otherwise returns
 * false.
 */
bool awFmReturnCodeSuccess(const enum AwFmReturnCode returnCode);

/*
 * Function:  getBlockIndexFromGlobalPosition
 * --------------------
 *  Computes the block index, given the full BWT query position.
 *  Inputs:
 *    globalQueryPosition: Position in the BWT that the occurrence function is
 * requesting Returns: Index of the block where the given query position
 * resides.
 */
size_t awFmGetBlockIndexFromGlobalPosition(const size_t globalQueryPosition);

/*
 * Function:  getBlockQueryPositionFromGlobalPosition
 * --------------------
 *  Computes bit position inside the block that represents the given full BWT
 * position. Inputs: globalQueryPosition: Position in the BWT that the
 * occurrence function is requesting Returns: Bit position into the block's AVX2
 * vectors where the query position lands.
 */
uint_fast8_t
awFmGetBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition);

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
 *    Number of positions in the given range if the range is valid (startPtr <
 * endPtr), or 0 otherwise, as that would imply that no instances of that kmer
 * were found.
 */
size_t
awFmSearchRangeLength(const struct AwFmSearchRange *_RESTRICT_ const range);

/*
 * Function:  awFmIndexIsVersionValid
 * --------------------
 * returns true if the given version number is one that is currently supported.
 *   The supported version numbers are enumerated at the top of this header
 * (AwFmIndexStruct.h)
 *
 *  Inputs:
 *    versionNumber: version number from the AwFmIndex struct's config.
 *
 *  Returns:
 *    true if the versionNumber is one of the supported versions.
 */
bool awFmIndexIsVersionValid(const uint16_t versionNumber);

/*
 * Function:  awFmIndexContainsFastaVector
 * --------------------
 * returns true if indices of the given version number contain a FastaVector
 * struct to catalog headers and the position of each sequence in the index.
 * This function should be updated every time a new version is added.
 *
 *  Inputs:
 *    index: allocated and constructed index to check if it contains fastaVector
 * information
 *
 *  Returns:
 *    True if the given version contains a FastaVector struct
 */
bool awFmIndexContainsFastaVector(
    const struct AwFmIndex *_RESTRICT_ const index);

#endif /* end of include guard: AW_FM_INDEX_STRUCT_H */
