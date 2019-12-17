#ifndef AW_FM_INDEX_STRUCTS_H
#define AW_FM_INDEX_STRUCTS_H

#include "AwFmGlobals.h"
#include <stdint.h>
#include <stdbool.h>
#include <immintrin.h>
#include <stdio.h>


#define AwFmAlphabetAminoAcid                   1
#define AwFmAlphabetAminoAcidVectorsPerWindow   5
#define AwFmAlphabetAminoAcidCardinality        20

#define AwFmAlphabetNucleotide                  2
#define AwFMAlphabetNucleotideVectorsPerWindow  2
#define AwFmAlphabetNucleotideCardinality       4


struct AwFmBackwardRange{
  uint64_t startPtr;
  uint64_t endPtr;
};

struct AwFmBiDirectionalRange{
  uint64_t startPtr;
  uint64_t endPtr;
  uint64_t startPrimePtr;
};

enum AwFmSearchDirection{
  AwFmSearchDirectionBackward,
  AwFmSearchDirectionForward
};

enum AwFmAlphabetType{
  AwFmAlphabetTypeAmino = 1, AwFmAlphabetTypeNucleotide = 2};

enum AwFmBwtType{
  AwFmBwtTypeBackwardOnly = 1, AwFmBwtTypeBiDirectional = 2};

/*Structs*/
struct AwFmAminoBlock{
  __m256i   letterBitVectors[AwFmAlphabetAminoAcidVectorsPerWindow];
  uint64_t  baseOccurrences[AwFmAlphabetAminoAcidCardinality];
};

struct AwFmNucleotideBlock{
  __m256i   letterBitVectors[AwFMAlphabetNucleotideVectorsPerWindow];
  uint64_t  baseOccurrences[AwFmAlphabetNucleotideCardinality];
};

union AwFmBwtBlockList{
  struct AwFmNucleotideBlock  *asNucleotide;
  struct AwFmAminoBlock       *asAmino;
};

/*Struct for the metadata in the AwFmIndex struct.
* This contains data that helps to build the index, or to determine how the index
* will function.*/
struct AwFmIndexMetadata{
  uint32_t              versionNumber;
  uint16_t              suffixArrayCompressionRatio;
  uint8_t               kmerLengthInSeedTable;
  enum AwFmAlphabetType alphabetType;
  enum AwFmBwtType      bwtType;
};

struct AwFmIndex{
  struct  AwFmIndexMetadata metadata;
          uint64_t          backwardSentinelCharacterPosition;
          uint64_t          forwardSentinelCharacterPosition;
          uint64_t          bwtLength;
  union   AwFmBwtBlockList  backwardBwtBlockList;
  union   AwFmBwtBlockList  forwardBwtBlockList;
          uint64_t          *prefixSums;
  struct  AwFmBackwardRange *kmerSeedTable;
          FILE              *fileHandle;
  };

//todo: remove unused return codes
enum AwFmReturnCode{
  AwFmSuccess             = 1,    AwFmFileReadOkay                = 2,    AwFmFileWriteOkay             = 3,
  AwFmGeneralFailure      = -1,   AwFmUnsupportedVersionError     = -2,   AwFmAllocationFailure         = -3,
  AwFmNullPtrError        = -4,   AwFmSuffixArrayCreationFailure  = -5,   AwFmIllegalPositionError      = -6,
  AwFmNoFileSrcGiven      = -7,   AwFmNoDatabaseSequenceGiven     = -8,   AwFmFileFormatError           = -9,
  AwFmFileOpenFail        = -10,  AwFmFileReadFail                = -11,  AwFmFileWriteFail             = -12,
  AwFmErrorDbSequenceNull = -13,  AwFmErrorSuffixArrayNull        = -14};


/*
 * Function:  awFmIndexAlloc
 * --------------------
 * Dynamically allocates memory for the AwFmIndex struct and all internally stored arrays.
 *
 *  Inputs:
 *    metadata:         metadata struct that describes the format and parameters of the index.
 *      The metadata struct will be memcpy'd directly into the index.
 *    sequenceLength:   Length of the sequence the index is built from.
 *
 *  Returns:
 *    Allocated AwFmIndex struct, or NULL on an allocation failure.
 *      If any dynamic allocation fails, all data used in the AwFmIndex will be deallocated, too.
 */
struct AwFmIndex *awFmIndexAlloc(const struct AwFmIndexMetadata *restrict const metadata,
  const size_t sequenceLength);


/*
 * Function:  awFmDeallocIndex
 * --------------------
 * Deallocates the given AwFmIndex struct, as well as all internal arrays.
 *   Performing this deallocation will also close the internally stored file handle.
 *
 *  Inputs:
 *    index:  Pointer to the AwFmIndex struct to deallocate
 */
void          awFmDeallocIndex(struct AwFmIndex *index);


/*
 * Function:  awFmGetAlphabetCardinality
 * --------------------
 * Returns the number of letters in the given alphabet.
 *  If AwFmAlphabetNucleotide is given, returns 4.
 *  If AwFmAlphabetAminoAcid  is given, returns 20.
 *  Inputs:
 *    alphabet: Alphabet to query.
 *
 *  Returns:
 *    Cardinality of the alphabet.
 */
uint_fast8_t  awFmGetAlphabetCardinality(const enum AwFmAlphabetType alphabet);


/*
 * Function:  awFmGetKmerTableLength
 * --------------------
 * Computes the number of AwFmBackwardRange structs in the kmerSeedTable.
 *  This value is equal to |A|^k, where |A| is the cardinalty of the alphabet as
 *    set in the given metadata, and k is the length of the kmers in the lookup table,
 *    also as set in the metadata.
 *  Inputs:
 *    metadata: Metadata struct corresponding to the AwFmIndex.
 *
 *  Returns:
 *    Number of AwFmBackwardRange structs in the table.
 */
size_t        awFmGetKmerTableLength(const struct AwFmIndexMetadata *restrict const metadata);


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
size_t        awFmNumBlocksFromBwtLength(const size_t suffixArrayLength);


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
bool          awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index, const uint64_t position);


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
uint64_t      awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index);


/*
 * Function:  awFmSearchRangeIsValid
 * --------------------
 * Determines if the given AwFmBackwardRange represents represents a range of positions, or
 *  if no positions are represented by the search range.
 *  A search range is valid if and only if the start ptr is less than or equal to the end pointer.
 *
 *  Inputs:
 *    searchRange: Pointer to the search range to query.
 *
 *  Returns:
 *    True if the search range represents a valid range of positions, or false if it represents no elements.
 */
bool          awFmSearchRangeIsValid(const struct AwFmBackwardRange *restrict const searchRange);


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
bool          awFmReturnCodeSuccess(const enum AwFmReturnCode returnCode);

#endif /* end of include guard: AW_FM_INDEX_STRUCTS_H */
