#ifndef AW_FM_INDEX_STRUCTS_H
#define AW_FM_INDEX_STRUCTS_H

#include "AwFmGlobals.h"
#include <stdint.h>
#include <stdbool.h>
#include <immintrin.h>
#include <stdio.h>


//TODO: move this to AwFmFile.c
#define AW_FM_INDEX_METADATA_BYTE_SIZE 64


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
  union   AwFmBwtBlockList  forwardBwtBlockList;
  union   AwFmBwtBlockList  backwardBwtBlockList;
          uint64_t          *prefixSums;
          uint64_t          *kmerSeedTable;
          FILE              *fileHandle;
  };


//TODO: update this block comment
/* enum AwFmReturnCode
* enum detailing the result a function that allocates data, performs disk I/O, or
*   otherwise might want to report an internal failure.
*
*   Successful actions will always generate positive numbers,
*   Failing actions will always generate negative numbers
*
*   AwFmFileReadOkay:         File was read sucessfully.
*   AwFmFileWriteOkay:        File was written successfully.
*   AwFmFileOpenFail:         File could not be opened.
*   AwFmFileReadFail:         An attempt to read the file failed.
*   AwFmFileWriteFail:        An attempt to write the file failed.
*     File may exist, but might be corrupted.
*   AwFmFileFormatError:      File was opened successfully, but the headed did not match.
*     This implies that the file provided was of an incorrect type, or was corrupted.
*   AwFmAllocationFailure:    Could not allocate memory to read in some portion of the file.
*   AwFmIllegalPositionError: The given position to read is outside the expected bounds of
*     the file's database sequence or compressed suffix array.
*   AwFmFileAccessAbandoned:  Something went wrong before attempting to open the file.
*   AwFmNoFileSrcGiven:       The fileSrc was null.
*/
enum AwFmReturnCode{
  AwFmSuccess             = 1,    AwFmFileReadOkay                = 2,    AwFmFileWriteOkay             = 3,
  AwFmGeneralFailure      = -1,   AwFmUnsupportedVersionError     = -2,   AwFmAllocationFailure         = -3,
  AwFmNullPtrError        = -4,   AwFmSuffixArrayCreationFailure  = -5,   AwFmIllegalPositionError      = -6,
  AwFmNoFileSrcGiven      = -7,   AwFmNoDatabaseSequenceGiven     = -8,   AwFmFileFormatError           = -9,
  AwFmFileOpenFail        = -10,  AwFmFileReadFail                = -11,  AwFmFileWriteFail             = -12,
  AwFmErrorDbSequenceNull = -13,  AwFmErrorSuffixArrayNull        = -14,  AwFmUnsupportedAlphabetError  = -15};

bool awFmReturnCodeSuccess(const enum AwFmReturnCode returnCode){
  return returnCode >= 0;
}

struct AwFmIndex *awFmIndexAlloc(const struct AwFmIndexMetadata *restrict const metadata,
  const size_t sequenceLength);
void awFmDeallocIndex(struct AwFmIndex *index);


uint_fast8_t  awFmGetAlphabetCardinality(const enum AwFmAlphabetType alphabet);
size_t        awFmGetKmerTableLength(const struct AwFmIndexMetadata *restrict const metadata);
size_t        awFmNumBlocksFromBwtLength(const size_t suffixArrayLength);
size_t        awFmNumBlocksFromSequenceLength(const size_t databaseSequenceLength);
bool          awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index, const uint64_t position);
uint64_t      awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index);
bool          awFmSearchRangeIsValid(const struct AwFmBackwardRange *restrict const searchRange);

#endif /* end of include guard: AW_FM_INDEX_STRUCTS_H */
