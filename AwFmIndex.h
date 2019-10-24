#ifndef AW_FM_INDEX_STRUCTS_H
#define AW_FM_INDEX_STRUCTS_H

#include "AwFmGlobals.h"
#include <stdint.h>
#include <stdbool.h>
#include <immintrin.h>


#define AW_FM_INDEX_METADATA_BYTE_SIZE 64

/*Structs*/
struct AwFmBlock{
  __m256i   letterBitVectors[5];
  uint64_t  baseOccupancies[20];
};

struct AwFmIndexMetadata{
  uint32_t  versionNumber;
};

union AwFmIndexPaddedMetadata{
  struct AwFmIndexMetadata data;
  uint8_t padding[AW_FM_INDEX_METADATA_BYTE_SIZE];
};

struct AwFmIndex{
  struct AwFmBlock *blockList;
  uint64_t  rankPrefixSums[AMINO_CARDINALITY + 1]; //last position acts as BWT length
  uint64_t  numBlocks;
  uint16_t  suffixArrayCompressionRatio;
  union AwFmIndexPaddedMetadata metadata;
  char            *fileSrc;
  const uint8_t   *databaseSequence;  //usually NULL, used in construction and saving to file
  uint64_t        *fullSuffixArray;   //usually NULL, used in construction and saving to file
};

struct AwFmSearchRange{
  uint64_t startPtr;
  uint64_t endPtr;
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
  AwFmSuccess         = 1,    AwFmFileReadOkay                = 2,    AwFmFileWriteOkay         = 3,
  AwFmGeneralFailure  = -1,   AwFmUnsupportedVersionError     = -2,   AwFmAllocationFailure     = -3,
  AwFmNullPtrError    = -4,   AwFmSuffixArrayCreationFailure  = -5,   AwFmIllegalPositionError  = -6,
  AwFmNoFileSrcGiven  = -7,   AwFmNoDatabaseSequenceGiven     = -8,   AwFmFileFormatError       = -9,
  AwFmFileOpenFail    = -10,  AwFmFileReadFail                = -11,  AwFmFileWriteFail         = -12};


enum AwFmReturnCode awFmIndexSetFileSrc(struct AwFmIndex *restrict const index, const char *restrict const fileSrc);
struct  AwFmIndex *awFmAlignedAllocAwFmIndex(void);
struct  AwFmBlock *awFmAlignedAllocBlockList(const size_t numBlocks);
        void      awFmDeallocateFmIndex(struct AwFmIndex *restrict index);
        void      awFmDeallocFullSuffixArray(struct AwFmIndex *const restrict index);
        size_t    awFmNumBlocksFromSuffixArrayLength(const size_t suffixArrayLength);
        size_t    awFmNumBlocksFromSequenceLength(const size_t databaseSequenceLength);
        bool      awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index, const uint64_t position);
        size_t    awFmSearchRangeLength(const struct AwFmSearchRange *restrict const range);
        uint64_t  awFmGetBwtLength(const struct AwFmIndex *restrict const index);
        uint64_t  awFmGetDbSequenceLength(const struct AwFmIndex *restrict const index);
        uint64_t  awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index);
        bool      awFmSearchRangeIsValid(const struct AwFmSearchRange *restrict const searchRange);

#endif /* end of include guard: AW_FM_INDEX_STRUCTS_H */
