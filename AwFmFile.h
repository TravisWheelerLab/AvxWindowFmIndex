#ifndef AW_FM_INDEX_FILE_H
#define AW_FM_INDEX_FILE_H

#include "AwFmIndex.h"
#include <stdint.h>

/*struct detailing the result of a file read/write request
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
enum AwFmFileAccessCode{AwFmFileReadOkay, AwFmFileWriteOkay, AwFmFileOpenFail,
                        AwFmFileReadFail, AwFmFileWriteFail, AwFmFileFormatError,
                        AwFmAllocationFailure, AwFmIllegalPositionError,
                        AwFmFileAccessAbandoned, AwFmNoFileSrcGiven};

/*Create an Awfmi file from the supplied AwFmIndex.*/
enum AwFmFileAccessCode awFmCreateIndexFile(const struct AwFmIndex *restrict const index,
  const size_t *restrict const fullSuffixArray, const char *restrict const databaseSequence,
  const bool allowOverwrite);

/*Allocate memory to the AwFmIndex pointer, and load it from the given fileSrc.*/
enum AwFmFileAccessCode awFmLoadIndexFromFile(const char *restrict const fileSrc,
  struct AwFmIndex *restrict *restrict index);

/*Loads the positions in the database sequence corresponding the the given array of BWT positions.
*   NOTE: positionArray is an array of BWT positions, not compressed suffix array indices.
*/
enum AwFmFileAccessCode awFmDbSequencePositionsFromSuffixArrayFile(const struct AwFmIndex *restrict const index,
  uint64_t *restrict positionArray, const uint64_t *restrict const offsetArray,
  const uint64_t positionArrayLength);

/*Loads a section of the database sequence from the Awmfi file.*/
enum AwFmFileAccessCode awFmLoadSequenceSectionFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankingSequenceLength,
  const size_t postFlankingSequenceLength, char **sequencePtr, size_t *charactersRead);

#endif /* end of include guard: AW_FM_INDEX_FILE_H */
