#ifndef AW_FM_INDEX_FILE_H
#define AW_FM_INDEX_FILE_H

#include "AwFmIndex.h"
#include <stdint.h>


enum AwFmFileAccessCode{AwFmFileReadOkay, AwFmFileWriteOkay, AwFmFileOpenFail,
                        AwFmFileReadFail, AwFmFileWriteFail, AwFmFileFormatError,
                        AwFmAllocationFailure, AwFmIllegalPositionError,
                        AwFmFileAccessAbandoned, AwFmNoFileSrcGiven};

enum AwFmFileAccessCode awFmCreateIndexFile(const struct AwFmIndex *restrict const index,
  const size_t *restrict const fullSuffixArray, const char *restrict const databaseSequence,
  const bool allowOverwrite);

enum AwFmFileAccessCode awFmLoadIndexFromFile(const char *restrict const fileSrc,
  struct AwFmIndex *restrict *restrict index);

enum AwFmFileAccessCode awFmDbSequencePositionsFromSuffixArrayFile(const struct AwFmIndex *restrict const index,
  uint64_t *restrict positionArray, const uint64_t *restrict const offsetArray,
  const uint64_t positionArrayLength);

enum AwFmFileAccessCode awFmLoadSequenceSectionFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankingSequenceLength,
  const size_t postFlankingSequenceLength, char **sequencePtr, size_t *charactersRead);

#endif /* end of include guard: AW_FM_INDEX_FILE_H */
