#ifndef AW_FM_INDEX_FILE_H
#define AW_FM_INDEX_FILE_H

#include "AwFmIndex.h"
#include <stdint.h>

enum AwFmFileAccessCode{AwFmFileReadOkay, AwFmFileWriteOkay, AwFmFileOpenFail,
                        AwFmFileReadFail, AwFmFileWriteFail, AwFmFileFormatError,
                        AwFmAllocationFailure, AwFmIllegalPositionError};

enum AwFmFileAccessCode awFmCreateIndexFile(const char *restrict const fileSrc,
  const struct AwFmIndex *restrict const index, const size_t *restrict const fullSuffixArray,
  const uint8_t *restrict const databaseSequence, const bool allowOverwrite);


enum AwFmFileAccessCode awFmLoadIndexFromFile(const char *restrict const fileSrc,
  const struct AwFmIndex *restrict *restrict index);

enum AwFmFileAccessCode awFmLoadSuffixArrayIndexFromFile(const char *restrict const fileSrc,
  const struct AwFmIndex *restrict const index, const size_t compressedSuffixArrayPosition,
  const uint64_t *restrict dbPositionOut);

enum AwFmFileAccessCode awFmLoadSequenceSectionFromFile(const char *restrict const fileSrc,
  const struct AwFmIndex *restrict const index, const size_t sequencePosition,
  const size_t priorFlankingSequenceLength, const size_t postFlankingSequenceLength,
  const char **sequencePtr, const size_t *charactersRead)


#endif /* end of include guard: AW_FM_INDEX_FILE_H */
