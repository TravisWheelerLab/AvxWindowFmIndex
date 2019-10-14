#ifndef AW_FM_INDEX_FILE_H
#define AW_FM_INDEX_FILE_H

#include "AwFmIndex.h"
#include <stdint.h>


/*Create an Awfmi file from the supplied AwFmIndex.*/
enum AwFmReturnCode awFmCreateIndexFile(const struct AwFmIndex *restrict const index,
  const size_t *restrict const fullSuffixArray, const char *restrict const databaseSequence,
  const bool allowOverwrite);

/*Allocate memory to the AwFmIndex pointer, and load it from the given fileSrc.*/
enum AwFmReturnCode awFmLoadIndexFromFile(const char *restrict const fileSrc,
  struct AwFmIndex *restrict *restrict index);

/*Loads the positions in the database sequence corresponding the the given array of BWT positions.
*   NOTE: positionArray is an array of BWT positions, not compressed suffix array indices.
*/
enum AwFmReturnCode awFmDbSequencePositionsFromSuffixArrayFile(const struct AwFmIndex *restrict const index,
  uint64_t *restrict positionArray, const uint64_t *restrict const offsetArray,
  const uint64_t positionArrayLength);

/*Loads a section of the database sequence from the Awmfi file.*/
enum AwFmReturnCode awFmLoadSequenceSectionFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankingSequenceLength,
  const size_t postFlankingSequenceLength, char **sequencePtr, size_t *charactersRead);

#endif /* end of include guard: AW_FM_INDEX_FILE_H */