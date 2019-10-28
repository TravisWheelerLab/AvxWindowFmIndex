#ifndef AW_FM_INDEX_CREATE_H
#define AW_FM_INDEX_CREATE_H

#include "AwFmIndex.h"
#include <stdint.h>

enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex **indexPtr, const uint8_t *restrict const databaseSequence,
  const size_t databaseSequenceLength, const uint16_t suffixArrayCompressionRatio, const char *restrict const fileSrc);



//includes for testing, not part of public API
enum AwFmReturnCode awFmCreateFullSuffixArray(const uint8_t *databaseSequence,
  const uint64_t databaseSequenceLength, uint64_t **fullSuffixArrayOut);

enum AwFmReturnCode awFmCreateBlockList(struct AwFmIndex *restrict const index,
  const uint8_t *restrict const databaseSequence, const uint64_t databaseSequenceLength,
  uint64_t *totalOccupancies);

void awFmSetRankPrefixSums(struct AwFmIndex *restrict const index, const uint64_t *restrict const totalOccupancies);

void awFmInitBlock(struct AwFmIndex *const restrict index, const uint64_t blockIndex,
  uint64_t *totalOccupanciesSoFar, const size_t suffixArrayLength);


#endif /* end of include guard: AW_FM_INDEX_CREATE_H */
