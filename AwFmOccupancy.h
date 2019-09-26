#ifndef AW_FM_OCCUPANCY_H
#define AW_FM_OCCUPANCY_H

#include "AwFmIndex.h"
#include <stdint.h>

uint64_t awFmGetOccupancy(const struct AwFmIndex *restrict const index,
  const size_t queryPosition, const uint8_t letter);

void  awFmOccupancyDataPrefetch(const struct AwFmIndex *restrict const index, const uint64_t queryPosition);

uint8_t awFmGetLetterAtBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition);

size_t awFmBackstepBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition);

#endif /* end of include guard: AW_FM_OCCUPANCY_H */
