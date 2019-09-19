#ifndef AW_FM_OCCUPANCY_H
#define AW_FM_OCCUPANCY_H

#include "AwFmIndex"
#include <stdint.h>

void awFmGetOccupancy(const struct AwFmBlock *restrict const blockList,
  const size_t queryPosition, const uint8_t letter);

inline void awFmRequestBlockDataPrefetch(const struct fm_block *const restrict blockList, const uint64_t blockIndex);

#endif /* end of include guard: AW_FM_OCCUPANCY_H */
