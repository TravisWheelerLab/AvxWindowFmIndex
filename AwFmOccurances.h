#ifndef AW_FM_OCCURANCES_H
#define AW_FM_OCCURANCES_H

#include "AwFmIndex.h"
#include <stdint.h>

uint64_t awFmGetOccurances(const struct AwFmIndex *restrict const index,
  const size_t queryPosition, const uint8_t letter);

void  awFmOccurancesDataPrefetch(const struct AwFmIndex *restrict const index, const uint64_t queryPosition);

uint8_t awFmGetLetterAtBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition);

size_t awFmBackstepBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition);

#endif /* end of include guard: AW_FM_OCCURANCES_H */
