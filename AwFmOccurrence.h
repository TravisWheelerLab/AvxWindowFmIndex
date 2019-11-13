#ifndef AW_FM_OCCURANCES_H
#define AW_FM_OCCURANCES_H

#include "AwFmIndex.h"
#include <stdint.h>

uint64_t awFmGetOccurrences(const struct AwFmIndex *restrict const index,
  const size_t queryPosition, const uint8_t letter);


  //deprecated
  // void  awFmOccurrencesDataPrefetch(const struct AwFmIndex *restrict const index, const uint64_t queryPosition);
void awFmAminoAcidDataPrefetch(const struct AwFmAminoBlock *restrict const blockList, const uint64_t queryPosition);


uint64_t awFmGetAminoAcidOccurrence(const struct AwFmAminoBlock *restrict const blockList, const size_t queryPosition, const uint8_t letter);

uint8_t awFmGetAminoAcidLetterAtBwtPosition(const struct AwFmAminoBlock *restrict const blockList, const uint64_t bwtPosition);

size_t awFmBackstepBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition);

#endif /* end of include guard: AW_FM_OCCURANCES_H */
