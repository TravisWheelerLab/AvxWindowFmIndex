#include <stdint.h>
#include "AwFmIndex.h"

inline struct AwFmBackwardRange awFmKmerTableLookup(const struct AwFmIndex *restrict const index,
  const uint8_t *restrict const kmerSeed);
