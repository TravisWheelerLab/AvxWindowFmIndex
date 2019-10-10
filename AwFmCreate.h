#ifndef AW_FM_INDEX_CREATE_H
#define AW_FM_INDEX_CREATE_H

#include <stdint.h>
#include "AwFmIndex.h"

enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex **indexPtr, const uint8_t *restrict const databaseSequence,
  const size_t databaseSequenceLength, const uint16_t suffixArrayCompressionRatio, const char *restrict const fileSrc);

#endif /* end of include guard: AW_FM_INDEX_CREATE_H */
