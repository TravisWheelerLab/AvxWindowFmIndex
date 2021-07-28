#ifndef AW_FM_SUFFIX_ARRAY_H
#define AW_FM_SUFFIX_ARRAY_H

#include "AwFmIndex.h"

struct AwFmCompressedSuffixArray{
  uint8_t   valueBitWidth;
  bool      responsibleForArrayDeallocation;
  uint8_t   *values;
  uint64_t  compressedByteLength;
};

enum AwFmReturnCode awFmInitSuffixArray(uint64_t *fullSa, size_t saLength,
  struct AwFmCompressedSuffixArray *compressedSuffixArray, bool buildSaInPlace);

size_t awFmGetValueFromCompressedSuffixArray(struct AwFmCompressedSuffixArray suffixArray, size_t valuePosition);



#endif
