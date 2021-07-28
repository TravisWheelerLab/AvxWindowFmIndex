#ifndef AW_FM_SUFFIX_ARRAY_H
#define AW_FM_SUFFIX_ARRAY_H

#include "AwFmIndex.h"

struct AwFmCompressedSuffixArray{
  uint8_t *values;
  size_t length;
  uint8_t valueBitWidth;
  uint64_t compressedByteLength;
};


struct AwFmCompressedSuffixArray awFmInitSuffixArray(uint64_t *fullSa, size_t saLength);

size_t awFmGetValueFromCompressedSuffixArray(struct AwFmCompressedSuffixArray suffixArray, size_t valuePosition);



#endif
