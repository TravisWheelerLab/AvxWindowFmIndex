#ifndef AWFM_SIMD_CONFIG_H
#define AWFM_SIMD_CONFIG_H

#include <stdint.h>
#include "AwFmIndex.h"



AwFmSimdVec256 AwFmSimdVecLoad(const AwFmSimdVec256* memAddr);
AwFmSimdVec256 AwFmSimdVecAnd(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2);
AwFmSimdVec256 AwFmSimdVecOr(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2);
AwFmSimdVec256 AwFmSimdVecAndNot(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2);
uint32_t AwFmMaskedVectorPopcount(const AwFmSimdVec256 vec, const uint8_t localQueryPosition);
void AwFmSimdPrefetch(const void *memAddr);

#endif
