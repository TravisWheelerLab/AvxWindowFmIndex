#ifndef AWFM_SIMD_CONFIG_H
#define AWFM_SIMD_CONFIG_H

#include <stdint.h>


// define the Simd vector type, which is determined by the architecture we're building for.
#ifdef AW_FM_SIMD_CONFIG_ARM_NEON
  #include <arm_neon.h>

  typedef struct AwFmSimdVec256{
    uint8x16_t lowVec;
    uint8x16_t highVec;
  } AwFmSimdVec256;

#else
  #include <immintrin.h>

  typedef __m256i AwFmSimdVec256;
#endif


AwFmSimdVec256 AwFmSimdVecLoad(const AwFmSimdVec256* mem_addr);
AwFmSimdVec256 AwFmSimdVecAnd(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2);
AwFmSimdVec256 AwFmSimdVecOr(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2);
AwFmSimdVec256 AwFmSimdVecAndNot(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2);
uint32_t AwFmMaskedVectorPopcount(const AwFmSimdVec256 vec, const uint8_t localQueryPosition);

#endif
