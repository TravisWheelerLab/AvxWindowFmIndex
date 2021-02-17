#include "AwFmSimdConfig.h"
#include "AwFmIndex.h"

#ifdef AW_FM_SIMD_CONFIG_ARM_NEON

  AwFmSimdVec256 AwFmSimdVecLoad(const AwFmSimdVec256* memAddr){
    AwFmSimdVec256 loadVector;
    loadVector.lowVec = vld1q_u8(memAddr);
    loadVector.highVec = vld1q_u8(memAddr+16); //each lane is 16B, so 2nd has offset of 16
    return loadVector;
  }

  AwFmSimdVec256 AwFmSimdVecAnd(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2){
    AwFmSimdVec256 resultVector;
    resultVector.lowVec   = vandq_u8(v1.lowVec,  v2.lowVec);
    resultVector.highVec  = vandq_u8(v1.highVec, v2.highVec);
    return resultVector;
  }

  AwFmSimdVec256 AwFmSimdVecOr(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2){
    AwFmSimdVec256 resultVector;
    resultVector.lowVec   = vorrq_u8(v1.lowVec,  v2.lowVec);
    resultVector.highVec  = vorrq_u8(v1.highVec, v2.highVec);
    return resultVector;
  }

  AwFmSimdVec256 AwFmSimdVecAndNot(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2){
    AwFmSimdVec256 resultVector;
    resultVector.lowVec   = vandq_u8(vmvnq_u8(v1.lowVec),  v2.lowVec);
    resultVector.highVec  = vandq_u8(vmvnq_u8(v1.highVec), v2.highVec);
    return resultVector;
  }

  uint32_t AwFmMaskedVectorPopcount(const AwFmSimdVec256 vec, const uint8_t localQueryPosition){
    uint64_t bitmasks[4] = {0};
    uint8_t  bitmaskedQuadWordIndex = localQueryPosition / 64;

    for(int8_t i = 0; i < bitmaskedQuadWordIndex; i++){
        bitmasks[i] = ~0UL;
    }
    bitmasks[bitmaskedQuadWordIndex] = ~0UL >> (63 - (localQueryPosition % 64));

    uint8x16_t lowVectorPopcnt  = vcntq_u8(vec.lowVec);
    uint8x16_t highVectorPopcnt = vcntq_u8(vec.highVec);
    uint8_t lowHorizontalSum    = vaddvq_u8(lowVectorPopcnt);
    uint8_t highHorizontalSum   = vaddvq_u8(highVectorPopcnt);

    return (uint32_t) lowHorizontalSum + (uint32_t) highHorizontalSum;
  }

  void AwFmSimdPrefetch(const void *memAddr){
    __builtin_prefetch(memAddr, 0, 0);
  }

#else

  AwFmSimdVec256 AwFmSimdVecLoad(const AwFmSimdVec256* memAddr){
      return _mm256_load_si256(memAddr);
  }

  AwFmSimdVec256 AwFmSimdVecAnd(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2){
    return _mm256_and_si256(v1, v2);
  }

  AwFmSimdVec256 AwFmSimdVecOr(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2){
    return _mm256_or_si256(v1, v2);
  }

  AwFmSimdVec256 AwFmSimdVecAndNot(const AwFmSimdVec256 v1, const AwFmSimdVec256 v2){
    return _mm256_andnot_si256(v1, v2);
  }

  uint32_t AwFmMaskedVectorPopcount(const AwFmSimdVec256 vec, const uint8_t localQueryPosition){
    uint64_t bitmasks[4] = {0};
    uint8_t  bitmaskedQuadWordIndex = localQueryPosition / 64;

    for(int8_t i = 0; i < bitmaskedQuadWordIndex; i++){
        bitmasks[i] = ~0UL;
    }
    bitmasks[bitmaskedQuadWordIndex] = ~0UL >> (63 - (localQueryPosition % 64));

    uint32_t popcount = _mm_popcnt_u64(_mm256_extract_epi64(vec, 0) & bitmasks[0]);
     popcount +=        _mm_popcnt_u64(_mm256_extract_epi64(vec, 1) & bitmasks[1]);
     popcount +=        _mm_popcnt_u64(_mm256_extract_epi64(vec, 2) & bitmasks[2]);
     popcount +=        _mm_popcnt_u64(_mm256_extract_epi64(vec, 3) & bitmasks[3]);

    return popcount;
  }


    void AwFmSimdPrefetch(const void *memAddr){
      _mm_prefetch(memAddr, _MM_HINT_NTA);  //prefetch with a non-temporal hint
    }

#endif
