#include <stdint.h>

#ifdef AW_FM_SIMD_CONFIG_M1
  //M1 config that uses AVX to emulate AVX2
  #include <emmintrin.h>
  #include <smmintrin.h>

  typedef struct __m256i{
    __m128i lowVector;
    __m128i highVector;
  } __m256i;

  //emulates the AVX2 load instruction
  __m256i _mm256_load_si256 (__m256i const * mem_addr){
    __m256i loadVector;
    loadVector.lowVector = _mm_load_si128(mem_addr);
    loadVector.highVector = _mm_load_si128(mem_addr + sizeof(__mm128i))
    return loadVector;
  }


  //emulates the AVX2 AND instruction
  __m256i _mm256_and_si256(const __m256i v1, const __m256i v2){
    __m256i andVector;
    andVector.lowVector   = _mm_and_si128(v1.lowVector, v2.lowVector);
    andVector.highVector  = _mm_and_si128(v1.highVector, v2.highVector);
    return andVector;
  }

  //emulates the AVX2 ANDNOT instruction
  __m256i _mm256_andnot_si256(const __m256i v1, const __m256i v2){
    __m256i andNotVector;
    andNotVector.lowVector   = _mm_andnot_si128(v1.lowVector, v2.lowVector);
    andNotVector.highVector  = _mm_andnot_si128(v1.highVector, v2.highVector);
    return andNotVector;
  }

  //emulates the AVX2 OR instruction
  __m256i _mm256_or_si256(const __m256i v1, const __m256i v2){
    __m256i orVector;
    orVector.lowVector   = _mm_andnot_si128(v1.lowVector, v2.lowVector);
    orVector.highVector  = _mm_andnot_si128(v1.highVector, v2.highVector);
    return orVector;
  }

  //emulates the AVX2 extract function
  inline __int64 _mm256_extract_epi64(const __m256i vec, const int index){
    if(index < 2){
      return _mm_extract_epi64(vec.lowVector, index & 1);
    }
    else{
      return _mm_extract_epi64(vec.highVector, index & 1);
    }
  }

#else
  //if no SIMD configuration was given, attempt to build as normal(AVX2 support)
  #include <immintrin.h>

#endif
