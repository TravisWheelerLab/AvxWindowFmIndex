#include "AwFmSuffixArray.h"

#include <assert.h>
#include <string.h>
#include "AwFmFile.h"

// adding padding bytes prevents buffer overflow problems recalling values from
// the suffix array.
#define AW_FM_SUFFIX_ARRAY_END_PADDING_BYTES 8

// simple log2 ceiling implementation, thanks builtin clzll!
uint8_t log2Floor(const uint64_t a) { return 64 - __builtin_clzll(a); }

// the number of bits required to store the suffix array values is just
// log2Floor of the highest val.
uint8_t awFmComputeSuffixArrayValueMinWidth(const size_t saLength) {
  return log2Floor(saLength - 1);
}

// TODO: use this to eliminate code redundency in figuring out which byte/bit a
// value should start on.
struct AwFmSuffixArrayOffset
awFmGetOffsetIntoSuffixArrayByteArray(const uint8_t compressedValueBitWidth,
                                      const size_t indexOfValueInCompressedSa) {
  struct AwFmSuffixArrayOffset offset;

  // if we group up the non-standard bit width values into groups of 8, the
  // lengths will always be a multiple of 8, thus aligning to a byte. Then, we
  // can add the num bytes of the remainder of values, and then find the bit
  // offset.
  size_t alignedByteOffset =
      (indexOfValueInCompressedSa / 8) * compressedValueBitWidth;
  size_t numEndingBits =
      (indexOfValueInCompressedSa % 8) * compressedValueBitWidth;
  offset.byteOffset = alignedByteOffset + (numEndingBits / 8);
  offset.bitOffset = numEndingBits % 8;

  return offset;
}

size_t awFmComputeCompressedSaSizeInBytes(size_t saLength,
                                          uint8_t samplingRatio) {
  uint64_t sampledSaLength =
      awFmGetSampledSuffixArrayLength(saLength, samplingRatio);
  uint8_t valueMinBitWidth = awFmComputeSuffixArrayValueMinWidth(saLength);
  struct AwFmSuffixArrayOffset offset =
      awFmGetOffsetIntoSuffixArrayByteArray(valueMinBitWidth, sampledSaLength);
  if (offset.bitOffset != 0) {
    offset.byteOffset++;
  }

  return offset.byteOffset + AW_FM_SUFFIX_ARRAY_END_PADDING_BYTES;
}

// note: fullSa MUST be dynamically allocated, this reallocs the array and
// claims ownership. the reallocated array will be freed when awFmDeallocIndex()
// is called.
enum AwFmReturnCode awFmInitCompressedSuffixArray(
    uint64_t *fullSa, size_t saLength,
    struct AwFmCompressedSuffixArray *compressedSuffixArray,
    uint8_t samplingRatio) {

  if (fullSa == NULL || compressedSuffixArray == NULL) {
    return AwFmNullPtrError;
  }
  uint8_t minimumBitWidth = awFmComputeSuffixArrayValueMinWidth(saLength);
  size_t compressedSaByteSize =
      awFmComputeCompressedSaSizeInBytes(saLength, samplingRatio);

  compressedSuffixArray->values = (uint8_t *)fullSa;
  compressedSuffixArray->valueBitWidth = minimumBitWidth;
  compressedSuffixArray->compressedByteLength = compressedSaByteSize;

  const size_t numSaSamples =
      awFmGetSampledSuffixArrayLength(saLength, samplingRatio);
  for (size_t i = 0; i < numSaSamples; i++) {
    uint64_t saValue = fullSa[i * samplingRatio];

    struct AwFmSuffixArrayOffset offset =
        awFmGetOffsetIntoSuffixArrayByteArray(minimumBitWidth, i);

    uint8_t byteToWrite = (saValue << offset.bitOffset) & 0xFF;

    uint8_t bitmask =
        (1ULL << offset.bitOffset) -
        1; // zero out the bits we're going to set, since it's in-place
    compressedSuffixArray->values[offset.byteOffset] &= bitmask;
    compressedSuffixArray->values[offset.byteOffset++] |= byteToWrite;
    int8_t bitsWrittenOnFirstWriteAction = 8 - offset.bitOffset;
    int8_t bitsRemaining =
        compressedSuffixArray->valueBitWidth - bitsWrittenOnFirstWriteAction;

    saValue >>= bitsWrittenOnFirstWriteAction;
    while (bitsRemaining > 0) {
      compressedSuffixArray->values[offset.byteOffset++] = saValue;
      saValue >>= 8;
      bitsRemaining -= 8;
    }
  }

  // reallocate the suffix array to save space.
  assert(fullSa != NULL);
  compressedSuffixArray->values =
      realloc(compressedSuffixArray->values,
              compressedSuffixArray->compressedByteLength);

  if (compressedSuffixArray->values == NULL) {
    return AwFmAllocationFailure;
  }

  return AwFmSuccess;
}

size_t awFmGetValueFromCompressedSuffixArray(
    const struct AwFmCompressedSuffixArray *suffixArray,
    size_t positionInArray) {
  struct AwFmSuffixArrayOffset offset = awFmGetOffsetIntoSuffixArrayByteArray(
      suffixArray->valueBitWidth, positionInArray);
  // memcpy the data containing the value into a uint64_t so we can shift the
  // data.
  uint64_t buffer;
  memcpy(&buffer, &suffixArray->values[offset.byteOffset], 8);
  const uint64_t bitmask = (1ULL << suffixArray->valueBitWidth) - 1;

  // shifts the 8-byte int into place. if this leaves some bits off the end,
  // we'll handle it in a bit.
  buffer >>= offset.bitOffset;

  // if we can  ignore the last byte of the SA memcpy'd buffer, thanks strategy
  // pattern!
  if (__builtin_expect(suffixArray->valueBitWidth > 57, 0)) {
    uint64_t lastByte = suffixArray->values[offset.byteOffset + 8];
    // shift the lastByte first by 1, then the rest of the way. This needs to be
    // done because shifting an int by its bitLength is undefined behavior in
    // GCC.
    lastByte <<= 1;
    lastByte <<= (63 - offset.bitOffset);
    buffer |= lastByte;
  }

  return buffer & bitmask;
}

inline size_t awFmGetSampledSuffixArrayLength(uint64_t bwtLength,
                                              uint64_t compressionRatio) {
  return (bwtLength + compressionRatio - 1) / compressionRatio;
}

enum AwFmReturnCode
awFmReadPositionsFromSuffixArray(const struct AwFmIndex *_RESTRICT_ const index,
                                 uint64_t *_RESTRICT_ const positionArray,
                                 const size_t positionArrayLength) {

  if (index->config.keepSuffixArrayInMemory) {
    for (size_t i = 0; i < positionArrayLength; i++) {
      size_t indexInSuffixArray =
          positionArray[i] / index->config.suffixArrayCompressionRatio;
      size_t position = awFmGetValueFromCompressedSuffixArray(
          &index->suffixArray, indexInSuffixArray);
      positionArray[i] = position;
    }

    return AwFmSuccess;
  } else {
    for (size_t i = 0; i < positionArrayLength; i++) {
      size_t indexInSuffixArray =
          positionArray[i] / index->config.suffixArrayCompressionRatio;
      enum AwFmReturnCode rc = awFmGetSuffixArrayValueFromFile(
          index, indexInSuffixArray, &positionArray[i]);
      if (rc != AwFmSuccess) {
        return rc;
      }
    }

    return AwFmFileReadOkay;
  }
}

enum AwFmReturnCode awFmSuffixArrayReadPositionParallel(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmBacktrace *_RESTRICT_ const backtracePtr) {

  if (__builtin_expect(index->config.keepSuffixArrayInMemory, 1)) {
    uint64_t suffixArrayPosition =
        backtracePtr->position / index->config.suffixArrayCompressionRatio;

    size_t saValue = awFmGetValueFromCompressedSuffixArray(&index->suffixArray,
                                                           suffixArrayPosition);
    backtracePtr->position =
        (saValue + backtracePtr->offset) % index->bwtLength;
    return AwFmSuccess;
  } else {
    uint64_t suffixArrayPosition =
        backtracePtr->position / index->config.suffixArrayCompressionRatio;
    size_t saValue;
    enum AwFmReturnCode rc =
        awFmGetSuffixArrayValueFromFile(index, suffixArrayPosition, &saValue);

    backtracePtr->position =
        (saValue + backtracePtr->offset) % index->bwtLength;
    return rc;
  }
}
