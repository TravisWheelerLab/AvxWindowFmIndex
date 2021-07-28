#include "AwFmSuffixArray.h"
#include "AwFmIndex.h"
#include <stdlib.h>
#include <string.h>

uint8_t log2Ceil(const uint64_t a){
    return 64 - __builtin_clzll(a);
}

size_t computeCompressedSaSizeInBytes(size_t saLength, uint8_t minWidth){
  uint8_t saLengthRemainder       = saLength % 8;
  //if each int is 'minWidth' wide, then each 8 values should require 'minWidth' bytes.
  //this division should also round down, so it automatically removes the remainder.
  size_t totalBytes               = (saLength/8) * minWidth;
  size_t remainingBits            = saLengthRemainder * minWidth;
  totalBytes += (remainingBits + 7) / 8;  //ceiling, since we need AT LEAST that many bits

  return totalBytes;
}


//need to check for null Sa
struct AwFmCompressedSuffixArray awFmInitSuffixArray(uint64_t *fullSa, size_t saLength){
  uint8_t minimumBitWidth = log2Ceil(saLength);
  size_t compressedSaByteSize = computeCompressedSaSizeInBytes(saLength, minimumBitWidth);
  struct AwFmCompressedSuffixArray compressedSuffixArray;

  compressedSuffixArray.values = (uint8_t*)fullSa;
  compressedSuffixArray.length = saLength;
  compressedSuffixArray.valueBitWidth = minimumBitWidth;
  compressedSuffixArray.compressedByteLength = compressedSaByteSize;

  for(size_t i = 0; i < saLength; i++){
    uint64_t saValueBuffer = fullSa[i];
    //the following finds the byte and bit that the element will start on.
    //this is written weirdly to avoid integer overflow for large saLengths near 2^64.
    //ever
    size_t alignedByteOffset = (i / 8) * minimumBitWidth; //every 8 positions consumes 'bitWidth' number of bytes.
    size_t numEndingBits = (i % 8) * minimumBitWidth;
    size_t bytePosition = alignedByteOffset + (numEndingBits / 8);
    size_t bitPosition = numEndingBits % 8;

    uint8_t byteToWrite = (saValueBuffer << bitPosition) & 0xFF;
    uint8_t bitmask = (1ULL << bitPosition) - 1; //zero out the bits we're going to set, since it's in-place
    compressedSuffixArray.values[bytePosition]    &=  bitmask;
    compressedSuffixArray.values[bytePosition++]  |= byteToWrite;

    int16_t bitsRemaining = compressedSuffixArray.valueBitWidth - (8- bitPosition);
    saValueBuffer >>= (8 - bitPosition);
    while(bitsRemaining > 0){
      compressedSuffixArray.values[bytePosition++] = saValueBuffer;
      saValueBuffer >>= 8;
      bitsRemaining -= 8;
    }
  }

  return compressedSuffixArray;
}


size_t awFmGetValueFromCompressedSuffixArray(struct AwFmCompressedSuffixArray suffixArray, size_t valuePosition){
  size_t alignedByteOffset = (valuePosition / 8) * suffixArray.valueBitWidth;
  size_t numEndingBits = (valuePosition % 8) * suffixArray.valueBitWidth;
  size_t bytePosition = alignedByteOffset + (numEndingBits / 8);
  size_t bitPosition = numEndingBits % 8;

  size_t value = 0;
  value = suffixArray.values[bytePosition++] >> bitPosition;
  uint8_t bitsFromFirstByte = 8 - bitPosition;
  // uint8_t bitsRemaining = suffixArray.valueBitWidth - bitsFromFirstByte;
  for(int_fast16_t bitOffset = bitsFromFirstByte;
    bitOffset < suffixArray.valueBitWidth; bitOffset += 8){
      value |= ((uint64_t)suffixArray.values[bytePosition++]) << bitOffset;
  }

  //create and apply a bitmask to remove any extra bits at the end.
  uint64_t bitmask =  (1 << suffixArray.valueBitWidth) - 1;
  return value & bitmask;

}
