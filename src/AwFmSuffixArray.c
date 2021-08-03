#include "AwFmSuffixArray.h"
#include <assert.h>

//simple log2 ceiling implementation, thanks builtin clzll!
uint8_t log2Ceil(const uint64_t a){
    return 64 - __builtin_clzll(a);
}

//the number of bits required to store the suffix array values is just log2ceil of the highest val.
uint8_t awFmComputeSuffixArrayValueMinWidth(const size_t saLength){
  return log2Ceil(saLength);
}

size_t awFmComputeCompressedSaSizeInBytes(size_t saLength, uint8_t samplingRatio){
  uint8_t minimumBitWidth = log2Ceil(saLength);
  uint64_t sampledSaLength = saLength / ((size_t)samplingRatio);
  uint8_t saLengthRemainder       = sampledSaLength % 8;
  //if each int is 'minWidth' wide, then each 8 values should require 'minWidth' bytes.
  //this division should also round down, so it automatically removes the remainder.
  size_t totalBytes               = (sampledSaLength/8) * minimumBitWidth;
  size_t remainingBits            = saLengthRemainder * minimumBitWidth;
  totalBytes += (remainingBits + 7) / 8;  //ceiling, since we need AT LEAST that many bits

  return totalBytes;
}


//note: fullSA MUST be dynamically allocated, this reallocs the array and claims ownership.
//the reallocated array will be freed when awFmDeallocIndex() is called.
 enum AwFmReturnCode awFmInitCompressedSuffixArray(uint64_t *fullSa, size_t saLength,
  struct AwFmCompressedSuffixArray *compressedSuffixArray, uint8_t samplingRatio){

  if(fullSa == NULL || compressedSuffixArray == NULL){
    return AwFmNullPtrError;
  }
  uint8_t minimumBitWidth = log2Ceil(saLength);
  size_t compressedSaByteSize = awFmComputeCompressedSaSizeInBytes(saLength, samplingRatio);

  compressedSuffixArray->values = (uint8_t*)fullSa;
  compressedSuffixArray->valueBitWidth = minimumBitWidth;
  compressedSuffixArray->compressedByteLength = compressedSaByteSize;

  const size_t numSaSamples = saLength / samplingRatio;
  for(size_t i = 0; i < numSaSamples; i++){
    uint64_t saValueBuffer = fullSa[i * samplingRatio];
    //the following finds the byte and bit that the element will start on.
    //this is written weirdly to avoid integer overflow for large saLengths near 2^64.
    size_t alignedByteOffset = (i / 8) * minimumBitWidth; //every 8 positions consumes 'bitWidth' number of bytes.
    size_t numEndingBits = (i % 8) * minimumBitWidth;
    size_t bytePosition = alignedByteOffset + (numEndingBits / 8);
    size_t bitPosition = numEndingBits % 8;

    uint8_t byteToWrite = (saValueBuffer << bitPosition) & 0xFF;
    uint8_t bitmask = (1ULL << bitPosition) - 1; //zero out the bits we're going to set, since it's in-place
    compressedSuffixArray->values[bytePosition]    &=  bitmask;
    compressedSuffixArray->values[bytePosition++]  |= byteToWrite;

    int16_t bitsRemaining = compressedSuffixArray->valueBitWidth - (8- bitPosition);
    saValueBuffer >>= (8 - bitPosition);
    while(bitsRemaining > 0){
      compressedSuffixArray->values[bytePosition++] = saValueBuffer;
      saValueBuffer >>= 8;
      bitsRemaining -= 8;
    }
  }

  //reallocate the suffix array to save space.
  assert(fullSa != NULL);
  compressedSuffixArray->values = realloc(fullSa, compressedSuffixArray->compressedByteLength);

  if(fullSa == NULL){
    printf("FULL SA WAS NULL after realloc\n");
  }
  if(compressedSuffixArray->values == NULL){
    return AwFmAllocationFailure;
  }
  return AwFmSuccess;
}


size_t awFmGetValueFromCompressedSuffixArray(struct AwFmCompressedSuffixArray suffixArray, size_t positionInArray){
  size_t alignedByteOffset = (positionInArray / 8) * suffixArray.valueBitWidth;
  size_t numEndingBits = (positionInArray % 8) * suffixArray.valueBitWidth;
  size_t bytePosition = alignedByteOffset + (numEndingBits / 8);
  size_t bitPosition = numEndingBits % 8;

  size_t value = 0;
  value = suffixArray.values[bytePosition++] >> bitPosition;
  uint8_t bitsFromFirstByte = 8 - bitPosition;
  for(int_fast16_t bitOffset = bitsFromFirstByte;
    bitOffset < suffixArray.valueBitWidth; bitOffset += 8){
      value |= ((uint64_t)suffixArray.values[bytePosition++]) << bitOffset;
  }

  //create and apply a bitmask to remove any extra bits at the end.
  uint64_t bitmask =  (1 << suffixArray.valueBitWidth) - 1;
  return value & bitmask;
}
