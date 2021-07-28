#include "../../src/AwFmSuffixArray.h"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

char messageBuffer[1024];

void generateFabricatedSuffixArray(uint64_t length, uint64_t *buffer){
  for(uint64_t i = 0; i < length; i++){
    buffer[i] = i;
  }
  for(uint64_t pivotPoint = length-1; pivotPoint > 0; pivotPoint--){
    //generate a 64 bit random number.
    uint64_t randChoice = rand();
    randChoice = ((randChoice << 32ULL) | rand()) % pivotPoint;
    uint64_t tmp = buffer[pivotPoint];
    buffer[pivotPoint] = buffer[randChoice];
    buffer[randChoice] = tmp;
  }
}


void testSuffixArrayCompressionStaticLengths(){
  uint64_t buffer[512];
  uint64_t referenceBuffer[512];
    for(uint64_t saLength = 4; saLength  < 512; saLength++){
    generateFabricatedSuffixArray(saLength, buffer);
    memcpy(referenceBuffer, buffer, saLength * sizeof(uint64_t));


    struct AwFmCompressedSuffixArray  compressedSa = awFmInitSuffixArray(buffer, saLength);



    uint64_t necessaryBitWidth = 0;
    while((1LL << (++necessaryBitWidth)) <= saLength);

    sprintf(messageBuffer, "sa length %zu expected bit width %zu, but SA said bit width %d.\n",
      saLength, necessaryBitWidth, compressedSa.valueBitWidth);
    testAssertString(necessaryBitWidth == compressedSa.valueBitWidth, messageBuffer);

    for(uint64_t i = 0; i < saLength; i++){
      uint64_t restoredValue = awFmGetValueFromCompressedSuffixArray(compressedSa, i);
      sprintf(messageBuffer, "at index %zu  of length %zu SA, compressed value %zu did not match reference val %zu.\n",
        i, saLength, restoredValue, referenceBuffer[i]);
      testAssertString(restoredValue == referenceBuffer[i], messageBuffer);

    }
  }
}


void testSuffixArrayCompressionRandomLengths(){
  uint64_t *buffer = NULL;
  uint64_t *referenceBuffer = NULL;
  for(uint64_t testNum = 0; testNum < 40; testNum++){
    uint64_t suffixArrayLength = (rand() >> 12);  //only testing to 2^20 because of memory limitations on my machine.
    buffer          = realloc(buffer, suffixArrayLength * sizeof(uint64_t));
    referenceBuffer = realloc(referenceBuffer, suffixArrayLength * sizeof(uint64_t));
    generateFabricatedSuffixArray(suffixArrayLength, buffer);
    memcpy(referenceBuffer, buffer, suffixArrayLength * sizeof(uint64_t));
    struct AwFmCompressedSuffixArray  compressedSa = awFmInitSuffixArray(buffer, suffixArrayLength);

    uint64_t necessaryBitWidth = 0;
    while((1LL << (++necessaryBitWidth)) <= suffixArrayLength);

    sprintf(messageBuffer, "sa length %zu expected bit width %zu, but SA said bit width %d.\n",
      suffixArrayLength, necessaryBitWidth, compressedSa.valueBitWidth);
    testAssertString(necessaryBitWidth == compressedSa.valueBitWidth, messageBuffer);

    for(uint64_t i = 0; i < suffixArrayLength; i++){
      uint64_t restoredValue = awFmGetValueFromCompressedSuffixArray(compressedSa, i);
      sprintf(messageBuffer, "at index %zu  of length %zu SA, compressed value %zu did not match reference val %zu.\n",
        i, suffixArrayLength,restoredValue, referenceBuffer[i]);
      testAssertString(restoredValue == referenceBuffer[i], messageBuffer);
    }

  }

  free(buffer);
  free(referenceBuffer);
}

int main(int argc, char **argv){
  srand(2);
  //srand(time(NULL));
  testSuffixArrayCompressionStaticLengths();
  testSuffixArrayCompressionRandomLengths();
}
