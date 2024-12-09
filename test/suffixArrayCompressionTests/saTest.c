#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../src/AwFmSuffixArray.h"
#include "../test.h"

char messageBuffer[1024];

void generateFabricatedSuffixArray(uint64_t length, uint64_t *buffer) {
  for (uint64_t i = 0; i < length; i++) {
    buffer[i] = i;
  }
  for (uint64_t pivotPoint = length - 1; pivotPoint > 0; pivotPoint--) {
    // generate a 64 bit random number.
    uint64_t randChoice = rand();
    randChoice = ((randChoice << 32ULL) | rand()) % pivotPoint;
    uint64_t tmp = buffer[pivotPoint];
    buffer[pivotPoint] = buffer[randChoice];
    buffer[randChoice] = tmp;
  }
}

void testSuffixArrayCompressionStaticLengths() {
  for (uint64_t saLength = 4; saLength < 1024; saLength++) {
    uint64_t *buffer = malloc(saLength * sizeof(uint64_t));
    uint64_t *referenceBuffer = malloc(saLength * sizeof(uint64_t));
    generateFabricatedSuffixArray(saLength, buffer);
    memcpy(referenceBuffer, buffer, saLength * sizeof(uint64_t));

    struct AwFmCompressedSuffixArray compressedSa;
    enum AwFmReturnCode returnCode =
        awFmInitCompressedSuffixArray(buffer, saLength, &compressedSa, true);
    if (returnCode != AwFmSuccess) {
      testAssertString(false,
                       "awFmInitSuffixArray did not return AwFmSuccess.\n");
      printf("initSuffixArray returned error code %d\n", returnCode);
    }

    uint64_t necessaryBitWidth = 0;
    while ((1LL << (++necessaryBitWidth)) < saLength)
      ;

    sprintf(messageBuffer,
            "sa length %zu expected bit width %zu, but SA said bit width %d.\n",
            saLength, necessaryBitWidth, compressedSa.valueBitWidth);
    testAssertString(necessaryBitWidth == compressedSa.valueBitWidth,
                     messageBuffer);

    for (uint64_t i = 1; i < saLength; i++) {
      uint64_t restoredValue =
          awFmGetValueFromCompressedSuffixArray(&compressedSa, i);
      uint64_t referenceValue = referenceBuffer[i];
      sprintf(messageBuffer,
              "at index %zu  of length %zu SA, compressed value %zu did not "
              "match reference val %zu.\n",
              i, saLength, restoredValue, referenceValue);
      testAssertString(restoredValue == referenceValue, messageBuffer);
    }
    free(compressedSa.values);
    free(referenceBuffer);
  }
}

void testSuffixArrayCompressionRandomLengths() {
  uint64_t *buffer = NULL;
  uint64_t *referenceBuffer = NULL;
  for (uint64_t testNum = 0; testNum < 80; testNum++) {
    uint64_t suffixArrayLength =
        (rand() >> 12); // only testing to 2^20 because of memory limitations on
                        // my machine.
    printf("test num %zu, SA len %zu.\n", testNum, suffixArrayLength);
    buffer = realloc(buffer, suffixArrayLength * sizeof(uint64_t));
    referenceBuffer =
        realloc(referenceBuffer, suffixArrayLength * sizeof(uint64_t));
    generateFabricatedSuffixArray(suffixArrayLength, buffer);
    memcpy(referenceBuffer, buffer, suffixArrayLength * sizeof(uint64_t));

    struct AwFmCompressedSuffixArray compressedSa;
    enum AwFmReturnCode returnCode = awFmInitCompressedSuffixArray(
        buffer, suffixArrayLength, &compressedSa, true);
    if (returnCode != AwFmSuccess) {
      testAssertString(false,
                       "awFmInitSuffixArray did not return AwFmSuccess.\n");
      printf("initSuffixArray returned error code %d\n", returnCode);
    }

    // re-set the buffer pointer, so we can realloc it again
    buffer = (uint64_t *)compressedSa.values;
    uint64_t necessaryBitWidth = 0;
    while ((1LL << (++necessaryBitWidth)) <= suffixArrayLength)
      ;

    sprintf(messageBuffer,
            "sa length %zu expected bit width %zu, but SA said bit width %d.\n",
            suffixArrayLength, necessaryBitWidth, compressedSa.valueBitWidth);
    testAssertString(necessaryBitWidth == compressedSa.valueBitWidth,
                     messageBuffer);

    for (uint64_t i = 0; i < suffixArrayLength; i++) {
      uint64_t restoredValue =
          awFmGetValueFromCompressedSuffixArray(&compressedSa, i);
      sprintf(messageBuffer,
              "at index %zu  of length %zu SA, compressed value %zu did not "
              "match reference val %zu.\n",
              i, suffixArrayLength, restoredValue, referenceBuffer[i]);
      testAssertString(restoredValue == referenceBuffer[i], messageBuffer);
    }
  }

  free(buffer);
  free(referenceBuffer);
}

int main(int argc, char **argv) {
  srand(time(NULL));
  testSuffixArrayCompressionStaticLengths();
  testSuffixArrayCompressionRandomLengths();
}
