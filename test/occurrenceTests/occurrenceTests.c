#include "../../AwFmIndex.h"
#include "../../AwFmOccurrence.h"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>

void popcountTestSuite(void);

char buffer[2048];

int main(int argc, char **argv){
  srand(time(NULL));
  popcountTestSuite();


  printf("occurrence function testing finished.\n");
}


void setVectorBytes(uint8_t *restrict const vector, const uint8_t byteValue){
  for(uint8_t i = 0; i < 32; i++){
    vector[i] = byteValue;
  }
}

uint16_t setVectorRandBits(uint8_t *restrict const vector){
  uint8_t bitsSet = 0;
  for(uint8_t byteIndex = 0; byteIndex < 32; byteIndex++){
    vector[byteIndex] = 0;

    for(uint8_t bitIndex = 0; bitIndex < 8; bitIndex++){
      uint8_t bitValue = rand() & 1;
      vector[byteIndex] |= bitValue << bitIndex;
      bitsSet += bitValue;
    }
  }
  return bitsSet;
}

//tests awFmVectorPopcount function
void popcountTestSuite(void){
  uint8_t vectorBytes[32];

  //test all 0s
  __m256i vector = vector = _mm256_set1_epi8(0x00);
  uint16_t popcount = awFmVectorPopcount(vector);
  sprintf(buffer, "vector popcount of all zeros should be zero, returned %d.", popcount);
  testAssertString(popcount == 0, buffer);

  //test all 1s
  setVectorBytes(vectorBytes, 0xFF);
  vector = _mm256_set1_epi8(0xff);
  popcount = awFmVectorPopcount(vector);
  sprintf(buffer, "vector popcount of all zeros should be all ones (256), returned %d.", popcount);
  testAssertString(popcount == 256, buffer);

  //test all 0x7fs
  setVectorBytes(vectorBytes, 0x7f);
  vector = _mm256_loadu_si256((__m256i*)vectorBytes);
  popcount = awFmVectorPopcount(vector);
  sprintf(buffer, "vector popcount of all 0x7f should be (224), returned %d.", popcount);
  testAssertString(popcount == 224, buffer);

  //test all 0xf7
  setVectorBytes(vectorBytes, 0xf7);
  vector = _mm256_loadu_si256((__m256i*)vectorBytes);
  popcount = awFmVectorPopcount(vector);
  sprintf(buffer, "vector popcount of all 0xf7 should be (224), returned %d.", popcount);
  testAssertString(popcount == 224, buffer);

  //test all 77s
  setVectorBytes(vectorBytes, 0x77);
  vector = _mm256_loadu_si256((__m256i*)vectorBytes);
  popcount = awFmVectorPopcount(vector);
  sprintf(buffer, "vector popcount of all 0x77 should be (192), returned %d.", popcount);
  testAssertString(popcount == 192, buffer);



//test all 77s
setVectorBytes(vectorBytes, 0x0f);
vector = _mm256_loadu_si256((__m256i*)vectorBytes);
popcount = awFmVectorPopcount(vector);
sprintf(buffer, "vector popcount of all 0x0f should be (128), returned %d.", popcount);
testAssertString(popcount == 128, buffer);

//test all 77s
setVectorBytes(vectorBytes, 0xf0);
vector = _mm256_loadu_si256((__m256i*)vectorBytes);
popcount = awFmVectorPopcount(vector);
sprintf(buffer, "vector popcount of all 0x0f should be (128), returned %d.", popcount);
testAssertString(popcount == 128, buffer);

  const uint16_t numRandTests = 20000;
  clock_t startTime = clock();
  for(uint16_t testNum = 0; testNum < numRandTests; testNum++){
    uint16_t bitsSet = setVectorRandBits(vectorBytes);
    vector = _mm256_loadu_si256((__m256i*)vectorBytes);
    popcount = awFmVectorPopcount(vector);
    sprintf(buffer, "vector popcount of all zeros should be %d, returned %d.", bitsSet, popcount);
    testAssertString(popcount == bitsSet, buffer);
  }
  clock_t endTime = clock();
  size_t timeTaken = endTime - startTime;
  printf("time: %zu\n", timeTaken);

}