#include "../../AwFmLetter.h"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <ctype.h>


struct ambiguityCode{
  char code;
  char *possibleChars;
};

const uint8_t numEncodingAminos = 20;
uint8_t encodingAminos[21] = {'a','c','d','e','f',
                              'g','h','i','k','l',
                              'm','n','p','q','r',
                              's','t','v','w','y',0};

uint8_t encodingAminosUpperCase[21] = {'A','C','D','E','F',
                              'G','H','I','K','L',
                              'M','N','P','Q','R',
                              'S','T','V','W','Y',0};


char buffer[2048];

void testNucleotideAscii(){
  testAssert(awFmAsciiNucleotideToLetterIndex('a') == 0);
  testAssert(awFmAsciiNucleotideToLetterIndex('A') == 0);
  testAssert(awFmAsciiNucleotideToLetterIndex('c') == 1);
  testAssert(awFmAsciiNucleotideToLetterIndex('C') == 1);
  testAssert(awFmAsciiNucleotideToLetterIndex('g') == 2);
  testAssert(awFmAsciiNucleotideToLetterIndex('G') == 2);
  testAssert(awFmAsciiNucleotideToLetterIndex('t') == 3);
  testAssert(awFmAsciiNucleotideToLetterIndex('T') == 3);


  testAssert(awFmNucleotideLetterIndexToAscii(0) == 'A');
  testAssert(awFmNucleotideLetterIndexToAscii(1) == 'C');
  testAssert(awFmNucleotideLetterIndexToAscii(2) == 'G');
  testAssert(awFmNucleotideLetterIndexToAscii(3) == 'T');
}


void testNucleotideAmbiguityCodes(){


  struct ambiguityCode codes[32] = {
    {'a', "A"},
    {'A', "A"},
    {'c', "C"},
    {'C', "C"},
    {'g', "G"},
    {'G', "G"},
    {'t', "T"},
    {'T', "T"},
    {'y', "CT"},
    {'Y', "CT"},
    {'r', "AG"},
    {'R', "AG"},
    {'w', "AT"},
    {'W', "AT"},
    {'s', "GC"},
    {'S', "GC"},
    {'k', "TG"},
    {'K', "TG"},
    {'m', "CA"},
    {'M', "CA"},
    {'d', "AGT"},
    {'D', "AGT"},
    {'v', "AGC"},
    {'V', "AGC"},
    {'h', "ACT"},
    {'H', "ACT"},
    {'b', "CGT"},
    {'B', "CGT"},
    {'x', "ACGT"},
    {'X', "ACGT"},
    {'n', "ACGT"},
    {'N', "ACGT"},
};

  for(uint8_t i = 0; i < 32; i++){
    for(uint8_t testNum = 0; testNum < 100; testNum++){
      uint8_t index = awFmAsciiNucleotideToLetterIndex(codes[i].code);
      char letterOut = awFmNucleotideLetterIndexToAscii(index);
      char charBuf[2] = {0,0};
      charBuf[0] = letterOut;
      sprintf(buffer, "test fail, %s not found in possible chars %s\n", charBuf, codes[i].possibleChars);
      testAssertString(strchr(codes[i].possibleChars, letterOut) != NULL, buffer);
    }
  }
}


void testAminoIndex(){
  for(uint8_t i = 0; i < 20; i++){
    char amino = encodingAminos[i];
    char upperCaseAmino = toupper(amino);
    uint8_t index = awFmAsciiAminoAcidToLetterIndex(amino);
    uint8_t upperIndex = awFmAsciiAminoAcidToLetterIndex(upperCaseAmino);
    sprintf(buffer, "test fail, amino %c expected to generate %d, but generated %d", amino, i, index);
    testAssertString(index == i, buffer);
    sprintf(buffer, "test fail, amino %c expected to generate %d, but generated %d", upperCaseAmino, i, upperIndex);
    testAssertString(index == i, buffer);
  }
}

void testAminoAmbiguityCodes(){
  struct ambiguityCode codes[6] = {
    {'b', "DN"},
    {'B', "DN"},
    {'z', "EQ"},
    {'Z', "EQ"},
    {'x', (char*)encodingAminosUpperCase},
    {'X', (char*)encodingAminosUpperCase}
  };

  for(uint8_t i = 0; i < 6; i++){
    for(uint8_t testNum = 0; testNum < 100; testNum++){
      uint8_t index = awFmAsciiAminoAcidToLetterIndex(codes[i].code);
      char letterOut = awFmAminoAcidLetterIndexToAscii(index);
      char charBuf[2] = {0,0};
      charBuf[0] = letterOut;
      sprintf(buffer, "test fail, %s not found in possible chars %s\n", charBuf, codes[i].possibleChars);
      testAssertString(strchr(codes[i].possibleChars, letterOut) != NULL, buffer);
    }
  }
}

void testCompressedAminos(){

  const uint8_t compressedEncodings[20] = {
    0x0C, 0x17, 0x03, 0x06, 0x1E, 0x1A, 0x1B, 0x19,
    0x15, 0x1C, 0x1D, 0x08, 0x09, 0x04, 0x13, 0x0A,
    0x05, 0x16, 0x01, 0x02};

  for(uint8_t i = 0; i < 20; i++){
    uint8_t asAscii = encodingAminos[i];
    uint8_t asAsciiUpper = encodingAminosUpperCase[i];

    uint8_t encoding = awFmAminoAcidAsciiLetterToCompressedVectorFormat(asAscii);
    uint8_t encodingUpper = awFmAminoAcidAsciiLetterToCompressedVectorFormat(asAsciiUpper);

    char charBuf[2] = {0,0};
    charBuf[0] = asAscii;
    sprintf(buffer, "test fail, lower case ascii %s should generate encoding 0x%.2X but generated 0x%.8X.\n",
      charBuf, compressedEncodings[i], encoding);
    testAssertString(compressedEncodings[i] == encoding, buffer);

    charBuf[0] = asAsciiUpper;
    sprintf(buffer, "test fail, upper case ascii %s should generate encoding 0x%.2X but generated 0x%.8X.\n",
      charBuf, compressedEncodings[i], encodingUpper);
    testAssertString(compressedEncodings[i] == encodingUpper, buffer);

    uint8_t backToIndex = awFmAminoAcidCompressedVectorToLetterIndex(encoding);
    sprintf(buffer, "test fail, encoding 0x%.2X for letter index %d generated incorrect index %d.\n",
      encodingUpper, i, backToIndex);
    testAssertString(backToIndex == i, buffer);


    uint8_t backToAscii = awFmAminoAcidCompressedVectorToAscii(encoding);
    sprintf(buffer, "test fail, encoding 0x%.2X for letter index %d created character %c instead of expected %c.\n",
      encodingUpper, i, asAsciiUpper, toupper(backToAscii));
    testAssertString(toupper(backToAscii) == asAsciiUpper, buffer);

  }
}



int main(int argc, char **argv){
  srand(time(NULL));
  testNucleotideAscii();
  testNucleotideAmbiguityCodes();
  testAminoIndex();
  testAminoAmbiguityCodes();
  testCompressedAminos();
  printf("letter tests finished\n");



}
