#include "AwFmLetter.h"
#include <stdlib.h>
#include <time.h>


/*
 * Function:  awFmAsciiLetterToLetterIndexSrand
 * --------------------
 * NOTICE: this function is exactly the same as awFmAsciiLetterToLetterIndex,
 * but calls srand to seed the RNG for ambiguity code resolution.
 *
 * Transforms an ascii-encoded character into a frequency-encoded amino acid value.
 * These frequency encoded values are used internally in AwFm to better optimize memory usage and aid performance.
 *  CAUTION: characters that are not ASCII letters may alias to represent an amino acid. Only letters a-z, A-Z should be given.
 *
 *  Inputs:
 *    asciiLetter: ascii-encoded letter representing an amino acid, or an amino acid ambiguity character (b, z, or x)
 *
 *  Returns:
 *    Frequency-encoded value representing the amino acid.
 *
 */
uint8_t awFmAsciiLetterToLetterIndexSrand(uint8_t asciiLetter, unsigned seed){
  srand(seed);
  return awFmAsciiLetterToLetterIndex(asciiLetter);
}


/*
 * Function:  awFmAsciiLetterToLetterIndex
 * --------------------
 * Transforms an ascii-encoded character into a frequency-encoded amino acid value.
 * These frequency encoded values are used internally in AwFm to better optimize memory usage and aid performance.
 *  CAUTION: characters that are not ASCII letters may alias to represent an amino acid. Only letters a-z, A-Z should be given.
 *  NOTICE: this function resolves ambiguity codes via the rand() function.
 *  In order to properly randomize the resolved amino acids, make sure rand is seeded with srand,
 *  or use the awFmAsciiLetterToLetterIndexSrand function for this purpose.
 *
 *  Inputs:
 *    asciiLetter: ascii-encoded letter representing an amino acid, or an amino acid ambiguity character (b, z, or x)
 *
 *  Returns:
 *    Frequency-encoded value representing the amino acid.
 *
 */
uint8_t awFmAsciiLetterToLetterIndex(uint8_t asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  const char charBitmask = 0x1F;
  asciiLetter = asciiLetter & charBitmask;

  //frequency indexed letter lookup table. values of 20 are illegal characters, and should not match to anything.
  //by frequency, L,A,G,V,E,S,I,K,R,D,T,P,N,Q,F,Y,M,H,C,W
  const uint8_t letterOffsets[32] = {
    //letters A - Z
    20, 1,20,18, 9, 4,14, 2,
    17, 6,20, 7, 0,16,12,20,
    11,13, 8, 5,10,20, 3,19,
    20,15,20,20,20,20,20,20
  };

  if(asciiLetter == ('B' & charBitmask)){
    asciiLetter = (rand() & 1)? 'D'& charBitmask: 'N'& charBitmask;
  }
  else if(asciiLetter == ('Z' & charBitmask)){
    asciiLetter = (rand() & 1)? 'E'& charBitmask: 'Q'& charBitmask;
  }

  uint8_t frequencyEncodedLetter = letterOffsets[asciiLetter];
  if(__builtin_expect(frequencyEncodedLetter == 20, 0)){
    return rand() % 20;
  }
  else{
    return frequencyEncodedLetter;

  }

}


/*
 * Function:  awFmAsciiLetterToCompressedVectorFormat
 * --------------------
 * Transforms an ascii-encoded character into a set of 5 bits that will represent the amino acid in the strided AVX vector.
 *  CAUTION: characters that are not ASCII letters may alias to represent an amino acid. Only letters a-z, A-Z should be given.
 *  Inputs:
 *    asciiLetter: ascii-encoded letter representing an amino acid, or an amino acid ambiguity character (b, z, or x)
 *
 *  Returns:
 *    Value representing the compressed vector format of the amino acid.
 *
 */
uint8_t awFmAsciiLetterToCompressedVectorFormat(const uint8_t asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  const char charBitmask = 0x1F;
  const uint8_t compressedVectorLetters[32] = {
    0x00, 0x1A, 0x00, 0x04, 0x06, 0x15, 0x1B, 0x16,
    0x02, 0x03, 0x00, 0x05, 0x1C, 0x01, 0x1E, 0x00,
    0x0C, 0x1D, 0x09, 0x13, 0x0A, 0x00, 0x19, 0x08,
    0x00, 0x17, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  };

  return compressedVectorLetters[asciiLetter & charBitmask];
}


uint8_t awFmCompressedVectorLetterToLetterIndex(const uint8_t compressedVectorLetter){
  const uint8_t frequencyIndexLetters[32] = {
    20, 16, 17,  6, 18,  7,  9, 20,
    19,  8, 10, 20, 11, 20, 20, 20,
    20, 20, 20,  5, 20,  4,  2, 15,
    20,  3,  1, 14,  0, 13, 12, 20
  };

  return frequencyIndexLetters[compressedVectorLetter];
}
