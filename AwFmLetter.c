#include "AwFmLetter.h"
#include <stdlib.h>



/*
 * Function:  awFmAsciiLetterToLetterIndex
 * --------------------
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
uint8_t awFmAsciiLetterToLetterIndex(char asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  const char charBitmask = 0x1F;
  asciiLetter = asciiLetter & charBitmask;

  //frequency indexed letter lookup table. values of 21 are illegal characters, and should not match to anything.
  //by frequency, L,A,G,V,E,S,I,K,R,D,T,P,N,Q,F,Y,M,H,C,W
  const letter_t letterOffsets[32] = {
    //letters A - Z
    21, 1,21,18, 9, 4,14, 2,
    17, 6,21, 7, 0,16,12,21,
    11,13, 8, 5,10,21, 3,19,
    21,15,21,21,21,21,21,21
  };

  if(asciiLetter)
  if(asciiLetter == ('B' & charBitmask)){
    asciiLetter = (rand() & 1)? 'D'& charBitmask: 'N'& charBitmask;
  }
  else if(asciiLetter == ('Z' & charBitmask)){
    asciiLetter = (rand() & 1)? 'E'& charBitmask: 'Q'& charBitmask;
  }
  else if(asciiLetter == 'X' & charBitmask){
    //X is the ambiguity character for completely unknown, randomize the letter.
    asciiLetter = (rand() % 19) + 1;
    //if we randomized to another illegal character, add 1 to get a legal char.
    if(letterOffsets[asciiLetter] == 21){
      asciiLetter += 1;
    }
  }

  return letterOffsets[asciiLetter];
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
uint8_t awFmAsciiLetterToCompressedVectorFormat(const char asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  const char charBitmask = 0x1F;
  const letter_t compressedVectorLetters[32] = {
    0x00, 0x1A, 0x00, 0x04, 0x06, 0x15, 0x1B, 0x16,
    0x02, 0x03, 0x00, 0x05, 0x1C, 0x01, 0x1E, 0x00,
    0x0C, 0x1D, 0x09, 0x13, 0x0A, 0x00, 0x19, 0x08,
    0x00, 0x17, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  };

  return compressedVectorLetters[asciiValue & charBitmask];
}
