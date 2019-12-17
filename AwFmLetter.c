#include "AwFmLetter.h"
#include "AwFmIndex.h"
#include <stdlib.h>
#include <time.h>


uint8_t awFmLetterToLetterIndex(const uint8_t asciiLetter, const enum AwFmAlphabetType alphabet){
  return (alphabet == AwFmAlphabetAminoAcid)?
    awFmAsciiAminoAcidToLetterIndex(asciiLetter):
    awFmAsciiNucleotideToLetterIndex(asciiLetter);
}


uint8_t awFmAsciiNucleotideToLetterIndex(const uint8_t asciiLetter){
  static const uint8_t charBitmask          = 0x1F;
  static const uint8_t lookupTableD[3]      = {0,2,3};
  static const uint8_t lookupTableH[3]      = {0,1,3};
  static const uint8_t nucleotideLookup[8]  = {0,0,0,1,3,0,0,2};

  uint8_t maskedLetter = asciiLetter & charBitmask;
  switch(maskedLetter){
    case 0x12:  //R or r
      //return A or G, represented by 0 or 2
      return (rand() & 1) << 1;
      break;
    case 0x19: //Y or y
      //return C or T, represented by 1 or 3
      return ((rand() & 1) << 1) | 1;
      break;
    case 0x13: //S or s
      //return C or G, represented by 1 or 2
      return 1 << (rand() & 1);
      break;
    case 0x17: //W or w
      //return A or T, represented by 0 or 3
     return (rand() & 1)? 0: 3;
     break;
    case 0x0b: //K or k
      //return G or T, represented by 0 or 3
      return 2 | (rand() & 1);
      break;
    case 0x0D: //M or m
      //return A or C, represented by 0 or 3
      return rand() & 1;
      break;
    case 0x02: //B or b
      //return C, G, or T
      return (rand() % 3) + 1;
      break;
    case 0x16: //V or v
      //return A, C, or G
      return (rand() % 3);
      break;
    case 0x04: //D or d
      //return A, C, or G
      return lookupTableD[rand()%3];
      break;
    case 0x08: //H or h
      return lookupTableH[rand()%3];
      break;
    default:
      return nucleotideLookup[maskedLetter & 0x07];
  }
}


uint8_t awFmAsciiAminoAcidToLetterIndex(const uint8_t asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  static const uint8_t charBitmask = 0x1F;

  //frequency indexed letter lookup table. values of 20 are illegal characters, and should not match to anything.
  static const uint8_t letterEncodings[26] = {
    //letters A - Z
    20,0,20, 1, 2, 3, 4, 5,
    6, 7, 20, 8, 9,10,11,20,
    12,13,14,15,16,20,17,18,
    20,19};

  uint8_t maskedLetter = asciiLetter & charBitmask;

  if(__builtin_expect(maskedLetter == ('B' & charBitmask), 0)){
    maskedLetter = (rand() & 1)? 'D'& charBitmask: 'N'& charBitmask;
  }
  else if(__builtin_expect(maskedLetter == ('Z' & charBitmask), 0)){
    maskedLetter = (rand() & 1)? 'E'& charBitmask: 'Q'& charBitmask;
  }

  return letterEncodings[maskedLetter];
}


uint8_t awFmNucleotideLetterIndexToAscii(const uint8_t letterIndex){
  static const uint8_t letterLookup[4] = {'A','C','G','T'};
  return letterLookup[letterIndex];
}


uint8_t awFmAminoAcidLetterIndexToAscii(const uint8_t letterIndex){
  static const uint8_t letterLookup[20] = {
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
    'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
    'T', 'V', 'W', 'Y'};
    return letterLookup[letterIndex];
}


uint8_t awFmAminoAcidAsciiLetterToCompressedVectorFormat(const uint8_t asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  static const char charBitmask = 0x1F;
  static const uint8_t compressedVectorLetters[32] = {
    0x00, 0x0C, 0x00, 0x17, 0x03, 0x06, 0x1E, 0x1A,
    0x1B, 0x19, 0x00, 0x15, 0x1C, 0x1D, 0x08, 0x00,
    0x09, 0x04, 0x13, 0x0A, 0x05, 0x00, 0x16, 0x01,
    0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  };

  return compressedVectorLetters[asciiLetter & charBitmask];
}


uint8_t awFmAminoAcidCompressedVectorToAscii(const uint8_t compressedVectorLetter){
  const uint8_t letterLookup[32] = {
    '$', 'W', 'Y', 'D', 'Q', 'T', 'E', '?',
    'N', 'P', 'S', '?', 'A', '?', '?', '?',
    '?', '?', '?', 'R', '?', 'K', 'U', 'C',
    '?', 'I', 'G', 'H', 'L', 'M', 'F', '?'};
  return letterLookup[compressedVectorLetter];
}


uint8_t awFmAminoAcidCompressedVectorToLetterIndex(const uint8_t compressedVectorLetter){
  static const uint8_t letterLookup[32] = {
    20, 18, 19, 2, 13, 16, 3,  20,
    11, 12, 15, 20, 0, 20, 20, 20,
    20, 20, 20, 14, 20, 8, 17, 1,
    20, 7,  5,  6,  9, 10, 4,  20};

  return letterLookup[compressedVectorLetter];
}
