#include "AwFmLetter.h"
#include "AwFmIndexStruct.h"
#include <stdlib.h>
#include <time.h>


uint8_t awFmAsciiNucleotideToLetterIndex(const uint8_t asciiLetter){
  uint8_t toLowerCase = asciiLetter | 0x20;
  switch(toLowerCase){
    case 'a': return 0;
    case 'c': return 1;
    case 'g': return 2;
    case 't': return 3;
    default:  return 4;
  }
}

uint8_t awFmAsciiNucleotideLetterSanitize(const uint8_t asciiLetter){
  uint8_t toLowerCase = asciiLetter | 0x20;
  switch(toLowerCase){
    case 'a': return 'a';
    case 'c': return 'c';
    case 'g': return 'g';
    case 't': return 't';
    default:  return '$';
  }
}


uint8_t awFmAsciiAminoAcidToLetterIndex(const uint8_t asciiLetter){
  if(__builtin_expect(asciiLetter < 'A', 0)){
    return 20;
  }
  else{
    static const uint8_t letterEncodings[32] = {
      20,0, 20, 1, 2, 3, 4, 5,
      6, 7, 20, 8, 9,10,11,20,
      12,13,14,15,16,20,17,18,
      20,19,20,20,20,20,20,20};
      //find the index and mod (to prevent array out of bounds for weird ascii inputs)
      const uint8_t lookupIndex = (asciiLetter & 0x1F);
      return letterEncodings[lookupIndex];
  }
}


uint8_t awFmAsciiAminoLetterSanitize(const uint8_t asciiLetter){
  uint8_t letterAsLowerCase = asciiLetter | 0x20;
  bool letterIsAmbiguityChar =
    letterAsLowerCase == 'b' ||
    letterAsLowerCase == 'z' ||
    letterAsLowerCase == 'x';

  if(__builtin_expect(letterIsAmbiguityChar, 0)){
    return '$';
  }
  else{
    return asciiLetter;
  }
}


uint8_t awFmAminoAcidLetterIndexToCompressedVector(const uint8_t letterIndex){
  static const uint8_t compressedVectorLookup[21] = {
    0x0C, 0x17, 0x03, 0x06, 0x1E, 0x1A, 0x1B,
    0x19, 0x15, 0x1C, 0x1D, 0x08, 0x09, 0x04,
    0x13, 0x0A, 0x05, 0x16, 0x01, 0x02, 0x00};

  return compressedVectorLookup[letterIndex];
}


uint8_t awFmAminoAcidCompressedVectorToLetterIndex(const uint8_t compressedVectorLetter){
  static const uint8_t letterLookup[31] = {
    20, 18, 19, 2, 13, 16, 3,  20,
    11, 12, 15, 20, 0, 20, 20, 20,
    20, 20, 20, 14, 20, 8, 17, 1,
    20, 7,  5,  6,  9, 10, 4};

  return letterLookup[compressedVectorLetter];
}
