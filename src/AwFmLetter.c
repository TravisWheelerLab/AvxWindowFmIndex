#include "AwFmLetter.h"
#include "AwFmIndexStruct.h"
#include <stdlib.h>
#include <time.h>


uint8_t awFmAsciiNucleotideToLetterIndex(const uint8_t asciiLetter){
  uint8_t toLowerCase = asciiLetter | 0x20;
  switch(toLowerCase){
    case 'a': return 1;
    case 'c': return 2;
    case 'g': return 3;
    case 't': return 4;
    default:  return 0;
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
    return 0;
  }
  else{
    static const uint8_t letterEncodings[32] = {
      0, 1, 0, 2, 3, 4, 5, 6,
      7, 8, 0, 9, 10,11,12,0,
      13,14,15,16,17, 0,18,19,
      0, 20, 0, 0, 0, 0, 0, 0};
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
    0x10, 0x0C, 0x17, 0x03, 0x06, 0x1E, 0x1A,
    0x1B, 0x19, 0x15, 0x1C, 0x1D, 0x08, 0x09,
    0x04, 0x13, 0x0A, 0x05, 0x16, 0x01, 0x02};

  return compressedVectorLookup[letterIndex];
}


uint8_t awFmAminoAcidCompressedVectorToLetterIndex(const uint8_t compressedVectorLetter){
  static const uint8_t letterLookup[31] = {
    0, 19,20, 3,14,17, 4, 0,
    12,13,16, 0, 1, 0, 0, 0,
    0,  0, 0,15, 0, 9,18, 2,
    0,  8, 6, 7,10,11, 5};

  return letterLookup[compressedVectorLetter];
}
