#include "AwFmLetter.h"
#include "AwFmIndexStruct.h"
#include <stdlib.h>
#include <time.h>


//TODO: force the sentinel to be 3 (0x03) for nucletides to make the occGte easier!
//otherwise, it will only show up in occGte when searching for A.

uint8_t awFmLetterToLetterIndex(const uint8_t asciiLetter, const enum AwFmAlphabetType alphabet){
  return (alphabet == AwFmAlphabetAmino)?
    awFmAsciiAminoAcidToLetterIndex(asciiLetter):
    awFmAsciiNucleotideToLetterIndex(asciiLetter);
}


uint8_t awFmAsciiNucleotideToLetterIndex(const uint8_t asciiLetter){
  static const uint8_t charBitmask          = 0x1F;
  static const uint8_t lookupTableD[3]      = {0,2,3};
  static const uint8_t lookupTableH[3]      = {0,1,3};
  static const uint8_t codeLookupTable[32]  = {14, 0, 4, 1, 5,14, 14, 2, 6, 14, 14, 7, 14, 8, 14, 14, 14,
                                                14, 9, 10, 3, 14, 11, 12, 14, 13, 14};

  uint8_t maskedLetter = asciiLetter & charBitmask;
  switch(codeLookupTable[maskedLetter]){
    case 0: return 0;                       break;//A
    case 1: return 1;                       break;//C
    case 2: return 2;                       break;//G
    case 3: return 3;                       break;//T
    case 4: return (rand() % 3) + 1;        break;//B
    case 5: return lookupTableD[rand()%3];  break;//D
    case 6: return lookupTableH[rand()%3];  break;//H
    case 7: return 2 | (rand() & 1);        break;//K
    case 8: return rand() & 1;              break;//M
    case 9: return (rand() & 1) << 1;       break;//R
    case 10:return 1 << (rand() & 1);       break;//S
    case 11:return (rand() % 3);            break;//V
    case 12:return (rand() & 1)? 0: 3;      break;//W
    case 13:return ((rand() & 1) << 1) | 1; break;//Y
    case 14:return rand()%4;//X or N
    default:
    __builtin_unreachable();
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
  else if(__builtin_expect(maskedLetter == ('X' & charBitmask),0)){
    return rand() % 20;
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
    '?', '?', '?', 'R', '?', 'K', 'V', 'C',
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
