#include "AwFmLetter.h"
#include "AwFmIndex.h"
#include <stdlib.h>
#include <time.h>

//TODO: write documentation to explain that there is no nucleotide compressed vector format since we just use
// the default index as the vector format.

inline uint8_t aminoAcidToLetterIndex(const uint8_t asciiLetter);
inline uint8_t nucleotideToLetterIndex(const uint8_t asciiLetter);


/*
 * Function:  awFmAsciiLetterToLetterIndex
 * --------------------
 * Transforms an ascii character into a bit-compressed index representation.
 *  CAUTION: characters that are not ASCII letters may alias to represent an amino acid. Only letters a-z, A-Z should be given.
 *  NOTICE: this function resolves ambiguity codes via the rand() function.
 *  In order to properly randomize the resolved amino acids, make sure rand is seeded with srand.
 *
 *  Inputs:
 *    asciiLetter: ascii-encoded letter representing an amino acid, nucleotide, or ambiguity code.
 *    alphabet: Alphabet of the characters being represented
 *
 *  Returns:
 *    index value representing the amino acid or nucleotide
 *
 */
uint8_t awFmLetterToLetterIndex(const uint8_t asciiLetter, const enum AwFmAlphabetType alphabet){
  return (alphabet == AwFmAlphabetAminoAcid)?
    aminoAcidToLetterIndex(asciiLetter):
    nucleotideToLetterIndex(asciiLetter);
}


/*
 * Function:  awFmNucleotideLetterIndexToAscii
 * --------------------
 * Transforms a nucleotide's letter index into it's ascii character.
 *  Inputs:
 *    letterIndex: index of the nucleotide.
 *
 *  Returns:
 *    Ascii representation of the nucleotide letter.
 *
 */
uint8_t awFmNucleotideLetterIndexToAscii(const uint8_t letterIndex){
  const uint8_t letterLookup[4] = {'A','C','G','T'};
  return letterLookup[letterIndex];
}

/*
 * Function:  awFmAminoAcidAsciiLetterToCompressedVectorFormat
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
uint8_t awFmAminoAcidAsciiLetterToCompressedVectorFormat(const uint8_t asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  const char charBitmask = 0x1F;
  const uint8_t compressedVectorLetters[32] = {
    0x00, 0x0C, 0x00, 0x17, 0x03, 0x06, 0x1E, 0x1A,
    0x1B, 0x19, 0x00, 0x15, 0x1C, 0x1D, 0x08, 0x00,
    0x09, 0x04, 0x13, 0x0A, 0x05, 0x00, 0x16, 0x01,
    0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  };

  return compressedVectorLetters[asciiLetter & charBitmask];
}


/*
 * Function:  awFmAminoAcidCompressedVectorToLetterIndex
 * --------------------
 * Transforms a compressed vector representation of an amino acid letter into its letter index.
 *  Inputs:
 *    compressedVectorLetter: format that the letter is stored in the AVX vectors.
 *
 *  Returns:
 *    Index of the amino acid between 0 and 19, or 20  as the sentinel.
 *
 */
uint8_t awFmAminoAcidCompressedVectorToLetterIndex(const uint8_t compressedVectorLetter){
  const uint8_t letterLookup[32] = {
    20, 18, 19, 2, 13, 16, 3,  20,
    11, 12, 15, 20, 0, 20, 20, 20,
    20, 20, 20, 14, 20, 8, 17, 1,
    20, 7,  5,  6,  9, 10, 4,  20};
  return letterLookup[compressedVectorLetter];
}


/*
 * Function:  awFmAminoAcidCompressedVectorToAscii
 * --------------------
 * Transforms a compressed vector representation of an amino acid letter into it's ascii character.
 *  Inputs:
 *    compressedVectorLetter: format that the letter is stored in the AVX vectors.
 *
 *  Returns:
 *    Ascii representation of the letter.
 *
 */
uint8_t awFmAminoAcidCompressedVectorToAscii(const uint8_t compressedVectorLetter){
  const uint8_t letterLookup[32] = {
    '$', 'W', 'Y', 'D', 'Q', 'T', 'E', '?',
    'N', 'P', 'S', '?', 'A', '?', '?', '?',
    '?', '?', '?', 'R', '?', 'K', 'U', 'C',
    '?', 'I', 'G', 'H', 'L', 'M', 'F', '?'};
  return letterLookup[compressedVectorLetter];
}

/*
 * Function:  awFmAminoAcidLetterIndexToAscii
 * --------------------
 * Transforms the letter index of an amino acid letter into it's ascii representation.
 *  Inputs:
 *    letterIndex: Index of the letter.
 *
 *  Returns:
 *    Ascii representation of the letter.
 *
 */
uint8_t awFmAminoAcidLetterIndexToAscii(const uint8_t letterIndex){
  const uint8_t letterLookup[21] = {
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
    'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
    'T', 'V', 'W', 'Y', '$'};
    return letterLookup[letterIndex];
}


//private functions
inline uint8_t nucleotideToLetterIndex(const uint8_t asciiLetter){
  const uint8_t charBitmask = 0x1F;
  const uint8_t lookupTableD[3] = {0,2,3};
  const uint8_t lookupTableH[3] = {0,1,3};
  const uint8_t nucleotideLookup[8] = {0,0,0,1,3,0,0,2};

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
      return nucleotideLookup[maskedLetter & 0x03];
  }
}

inline uint8_t aminoAcidToLetterIndex(const uint8_t asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  const uint8_t charBitmask = 0x1F;
  uint8_t maskedLetter = asciiLetter & charBitmask;

  //frequency indexed letter lookup table. values of 20 are illegal characters, and should not match to anything.
  //by frequency, L,A,G,V,E,S,I,K,R,D,T,P,N,Q,F,Y,M,H,C,W
  const uint8_t letterOffsets[32] = {
    //letters A - Z
    20,0,20, 1, 2, 3, 4, 5,
    6, 7, 20, 8, 9,10,11,20,
    12,13,14,15,16,20,17,18,
    20,19,20,20,20,20,20,20
  };

  if(maskedLetter == ('B' & charBitmask)){
    maskedLetter = (rand() & 1)? 'D'& charBitmask: 'N'& charBitmask;
  }
  else if(maskedLetter == ('Z' & charBitmask)){
    maskedLetter = (rand() & 1)? 'E'& charBitmask: 'Q'& charBitmask;
  }

  uint8_t frequencyEncodedLetter = letterOffsets[maskedLetter];
  if(__builtin_expect(frequencyEncodedLetter == 20, 0)){
    return rand() % 20;
  }
  else{
    return frequencyEncodedLetter;
  }
}
