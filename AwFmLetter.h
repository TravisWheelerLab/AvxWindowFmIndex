#ifndef AW_FM_LETTER_H
#define AW_FM_LETTER_H

#include "AwFmIndex.h"
#include <stdint.h>


/*
 * Function:  awFmAsciiLetterToLetterIndex
 * --------------------
 * Transforms an ascii character into a bit-compressed index representation.
 *  This function acts as a wrapper around calls to the ascii nucleotide to index and ascii amino to index calls.
 *  CAUTION: characters that are not ASCII letters may alias to represent another letter.
 *    Only letters a-z, A-Z should be given.
 *  NOTE: this function resolves ambiguity codes via the rand() function.
 *    As such, srand should be used if true random generation is desired.
 *
 *  Inputs:
 *    asciiLetter: ascii-encoded letter representing an amino acid, nucleotide, or ambiguity code.
 *    alphabet: Alphabet of the characters being represented
 *
 *  Returns:
 *    index value representing the amino acid or nucleotide
 *
 */
uint8_t awFmLetterToLetterIndex(const uint8_t asciiLetter, const enum AwFmAlphabetType alphabet);


/*
 * Function:  awFmAsciiNucleotideToLetterIndex
 * --------------------
 * Transforms a nucleotide's ascii letter into a bit compressed index.
 *  The letter can be either lower case or upper case.
 *
 *  Inputs:
 *    asciiLetter: ascii representation of the nucleotide.
 *
 *  Returns:
 *    Bit compressed index of the nucleotide, between 0 and 3 inclusive.
 */
uint8_t awFmAsciiNucleotideToLetterIndex(const uint8_t asciiLetter);


/*
 * Function:  awFmAsciiAminoAcidToLetterIndex
 * --------------------
 * Transforms an amino acid's ascii letter into a bit compressed index.
 *  The letter can be either lower case or upper case.
 *  Note: Giving this function anything besides amino acid letters or ambiguity characters
 *    will result in undefined behavior.
 *    As such, the sentinel character $ should not be given to this function.
 *
 *    List of valid letters: a,b,c,d,e,f,g,h,i,k,l,m,n,p,q,r,s,t,v,w,x,y,z
 *      j, o, u, and any non-alphabet letter are unsupported.
 *
 *  Inputs:
 *    asciiLetter: ascii representation of the amino acid.
 *
 *  Returns:
 *    Bit compressed index of the amino acid, between 0 and 19 inclusive.
 */
uint8_t awFmAsciiAminoAcidToLetterIndex(const uint8_t asciiLetter);


/*
 * Function:  awFmNucleotideLetterIndexToAscii
 * --------------------
 * Transforms a nucleotide's letter index into it's ascii character.
 *  Inputs:
 *    letterIndex: index of the nucleotide.
 *
 *  Returns:
 *    Ascii representation of the nucleotide letter.
 */
uint8_t awFmNucleotideLetterIndexToAscii(const uint8_t letterIndex);


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
uint8_t awFmAminoAcidLetterIndexToAscii(const uint8_t letterIndex);


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
uint8_t awFmAminoAcidAsciiLetterToCompressedVectorFormat(const uint8_t asciiLetter);


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
uint8_t awFmAminoAcidCompressedVectorToAscii(const uint8_t compressedVectorLetter);


/*
* Function:  awFmAminoAcidCompressedVectorToLetterIndex
* --------------------
* Transforms a compressed vector representation of an amino acid letter into its letter index.
*  Inputs:
*    compressedVectorLetter: format that the letter is stored in the AVX vectors.
*
*  Returns:
*    Index of the amino acid between 0 and 19, or 20  as the sentinel.
*/
uint8_t awFmAminoAcidCompressedVectorToLetterIndex(const uint8_t compressedVectorLetter);


#endif /* end of include guard: AW_FM_LETTER_H */
