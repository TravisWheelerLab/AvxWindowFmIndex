#ifndef AW_FM_LETTER_H
#define AW_FM_LETTER_H

#include <stdint.h>
#include "AwFmIndexStruct.h"

/*
 * Function:  awFmAsciiNucleotideToLetterIndex
 * --------------------
 * Transforms an ascii nucleotide character into a bit-compressed index
 * representation. Any non-nucleotide (a,c,g,t,u, A,C,G,T, or U) character will
 * be converted to index 4, signifying a sentinel character.
 *
 *  This function should only be given nucleotides or sentinel '$' characters.
 * Passing ambiguity codes will result in undefined behavior. To sanitize out
 * any ambiguity codes, use the function awFmAsciiNucleotideLetterSanitize().
 *
 *  Inputs:
 *    asciiLetter: ascii-encoded letter representing a nucleotide, ambiguity
 * code, or sentinel '$'.
 *
 *  Returns:
 *    index value representing the nucleotide or a sentinel character.
 *
 */
uint8_t awFmAsciiNucleotideToLetterIndex(const uint8_t asciiLetter);

/*
 * Function:  awFmAsciiNucleotideLetterSanitize
 * --------------------
 * With a given ascii letter, returns the input letter if it is a valid
 * nucleotide. Otherwise, this function returns a sentinel '$' character.
 *  Inputs:
 *    asciiLetter: ascii-encoded letter representing a nucleotide, ambiguity
 * code, or sentinel.
 *
 *  Returns:
 *    The input nucleotide letter or a sentinel '$' character.
 *
 */
uint8_t awFmAsciiNucleotideLetterSanitize(const uint8_t asciiLetter);

/*
 * Function:  awNucleotideLetterIndexToCompressedVector
 * --------------------
 * Transforms an nucleotide's letter-index form into a compressed vector
 * representation.
 *
 *  Inputs:
 *    letterIndex: letter index of the nucleotide.
 *
 *  Returns:
 *    Compressed vector representation of the nucleotide.
 */
uint8_t awFmNucleotideLetterIndexToCompressedVector(const uint8_t asciiLetter);

/*
 * Function:  awFmNucleotideCompressedVectorToLetterIndex
 * --------------------
 * Transforms a compressed vector representation of a nucleotide letter into its
 * letter index. Inputs: compressedVectorLetter: format that the letter is
 * stored in the AVX vectors.
 *
 *  Returns:
 *    Letter index of the nucleotide.
 */
uint8_t awFmNucleotideCompressedVectorToLetterIndex(
    const uint8_t compressedVectorLetter);

/*
 * Function:  awFmAsciiAminoAcidToLetterIndex
 * --------------------
 * Transforms an ascii amino acid character into a bit-compressed index
 * representation. Any ambiguity character or sentinel '$' character will be
 * converted into a sentinel index of 20.
 *
 *  This function should only be given amino acid or sentinel '$' characters.
 * Passing ambiguity codes will result in undefined behavior. To sanitize out
 * any ambiguity codes, use the function awFmAsciiAminoLetterSanitize(). Inputs:
 *    asciiLetter: ascii-encoded letter representing a nucleotide, ambiguity
 * code, or sentinel '$'.
 *
 *  Returns:
 *    index value representing the amino, or the sentinel index of 20.
 *
 */
uint8_t awFmAsciiAminoAcidToLetterIndex(const uint8_t asciiLetter);

/*
 * Function:  awFmAsciiAminoLetterSanitize
 * --------------------
 * With a given ascii letter, returns the input letter if it is a valid amino
 * acid. Otherwise, this function returns a sentinel '$' character. Inputs:
 *    asciiLetter: ascii-encoded letter representing an amino acid, ambiguity
 * code, or sentinel.
 *
 *  Returns:
 *    The input amino letter or a sentinel '$' character.
 *
 */
uint8_t awFmAsciiAminoLetterSanitize(const uint8_t asciiLetter);

/*
 * Function:  awFmAminoAcidLetterIndexToCompressedVector
 * --------------------
 * Transforms an amino acid's in letter-index form into a compressed vector
 * representation.
 *
 *  Inputs:
 *    letterIndex: letter index of the amino or the sentinel '$' index (value
 * 20).
 *
 *  Returns:
 *    Compressed vector representation of the amino or sentinel.
 */
uint8_t awFmAminoAcidLetterIndexToCompressedVector(const uint8_t letterIndex);

/*
 * Function:  awFmAminoAcidCompressedVectorToLetterIndex
 * --------------------
 * Transforms a compressed vector representation of an amino acid letter into
 * its letter index. Inputs: compressedVectorLetter: format that the letter is
 * stored in the AVX vectors.
 *
 *  Returns:
 *    Index of the amino acid between 0 and 19, or 20  as the sentinel.
 */
uint8_t awFmAminoAcidCompressedVectorToLetterIndex(
    const uint8_t compressedVectorLetter);

/*
 * Function:  awFmLetterIsAmbiguous
 * --------------------
 * Determines if the given character representation is ambiguous.
 *  Inputs:
 *      letter: ascii character to check for ambiguity
 *      alphabet: alphabet the character is from.
 *
 *  Returns:
 *    true if the letter is an ambiguity code, false if it
 *      represents a specific nucleotide or amino acid.
 */
bool awFmLetterIsAmbiguous(const char letter,
                           const enum AwFmAlphabetType alphabet);

#endif /* end of include guard: AW_FM_LETTER_H */
