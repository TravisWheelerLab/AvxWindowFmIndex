#ifndef AW_FM_LETTER_H
#define AW_FM_LETTER_H

#include "AwFmIndex.h"
#include <stdint.h>

uint8_t awFmLetterToLetterIndex(const uint8_t asciiLetter, const enum AwFmAlphabetType alphabet);

uint8_t awFmNucleotideLetterIndexToAscii(const uint8_t letterIndex);

uint8_t awFmAminoAcidAsciiLetterToCompressedVectorFormat(const uint8_t asciiLetter);

uint8_t awFmAminoAcidCompressedVectorToLetterIndex(const uint8_t compressedVectorLetter);

uint8_t awFmAminoAcidCompressedVectorToAscii(const uint8_t compressedVectorLetter);

uint8_t awFmAminoAcidLetterIndexToAscii(const uint8_t letterIndex);

#endif /* end of include guard: AW_FM_LETTER_H */
