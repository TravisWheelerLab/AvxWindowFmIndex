#ifndef AW_FM_LETTER_H
#define AW_FM_LETTER_H

#include <stdint.h>


uint8_t awFmAsciiLetterToLetterIndex(const char asciiLetter);
uint8_t awFmAsciiLetterToCompressedVectorFormat(const char asciiLetter);

void awFmAsciiLettersToBitCompressedAvxVectors(const char *const restrict letters, const uint8_t *letterData);

void

#endif /* end of include guard: AW_FM_LETTER_H */
