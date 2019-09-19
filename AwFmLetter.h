#ifndef AW_FM_LETTER_H
#define AW_FM_LETTER_H

#include <stdint.h>

typedef frequencyIndexedLetter_t  uint8_t;
typedef compressedVectorLetter_t  uint8_t;

frequencyIndexedLetter_t awFmAsciiLetterToLetterIndex(const char asciiLetter);
compressedVectorLetter_t awFmAsciiLetterToCompressedVectorFormat(const char asciiLetter);

void awFmAsciiLettersToBitCompressedAvxVectors(const char *const restrict letters, const uint8_t *letterData);

void 

#endif /* end of include guard: AW_FM_LETTER_H */
