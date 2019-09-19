#ifndef AW_FM_LETTER_H
#define AW_FM_LETTER_H

#include <stdint.h>

/*converts ASCII characters to frequency encoded values for query sequence*/
uint8_t awFmAsciiLetterToLetterIndex(const char asciiLetter);

/*converts ASCII characters to stridex, compressed vector format for database sequence*/
uint8_t awFmAsciiLetterToCompressedVectorFormat(const char asciiLetter);

#endif /* end of include guard: AW_FM_LETTER_H */
