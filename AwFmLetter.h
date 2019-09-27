#ifndef AW_FM_LETTER_H
#define AW_FM_LETTER_H

#include <stdint.h>


/*converts ASCII characters to frequency encoded values for query sequence*/
uint8_t awFmAsciiLetterToLetterIndex(const uint8_t asciiLetter);

/*converts ASCII characters to stridex, compressed vector format for database sequence*/
uint8_t awFmAsciiLetterToCompressedVectorFormat(const uint8_t asciiLetter);

/*converts a compressed vector format letter into a frequency encoded letter.*/
uint8_t awFmCompressedVectorLetterToLetterIndex(const uint8_t compressedVectorLetter);

#endif /* end of include guard: AW_FM_LETTER_H */
