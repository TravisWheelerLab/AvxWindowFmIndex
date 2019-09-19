#include "AwFmLetter.h"
#include <stdlib.h>




uint8_t awFmAsciiLetterToLetterIndex(char asciiLetter){
  //bitmask to alias lowercase letter to uppercase letters
  const char charBitmask = 0x1F;
  asciiLetter = asciiLetter & charBitmask;

  //frequency indexed letter lookup table. values of 21 are illegal characters, and should not match to anything.
  //by frequency, L,A,G,V,E,S,I,K,R,D,T,P,N,Q,F,Y,M,H,C,W
  const letter_t letterOffsets[32] = {
    //letters A - Z
    21, 1,21,18, 9, 4,14, 2,
    17, 6,21, 7, 0,16,12,21,
    11,13, 8, 5,10,21, 3,19,
    21,15,21,21,21,21,21,21
  };

  if(asciiLetter)
  if(asciiLetter == ('B' & charBitmask)){
    asciiLetter = (rand() & 1)? 'D'& charBitmask: 'N'& charBitmask;
  }
  else if(asciiLetter == ('Z' & charBitmask)){
    asciiLetter = (rand() & 1)? 'E'& charBitmask: 'Q'& charBitmask;
  }
  else if(asciiLetter == 'X' & charBitmask){
    //X is the ambiguity character for completely unknown, randomize the letter.
    asciiLetter = (rand() % 19) + 1;
    //if we randomized to another illegal character, add 1 to get a legal char.
    if(letterOffsets[asciiLetter] == 21){
      asciiLetter += 1;
    }
  }

  return letterOffsets[asciiLetter];
}
