#include "../../AwFmLetter.h"
#include "../../AwFmGlobals.h"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
const uint8_t numEncodingAminos = 20;
uint8_t encodingAminos[20] = {'a','c','d','e','f',
                              'g','h','i','k','l',
                              'm','n','p','q','r',
                              's','t','v','w','y'};

uint8_t frequencyEncodingToAscii[20] = {'l','a','g','v','e',
                                        's','i','k','r','d',
                                        't','p','n','q','f',
                                        'y','m','h','c','w'};


void testAssertAsciiToFrequency(){
  resetAssertionNumber();
  for(uint8_t aminoIndex = 0; aminoIndex < numEncodingAminos; aminoIndex++){
    uint8_t asciiLetter = encodingAminos[aminoIndex];
    uint8_t frequencyEncodedLetter = awFmAsciiLetterToLetterIndex(asciiLetter);
    testAssert(frequencyEncodedLetter < 21);

    uint8_t encodingAminoUpperCase = ('a' - 'A') + asciiLetter;
    uint8_t frequencyEncodedLetterUpperCase =  awFmAsciiLetterToLetterIndex(encodingAminoUpperCase);
    testAssert(frequencyEncodedLetterUpperCase < 21);
  }
}


void testAssertAsciiToFrequencyUnique(){
resetAssertionNumber();
  bool aminoWasEncoded[20] = {false};
  for(uint8_t aminoIndex = 0; aminoIndex < numEncodingAminos; aminoIndex++){
    uint8_t asciiLetter = encodingAminos[aminoIndex];
    uint8_t frequencyEncodedLetter = awFmAsciiLetterToLetterIndex(asciiLetter);
    testAssert(aminoWasEncoded[frequencyEncodedLetter] == false);
    aminoWasEncoded[frequencyEncodedLetter] = true;
  }

  //check that every one was set true
  for(uint8_t aminoIndex = 0; aminoIndex < numEncodingAminos; aminoIndex++){
    testAssert(aminoWasEncoded[aminoIndex] == true);
  }
}


void testAssertAmbiguityCharacterToFrequency(){
resetAssertionNumber();
  uint8_t ambiguityCharBResult = awFmAsciiLetterToLetterIndex('b');
  testAssert(ambiguityCharBResult == awFmAsciiLetterToLetterIndex('d') || awFmAsciiLetterToLetterIndex('n'));
  testAssert(ambiguityCharBResult == awFmAsciiLetterToLetterIndex('D') || awFmAsciiLetterToLetterIndex('N'));
  ambiguityCharBResult = awFmAsciiLetterToLetterIndex('B');
  testAssert(ambiguityCharBResult == awFmAsciiLetterToLetterIndex('d') || awFmAsciiLetterToLetterIndex('n'));
  testAssert(ambiguityCharBResult == awFmAsciiLetterToLetterIndex('D') || awFmAsciiLetterToLetterIndex('N'));


  ambiguityCharBResult = awFmAsciiLetterToLetterIndex('z');
  testAssert(ambiguityCharBResult == awFmAsciiLetterToLetterIndex('e') || awFmAsciiLetterToLetterIndex('q'));
  testAssert(ambiguityCharBResult == awFmAsciiLetterToLetterIndex('E') || awFmAsciiLetterToLetterIndex('Q'));
  ambiguityCharBResult = awFmAsciiLetterToLetterIndex('Z');
  testAssert(ambiguityCharBResult == awFmAsciiLetterToLetterIndex('e') || awFmAsciiLetterToLetterIndex('q'));
  testAssert(ambiguityCharBResult == awFmAsciiLetterToLetterIndex('E') || awFmAsciiLetterToLetterIndex('Q'));

  for(int testNum = 0; testNum < 1000; testNum++){
    testAssert(awFmAsciiLetterToLetterIndex('x') < 20);
    testAssert(awFmAsciiLetterToLetterIndex('X') < 20);
  }
}


void testAssertAsciiLetterToCompressedVectorFormat(){
resetAssertionNumber();
  char assertMessage[256];
  for(uint8_t group1AminoCode = 0; group1AminoCode < 12; group1AminoCode++){
    uint8_t asciiEncodedLetter = frequencyEncodingToAscii[group1AminoCode];
    uint8_t compressedVectorFormat = awFmAsciiLetterToCompressedVectorFormat(asciiEncodedLetter);
    uint8_t lastBit = (compressedVectorFormat >> 4) & 1;
    uint8_t numBitsDifferFromLast = 0;
    for(uint8_t bit = 0; bit < 4; bit++){
      if(((compressedVectorFormat >> bit) & 1) != lastBit){
        numBitsDifferFromLast++;
      }
    }
    sprintf(assertMessage, "group1 amino code %c with compressedFormat %d was supposed to have 2 bits \
differing from sign bit, but had %d instead. the encoding the function gave was %d", asciiEncodedLetter, compressedVectorFormat, numBitsDifferFromLast, compressedVectorFormat);
    testAssertString(numBitsDifferFromLast == 2, assertMessage);
  }

  for(uint8_t group2AminoCode = 12; group2AminoCode < 20; group2AminoCode++){
    uint8_t asciiEncodedLetter = frequencyEncodingToAscii[group2AminoCode];
    uint8_t compressedVectorFormat = awFmAsciiLetterToCompressedVectorFormat(asciiEncodedLetter);
    uint8_t lastBit = (compressedVectorFormat >> 4) & 1;
    uint8_t numBitsDifferFromLast = 0;
    for(uint8_t bit = 0; bit < 4; bit++){
      if(((compressedVectorFormat >> bit) & 1) != lastBit){
        numBitsDifferFromLast++;
      }
    }
    sprintf(assertMessage, "group1 amino code %c with compressed vec format %d was supposed to have 1 bit \
differing from sign bit, but had %d instead. the encoding the function gave was %d", asciiEncodedLetter, compressedVectorFormat, numBitsDifferFromLast, compressedVectorFormat);
    testAssertString(numBitsDifferFromLast == 1, assertMessage);
  }
}


void testAssertCompressedVectorToLetterIndex(){
resetAssertionNumber();
  for(uint8_t frequencyEncodedLetter = 0; frequencyEncodedLetter < 20; frequencyEncodedLetter++){
    uint8_t asciiLetter = frequencyEncodingToAscii[frequencyEncodedLetter];
    uint8_t compressedVectorFormat = awFmAsciiLetterToCompressedVectorFormat(asciiLetter);
    testAssert(frequencyEncodedLetter == awFmCompressedVectorLetterToLetterIndex(compressedVectorFormat));
  }
}


int main(int argc, char **argv){
  srand(time(NULL));
  testAssertAsciiToFrequency();
  testAssertAsciiToFrequencyUnique();
  testAssertAmbiguityCharacterToFrequency();
  testAssertAsciiLetterToCompressedVectorFormat();
  testAssertCompressedVectorToLetterIndex();





}
