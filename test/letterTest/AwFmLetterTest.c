#include "../../src/AwFmIndexStruct.h"
#include "../../src/AwFmLetter.h"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <ctype.h>



char buffer[2048];

void testNucleotideAscii(){
  for(char c = 'a'; c <= 'z'; c++){
    switch(c){
      case 'a':
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 1);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 1);
        break;
      case 'c':
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 2);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 2);
        break;
      case 'g':
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 3);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 3);
        break;
      case 't':
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 4);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 4);
        break;
      default:
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 0);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 0);
    }
  }
}


void testAminoIndex(){
  for(char c = 'a'; c <= 'z'; c++){
    switch(c){
      case 'a':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 1);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 1);
        break;
      case 'c':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 2);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 2);
        break;
      case 'd':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 3);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 3);
        break;
      case 'e':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 4);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 4);
        break;
      case 'f':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 5);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 5);
        break;
      case 'g':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 6);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 6);
        break;
      case 'h':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 7);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 7);
        break;
      case 'i':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 8);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 8);
        break;
      case 'k':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 9);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 9);
        break;
      case 'l':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 10);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 10);
        break;
      case 'm':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 11);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 11);
        break;
      case 'n':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 12);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 12);
        break;
      case 'p':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 13);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 13);
        break;
      case 'q':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 14);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 14);
        break;
      case 'r':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 15);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 15);
        break;
      case 's':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 16);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 16);
        break;
      case 't':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 17);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 17);
        break;
      case 'v':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 18);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 18);
        break;
      case 'w':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 19);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 19);
        break;
      case 'y':
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 20);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 20);
        break;
      default:
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(c)) == 0);
        testAssert(awFmAsciiAminoLetterSanitize(
          awFmAsciiAminoAcidToLetterIndex(toupper(c))) == 0);
    }
  }
}

void testCompressedAminos(){
  for(char c = 'a'; c <= 'z'; c++){
    switch(c){
      case 'a':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 1);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 1);
        break;
      case 'c':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 2);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 2);
        break;
      case 'd':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 3);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 3);
        break;
      case 'e':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 4);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 4);
        break;
      case 'f':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 5);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 5);
        break;
      case 'g':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 6);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 6);
        break;
      case 'h':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 7);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 7);
        break;
      case 'i':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 8);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 8);
        break;
      case 'k':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 9);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 9);
        break;
      case 'l':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 10);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 10);
        break;
      case 'm':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 11);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 11);
        break;
      case 'n':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 12);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 12);
        break;
      case 'p':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 13);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 13);
        break;
      case 'q':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 14);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 14);
        break;
      case 'r':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 15);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 15);
        break;
      case 's':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 16);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 16);
        break;
      case 't':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 17);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 17);
        break;
      case 'v':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 18);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 18);
        break;
      case 'w':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 19);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 19);
        break;
      case 'y':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 20);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 20);
        break;
      default:
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 0);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 0);
    }
  }
}



int main(int argc, char **argv){
  srand(time(NULL));
  testNucleotideAscii();
  testAminoIndex();
  testCompressedAminos();
  printf("letter tests finished\n");

}
