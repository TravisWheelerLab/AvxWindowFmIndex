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
          awFmAsciiNucleotideLetterSanitize(c)) == 0);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 0);
        break;
      case 'c':
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 1);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 1);
        break;
      case 'g':
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 2);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 2);
        break;
      case 't':
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 3);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 3);
        break;
      case '$':
      testAssert(awFmAsciiNucleotideToLetterIndex(
        awFmAsciiNucleotideLetterSanitize(c)) == 5);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 5);
          break;
      default:
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(c)) == 4);
        testAssert(awFmAsciiNucleotideToLetterIndex(
          awFmAsciiNucleotideLetterSanitize(toupper(c))) == 4);
    }
  }
}


void testAminoIndex(){

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('a')) == 0);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('A'))) == 0);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('b')) == 19);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('B'))) == 19);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('c')) == 1);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('C'))) == 1);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('d')) == 2);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('D'))) == 2);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('e')) == 3);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('E'))) == 3);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('f')) == 4);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('F'))) == 4);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('g')) == 5);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('G'))) == 5);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('h')) == 6);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('H'))) == 6);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('i')) == 7);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('I'))) == 7);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('j')) == 19);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('J'))) == 19);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('k')) == 8);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('K'))) == 8);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('l')) == 9);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('L'))) == 9);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('m')) == 10);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('M'))) == 10);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('n')) == 11);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('N'))) == 11);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('o')) == 19);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('O'))) == 19);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('p')) == 12);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('P'))) == 12);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('q')) == 13);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('Q'))) == 13);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('r')) == 14);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('R'))) == 14);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('s')) == 15);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('S'))) == 15);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('t')) == 16);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('T'))) == 16);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('u')) == 19);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('u'))) == 19);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('v')) == 17);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('v'))) == 17);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('w')) == 18);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('W'))) == 18);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('x')) == 19);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('X'))) == 19);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('y')) == 20);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('Y'))) == 20);

testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex('z')) == 19);
testAssert(awFmAsciiAminoLetterSanitize(
  awFmAsciiAminoAcidToLetterIndex(toupper('Z'))) == 19);


  testAssert(awFmAsciiAminoLetterSanitize(
    awFmAsciiAminoAcidToLetterIndex('$')) == 21);

}

void testCompressedAminos(){
  for(char c = 'a'; c <= 'z'; c++){
    switch(c){
      case 'a':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 0);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 0);
        break;
      case 'c':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 1);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 1);
        break;
      case 'd':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 2);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 2);
        break;
      case 'e':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 3);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 3);
        break;
      case 'f':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 4);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 4);
        break;
      case 'g':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 5);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 5);
        break;
      case 'h':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 6);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 6);
        break;
      case 'i':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 7);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 7);
        break;
      case 'k':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 8);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 8);
        break;
      case 'l':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 9);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 9);
        break;
      case 'm':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 10);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 10);
        break;
      case 'n':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 11);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 11);
        break;
      case 'p':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 12);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 12);
        break;
      case 'q':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 13);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 13);
        break;
      case 'r':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 14);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 14);
        break;
      case 's':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 15);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 15);
        break;
      case 't':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 16);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 16);
        break;
      case 'v':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 17);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 17);
        break;
      case 'w':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 18);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 18);
        break;
      case 'y':
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 20);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 20);
        break;
      default:
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(c)) == 19);
        testAssert(awFmAsciiAminoAcidToLetterIndex(
          awFmAsciiAminoLetterSanitize(toupper(c))) == 19);
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
