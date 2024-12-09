#include "../../src/AwFmLetter.h"

#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../test.h"

char buffer[2048];

void testNucleotideAscii() {
  for (char c = 'a'; c <= 'z'; c++) {
    switch (c) {
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
    case 'u':
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

void testAminoIndex() {
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('a')) == 0);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('A')) == 0);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('b')) == 20);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('B')) == 20);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('c')) == 1);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('C')) == 1);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('d')) == 2);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('D')) == 2);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('e')) == 3);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('E')) == 3);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('f')) == 4);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('F')) == 4);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('g')) == 5);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('G')) == 5);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('h')) == 6);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('H')) == 6);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('i')) == 7);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('I')) == 7);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('j')) == 20);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('J')) == 20);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('k')) == 8);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('K')) == 8);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('l')) == 9);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('L')) == 9);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('m')) == 10);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('M')) == 10);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('n')) == 11);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('N')) == 11);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('o')) == 20);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('O')) == 20);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('p')) == 12);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('P')) == 12);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('q')) == 13);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('Q')) == 13);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('r')) == 14);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('R')) == 14);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('s')) == 15);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('S')) == 15);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('t')) == 16);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('T')) == 16);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('u')) == 20);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('u')) == 20);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('v')) == 17);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('v')) == 17);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('w')) == 18);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('W')) == 18);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('x')) == 20);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('X')) == 20);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('y')) == 19);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('Y')) == 19);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('z')) == 20);
  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('Z')) == 20);

  testAssert(
      awFmAsciiAminoAcidToLetterIndex(awFmAsciiAminoLetterSanitize('$')) == 21);
}

void testCompressedAminos() {
  for (char c = 'a'; c <= 'z'; c++) {
    switch (c) {
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
                     awFmAsciiAminoLetterSanitize(c)) == 19);
      testAssert(awFmAsciiAminoAcidToLetterIndex(
                     awFmAsciiAminoLetterSanitize(toupper(c))) == 19);
      break;
    default:
      testAssert(awFmAsciiAminoAcidToLetterIndex(
                     awFmAsciiAminoLetterSanitize(c)) == 20);
      testAssert(awFmAsciiAminoAcidToLetterIndex(
                     awFmAsciiAminoLetterSanitize(toupper(c))) == 20);
    }
  }
}

int main(int argc, char **argv) {
  srand(time(NULL));
  testNucleotideAscii();
  testAminoIndex();
  testCompressedAminos();
  printf("letter tests finished\n");
}
