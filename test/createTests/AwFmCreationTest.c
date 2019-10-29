#include <divsufsort64.h>
#include "../test.h"
#include "../../AwFmIndex.h"
#include "../../AwFmLetter.h"
#include "../../AwFmOccupancy.h"
#include "../../AwFmGlobals.h"
#include "../../AwFmCreate.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <time.h>

uint8_t *allocAndAssignDbSequence(const size_t dbSequenceLength);
void suffixArrayCreationTest(const size_t dbSequenceLength);
void blockInitTest(const size_t dbSequenceLength);
void testBlockBaseOccupancies(const struct AwFmIndex *restrict const index, const size_t blockIndex,
  const size_t *restrict const expectedOccupancies);
void testBlockForCorrectValues(const struct AwFmIndex *restrict const index,
  const size_t blockIndex, const size_t suffixArrayLength);
void setRankPrefixSumTest(const size_t dbSequenceLength);

size_t dbSequenceLength = 0;
size_t numTests = 0;
int testNumber = -1;
bool verbose = false;

const uint8_t asciiAminoLookupTableLength = 22;
const uint8_t asciiAminoLookupTable[] ={
  'a','c','d','e','f',
  'g','h','i','k','l',
  'm','n','p','q','r',
  's','t','v','w','y','b','x'
};

const uint8_t frequencyIndexToAsciiTable[] = {
  'l','a','g','v','e','s','i','k','r','d','t','p','n','q','f','y','m','h','c','w'
};

void parseArgs(int argc, char **argv){
  int option = 0;
  while((option = getopt(argc, argv, "vVl:n:t:")) != -1){
    printf("option: %c ", option);
    switch(option){
      case 'V':
      case 'v':
        verbose = true;
        printf("set to verbose output\n");
      break;
      case 'l':
        sscanf(optarg , "%zu" , &dbSequenceLength);
        printf("sequence length: dbSequenceLength\n");
        break;
      case 'n':
        sscanf(optarg , "%zu" , &numTests);
        printf("num tests: %zu\n", numTests);
        break;
      case 't':
        sscanf(optarg , "%d" , &testNumber);
        printf("test number: %d\n", testNumber);
        break;
      case '?':
      printf("unknown option flag given: %c is not a supported flag.\n", optopt);
      break;
    }
  }
}



int main (int argc, char **argv){
  srand(time(NULL));

  parseArgs(argc, argv);
  if(dbSequenceLength == 0){
    printf("error: db sequence length is required, and must not be 0.\n");
    exit(EXIT_FAILURE);
  }

  setRankPrefixSumTest(dbSequenceLength);
    suffixArrayCreationTest(dbSequenceLength);
    blockInitTest(dbSequenceLength);
}

uint8_t *allocAndAssignDbSequence(const size_t dbSequenceLength){
  uint8_t *dbSequence = malloc(dbSequenceLength * sizeof(uint8_t));
  for(size_t i = 0; i < dbSequenceLength; i++){
    dbSequence[i] = asciiAminoLookupTable[rand() % 20];
  }

  return dbSequence;
}


void suffixArrayCreationTest(const size_t dbSequenceLength){
  if(verbose){
    printf("suffix array creation test...\n");
  }

  char stringBuffer[2048];
  uint8_t *dbSequence = allocAndAssignDbSequence(dbSequenceLength);
  uint64_t *suffixArray;

  if(verbose){
    printf("database sequence generated:\n");
    for(size_t i = 0; i < dbSequenceLength; i++){
      printf("%c", dbSequence[i]);
    }
    printf("\n\n");
  }

  enum AwFmReturnCode returnCode =  awFmCreateFullSuffixArray(dbSequence,
    dbSequenceLength, &suffixArray);

  if(verbose){
    printf("testing suffix array creation return code.\n");
  }
  sprintf(stringBuffer, "return value was %i", returnCode);
  testAssertString(returnCode >= 0, stringBuffer);

  testAssertString(suffixArray != NULL, "The suffix array was null.");

  if(verbose){
    printf("testing null terminator location\n");
  }

  //test the null terminator, this should be in the first index of the suffix array.
  sprintf(stringBuffer,
    "the first position in the suffix array should have pointed past the end of the database, (len %zu), but had value %zu instead.",
  dbSequenceLength, suffixArray[0]);
  testAssertString(suffixArray[0] == dbSequenceLength,stringBuffer);

  if(verbose){
    printf("testing suffix array by looking back into sequence.\n");
  }

  uint8_t currentlySeenLetter = 0;
  //test the suffix array, making sure each letter corresponds to a letter that's not
  // lexicographically less than one we've seen.
  for(size_t i = 1; i < dbSequenceLength + 1; i++){
    size_t valueInSuffixArray = suffixArray[i];
    sprintf(stringBuffer,
      "position %zu in suffix array should be less than db length %zu, but was %zu instead", i, dbSequenceLength, valueInSuffixArray);
    testAssertString(valueInSuffixArray < dbSequenceLength, stringBuffer);

    uint8_t letterAtGivenPosition = dbSequence[valueInSuffixArray];
    sprintf(stringBuffer, "character in db should be between 'a'(%i) and 'z'(%i), but was value %i.", 'a','z',letterAtGivenPosition);
    testAssertString(letterAtGivenPosition >= 'a' && letterAtGivenPosition <= 'y', stringBuffer);

    sprintf(stringBuffer, "letter at position %zu was expected to be at least equal to the last seen letter %c, but was %c instead.",
      i, currentlySeenLetter, letterAtGivenPosition);
    testAssertString(currentlySeenLetter <= letterAtGivenPosition, stringBuffer);

    currentlySeenLetter = letterAtGivenPosition;
  }

  if(verbose){
    printf("printing suffix array (map to letter)\n");
    if(suffixArray[0] == dbSequenceLength){
      printf("$");
    }
    else{
      printf("COLOSSAL ERROR: the first position of the suffix array isn't a null terminator, what the hell happened? (it's %c)\n", dbSequence[0]);
    }

    for( size_t i = 1; i < dbSequenceLength + 1; i++){
      printf("%c", dbSequence[suffixArray[i]]);
    }
    printf("\n\n");
  }
  //deallocate the memory for this test.
  free(dbSequence);
  free(suffixArray);

}


void testBlockBaseOccupancies(const struct AwFmIndex *restrict const index, const size_t blockIndex,
  const size_t *restrict const expectedOccupancies){
  char stringBuffer[2048];
  for(uint8_t i = 0; i < 20; i++){
    const size_t expectedOcc = expectedOccupancies[i];
    const size_t storedOcc = index->blockList[blockIndex].baseOccupancies[i];
    sprintf(stringBuffer, "occupancy for frequency index %i, expected %zu, got %zu.", i, expectedOcc, storedOcc);
    testAssertString(expectedOcc == storedOcc, stringBuffer);
  }
}


void testBlockForCorrectLetters(const struct AwFmIndex *restrict const index,
  const size_t blockIndex, const size_t suffixArrayLength, size_t *baseOccupancies){
  char stringBuffer[2048];

  bool isLastBlock  = blockIndex == (index->numBlocks -1);
  //skip the null terminator on the first block, and stop at the end of the suffix array for the last block
  uint_fast16_t blockEndingPosition = isLastBlock? suffixArrayLength - (blockIndex * POSITIONS_PER_FM_BLOCK): POSITIONS_PER_FM_BLOCK;

  for(uint_fast16_t i = 0; i < blockEndingPosition; i++){
    uint8_t assembledVectorLetter  = 0;
    uint8_t byteInLetterVectors    = i / 8;
    uint8_t bitInLetterVectorByte  = i % 8;

    for(int_fast8_t bitIndex = 4; bitIndex >= 0; bitIndex--){
      uint8_t *thisBitsLetterVectorAsBytePtr = (uint8_t *)(&index->blockList[blockIndex].letterBitVectors[bitIndex]);
      uint8_t thisBit = (thisBitsLetterVectorAsBytePtr[byteInLetterVectors] >> bitInLetterVectorByte) & 1;
      assembledVectorLetter = (assembledVectorLetter << 1) | thisBit;
    }

    uint8_t thisLetterAsFrequencyFormat = awFmCompressedVectorLetterToLetterIndex(assembledVectorLetter);
    size_t positionInSuffixArray = blockIndex * POSITIONS_PER_FM_BLOCK + i;
    size_t suffixArrayValue = index->fullSuffixArray[positionInSuffixArray];
    if(suffixArrayValue == 0){
      testAssertString(assembledVectorLetter == 0, stringBuffer);

    }
    else{
      size_t positionInBwt = suffixArrayValue == 0? suffixArrayLength - 1: suffixArrayValue - 1;
      uint8_t asciiLetterAtBwtPosition = index->databaseSequence[positionInBwt];
      uint8_t databaseLetterAsFrequencyFormat = awFmAsciiLetterToLetterIndex(asciiLetterAtBwtPosition);

      sprintf(stringBuffer,
        "letter as SA position %zu (database position %zu) found a frequency index %i in database, but found a %i in the letter bit vectors. Assmebled letter = 0x%x.",
        positionInSuffixArray, suffixArrayValue, databaseLetterAsFrequencyFormat, thisLetterAsFrequencyFormat, assembledVectorLetter);

      testAssertString(thisLetterAsFrequencyFormat == databaseLetterAsFrequencyFormat, stringBuffer);
    }
  }
}


void blockInitTest(const size_t dbSequenceLength){
  if(verbose){
    printf("starting block init test...\n");
  }
  const size_t suffixArrayLength = dbSequenceLength + 1;
  char stringBuffer[2048];
  struct AwFmIndex *index = awFmAlignedAllocAwFmIndex();
  size_t numBlocks = awFmNumBlocksFromSequenceLength(dbSequenceLength);

  if(index == NULL){
    testAssertString(false, "test failure: index was null");
    exit(1);
  }

  if(verbose){
    printf("creating the blockList.\n");
  }

  index->blockList = awFmAlignedAllocBlockList(numBlocks);
  if(index->blockList == NULL){
    testAssertString(false, "test failure, block list returned null.");
    exit(1);
  }

  uint8_t *dbSequence = allocAndAssignDbSequence(dbSequenceLength);


  size_t *suffixArray;
  enum AwFmReturnCode returnCode = awFmCreateFullSuffixArray(dbSequence, dbSequenceLength, &suffixArray);

  if(verbose){

    printf("db sequence:\n");
    for(size_t i = 0; i < dbSequenceLength; i++){
      printf("%c", dbSequence[i]);
    }
    printf("\n\n");

    printf("suffix array (map to char)\n");
    printf(suffixArray[0] == dbSequenceLength? "$": "Colossal error! first element in SA was not null terminator.\n");
    for(size_t i = 1; i < dbSequenceLength + 1; i++){
      printf("%c", dbSequence[suffixArray[i]]);
    }
    printf("\n\n");

    printf("bwt:\n");
    for(size_t i = 0; i < dbSequenceLength + 1; i++){
      size_t suffixArrayValue = suffixArray[i];
      printf("%c", suffixArrayValue == 0? '$': dbSequence[suffixArray[i] - 1]);
    }
    printf("\n\n");
  }

  index->numBlocks        = numBlocks;
  index->databaseSequence = dbSequence;
  index->fullSuffixArray  = suffixArray;

  if(verbose){
    printf("creating the suffix array\n");
  }

  sprintf(stringBuffer, "suffix array creation returned failure code %i", returnCode);
  testAssertString(returnCode >= 0, stringBuffer);

  size_t totalOccupancies[21];
  memset(totalOccupancies, 0, 21 * sizeof(size_t));

  if(verbose){
    printf("testing the block list...\n");
  }

  for(size_t blockIndex = 0; blockIndex < numBlocks; blockIndex++){
    size_t baseOccupancies[21];
    memcpy(baseOccupancies, totalOccupancies, 21 * sizeof(size_t));
    awFmInitBlock(index, blockIndex, totalOccupancies, suffixArrayLength);

    struct AwFmBlock * blockPtr = &index->blockList[blockIndex];
    if(verbose){
      printf("block contents:\n");
      printf("occupancy: ");
      for( uint8_t i = 0; i < 20; i++){
        printf("%zu, ", blockPtr->baseOccupancies[i]);
      }
      printf("\nvectors:");
      for(uint8_t vec = 0; vec < 5; vec++){
        uint8_t *vecAsBytePtr = (uint8_t *) (&blockPtr->letterBitVectors[vec]);
        for(uint8_t byteInVec = 0; byteInVec < BYTES_PER_AVX2_REGISTER; byteInVec++){
          printf("%x ", vecAsBytePtr[byteInVec]);
        }
        printf("\n");
      }
      printf("\n");
    }

    testBlockBaseOccupancies(index, blockIndex, baseOccupancies);
    testBlockForCorrectLetters(index, blockIndex,suffixArrayLength, totalOccupancies);

  }
}


void setRankPrefixSumTest(const size_t dbSequenceLength){
  char stringBuffer[2048];
  size_t fullOccupanciesExpected[20];
  size_t rankPrefixSumsExpected[21];

  memset(fullOccupanciesExpected, 0, 20 * sizeof(size_t));

  struct AwFmIndex *index;
  uint8_t *dbSequence = malloc(dbSequenceLength * sizeof(uint8_t));


  for(size_t i = 0; i < dbSequenceLength; i++){
    const uint8_t letterAtThisPosition = rand() %20;
    dbSequence[i] = frequencyIndexToAsciiTable[letterAtThisPosition];
    fullOccupanciesExpected[letterAtThisPosition]++;
  }

  if(verbose){
    printf("db sequence: ");
    for(size_t i = 0; i < dbSequenceLength; i++){
      printf("%c",dbSequence[i]);
    }
    printf("\n");
  }

  //generate the expected rankPrefixSums
  //start with a 1, since the null terminationg $ is below 'a'
  size_t lettersBelowThisValue = 1;
  for(uint8_t i = 0; i < 20; i++){
    rankPrefixSumsExpected[i] = lettersBelowThisValue;
    lettersBelowThisValue += fullOccupanciesExpected[i];
  }
  rankPrefixSumsExpected[20] = lettersBelowThisValue;

  if(verbose){
    printf("prefix sums:\t\t\t");
    for(uint8_t i = 0; i < 21; i++){
      printf("(%i, %zu), ", i, rankPrefixSumsExpected[i]);
    }
    printf("\n");
  }
  //run the creation code
  enum AwFmReturnCode creationReturnCode = awFmCreateIndex(&index, dbSequence,
    dbSequenceLength, 1);

  if(verbose){
    printf("suffix array: ");
    for( size_t i = 0; i < dbSequenceLength + 1; i++){
      printf("%zu, ", index->fullSuffixArray[i]);
    }

    printf("\nbwt:\t");
    for(size_t i = 0; i < dbSequenceLength + 1; i++){
      size_t suffixArrayValue = index->fullSuffixArray[i];
      printf("%c",suffixArrayValue == 0? '$': dbSequence[suffixArrayValue-1]);
    }

    printf("\ncomputed rank prefix sums:\t");
    for(uint8_t i = 0; i < 21; i++){
      printf("(%i, %zu), ", i, index->rankPrefixSums[i]);
    }
    printf("\n");
  }

  sprintf(stringBuffer, "awfmi creation returned failure code %i", creationReturnCode);
  testAssertString(creationReturnCode >= 0, stringBuffer);

  for(uint8_t i = 0; i < 21; i++){
    size_t prefixSumResult  = index->rankPrefixSums[i];
    size_t expectedPrefixSum = rankPrefixSumsExpected[i];
    sprintf(stringBuffer, "rank prefix sum value for letter @ frequency %i expected value %zu, but got %zu.", i, expectedPrefixSum, prefixSumResult);
    testAssertString(prefixSumResult == expectedPrefixSum, stringBuffer);
  }
}
