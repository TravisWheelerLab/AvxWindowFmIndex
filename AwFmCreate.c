#include "AwFmIndex.h"
#include "AwFmCreate.h"
#include "AwFmLetter.h"
#include "AwFmSearch.h"
#include "AwFmFile.h"
#include <stdlib.h>
#include <string.h>
#include <divsufsort64.h>

void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t sequenceLength,
  const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray);


void populateKmerSeedTable(struct AwFmIndex *restrict const index, struct AwFmSearchRange searchRange,
  uint8_t currentKmerLength, uint64_t currentKmerIndex);


enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex *restrict *index,
  const struct AwFmIndexMetadata *restrict const metadata, const uint8_t *restrict const sequence, const size_t sequenceLength,
  const char *restrict const fileSrc, const bool allowFileOverwrite){

  //first, do a sanity check on inputs
  if(metadata == NULL){
    return AwFmNullPtrError;
  }
  if(metadata->suffixArrayCompressionRatio <= 0){
    return AwFmFileFormatError;
  }

  //allocate the index and all internal arrays.
  struct AwFmIndex *restrict indexData = awFmIndexAlloc(metadata, sequenceLength + 1);
  if(indexData == NULL){
    return AwFmAllocationFailure;
  }

  memcpy(&indexData->metadata, metadata, sizeof(struct AwFmIndexMetadata));


  //set the bwtLength
  indexData->bwtLength = sequenceLength + 1;

  uint64_t *suffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
  //hard code in the sentinel character, since it will always be in the first element, we know
  //the value, and libdivsufsort doesn't add the sentinel to the beginning when it creates
  //the suffix array.
  suffixArray[0] = sequenceLength;

  if(suffixArray == NULL){
    awFmDeallocIndex(indexData);
    return AwFmAllocationFailure;
  }
  //create the suffix array, storing it starting in the second element of the suffix array we allocated.
  //this doesn't clobber the sentinel we added earlier, and makes for easier bwt creation.
  uint64_t divSufSortReturnCode = divsufsort64(sequence, (int64_t*)(suffixArray + 1), sequenceLength);
  if(divSufSortReturnCode < 0){
    free(suffixArray);
    awFmDeallocIndex(indexData);
    return AwFmSuffixArrayCreationFailure;
  }

  //set the bwt and prefix sums
  setBwtAndPrefixSums(indexData, indexData->bwtLength, sequence, suffixArray);

  struct AwFmSearchRange searchRange;
  populateKmerSeedTable(indexData, searchRange, 0, 0);


  //set the index as an out argument.
  *index = indexData;

  indexData->suffixArrayFileOffset = awFmGetSuffixArrayFileOffset(*index);
  indexData->sequenceFileOffset    = awFmGetSequenceFileOffset(*index);
  //file descriptor will be set in awFmWriteIndexToFile

  //create the file
  enum AwFmReturnCode returnCode = awFmWriteIndexToFile(indexData, suffixArray, sequence, sequenceLength,
    fileSrc, allowFileOverwrite);

  //since the suffix array is allocated in this function, we need to clean it up before returning.
  free(suffixArray);

  return returnCode;
}



void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t bwtLength,
  const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray){
  if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
    uint64_t baseOccurrences[4] = {0};

    for(uint64_t suffixArrayPosition = 0; suffixArrayPosition < bwtLength; suffixArrayPosition++){
      const size_t  blockIndex      = suffixArrayPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t positionInBlock = suffixArrayPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t byteInVector    = positionInBlock / 8;
      const uint8_t bitInVectorByte = positionInBlock % 8;
      struct AwFmNucleotideBlock *nucleotideBlockPtr = &index->bwtBlockList.asNucleotide[blockIndex];
      uint8_t *restrict const letterBitVectorBytes = (uint8_t*)nucleotideBlockPtr->letterBitVectors;

      if(__builtin_expect(positionInBlock == 0, 0)){
        //when we start a new block, copy over the base occurrences, and initialize the bit vectors
        memcpy(nucleotideBlockPtr->baseOccurrences, baseOccurrences, AW_FM_NUCLEOTIDE_CARDINALITY * sizeof(uint64_t));
        memset(nucleotideBlockPtr->letterBitVectors, 0, sizeof(__m256i) * AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW);
      }

      uint64_t sequencePositionInSuffixArray = suffixArray[suffixArrayPosition];
      if(__builtin_expect(sequencePositionInSuffixArray != 0, 1)){
        // printf("\tsequence position %zu, ", sequencePositionInSuffixArray);
        // printf("bwt letter: %c.\n", sequence[sequencePositionInSuffixArray - 1]);
        uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
        uint8_t letterIndex = awFmAsciiNucleotideToLetterIndex(sequence[positionInBwt]);
        letterBitVectorBytes[byteInVector]      = letterBitVectorBytes[byteInVector] | (((letterIndex >> 0) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 32] = letterBitVectorBytes[byteInVector] | (((letterIndex >> 1) & 0x1) << bitInVectorByte);
      }
    }

    //set the prefix sums
    uint64_t tempBaseOccurrence;
    uint64_t prevBaseOccurrence = 1;  //1 is for the sentinel character
    for(uint8_t i = 0; i < AW_FM_NUCLEOTIDE_CARDINALITY; i++){
      tempBaseOccurrence = baseOccurrences[i];
      baseOccurrences[i] = prevBaseOccurrence;
      prevBaseOccurrence += tempBaseOccurrence;
    }
    memcpy(index->prefixSums, baseOccurrences, AW_FM_NUCLEOTIDE_CARDINALITY * sizeof(uint64_t));
  }
  else{
    uint64_t baseOccurrences[20] = {0};

    for(uint64_t suffixArrayPosition = 0; suffixArrayPosition < bwtLength; suffixArrayPosition++){
      const size_t  blockIndex      = suffixArrayPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t positionInBlock = suffixArrayPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t byteInVector    = positionInBlock / 8;
      const uint8_t bitInVectorByte = positionInBlock % 8;
      struct AwFmAminoBlock *restrict const aminoBlockPointer = &index->bwtBlockList.asAmino[blockIndex];
      uint8_t *restrict const letterBitVectorBytes = (uint8_t*)aminoBlockPointer->letterBitVectors;

      if(__builtin_expect(positionInBlock == 0, 0)){
        //when we start a new block, copy over the base occurrences, and initialize the bit vectors
        memcpy(aminoBlockPointer->baseOccurrences, baseOccurrences, AW_FM_AMINO_CARDINALITY * sizeof(uint64_t));
        memset(aminoBlockPointer->letterBitVectors, 0, sizeof(__m256i) * AW_FM_AMINO_VECTORS_PER_WINDOW);
      }

      uint64_t sequencePositionInSuffixArray = suffixArray[suffixArrayPosition];
      if(__builtin_expect(sequencePositionInSuffixArray != 0, 1)){
        // printf("\tsequence position %zu, ", sequencePositionInSuffixArray);
        // printf("bwt letter: %c.\n", sequence[sequencePositionInSuffixArray - 1]);
        uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
        uint8_t letterIndex = awFmAsciiAminoAcidToLetterIndex(sequence[positionInBwt]);
        letterBitVectorBytes[byteInVector]        = letterBitVectorBytes[byteInVector] | (((letterIndex >> 0) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 32]   = letterBitVectorBytes[byteInVector] | (((letterIndex >> 1) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 64]   = letterBitVectorBytes[byteInVector] | (((letterIndex >> 2) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 96]   = letterBitVectorBytes[byteInVector] | (((letterIndex >> 3) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 128]  = letterBitVectorBytes[byteInVector] | (((letterIndex >> 4) & 0x1) << bitInVectorByte);
      }
    }

    //set the prefix sums
    uint64_t tempBaseOccurrence;
    uint64_t prevBaseOccurrence = 1;  //1 is for the sentinel character
    for(uint8_t i = 0; i < AW_FM_AMINO_CARDINALITY; i++){
      tempBaseOccurrence = baseOccurrences[i];
      baseOccurrences[i] = prevBaseOccurrence;
      prevBaseOccurrence += tempBaseOccurrence;
    }
    memcpy(index->prefixSums, baseOccurrences, AW_FM_AMINO_CARDINALITY * sizeof(uint64_t));
  }
}



void populateKmerSeedTable(struct AwFmIndex *restrict const index, struct AwFmSearchRange searchRange,
  uint8_t currentKmerLength, uint64_t currentKmerIndex){
  const uint8_t alphabetSize = awFmGetAlphabetCardinality(index->metadata.alphabetType);

  const uint8_t kmerLength  = index->metadata.kmerLengthInSeedTable;
    // printf("populating kmer seed table, card %u, kmer length %u, current kmerLength %u, cur ind %zu\n", alphabetSize, kmerLength, currentKmerLength, currentKmerIndex);

  //base case
  if(kmerLength == currentKmerLength){
    // printf("memoizing kmer range %zu-%zu into index %zu\n", searchRange.startPtr, searchRange.endPtr, currentKmerIndex);
    memcpy(&index->kmerSeedTable[currentKmerIndex], &searchRange, sizeof(struct AwFmSearchRange));
    return;
  }
  else if(currentKmerLength == 0){
    //initial case
    // printf("init case!\n");
    searchRange.startPtr  = 1;
    searchRange.endPtr    = index->bwtLength - 1;
  }

  //recursive case
  // printf("recursive case\n");
  for(uint8_t extendedLetter = 0; extendedLetter < alphabetSize; extendedLetter++){
    //... do update search range
    if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
      //update the searchRange
      // printf("updating range...\n");
      awFmNucleotideIterativeStepBackwardSearch(index, &searchRange, extendedLetter);
      // printf("added ascii letter %c, new range is %zu-%zu\n", extendedLetter, searchRange.startPtr, searchRange.endPtr);
    }
    else{
      awFmAminoIterativeStepBackwardSearch(index, &searchRange, extendedLetter);
    }

    uint64_t newKmerIndex = (currentKmerIndex * alphabetSize)+ extendedLetter;
    populateKmerSeedTable(index, searchRange, currentKmerLength + 1, newKmerIndex);
  }

  // printf("leaving function!\n");
}
