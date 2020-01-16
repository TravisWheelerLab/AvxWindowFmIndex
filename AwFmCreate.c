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


enum AwFmReturnCode awFmCreateIndex(const struct AwFmIndex *restrict *index,
  const struct AwFmIndexMetadata *restrict const metadata, const uint8_t *restrict const sequence, const size_t sequenceLength,
  const char *restrict const fileSrc, const bool allowFileOverwrite){

  //allocate the index and all internal arrays.
  struct AwFmIndex *restrict indexData = awFmIndexAlloc(metadata, sequenceLength);
  if(indexData == NULL){
    return AwFmAllocationFailure;
  }

  //set the bwtLength
  indexData->bwtLength = sequenceLength + 1;

  uint64_t *backwardSuffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
  if(backwardSuffixArray == NULL){
    awFmDeallocIndex(indexData);
    return AwFmAllocationFailure;
  }

  uint64_t divSufSortReturnCode = divsufsort64(sequence, (int64_t*)backwardSuffixArray, sequenceLength);
  if(divSufSortReturnCode < 0){
    free(backwardSuffixArray);
    awFmDeallocIndex(indexData);
    return AwFmSuffixArrayCreationFailure;
  }

  //set the bwt and prefix sums
  setBwtAndPrefixSums(indexData, sequenceLength + 1, sequence, suffixArray);

  struct AwFmSearchRange searchRange;
  populateKmerSeedTable(indexData, searchRange, 0, 0);

  //set the index as an out argument.
  *index = indexData;

  indexData->suffixArrayFileOffset = awFmGetSuffixArrayFileOffset(*index);
  indexData->sequenceFileOffset    = awFmGetSequenceFileOffset(*index);
  //file descriptor will be set in awFmWriteIndexToFile

  //create the file and return
  return awFmWriteIndexToFile(indexData, backwardSuffixArray, sequence, sequenceLength,
    fileSrc, allowFileOverwrite);
}



void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t sequenceLength,
  const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray){
  if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
    uint64_t baseOccurrences[4] = {0};

    for(uint64_t suffixArrayPosition = 0; suffixArrayPosition < sequenceLength; suffixArrayPosition++){
      const size_t  blockIndex      = suffixArrayPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t positionInBlock = suffixArrayPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t byteInVector    = positionInBlock / 8;
      const uint8_t bitInVectorByte = positionInBlock % 8;
      struct AwFmNucleotideBlock *nucleotideBlockPtr = &index->bwtBlockList.asNucleotide[blockIndex];
      uint8_t *restrict const letterBitVectorBytes = (uint8_t*)nucleotideBlockPtr->letterBitVectors;

      if(positionInBlock == 0){
        //when we start a new block, copy over the base occurrences, and initialize the bit vectors
        memcpy(nucleotideBlockPtr->baseOccurrences, baseOccurrences, AW_FM_NUCLEOTIDE_CARDINALITY * sizeof(uint64_t));
        memset(nucleotideBlockPtr->letterBitVectors, 0, sizeof(__m256i) * AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW);
      }

      uint64_t sequencePosition;
      if(__builtin_expect(suffixArrayPosition == 0, 0)){
        //the first position in the bwt stores the char before the sentinel, aka the last letter in the sequence.
        sequencePosition = sequenceLength - 1;
      }
      else{
        sequencePosition = suffixArray[suffixArrayPosition] - 1;
        if(__builtin_expect(sequencePosition == 0, 0)){
          //sentinel character at this position, don't shift anything in
          continue;
        }
        uint8_t letterIndex = awFmAsciiNucleotideToLetterIndex(sequence[suffixArrayPosition]);
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
  }
  else{
    uint64_t baseOccurrences[20] = {0};

    for(uint64_t suffixArrayPosition = 0; suffixArrayPosition < sequenceLength; suffixArrayPosition++){
      const size_t  blockIndex      = suffixArrayPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t positionInBlock = suffixArrayPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t byteInVector    = positionInBlock / 8;
      const uint8_t bitInVectorByte = positionInBlock % 8;
      struct AwFmAminoBlock *restrict const aminoBlockPointer = &index->bwtBlockList.asAmino[blockIndex];
      uint8_t *restrict const letterBitVectorBytes = (uint8_t*)aminoBlockPointer->letterBitVectors;

      if(positionInBlock == 0){
        //when we start a new block, copy over the base occurrences, and initialize the bit vectors
        memcpy(aminoBlockPointer->baseOccurrences, baseOccurrences, AW_FM_AMINO_CARDINALITY * sizeof(uint64_t));
        memset(aminoBlockPointer->letterBitVectors, 0, sizeof(__m256i) * AW_FM_AMINO_VECTORS_PER_WINDOW);
      }

      uint64_t sequencePosition;
      if(__builtin_expect(suffixArrayPosition == 0, 0)){
        //the first position in the bwt stores the char before the sentinel, aka the last letter in the sequence.
         sequencePosition = sequenceLength-1;
      }
      else{
        uint64_t sequencePosition = suffixArray[suffixArrayPosition] - 1;
        if(__builtin_expect(sequencePosition == 0, 0)){
          //sentinel character at this position
          continue;
        }
      }
      uint8_t letterIndex = awFmAsciiAminoAcidToLetterIndex(sequence[sequencePosition]);
      letterBitVectorBytes[byteInVector]        = letterBitVectorBytes[byteInVector] | (((letterIndex >> 0) & 0x1) << bitInVectorByte);
      letterBitVectorBytes[byteInVector + 32]   = letterBitVectorBytes[byteInVector] | (((letterIndex >> 1) & 0x1) << bitInVectorByte);
      letterBitVectorBytes[byteInVector + 64]   = letterBitVectorBytes[byteInVector] | (((letterIndex >> 2) & 0x1) << bitInVectorByte);
      letterBitVectorBytes[byteInVector + 96]   = letterBitVectorBytes[byteInVector] | (((letterIndex >> 3) & 0x1) << bitInVectorByte);
      letterBitVectorBytes[byteInVector + 128]  = letterBitVectorBytes[byteInVector] | (((letterIndex >> 4) & 0x1) << bitInVectorByte);
    }

    //set the prefix sums
    uint64_t tempBaseOccurrence;
    uint64_t prevBaseOccurrence = 1;  //1 is for the sentinel character
    for(uint8_t i = 0; i < AW_FM_AMINO_CARDINALITY; i++){
      tempBaseOccurrence = baseOccurrences[i];
      baseOccurrences[i] = prevBaseOccurrence;
      prevBaseOccurrence += tempBaseOccurrence;
    }
  }
}



void populateKmerSeedTable(struct AwFmIndex *restrict const index, struct AwFmSearchRange searchRange,
  uint8_t currentKmerLength, uint64_t currentKmerIndex){
  const uint8_t alphabetSize = awFmGetAlphabetCardinality(index->metadata.alphabetType);

  const uint8_t kmerLength  = index->metadata.kmerLengthInSeedTable;

  //init case
  if(currentKmerLength == 0){
    searchRange.startPtr = 0;
    searchRange.endPtr = index->bwtLength;
  }
  //base case
  else if(kmerLength == (currentKmerLength+1)){
    memcpy(&index->kmerSeedTable[currentKmerIndex], &searchRange, sizeof(struct AwFmSearchRange));
    return;
  }
  //recursive case
  else{
    for(uint8_t extendedLetter = 0; extendedLetter < alphabetSize; extendedLetter++){
      //... do update search range
      uint8_t extendedLetterAsAscii;
      if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
          extendedLetterAsAscii = awFmNucleotideLetterIndexToAscii(extendedLetter);
        //update the searchRange
        awFmNucleotideIterativeStepBackwardSearch(index, &searchRange, extendedLetterAsAscii);
      }
      else{
        extendedLetterAsAscii = awFmAminoAcidLetterIndexToAscii(extendedLetter);
        awFmAminoIterativeStepBackwardSearch(index, &searchRange, extendedLetterAsAscii);
      }

      uint64_t newKmerIndex = (currentKmerIndex * alphabetSize)+ extendedLetter;
      populateKmerSeedTable(index, searchRange, currentKmerLength + 1, newKmerIndex);
    }
  }
}
