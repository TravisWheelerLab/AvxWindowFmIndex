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


void populateKmerSeedTable(struct AwFmIndex *restrict const index);

void populateKmerSeedTableRecursive(struct AwFmIndex *restrict const index, struct AwFmSearchRange range,
  uint8_t currentKmerLength, uint64_t currentKmerIndex, uint64_t letterIndexMultiplier);


void createSequenceEndKmerEncodings(struct AwFmIndex *restrict const index,
  const uint8_t *restrict const sequence, const uint64_t sequenceLength);


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

  populateKmerSeedTable(indexData);

  createSequenceEndKmerEncodings(indexData, sequence, sequenceLength);


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
        uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
        uint8_t letterIndex = awFmAsciiNucleotideToLetterIndex(sequence[positionInBwt]);
        baseOccurrences[letterIndex]++;

        letterBitVectorBytes[byteInVector]      = letterBitVectorBytes[byteInVector] | (((letterIndex >> 0) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 32] = letterBitVectorBytes[byteInVector+ 32] | (((letterIndex >> 1) & 0x1) << bitInVectorByte);
      }
      else{
        index->sentinelCharacterPosition = suffixArrayPosition;
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
    index->prefixSums[AW_FM_NUCLEOTIDE_CARDINALITY] = index->bwtLength;
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
        uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
        uint8_t letterIndex = awFmAsciiAminoAcidToLetterIndex(sequence[positionInBwt]);
        uint8_t letterAsVectorFormat = awFmAminoAcidAsciiLetterToCompressedVectorFormat(sequence[positionInBwt]);
        baseOccurrences[letterIndex]++;
        letterBitVectorBytes[byteInVector]        = letterBitVectorBytes[byteInVector]        | (((letterAsVectorFormat >> 0) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 32]   = letterBitVectorBytes[byteInVector + 32]   | (((letterAsVectorFormat >> 1) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 64]   = letterBitVectorBytes[byteInVector + 64]   | (((letterAsVectorFormat >> 2) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 96]   = letterBitVectorBytes[byteInVector + 96]   | (((letterAsVectorFormat >> 3) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 128]  = letterBitVectorBytes[byteInVector + 128]  | (((letterAsVectorFormat >> 4) & 0x1) << bitInVectorByte);
      }
    }

    //set the prefix sums
    uint64_t tempBaseOccurrence;
    uint64_t prevBaseOccurrence = 1;  //1 is for the sentinel character
    for(uint8_t i = 0; i < awFmGetPrefixSumsLength(index->metadata.alphabetType); i++){
      tempBaseOccurrence = baseOccurrences[i];
      baseOccurrences[i] = prevBaseOccurrence;
      prevBaseOccurrence += tempBaseOccurrence;
    }
    memcpy(index->prefixSums, baseOccurrences, AW_FM_AMINO_CARDINALITY * sizeof(uint64_t));
    index->prefixSums[AW_FM_AMINO_CARDINALITY] = index->bwtLength;
  }
}

void populateKmerSeedTable(struct AwFmIndex *restrict const index){
  const uint8_t alphabetCardinality = awFmGetAlphabetCardinality(index->metadata.alphabetType);

  for(uint8_t i = 0; i < alphabetCardinality; i++){
    struct AwFmSearchRange range = {
      .startPtr= index->prefixSums[i],
      .endPtr= (i == (alphabetCardinality-1)? index->bwtLength: index->prefixSums[i+1])-1
    };
    populateKmerSeedTableRecursive(index, range, 1, i, alphabetCardinality);
  }
}


void populateKmerSeedTableRecursive(struct AwFmIndex *restrict const index, struct AwFmSearchRange range,
  uint8_t currentKmerLength, uint64_t currentKmerIndex, uint64_t letterIndexMultiplier){
  const uint8_t alphabetSize = awFmGetAlphabetCardinality(index->metadata.alphabetType);

  const uint8_t kmerLength  = index->metadata.kmerLengthInSeedTable;

  //base case
  if(kmerLength == currentKmerLength){
    index->kmerSeedTable.table[currentKmerIndex] = range;
    return;
  }

  //recursive case
  for(uint8_t extendedLetter = 0; extendedLetter < alphabetSize; extendedLetter++){
    struct AwFmSearchRange newRange = range;
    if(__builtin_expect(awFmSearchRangeIsValid(&range), 1)){
      if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
        awFmNucleotideIterativeStepBackwardSearch(index, &newRange, extendedLetter);
      }
      else{
        awFmAminoIterativeStepBackwardSearch(index, &newRange, extendedLetter);
      }
    }

    uint64_t newKmerIndex = currentKmerIndex + (extendedLetter * letterIndexMultiplier);
    populateKmerSeedTableRecursive(index, newRange, currentKmerLength + 1, newKmerIndex, letterIndexMultiplier * alphabetSize);
  }
}


void createSequenceEndKmerEncodings(struct AwFmIndex *restrict const index,
  const uint8_t *restrict const sequence, const uint64_t sequenceLength){

  const uint64_t alphabetCardinality = awFmGetAlphabetCardinality(index->metadata.alphabetType);
  const uint8_t sequenceEndingLength = index->metadata.kmerLengthInSeedTable - 1;

  for(uint8_t kmerPosition = 0; kmerPosition < sequenceEndingLength; kmerPosition++){
    uint64_t kmerIndexEncoding = 0;
    for(uint8_t letterInKmer = kmerPosition; letterInKmer < sequenceEndingLength; letterInKmer++){
      uint64_t sequencePosition = (sequenceLength - sequenceEndingLength) + letterInKmer;
      kmerIndexEncoding  = (kmerIndexEncoding * alphabetCardinality) +
      (index->metadata.alphabetType == AwFmAlphabetNucleotide?
        awFmAsciiAminoAcidToLetterIndex(sequence[sequencePosition]):
        awFmAsciiNucleotideToLetterIndex(sequence[sequencePosition]));
    }

    //add the implicit 'a' character to the end to extend the kmer to the required length
    //in order to get the correct index into the kmerTable
    for(uint8_t appendedLetterIndex = 0; appendedLetterIndex < kmerPosition; appendedLetterIndex++){
      kmerIndexEncoding *= alphabetCardinality;
    }

    index->kmerSeedTable.sequenceEndingKmerEncodings[kmerPosition] = kmerIndexEncoding;
  }
}
