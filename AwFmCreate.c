#include "AwFmIndex.h"
#include "AwFmCreate.h"
#include "AwFmLetter.h"
#include "AwFmSearch.h"
#include "AwFmFile.h"
#include <stdlib.h>
#include <string.h>
#include <divsufsort64.h>

void createBwt(struct AwFmIndex *restrict const index, const enum AwFmSearchDirection direction,
  const size_t suffixArrayLength, const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray);

void createPrefixSums(uint64_t *restrict const prefixSums, const uint8_t *restrict const sequence,
  const uint64_t sequenceLength, const enum AwFmAlphabetType alphabet);

void populateKmerSeedTable(struct AwFmIndex *restrict const index, struct AwFmBackwardRange searchRange,
  uint8_t currentKmerLength, uint64_t currentKmerIndex);


enum AwFmReturnCode awFmCreateIndex(const struct AwFmIndex *restrict *index,
  const struct AwFmIndexMetadata *restrict const metadata, const uint8_t *restrict const sequence, const size_t sequenceLength,
  const char *restrict const fileSrc, const bool allowFileOverwrite){

  //allocate the index and all internal arrays.
  struct AwFmIndex *restrict indexData = awFmIndexAlloc(metadata, sequenceLength);
  if(indexData == NULL){
    return AwFmAllocationFailure;
  }

  createPrefixSums(indexData->prefixSums, sequence, sequenceLength, indexData->metadata.alphabetType);
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

  //create the Bwt corresponding to the created suffixArray
  createBwt(indexData, AwFmSearchDirectionBackward,
    sequenceLength + 1, sequence, backwardSuffixArray);


  uint64_t *forwardSuffixArray = NULL;
  //if the index is bidirectional, create the forward suffix array and bwt
  if(metadata->bwtType == AwFmBwtTypeBiDirectional){

    //create the reversed sequence
    uint8_t *restrict const reverseSequence = malloc(sequenceLength * sizeof(uint8_t));
    if(reverseSequence == NULL){
      free(backwardSuffixArray);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }

    for(size_t i = 0; i < sequenceLength; i++){
      reverseSequence[sequenceLength - i - 1] = sequence[i];
    }

    forwardSuffixArray = malloc((sequenceLength + 1) * sizeof(uint64_t));
    if(forwardSuffixArray == NULL){
      free(reverseSequence);
      free(backwardSuffixArray);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }

    //create the forward suffix array
    uint64_t divSufSortReturnCode = divsufsort64(reverseSequence, (int64_t*)forwardSuffixArray, sequenceLength);
    if(divSufSortReturnCode < 0){
      free(reverseSequence);
      free(backwardSuffixArray);
      awFmDeallocIndex(indexData);
      return AwFmSuffixArrayCreationFailure;
    }

    createBwt(indexData, AwFmSearchDirectionForward,
      sequenceLength + 1, reverseSequence, forwardSuffixArray);
    free(reverseSequence);
    free(forwardSuffixArray);
  }

  struct AwFmBackwardRange searchRange;
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


void createBwt(struct AwFmIndex *restrict const index, const enum AwFmSearchDirection direction,
  const size_t suffixArrayLength, const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray){

  //get the correct block list
  union AwFmBwtBlockList blockList;
  uint64_t *sentinelCharacterPositionPtr;
  if(direction == AwFmSearchDirectionBackward){
    blockList = index->backwardBwtBlockList;
    sentinelCharacterPositionPtr = &(index->backwardSentinelCharacterPosition);
  }
  else{
    blockList = index->forwardBwtBlockList;
    sentinelCharacterPositionPtr = &(index->forwardSentinelCharacterPosition);
  }

  if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
    uint64_t baseOccurrences[5] = {0};
    struct AwFmNucleotideBlock workingBlock = {0};
    for(size_t position = 0; position < suffixArrayLength; position++){
      const size_t blockIndex = position / AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t positionInBlock = position % AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t byteInVector    = positionInBlock / 8;
      const uint8_t bitInVectorByte = positionInBlock % 8;

      const uint8_t asciiLetterFromSequence = sequence[suffixArray[position]];
      //if we encounter the sentinel character, notate it, but don't try to put it inside the letterBitVectors
      if(__builtin_expect(asciiLetterFromSequence == '$', 0)){
        *sentinelCharacterPositionPtr = position;
        continue;
      }
      else{
        const uint8_t encodedLetter = awFmNucleotideLetterIndexToAscii(asciiLetterFromSequence);
        baseOccurrences[encodedLetter]++;

        uint8_t * letterBitVectorsAsBytePtr = (uint8_t*) workingBlock.letterBitVectors;
        *(letterBitVectorsAsBytePtr + byteInVector) |= (encodedLetter & 1) << bitInVectorByte;
        *(letterBitVectorsAsBytePtr + byteInVector + sizeof(__m256i)) |= ((encodedLetter >> 1) & 1) << bitInVectorByte;
      }

      //if this was the last position of the block, memcpy the block over and initialize the next block
      if(positionInBlock == (AW_FM_POSITIONS_PER_FM_BLOCK - 1)){
        memcpy(&blockList.asNucleotide[blockIndex], &workingBlock, sizeof(struct AwFmNucleotideBlock));
        memcpy(workingBlock.baseOccurrences, baseOccurrences, AW_FM_NUCLEOTIDE_CARDINALITY * sizeof(uint64_t));
        memset(workingBlock.letterBitVectors, 0, 2 * sizeof(__m256i));
      }
    }

    //transfer whatever is left over in the final block to the block list
    size_t finalBlockIndex = suffixArrayLength / AW_FM_POSITIONS_PER_FM_BLOCK;
    memcpy(&blockList.asNucleotide[finalBlockIndex], &workingBlock, sizeof(struct AwFmNucleotideBlock));

  }
  else{
      // +1 on the array length is to allow illegal characters to collect in the last element.
    uint64_t baseOccurrences[AW_FM_AMINO_CARDINALITY + 1] = {0};
    struct AwFmAminoBlock workingBlock = {0};
    for(size_t position = 0; position < suffixArrayLength; position++){
      const size_t blockIndex = position / AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t positionInBlock = position % AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t byteInVector    = positionInBlock / 8;
      const uint8_t bitInVectorByte = positionInBlock % 8;

      const uint8_t asciiLetterFromSequence = sequence[suffixArray[position]];
      //check for the null terminator
      if(__builtin_expect(asciiLetterFromSequence == '$', 0)){
        *sentinelCharacterPositionPtr = position;
      }
      else{
        const uint8_t encodedLetter           = awFmAsciiAminoAcidToLetterIndex(asciiLetterFromSequence);
        const uint8_t vectorFormatLetter      = awFmAminoAcidAsciiLetterToCompressedVectorFormat(asciiLetterFromSequence);
        baseOccurrences[encodedLetter]++;

        uint8_t * letterBitVectorsAsBytePtr = (uint8_t*) workingBlock.letterBitVectors;
        for(uint8_t bit = 0; bit < 5; bit++){
          const uint8_t bitInVectorFormatLetter = (vectorFormatLetter >> bit) & 1;
          *(letterBitVectorsAsBytePtr + byteInVector + (bit * sizeof(__m256i))) |= bitInVectorFormatLetter << bitInVectorByte;
        }
      }

      //if this was the last position of the block, memcpy the block over and initialize the next block
      if(positionInBlock == (AW_FM_POSITIONS_PER_FM_BLOCK - 1)){
        memcpy(&blockList.asAmino[blockIndex], &workingBlock, sizeof(struct AwFmAminoBlock));
        memcpy(workingBlock.baseOccurrences, baseOccurrences, AW_FM_AMINO_CARDINALITY * sizeof(uint64_t));
        memset(workingBlock.letterBitVectors, 0, 5 * sizeof(__m256i));
      }
    }

    //transfer whatever is left over in the final block to the block list
    size_t finalBlockIndex = suffixArrayLength / AW_FM_POSITIONS_PER_FM_BLOCK;
    memcpy(&blockList.asAmino[finalBlockIndex], &workingBlock, sizeof(struct AwFmAminoBlock));
  }

}


void createPrefixSums(uint64_t *restrict const prefixSums, const uint8_t *restrict const sequence,
  const uint64_t sequenceLength, const enum AwFmAlphabetType alphabet){
  if(alphabet == AwFmAlphabetNucleotide){
    uint64_t tempPrefixSums[5] = {0};
    //count how many times each letter occurrs.
    for(size_t position = 0; position < sequenceLength; position++){
      const uint8_t asciiLetter = sequence[position];
      const uint8_t letterIndex = awFmAsciiNucleotideToLetterIndex(asciiLetter);
      tempPrefixSums[letterIndex]++;
    }

    //perform a scan over the letter occurrences, creating the prefix sums.
    for(uint8_t i = 0; i < AW_FM_NUCLEOTIDE_CARDINALITY; i++){
      tempPrefixSums[i+1]+= tempPrefixSums[i];
    }
    memcpy(prefixSums, tempPrefixSums, AW_FM_NUCLEOTIDE_CARDINALITY * sizeof(uint64_t));

  }
  else{
    // +1 on the array length is to allow illegal characters to collect in the last element.
    uint64_t tempPrefixSums[AW_FM_AMINO_CARDINALITY + 1] = {0};
    //count how many times each letter occurrs.
    for(size_t position = 0; position < sequenceLength; position++){
      const uint8_t asciiLetter = sequence[position];
      const uint8_t letterIndex = awFmAsciiAminoAcidToLetterIndex(asciiLetter);
      tempPrefixSums[letterIndex]++;
    }

    //perform a scan over the letter occurrences, creating the prefix sums.
    for(uint8_t i = 0; i < AW_FM_AMINO_CARDINALITY; i++){
      tempPrefixSums[i+1]+= tempPrefixSums[i];
    }
    memcpy(prefixSums, tempPrefixSums, AW_FM_AMINO_CARDINALITY * sizeof(uint64_t));
  }
}


void populateKmerSeedTable(struct AwFmIndex *restrict const index, struct AwFmBackwardRange searchRange,
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
    memcpy(&index->kmerSeedTable[currentKmerIndex], &searchRange, sizeof(struct AwFmBackwardRange));
    return;
  }
  //recursive case
  else{
    for(uint8_t extendedLetter = 0; extendedLetter < alphabetSize; extendedLetter++){
      //... do update search range
      uint8_t extendedLetterAsAscii;
      if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
          extendedLetterAsAscii = awFmNucleotideLetterIndexToAscii(extendedLetter);
      }
      else{
        extendedLetterAsAscii = awFmAminoAcidLetterIndexToAscii(extendedLetter);
      }

      //update the searchRange
      awFmIterativeStepBackwardSearch(index, &searchRange, extendedLetterAsAscii);
      uint64_t newKmerIndex = (currentKmerIndex * alphabetSize)+ extendedLetter;
      populateKmerSeedTable(index, searchRange, currentKmerLength + 1, newKmerIndex);
    }
  }
}
