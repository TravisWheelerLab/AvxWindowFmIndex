#include "AwFmIndexStruct.h"
#include "AwFmCreate.h"
#include "AwFmLetter.h"
#include "AwFmSearch.h"
#include "AwFmFile.h"
#include <stdlib.h>
#include <string.h>
#include "divsufsort64.h"

void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t sequenceLength,
  const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray);


void populateKmerSeedTable(struct AwFmIndex *restrict const index);

void populateKmerSeedTableRecursive(struct AwFmIndex *restrict const index, struct AwFmSearchRange range,
  uint8_t currentKmerLength, uint64_t currentKmerIndex, uint64_t letterIndexMultiplier);


void compressSuffixArrayInPlace(uint64_t *const suffixArray, uint64_t suffixArrayLength, const uint8_t compressionRatio);


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

  //set the index out arg initally to NULL, if this function fully completes this will get overwritten
  *index = NULL;

  const size_t suffixArrayLength = sequenceLength + 1;
  //create a sanitized copy of the input sequence
  uint8_t *sanitizedSequenceCopy = malloc(suffixArrayLength);
  if(sanitizedSequenceCopy == NULL){
    return AwFmAllocationFailure;
  }
  if(metadata->alphabetType == AwFmAlphabetNucleotide){
    for(size_t i = 0; i < sequenceLength; i++){
      sanitizedSequenceCopy[i] = awFmAsciiNucleotideLetterSanitize(sequence[i]);
    }
  }
  else{
    for(size_t i = 0; i < sequenceLength; i++){
      sanitizedSequenceCopy[i] = awFmAsciiAminoLetterSanitize(sequence[i]);
    }
  }
  // append the final sentinel character as a terminator.
  sanitizedSequenceCopy[suffixArrayLength - 1] = '$';

  //allocate the index and all internal arrays.
  struct AwFmIndex *restrict indexData = awFmIndexAlloc(metadata, suffixArrayLength);
  if(indexData == NULL){
    return AwFmAllocationFailure;
  }

  memcpy(&indexData->metadata, metadata, sizeof(struct AwFmIndexMetadata));

  //init the in memory suffix array to NULL, to be safe. this will get overwritten on success,
  //if the metadata demands in memory SA. If not, this will be left NULL.
  indexData->inMemorySuffixArray = NULL;

  //set the bwtLength
  indexData->bwtLength = suffixArrayLength;

  uint64_t *suffixArray = malloc((suffixArrayLength) * sizeof(uint64_t));

  if(suffixArray == NULL){
    free(sanitizedSequenceCopy);
    awFmDeallocIndex(indexData);
    return AwFmAllocationFailure;
  }

  //create the suffix array, storing it starting in the second element of the suffix array we allocated.
  //this doesn't clobber the sentinel we added earlier, and makes for easier bwt creation.
  int64_t divSufSortReturnCode = divsufsort64(sanitizedSequenceCopy, (int64_t*)(suffixArray), suffixArrayLength);
  if(divSufSortReturnCode < 0){
    free(sanitizedSequenceCopy);
    free(suffixArray);
    awFmDeallocIndex(indexData);
    return AwFmSuffixArrayCreationFailure;
  }

  //set the bwt and prefix sums
  setBwtAndPrefixSums(indexData, indexData->bwtLength, sanitizedSequenceCopy, suffixArray);
  //after generating the bwt, the sequence copy is no longer needed.
  free(sanitizedSequenceCopy);

  populateKmerSeedTable(indexData);


  indexData->suffixArrayFileOffset = awFmGetSuffixArrayFileOffset(indexData);
  indexData->sequenceFileOffset    = awFmGetSequenceFileOffset(indexData);
  //file descriptor will be set in awFmWriteIndexToFile

  //create the file
  enum AwFmReturnCode returnCode = awFmWriteIndexToFile(indexData, suffixArray, sequence, sequenceLength,
    fileSrc, allowFileOverwrite);
    //if suffix array was requested to be kept in memory, realloc it to it's compressed shape
    if(metadata->keepSuffixArrayInMemory){
      //if the suffix array is uncompressed, we get to skip the compression and realloc
      if(metadata->suffixArrayCompressionRatio != 1){
        compressSuffixArrayInPlace(suffixArray, suffixArrayLength,indexData->metadata.suffixArrayCompressionRatio);
        uint64_t *compressedSuffixArray = realloc(suffixArray, awFmGetCompressedSuffixArrayLength(indexData) * sizeof(uint64_t));

        //check for allocation failure in the realloc
        if(compressedSuffixArray == NULL){
          free(suffixArray);
          awFmDeallocIndex(indexData);
          return AwFmAllocationFailure;
        }
        suffixArray = compressedSuffixArray;
      }

      indexData->inMemorySuffixArray = suffixArray;
    }
    else{
      //if the suffix array isn't supposed to be kept in memory, free it to free up memory.
      indexData->inMemorySuffixArray = NULL;
      free(suffixArray);
    }
    //set the index as an out argument.
    *index = indexData;

  return returnCode;
}



void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t bwtLength,
  const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray){
  if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
    uint64_t baseOccurrences[5] = {0};

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
      uint8_t letterIndex;
      if(__builtin_expect(sequencePositionInSuffixArray != 0, 1)){
        uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
        letterIndex = awFmAsciiNucleotideToLetterIndex(sequence[positionInBwt]);
      }
      else{
        letterIndex = AW_FM_NUCLEOTIDE_CARDINALITY;
      }
        baseOccurrences[letterIndex]++;

        letterBitVectorBytes[byteInVector]      = letterBitVectorBytes[byteInVector] | (((letterIndex >> 0) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 32] = letterBitVectorBytes[byteInVector+ 32] | (((letterIndex >> 1) & 0x1) << bitInVectorByte);
        //for the last bit, flip the bit so all nucleotides have a 1, and sentinels have a zero.
        letterBitVectorBytes[byteInVector + 64] = letterBitVectorBytes[byteInVector+64] | (((letterIndex >>2) ^ 1) & 0x01) << bitInVectorByte;
    }

    //set the prefix sums
    uint64_t tempBaseOccurrence;

    //set the base occurrence to the number of sentinels
    uint64_t prevBaseOccurrence = baseOccurrences[AW_FM_NUCLEOTIDE_CARDINALITY];
    for(uint8_t i = 0; i < AW_FM_NUCLEOTIDE_CARDINALITY; i++){
      tempBaseOccurrence = baseOccurrences[i];
      baseOccurrences[i] = prevBaseOccurrence;
      prevBaseOccurrence += tempBaseOccurrence;
    }
    memcpy(index->prefixSums, baseOccurrences, AW_FM_NUCLEOTIDE_CARDINALITY * sizeof(uint64_t));
    index->prefixSums[AW_FM_NUCLEOTIDE_CARDINALITY] = index->bwtLength;
  }
  else{
    uint64_t baseOccurrences[21] = {0};

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
      uint8_t letterIndex;
      if(__builtin_expect(sequencePositionInSuffixArray != 0, 1)){
        uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
        letterIndex = awFmAsciiAminoAcidToLetterIndex(sequence[positionInBwt]);
      }
      else{
        letterIndex = AW_FM_AMINO_CARDINALITY;
      }
        uint8_t letterAsVectorFormat = awFmAminoAcidLetterIndexToCompressedVector(letterIndex);
        baseOccurrences[letterIndex]++;
        letterBitVectorBytes[byteInVector]        = letterBitVectorBytes[byteInVector]        | (((letterAsVectorFormat >> 0) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 32]   = letterBitVectorBytes[byteInVector + 32]   | (((letterAsVectorFormat >> 1) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 64]   = letterBitVectorBytes[byteInVector + 64]   | (((letterAsVectorFormat >> 2) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 96]   = letterBitVectorBytes[byteInVector + 96]   | (((letterAsVectorFormat >> 3) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 128]  = letterBitVectorBytes[byteInVector + 128]  | (((letterAsVectorFormat >> 4) & 0x1) << bitInVectorByte);

    }

    //set the prefix sums
    uint64_t tempBaseOccurrence;
    //set the base occurrence to the number of sentinels
    uint64_t prevBaseOccurrence = baseOccurrences[AW_FM_AMINO_CARDINALITY];
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
    index->kmerSeedTable[currentKmerIndex] = range;
    return;
  }

  //recursive case
  for(uint8_t extendedLetter = 0; extendedLetter < alphabetSize; extendedLetter++){
    struct AwFmSearchRange newRange = range;
    if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
      awFmNucleotideIterativeStepBackwardSearch(index, &newRange, extendedLetter);
    }
    else{
      awFmAminoIterativeStepBackwardSearch(index, &newRange, extendedLetter);
    }

    uint64_t newKmerIndex = currentKmerIndex + (extendedLetter * letterIndexMultiplier);
    populateKmerSeedTableRecursive(index, newRange, currentKmerLength + 1, newKmerIndex, letterIndexMultiplier * alphabetSize);
  }
}


//copies the compressed suffix array over the full suffix array.
void compressSuffixArrayInPlace(uint64_t *const suffixArray, uint64_t suffixArrayLength, const uint8_t compressionRatio){
  for(size_t i = 1; i * compressionRatio < suffixArrayLength; i++){
    suffixArray[i] = suffixArray[i*compressionRatio];
  }
}
