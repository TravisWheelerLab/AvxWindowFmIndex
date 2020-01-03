//_XOPEN_SOURCE needs to be defined before including unistd to make it include pread.
//I have absolutely no idea why this is needed, but the linter seems to think so...
#define _XOPEN_SOURCE 500

#include "AwFmIndex.h"
#include "AwFmFile.h"
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include <unistd.h>


static const uint8_t IndexFileFormatIdHeaderLength  = 10;
static const char    IndexFileFormatIdHeader[10]    = "AwFmIndex\n";


enum AwFmReturnCode awFmWriteIndexToFile(struct AwFmIndex *restrict const index,
  const uint64_t *restrict const suffixArray, const uint8_t *restrict const sequence,
  const uint64_t sequenceLength, const char *restrict const fileSrc, const bool allowOverwrite){

  if(__builtin_expect(fileSrc == NULL, 0)){
    return AwFmNoFileSrcGiven;
  }

  if(__builtin_expect(index == NULL, 0)){
    return AwFmNullPtrError;
  }

  if(__builtin_expect(suffixArray == NULL, 0)){
    return AwFmNullPtrError;
  }

  if(__builtin_expect(sequence == NULL, 0)){
    return AwFmNullPtrError;
  }

  char fileOpenMode[4] = {'w','b', (allowOverwrite? 0:'x'), 0};
  index->fileHandle = fopen(fileSrc, fileOpenMode);

  if(index->fileHandle == NULL){
    return AwFmFileAlreadyExists;
  }

  const size_t bytesPerBwtBlock = index->metadata.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);
  const size_t bytesPerBwt = awFmNumBlocksFromBwtLength(index->bwtLength) * bytesPerBwtBlock;
  const uint64_t prefixSumsByteSize = awFmGetAlphabetCardinality(index->metadata.alphabetType) * sizeof(uint64_t);
  const uint64_t kmerTableByteSize  = awFmGetKmerTableLength(index) * sizeof(struct AwFmBackwardRange);

  //write the metadata, bwt length, sentinel locations, and backward bwt.
  struct FileWriteData{const void *dataPtr; const size_t numBytes;};
  struct FileWriteData fileWriteDataArray[10]  = {
    {IndexFileFormatIdHeader,                   IndexFileFormatIdHeaderLength},
    {&index->metadata,                          sizeof(struct AwFmIndexMetadata)},
    {&index->bwtLength,                         sizeof(uint64_t)},
    {&index->backwardSentinelCharacterPosition, sizeof(uint64_t)},
    {&index->forwardSentinelCharacterPosition,  sizeof(uint64_t)},
    {index->backwardBwtBlockList.asNucleotide,  bytesPerBwt},
    {index->forwardBwtBlockList.asNucleotide,   bytesPerBwt},
    {index->prefixSums,                         prefixSumsByteSize},
    {index->kmerSeedTable,                      kmerTableByteSize},
    {sequence,                                  sequenceLength}
  };

  for(uint8_t i = 0; i < 6; i++){
    size_t elementsWritten = fwrite(fileWriteDataArray[i].dataPtr, fileWriteDataArray[i].numBytes, 1, index->fileHandle);
    if(elementsWritten != 1){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
  }

  //only write the forwardBwt (element index 6) if the bwt is bidirectional
  const uint8_t secondPassFileWriteStart = index->metadata.bwtType == AwFmBwtTypeBiDirectional? 6: 7;

  for(uint8_t i = secondPassFileWriteStart; i < 10; i++){
    size_t elementsWritten = fwrite(fileWriteDataArray[i].dataPtr, fileWriteDataArray[i].numBytes, 1, index->fileHandle);
    if(elementsWritten != 1){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
  }

  //write the compressed suffix array.
  for(size_t i = 0; i < index->bwtLength; i+= index->metadata.suffixArrayCompressionRatio){
    size_t elementsWritten = fwrite(suffixArray + i, sizeof(uint64_t), 1, index->fileHandle);
    if(elementsWritten != 1){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
  }

  index->fileDescriptor        = fileno(index->fileHandle);

  fflush(index->fileHandle);
  return AwFmFileWriteOkay;
}


enum AwFmReturnCode awFmReadIndexFromFile(struct AwFmIndex *restrict *restrict index, const char *fileSrc){
  FILE *fileHandle = fopen(fileSrc, "r");
  if(!fileHandle){
    return AwFmFileOpenFail;
  }

  //read the header
  char headerBuffer[IndexFileFormatIdHeaderLength + 1];
  size_t elementsRead = fread(headerBuffer, sizeof(char), IndexFileFormatIdHeaderLength, fileHandle);
  if(elementsRead != IndexFileFormatIdHeaderLength){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  if(strncmp(IndexFileFormatIdHeader, headerBuffer, IndexFileFormatIdHeaderLength) != 0){
    fclose(fileHandle);
    return AwFmFileFormatError;
  }

  //read the metadata
  struct AwFmIndexMetadata metadata;
  elementsRead = fread(&metadata, sizeof(struct AwFmIndexMetadata), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }

  //read the bwt length
  uint64_t bwtLength;
  elementsRead = fread(&bwtLength, sizeof(uint64_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }

  //allocate the index
  (*index) = awFmIndexAlloc(&metadata, bwtLength -1);
  if(*index == NULL){
    fclose(fileHandle);
    return AwFmAllocationFailure;
  }

  elementsRead = fread(&(*index)->backwardSentinelCharacterPosition, sizeof(uint64_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    awFmDeallocIndex(*index);
    return AwFmFileReadFail;
  }

  elementsRead = fread(&(*index)->forwardSentinelCharacterPosition, sizeof(uint64_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    awFmDeallocIndex(*index);
    return AwFmFileReadFail;
  }

  //read the backward bwt block list
  const size_t numBwtBlocks = (*index)->bwtLength / AW_FM_POSITIONS_PER_FM_BLOCK;
  const size_t blockSizeInBytes = (*index)->metadata.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock):
    sizeof(struct AwFmAminoBlock);

  elementsRead = fread((*index)->backwardBwtBlockList.asNucleotide, blockSizeInBytes, numBwtBlocks, fileHandle);
  if(elementsRead != numBwtBlocks){
    fclose(fileHandle);
    awFmDeallocIndex(*index);
    return AwFmFileReadFail;
  }

  //read the forward bwt block list
  if((*index)->metadata.bwtType == AwFmBwtTypeBiDirectional){
    elementsRead = fread((*index)->forwardBwtBlockList.asNucleotide, blockSizeInBytes, numBwtBlocks, fileHandle);
    if(elementsRead != numBwtBlocks){
      fclose(fileHandle);
      awFmDeallocIndex(*index);
      return AwFmFileReadFail;
    }
  }

  const size_t prefixSumsLength     = awFmGetAlphabetCardinality((*index)->metadata.alphabetType);
  const size_t kmerSeedTableLength  = awFmGetKmerTableLength(*index);

  elementsRead = fread((*index)->prefixSums, sizeof(uint64_t), prefixSumsLength, fileHandle);
  if(elementsRead != prefixSumsLength){
    fclose(fileHandle);
    awFmDeallocIndex(*index);
    return AwFmFileReadFail;
  }

  elementsRead = fread((*index)->kmerSeedTable, sizeof(struct AwFmBackwardRange), kmerSeedTableLength, fileHandle);
  if(elementsRead != kmerSeedTableLength){
    fclose(fileHandle);
    awFmDeallocIndex(*index);
    return AwFmFileReadFail;
  }

  (*index)->suffixArrayFileOffset = awFmGetSuffixArrayFileOffset(*index);
  (*index)->sequenceFileOffset    = awFmGetSequenceFileOffset(*index);
  (*index)->fileDescriptor        = fileno(fileHandle);

  return AwFmFileReadOkay;
}


enum AwFmReturnCode awFmReadPositionsFromSuffixArray(const struct AwFmIndex *restrict const index,
  uint64_t *restrict const positionArray, const size_t positionArrayLength){

  for(size_t i = 0; i < positionArrayLength; i++){
    uint64_t seekPosition = index->suffixArrayFileOffset + ( positionArray[i] * sizeof(uint64_t));
    int returnValue = fseek(index->fileHandle, seekPosition, SEEK_SET);
    if(returnValue != 0){
      return AwFmFileReadFail;
    }

    size_t elementsRead = fread(&(positionArray[i]), sizeof(uint64_t), 1, index->fileHandle);
    if(elementsRead != 1){
      return AwFmFileReadFail;
    }
  }
  return AwFmFileReadOkay;
}


enum AwFmReturnCode awFmReadSequenceFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankLength, const size_t postFlankLength,
  char *const sequenceBuffer){

    //check to make sure that the prior flank length won't underflow the sequence.
  bool bufferWouldUnderflowSequence = priorFlankLength > sequencePosition;
  const size_t actualPriorFlankLength = (bufferWouldUnderflowSequence?
    sequencePosition: priorFlankLength);

  const size_t seekPosition = index->sequenceFileOffset + sequencePosition - priorFlankLength;
  int seekResult = fseek(index->fileHandle, seekPosition, SEEK_SET);
  if(seekResult != 0){
    return AwFmFileReadFail;
  }

  const size_t sequenceSegmentLength = actualPriorFlankLength + postFlankLength;
  size_t elementsRead = fread(sequenceBuffer, sizeof(char), sequenceSegmentLength, index->fileHandle);
  if(elementsRead != sequenceSegmentLength){
    return AwFmFileReadFail;
  }

  //add the null terminator
  sequenceBuffer[sequenceSegmentLength] = 0;

  return AwFmFileReadOkay;
}


enum AwFmReturnCode awFmSuffixArrayReadPositionParallel(const struct AwFmIndex *restrict const index,
  struct AwFmBacktrace *restrict const backtracePtr){

  uint64_t suffixArrayPosition = backtracePtr->position / index->metadata.suffixArrayCompressionRatio;
  const size_t suffixArrayFileOffset = index->suffixArrayFileOffset + ( backtracePtr->position * sizeof(uint64_t));
  ssize_t numBytesRead = pread(index->fileDescriptor, &suffixArrayPosition, sizeof(uint64_t), suffixArrayFileOffset);

  if(numBytesRead == sizeof(uint64_t)){
    backtracePtr->position += backtracePtr->offset;
    return AwFmFileReadOkay;
  }
  else{
    backtracePtr->position = -1UL;
    return AwFmFileReadFail;
  }
}


size_t awFmGetSequenceFileOffset(const struct AwFmIndex *restrict const index){
  const size_t prefixSumLength = awFmGetAlphabetCardinality(index->metadata.alphabetType);
  const size_t kmerSeedTableLength = awFmGetKmerTableLength(index);
  return IndexFileFormatIdHeaderLength + sizeof(struct AwFmIndexMetadata) +
    (3 * sizeof(uint64_t)) + (prefixSumLength * sizeof(uint64_t)) +
    (kmerSeedTableLength * sizeof(struct AwFmBackwardRange));
}


size_t awFmGetSuffixArrayFileOffset(const struct AwFmIndex *restrict const index){
  return awFmGetSequenceFileOffset(index) + ((index->bwtLength - 1) * sizeof(char));
}
