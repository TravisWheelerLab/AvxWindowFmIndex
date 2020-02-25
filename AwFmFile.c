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

  //open the file
  char fileOpenMode[5] = {'w','+','b', (allowOverwrite? 0:'x'), 0};
  index->fileHandle = fopen(fileSrc, fileOpenMode);

  if(index->fileHandle == NULL){
    return AwFmFileAlreadyExists;
  }

  //write the header
  size_t elementsWritten = fwrite(IndexFileFormatIdHeader, sizeof(char), IndexFileFormatIdHeaderLength, index->fileHandle);
  if(elementsWritten != IndexFileFormatIdHeaderLength){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the metadata, starting with the version number
  elementsWritten = fwrite(&index->metadata.versionNumber, sizeof(uint16_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  elementsWritten = fwrite(&index->metadata.suffixArrayCompressionRatio, sizeof(uint8_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  elementsWritten = fwrite(&index->metadata.kmerLengthInSeedTable, sizeof(uint8_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  uint8_t alphabetType = index->metadata.alphabetType;
  elementsWritten = fwrite(&alphabetType, sizeof(uint8_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the bwt length
  elementsWritten = fwrite(&index->bwtLength, sizeof(uint64_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the sentinel character position
  elementsWritten = fwrite(&index->sentinelCharacterPosition, sizeof(uint64_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  const size_t numBlockInBwt = awFmNumBlocksFromBwtLength(index->bwtLength);
  const size_t bytesPerBwtBlock = index->metadata.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

  elementsWritten = fwrite(index->bwtBlockList.asNucleotide, bytesPerBwtBlock, numBlockInBwt, index->fileHandle);
  if(elementsWritten != numBlockInBwt){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the prefix sums table
  const size_t prefixSumsLength = awFmGetPrefixSumsLength(index->metadata.alphabetType);
  elementsWritten = fwrite(index->prefixSums, sizeof(uint64_t), prefixSumsLength, index->fileHandle);
  if(elementsWritten != prefixSumsLength){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the kmer seed table
  const size_t numElementsInKmerSeedTable = awFmGetKmerTableLength(index);
  elementsWritten = fwrite(index->kmerSeedTable.table, sizeof(struct AwFmSearchRange), numElementsInKmerSeedTable, index->fileHandle);
  if(elementsWritten != numElementsInKmerSeedTable){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the sequence ending
  const size_t sequenceEndingLength = index->metadata.kmerLengthInSeedTable - 1;
  elementsWritten = fwrite(index->kmerSeedTable.sequenceEndingKmerEncodings, sizeof(uint64_t), sequenceEndingLength, index->fileHandle);
  if(elementsWritten != sequenceEndingLength){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the sequence
  elementsWritten = fwrite(sequence, sizeof(char), sequenceLength, index->fileHandle);
  if(elementsWritten != sequenceLength){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the compressed suffix array
  for(size_t i = 0; i < index->bwtLength; i+= index->metadata.suffixArrayCompressionRatio){
    //write the 'implicit' position of the sentinel character at the end in the first element of the SA
    size_t valueToWrite = suffixArray[i];
    size_t elementsWritten = fwrite(&valueToWrite, sizeof(uint64_t), 1, index->fileHandle);
    if(elementsWritten != 1){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
  }

  fflush(index->fileHandle);
  index->fileDescriptor = fileno(index->fileHandle);

  return AwFmFileWriteOkay;
}


enum AwFmReturnCode awFmReadIndexFromFile(struct AwFmIndex *restrict *restrict index, const char *fileSrc){
  if(__builtin_expect(fileSrc == NULL, 0)){
    return AwFmNoFileSrcGiven;
  }

  FILE *fileHandle = fopen(fileSrc, "r");
  if(!fileHandle){
    return AwFmFileOpenFail;
  }

  //create a local-scope pointer for the index, when we're done we'll set the index out-arg to this pointer.
  struct AwFmIndex *restrict indexData;

  //read the header, and check to make sure it matches
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

  //read in the metadata
  struct AwFmIndexMetadata metadata;
  elementsRead = fread(&metadata.versionNumber, sizeof(uint16_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  elementsRead = fread(&metadata.suffixArrayCompressionRatio, sizeof(uint8_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  elementsRead = fread(&metadata.kmerLengthInSeedTable, sizeof(uint8_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  uint8_t alphabetType;
  elementsRead = fread(&alphabetType, sizeof(uint8_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  metadata.alphabetType = alphabetType;


  //read the bwt length
  uint64_t bwtLength;
  elementsRead = fread(&bwtLength, sizeof(uint64_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }

  //allocate the index
  indexData = awFmIndexAlloc(&metadata, bwtLength);
  if(indexData == NULL){
    fclose(fileHandle);
    return AwFmAllocationFailure;
  }

  //read the sentinel character position
  elementsRead = fread(&indexData->sentinelCharacterPosition, sizeof(uint64_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }

  //read the bwt block list
  const size_t numBlockInBwt = awFmNumBlocksFromBwtLength(indexData->bwtLength);
  const size_t bytesPerBwtBlock = indexData->metadata.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);
  elementsRead = fread(indexData->bwtBlockList.asNucleotide, bytesPerBwtBlock, numBlockInBwt, fileHandle);
  if(elementsRead != numBlockInBwt){
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }
  //read the prefix sums array
  const size_t prefixSumsLength = awFmGetPrefixSumsLength(indexData->metadata.alphabetType);
  elementsRead = fread(indexData->prefixSums, sizeof(uint64_t), prefixSumsLength, fileHandle);
  if(elementsRead != prefixSumsLength){
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }

  //read the kmer seed table
  const size_t kmerSeedTableLength  = awFmGetKmerTableLength(indexData);
  elementsRead = fread(indexData->kmerSeedTable.table, sizeof(struct AwFmSearchRange), kmerSeedTableLength, fileHandle);
  if(elementsRead != kmerSeedTableLength){
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }

  const size_t sequenceEndingLength = indexData->metadata.kmerLengthInSeedTable - 1;
  elementsRead = fread(indexData->kmerSeedTable.sequenceEndingKmerEncodings, sizeof(uint64_t), sequenceEndingLength, fileHandle);
  if(elementsRead != sequenceEndingLength){
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }

  indexData->suffixArrayFileOffset  = awFmGetSuffixArrayFileOffset(indexData);
  indexData->sequenceFileOffset     = awFmGetSequenceFileOffset(indexData);
  indexData->fileHandle             = fileHandle;
  indexData->fileDescriptor         = fileno(fileHandle);

  *index = indexData;
  return AwFmFileReadOkay;
}


enum AwFmReturnCode awFmReadPositionsFromSuffixArray(const struct AwFmIndex *restrict const index,
  uint64_t *restrict const positionArray, const size_t positionArrayLength){

  for(size_t i = 0; i < positionArrayLength; i++){
    size_t indexInSuffixArray = positionArray[i] / index->metadata.suffixArrayCompressionRatio;
    uint64_t seekPosition = index->suffixArrayFileOffset + (indexInSuffixArray * sizeof(uint64_t));

    size_t bytesRead = pread(index->fileDescriptor, &(positionArray[i]), sizeof(uint64_t), seekPosition);
    if(bytesRead != sizeof(uint64_t)){
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
  const size_t sequenceSegmentLength = actualPriorFlankLength + postFlankLength;
  size_t bytesRead = pread(index->fileDescriptor, sequenceBuffer, sequenceSegmentLength * sizeof(char), seekPosition);

  if(bytesRead != sequenceSegmentLength * sizeof(char)){
    return AwFmFileReadFail;
  }

  //null terminate the string
  sequenceBuffer[sequenceSegmentLength] = 0;
  return AwFmFileReadOkay;
}


enum AwFmReturnCode awFmSuffixArrayReadPositionParallel(const struct AwFmIndex *restrict const index,
  struct AwFmBacktrace *restrict const backtracePtr){
  uint64_t suffixArrayPosition = backtracePtr->position / index->metadata.suffixArrayCompressionRatio;
  const size_t suffixArrayFileOffset = index->suffixArrayFileOffset + (suffixArrayPosition * sizeof(uint64_t));

  ssize_t numBytesRead = pread(index->fileDescriptor, &suffixArrayPosition, sizeof(uint64_t), suffixArrayFileOffset);
  if(numBytesRead == sizeof(uint64_t)){
    backtracePtr->position = suffixArrayPosition + backtracePtr->offset;
    return AwFmFileReadOkay;
  }
  else{
    backtracePtr->position = -1ULL;
    return AwFmFileReadFail;
  }
}


size_t awFmGetSequenceFileOffset(const struct AwFmIndex *restrict const index){
  const size_t metadataLength               = 5 * sizeof(uint8_t);
  const size_t bytesPerBwtBlock             = index->metadata.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);
  const size_t bwtLengthDataLength          = sizeof(uint64_t);
  const size_t sentinelPositionDataLength   = sizeof(uint64_t);
  const size_t bwtLengthInBytes             = awFmNumBlocksFromBwtLength(index->bwtLength) * bytesPerBwtBlock;
  const size_t prefixSumLengthInBytes       = awFmGetPrefixSumsLength(index->metadata.alphabetType) * sizeof(uint64_t);
  const size_t kmerSeedTableLength          = awFmGetKmerTableLength(index);
  const size_t endingSequenceLenghtInBytes  = (index->metadata.kmerLengthInSeedTable - 1) * sizeof(uint64_t);

  return IndexFileFormatIdHeaderLength + metadataLength +
    bwtLengthDataLength + sentinelPositionDataLength + bwtLengthInBytes +
    prefixSumLengthInBytes + (kmerSeedTableLength * sizeof(struct AwFmSearchRange) +
    endingSequenceLenghtInBytes);
}


size_t awFmGetSuffixArrayFileOffset(const struct AwFmIndex *restrict const index){
  return awFmGetSequenceFileOffset(index) + ((index->bwtLength - 1) * sizeof(char));
}
