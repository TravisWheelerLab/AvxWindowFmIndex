//_XOPEN_SOURCE needs to be defined before including unistd to make it include pread.
//I have absolutely no idea why this is needed, but the linter seems to think so...
#define _XOPEN_SOURCE 500

#include "AwFmFile.h"
#include "AwFmIndexStruct.h"
#include "AwFmIndex.h"
#include "AwFmSuffixArray.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

static const uint8_t IndexFileFormatIdHeaderLength  = 10;
static const char    IndexFileFormatIdHeader[11]    = "AwFmIndex\n\0";

enum AwFmReturnCode awFmWriteIndexToFile(struct AwFmIndex *restrict const index, const uint8_t *restrict const sequence,
  const uint64_t sequenceLength, const char *restrict const fileSrc, const bool allowOverwrite){
  if(__builtin_expect(fileSrc == NULL, 0)){
    return AwFmNoFileSrcGiven;
  }

  if(__builtin_expect(index == NULL, 0)){
    return AwFmNullPtrError;
  }

  if(__builtin_expect(sequence == NULL, 0)){
    return AwFmNullPtrError;
  }

  const bool indexContainsFastaVector = awFmIndexContainsFastaVector(index);
  const bool storeOriginalSequence = index->config.storeOriginalSequence;

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

  //write the config, starting with the version number
  elementsWritten = fwrite(&index->versionNumber, sizeof(uint32_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  //write the config, starting with the version number
  elementsWritten = fwrite(&index->featureFlags, sizeof(uint32_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  elementsWritten = fwrite(&index->config.suffixArrayCompressionRatio, sizeof(uint8_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  elementsWritten = fwrite(&index->config.kmerLengthInSeedTable, sizeof(uint8_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  uint8_t alphabetType = index->config.alphabetType;
  elementsWritten = fwrite(&alphabetType, sizeof(uint8_t), 1, index->fileHandle);
  if(elementsWritten != 1){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  elementsWritten = fwrite(&storeOriginalSequence, sizeof(uint8_t), 1, index->fileHandle);
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

  const size_t numBlockInBwt = awFmNumBlocksFromBwtLength(index->bwtLength);
  const size_t bytesPerBwtBlock = index->config.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

  elementsWritten = fwrite(index->bwtBlockList.asNucleotide, bytesPerBwtBlock, numBlockInBwt, index->fileHandle);
  if(elementsWritten != numBlockInBwt){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the prefix sums table
  const size_t prefixSumsLength = awFmGetPrefixSumsLength(index->config.alphabetType);
  elementsWritten = fwrite(index->prefixSums, sizeof(uint64_t), prefixSumsLength, index->fileHandle);
  if(elementsWritten != prefixSumsLength){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  //write the kmer seed table
  const size_t numElementsInKmerSeedTable = awFmGetKmerTableLength(index);
  elementsWritten = fwrite(index->kmerSeedTable, sizeof(struct AwFmSearchRange), numElementsInKmerSeedTable, index->fileHandle);
  if(elementsWritten != numElementsInKmerSeedTable){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  if(storeOriginalSequence){
    //write the sequence
    elementsWritten = fwrite(sequence, sizeof(char), sequenceLength, index->fileHandle);
    if(elementsWritten != sequenceLength){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
  }

  size_t downsampledSuffixArrayLengthInBytes = index->suffixArray.compressedByteLength;
  elementsWritten = fwrite(index->suffixArray.values, sizeof(uint8_t), downsampledSuffixArrayLengthInBytes, index->fileHandle);
  if(elementsWritten != downsampledSuffixArrayLengthInBytes){
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  if(indexContainsFastaVector){
    //write the lengths for the header string and the metadata vector
    const size_t headerStringLength = index->fastaVector->header.count;
    const size_t fastaVectorMetadataLength = index->fastaVector->metadata.count;
    size_t headerStringLengthWritten = fwrite(&headerStringLength, sizeof(size_t), 1, index->fileHandle);
    if(headerStringLengthWritten != 1){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
    size_t fastaVectorMetadataLengthWritten = fwrite(&fastaVectorMetadataLength, sizeof(size_t), 1, index->fileHandle);
    if(fastaVectorMetadataLengthWritten != 1){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
    size_t headerStringBytesWritten = fwrite(index->fastaVector->header.charData, sizeof(char), headerStringLength, index->fileHandle);
    if(headerStringBytesWritten != headerStringLength){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
    size_t fastaVectorMetadataElementsWritten = fwrite(index->fastaVector->metadata.data, sizeof(struct FastaVectorMetadata), fastaVectorMetadataLength, index->fileHandle);
    if(fastaVectorMetadataElementsWritten != fastaVectorMetadataLength){
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
  }

  fflush(index->fileHandle);
  index->fileDescriptor = fileno(index->fileHandle);

  return AwFmFileWriteOkay;
}


enum AwFmReturnCode awFmReadIndexFromFile(struct AwFmIndex *restrict *restrict index,
  const char *fileSrc, const bool keepSuffixArrayInMemory){

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
  struct AwFmIndexConfiguration config;
  uint32_t versionNumber;
  elementsRead = fread(&versionNumber, sizeof(uint32_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }

  //check to make sure the version is currently supported.
  if(!awFmIndexIsVersionValid(versionNumber)){
    fclose(fileHandle);
    return AwFmUnsupportedVersionError;
  }
  uint32_t featureFlags;
  elementsRead = fread(&featureFlags, sizeof(uint32_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  elementsRead = fread(&config.suffixArrayCompressionRatio, sizeof(uint8_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  elementsRead = fread(&config.kmerLengthInSeedTable, sizeof(uint8_t), 1, fileHandle);
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
  config.alphabetType = alphabetType;


  uint8_t storeOriginalSequence;

  elementsRead = fread(&storeOriginalSequence, sizeof(uint8_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  //boolean-ify the byte (should be 0 or 1 already, but just in case)
  config.storeOriginalSequence = !!storeOriginalSequence;


  //read the bwt length
  uint64_t bwtLength;
  elementsRead = fread(&bwtLength, sizeof(uint64_t), 1, fileHandle);
  if(elementsRead != 1){
    fclose(fileHandle);
    return AwFmFileReadFail;
  }

  //allocate the index
  indexData = awFmIndexAlloc(&config, bwtLength);
  if(indexData == NULL){
    fclose(fileHandle);
    return AwFmAllocationFailure;
  }
  indexData->versionNumber = versionNumber;

  indexData->featureFlags = featureFlags;
  //from the version number, determine if we expect to find FastaVector data.
  const bool indexContainsFastaVector = awFmIndexContainsFastaVector(indexData);

  //read the bwt block list
  const size_t numBlockInBwt = awFmNumBlocksFromBwtLength(indexData->bwtLength);
  const size_t bytesPerBwtBlock = indexData->config.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);
  elementsRead = fread(indexData->bwtBlockList.asNucleotide, bytesPerBwtBlock, numBlockInBwt, fileHandle);
  if(elementsRead != numBlockInBwt){
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }
  //read the prefix sums array
  const size_t prefixSumsLength = awFmGetPrefixSumsLength(indexData->config.alphabetType);
  elementsRead = fread(indexData->prefixSums, sizeof(uint64_t), prefixSumsLength, fileHandle);
  if(elementsRead != prefixSumsLength){
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }

  //read the kmer seed table
  const size_t kmerSeedTableLength  = awFmGetKmerTableLength(indexData);
  elementsRead = fread(indexData->kmerSeedTable, sizeof(struct AwFmSearchRange), kmerSeedTableLength, fileHandle);
  if(elementsRead != kmerSeedTableLength){
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }

  //handle the in memory suffix array, if requested.
  indexData->config.keepSuffixArrayInMemory = keepSuffixArrayInMemory;
  indexData->suffixArray.values = NULL;  //probably unnecessary, but here for safety.

  size_t compressedSuffixArrayByteLength = awFmComputeCompressedSaSizeInBytes(bwtLength, config.suffixArrayCompressionRatio);
  indexData->suffixArray.compressedByteLength = compressedSuffixArrayByteLength;
  indexData->suffixArray.valueBitWidth = awFmComputeSuffixArrayValueMinWidth(bwtLength);
  if(keepSuffixArrayInMemory){

    //seek to the suffix array
    fseek(fileHandle, awFmGetSuffixArrayFileOffset(indexData), SEEK_SET);
    indexData->suffixArray.values = malloc(compressedSuffixArrayByteLength);

    if(indexData->suffixArray.values == NULL){
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }
    else{
      elementsRead = fread(indexData->suffixArray.values, 1, compressedSuffixArrayByteLength, fileHandle);
      if(elementsRead != compressedSuffixArrayByteLength){
        fclose(fileHandle);
        awFmDeallocIndex(indexData);
        return AwFmFileReadFail;
      }
    }
  }
  if(indexContainsFastaVector){
    fseek(fileHandle, awFmGetFastaVectorFileOffset(indexData), SEEK_SET);
    //allocate and init the fastaVector struct
    struct FastaVector *fastaVector = malloc(sizeof(fastaVector));
    if(!fastaVector){
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
    }
    enum FastaVectorReturnCode fastaVectorReturnCode = fastaVectorInit(fastaVector);
    if(fastaVectorReturnCode == FASTA_VECTOR_ALLOCATION_FAIL){
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }

    //free the sequence buffer in the fastaVector, since it won't be used here
    free(fastaVector->sequence.charData);
    fastaVector->sequence.charData = NULL;

    indexData->fastaVector = fastaVector;
    size_t fastaVectorHeaderLength;
    size_t fastaVectorMetadataLength;
    elementsRead = fread(&fastaVectorHeaderLength, sizeof(size_t), 1, fileHandle);
    if(elementsRead != 1){
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }
    elementsRead = fread(&fastaVectorMetadataLength, sizeof(size_t), 1, fileHandle);
    if(elementsRead != 1){
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }

    //now, to do some hacking to the FastaVector struct
    fastaVector->header.charData = realloc(fastaVector->header.charData, fastaVectorHeaderLength * sizeof(char));
    if(!fastaVector->header.charData){
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }
    fastaVector->metadata.data = realloc(fastaVector->metadata.data, fastaVectorMetadataLength * sizeof(struct FastaVectorMetadata));
    if(!fastaVector->metadata.data){
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }
    fastaVector->header.count       = fastaVectorHeaderLength;
    fastaVector->header.capacity    = fastaVectorHeaderLength;
    fastaVector->metadata.count     = fastaVectorMetadataLength;
    fastaVector->metadata.capacity  = fastaVectorMetadataLength;
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

  if(index->config.keepSuffixArrayInMemory){
    for(size_t i = 0; i < positionArrayLength; i++){
      size_t indexInSuffixArray = positionArray[i] / index->config.suffixArrayCompressionRatio;
      size_t  position = awFmGetValueFromCompressedSuffixArray(index->suffixArray, indexInSuffixArray);
      positionArray[i] = position;
    }

    return AwFmSuccess;
  }
  else{
    for(size_t i = 0; i < positionArrayLength; i++){
      size_t indexInSuffixArray = positionArray[i] / index->config.suffixArrayCompressionRatio;
      enum AwFmReturnCode rc = awFmGetSuffixArrayValueFromFile(index, indexInSuffixArray, &positionArray[i]);
      if(rc != AwFmSuccess){
        return rc;
      }
    }

    return AwFmFileReadOkay;
  }
}


enum AwFmReturnCode awFmReadSequenceFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequenceStartPosition, const size_t sequenceSegmentLength,
  char *const sequenceBuffer){

  if(index->config.storeOriginalSequence){
    if(__builtin_expect((sequenceStartPosition + sequenceSegmentLength) > index->bwtLength, 0)){
      return AwFmIllegalPositionError;
    }

    const size_t seekPosition = index->sequenceFileOffset + sequenceStartPosition;

    size_t bytesRead = pread(index->fileDescriptor, sequenceBuffer, sequenceSegmentLength * sizeof(char), seekPosition);

    if(bytesRead != sequenceSegmentLength * sizeof(char)){
      return AwFmFileReadFail;
    }

    //null terminate the string
    sequenceBuffer[sequenceSegmentLength] = 0;
    return AwFmFileReadOkay;

  }
  else{

   return AwFmUnsupportedVersionError;
  }
}


enum AwFmReturnCode awFmSuffixArrayReadPositionParallel(const struct AwFmIndex *restrict const index,
  struct AwFmBacktrace *restrict const backtracePtr){

  if(index->config.keepSuffixArrayInMemory){
    uint64_t suffixArrayPosition = backtracePtr->position / index->config.suffixArrayCompressionRatio;

    size_t  saValue = awFmGetValueFromCompressedSuffixArray(index->suffixArray, suffixArrayPosition);
    backtracePtr->position = saValue + backtracePtr->offset;

    return AwFmSuccess;

  }
  else{
    uint64_t suffixArrayPosition = backtracePtr->position / index->config.suffixArrayCompressionRatio;
    size_t saValue;
    enum AwFmReturnCode rc = awFmGetSuffixArrayValueFromFile(index,suffixArrayPosition, &saValue);

    backtracePtr->position = saValue + backtracePtr->offset;
    return rc;
  }
}


enum AwFmReturnCode awFmGetSuffixArrayValueFromFile(const struct AwFmIndex *restrict const index, const size_t positionInArray, size_t *valueOut){

  //if we group up the non-standard bit width values into groups of 8, the lengths will always be a multiple of 8,
  // thus aligning to a byte. Then, we can add the num bytes of the remainder of values, and then find the bit offset.
  size_t alignedByteOffset = (positionInArray / 8) * index->suffixArray.valueBitWidth;
  size_t numEndingBits = (positionInArray % 8) * index->suffixArray.valueBitWidth;
  size_t bytePosition = alignedByteOffset + (numEndingBits / 8);
  size_t bitPosition = numEndingBits % 8;

  uint64_t buffer[2];
  size_t fileOffset = index->suffixArrayFileOffset + bytePosition;
  size_t numReadBytesRequired = (bitPosition + index->suffixArray.valueBitWidth + 7) / 8; //ceil to not round down
  memset(buffer, 0, sizeof(uint64_t)*2);

  ssize_t totalBytesRead = 0;
  do{
    //pread does not guarantee that it can read the entire pread request in a single transaction.
    //thus, we may need to loop until all bytes are read, or until we hit an error (negative return value)

    size_t currentReadFileOffset = fileOffset + totalBytesRead;
    size_t bytesInThisReadRequest = numReadBytesRequired - totalBytesRead;
    ssize_t bytesRead = pread(index->fileDescriptor, ((uint8_t*)buffer)+totalBytesRead, bytesInThisReadRequest, currentReadFileOffset);

    if(bytesRead <= 0){
      fclose(index->fileHandle);
      *valueOut = 0;
      return AwFmFileReadFail;
    }
    totalBytesRead += bytesRead;
  }while(totalBytesRead < numReadBytesRequired);

  //shift the bits of each element into place, and OR them together.
  buffer[0] = buffer[0] >> bitPosition;
  buffer[1] = buffer[1] << (64-bitPosition);
  buffer[0] = buffer[0] | buffer[1];
  uint64_t bitmask = (1 << index->suffixArray.valueBitWidth) - 1;

  *valueOut = buffer[0] & bitmask;
  return AwFmSuccess;
}


size_t awFmGetSequenceFileOffset(const struct AwFmIndex *restrict const index){
  const size_t configLength                 = 12 * sizeof(uint8_t);
  const size_t bytesPerBwtBlock             = index->config.alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);
  const size_t bwtLengthDataLength          = sizeof(uint64_t);
  const size_t bwtLengthInBytes             = awFmNumBlocksFromBwtLength(index->bwtLength) * bytesPerBwtBlock;
  const size_t prefixSumLengthInBytes       = awFmGetPrefixSumsLength(index->config.alphabetType) * sizeof(uint64_t);
  const size_t kmerSeedTableLength          = awFmGetKmerTableLength(index);

  return IndexFileFormatIdHeaderLength + configLength +
    bwtLengthDataLength + bwtLengthInBytes + prefixSumLengthInBytes +
    (kmerSeedTableLength * sizeof(struct AwFmSearchRange));
}


size_t awFmGetSuffixArrayFileOffset(const struct AwFmIndex *restrict const index){
  if(index->config.storeOriginalSequence){
    return awFmGetSequenceFileOffset(index) + ((index->bwtLength - 1) * sizeof(char));
  }
  else{
    return awFmGetSequenceFileOffset(index);
  }
}

size_t awFmGetFastaVectorFileOffset(const struct AwFmIndex *restrict const index){
  size_t compressedSuffixArrayByteLength = index->suffixArray.compressedByteLength;
  return awFmGetSuffixArrayFileOffset(index) + compressedSuffixArrayByteLength;
}
