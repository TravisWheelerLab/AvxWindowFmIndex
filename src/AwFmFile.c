//_XOPEN_SOURCE needs to be defined before including unistd to make it include
//pread.
// I have absolutely no idea why this is needed, but the linter seems to think
// so...
#define _XOPEN_SOURCE 500

#include "AwFmFile.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "AwFmIndex.h"
#include "AwFmIndexStruct.h"
#include "AwFmSuffixArray.h"

static const uint8_t IndexFileFormatIdHeaderLength = 10;
static const char IndexFileFormatIdHeader[11] = "AwFmIndex\n\0";

enum AwFmReturnCode
awFmWriteIndexToFile(struct AwFmIndex *_RESTRICT_ const index,
                     const uint8_t *_RESTRICT_ const sequence,
                     const uint64_t sequenceLength,
                     const char *_RESTRICT_ const fileSrc) {
  if (__builtin_expect(fileSrc == NULL, 0)) {
    return AwFmNoFileSrcGiven;
  }

  if (__builtin_expect(index == NULL, 0)) {
    return AwFmNullPtrError;
  }

  if (__builtin_expect(sequence == NULL, 0)) {
    return AwFmNullPtrError;
  }

  const bool indexContainsFastaVector = awFmIndexContainsFastaVector(index);
  const bool storeOriginalSequence = index->config.storeOriginalSequence;

  // open the file
  char fileOpenMode[4] = "w+b";
  index->fileHandle = fopen(fileSrc, fileOpenMode);

  if (index->fileHandle == NULL) {
    return AwFmFileAlreadyExists;
  }

  // write the header
  size_t elementsWritten =
      fwrite(IndexFileFormatIdHeader, sizeof(char),
             IndexFileFormatIdHeaderLength, index->fileHandle);
  if (elementsWritten != IndexFileFormatIdHeaderLength) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  // write the config, starting with the version number
  elementsWritten =
      fwrite(&index->versionNumber, sizeof(uint32_t), 1, index->fileHandle);
  if (elementsWritten != 1) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  // write the config, starting with the version number
  elementsWritten =
      fwrite(&index->featureFlags, sizeof(uint32_t), 1, index->fileHandle);
  if (elementsWritten != 1) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  elementsWritten = fwrite(&index->config.suffixArrayCompressionRatio,
                           sizeof(uint8_t), 1, index->fileHandle);
  if (elementsWritten != 1) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  elementsWritten = fwrite(&index->config.kmerLengthInSeedTable,
                           sizeof(uint8_t), 1, index->fileHandle);
  if (elementsWritten != 1) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  uint8_t alphabetType = index->config.alphabetType;
  elementsWritten =
      fwrite(&alphabetType, sizeof(uint8_t), 1, index->fileHandle);
  if (elementsWritten != 1) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }
  elementsWritten =
      fwrite(&storeOriginalSequence, sizeof(uint8_t), 1, index->fileHandle);
  if (elementsWritten != 1) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  // write the bwt length
  elementsWritten =
      fwrite(&index->bwtLength, sizeof(uint64_t), 1, index->fileHandle);
  if (elementsWritten != 1) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  const size_t numBlockInBwt = awFmNumBlocksFromBwtLength(index->bwtLength);
  const size_t bytesPerBwtBlock =
      index->config.alphabetType == AwFmAlphabetAmino
          ? sizeof(struct AwFmAminoBlock)
          : sizeof(struct AwFmNucleotideBlock);

  elementsWritten = fwrite(index->bwtBlockList.asNucleotide, bytesPerBwtBlock,
                           numBlockInBwt, index->fileHandle);
  if (elementsWritten != numBlockInBwt) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  // write the prefix sums table
  const size_t prefixSumsLength =
      awFmGetPrefixSumsLength(index->config.alphabetType);
  elementsWritten = fwrite(index->prefixSums, sizeof(uint64_t),
                           prefixSumsLength, index->fileHandle);
  if (elementsWritten != prefixSumsLength) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  // write the kmer seed table
  const size_t numElementsInKmerSeedTable = awFmGetKmerTableLength(index);
  elementsWritten = fwrite(index->kmerSeedTable, sizeof(struct AwFmSearchRange),
                           numElementsInKmerSeedTable, index->fileHandle);
  if (elementsWritten != numElementsInKmerSeedTable) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  if (storeOriginalSequence) {
    // write the sequence
    elementsWritten =
        fwrite(sequence, sizeof(char), sequenceLength, index->fileHandle);
    if (elementsWritten != sequenceLength) {
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
  }

  size_t downsampledSuffixArrayLengthInBytes =
      index->suffixArray.compressedByteLength;
  elementsWritten =
      fwrite(index->suffixArray.values, sizeof(uint8_t),
             downsampledSuffixArrayLengthInBytes, index->fileHandle);
  if (elementsWritten != downsampledSuffixArrayLengthInBytes) {
    fclose(index->fileHandle);
    return AwFmFileWriteFail;
  }

  if (indexContainsFastaVector) {
    // write the lengths for the header string and the metadata vector
    const size_t headerStringLength = index->fastaVector->header.count;
    const size_t fastaVectorMetadataLength = index->fastaVector->metadata.count;
    size_t headerStringLengthWritten =
        fwrite(&headerStringLength, sizeof(size_t), 1, index->fileHandle);
    if (headerStringLengthWritten != 1) {
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
    size_t fastaVectorMetadataLengthWritten = fwrite(
        &fastaVectorMetadataLength, sizeof(size_t), 1, index->fileHandle);
    if (fastaVectorMetadataLengthWritten != 1) {
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
    size_t headerStringBytesWritten =
        fwrite(index->fastaVector->header.charData, sizeof(char),
               headerStringLength, index->fileHandle);
    if (headerStringBytesWritten != headerStringLength) {
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
    size_t fastaVectorMetadataElementsWritten = fwrite(
        index->fastaVector->metadata.data, sizeof(struct FastaVectorMetadata),
        fastaVectorMetadataLength, index->fileHandle);
    if (fastaVectorMetadataElementsWritten != fastaVectorMetadataLength) {
      fclose(index->fileHandle);
      return AwFmFileWriteFail;
    }
  }

  fflush(index->fileHandle);
  index->fileDescriptor = fileno(index->fileHandle);

  return AwFmFileWriteOkay;
}

enum AwFmReturnCode
awFmReadIndexFromFile(struct AwFmIndex *_RESTRICT_ *_RESTRICT_ index,
                      const char *fileSrc, const bool keepSuffixArrayInMemory) {

  if (__builtin_expect(fileSrc == NULL, 0)) {
    return AwFmNoFileSrcGiven;
  }

  FILE *fileHandle = fopen(fileSrc, "r");
  if (!fileHandle) {
    return AwFmFileOpenFail;
  }

  // create a local-scope pointer for the index, when we're done we'll set the
  // index out-arg to this pointer.
  struct AwFmIndex *_RESTRICT_ indexData;

  // read the header, and check to make sure it matches
  char headerBuffer[IndexFileFormatIdHeaderLength + 1];
  size_t elementsRead = fread(headerBuffer, sizeof(char),
                              IndexFileFormatIdHeaderLength, fileHandle);
  if (elementsRead != IndexFileFormatIdHeaderLength) {
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  if (strncmp(IndexFileFormatIdHeader, headerBuffer,
              IndexFileFormatIdHeaderLength) != 0) {
    fclose(fileHandle);
    return AwFmFileFormatError;
  }

  // read in the metadata
  struct AwFmIndexConfiguration config;
  uint32_t versionNumber;
  elementsRead = fread(&versionNumber, sizeof(uint32_t), 1, fileHandle);
  if (elementsRead != 1) {
    fclose(fileHandle);
    return AwFmFileReadFail;
  }

  // check to make sure the version is currently supported.
  if (!awFmIndexIsVersionValid(versionNumber)) {
    fclose(fileHandle);
    return AwFmUnsupportedVersionError;
  }
  uint32_t featureFlags;
  elementsRead = fread(&featureFlags, sizeof(uint32_t), 1, fileHandle);
  if (elementsRead != 1) {
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  elementsRead = fread(&config.suffixArrayCompressionRatio, sizeof(uint8_t), 1,
                       fileHandle);
  if (elementsRead != 1) {
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  elementsRead =
      fread(&config.kmerLengthInSeedTable, sizeof(uint8_t), 1, fileHandle);
  if (elementsRead != 1) {
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  uint8_t alphabetType;
  elementsRead = fread(&alphabetType, sizeof(uint8_t), 1, fileHandle);
  if (elementsRead != 1) {
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  config.alphabetType = alphabetType;

  uint8_t storeOriginalSequence;

  elementsRead = fread(&storeOriginalSequence, sizeof(uint8_t), 1, fileHandle);
  if (elementsRead != 1) {
    fclose(fileHandle);
    return AwFmFileReadFail;
  }
  // boolean-ify the byte (should be 0 or 1 already, but just in case)
  config.storeOriginalSequence = !!storeOriginalSequence;

  // read the bwt length
  uint64_t bwtLength;
  elementsRead = fread(&bwtLength, sizeof(uint64_t), 1, fileHandle);
  if (elementsRead != 1) {
    fclose(fileHandle);
    return AwFmFileReadFail;
  }

  // allocate the index
  indexData = awFmIndexAlloc(&config, bwtLength);
  if (indexData == NULL) {
    fclose(fileHandle);
    return AwFmAllocationFailure;
  }
  indexData->versionNumber = versionNumber;

  indexData->featureFlags = featureFlags;
  // from the version number, determine if we expect to find FastaVector data.
  const bool indexContainsFastaVector = awFmIndexContainsFastaVector(indexData);

  // read the bwt block list
  const size_t numBlockInBwt = awFmNumBlocksFromBwtLength(indexData->bwtLength);
  const size_t bytesPerBwtBlock =
      indexData->config.alphabetType == AwFmAlphabetAmino
          ? sizeof(struct AwFmAminoBlock)
          : sizeof(struct AwFmNucleotideBlock);
  elementsRead = fread(indexData->bwtBlockList.asNucleotide, bytesPerBwtBlock,
                       numBlockInBwt, fileHandle);
  if (elementsRead != numBlockInBwt) {
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }
  // read the prefix sums array
  const size_t prefixSumsLength =
      awFmGetPrefixSumsLength(indexData->config.alphabetType);
  elementsRead = fread(indexData->prefixSums, sizeof(uint64_t),
                       prefixSumsLength, fileHandle);
  if (elementsRead != prefixSumsLength) {
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }

  // read the kmer seed table
  const size_t kmerSeedTableLength = awFmGetKmerTableLength(indexData);
  elementsRead = fread(indexData->kmerSeedTable, sizeof(struct AwFmSearchRange),
                       kmerSeedTableLength, fileHandle);
  if (elementsRead != kmerSeedTableLength) {
    fclose(fileHandle);
    awFmDeallocIndex(indexData);
    return AwFmFileReadFail;
  }

  // handle the in memory suffix array, if requested.
  indexData->config.keepSuffixArrayInMemory = keepSuffixArrayInMemory;
  indexData->suffixArray.values =
      NULL; // probably unnecessary, but here for safety.

  size_t compressedSuffixArrayByteLength = awFmComputeCompressedSaSizeInBytes(
      bwtLength, config.suffixArrayCompressionRatio);
  indexData->suffixArray.compressedByteLength = compressedSuffixArrayByteLength;
  indexData->suffixArray.valueBitWidth =
      awFmComputeSuffixArrayValueMinWidth(bwtLength);
  if (keepSuffixArrayInMemory) {

    // seek to the suffix array
    fseek(fileHandle, awFmGetSuffixArrayFileOffset(indexData), SEEK_SET);
    indexData->suffixArray.values = malloc(compressedSuffixArrayByteLength);

    if (indexData->suffixArray.values == NULL) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    } else {
      elementsRead = fread(indexData->suffixArray.values, 1,
                           compressedSuffixArrayByteLength, fileHandle);
      if (elementsRead != compressedSuffixArrayByteLength) {
        fclose(fileHandle);
        awFmDeallocIndex(indexData);
        return AwFmFileReadFail;
      }
    }
  }
  if (indexContainsFastaVector) {
    fseek(fileHandle, awFmGetFastaVectorFileOffset(indexData), SEEK_SET);
    // allocate and init the fastaVector struct
    struct FastaVector *fastaVector = malloc(sizeof(struct FastaVector));
    if (!fastaVector) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }
    enum FastaVectorReturnCode fastaVectorReturnCode =
        fastaVectorInit(fastaVector);
    if (fastaVectorReturnCode == FASTA_VECTOR_ALLOCATION_FAIL) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }

    // free the sequence buffer in the fastaVector, since it won't be used here
    fastaVectorStringDealloc(&fastaVector->sequence);
    fastaVector->sequence.charData = NULL;
    fastaVector->sequence.capacity = 0;
    fastaVector->sequence.count = 0;

    indexData->fastaVector = fastaVector;
    size_t fastaVectorHeaderLength;
    size_t fastaVectorMetadataLength;
    elementsRead =
        fread(&fastaVectorHeaderLength, sizeof(size_t), 1, fileHandle);
    if (elementsRead != 1) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmFileReadFail;
    }
    elementsRead =
        fread(&fastaVectorMetadataLength, sizeof(size_t), 1, fileHandle);
    if (elementsRead != 1) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmFileReadFail;
    }

    // now, to do some hacking to the FastaVector struct
    fastaVector->header.charData = realloc(
        fastaVector->header.charData, fastaVectorHeaderLength * sizeof(char));
    if (!fastaVector->header.charData) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }

    elementsRead = fread(fastaVector->header.charData, sizeof(char),
                         fastaVectorHeaderLength, fileHandle);
    if (elementsRead != fastaVectorHeaderLength) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmFileReadFail;
    }

    fastaVector->metadata.data =
        realloc(fastaVector->metadata.data,
                fastaVectorMetadataLength * sizeof(struct FastaVectorMetadata));
    if (!fastaVector->metadata.data) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmAllocationFailure;
    }

    elementsRead =
        fread(fastaVector->metadata.data, sizeof(struct FastaVectorMetadata),
              fastaVectorMetadataLength, fileHandle);
    if (elementsRead != fastaVectorMetadataLength) {
      fclose(fileHandle);
      awFmDeallocIndex(indexData);
      return AwFmFileReadFail;
    }

    fastaVector->header.count = fastaVectorHeaderLength;
    fastaVector->header.capacity = fastaVectorHeaderLength;
    fastaVector->metadata.count = fastaVectorMetadataLength;
    fastaVector->metadata.capacity = fastaVectorMetadataLength;
  }

  indexData->suffixArrayFileOffset = awFmGetSuffixArrayFileOffset(indexData);
  indexData->sequenceFileOffset = awFmGetSequenceFileOffset(indexData);
  indexData->fileHandle = fileHandle;
  indexData->fileDescriptor = fileno(fileHandle);

  *index = indexData;
  return AwFmFileReadOkay;
}

enum AwFmReturnCode
awFmReadSequenceFromFile(const struct AwFmIndex *_RESTRICT_ const index,
                         const size_t sequenceStartPosition,
                         const size_t sequenceSegmentLength,
                         char *const sequenceBuffer) {

  if (index->config.storeOriginalSequence) {
    if (__builtin_expect((sequenceStartPosition + sequenceSegmentLength) >
                             index->bwtLength,
                         0)) {
      return AwFmIllegalPositionError;
    }

    const size_t seekPosition =
        index->sequenceFileOffset + sequenceStartPosition;

    size_t bytesRead =
        pread(index->fileDescriptor, sequenceBuffer,
              sequenceSegmentLength * sizeof(char), seekPosition);

    if (bytesRead != sequenceSegmentLength * sizeof(char)) {
      return AwFmFileReadFail;
    }

    // null terminate the string
    sequenceBuffer[sequenceSegmentLength] = 0;
    return AwFmFileReadOkay;
  } else {

    return AwFmUnsupportedVersionError;
  }
}

enum AwFmReturnCode
awFmGetSuffixArrayValueFromFile(const struct AwFmIndex *_RESTRICT_ const index,
                                const size_t positionInArray,
                                size_t *valueOut) {

  struct AwFmSuffixArrayOffset offset = awFmGetOffsetIntoSuffixArrayByteArray(
      index->suffixArray.valueBitWidth, positionInArray);

  const size_t suffixArrayFileOffset = index->suffixArrayFileOffset;
  uint8_t valueBuffer[9] = {
      0}; // 9 bytes are required for this buffer, since ceil(((bitOffset=7) +
          // (valBitWidth=64)) / 8) = 9
  int8_t bytesToRead =
      (offset.bitOffset + index->suffixArray.valueBitWidth + 7) /
      8; // rounded up.
  int8_t totalBytesRead = 0;

  while (bytesToRead > totalBytesRead) {
    int8_t bytesLeft = bytesToRead - totalBytesRead;
    size_t readPosition =
        suffixArrayFileOffset + offset.byteOffset + totalBytesRead;
    ssize_t bytesRead =
        pread(index->fileDescriptor, valueBuffer + totalBytesRead, bytesLeft,
              readPosition);
    if (bytesRead <= 0) {
      return AwFmFileWriteFail;
    }
    totalBytesRead += bytesRead;
  }

  // reconstruct the size_t valueOut from the valueBuffer.
  // first memcpy the first 8 bytes into valueOut
  memcpy(valueOut, valueBuffer, 8);
  *valueOut >>= offset.bitOffset;
  *valueOut |= ((uint64_t)valueBuffer[8]) << (64 - offset.bitOffset);
  size_t bitmask = (((size_t)1) << index->suffixArray.valueBitWidth) - 1;
  *valueOut &= bitmask;
  return AwFmSuccess;
}

size_t
awFmGetSequenceFileOffset(const struct AwFmIndex *_RESTRICT_ const index) {
  const size_t configLength = 12 * sizeof(uint8_t);
  const size_t bytesPerBwtBlock =
      index->config.alphabetType == AwFmAlphabetAmino
          ? sizeof(struct AwFmAminoBlock)
          : sizeof(struct AwFmNucleotideBlock);
  const size_t bwtLengthDataLength = sizeof(uint64_t);
  const size_t bwtLengthInBytes =
      awFmNumBlocksFromBwtLength(index->bwtLength) * bytesPerBwtBlock;
  const size_t prefixSumLengthInBytes =
      awFmGetPrefixSumsLength(index->config.alphabetType) * sizeof(uint64_t);
  const size_t kmerSeedTableLength = awFmGetKmerTableLength(index);

  return IndexFileFormatIdHeaderLength + configLength + bwtLengthDataLength +
         bwtLengthInBytes + prefixSumLengthInBytes +
         (kmerSeedTableLength * sizeof(struct AwFmSearchRange));
}

size_t
awFmGetSuffixArrayFileOffset(const struct AwFmIndex *_RESTRICT_ const index) {
  if (index->config.storeOriginalSequence) {
    return awFmGetSequenceFileOffset(index) +
           ((index->bwtLength - 1) * sizeof(char));
  } else {
    return awFmGetSequenceFileOffset(index);
  }
}

size_t
awFmGetFastaVectorFileOffset(const struct AwFmIndex *_RESTRICT_ const index) {
  size_t compressedSuffixArrayByteLength =
      index->suffixArray.compressedByteLength;
  return awFmGetSuffixArrayFileOffset(index) + compressedSuffixArrayByteLength;
}
