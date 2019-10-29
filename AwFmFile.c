#include "AwFmIndex.h"
#include "AwFmFile.h"
#include <stdbool.h>
#include <stdio.h>
#include <string.h>


/*Private Function Prototypes*/

/*Returns the start position of the database sequence segment, or zero if the position would underflow*/
static inline size_t getSequenceSegmentStartPositionWithoutUnderflow(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankingSequenceLength);

/*Returns the end position of the database sequence segment, or sequence end position if the position would overflow*/
static inline size_t getSequenceSegmentEndPositionWithoutOverflow(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t postFlankingSequenceLength, const size_t dbSequenceLength);

/*Returns the number of bytes between the start of the AwFmIndex file and the start of the database sequence.*/
static inline size_t fileOffsetForSequenceBeginning(const struct AwFmIndex *restrict const index);

/*Returns the number of bytes between the start of the AwFmIndex file and the given beginning of the requested sequence segment.*/
static inline size_t fileOffsetForSequenceSegment(const struct AwFmIndex *restrict const index,
  const size_t sequenceSegmentStartPosition);

/*Returns the number of bytes between the start of the AwFmIndex file and the given position in the compressed suffix array.*/
static inline size_t fileOffsetForSuffixArrayPosition(const struct AwFmIndex *restrict const index,
  const size_t compressedSuffixArrayPosition);


static const uint8_t IndexFileFormatIdHeaderLength = 10;
static const char    IndexFileFormatIdHeader[10]  = "AwFmIndex\n";

/*
 * Function:  awFmCreateIndexFile
 * --------------------
 * Creates an AwFm-Index file from the given data. The file contains the index data structure,
 *  as well as the original database sequence and compressed suffix array.
 *  NOTE: This function requires that the database sequence and full suffix array
 *  be set in the AwFmIndex struct. This is done automatically when the index is
 *  constructed from a database sequence, but won't be done when loaded from a file.
 *
 *  Inputs:
 *    index:            prebuilt AwFmIndex struct to be saved to file.
 *    fullSuffixArray:  uncompresseed suffix array, that will be compressed by the
 *      ratio stored in the index's metadata
 *    databaseSequence: sequence used to build the given AwFmIndex and suffix array
 *    allowOverwrite:   flag that determines behavior if the index's fileSrc member data
 *      points to an existing file.
 *      If allowOverwrite is true, the file will be overwritten with the new index.
 *      If allowOverwrite is false, this function will return AwFmFileOpenFail.
 *
 *  Returns:
 *    AwFmReturnCode detailing the result of the write attempt.
 *      returns AwFmFileWriteOkay on success, AwFmFileOpenFail or AwFmFileWriteFail on failure.
 */
enum AwFmReturnCode awFmCreateIndexFile(const struct AwFmIndex *restrict const index, const bool allowOverwrite){

  //check to make sure the fileSrc seems valid
  if(index->fileSrc == NULL){
    return AwFmNoFileSrcGiven;
  }

  if(index->databaseSequence == NULL){
    return AwFmErrorDbSequenceNull;
  }

  if(index->fullSuffixArray == NULL){
    return AwFmErrorSuffixArrayNull;
  }
  const char *const fileOpenMode = allowOverwrite? "w": "wx";
  //open file to save to, and check for
  FILE *datafile = fopen(index->fileSrc, fileOpenMode);
  if(datafile == NULL){
    fclose(datafile);
    return AwFmFileOpenFail;
  }

  size_t elementsWritten;
  elementsWritten = fwrite(IndexFileFormatIdHeader, sizeof(char), IndexFileFormatIdHeaderLength, datafile);
  if(elementsWritten != IndexFileFormatIdHeaderLength){
    return AwFmFileWriteFail;
  }

  elementsWritten = fwrite(&index->metadata, sizeof(union AwFmIndexPaddedMetadata), 1, datafile);
  if(elementsWritten != 1){
    fclose(datafile);
    return AwFmFileWriteFail;
  }

  elementsWritten = fwrite(&index->numBlocks, sizeof(uint64_t), 1, datafile);
  if(elementsWritten != 1){
    fclose(datafile);
    return AwFmFileWriteFail;
  }

  elementsWritten = fwrite(index->rankPrefixSums, sizeof(uint64_t), AMINO_CARDINALITY + 1, datafile);
  if(elementsWritten != (AMINO_CARDINALITY + 1)){
    fclose(datafile);
    return AwFmFileWriteFail;
  }

  //write the compressed suffix array
  const uint8_t suffixArrayCompressionRatio = index->suffixArrayCompressionRatio;
  const uint64_t bwtLength = awFmGetBwtLength(index);

  for(size_t suffixArrayIndex = 0; suffixArrayIndex < bwtLength; suffixArrayIndex += suffixArrayCompressionRatio){
    elementsWritten = fwrite(index->fullSuffixArray + suffixArrayIndex, sizeof(uint64_t), 1, datafile);
    if(elementsWritten != 1){
      fclose(datafile);
      return AwFmFileWriteFail;
    }
  }

  //write the database sequence.
  elementsWritten = fwrite(index->databaseSequence, sizeof(char), awFmGetDbSequenceLength(index), datafile);
  if(elementsWritten != 1){
    fclose(datafile);
    return AwFmFileWriteFail;
  }

  fclose(datafile);
  return AwFmFileWriteOkay;
}


/*
 * Function:  awFmLoadIndexFromFile
 * --------------------
 * Loads an AwFmIndex struct from the given fileSrc. This function dynamically allocates memory for the index,
 *  so the index argument should not point to allocated memory before calling this function.
 *
 *  While the AwFmIndex file contains index, database sequence, and compressed suffix array, this function
 *    only reads the index struct, and does not read the database sequence or compressed suffix array.
 *
 *  Inputs:
 *    fileSrc:          Fully qualified file path to load the index from. the index's fileSrc member data
 *      will be set to this value. The fileSrc is not copied,
 *    index:            Pointer to the AwFmIndex pointer to be dynamically allocated and
 *      populated with the index from the given file.
 *
 *  Returns:
 *    AwFmReturnCode detailing the result of the read attempt. Possible return values:
 *      AwFmFileReadOkay on success,
 *      AwFmAllocationFailure on failure to allocate memory for the index struct
 *      AwFmFileFormatError on error caused by reading a file that does not look like the correct format,
 *      AwFmFileOpenFail on failure opening the file,
 *      AwFmFileWriteFail on failure writing to the file.
 */
enum AwFmReturnCode awFmLoadIndexFromFile(const char *restrict const fileSrc,
  struct AwFmIndex *restrict *restrict index){

  //make sure that the fileSrc was actually supplied
  if(fileSrc == NULL){
    return AwFmNoFileSrcGiven;
  }
  //allocate the index struct
  *index = awFmAlignedAllocAwFmIndex(); //TODO: also set the fileSrc here
  if(index == NULL){
    return AwFmAllocationFailure;
  }



  FILE *datafile = fopen(fileSrc, "r");
  if(datafile == NULL){
    return AwFmFileOpenFail;
  }

  size_t elementsRead;
  char headerLabelBuffer[IndexFileFormatIdHeaderLength];
  elementsRead = fread(headerLabelBuffer, sizeof(char), IndexFileFormatIdHeaderLength, datafile);
  if(elementsRead != sizeof(char) * IndexFileFormatIdHeaderLength){
    free(*index);
    fclose(datafile);
    return AwFmFileReadFail;
  }

  //check to make sure the first bytes are the expected label for this format
  if(strncmp(headerLabelBuffer, IndexFileFormatIdHeader, IndexFileFormatIdHeaderLength) != 0){
    free(*index);
    fclose(datafile);
    return AwFmFileFormatError;
  }

  elementsRead = fread(&(*index)->metadata, sizeof(union AwFmIndexPaddedMetadata), 1, datafile);
  if(elementsRead != 1){
    free(*index);
    fclose(datafile);
    return AwFmFileReadFail;
  }

  elementsRead = fread(&(*index)->numBlocks, sizeof(uint64_t), 1, datafile);
  if(elementsRead != 1){
    free(*index);
    fclose(datafile);
    return AwFmFileReadFail;
  }

  elementsRead = fread((*index)->rankPrefixSums, sizeof(uint64_t), AMINO_CARDINALITY + 1, datafile);
  if(elementsRead != (AMINO_CARDINALITY + 1)){
    free(*index);
    fclose(datafile);
    return AwFmFileReadFail;
  }

  //allocate the block list
  (*index)->blockList = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, sizeof(struct AwFmBlock) * (*index)->numBlocks);
  if((*index)->blockList ==  NULL){
    free(*index);
    fclose(datafile);
    return AwFmAllocationFailure;
  }

  elementsRead = fread((*index)->blockList, sizeof(struct AwFmBlock), (*index)->numBlocks, datafile);
  if(elementsRead != (*index)->numBlocks){
    free((*index)->blockList);
    free(*index);
    fclose(datafile);
    return AwFmFileReadFail;
  }

  fclose(datafile);
  return AwFmFileReadOkay;
}


/*
 * Function:  awFmDbSequencePositionsFromSuffixArrayFile
 * --------------------
 *  Loads the database sequence positions from the compressed suffix array for
 *    each element in the given positionArray.
 *
 *  This function will overwrite the data in the positionArray, returning the
 *    database sequence positions in the corresponding elements of the positionArray.
 *
 *  Note: positionArray and offsetArray should be of equal length. Having either shorter
 *    than positionArrayLength will result in undefined behavior.
 *
 *  Inputs:
 *    index:                  Pointer to the valid AwFmIndex struct so that we can determine the location in the file.
 *    positionArray:          Array of positions in the implied full suffix array to load.
 *      This function will convert these positions to indices in the compressed suffix array,
 *      as long as the positions are multiples of the suffix array compression ratio.
 *    offsetArray:  array of offsets to be added to the database sequence positions.
 *      This is needed because we can only query the suffix array on positions that are sampled.
 *    positionArrayLength:    Length of the positionArray and offsetArray.
 *
 *  Returns:
 *    AwFmReturnCode detailing the result of the read attempt. Possible return values:
 *      AwFmFileReadOkay on success,
 *      AwFmFileOpenFail on failure to open the AwFm Index file
 *      AwFmFileReadFail on failure to read as many characters as was expected by the sequence.
 *      AwFmIllegalPositionError on a suffix array position being out of bounds of
 *        the file's compressed suffix array.
 */
enum AwFmReturnCode awFmDbSequencePositionsFromSuffixArrayFile(const struct AwFmIndex *restrict const index,
  uint64_t *restrict positionArray, const uint64_t *restrict const offsetArray,
  const uint64_t positionArrayLength){

  if(index->fileSrc == NULL){
    return AwFmNoFileSrcGiven;
  }
  //open the file
  FILE *datafile = fopen( index->fileSrc, "r");
  if(datafile == NULL){
    return AwFmFileOpenFail;
  }

  const uint64_t compressedSuffixArrayLength = awFmGetCompressedSuffixArrayLength(index);
  for(uint64_t i = 0; i < positionArrayLength; i++){

    //for each position, find the index in the compressed suffix array.
    uint64_t compressedSuffixArrayPosition = positionArray[i] / index->suffixArrayCompressionRatio;
    if(compressedSuffixArrayPosition > compressedSuffixArrayLength){
      return AwFmIllegalPositionError;
    }

    //seek to the correct value in the compressed suffix array
    const size_t suffixFileOffset = fileOffsetForSuffixArrayPosition(index, compressedSuffixArrayPosition);
    fseek(datafile, suffixFileOffset, SEEK_SET);

    size_t elementsRead = fread(positionArray + i, sizeof(uint64_t), 1 , datafile);

    if(elementsRead != 1){
      fclose(datafile);
      return AwFmFileReadFail;
    }
    else{
      fclose(datafile);
      return AwFmFileReadOkay;
    }

    //add the offset to the position read from file
    positionArray[i] += offsetArray[i];
  }

  fclose(datafile);
  return AwFmFileReadOkay;

}


/*
 * Function:  awFmLoadSequenceSectionFromFile
 * --------------------
 * Loads a section of the database sequence from the AwFm Index file. If the specified range
 *  would start before the sequence's beginning, or would end after the sequence end, the range is trimmed
 *  to the sequence beginning or end
 *
 *
 *  Inputs:
 *    index:                        Pointer to the valid AwFmIndex struct so that we can determine the location in the file.
 *    sequencePosition:             position in the sequence segment to load, usually a hit from the FmIndex.
 *    priorFlankingSequenceLength:  How many characters to load before the given position.
 *    postFlankingSequenceLength:   How many characters to load after the given position.
 *    sequencePtr:                  pointer to the array to load the sequence into.
 *      NOTE: The buffer located at this pointer must be allocated by the caller, and
 *      the lifetime of the given buffer is the caller's responsibility.
 *    charactersRead:               output argument that returns the number of characters actually read.
 *      This may differ from the post-prior if either would violate the bounds of the sequence.
 *
 *  Returns:
 *    AwFmReturnCode detailing the result of the read attempt. Possible return values:
 *      AwFmFileReadOkay on success,
 *      AwFmFileOpenFail on failure to open the AwFm Index file
 *      AwFmFileReadFail on failure to read as many characters as was expected by the sequence.
 */
enum AwFmReturnCode awFmLoadSequenceSectionFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankingSequenceLength, const size_t postFlankingSequenceLength,
  char **sequencePtr, size_t *charactersRead){

  //check to make sure that the fileSrc is valid
  if(index->fileSrc == NULL){
    return AwFmNoFileSrcGiven;
  }

  //calculate the start position in the sequence, and how many characters to read.
  const size_t databaseSequenceLength       = awFmGetDbSequenceLength(index);
  const size_t sequenceSegmentStartPosition = getSequenceSegmentStartPositionWithoutUnderflow(index,
    sequencePosition, priorFlankingSequenceLength);
  const size_t sequenceSegmentEndPosition   = getSequenceSegmentEndPositionWithoutOverflow(index,
    sequencePosition, postFlankingSequenceLength, databaseSequenceLength);
  const size_t sequenceSegmentLength        = sequenceSegmentEndPosition - sequenceSegmentStartPosition;

  //find the byte offset from the start of the file where the sequence segment starts
  const size_t fileOffsetForSequenceSegmentStartPosition = fileOffsetForSequenceSegment(index,
    sequenceSegmentStartPosition);


  //open the file

  FILE *datafile = fopen(index->fileSrc, "r");
  if(datafile == NULL){
    return AwFmFileOpenFail;
  }

  //in the AwFmIndex file, seek to the beginning of the sequence segment
  int fseekResult = fseek(datafile, fileOffsetForSequenceSegmentStartPosition, SEEK_SET);

  if(fseekResult != 0){
    fclose(datafile);
    return AwFmFileReadFail;
  }

  *charactersRead = fread(*sequencePtr, sizeof(char), sequenceSegmentLength, datafile);
  if(*charactersRead != sequenceSegmentLength){
    fclose(datafile);
    return AwFmFileReadFail;
  }

  return AwFmFileReadOkay;
}


/*Private Function Implementations*/

/*
 * Function:  getSequenceSegmentStartPositionWithoutUnderflow
 * --------------------
 *  Returns the start position of the database sequence segment, or zero if the position would underflow.
 *    This function acts as a data-access safeguard for when a sequence segment is requested that would
 *    start before the actual start of the database sequence.
 *
 *  Inputs:
 *    index:                        Pointer to the AwFmIndex to be queried.
 *    sequencePosition:             Position in the sequence to build the segment around
 *    priorFlankingSequenceLength:  How many characters before the sequencePosition to include in the segment.
 *  Returns:
 *    priorFlankingSequenceLength if no underflow would occur,
 *      or sequencePosition - priorFlankingSequenceLength otherwise.
 */
size_t getSequenceSegmentStartPositionWithoutUnderflow(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankingSequenceLength){
    if(priorFlankingSequenceLength > sequencePosition)
      return 0;
    else
      return sequencePosition - priorFlankingSequenceLength;
}


/*
 * Function:  getSequenceSegmentEndPositionWithoutOverflow
 * --------------------
 *  Returns the end position of the database sequence segment, or sequence end position if the position would overflow.
 *    This function acts as a data-access safeguard for when a sequence segment is requested that would
 *    end after the actual end of the database sequence.
 *
 *  Inputs:
 *    index:                        Pointer to the AwFmIndex to be queried.
 *    sequencePosition:             Position in the sequence to build the segment around
 *    postFlankingSequenceLength:   How many characters after the sequencePosition to include in the segment.
 *  Returns:
 *    postFlankingSequenceLength if no overflow would occur,
 *      or postFlankingSequenceLength  - sequencePosition otherwise.
 */
size_t getSequenceSegmentEndPositionWithoutOverflow(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t postFlankingSequenceLength, const size_t dbSequenceLength){
    if((sequencePosition + postFlankingSequenceLength > dbSequenceLength))
      return dbSequenceLength - sequencePosition;
    else
      return postFlankingSequenceLength - sequencePosition;
}

/*
 * Function:  fileOffsetForSequenceBeginning
 * --------------------
 *  Returns the number of bytes between the start of the AwFmIndex file and the start of the database sequence.
 *
 *  Inputs:
 *    index:  Pointer to the AwFmIndex to be queried.
 *
 *  Returns:
 *    Number of bytes between the start of the AwFmIndex file and the start of the database sequence.
 *      While most of the AwFmIndex file is of static length, the number of blocks is variable,
 *      so this function accounts for the size of the block list.
 */
inline size_t fileOffsetForSequenceBeginning(const struct AwFmIndex *restrict const index){
  const uint_fast8_t  headerLength          = sizeof(char) * IndexFileFormatIdHeaderLength;
  const uint_fast16_t metadataLength        = sizeof(union AwFmIndexPaddedMetadata);
  const uint_fast8_t  numBlocksLength       = sizeof(uint64_t);
  const uint_fast16_t rankPrefixSumsLength  = sizeof(uint64_t) * (AMINO_CARDINALITY + 1);
  const uint64_t      blockListLength       = sizeof(struct AwFmBlock) * index->numBlocks;

  return headerLength + metadataLength + numBlocksLength + rankPrefixSumsLength + blockListLength;
}


/*
 * Function:  fileOffsetForSequenceSegment
 * --------------------
 *  Returns the number of bytes between the start of the AwFmIndex file and the given beginning of the requested sequence segment.
 *
 *  Inputs:
 *    index:                          Pointer to the AwFmIndex to be queried.
 *    sequenceSegmentStartPosition:   Position in the database sequence to start the segment.
 *
 *  Returns:
 *    Number of bytes between the start of the AwFmIndex file and the start of the sequence segment.
 *      While most of the AwFmIndex file is of static length, the number of blocks is variable,
 *      so this function accounts for the size of the block list (by calling fileOffsetForSequenceBeginning).
 */
inline size_t fileOffsetForSequenceSegment(const struct AwFmIndex *restrict const index, const size_t sequenceSegmentStartPosition){
  return fileOffsetForSequenceBeginning(index) + sequenceSegmentStartPosition;
}


/*
 * Function:  fileOffsetForSuffixArrayPosition
 * --------------------
 *  Returns the number of bytes between the start of the AwFmIndex file and the given position in the compressed suffix array.
 *
 *  Inputs:
 *    index:                          Pointer to the AwFmIndex to be queried.
 *    compressedSuffixArrayPosition:  Position in the compressed suffix array to load.
 *      Note: this is the index in the -compressed- suffix array, not the BWT position.
 *      To determine this index, use the awFmGetCompressedSuffixArrayLength function.
 *
 *  Returns:
 *    Number of bytes between the start of the AwFmIndex file and the requested position in the compressed suffix array.
 *      While most of the AwFmIndex file is of static length, the number of blocks is variable,
 *      so this function accounts for the size of the block list (by calling fileOffsetForSequenceBeginning).
 */
inline size_t fileOffsetForSuffixArrayPosition(const struct AwFmIndex *restrict const index,
  const size_t compressedSuffixArrayPosition){
  const size_t databaseSequenceSize = awFmGetDbSequenceLength(index) * sizeof(char);
  const size_t byteOffsetInSuffixArray = compressedSuffixArrayPosition * sizeof(uint64_t);

  return fileOffsetForSequenceBeginning(index) + databaseSequenceSize + byteOffsetInSuffixArray;

}
