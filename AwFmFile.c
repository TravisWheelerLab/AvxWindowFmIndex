#include "AwFmFile.h"
#include "AwFmIndex.h"
#include <stdbool.h>
#include <string.h>


/*Private Function Prototypes*/
static inline size_t getSequenceSegmentStartPositionWithoutUnderflow(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankingSequenceLength, const size_t dbSequenceLength);

static inline size_t getSequenceSegmentEndPositionWithoutOverflow(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t postFlankingSequenceLength, const size_t dbSequenceLength);

static inline size_t fileOffsetForSequenceBeginning(const struct AwFmIndex *restrict const index);

static inline size_t fileOffsetForSequenceSegment(const struct AwFmIndex *restrict const index,
  const size_t sequenceSegmentStartPosition);

static inline size_t fileOffsetForSuffixArrayPosition(const struct AwFmIndex *restrict const index,
  const size_t compressedSuffixArrayPosition);



static const char    IndexFileExtension[5]         = "awfmi";
static const uint8_t IndexFileFormatIdHeaderLength = 10
static const char    IndexFileFormatIdHeader[IndexFileFormatIdHeaderLength]  = "AwFmIndex\n";

/*
 * Function:  awFmCreateIndexFile
 * --------------------
 * Creates an AwFm-Index file from the given data. The file contains the index data structure,
 *  as well as the original database sequence and compressed suffix array.
 *
 *  Inputs:
 *    fileSrc:          Fully qualified file path to save the index to.
 *    index:            prebuilt AwFmIndex struct to be saved to file.
 *    fullSuffixArray:  uncompresseed suffix array, that will be compressed by the
 *      ratio stored in the index's metadata
 *    databaseSequence: sequence used to build the given AwFmIndex and suffix array
 *    allowOverwrite:   flag that determines behavior if the given fileSrc points to an existing file.
 *      If allowOverwrite is true, the file will be overwritten with the new index.
 *      If allowOverwrite is false, this function will return AwFmFileOpenFail.
 *
 *  Returns:
 *    AwFmFileAccessCode detailing the result of the write attempt.
 *      returns AwFmFileWriteOkay on success, AwFmFileOpenFail or AwFmFileWriteFail on failure.
 */
enum AwFmFileAccessCode awFmCreateIndexFile(const char *restrict const fileSrc,
  const struct AwFmIndex *restrict const index, const size_t *restrict const fullSuffixArray,
  const char *restrict const databaseSequence, const bool allowOverwrite){

  FILE *datafile;

  const char *const fileOpenMode = allowOverwrite? "w": "wx";
  //open file to save to, and check for
  errno_t errorCode = fopen_s(&datafile, fileSrc, fileOpenMode);
  if(res != 0){
    fclose(datafile);
    return AwFmFileOpenFail;
  }

  size_t elementsWritten;
  elementsWritten = fwrite(IndexFileFormatIdHeader, sizeof(char), IndexFileFormatIdHeaderLength, datafile);
  if(elementsWritten != IndexFileFormatIdHeaderLength){
    return AwFmFileWriteFail;
  }

  elementsWritten = fwrite(index->metadata, sizeof(struct AwFmIndexPaddedMetadata), 1, datafile);
  if(elementsWritten != 1){
    fclose(datafile);
    return AwFmFileWriteFail;
  }

  elementsWritten = fwrite(index->numBlocks, sizeof(uint64_t), 1, datafile);
  if(elementsWritten != 1){
    fclose(datafile);
    return AwFmFileWriteFail;
  }

  elementsWritten = fwrite(index->rankPrefixSums, sizeof(uint64_t), AMINO_CARDINALITY + 1, datafile);
  if(elementsWritten != (AMINO_CARDINALITY + 1)){
    fclose(datafile);
    return AwFmFileWriteFail;
  }

  //write the database sequence
  const uint64_t dbSequenceLength = bwtLength - 1; //subtract 1 since the original sequence doesn't have a null terminator.
  fwrite(databaseSequence, sizeof(char), dbSequenceLength, datafile);
  if(elementsWritten != dbSequenceLength){
    fclose(datafile);
    return AwFmFileWriteFail;
  }

  //write the suffix array
  const uint8_t suffixArrayCompressionRatio = index->metadata.data.suffixArrayCompressionRatio;
  const uint64_t bwtLength = awFmGetBwtLength(index);

  for(size_t suffixArrayIndex = 0; suffixArrayIndex < bwtLength; suffixArrayIndex += suffixArrayCompressionRatio){
    fwrite(fullSuffixArray + suffixArrayIndex, sizeof(uint64_t), 1, datafile);
    if(elementsWritten != 1){
      fclose(datafile);
      return AwFmFileWriteFail;
    }
  }

  fclose();
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
 *    fileSrc:          Fully qualified file path to load the index from.
 *    index:            Pointer to the AwFmIndex pointer to be dynamically allocated and
 *      populated with the index from the given file.
 *
 *  Returns:
 *    AwFmFileAccessCode detailing the result of the read attempt. Possible return values:
 *      AwFmFileReadOkay on success,
 *      AwFmAllocationFailure on failure to allocate memory for the index struct
 *      AwFmFileFormatError on error caused by reading a file that does not look like the correct format,
 *      AwFmFileOpenFail on failure opening the file,
 *      AwFmFileWriteFail on failure writing to the file.
 */
enum AwFmFileAccessCode awFmLoadIndexFromFile(const char *restrict const fileSrc,
  const struct AwFmIndex *restrict *restrict index){
  //allocate the index struct
  *index = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, sizeof(index));
  if(index == NULL){
    return AwFmAllocationFailure;
  }

  FILE *datafile;
  errno_t errorCode = fopen_s(&datafile, fileSrc, "r");
  if(errorCode != 0){
    free(*index);
    return AwFmFileOpenFail;
  }

  size_t elementsRead;
  char headerLabelBuffer[IndexFileFormatIdHeaderLength];
  elementsRead = fread(headerLabelBuffer, sizeof(char), IndexFileFormatIdHeaderLength, datafile);
  if(elementsRead != sizeof(char) * IndexFileFormatLabelLength){
    free(*index);
    fclose(datafile);
    retrun AwFmFileReadFail;
  }

  //check to make sure the first bytes are the expected label for this format
  if(strncmp(headerLabelBuffer, IndexFileFormatIdHeader, fmIndexFileFormatHeadLength) != 0){
    free(*index);
    fclose(datafile);
    return AwFmFileFormatError;
  }

  elementsRead = fread((*index)->metadata, sizeof(struct AwFmIndexPaddedMetadata), 1, datafile);
  if(elementsRead != 1){
    free(*index);
    fclose(datafile);
    retrun AwFmFileReadFail;
  }

  elementsRead = fread((*index)->numBlocks, sizeof(uint64_t), 1, datafile);
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

  elementsRead = fread((*index->blockList), sizeof(struct AwFmBlock), (*index)->numBlocks);
  if(elementsRead != (*index)->numBlocks){
    free((*index)->blockList);
    free(*index);
    fclose(datafile);
    retrun AwFmFileReadFail;
  }

  fclose(datafile);
  return AwFmFileReadOkay;
}


/*
 * Function:  awFmLoadSuffixArrayIndexFromFile
 * --------------------
 * Loads the position in the database sequence for the given compressed suffix array index from the AwFm Index file.
 *
 *
 *  Inputs:
 *    fileSrc:                        Fully qualified file path to load the database sequence from.
 *    index:                          Pointer to the valid AwFmIndex struct so that we can determine the location in the file.
 *    compressedSuffixArrayPosition:  Index into the compressed suffix array to load. Note, this is not the same as the BWT position.
 *      Since the SA is likely compressed, the BWT positon is the compressedSuffixArrayPosition multiplied by the compression factor.
 *    dbPositionOut:                  Output argument for the position in the database sequence.
 *
 *  Returns:
 *    AwFmFileAccessCode detailing the result of the read attempt. Possible return values:
 *      AwFmFileReadOkay on success,
 *      AwFmFileOpenFail on failure to open the AwFm Index file.
 *      AwFmFileReadFail on failure to read the data at the expected file position.
 *        This would be likely due to the file somehow being smaller than the index expects.
 *      AwFmIllegalPositionError on the compressedSuffixArrayPosition being outisde the range that
 *        would be stored in the AwFm Index file.
 */
enum AwFmFileAccessCode awFmLoadSuffixArrayIndexFromFile(const char *restrict const fileSrc,
  const struct AwFmIndex *restrict const index, const size_t compressedSuffixArrayPosition,
  const uint64_t *restrict dbPositionOut){

  //sanity check on the compressed SA index
  if((compressedSuffixArrayPosition * index->metadata.data.suffixArrayCompressionRatio) > awFmGetDbSequenceLength(index)){
    return AwFmIllegalPositionError;
  }

  //open the file
  FILE *datafile;
  errno_t errorCode = fopen_s(&datafile, "r");
  if(errorCode != 0){
    return AwFmFileOpenFail;
  }

  //seek to the correct value in the compressed suffix array
  const size_t suffixFileOffset = fileOffsetForSuffixArrayPosition(index, compressedSuffixArrayPosition);
  fseek(datafile, suffixFileOffset, SEEK_SET);

  size_t elementsRead = fread(dbPositionOut, sizeof(uint64_t), 1 , datafile);
  fclose(datafile);

  if(elementsRead != 1)
    return AwFmFileReadFail;
  else
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
 *    fileSrc:                      Fully qualified file path to load the database sequence from.
 *    index:                        Pointer to the valid AwFmIndex struct so that we can determine the location in the file.
 *    sequencePosition:             position in the sequence segment to load, usually a hit from the FmIndex.
 *    priorFlankingSequenceLength:  How many characters to load before the given position.
 *    postFlankingSequenceLength:   How many characters to load after the given position.
 *    sequencePtr:                  pointer to the array to load the sequence into.
 *    charactersRead:               output argument that returns the number of characters actually read.
 *      This may differ from the post-prior if either would violate the bounds of the sequence.
 *
 *  Returns:
 *    AwFmFileAccessCode detailing the result of the read attempt. Possible return values:
 *      AwFmFileReadOkay on success,
 *      AwFmFileOpenFail on failure to open the AwFm Index file
 *      AwFmFileReadFail on failure to read as many characters as was expected by the sequence.
 */
enum AwFmFileAccessCode awFmLoadSequenceSectionFromFile(const char *restrict const fileSrc,
  const struct AwFmIndex *restrict const index, const size_t sequencePosition,
  const size_t priorFlankingSequenceLength, const size_t postFlankingSequenceLength,
  const char **sequencePtr, const size_t *charactersRead){

  //calculate the start position in the sequence, and how many characters to read.
  const size_t databaseSequenceLength = awFmGetDbSequenceLength(index);
  const size_t sequenceSegmentStartPosition  = getSequenceSegmentStartPositionWithoutUnderflow(index,
    sequencePosition, priorFlankingSequenceLength, databaseSequenceLength);
  const size_t sequenceSegmentEndPosition    = getSequenceSegmentEndPositionWithoutOverflow(index,
    sequencePosition, postFlankingSequenceLength, databaseSequenceLength);
  const size_t sequenceSegmentLength = sequenceEndPosition - sequenceStartPosition;

  //find the byte offset from the start of the file where the sequence segment starts
  const size_t fileOffsetForSequenceSegmentStartPosition = fileOffsetForSequenceSegment(index->numBlocks,
    sequenceSegmentStartPosition);


  //open the file
  FILE *datafile;
  errno_t errorCode = fopen_s(&datafile, "r");
  if(errorCode != 0){
    return AwFmFileOpenFail;
  }

  //in the AwFmIndex file, seek to the beginning of the sequence segment
  int fseekResult = fseek(datafile, fileOffsetForSequenceSegmentStartPosition, SEEK_SET);

  if(fseekResult != 0){
    fclose(datafile);
    return AwFmFileReadFail;
  }

  *charactersRead = fread(*sequencePtr, sizeof(char), charactersToRead, datafile);
  if(*charactersRead != sequenceSegmentLength){
    fclose(datafile);
    return AwFmFileReadFail;
  }

  return AwFmFileReadOkay;
}



/*private function implementations*/
size_t getSequenceSegmentStartPositionWithoutUnderflow(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankingSequenceLength, const size_t dbSequenceLength){
    if(priorFlankingSequenceLength > sequencePosition)
      return 0;
    else
      return sequencePosition - priorFlankingSequenceLength
}

size_t getSequenceSegmentEndPositionWithoutOverflow(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t postFlankingSequenceLength, const size_t dbSequenceLength){
    if((sequencePosition + postFlankingSequenceLength > dbSequenceLength))
      return dbSequenceLength - sequencePosition;
    else
      return postFlankingSequenceLength - sequencePosition;
}

inline size_t fileOffsetForSequenceBeginning(const struct AwFmIndex *restrict const index){
  const uint_fast8_t  headerLength          = sizeof(char) * IndexFileFormatIdHeaderLength;
  const uint_fast16_t metadataLength        = sizeof(struct AwFmIndexPaddedMetadata);
  const uint_fast8_t  numBlocksLength       = sizeof(uint64_t);
  const uint_fast16_t rankPrefixSumsLength  = sizeof(uint64_t) * (AMINO_CARDINALITY + 1);
  const uint_64_t     blockListLength       = sizeof(struct AwFmBlock) * index->numBlocks

  return headerLength + metadataLength + numBlocksLength + rankPrefixSumsLength + blockListLength;
}

inline size_t fileOffsetForSequenceSegment(const struct AwFmIndex *restrict const index, const size_t sequenceSegmentStartPosition){
  return fileOffsetForSequenceBeginning(index) + sequenceSegmentStartPosition;
}

inline size_t fileOffsetForSuffixArrayPosition(const struct AwFmIndex *restrict const index,
  const size_t compressedSuffixArrayPosition){
  const size_t databaseSequenceSize = awFmGetDbSequenceLength(index) * sizeof(char);
  const size_t byteOffsetInSuffixArray = suffixArrayPosition * sizeof(uint64_t);

  return fileOffsetForSequenceBeginning(index) + databaseSequenceSize + byteOffsetInSuffixArray;

}
