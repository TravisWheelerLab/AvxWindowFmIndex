#include "AwFmIndex.h"
#include "AwFmIndexStruct.h"
#include "AwFmCreate.h"
#include "AwFmLetter.h"
#include "AwFmSearch.h"
#include "AwFmSimdConfig.h"
#include "AwFmFile.h"
#include "AwFmIndexStruct.h"
#include <stdlib.h>
#include <string.h>
#include "divsufsort64.h"
#include "FastaVector.h"



/*private function prototypes*/
void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t sequenceLength,
  const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray);

void populateKmerSeedTable(struct AwFmIndex *restrict const index);

void populateKmerSeedTableRecursive(struct AwFmIndex *restrict const index, struct AwFmSearchRange range,
  uint8_t currentKmerLength, uint64_t currentKmerIndex, uint64_t letterIndexMultiplier);

void compressSuffixArrayInPlace(uint64_t *const suffixArray, uint64_t suffixArrayLength, const uint8_t compressionRatio);

void fullSequenceSanitize(const uint8_t *const sequence, uint8_t *const sanitizedSequenceCopy,
  const size_t sequenceLength, const enum AwFmAlphabetType alphabetType);


/*function implementations*/
enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex *restrict *index,
  struct AwFmIndexMetadata *restrict const metadata, const uint8_t *restrict const sequence, const size_t sequenceLength,
  const char *restrict const fileSrc, const bool allowFileOverwrite){

  //first, do a sanity check on inputs
  if(metadata == NULL){
    return AwFmNullPtrError;
  }
  if(sequence == NULL){
    return AwFmNullPtrError;
  }
  if(fileSrc == NULL){
    return AwFmNullPtrError;
  }

  metadata->versionNumber = 1;

  //set the index out arg initally to NULL, if this function fully completes this will get overwritten
  *index = NULL;

  const size_t suffixArrayLength = sequenceLength + 1;
  //create a sanitized copy of the input sequence
  uint8_t *sanitizedSequenceCopy = malloc(suffixArrayLength);
  if(sanitizedSequenceCopy == NULL){
    return AwFmAllocationFailure;
  }

  //sanitize the sequence, turning ambiguity characters into the singular ambiguity character
  //character (x for nucleotide, z for amino)
  fullSequenceSanitize(sequence, sanitizedSequenceCopy, sequenceLength, metadata->alphabetType);

  // append the final sentinel character as a terminator.
  sanitizedSequenceCopy[suffixArrayLength - 1] = '$';

  //allocate the index and all internal arrays.
  struct AwFmIndex *restrict indexData = awFmIndexAlloc(metadata, suffixArrayLength);
  if(indexData == NULL){
    return AwFmAllocationFailure;
  }
  indexData->fastaVector = NULL;  //set the fastaVector struct to null, since we aren't using it for this version.
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
    //if the suffix array isn't supposed to be kept in memory, there's no need to keep it, so it gets freed here.
    indexData->inMemorySuffixArray = NULL;
    free(suffixArray);
  }
  //set the index as an out argument.
  *index = indexData;

  return returnCode;
}

/*function implementations*/
enum AwFmReturnCode awFmCreateIndexFromFasta(struct AwFmIndex *restrict *index,
  struct AwFmIndexMetadata *restrict const metadata, const char *fastaSrc,
  const char *restrict const indexFileSrc, const bool allowFileOverwrite){

  //first, do a sanity check on inputs
  if(metadata == NULL){
    return AwFmNullPtrError;
  }
  if(fastaSrc == NULL){
    return AwFmNullPtrError;
  }
  if(indexFileSrc == NULL){
    return AwFmNullPtrError;
  }

  metadata->versionNumber = 2;
  //set the index out arg initally to NULL, if this function fully completes this will get overwritten
  *index = NULL;

  //read in the sequence file with FastaVector
  struct FastaVector *fastaVector = malloc(sizeof(struct FastaVector));
  if(fastaVector == NULL){
    return AwFmAllocationFailure;
  }
  enum FastaVectorReturnCode fastaVectorReturnCode = fastaVectorInit(fastaVector);
  if(fastaVectorReturnCode == FASTA_VECTOR_ALLOCATION_FAIL){
    return AwFmAllocationFailure;
  }
  fastaVectorReturnCode = fastaVectorReadFasta(fastaSrc, fastaVector);
  if(fastaVectorReturnCode == FASTA_VECTOR_FILE_OPEN_FAIL){
    return AwFmFileOpenFail;
  }
  else if (fastaVectorReturnCode == FASTA_VECTOR_ALLOCATION_FAIL){
    return AwFmAllocationFailure;
  }

  const size_t fullSequenceLength = fastaVector->sequence.count;
  const char *fullSequencePtr     = fastaVector->sequence.charData;

  const size_t suffixArrayLength = fullSequenceLength + 1;
  //create a sanitized copy of the input sequence
  uint8_t *sanitizedSequenceCopy = malloc(suffixArrayLength);
  if(sanitizedSequenceCopy == NULL){
    return AwFmAllocationFailure;
  }

  //sanitize the sequence, turning ambiguity characters into the singular ambiguity character
  //character (x for nucleotide, z for amino)
  //giving the sequencePtr as the in and out arrays makes this an in-place algorithm
  fullSequenceSanitize((uint8_t*)fullSequencePtr, sanitizedSequenceCopy, fullSequenceLength, metadata->alphabetType);

  // FastaVector will always have one extra storage at the end of the string, so this shouldn't ever
  // cause a buffer overflow
  sanitizedSequenceCopy[suffixArrayLength - 1] = '$';

  //allocate the index and all internal arrays.
  struct AwFmIndex *restrict indexData = awFmIndexAlloc(metadata, suffixArrayLength);
  if(indexData == NULL){
    return AwFmAllocationFailure;
  }
  indexData->fastaVector = fastaVector;
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
    enum AwFmReturnCode returnCode = awFmWriteIndexToFile(indexData, suffixArray, (uint8_t*)fullSequencePtr, fullSequenceLength,
      fastaSrc, allowFileOverwrite);
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
      //if the suffix array isn't supposed to be kept in memory, there's no need to keep it, so it gets freed here.
      indexData->inMemorySuffixArray = NULL;
      free(suffixArray);
    }



//TODO: populate the member data of the AwFmIndex struct with the FastaVector header and metadata.

    //set the index as an out argument.
    *index = indexData;

    return returnCode;
}

void setBwtAndPrefixSums(struct AwFmIndex *restrict const index, const size_t bwtLength,
  const uint8_t *restrict const sequence, const uint64_t *restrict const suffixArray){
  if(index->metadata.alphabetType == AwFmAlphabetNucleotide){
    uint64_t baseOccurrences[8] = {0};
    //baseOccurrences is length 8 because that's how long the signpost baseOccurrences in
    //each window need to be to keep alignment to 32B AVX2 boundries.

    for(uint64_t suffixArrayPosition = 0; suffixArrayPosition < bwtLength; suffixArrayPosition++){
      const size_t  blockIndex      = suffixArrayPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t positionInBlock = suffixArrayPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t byteInVector    = positionInBlock / 8;
      const uint8_t bitInVectorByte = positionInBlock % 8;
      struct AwFmNucleotideBlock *nucleotideBlockPtr = &index->bwtBlockList.asNucleotide[blockIndex];
      uint8_t *restrict const letterBitVectorBytes = (uint8_t*)nucleotideBlockPtr->letterBitVectors;

      if(__builtin_expect(positionInBlock == 0, 0)){
        //when we start a new block, copy over the base occurrences, and initialize the bit vectors
        //while we only use 5 elements, copy over all 8 (to preserve padding and so valgrind
        //doesn't complain about invalid writes)
        memcpy(nucleotideBlockPtr->baseOccurrences, baseOccurrences, 8 * sizeof(uint64_t));
        memset(nucleotideBlockPtr->letterBitVectors, 0, sizeof(AwFmSimdVec256) * AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW);
      }

      uint64_t sequencePositionInSuffixArray = suffixArray[suffixArrayPosition];
      uint8_t letterIndex;
      if(__builtin_expect(sequencePositionInSuffixArray == 0, 0)){
        //set to the sentinel value (index 5) if we're looking at the letter before the first character.
        letterIndex = 5;
      }
      else{
        uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
        letterIndex = awFmAsciiNucleotideToLetterIndex(sequence[positionInBwt]);
      }

      uint8_t letterAsCompressedVector = awFmNucleotideLetterIndexToCompressedVector(letterIndex);
      baseOccurrences[letterIndex]++;
      letterBitVectorBytes[byteInVector]      |= ((letterAsCompressedVector >> 0) & 0x1) << bitInVectorByte;
      letterBitVectorBytes[byteInVector + 32] |= ((letterAsCompressedVector >> 1) & 0x1) << bitInVectorByte;
      letterBitVectorBytes[byteInVector + 64] |= ((letterAsCompressedVector >> 2) & 0x1) << bitInVectorByte;
    }

    //set the prefix sums
    index->prefixSums[0] = 1;   //1 is for the sentinel
    baseOccurrences[0]++;  //add the sentinel to the count of a's
    for(uint8_t i = 1; i < AW_FM_NUCLEOTIDE_CARDINALITY + 2; i++){
      index->prefixSums[i] = baseOccurrences[i-1];
      baseOccurrences[i]  += baseOccurrences[i-1];
    }

  }
  else{
    uint64_t baseOccurrences[24] = {0};

    for(uint64_t suffixArrayPosition = 0; suffixArrayPosition < bwtLength; suffixArrayPosition++){
      const size_t  blockIndex      = suffixArrayPosition / AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t positionInBlock = suffixArrayPosition % AW_FM_POSITIONS_PER_FM_BLOCK;
      const uint8_t byteInVector    = positionInBlock / 8;
      const uint8_t bitInVectorByte = positionInBlock % 8;
      struct AwFmAminoBlock *restrict const aminoBlockPointer = &index->bwtBlockList.asAmino[blockIndex];
      uint8_t *restrict const letterBitVectorBytes = (uint8_t*)aminoBlockPointer->letterBitVectors;

      if(__builtin_expect(positionInBlock == 0, 0)){
        //when we start a new block, copy over the base occurrences, and initialize the bit vectors
        //while we only use 21 elements, copy over all 24 (to preserve padding and so valgrind
        //doesn't complain about invalid writes)
        memcpy(aminoBlockPointer->baseOccurrences, baseOccurrences, 24 * sizeof(uint64_t));
        memset(aminoBlockPointer->letterBitVectors, 0, sizeof(AwFmSimdVec256) * AW_FM_AMINO_VECTORS_PER_WINDOW);
      }

      uint64_t sequencePositionInSuffixArray = suffixArray[suffixArrayPosition];
      uint8_t letterIndex;
      if(__builtin_expect(sequencePositionInSuffixArray != 0, 1)){
        uint64_t positionInBwt = sequencePositionInSuffixArray - 1;
        letterIndex = awFmAsciiAminoAcidToLetterIndex(sequence[positionInBwt]);
      }
      else{
        //set to the sentinel value (index 21) if we're looking at the letter before the first character.
        letterIndex = 21;
      }
        uint8_t letterAsVectorFormat = awFmAminoAcidLetterIndexToCompressedVector(letterIndex);
        baseOccurrences[letterIndex]++;
        letterBitVectorBytes[byteInVector]       |= (((letterAsVectorFormat >> 0) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 32]  |= (((letterAsVectorFormat >> 1) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 64]  |= (((letterAsVectorFormat >> 2) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 96]  |= (((letterAsVectorFormat >> 3) & 0x1) << bitInVectorByte);
        letterBitVectorBytes[byteInVector + 128] |= (((letterAsVectorFormat >> 4) & 0x1) << bitInVectorByte);

    }
    //set the prefix sums
    index->prefixSums[0] = 1;
    baseOccurrences[0]++;   //add the sentinel into the count of a's
    for(uint8_t i = 1; i < AW_FM_AMINO_CARDINALITY + 2; i++){
      index->prefixSums[i] = baseOccurrences[i-1];
      baseOccurrences[i]  += baseOccurrences[i-1];
    }
  }
}

void populateKmerSeedTable( struct AwFmIndex *restrict const index){
  const uint8_t alphabetCardinality = awFmGetAlphabetCardinality(index->metadata.alphabetType);
  for(uint8_t i = 0; i < alphabetCardinality; i++){
    struct AwFmSearchRange range = {
      .startPtr=  index->prefixSums[i],
      .endPtr=    index->prefixSums[i+1] -1
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

inline void fullSequenceSanitize(const uint8_t *const sequence, uint8_t *const sanitizedSequenceCopy,
  const size_t sequenceLength, const enum AwFmAlphabetType alphabetType){

  if(alphabetType == AwFmAlphabetNucleotide){
    for(size_t i = 0; i < sequenceLength; i++){
      sanitizedSequenceCopy[i] = awFmAsciiNucleotideLetterSanitize(sequence[i]);
    }
  }
  else{
    for(size_t i = 0; i < sequenceLength; i++){
      sanitizedSequenceCopy[i] = awFmAsciiAminoLetterSanitize(sequence[i]);
    }
  }
}
