#include "AwFmIndex.h"
#include "AwFmCreate.h"
#include "AwFmLetter.h"
#include <stdlib.h>
#include <string.h>
#include <divsufsort64.h>


//TODO: rewrite and check AwFmCreateBlockList



struct AwFmIndex *awFmCreateIndex(const struct AwFmIndexMetadata *restrict const metadata,
  const uint8_t *const sequence, const size_t sequenceLength,
  uint64_t *restrict backwardSuffixArray, uint64_t *forwardSuffixArray){
  //allocate the index and all internal arrays.
  struct AwFmIndex *restrict index = awFmIndexAlloc(metadata, sequenceLength);
  if(index == NULL){
    return NULL;
  }

  //set the bwtLength
  index->bwtLength = sequenceLength + 1;

  //create the backwards suffix array. this function call also populates the prefix sums
  enum AwFmReturnCode returnCode = awFmCreateFullSuffixArray(sequence,
    sequenceLength, &backwardSuffixArray);
  if(!awFmReturnCodeSuccess(returnCode)){
    awFmDeleteIndex(index);
    return NULL;
  }

  //create the Bwt corresponding to the created suffixArray
  createBwt(index->backwardBwtBlockList, metadata, sequence, backwardSuffixArray);


  //if the index is bidirectional, create the forward suffix array and bwt
  if(metadata->bwtType == AwFmBwtTypeBiDirectional){

    //create the reversed sequence
    uint8_t *restrict const reverseSequence = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, sequenceLength * sizeof(uint8_t));
    if(reverseSequence == NULL){
      free(backwardSuffixArray);
      awFmDeleteIndex(index);
      return NULL;
    }
    for(size_t i = 0; i < sequenceLength; i++){
      reverseSequence[sequenceLength - i - 1] = sequence[i];
    }

    //create the forward suffix array
    returnCode = awFmCreateFullSuffixArray(reverseSequence,
      sequenceLength, &forwardSuffixArray);
    if(!awFmReturnCodeSuccess(returnCode)){
      free(reverseSequence);
      free(backwardSuffixArray);
      awFmDeleteIndex(index);
      return NULL;
    }

    //TODO: check to confirm this uses reverseSequence here, not sequence
    createBwt(index->forwardBwtBlockList, metadata, reverseSequence, forwardSuffixArray, &index->forwardSentinelCharacterPosition);
    free(reverseSequence);
  }

  createKmerSeedTable(index);

  return index;
}



void createBwt(const union AwFmBwtBlockList blockList, const struct AwFmIndexMetadata *restrict const metadata,
  const uint8_t *sequence, const uint64_t *suffixArray, const size_t suffixArrayLength){

  for(size_t positionInBwt = 0; positionInBwt < suffixArrayLength; positionInBwt++){
    const size_t blockIndex = positionInBwt % POSITIONS_PER_FM_BLOCK;
    const uint8_t letterForPosition = sequence[suffixArray[positionInBwt]];

    const uint8_t letterIndex = awFmLetterToLetterIndex(letterForPosition, metadata->alphabetType);

    uint8_t vectorEncodedLetter;
    if(metadata->alphabetType == AwFmAlphabetNucleotide){
      blockList.asNucleotide[blockIndex].baseOccurrences[letterIndex]++;
      vectorEncodedLetter = letterIndex;
    }
    else{
      vectorEncodedLetter = awFmAminoAcidAsciiLetterToCompressedVectorFormat(letterForPosition);
    }




  }


  size_t blocksInBwt = suffixArrayLength % POSITIONS_PER_FM_BLOCK;
  for(size_t blockIndex = 0; blockIndex < blocksInBwt; blockIndex++){
    for(uint16_t positionInBlock = 0; positionInBlock < POSITIONS_PER_FM_BLOCK; positionInBlock++){

      uint8_t letterForPosition = sequence[suffixArray[positionInBlock]]
      const uint8_t vectorFormatLetter =
    }
  }

}



uint64_t *createSuffixArray(const uint8_t *restrict const sequence,
  const uint64_t sequenceLength, uint64_t *restrict const prefixSums){


}




//to sort...
// at each position, grab 8 char windows (64 bit)

//16 windows at a time, use avx2 shuffle?

enum AwFmReturnCode awFmCreateReverseBwtBlockList(const char *restrict const originalString,
  const size_t textLength, struct AwFmIndex *restrict const index){
  uint8_t *restrict const reverseText = malloc(textLength *sizeof(char));
  if(reverseText == NULL){
    return AwFmAllocationFailure;
  }

  for(size_t i = 0; i < textLength; i++){
    const size_t reverseTextPosition = textLength - i;
    reverseText[reverseTextPosition] = originalString[i];
  }

  const size_t suffixArrayLength = textLength + 1;
  int64_t *restrict const reverseTextSuffixArray = malloc(suffixArrayLength * sizeof(int64_t));
  if(reverseTextSuffixArray == NULL){
    free(reverseText);
    return AwFmAllocationFailure;
  }

  uint64_t divSufSortReturnCode = divsufsort64(reverseText, reverseTextSuffixArray, (int64_t)textLength);
  if(divSufSortReturnCode < 0){
    free(reverseText);
    free(reverseTextSuffixArray);
    return AwFmSuffixArrayCreationFailure;
  }


  //TODO: use the now-created reverseTextSuffixArray to make the reverse bwt blockList
  //and assign it to the ptr in the index.


  return AwFmSuccess;
}

/*Creates the fully-sampled suffix array from the given database sequence.*/
// enum AwFmReturnCode awFmCreateFullSuffixArray(const uint8_t *databaseSequence,
//   const uint64_t databaseSequenceLength, uint64_t **fullSuffixArrayOut);

/*Creates the AwFmIndex blockList from the given database sequence.*/
// enum AwFmReturnCode awFmCreateBlockList(struct AwFmIndex *restrict const index,
//   const uint8_t *restrict const databaseSequence, const uint64_t databaseSequenceLength,
//   uint64_t *totalOccurances);

/*sets the rankPrefixSums array. This function requires that the blockList and suffix array already be set..*/
// void AwFmSetRankPrefixSums(struct AwFmIndex *restrict const index, const uint64_t *restrict const totalOccurances);

// void awFmInitBlock(struct AwFmIndex *const restrict index, const uint64_t blockIndex, uint64_t *totalOccurancesSoFar);

/*
 * Function:  awFmCreateIndex
 * --------------------
 * Create an AwFmIndex struct from a given database sequence. This function will deallocate
 *  all data it allocated if any part fails.
 *  Note: databaseSequence is not copied, instead, the AwFmIndex copies the pointer to the string.
 *  As such, the caller is responsible for the lifetime of this string, and it should not be
 *  deallocated until the AwFmIndex is written to file.
 *
 *  Inputs:
 *    indexPtr:                     Pointer to the AwFmIndex ptr to be allocated.
 *    databaseSequence:             ASCII-encoded string representing the
 *      amino-acid sequence to build the database from.
 *      NOTE: this is passes by reference (pointer copy). The caller is expected to maintain
 *      the lifetime of the databaseSequence string, and assuming that the AwFmIndex struct
 *      is going to be written to a file, the caller is expected to not deallocate this string until
 *      the AwFmIndex is written to file.
 *    databaseSequenceLength:       length of the databaseSequence string.
 *    suffixArrayCompressionRatio:  how often to sample the suffix array. Larger values
 *      reduce the size of the Fm-index data file, but increase the time it takes to reconstruct
 *      the database sequence positions when hits are found.
 *
 *  Returns:
 *    AwFmReturnCode showing the result of this action. Possible returns are:
 *      AwFmSuccess on success,
 *      AwFmNullPtrError if the indexPtr or databaseSequence were null.
 *      AwFmAllocationFailure on failure to allocate the AwFmIndex, suffix array copy,
 *        or block list.
 *      AwFmSuffixArrayCreationFailure on failure in building the suffix array (this error
 *        is caused by internal failure in the libdivsufsort library).
*/
// enum AwFmReturnCode awFmCreateIndex(const uint8_t *restrict const databaseSequence,
//   const size_t databaseSequenceLength, const struct AwFmIndexMetadata *const restrict indexMetadata,
//   struct AwFmIndex *restrict *restrict const indexPtr, uint64_t *restrict *restrict const fullSuffixArray){
//
//   if(indexPtr == NULL || databaseSequence == NULL){
//     return AwFmNullPtrError;
//   }
//
//   struct AwFmIndex *index = awFmAlignedAllocAwFmIndex();
//   if(index == NULL){
//     return AwFmAllocationFailure;
//   }
//
//   //copy over the given metadata
//   memcpy(&index->metadata, indexMetadata, sizeof(struct AwFmIndexMetadata));
//
//   //create the full suffix array
//   enum AwFmReturnCode suffixArrayCreationReturnCode = awFmCreateFullSuffixArray(databaseSequence,
//     databaseSequenceLength, fullSuffixArray);
//   if(__builtin_expect(suffixArrayCreationReturnCode < 0, 0)){
//     awFmDeallocateFmIndex(index);
//     return suffixArrayCreationReturnCode;
//   }
//
//   //the occurances values will accumulate into this array.
//   // the +1 leaves one slot at the end for unrecognized letters to accumulate into,
//   // and will be treated as a garbage data slot so that illegal characters don't crash the application.
//   uint64_t totalOccurances[AMINO_CARDINALITY + 1];
//   //create the block list
//   enum AwFmReturnCode blockListCreationReturnCode = awFmCreateBlockList(index,
//     databaseSequence, databaseSequenceLength, totalOccurances);
//
//   if(__builtin_expect(blockListCreationReturnCode < 0, 0)){
//     awFmDeallocateFmIndex(index);
//     return blockListCreationReturnCode;
//   }
//
//   awFmSetRankPrefixSums(index, totalOccurances);
//
//   //set the out pointer for the index.
//   *indexPtr = index;
//
//   return AwFmSuccess;
// }


void awFmSetRankPrefixSums(struct AwFmIndex *restrict const index, const uint64_t *restrict const totalOccurances){
  //set the rank prefix sums
  //accumulator starts at 1 because of the null terminator $ being before all other letters
  uint64_t rankPrefixAccumulator = 1;

  index->rankPrefixSums[0] = rankPrefixAccumulator;
  for(uint_fast8_t i = 0; i < AMINO_CARDINALITY; i++){
    rankPrefixAccumulator += totalOccurances[i];
    index->rankPrefixSums[i + 1] = rankPrefixAccumulator;
  }
}


//if failure, the caller may need to deallocate index.
/*
 * Function:  awFmCreateBlockList
 * --------------------
 * Allocates and initializes the BlockList with the suffix array created by
 *  the given database sequence. If this function returns a failure value (any
 *  value < 0), the given index should be deallocated.
 *
 *  Inputs:
 *    index              AwFmIndex struct that's block list will be set.
 *    databaseSequence:             ASCII-encoded string representing the
 *      amino-acid sequence to build the blockList from.
 *    databaseSequenceLength:       length of the databaseSequence string to be copied
 *      into the AwFmIndex struct.
 *    totalOccurances:   Array of at least length 21 (AMINO_CARDINALITY + 1) whose
 *      values will be assigned the full occurances for each amino over the entire database.
 *      the array has one final element at the end so that the null terminator has
 *      somewhere to accumulate to, which will later be thrown away.
 *
 *  Returns:
 *    AwFmReturnCode showing the result of this action. Possible returns are:
 *      AwFmSuccess on success,
 *      AwFmNullPtrError if the indexPtr, databaseSequence, or totalOccurances were null.
 *      AwFmAllocationFailure on failure to allocate the block list.
*/
enum AwFmReturnCode awFmCreateBlockList(struct AwFmIndex *const restrict index,
  const uint8_t *restrict const databaseSequence, const uint64_t databaseSequenceLength,
  uint64_t *totalOccurances){

  if(index == NULL || databaseSequence == NULL || totalOccurances == NULL){
    return AwFmNullPtrError;
  }

  const uint64_t suffixArrayLength = databaseSequenceLength + 1;

  //ceiling function to compute the number of blocks needed
  //I would normally write this in the form (x+(y-1)/y), but written like this
  // avoids the potential of overflow.
  index->numBlocks = awFmNumBlocksFromSuffixArrayLength(suffixArrayLength);
  index->blockList = awFmAlignedAllocBlockList(index->numBlocks);
  if(index->blockList == NULL){
    return AwFmAllocationFailure;
  }

  //keep an accumulator for the baseOccurances and current position in suffix array
  //by adding 1 to the length, anything that parses incorrectly won't crash the program,
  //but will accumulate into an unused value
  memset(totalOccurances, 0,  AMINO_CARDINALITY * sizeof(uint64_t));

  for(uint64_t blockIndex = 0; blockIndex < index->numBlocks; blockIndex++){
    awFmInitBlock(index, blockIndex, totalOccurances, suffixArrayLength);
  }

  return AwFmSuccess;
}


enum AwFmReturnCode AwFmCreateBlockList(const uint8_t *restrict const text, const size_t suffixArrayLength,
  const uint64_t *restrict const suffixArray, const enum AwFmAlphabetType alphabetType,
  union AwFmBwtBlockList *restrict const blockList, uint64_t *totalOccurances){

  size_t numBlocks = awFmNumBlocksFromSuffixArrayLength(suffixArrayLength);

  if(alphabetType == AwFmAlphabetNucleotide){
    struct AwFmNucleotideBlock *nucleotideBlockList = calloc(numBlocks, sizeof(struct AwFmNucleotideBlock));
    if(nucleotideBlockList == NULL){
      return AwFmAllocationFailure;
    }

    uint64_t *totalOccurancesSoFar = malloc(AwFmAlphabetNucleotideCardinality * sizeof(uint64_t));
    if(totalOccurancesSoFar == NULL){
      free(nucleotideBlockList);
      return AwFmAllocationFailure;
    }

    for(size_t blockIndex = 0; blockIndex < numBlocks; blockIndex++){
      AwFmInitNucleotideBlock(&nucleotideBlockList[blockIndex],
      totalOccurancesSoFar, suffixArrayLength);
    }

    totalOccurances         = totalOccurancesSoFar;
    blockList->asNucleotide = nucleotideBlockList;
    return AwFmSuccess;

  }
  else if(alphabetType == AwFmAlphabetAminoAcid){
    struct AwFmAminoBlock *aminoBlockList = malloc(numBlocks * sizeof(struct AwFmAminoBlock));
    if(aminoBlockList == NULL){
      return AwFmAllocationFailure;
    }

    uint64_t *totalOccurancesSoFar = calloc(AwFmAlphabetAminoAcidCardinality, sizeof(uint64_t));
    if(totalOccurancesSoFar == NULL){
      free(aminoBlockList);
      return AwFmAllocationFailure;
    }

    for(size_t blockIndex = 0; blockIndex < numBlocks; blockIndex++){
      AwFmInitAminoBlock(&aminoBlockList[blockIndex],
      totalOccurancesSoFar, suffixArrayLength);
    }

    totalOccurances     = totalOccurancesSoFar;
    blockList->asAmino  = aminoBlockList;
    return AwFmSuccess;

  }
  else{
    //return an error code if the alphabet type isn't of one of the supported types.
    //currently, only nucleotide and amino acid data are supported.
    return AwFmUnsupportedAlphabetError;
  }

}



void awFmInitBlock(union AwFmBwtBlockList blockList, const enum AwFmAlphabetType alphabet, const uint64_t blockIndex,
  const uint64_t numBlocksInList, const uint64_t *restrict const suffixArray, const size_t suffixArrayLength,
  const uint8_t *sequence, uint64_t *totalOccurancesSoFar){


    uint8_t *restrict letterBitVectorsAsBytePtr;
    //copy the base occurances over, and set the bit vectors to all zeros.
  if(alphabet ==AwFmAlphabetNucleotide){
    memcpy(blockList.asNucleotide[blockIndex].baseOccurrences, totalOccurancesSoFar, 4 * sizeof(uint64_t));
    memset(blockList.asNucleotide[blockIndex].letterBitVectors, 0, sizeof(__m256i) * 2);
    letterBitVectorsAsBytePtr = (uint8_t*)blockList.asNucleotide[blockIndex].letterBitVectors;
  }
  else{
    memcpy(blockList.asAmino[blockIndex].baseOccurrences, totalOccurancesSoFar, 20 * sizeof(uint64_t));
    memset(blockList.asAmino[blockIndex].letterBitVectors, 0, sizeof(__m256i) * 5);
    letterBitVectorsAsBytePtr = (uint8_t*)blockList.asAmino[blockIndex].letterBitVectors;
  }


  bool isLastBlock = blockIndex == numBlocksInList - 1;
  //if we're in the last block, we need to stop at the actual end of the suffix array.
  uint_fast16_t lengthOfThisBlock = isLastBlock?
    (suffixArrayLength - (blockIndex * POSITIONS_PER_FM_BLOCK)):
    POSITIONS_PER_FM_BLOCK;


  for(uint_fast16_t i = 0; i < lengthOfThisBlock; i++){

    const uint8_t byteInBlock             = i / 8;
    const uint8_t bitInBlockByte          = i % 8;
    const uint64_t positionInSuffixArray  = (blockIndex * POSITIONS_PER_FM_BLOCK) + i;

    const uint64_t suffixArrayValue       = suffixArray[positionInSuffixArray];
    //the character in the BWT is the character immediately before the one represented
    // in the suffix array. But, since we can point to the first character (at index 0),
    // we need to check for 0 so we actually point to the last character, rather
    //than overflowing.


    if(__builtin_expect(suffixArrayValue != 0, 1)){
      const uint8_t letterAsAscii           = sequence[suffixArrayValue - 1];
      const uint8_t letterIndex             = awFmLetterToLetterIndex(letterAsAscii, metadata->alphabetType);
      const uint8_t letterAsVectorFormat    = metadata->alphabetType == AwFmAlphabetNucleotide?
        letterIndex:
        awFmAminoAcidAsciiLetterToCompressedVectorFormat(letterAsAscii);

      //accumulate the letter's occurances.
      totalOccurancesSoFar[letterIndex]++;

      const uint8_t numberOfBitVectors = metadata->alphabetType == AwFmAlphabetNucleotide? 2: 5;
      //set the correct bits in the letterBitVectors
      for(uint_fast8_t bitInVectorLetter = 0; bitInVectorLetter < numberOfBitVectors; bitInVectorLetter++){
        const uint8_t letterBit = (letterAsVectorFormat >> bitInVectorLetter) & 1;
        const uint8_t byteInBlockVectors = (bitInVectorLetter * BYTES_PER_AVX2_REGISTER) + byteInBlock;
        //or in the shifted bit.
        letterBitVectorsAsBytePtr[byteInBlockVectors] |= (letterBit << bitInBlockByte);
      }
    }
  }
}



/*
 * Function:  awFmCreateFullSuffixArray
 * --------------------
 * Creates a fully-sampled suffix array from the given database sequence.
 *
 *  Inputs:
 *    databaseSequence:       ASCII-encoded string representing the
 *      amino-acid sequence to build the blockList from.
 *    databaseSequenceLength: Length of the databaseSequence string to be copied
 *      into the AwFmIndex struct.
 *    fullSuffixArrayOut:     Pointer to the suffixArray ptr to be allocated and set
 *      to the computed suffix array.
 *
 *  Returns:
 *    AwFmReturnCode showing the result of this action. Possible returns are:
 *      AwFmSuccess on success,
 *      AwFmAllocationFailure on failure to allocate the suffix array.
 *      AwFmNullPtrError if the indexPtr or databaseSequence were null.
 *      AwFmSuffixArrayCreationFailure on failure in building the suffix array (this error
 *        is caused by internal failure in the libdivsufsort library).
*/
enum AwFmReturnCode awFmCreateFullSuffixArray(const uint8_t *databaseSequence,
  const uint64_t databaseSequenceLength, uint64_t *restrict *restrict const fullSuffixArrayOut){

  //suffix array is one longer than database sequence, since it needs to store the null terminating $.
  const uint64_t suffixArrayLength = databaseSequenceLength + 1;
  uint64_t *fullSuffixArray = malloc(suffixArrayLength * sizeof(uint64_t));
  if(fullSuffixArray == NULL){
    return AwFmAllocationFailure;
  }
  int64_t *suffixArrayAfterTerminator = (int64_t*)(fullSuffixArray + 1);
  uint64_t divSufSortReturnCode = divsufsort64(databaseSequence, (int64_t *)suffixArrayAfterTerminator, (int64_t)databaseSequenceLength);

  //the first position will be the location of the null terminator.
  fullSuffixArray[0] = databaseSequenceLength;

  *fullSuffixArrayOut = fullSuffixArray;

  if(divSufSortReturnCode < 0){
    free(fullSuffixArray);
    return AwFmSuffixArrayCreationFailure;
  }
  else{
    return AwFmSuccess;
  }
}
