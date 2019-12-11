#include "AwFmIndex.h"
#include "AwFmGlobals.h"
#include "AwFmFile.h"
#include <string.h>
#include <stdlib.h>



struct AwFmIndex *awFmIndexAlloc(const struct AwFmIndexMetadata *restrict const metadata,
  const size_t sequenceLength){

  //allocate the index
  struct AwFmIndex *index = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, sizeof(struct AwFmIndex));
  if(index == NULL){
    return NULL;
  }
  //initialize all bytes in the index to 0.
  memset(index, 0, sizeof(struct AwFmIndex));

  //allocate the prefixSums
  size_t alphabetSize = metadata->alphabetType == AwFmAlphabetNucleotide? 4: 20;
  index->prefixSums = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, alphabetSize * sizeof(uint64_t));
  if(index->prefixSums == NULL){
    awFmDeallocIndex(index);
    return NULL;
  }

  //allocate the blockLists
  size_t numBlocksInBwt = (sequenceLength + 1) % POSITIONS_PER_FM_BLOCK;
  size_t sizeOfBwtBlock = metadata->alphabetType == AwFmAlphabetNucleotide?
    sizeof(struct AwFmNucleotideBlock): sizeof(struct AwFmAminoBlock);

  //alloc the backward bwt
  index->backwardBwtBlockList.asNucleotide = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, numBlocksInBwt * sizeOfBwtBlock);
  if(index->backwardBwtBlockList.asNucleotide == NULL){
    awFmDeallocIndex(index);
    return NULL;
  }

  //alloc the forward bwt (if bidirectional)
  if(metadata->bwtType == AwFmBwtTypeBiDirectional){
    index->forwardBwtBlockList.asNucleotide = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, numBlocksInBwt * sizeOfBwtBlock);
    if(index->forwardBwtBlockList.asNucleotide == NULL){
      awFmDeallocIndex(index);
      return NULL;
    }
  }
  else{
    index->forwardBwtBlockList.asNucleotide = NULL;
  }

  //compute the size of the kmer seed table (essentially a power function)
  size_t kmerSeedTableSize = 1;
  for(size_t i = 0; i < metadata->kmerLengthInSeedTable; i++){
    kmerSeedTableSize *= alphabetSize;
  }

  //allocate the kmerSeedTable
  index->kmerSeedTable = aligned_alloc(CACHE_LINE_SIZE_IN_BYTES, kmerSeedTableSize * sizeof(uint64_t));
  if(index->kmerSeedTable == NULL){
    awFmDeallocIndex(index);
    return NULL;
  }

  return index;
}


void awFmDeallocIndex(struct AwFmIndex *index){
  if(index != NULL){
    fclose(index->fileHandle);
    free(index->backwardBwtBlockList.asNucleotide);
    free(index->forwardBwtBlockList.asNucleotide);
    free(index->prefixSums);
    free(index->kmerSeedTable);
    free(index);
  }
}


uint_fast8_t awFmGetAlphabetCardinality(const enum AwFmAlphabetType alphabet){
  return (alphabet == AwFmAlphabetNucleotide)? 4: 20;
}

size_t awFmGetKmerTableLength(const struct AwFmIndexMetadata *restrict const metadata){
  const size_t multiplier = awFmGetAlphabetCardinality(metadata->alphabetType);
  size_t length = 1;
  for(size_t i = 0; i < metadata->kmerLengthInSeedTable; i++){
    length *= multiplier;
  }

  return length;
}

/*
 * Function:  awFmBwtPositionIsSampled
 * --------------------
 * Determines if the given BWT position is sampled under the given AwFmIndex's
 *  suffix array compression ratio.
 *
 *  Inputs:
 *    index:      Pointer to the BWT's AwFmIndex.
 *    position:   Position in the BWT.
 *
 *  Outputs:
 *    True if the position is sampled in the compressedSuffixArray, or false otherwise.
 */
 bool awFmBwtPositionIsSampled(const struct AwFmIndex *restrict const index, const uint64_t position){
  return position % index->metadata.suffixArrayCompressionRatio == 0;
}


 uint64_t awFmGetCompressedSuffixArrayLength(const struct AwFmIndex *restrict const index){
  return index->bwtLength / index->metadata.suffixArrayCompressionRatio;
}


/*
 * Function:  awFmSearchRangeIsValid
 * --------------------
 * Compares the start pointer and end pointer to determine if this range
 *  contains a valid range of the relevant kmer suffix.
 *
 *  Inputs:
 *    searchRange: AwFmSearchRange struct to compare.
 *
 *  Outputs:
 *    True if this range represents a range of positions for the given kmer suffix in the database,
 *      or false if the database does not contain the given kmer suffix.
 */
 bool awFmSearchRangeIsValid(const struct AwFmBackwardRange *restrict const searchRange){
  return searchRange->startPtr <= searchRange->endPtr;
}


size_t awFmNumBlocksFromBwtLength(const size_t suffixArrayLength){
  return  1 + ((suffixArrayLength -1) / POSITIONS_PER_FM_BLOCK);
}

size_t awFmNumBlocksFromSequenceLength(const size_t databaseSequenceLength){
  return awFmNumBlocksFromBwtLength(databaseSequenceLength + 1);
}
