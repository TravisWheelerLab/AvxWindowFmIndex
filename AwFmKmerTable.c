#include "AwFmKmerTable.h"
#include "AwFmLetter.h"

/*
 * Function:  awFmKmerTableLookup
 * --------------------
 * Given a kmer seed and the AwFmIndex, this function gets the backward BWT searchRange for the seed.
 * The length of the kmer seed is sequal to the kmerLengthInSeedTable member data in the index metadata.
 * As such, the kmerSeed argument must be this long.
 *  Inputs:
 *    index: AwFmIndex structure to query for the range.
 *    kmerSeed: string representing the kmer to query for.
 *
 *  Returns:
 *    Range in the backward BWT corresponding to this seed.
 *
 */
inline struct AwFmBackwardRange awFmKmerTableLookup(const struct AwFmIndex *restrict const index,
  const uint8_t *restrict const kmerSeed){
  const uint64_t sentinelCharacterPosition  = index->sentinelCharacterPosition;
  const uint8_t kmerLengthInSeedTable       = index->metadata.kmerLengthInSeedTable;

  //build the position in the table that stores the startPtr
  uint64_t kmerPositionInTable = 0;
  for( uint8_t i = 0; i < kmerLengthInSeedTable; i++){
    uint8_t letterAsCompressedEncoding = awFmLetterToLetterIndex(kmerSeed[i], index->metadata.alphabetType);
    kmerPositionInTable *= (index->metadata.alphabetType == AwFmAlphabetNucleotide)? 4: 20;

    kmerPositionInTable += letterAsCompressedEncoding;
  }

  const uint64_t positionBeforeNextKmer = index->kmerSeedTable[kmerPositionInTable + 1] - 1;

  struct AwFmBackwardRange range;
  range.startPtr  = index->kmerSeedTable[kmerPositionInTable];
  range.endPtr    = (positionBeforeNextKmer == sentinelCharacterPosition)? positionBeforeNextKmer - 1: positionBeforeNextKmer;

  return range;
}
