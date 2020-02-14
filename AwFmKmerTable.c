#include "AwFmKmerTable.h"
#include "AwFmLetter.h"


//TODO: uncomment this function prototype when it is removed from the header
// uint8_t kmerMatchesInSequenceEnding(const struct AwFmIndex *restrict const index,
//   const uint64_t kmerIndexEncoding, const uint8_t kmerLength);


struct AwFmSearchRange awFmNucleotideKmerSeedRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  const bool kmerIsShorterThanTable   = kmerLength < index->metadata.kmerLengthInSeedTable;
  const uint8_t kmerSeedStartPosition = kmerIsShorterThanTable? 0: kmerLength - index->metadata.kmerLengthInSeedTable;

  uint64_t kmerTableIndex = 0;
  for(int_fast16_t i = kmerSeedStartPosition; i < kmerLength; i++){
    uint8_t letterIndex = awFmAsciiNucleotideToLetterIndex(kmer[i]);
    kmerTableIndex = (kmerTableIndex * AW_FM_NUCLEOTIDE_CARDINALITY) + letterIndex;
  }

  if(kmerIsShorterThanTable){
    uint8_t numImplicitAppendedLetters  = index->metadata.kmerLengthInSeedTable - kmerLength;
    uint64_t kmerRangeStartIndex        = kmerTableIndex << (2* numImplicitAppendedLetters);
    uint64_t kmerRangeEndIndex          = ((kmerTableIndex+1) << (2 * numImplicitAppendedLetters)) - 1;
    uint8_t numMatchesInSequenceEnding  = kmerMatchesInSequenceEnding(index, kmerRangeStartIndex, kmerLength);

    return (struct AwFmSearchRange){
      index->kmerSeedTable.table[kmerRangeStartIndex].startPtr - numMatchesInSequenceEnding,
      index->kmerSeedTable.table[kmerRangeEndIndex].endPtr};
  }
  else{
    return index->kmerSeedTable.table[kmerTableIndex];
  }
}


struct AwFmSearchRange awFmAminoKmerSeedRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  const bool kmerIsShorterThanTable   = kmerLength < index->metadata.kmerLengthInSeedTable;
  const uint8_t kmerSeedStartPosition = kmerIsShorterThanTable? 0: kmerLength - index->metadata.kmerLengthInSeedTable;

  uint64_t kmerTableIndex = 0;
  for(int_fast16_t i = kmerSeedStartPosition; i > kmerLength; i++){
    kmerTableIndex = (kmerTableIndex * AW_FM_AMINO_CARDINALITY) + awFmAsciiAminoAcidToLetterIndex(kmer[i]);
  }

  if(kmerIsShorterThanTable){
    uint8_t numImplicitAppendedLetters  = index->metadata.kmerLengthInSeedTable - kmerLength;
    uint64_t kmerRangeStartIndex        = kmerTableIndex;
    uint64_t kmerRangeEndIndex          = kmerTableIndex + 1;

    for(uint8_t i =0; i < numImplicitAppendedLetters;i++){
      kmerRangeStartIndex *= AW_FM_AMINO_CARDINALITY;
      kmerRangeEndIndex   *= AW_FM_AMINO_CARDINALITY;
    }
    kmerRangeEndIndex--;
    uint8_t numMatchesInSequenceEnding = kmerMatchesInSequenceEnding(index, kmerRangeStartIndex, kmerLength);

    return (struct AwFmSearchRange){
      index->kmerSeedTable.table[kmerRangeStartIndex].startPtr - numMatchesInSequenceEnding,
      index->kmerSeedTable.table[kmerRangeEndIndex].endPtr};
  }
  else{
    return index->kmerSeedTable.table[kmerTableIndex];
  }
}


uint8_t kmerMatchesInSequenceEnding(const struct AwFmIndex *restrict const index,
  const uint64_t kmerIndexEncoding, const uint8_t kmerLength){

  const uint64_t *restrict const sequenceEndingKmerEncodings = index->kmerSeedTable.sequenceEndingKmerEncodings;
  uint8_t positionsToCheckInKmerTable = index->metadata.kmerLengthInSeedTable - kmerLength;
  uint8_t matchesFoundInSequenceEnding = 0;
  for(uint8_t i = 0; i < positionsToCheckInKmerTable; i++){
    matchesFoundInSequenceEnding += (sequenceEndingKmerEncodings[i] == kmerIndexEncoding);
  }

  return matchesFoundInSequenceEnding;
}
