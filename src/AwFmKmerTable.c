#include "AwFmKmerTable.h"
#include "AwFmLetter.h"
#include "AwFmSearch.h"


struct AwFmSearchRange awFmNucleotideKmerSeedRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  const uint8_t kmerSeedStartPosition = kmerLength - index->metadata.kmerLengthInSeedTable;
  uint64_t kmerTableIndex = 0;
  for(int_fast16_t i = kmerSeedStartPosition; i < kmerLength; i++){
    uint8_t letterIndex = awFmAsciiNucleotideToLetterIndex(kmer[i]);
    kmerTableIndex = (kmerTableIndex * AW_FM_NUCLEOTIDE_CARDINALITY) + letterIndex;
  }

  return index->kmerSeedTable[kmerTableIndex];
}


struct AwFmSearchRange awFmAminoKmerSeedRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  const uint8_t kmerSeedStartPosition = kmerLength - index->metadata.kmerLengthInSeedTable;
  uint64_t kmerTableIndex = 0;
  for(int_fast16_t i = kmerSeedStartPosition; i < kmerLength; i++){
    uint8_t letterIndex = awFmAsciiAminoAcidToLetterIndex(kmer[i]);
    kmerTableIndex = (kmerTableIndex * AW_FM_AMINO_CARDINALITY) + letterIndex;
  }

  return index->kmerSeedTable[kmerTableIndex];
}


// definition check for an intentionally undefined variable so this code doesn't
// get implemented, get used, or throw
#ifdef AW_FM_PARTIAL_SEED_CODE_IMPLEMENTATION_PROVIDED
static inline struct AwFmSearchRange awFmNucleotidePartialKmerSeedRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  assert(false, "this function is unfinished, untested, and left only to be finished for future releases");
  struct AwFmSearchRange range = {0,0};

  if(__builtin_expect((kmer[kmerLength-1] | 0x20) == 't', 0)){
    //lookup kmer from scratch
    awFmNucleotideNonSeededSearch(index, kmer, kmerLength, &range);
  }
  else{

    const uint64_t endSequenceKmerTableIndex = 0;// index->endSequenceKmerTableIndex;
    const uint8_t kmerLengthInSeedTable = index->metadata.kmerLengthInSeedTable;
    const uint8_t nucleotideBitOffset = 2;
    uint64_t queryKmerStartTableIndex = 0;
    uint64_t queryKmerEndTableIndex = 0;

    for(uint8_t kmerLetterIndex = 0; kmerLetterIndex < kmerLength; kmerLetterIndex++){
      queryKmerStartTableIndex <<= nucleotideBitOffset;
      queryKmerStartTableIndex |= awFmAsciiNucleotideToLetterIndex(kmer[kmerLetterIndex]);
    }

    //adds 1 to the last character, this will never overflow the last character,
    //since kmers that end with a 't' must use non-seeded search.
    queryKmerEndTableIndex = queryKmerStartTableIndex + 1;
    uint64_t indexBitmask = (1 << (kmerLength * nucleotideBitOffset)) - 1;

    for(uint8_t i = kmerLength; i < kmerLengthInSeedTable; i++){
      bool baseMatchesEndingIndex = ((endSequenceKmerTableIndex & indexBitmask) ^ queryKmerStartTableIndex) == 0;
      bool endMatchesEndingIndex  = ((endSequenceKmerTableIndex & indexBitmask) ^ queryKmerEndTableIndex) == 0;

      if(baseMatchesEndingIndex){
        range.startPtr--;
      }
      if(endMatchesEndingIndex){
        range.endPtr--;
      }
      indexBitmask              = (indexBitmask << nucleotideBitOffset) | 0x3;
      queryKmerStartTableIndex  <<= nucleotideBitOffset;
      queryKmerEndTableIndex    <<= nucleotideBitOffset;
    }

    //add in the start pointers from the ranges extended with a's.
    range.startPtr  += index->kmerSeedTable[queryKmerStartTableIndex].startPtr;
    range.endPtr    += index-> kmerSeedTable[queryKmerEndTableIndex].startPtr - 1;
  }

  return range;
}


inline struct AwFmSearchRange awFmAminoPartialKmerSeedRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){
    //todo: use bitpacked 5 bit encodings for partial
  struct AwFmSearchRange range = {0,0};

  if(__builtin_expect((kmer[kmerLength-1] | 0x20) == 'y', 0)){
    //lookup kmer from scratch
    awFmAminoNonSeededSearch(index, kmer, kmerLength, &range);
  }
  else{

    const uint64_t endSequenceKmerTableIndex  = 0;//index->endSequenceKmerTableIndex;
    const uint8_t kmerLengthInSeedTable       = index->metadata.kmerLengthInSeedTable;
    uint64_t queryKmerStartTableIndex         = 0;
    uint64_t queryKmerEndTableIndex           = 0;
    uint64_t indexModFactor                   = 1;

    for(uint8_t kmerLetterIndex = 0; kmerLetterIndex < kmerLength; kmerLetterIndex++){
      indexModFactor            *= AW_FM_AMINO_CARDINALITY;
      queryKmerStartTableIndex  *= AW_FM_AMINO_CARDINALITY;
      queryKmerStartTableIndex  |= awFmAsciiAminoAcidToLetterIndex(kmer[kmerLetterIndex]);
    }

    //adds 1 to the last character, this will never overflow the last character,
    //since kmers that end with a 'y' must use non-seeded search.
    queryKmerEndTableIndex = queryKmerStartTableIndex + 1;

    for(uint8_t i = kmerLength; i < kmerLengthInSeedTable; i++){
      uint64_t endSequenceKmerTableIndexMod = endSequenceKmerTableIndex % indexModFactor;
      bool baseMatchesEndingIndex = (queryKmerStartTableIndex  == endSequenceKmerTableIndexMod);
      bool endMatchesEndingIndex  = (queryKmerEndTableIndex    == endSequenceKmerTableIndexMod);

      if(baseMatchesEndingIndex){
        range.startPtr--;
      }
      if(endMatchesEndingIndex){
        range.endPtr--;
      }
      indexModFactor            *= AW_FM_AMINO_CARDINALITY;
      queryKmerStartTableIndex  *= AW_FM_AMINO_CARDINALITY;
      queryKmerEndTableIndex    *= AW_FM_AMINO_CARDINALITY;
    }

    //add in the start pointers from the ranges extended with a's.
    range.startPtr  += index->kmerSeedTable[queryKmerStartTableIndex].startPtr;
    range.endPtr    += index->kmerSeedTable[queryKmerEndTableIndex].startPtr - 1;
  }

  return range;
}

#endif
