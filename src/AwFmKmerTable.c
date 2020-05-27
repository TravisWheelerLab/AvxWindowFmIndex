#include "AwFmKmerTable.h"
#include "AwFmLetter.h"


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
