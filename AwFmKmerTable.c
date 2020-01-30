#include "AwFmKmerTable.h"
#include "AwFmLetter.h"
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>


struct AwFmSearchRange awFmNucleotideSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  const int_fast16_t finalSuffixLetterPosition = kmerLength - index->metadata.kmerLengthInSeedTable;

  size_t tableIndex = 0;

  for(int_fast16_t i = kmerLength - 1; i >= finalSuffixLetterPosition; i--){
    tableIndex = (tableIndex * AW_FM_NUCLEOTIDE_CARDINALITY) + awFmAsciiNucleotideToLetterIndex(kmer[i]);
  }
  // struct AwFmSearchRange range = {.startPtr= index->kmerSeedTable[tableIndex],
  //                                 .endPtr= index->kmerSeedTable[tableIndex+1]-1};
  return index->kmerSeedTable[tableIndex];
}


struct AwFmSearchRange awFmAminoSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  const int_fast16_t finalSuffixLetterPosition = kmerLength - index->metadata.kmerLengthInSeedTable;
  size_t tableIndex = 0;
  for(int_fast16_t i = kmerLength - 1; i >= finalSuffixLetterPosition; i--){
    tableIndex = (tableIndex * AW_FM_AMINO_CARDINALITY) + awFmAsciiAminoAcidToLetterIndex(kmer[i]);
  }

  return index->kmerSeedTable[tableIndex];
}
