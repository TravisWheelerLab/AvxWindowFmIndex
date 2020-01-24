#include "AwFmKmerTable.h"
#include "AwFmLetter.h"
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>


struct AwFmSearchRange *awFmNucleotideSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){
  //find the index of the last prefix letter, either 0 or the index of the last letter in the kmer seed.
  const int_fast16_t finalSuffixLetterPosition = kmerLength - index->metadata.kmerLengthInSeedTable;

  size_t tableIndex = 0;

  for(int_fast16_t i = kmerLength - 1; i >= finalSuffixLetterPosition; i--){
    tableIndex = (tableIndex * AW_FM_NUCLEOTIDE_CARDINALITY) + awFmAsciiNucleotideToLetterIndex(kmer[i]);
  }

  struct AwFmSearchRange *restrict const range = &index->kmerSeedTable[tableIndex];
  //prefetch the range
  _mm_prefetch(range, _MM_HINT_NTA);
  return range;
}


struct AwFmSearchRange *awFmAminoSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  //find the index of the last prefix letter, either 0 or the index of the last letter in the kmer seed.
  const int_fast16_t finalSuffixLetterPosition = kmerLength - index->metadata.kmerLengthInSeedTable;

  size_t tableIndex = 0;
  for(int_fast16_t i = kmerLength - 1; i >= finalSuffixLetterPosition; i--){
    tableIndex = (tableIndex * AW_FM_AMINO_CARDINALITY) + awFmAsciiAminoAcidToLetterIndex(kmer[i]);
  }

  struct AwFmSearchRange *restrict const range = &index->kmerSeedTable[tableIndex];
  //prefetch the range
  _mm_prefetch(range, _MM_HINT_NTA);
  return range;
}
