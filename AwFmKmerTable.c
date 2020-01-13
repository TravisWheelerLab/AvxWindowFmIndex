#include "AwFmKmerTable.h"
#include "AwFmLetter.h"
#include <stdlib.h>
#include <string.h>


struct AwFmSearchRange awFmSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  if(kmerLength < index->metadata.kmerLengthInSeedTable){
    char extendedKmerBuffer[index->metadata.kmerLengthInSeedTable];
    memcpy(extendedKmerBuffer, kmer, kmerLength);
    memset(extendedKmerBuffer + kmerLength, 'a', index->metadata.kmerLengthInSeedTable - kmerLength);
    const uint64_t lowerRange = awFmSeedKmerRangeFromTableExactLength(index, extendedKmerBuffer).startPtr;

    const char lastChar = index->metadata.alphabetType == AwFmAlphabetNucleotide? 't': 'y';
    memset(extendedKmerBuffer + kmerLength, lastChar, index->metadata.kmerLengthInSeedTable - kmerLength);
    const uint64_t upperRange = awFmSeedKmerRangeFromTableExactLength(index, extendedKmerBuffer).startPtr;

    return (struct AwFmSearchRange){lowerRange, upperRange};
  }
  else{
    return awFmSeedKmerRangeFromTableExactLength(index, kmer);
  }
}


struct AwFmSearchRange awFmSeedKmerRangeFromTableExactLength(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer){

  const uint8_t kmerSeedLength = index->metadata.kmerLengthInSeedTable;
  size_t tableIndex = 0;

  for(int_fast16_t i = kmerSeedLength - 1; i >= 0; i--){
    if(index->metadata.alphabetType == AwFmAlphabetNucleotide)
      tableIndex = (tableIndex * AW_FM_NUCLEOTIDE_CARDINALITY) + awFmNucleotideLetterIndexToAscii(kmer[i]);
    else
      tableIndex = (tableIndex * AW_FM_AMINO_CARDINALITY) + awFmAminoAcidLetterIndexToAscii(kmer[i]);
  }

  return (index->kmerSeedTable[tableIndex]);
}
