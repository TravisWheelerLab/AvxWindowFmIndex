#include "AwFmKmerTable.h"
#include "AwFmLetter.h"

 struct AwFmBackwardRange *awFmSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
   const char *restrict const kmer){

   const uint8_t kmerSeedLength = index->metadata.kmerLengthInSeedTable;
   size_t tableIndex = 0;

   for(int_fast16_t i = kmerSeedLength - 1; i >= 0; i--){
     if(index->metadata.alphabetType == AwFmAlphabetNucleotide)
       tableIndex = (tableIndex * 4) + awFmNucleotideLetterIndexToAscii(kmer[i]);
     else
       tableIndex = (tableIndex * 20) + awFmAminoAcidLetterIndexToAscii(kmer[i]);
   }

   return &(index->kmerSeedTable[tableIndex]);
 }
