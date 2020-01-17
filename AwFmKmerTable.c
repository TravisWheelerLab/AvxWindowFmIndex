#include "AwFmKmerTable.h"
#include "AwFmLetter.h"
#include <stdlib.h>
#include <string.h>


struct AwFmSearchRange awFmNucleotideSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){
  //find the index of the last prefix letter, either 0 or the index of the last letter in the kmer seed.
  const int_fast16_t finalSuffixLetterPosition = kmerLength >= index->metadata.kmerLengthInSeedTable?
    kmerLength - index->metadata.kmerLengthInSeedTable: 0;

  size_t tableIndex = 0;

  for(int_fast16_t i = kmerLength - 1; i >= finalSuffixLetterPosition; i--){
    tableIndex = (tableIndex * AW_FM_NUCLEOTIDE_CARDINALITY) + awFmAsciiNucleotideToLetterIndex(kmer[i]);
  }

  if(kmerLength <= index->metadata.kmerLengthInSeedTable){
    return index->kmerSeedTable[tableIndex];
  }
  else{
    const uint64_t lowerRange = (index->kmerSeedTable[tableIndex]).startPtr;
    const uint_fast16_t lettersLeftInSeed = index->metadata.kmerLengthInSeedTable - kmerLength;

    //find the index of the kmer seed in the end of the range of valid kmers for the given one.
    //this is akin to adding 't's to the beginning.
    //we could do this with tableIndex = tableIndex * 4 + 3 in a for loop, but this should be equivalent
    //and not require a for loop.
    tableIndex = (tableIndex << (2* lettersLeftInSeed)) - 1;

    const uint64_t upperRange = (index->kmerSeedTable[tableIndex]).endPtr;
    return (struct AwFmSearchRange){lowerRange, upperRange};
  }

}


struct AwFmSearchRange awFmAminoSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength){

  //find the index of the last prefix letter, either 0 or the index of the last letter in the kmer seed.
  const int_fast16_t finalSuffixLetterPosition = kmerLength >= index->metadata.kmerLengthInSeedTable?
    kmerLength - index->metadata.kmerLengthInSeedTable: 0;

  size_t tableIndex = 0;
  for(int_fast16_t i = kmerLength - 1; i >= finalSuffixLetterPosition; i--){
    tableIndex = (tableIndex * AW_FM_AMINO_CARDINALITY) + awFmAsciiAminoAcidToLetterIndex(kmer[i]);
  }

  if(kmerLength <= index->metadata.kmerLengthInSeedTable){
    return index->kmerSeedTable[tableIndex];
  }
  else{
    const uint64_t lowerRange = (index->kmerSeedTable[tableIndex]).startPtr;
    const uint_fast16_t lettersLeftInSeed = index->metadata.kmerLengthInSeedTable - kmerLength;

    //loop through the last letters, conceptually adding 'y' aminos to the beginning
    //until we've reached the length of seeds stored in the kmer table.
    for(uint8_t i = 0; i < lettersLeftInSeed; i++){
      tableIndex = tableIndex * AW_FM_AMINO_CARDINALITY;
    }
    tableIndex--;

    const uint64_t upperRange = (index->kmerSeedTable[tableIndex]).endPtr;
    return (struct AwFmSearchRange){lowerRange, upperRange};
  }
}
