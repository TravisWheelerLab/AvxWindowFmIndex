#ifndef AW_FM_KMER_TABLE_H
#define AW_FM_KMER_TABLE_H

#include "AwFmIndex.h"
#include <stdint.h>


/*
 * Function:  awFmNucleotideSeedKmerRangeFromTable
 * --------------------
 * Given an ascii nucleotide kmer, queries the kmerSeedTable for the partially completed range in the backward suffix array.
 * This function checks the kmerLength, and if the kmer is smaller than those memoized in the seed table,
 * it extends the kmer to find the range the kmer comprises.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: ascii nucleotide character string to search for in the kmerSeedTable.
 *    kmerLength: length of the given kmer
 *
 *  Returns:
 *    AwFmSearchRange containing the startPtr and endPtr for the kmer seed.
 */
struct AwFmSearchRange awFmNucleotideSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength);

/*
 * Function:  awFmAminoSeedKmerRangeFromTable
 * --------------------
 * Given an ascii amino acid kmer, queries the kmerSeedTable for the partially completed range in the backward suffix array.
 * This function checks the kmerLength, and if the kmer is smaller than those memoized in the seed table,
 * it extends the kmer to find the range the kmer comprises.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: ascii amino acid character string to search for in the kmerSeedTable.
 *    kmerLength: length of the given kmer
 *
 *  Returns:
 *    AwFmSearchRange containing the startPtr and endPtr for the kmer seed.
 */
struct AwFmSearchRange awFmAminoSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength);

#endif /* end of include guard: AW_FM_KMER_TABLE_H */
