#ifndef AW_FM_KMER_TABLE_H
#define AW_FM_KMER_TABLE_H

#include "AwFmIndex.h"
#include <stdint.h>

//TODO: change kmer table to return a pointer (with prefetch!), not a value.


/*
 * Function:  awFmNucleotideKmerSeedRangeFromTable
 * --------------------
 * Given an ascii nucleotide kmer, queries the kmerSeedTable for the partially completed range in the backward suffix array,
 *  and returns a pointer to the range of that kmer. This function looks up the range for the
 *  suffix of length equal to the kmers stored in the seed table. As such, using this to query for kmers
 *  that are shorter than those in the table is undefined behavior.
 *
 *  While this function doesn't return the actual value, it does prefetch the value.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: ascii nucleotide character string to search for in the kmerSeedTable.
 *    kmerLength: length, in characters of the kmer. This value is used to determine
 *      where the suffix occurrs in the kmer string.
 *
 *  Returns:
 *    pointer to an AwFmSearchRange containing the startPtr and endPtr for the kmer seed.
 */
struct AwFmSearchRange awFmNucleotideKmerSeedRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength);


/*
 * Function:  awFmAminoKmerSeedRangeFromTable
 * --------------------
 * Given an ascii amino acid kmer, queries the kmerSeedTable for the partially completed range in the backward suffix array,
 *  and returns a pointer to the range of that kmer. This function looks up the range for the
 *  suffix of length equal to the kmers stored in the seed table. As such, using this to query for kmers
 *  that are shorter than those in the table is undefined behavior.
 *
 *  While this function doesn't return the actual value, it does prefetch the value.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: ascii amino acid character string to search for in the kmerSeedTable.
 *    kmerLength: length, in characters of the kmer. This value is used to determine
 *      where the suffix occurrs in the kmer string.
 *
 *  Returns:
 *    pointer to an AwFmSearchRange containing the startPtr and endPtr for the kmer seed.
 */
struct AwFmSearchRange awFmAminoKmerSeedRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength);


//Added for testing ONLY, not to be included in final release
//TODO: remove this function prototype
uint8_t kmerMatchesInSequenceEnding(const struct AwFmIndex *restrict const index,
  const uint64_t kmerIndexEncoding, const uint8_t kmerLength);

#endif /* end of include guard: AW_FM_KMER_TABLE_H */
