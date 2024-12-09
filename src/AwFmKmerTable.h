#ifndef AW_FM_KMER_TABLE_H
#define AW_FM_KMER_TABLE_H

#include <stdint.h>
#include "AwFmIndexStruct.h"

/*
 * Function:  awFmNucleotideKmerSeedRangeFromTable
 * --------------------
 * Given an ascii nucleotide kmer, queries the kmerSeedTable for the partially
 * completed range in the backward suffix array, and returns a copy of that
 * range. This function looks up the range for the suffix of length equal to the
 * kmers stored in the seed table. As such, using this to query for kmers that
 * are shorter than those in the table is undefined behavior.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: ascii nucleotide character string to search for in the
 * kmerSeedTable. kmerLength: length, in characters of the kmer. This value is
 * used to determine where the suffix occurrs in the kmer string.
 *
 *  Returns:
 *    Copy of the AwFmSearchRange containing the startPtr and endPtr for the
 * kmer seed.
 */
struct AwFmSearchRange awFmNucleotideKmerSeedRangeFromTable(
    const struct AwFmIndex *_RESTRICT_ const index,
    const char *_RESTRICT_ const kmer, const size_t kmerLength);

/*
 * Function:  awFmAminoKmerSeedRangeFromTable
 * --------------------
 * Given an ascii amino acid kmer, queries the kmerSeedTable for the partially
 * completed range in the backward suffix array, and returns a copy of that
 * range. This function looks up the range for the suffix of length equal to the
 * kmers stored in the seed table. As such, using this to query for kmers that
 * are shorter than those in the table is undefined behavior.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: ascii amino acid character string to search for in the
 * kmerSeedTable. kmerLength: length, in characters of the kmer. This value is
 * used to determine where the suffix occurrs in the kmer string.
 *
 *  Returns:
 *    Copy of the AwFmSearchRange containing the startPtr and endPtr for the
 * kmer seed.
 */
struct AwFmSearchRange
awFmAminoKmerSeedRangeFromTable(const struct AwFmIndex *_RESTRICT_ const index,
                                const char *_RESTRICT_ const kmer,
                                const size_t kmerLength);

/*
 * Function:  awFmQueryCanUseKmerTable
 * --------------------
 * Determines if a given kmer can use the kmer seed table. A kmer is ineligible
 * for using the kmer seed table if it is too short or contains any ambiguity
 * characters in the suffix characters that would be used to query the table.
 *  Inputs:
 *    	index: AwFmIndex that contains the table to be used
 * 		kmer: pointer to the start of the kmer
 * 		kmerLength: length of the kmer in question.
 *
 *  Returns:
 *    true if the kmer is eligible for using the kmer seed table, or false
 * otherwise.
 */
bool awFmQueryCanUseKmerTable(const struct AwFmIndex *_RESTRICT_ const index,
                              const char *_RESTRICT_ const kmer,
                              const size_t kmerLength);

#endif /* end of include guard: AW_FM_KMER_TABLE_H */
