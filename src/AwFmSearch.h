#ifndef AW_FM_INDEX_SEARCH_H
#define AW_FM_INDEX_SEARCH_H

#include <stdbool.h>
#include <stdint.h>
#include "AwFmFile.h"
#include "AwFmIndex.h"
#include "AwFmIndexStruct.h"

/*
 * Function:  awFmBacktraceBwtPosition
 * --------------------
 * Given a specified Bwt position, backsteps to find the position one before in
 *  original sequence.
 *
 *  Inputs:
 *    index: Index to backstep
 *    alphabet: alphabet of the index, either Dna, Rna, or Amino
 *    bwtPosition: Position of the character to be returned.
 *
 *  Returns:
 *    Position in the suffix array of the character in the sequence immediately
 * preceeding the one found at the given bwtPosition.
 */
size_t awFmNucleotideBacktraceBwtPosition(
    const struct AwFmIndex *_RESTRICT_ const index, const uint64_t bwtPosition);

/*
 * Function:  awFmBacktraceBwtPosition
 * --------------------
 * Given a specified Bwt position, backsteps to find the position one before in
 *  original sequence.
 *
 *  Inputs:
 *    index: Index to backstep
 *    alphabet: alphabet of the index, either Dna, Rna, or Amino
 *    bwtPosition: Position of the character to be returned.
 *
 *  Returns:
 *    Position in the suffix array of the character in the sequence immediately
 * preceeding the one found at the given bwtPosition.
 */
size_t
awFmAminoBacktraceBwtPosition(const struct AwFmIndex *_RESTRICT_ const index,
                              const uint64_t bwtPosition);

/*
 * Function:  awFmSingleKmerExists
 * --------------------
 *  Queries the FM-Index to determine if the database sequence contains any
 * instances of the given kmer
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur,
 * including creating potential segfauts. kmerLength:   Length of the kmer to be
 * queried. Undefined behavior may occur if the function is given a kmerLength
 * of 0.
 *
 *  Returns:
 *    True if the kmer exists in the database sequence, or false if it
 *      cannot be found in the database sequence.
 */
bool awFmSingleKmerExists(const struct AwFmIndex *_RESTRICT_ const index,
                          const char *_RESTRICT_ const kmer,
                          const size_t kmerLength);

/*
 * Function:  awFmNucleotideNonSeededSearch
 * --------------------
 *  Finds the range in the bwt corresponding to the entire given nucleotide
 * kmer. While this function will function on any input kmer, this function
 * explicitly exists to query for kmers that are too short to be memoized in the
 * kmerTable. If the kmer you are searching for is at least as long as those in
 * the kmerTable, do not use this function. Instead, start with the memozied
 * range for the kmer suffix and extend from there.
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur,
 * including creating potential segfauts. kmerLength:   Length of the kmer to be
 * queried. Undefined behavior may occur if the function is given a kmerLength
 * of 0. range:        Pointer to a search range where the output will be
 * written. Since gcc doesn't seem to perform NRVO correctly, this performs
 * slightly better.
 *
 */
void awFmNucleotideNonSeededSearch(
    const struct AwFmIndex *_RESTRICT_ const index,
    const char *_RESTRICT_ const kmer, const size_t kmerLength,
    struct AwFmSearchRange *range);

/*
 * Function:  awFmAminoNonSeededSearch
 * --------------------
 *  Finds the range in the bwt corresponding to the entire given amino kmer.
 * While this function will function on any input kmer, this function explicitly
 * exists to query for kmers that are too short to be memoized in the kmerTable.
 *    If the kmer you are searching for is at least as long as those in the
 * kmerTable, do not use this function. Instead, start with the memozied range
 * for the kmer suffix and extend from there.
 *
 *  Inputs:
 *    index:        Pointer to the valid AwFmIndex struct.
 *    kmer:         Pointer to the kmer character string.
 *      kmer MUST point to valid data, otherwise, undefined behavior may occur,
 * including creating potential segfauts. kmerLength:   Length of the kmer to be
 * queried. Undefined behavior may occur if the function is given a kmerLength
 * of 0. range:        Pointer to a search range where the output will be
 * written. Since gcc doesn't seem to perform NRVO correctly, this performs
 * slightly better.
 *
 */
void awFmAminoNonSeededSearch(const struct AwFmIndex *_RESTRICT_ const index,
                              const char *_RESTRICT_ const kmer,
                              const size_t kmerLength,
                              struct AwFmSearchRange *range);

#endif /* end of include guard: AW_FM_INDEX_SEARCH_H */
