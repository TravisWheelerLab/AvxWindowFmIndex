#ifndef AW_FM_KMER_TABLE_H
#define AW_FM_KMER_TABLE_H

#include <stdint.h>
#include "AwFmIndex.h"


/*
 * Function:  awFmSeedKmerRangeFromTable
 * --------------------
 * Given a kmer, queries the kmerSeedTable  for the partially completed range in the backward suffix array.
 * This function checks the kmerLength, and if the kmer is smaller than those memoized in the seed table,
 * it extends the kmer to find the range the kmer comprises.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: character string to search for in the kmerSeedTable.
 *    kmerLength: length of the given kmer
 *
 *  Returns:
 *    struct containing the startPtr and endPtr for the kmer.
 */
struct AwFmBackwardRange awFmSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer, const uint8_t kmerLength);


/*
 * Function:  awFmSeedKmerRangeFromTableExactLength
 * --------------------
 * Given a kmer, queries the kmerSeedTable for the partially completed range in the backward suffix array.
 * The length of the kmer's stored in the seed table can be found in the index metadata's kmerLengthInSeedTable attribute.
 * The given kmer must be at least this long.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: character string to search for in the kmerSeedTable.
 *
 *  Returns:
 *    struct containing the backwardsRange struct in the kmerSeedTable.
 */
struct AwFmBackwardRange awFmSeedKmerRangeFromTableExactLength(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer);

#endif /* end of include guard: AW_FM_KMER_TABLE_H */
