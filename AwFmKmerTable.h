#include <stdint.h>
#include "AwFmIndex.h"

/*
 * Function:  awFmSeedKmerRangeFromTable
 * --------------------
 * Given a kmer, queries the kmerSeedTable  for the partially completed range in the backward suffix array.
 * The length of the kmer's stored in the seed table can be found in the index metadata's kmerLengthInSeedTable attribute.
 * The given kmer must be at least this long.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    kmer: character string to search for in the kmerSeedTable.
 *
 *  Returns:
 *    pointer to the backwardsRange struct in the kmerSeedTable.
 */
struct AwFmBackwardRange *awFmSeedKmerRangeFromTable(const struct AwFmIndex *restrict const index,
  const char *restrict const kmer);
