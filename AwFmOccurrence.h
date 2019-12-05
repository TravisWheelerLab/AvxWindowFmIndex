#ifndef AW_FM_OCCURANCES_H
#define AW_FM_OCCURANCES_H

#include "AwFmIndex.h"
#include "AwFmSearch.h"
#include <stdint.h>

//deprecated?
struct AwFmOccurrencePair{
  uint64_t occurrence;
  uint64_t occurrenceLessThan;
};
struct AwFmOccurrenceVectorPair{
  __m256i occurrenceVector;
  __m256i occurrenceGteVector;
};

inline void awFmMakeNucleotideOccurrenceVectorPair(const struct AwFmNucleotideBlock *restrict const blockPtr,
  const uint64_t queryPosition, const uint8_t letter, const uint64_t sentinelCharacterPosition,
  struct AwFmOccurrenceVectorPair *occurrenceVectors);

inline void awFmMakeAminoAcidOccurrenceVectorPair(const struct AwFmAminoBlock *restrict const blockPtr,
  const uint64_t queryPosition, const uint8_t letter, struct AwFmOccurrenceVectorPair *vectorPair);

inline uint_fast8_t awFmVectorPopcount(const __m256i occurrenceVector);

inline void awFmBlockPrefetch(const uint8_t *restrict const baseBlockListPtr, const uint64_t blockByteWidth,
  const uint64_t nextQueryPosition);

inline uint8_t awFmGetLetterAtBwtPosition(const union AwFmBwtBlockList blockList, const enum AwFmAlphabetType alphabet, const uint64_t bwtPosition);


// uint64_t awFmGetOccurrences(const struct AwFmIndex *restrict const index,
//   const size_t queryPosition, const uint8_t letter);
//
// struct AwFmOccurrencePair awFmGetAminoAcidOccurrencePair(const struct AwFmAminoBlock *restrict const blockList,
//   const uint64_t *restrict const rankPrefixSums, const size_t queryPosition, const uint8_t letter);
//
//
//   //deprecated
//   // void  awFmOccurrencesDataPrefetch(const struct AwFmIndex *restrict const index, const uint64_t queryPosition);
// void awFmAminoAcidDataPrefetch(const struct AwFmAminoBlock *restrict const blockList, const uint64_t queryPosition);
//
//
// uint64_t awFmGetAminoAcidOccurrence(const struct AwFmAminoBlock *restrict const blockList, const size_t queryPosition, const uint8_t letter);
//
uint8_t awFmGetAminoAcidLetterAtBwtPosition(const struct AwFmAminoBlock *restrict const blockList, const uint64_t bwtPosition);
//
size_t awFmBackstepBwtPosition(const struct AwFmIndex *restrict const index, const uint64_t bwtPosition);

#endif /* end of include guard: AW_FM_OCCURANCES_H */
