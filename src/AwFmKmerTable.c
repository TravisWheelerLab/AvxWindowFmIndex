#include "AwFmKmerTable.h"
#include "AwFmLetter.h"

bool awFmQueryCanUseKmerTable(const struct AwFmIndex *_RESTRICT_ const index,
                              const char *_RESTRICT_ const kmer,
                              const size_t kmerLength) {

  if (kmerLength < index->config.kmerLengthInSeedTable) {
    return false;
  }
  for (size_t letterIdx = kmerLength - index->config.kmerLengthInSeedTable;
       letterIdx < kmerLength; letterIdx++) {
    if (awFmLetterIsAmbiguous(kmer[letterIdx], index->config.alphabetType)) {
      return false;
    }
  }

  return true;
}

struct AwFmSearchRange awFmNucleotideKmerSeedRangeFromTable(
    const struct AwFmIndex *_RESTRICT_ const index,
    const char *_RESTRICT_ const kmer, const size_t kmerLength) {

  const size_t kmerSeedStartPosition =
      kmerLength - index->config.kmerLengthInSeedTable;
  size_t kmerTableIndex = 0;
  for (size_t i = kmerSeedStartPosition; i < kmerLength; i++) {
    uint8_t letterIndex = awFmAsciiNucleotideToLetterIndex(kmer[i]);
    kmerTableIndex =
        (kmerTableIndex * AW_FM_NUCLEOTIDE_CARDINALITY) + letterIndex;
  }

  return index->kmerSeedTable[kmerTableIndex];
}

struct AwFmSearchRange
awFmAminoKmerSeedRangeFromTable(const struct AwFmIndex *_RESTRICT_ const index,
                                const char *_RESTRICT_ const kmer,
                                const size_t kmerLength) {

  const size_t kmerSeedStartPosition =
      kmerLength - index->config.kmerLengthInSeedTable;
  size_t kmerTableIndex = 0;
  for (size_t i = kmerSeedStartPosition; i < kmerLength; i++) {
    uint8_t letterIndex = awFmAsciiAminoAcidToLetterIndex(kmer[i]);
    kmerTableIndex = (kmerTableIndex * AW_FM_AMINO_CARDINALITY) + letterIndex;
  }

  return index->kmerSeedTable[kmerTableIndex];
}

// definition check for an intentionally undefined variable so this code doesn't
// get implemented, get used, or throw
#ifdef AW_FM_PARTIAL_SEED_CODE_IMPLEMENTATION_PROVIDED

/*

        This code is a snippit of code that may be used in future releases to
search for kmers smaller than what are memoized in the kmerSeedTable.

        The main gist is that if your seed table memoizes length 4 kmers, and
you have a query "gc", you can use the memoized ranges for "gcaa" and "gctt" and
stitch them together to make the range for the smaller query kmer.

        This gets more complicated, since the sentinel at the end can result in
some smaller kmers at the very end of the sequence, for instance if the sequence
ends with "gc", it'll be after the endPtr for the "gctt" range. You can keep
track of some representation of the last k-1 letters in the sequence, where k is
the lengh of kmers memoized, and just check in there. This could be a char *
string, and use strncmp, but a more efficient option is likely to just create
the table index of the last letters, and bit shift it to remove letters to
check.

        This all gets harder when ambiguity characters are added, but not
memoized. All of a sudden, kmers that end with the greatest nucleotide ('t') or
amino ("y"), can't easily tell where "gttt" ends and "gxaaa" begins.

        This problem was hard, and a better solution than "just ignore all that
don't end in "t" or "y" was just too hacky to implement in the first official
open beta release.


Ideas as to how to fix this and get it working:

        So, if the query kmer ends with the largest valid letter ('t' or 'y'),
this current technique will fail. This is because it's impossible to separate
out the endPtr for these kmers, since we're not storing the ranges for things
that include ambiguity characters.

        One suggested solution is to make the kmer seed table bigger by also
memoizing kmers that contain ambiguity characters

        The size of the tables will increase, and as a result, the size of the
seeds will likely have to be 1 smaller in practice. instead of 4^n * 16 bytes
for a nucleotide table with length n, the new size would be 5^n * 16.  However,
at that point, searching for the smaller kmers becomes much easier.

        For Aminos, the same thing would happen (21^k * 16, instead of 20).



        Another future update: reverse the direction that the kmerSeedTable
index is generated. instead of having the left-most character contribute the
largest bits, the right-most have the largest bits. I think this would result in
more contiguous memory writes during index creation, but might also help with
the indexing for ideas listed above.

*/
static inline struct AwFmSearchRange
awFmNucleotidePartialKmerSeedRangeFromTable(
    const struct AwFmIndex *_RESTRICT_ const index,
    const char *_RESTRICT_ const kmer, const uint8_t kmerLength) {

  assert(false, "this function is unfinished, untested, and left only to be "
                "finished for future releases");
  struct AwFmSearchRange range = {0, 0};

  // this check is only for the current assumption that the kmers can't end with
  // the last letter
  //(due to being unable to tell) where the partial kmer range ends and the
  //range containing ambig. chararacters begins.
  if (__builtin_expect((kmer[kmerLength - 1] | 0x20) == 't', 0)) {
    // lookup kmer from scratch
    awFmNucleotideNonSeededSearch(index, kmer, kmerLength, &range);
  } else {

    const uint64_t endSequenceKmerTableIndex =
        0; // index->endSequenceKmerTableIndex;
    const uint8_t kmerLengthInSeedTable = index->config.kmerLengthInSeedTable;
    const uint8_t nucleotideBitOffset = 2;
    uint64_t queryKmerStartTableIndex = 0;
    uint64_t queryKmerEndTableIndex = 0;

    for (uint8_t kmerLetterIndex = 0; kmerLetterIndex < kmerLength;
         kmerLetterIndex++) {
      queryKmerStartTableIndex <<= nucleotideBitOffset;
      queryKmerStartTableIndex |=
          awFmAsciiNucleotideToLetterIndex(kmer[kmerLetterIndex]);
    }

    // adds 1 to the last character, this will never overflow the last
    // character, since kmers that end with a 't' must use non-seeded search.
    queryKmerEndTableIndex = queryKmerStartTableIndex + 1;
    uint64_t indexBitmask = (1 << (kmerLength * nucleotideBitOffset)) - 1;

    for (uint8_t i = kmerLength; i < kmerLengthInSeedTable; i++) {
      bool baseMatchesEndingIndex =
          ((endSequenceKmerTableIndex & indexBitmask) ^
           queryKmerStartTableIndex) == 0;
      bool endMatchesEndingIndex = ((endSequenceKmerTableIndex & indexBitmask) ^
                                    queryKmerEndTableIndex) == 0;

      if (baseMatchesEndingIndex) {
        range.startPtr--;
      }
      if (endMatchesEndingIndex) {
        range.endPtr--;
      }
      indexBitmask = (indexBitmask << nucleotideBitOffset) | 0x3;
      queryKmerStartTableIndex <<= nucleotideBitOffset;
      queryKmerEndTableIndex <<= nucleotideBitOffset;
    }

    // add in the start pointers from the ranges extended with a's.
    range.startPtr += index->kmerSeedTable[queryKmerStartTableIndex].startPtr;
    range.endPtr += index->kmerSeedTable[queryKmerEndTableIndex].startPtr - 1;
  }

  return range;
}

inline struct AwFmSearchRange awFmAminoPartialKmerSeedRangeFromTable(
    const struct AwFmIndex *_RESTRICT_ const index,
    const char *_RESTRICT_ const kmer, const uint8_t kmerLength) {
  struct AwFmSearchRange range = {0, 0};

  if (__builtin_expect((kmer[kmerLength - 1] | 0x20) == 'y', 0)) {
    // lookup kmer from scratch
    awFmAminoNonSeededSearch(index, kmer, kmerLength, &range);
  } else {

    const uint64_t endSequenceKmerTableIndex =
        0; // index->endSequenceKmerTableIndex;
    const uint8_t kmerLengthInSeedTable = index->config.kmerLengthInSeedTable;
    uint64_t queryKmerStartTableIndex = 0;
    uint64_t queryKmerEndTableIndex = 0;
    uint64_t indexModFactor = 1;

    for (uint8_t kmerLetterIndex = 0; kmerLetterIndex < kmerLength;
         kmerLetterIndex++) {
      indexModFactor *= AW_FM_AMINO_CARDINALITY;
      queryKmerStartTableIndex *= AW_FM_AMINO_CARDINALITY;
      queryKmerStartTableIndex |=
          awFmAsciiAminoAcidToLetterIndex(kmer[kmerLetterIndex]);
    }

    // adds 1 to the last character, this will never overflow the last
    // character, since kmers that end with a 'y' must use non-seeded search.
    queryKmerEndTableIndex = queryKmerStartTableIndex + 1;

    for (uint8_t i = kmerLength; i < kmerLengthInSeedTable; i++) {
      uint64_t endSequenceKmerTableIndexMod =
          endSequenceKmerTableIndex % indexModFactor;
      bool baseMatchesEndingIndex =
          (queryKmerStartTableIndex == endSequenceKmerTableIndexMod);
      bool endMatchesEndingIndex =
          (queryKmerEndTableIndex == endSequenceKmerTableIndexMod);

      if (baseMatchesEndingIndex) {
        range.startPtr--;
      }
      if (endMatchesEndingIndex) {
        range.endPtr--;
      }
      indexModFactor *= AW_FM_AMINO_CARDINALITY;
      queryKmerStartTableIndex *= AW_FM_AMINO_CARDINALITY;
      queryKmerEndTableIndex *= AW_FM_AMINO_CARDINALITY;
    }

    // add in the start pointers from the ranges extended with a's.
    range.startPtr += index->kmerSeedTable[queryKmerStartTableIndex].startPtr;
    range.endPtr += index->kmerSeedTable[queryKmerEndTableIndex].startPtr - 1;
  }

  return range;
}

#endif
