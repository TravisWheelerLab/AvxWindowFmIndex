#ifndef AW_FM_INDEX_STRUCTS_H
#define AW_FM_INDEX_STRUCTS_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "FastaVector.h"

#ifdef __cplusplus
#define _RESTRICT_ __restrict__
#else
#define _RESTRICT_ restrict
#endif

#ifndef AW_FM_NUM_CONCURRENT_QUERIES
#define AW_FM_NUM_CONCURRENT_QUERIES 8
#endif

#define AW_FM_POSITIONS_PER_FM_BLOCK 256
#define AW_FM_CACHE_LINE_SIZE_IN_BYTES 64

#define AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW 3
#define AW_FM_NUCLEOTIDE_CARDINALITY 4

#define AW_FM_AMINO_VECTORS_PER_WINDOW 5
#define AW_FM_AMINO_CARDINALITY 20

enum AwFmAlphabetType {
  AwFmAlphabetAmino = 1,
  AwFmAlphabetDna = 2,
  AwFmAlphabetRna = 3
};

// NOTE: not currently used, but this enum is kept for future use.
enum AwFmBwtType { AwFmBwtTypeBackwardOnly = 1, AwFmBwtTypeBiDirectional = 2 };

// define the Simd vector type, which is determined by the architecture we're
// building for.
#ifdef __aarch64__
#include <arm_neon.h>

typedef struct AwFmSimdVec256 {
  uint8x16_t lowVec;
  uint8x16_t highVec;
} AwFmSimdVec256;

#else
#include <immintrin.h>

typedef __m256i AwFmSimdVec256;
#endif

// Types for the actual FM index structs
struct AwFmAminoBlock {
  AwFmSimdVec256 letterBitVectors[AW_FM_AMINO_VECTORS_PER_WINDOW];
  uint64_t baseOccurrences[AW_FM_AMINO_CARDINALITY +
                           4]; //+4 is for sentinel count and 32B padding
};

struct AwFmNucleotideBlock {
  AwFmSimdVec256 letterBitVectors[AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW];
  uint64_t baseOccurrences[AW_FM_NUCLEOTIDE_CARDINALITY +
                           4]; //+4 is for sentinel count and 32B padding
};

union AwFmBwtBlockList {
  struct AwFmNucleotideBlock *asNucleotide;
  struct AwFmAminoBlock *asAmino;
};

/*Struct for the configuration in the AwFmIndex struct.
 * This contains data that may be set by the user toconfigure the index.*/
struct AwFmIndexConfiguration {
  uint8_t suffixArrayCompressionRatio;
  uint8_t kmerLengthInSeedTable;
  enum AwFmAlphabetType alphabetType;
  bool keepSuffixArrayInMemory;
  bool storeOriginalSequence;
};

struct AwFmCompressedSuffixArray {
  uint8_t valueBitWidth;
  uint8_t *values;
  uint64_t compressedByteLength;
};

struct AwFmSearchRange {
  uint64_t startPtr;
  uint64_t endPtr;
};

// feature flags, hardcode version
struct AwFmIndex {
  uint32_t versionNumber;
  uint32_t featureFlags; // for non user-customizable options.
  uint64_t bwtLength;
  union AwFmBwtBlockList bwtBlockList;
  uint64_t *prefixSums;
  struct AwFmSearchRange *kmerSeedTable;
  FILE *fileHandle;
  struct AwFmIndexConfiguration config;
  int fileDescriptor;
  size_t suffixArrayFileOffset;
  size_t sequenceFileOffset;
  // optional member data, dependant on the index version.
  struct FastaVector *fastaVector; // ptr should be null if not in use.
  struct AwFmCompressedSuffixArray suffixArray;
};

struct AwFmKmerSearchData {
  char *kmerString;
  uint64_t kmerLength;
  uint64_t *positionList;
  uint32_t count;
  uint32_t capacity;
};

struct AwFmKmerSearchList {
  size_t capacity;
  size_t count;
  struct AwFmKmerSearchData *kmerSearchData;
};

// for internal use during backtrace, you can likely ignore this
struct AwFmBacktrace {
  uint64_t position;
  uint64_t offset;
};

/* clang-format off */
enum AwFmReturnCode{
  AwFmSuccess             = 1,    AwFmFileReadOkay                = 2,    AwFmFileWriteOkay         = 3,
  AwFmGeneralFailure      = -1,   AwFmUnsupportedVersionError     = -2,   AwFmAllocationFailure     = -3,
  AwFmNullPtrError        = -4,   AwFmSuffixArrayCreationFailure  = -5,   AwFmIllegalPositionError  = -6,
  AwFmNoFileSrcGiven      = -7,   AwFmNoDatabaseSequenceGiven     = -8,   AwFmFileFormatError       = -9,
  AwFmFileOpenFail        = -10,  AwFmFileReadFail                = -11,  AwFmFileWriteFail         = -12,
  AwFmErrorDbSequenceNull = -13,  AwFmErrorSuffixArrayNull        = -14,  AwFmFileAlreadyExists     = -15};
/* clang-format on */

/*
 * Function:  awFmCreateIndex
 * --------------------
 * Allocates a new AwFmIndex from the sequence using the given configuration.
 *
 *  Inputs:
 *    index:          Double pointer to a AwFmIndex struct to be allocated and
 * constructed. configuration:       Fully initialized index config to construct
 * the index with. This configuration will be memcpy'd into the created index.
 *    sequence:       Database sequence that the AwFmIndex is built from.
 *    sequenceLength: Length of the sequence.
 *    fileSrc:        File path to write the Index file to.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the write. Possible returns are:
 *      AwFmFileWriteOkay on success.
 *      AwFmNullPtrError on passing an argument as a null ptr
 *      AwFmAllocationFailure if memory could not be allocated during the
 * creation process. AwFmFileAlreadyExists if a file exists at the given
 * fileSrc, but allowOverwite was false. AwFmSuffixArrayCreationFailure if an
 * error was caused by divsufsort64 in suffix array creation. AwFmFileWriteFail
 * if a file write failed.
 */
enum AwFmReturnCode
awFmCreateIndex(struct AwFmIndex *_RESTRICT_ *index,
                struct AwFmIndexConfiguration *_RESTRICT_ const config,
                const uint8_t *_RESTRICT_ const sequence,
                const size_t sequenceLength,
                const char *_RESTRICT_ const fileSrc);

/*
 * Function:  awFmCreateIndexFromFasta
 * --------------------
 * Loads the sequence and header data from the given fasta, and a allocates a
 * new AwFmIndex from the sequence using the given configuration.
 *
 *  Inputs:
 *    index:          Double pointer to a AwFmIndex struct to be allocated and
 * constructed. configuration:       Fully initialized config struct to
 * construct the index with. This config will be memcpy'd into the created
 * index. fastaSrc:       File source of the fasta to use to generate the index.
 *      Every sequence in the fasta file will be included in the index.
 *    indexFileSrc:        File path to write the Index file to.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the write. Possible returns are:
 *      AwFmFileWriteOkay on success.
 *      AwFmNullPtrError on passing an argument as a null ptr
 *      AwFmFileOpenFail if the fasta cannot be opened for reading.
 *      AwFmAllocationFailure if memory could not be allocated during the
 * creation process. AwFmFileAlreadyExists if a file exists at the given
 * fileSrc, but allowOverwite was false. AwFmSuffixArrayCreationFailure if an
 * error was caused by divsufsort64 in suffix array creation. AwFmFileWriteFail
 * if a file write failed.
 */
enum AwFmReturnCode
awFmCreateIndexFromFasta(struct AwFmIndex *_RESTRICT_ *index,
                         struct AwFmIndexConfiguration *_RESTRICT_ const config,
                         const char *fastaSrc,
                         const char *_RESTRICT_ const indexFileSrc);

/*
 * Function:  awFmDeallocIndex
 * --------------------
 * Deallocates the given AwFmIndex struct, as well as all internal arrays.
 *   Performing this deallocation will also close the internally stored file
 * handle.
 *
 *  Inputs:
 *    index:  Pointer to the AwFmIndex struct to deallocate
 */
void awFmDeallocIndex(struct AwFmIndex *index);

/*
 * Function:  awFmWriteIndexToFile
 * --------------------
 * With a given AwFmIndex and associated data, writes the index to the file.
 *  If you, the user, want to create a new index file, use the awFmCreateIndex
 *  function in AwFmCreate.h
 *
 *  Inputs:
 *    index:          AwFmIndex struct to be written.
 *    suffixArray:    Full (uncompressed) suffix array that the AwFmIndex is
 * built from. sequence:       Database sequence that the AwFmIndex is built
 * from. sequenceLength: Length of the sequence. fileSrc:        File path to
 * write the Index file to.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the write. Possible returns are:
 *      AwFmFileWriteOkay on success.
 *      AwFmFileAlreadyExists if a file exists at the given fileSrc, but
 * allowOverwite was false. AwFmFileWriteFail if a file write failed.
 */
enum AwFmReturnCode
awFmWriteIndexToFile(struct AwFmIndex *_RESTRICT_ const index,
                     const uint8_t *_RESTRICT_ const sequence,
                     const uint64_t sequenceLength,
                     const char *_RESTRICT_ const fileSrc);

/*
 * Function:  awFmReadIndexFromFile
 * --------------------
 * Reads the AwFmIndex file from the given fileSrc
 *
 *  Inputs:
 *    index:                    Double pointer to an unallocated AwFmIndex to be
 *        allocated and populated by this function.
 *    fileSrc:                  Path to the file containing the AwFmIndex.
 *    keepSuffixArrayInMemory:  flag to determine whether to read the compressed
 *        suffix array into memory, or leave it on disk.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the write. Possible returns are:
 *      AwFmFileReadOkay on success.
 *      AwFmFileAlreadyExists if no file could be opened at the given fileSrc.
 *      AwFmFileFormatError if the header was not correct. This suggests that
 * the file at this location is not the correct format. AwFmAllocationFailure on
 * failure to allocated the necessary memory for the index.
 */
enum AwFmReturnCode
awFmReadIndexFromFile(struct AwFmIndex *_RESTRICT_ *_RESTRICT_ index,
                      const char *fileSrc, const bool keepSuffixArrayInMemory);

/*
 * Function:  awFmFindSearchRangeForString
 * --------------------
 *  Queries the FM-Index for the range of BWT positions that represent instances
 *    of the given Kmer found in the database.
 *
 *  If the given kmer is not found, the AwFmSearch Range will result in a false
 * value when given to the awFmSearchRangeIsValid function.
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
 *    AwFmSearchRange representing the range of BWT positions where the given
 *    kmer may be found, as long as startPtr < endPtr. Otherwise (startPtr >=
 * endPtr), the given kmer does not exist in the database sequence.
 */
struct AwFmSearchRange
awFmFindSearchRangeForString(const struct AwFmIndex *_RESTRICT_ const index,
                             const char *_RESTRICT_ const kmer,
                             const size_t kmerLength);

/*
 * Function:  awFmCreateKmerSearchList
 * --------------------
 *  Allocates and initializes an AwFmKmerSearchList struct to be used to search
 *  for groups of kmers in a thread-parallel manner that also hides memory read
 * latency through multiple concurrent queries per thread. This struct can be
 * used for both count and locate functions.
 *
 *  Note that the kmers inside the searchList are not allocated, and only
 * contain char pointers that can be set to the kmers you want to query for.
 *
 *  Inputs:
 *    capacity:     How many kmers the searchData struct can hold.
 *
 *  Returns:
 *    Pointer to the allocated searchData struct, or null on failure.
 */
struct AwFmKmerSearchList *awFmCreateKmerSearchList(const size_t capacity);

/*
 * Function:  awFmDeallocKmerSearchList
 * --------------------
 *  Deallocates the given search list struct, and both the internal searchData
 * list and positionLists inside the searchData.
 *
 *  Note that, since the searchData doesn't own the kmer char strings, it will
 * not try to deallocate them. If kmers were dynamically allocated externally,
 * it is the caller's responsibility to deallocate them. This means that if
 * these pointers are the only pointer to the data, and if they were dynamically
 * allocated, forgetting to deallocate them before calling this function will
 * leak the data.
 *
 *  Inputs:
 *    searchData:   pointer to the searchData struct to deallocate
 */
void awFmDeallocKmerSearchList(
    struct AwFmKmerSearchList *_RESTRICT_ const searchList);

/*
 * Function:  awFmParallelSearchLocate
 * --------------------
 *  Using the given index and a searchData struct preloaded with kmers, query
 *the kmers in a concurrent, thread-parallel manner to find the positions of
 *each occurrence of each respective kmer. The suggested use case for this
 *function is as follows:
 *
 *    1. Allocate a searchList struct with awFmCreateKmerSearchList().
 *    2. For searching for n kmers, set count to n (must be smaller than
 *capacity!), and set the first n kmers with the correct char pointer and
 *length.
 *    3. Call this function using the index to search. The positions of each
 *kmer will be loaded into the positionList member variable of the corresponding
 *AwFmKmerSearchData struct.
 *    4. To query for additional kmers, reuse the searchList struct, starting
 *with step (2).
 *    5. Deallocate the searchData with awFmDeallocKmerSearchList when finished.
 *
 *  IMPORTANT NOTE:
 *    The locate and count functions cannot search for ambiguity characters, and
 *      providing ambiguity characteres in the kmers will result in undefined
 *behavior. ensure that the query kmers contain only nucleotide or amino acid
 *characters, depending on the alphabet being used.
 *
 *  Inputs:
 *    index:        pointer to the index to search.
 *    searchList:   pointer to the searchList struct loaded with kmers to search
 *for. numThreads:   How many threads to direct OpenMP to use. The best value
 *for this argument will likely vary from system to system. Suggested default
 *value is 4 Returns: AwFmReturnCode represnting the result of the read.
 *Possible returns are: AwFmFileReadOkay on success. AwFmFileReadFail if the
 *file could not be read sucessfully (If suffix array is stored on file, not in
 *memory)
 */
enum AwFmReturnCode
awFmParallelSearchLocate(const struct AwFmIndex *_RESTRICT_ const index,
                         struct AwFmKmerSearchList *_RESTRICT_ const searchList,
                         uint32_t numThreads);

/*
 * Function:  awFmParallelSearchCount
 * --------------------
 *  Using the given index and a searchData struct preloaded with kmers, query
 * the kmers in a concurrent, thread-parallel manner to find the count of
 * occurrences of each respective kmer. The suggested use case for this function
 * is as follows:
 *
 *    1. Allocate a searchList struct with awFmCreateKmerSearchList().
 *    2. For searching for n kmers, set count to n (must be smaller than
 * capacity!), and set the first n kmers with the correct char pointer and
 * length.
 *    3. Call this function using the index to search. The count of matching
 * kmers will be loaded into the AwFmKmerSearchData's count member variable.
 *    4. To query for additional kmers, reuse the searchList struct, starting
 * with step (2).
 *    5. Deallocate the searchData with awFmDeallocKmerSearchList when finished.
 *
 *  IMPORTANT NOTE:
 *    The locate and count functions cannot search for ambiguity characters, and
 *      providing ambiguity characteres in the kmers will result in undefined
 * behavior. ensure that the query kmers contain only nucleotide or amino acid
 * characters, depending on the alphabet being used.
 *
 *  Inputs:
 *    index:        pointer to the index to search.
 *    searchList:   pointer to the searchList struct loaded with kmers to search
 * for. numThreads:   How many threads to direct OpenMP to use. The best value
 * for this argument will likely vary from system to system. Suggested default
 * value is 4
 */
void awFmParallelSearchCount(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmKmerSearchList *_RESTRICT_ const searchList,
    uint32_t numThreads);

/*
 * Function:  awFmReadSequenceFromFile
 * --------------------
 * Given a sequence position, reads a section of sequence surrounding that
 *position from the corresponding index file.
 *
 *  Inputs:
 *    index:          Pointer to the AwFmIndex struct that contains the file
 *handle to read. index file that stores the compressed suffix array.
 *    sequenceStartPosition:  Position in the sequence to begin reading to fill
 *the provided buffer. sequenceSegmentLength:  Length of the sequence segment to
 *read. sequenceBuffer: Pointer to the buffer to read the sequence segment into.
 *This buffer  must be large enough to hold (sequenceStartPosition -
 *sequenceEndPosition +1) characters. the + 1 on this length is for the null
 *terminator added to the end of the string.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the read. Possible returns are:
 *    AwFmFileReadOkay on success.
 *    AwFmFileReadFail if the file could not be read sucessfully.
 *    AwFmIllegalPositionError if the start position is not less than the end
 *position AwFmUnsupportedVersionError if the index was configured to not store
 *the original sequence.
 */
enum AwFmReturnCode
awFmReadSequenceFromFile(const struct AwFmIndex *_RESTRICT_ const index,
                         const size_t sequenceStartPosition,
                         const size_t sequenceSegmentLength,
                         char *const sequenceBuffer);

/*
 * Function:  awFmCreateInitialQueryRange
 * --------------------
 * Creates the initial Start Pointer - End Pointer range for the given query.
 *   This range represents all positions in the suffix array that correspond
 *   to the last character in the given query.
 *   The alphabet (nucleotide or amino acid) is determined by the alphabet of
 * the index.
 *
 *   NOTE: This function is likely only useful if you need to query on a
 * letter-by-letter basis, e.g., to implement inexact matching. If you just want
 * to query for exact matches, especially for many queries, use
 * awFmParallelSearchLocate or awFmParallelSearchCount instead.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    query: ASCII string of the query. The last character of this string will
 * be used to generate the initial range. queryLength: length of the query
 * argument. This value is used to find the last character in the input query.
 */
struct AwFmSearchRange
awFmCreateInitialQueryRange(const struct AwFmIndex *_RESTRICT_ const index,
                            const char *_RESTRICT_ const query,
                            const uint64_t queryLength);

/*
 * Function:  awFmCreateInitialQueryRangeFromChar
 * --------------------
 * Creates the initial Start Pointer - End Pointer range for the given query.
 *   This is a simpler version of the awFmCreateInitialQueryRange function for
 * when you just want to start the range based off a single character
 *
 *   NOTE: This function is likely only useful if you need to query on a
 * letter-by-letter basis, e.g., to implement inexact matching. If you just want
 * to query for exact matches, especially for many queries, use
 * awFmParallelSearchLocate or awFmParallelSearchCount instead.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    letter: letter to begin the search with. This will represent the final
 * chacter in a searched kmer.
 */
struct AwFmSearchRange awFmCreateInitialQueryRangeFromChar(
    const struct AwFmIndex *_RESTRICT_ const index, const char letter);

/*
 * Function:  awFmIterativeStepBackwardSearch
 * --------------------
 * Performs a single backward search step on the given index.
 *  In lieu of returning an additional value, this function updates the data
 * pointed to by the range ptr.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    range: range in the BWT that corresponds to the implicit kmer that is
 * about to be extended. this acts as an out-parameter, and will update to the
 * newly extended range once finished. letterIndex: letter index of the prefix
 * or suffix character, between 0 and 3.
 */
void awFmNucleotideIterativeStepBackwardSearch(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmSearchRange *_RESTRICT_ const range, const uint8_t letterIndex);

/*
 * Function:  awFmIterativeStepBackwardSearch
 * --------------------
 * Performs a single backward search step on the given index.
 *  In lieu of returning an additional value, this function updates the data
 * pointed to by the range ptr.
 *
 *  Inputs:
 *    index: AwFmIndex struct to search
 *    range: range in the BWT that corresponds to the implicit kmer that is
 * about to be extended. this acts as an out-parameter, and will update to the
 * newly extended range once finished. letterIndex: letter index of the suffix
 * character, between 0 and 19.
 */
void awFmAminoIterativeStepBackwardSearch(
    const struct AwFmIndex *_RESTRICT_ const index,
    struct AwFmSearchRange *_RESTRICT_ const range, const uint8_t letterIndex);

/*
 * Function:  awFmFindDatabaseHitPositions
 * --------------------
 *  Takes a range of BWT positions, backtraces each position to find the nearest
 * sample in the compressed suffix array, and looks up those suffix array
 * positions on disk to determine the corresponding database sequence position
 * for each BWT position between the searchRange's pointers (inclusive startPtr,
 * exclusive endPtr).
 *
 *
 *  It is the caller's responsibility to free() the returned sequence position
 * array and offset array.
 *
 *
 *
 *  Inputs:
 *    index:              Pointer to the valid AwFmIndex struct.
 *    searchRange:        Range in the index. These BWT positions will be
 * converted to real global sequence positions fileAssessResult:
 * Returns the result of this action as an out-variable. possible return values
 * are: AwFmFileReadOkay on success AwFmGeneralFailure if the search range did
 * not represent any values (aka SP > EP) AwFmFileReadFail on failure to read
 * from the position array (if left on file and not in memory)
 * 			AwFmAllocationFailure if the function could not allocate
 * memory for the return array.
 *
 *  Returns:
 * 		Dynamically allocated array of hit positions. the length of the
 * array will be the same length as the search range. This length can be easily
 * determined with the awFmSearchRangeLength() function
 */
uint64_t *awFmFindDatabaseHitPositions(
    const struct AwFmIndex *_RESTRICT_ const index,
    const struct AwFmSearchRange *_RESTRICT_ const searchRange,
    enum AwFmReturnCode *_RESTRICT_ fileAccessResult);

/*
 * Function:  awFmFindDatabaseHitPositionSingle
 * --------------------
 *  Backtraces a single BWT position to find its position in the sequence.
 * 		If the Suffix Array is stored on disk, this will occur a single
 * disk read. If your goal is to query many BWT positions, use
 * awFmFindDatabaseHitPositions instead.
 *
 *
 *
 *  Inputs:
 *    index:              	Pointer to the valid AwFmIndex struct.
 *    bwtPosition:			Position in the bwt to query
 *    fileAssessResult:		Returns the result of this action as an
 * out-variable. possible return values are: AwFmFileReadOkay on success
 * 			AwFmFileReadFail on failure to read from the position
 * array (if left on file and not in memory)
 *
 *  Returns:
 * 		Position in the original sequence that corresponds to the bwt
 * position
 */
uint64_t awFmFindDatabaseHitPositionSingle(
    const struct AwFmIndex *_RESTRICT_ const index, const uint64_t bwtPosition,
    enum AwFmReturnCode *_RESTRICT_ fileAccessResult);

/*
 * Function:  awFmGetLocalSequencePositionFromIndexPosition
 * --------------------
 *  For indices that are built from a fasta file (indices that use an internal
 * FastaVector struct), this function takes the global position into the full
 * indexed sequence collection, and returns the index of the individual sequence
 * that the global position lands in and the local position in that sequence
 * that corresponds to the global position.
 *
 *  Inputs:
 *     index:                  Pointer to the valid AwFmIndex struct.
 *     globalPosition:         Position into the sequence collection, as is
 * stored in the suffix array. sequenceNumber:         Out-argument that returns
 * the index of the sequence the global position falls into, or 0 on error.
 *     localSequencePosition:  Out-argument that returns the local position in
 * the relevant sequence, or 0 on error.
 *
 *  Returns:
 *    AwFmReturnCode detailing the result of the read attempt. Possible return
 * values: AwFmSuccess on success, AwFmIllegalPositionError if the
 * globalPosition is greater than the length of the index.
 *      AwFmUnsupportedVersionError if the version of the AwFmIndex does not
 * support FastaVector.
 */
enum AwFmReturnCode awFmGetLocalSequencePositionFromIndexPosition(
    const struct AwFmIndex *_RESTRICT_ const index, size_t globalPosition,
    size_t *sequenceNumber, size_t *localSequencePosition);

/*
 * Function:  awFmNucleotideBacktraceReturnPreviousLetterIndex
 * --------------------
 *  Backtraces the given position in the bwt, and returns the the previous
 * character to that position. The given BWT position will then be updated to
 * refer to the position of the returned previous character.
 *
 *  Inputs:
 *     index:			Pointer to the valid AwFmIndex struct.
 *     bwtPosition:     A valid position in the fm index
 *
 *  Returns:
 *    Letter index of the character found at the given position in the bwt. In
 * other words, the letter previous in the sequence to the given position.
 */
uint8_t awFmNucleotideBacktraceReturnPreviousLetterIndex(
    const struct AwFmIndex *_RESTRICT_ const index, uint64_t *bwtPosition);

/*
 * Function:  awFmAminoBacktraceReturnPreviousLetterIndex
 * --------------------
 *  Backtraces the given position in the bwt, and returns the the previous
 * character to that position. The given BWT position will then be updated to
 * refer to the position of the returned previous character.
 *
 *  Inputs:
 *     index:			Pointer to the valid AwFmIndex struct.
 *     bwtPosition:     A valid position in the fm index
 *
 *  Returns:
 *    Letter index of the character found at the given position in the bwt. In
 * other words, the letter previous in the sequence to the given position.
 */
uint8_t awFmAminoBacktraceReturnPreviousLetterIndex(
    const struct AwFmIndex *_RESTRICT_ const index, uint64_t *bwtPosition);

/*
 * Function:  awFmGetHeaderStringFromSequenceNumber
 * --------------------
 *  For indices that are built from a fasta file (indices that use an internal
 * FastaVector struct), this function returns a pointer to the header and the
 * length of said for for the given sequenceNumber as out-arguments.
 *
 *  Inputs:
 *     index:            Pointer to the valid AwFmIndex struct.
 *     sequenceNumber:   Index of the sequence whose corresponding header the
 * function will return. headerBuffer:     Upon success, this function returns
 * the pointer to the (non null-terminated) header in this out-argument.
 * headerLength:     Upon success, this fucntion reuturns the length of the
 * header in this out-argument.
 *
 *  Returns:
 *    AwFmReturnCode detailing the result of the read attempt. Possible return
 * values: AwFmFileReadOkay on success, AwFmUnsupportedVersionError if the
 * version of the AwFmIndex does not support FastaVector.
 */
enum AwFmReturnCode awFmGetHeaderStringFromSequenceNumber(
    const struct AwFmIndex *_RESTRICT_ const index, size_t sequenceNumber,
    char **headerBuffer, size_t *headerLength);

/*
 * Function:  awFmSearchRangeLength
 * --------------------
 * Gets the number of positions included in the given AwFmSearchRange
 *
 *  Inputs:
 *    range: Range of positions in the BWT that corresponds to some number of
 *      instances of a given kmer.
 *
 *  Returns:
 *    Number of positions in the given range if the range is valid (startPtr <
 * endPtr), or 0 otherwise, as that would imply that no instances of that kmer
 * were found.
 */
size_t
awFmSearchRangeLength(const struct AwFmSearchRange *_RESTRICT_ const range);

/*
 * Function:  awFmReturnCodeIsFailure
 * --------------------
 *  For a given return code, returns true if the code describes an error
 * condition, or true on a successful condititon.
 *
 *  Inputs:
 *     rc:            AwFmReturnCode returned by an AwFm function.
 *  Returns:
 *    True if return code describes an error of some kind.
 */
bool awFmReturnCodeIsFailure(const enum AwFmReturnCode rc);

/*
 * Function:  awFmReturnCodeIsFailure
 * --------------------
 *  returns the negation of awFmReturnCodeIsFailure().
 *
 *  Inputs:
 *     rc:            AwFmReturnCode returned by an AwFm function.
 *  Returns:
 *    True if return code describes a successful operation.
 */
bool awFmReturnCodeIsSuccess(const enum AwFmReturnCode rc);

/*
 * Function:  awFmGetNumSequences
 * --------------------
 * Returns the number of sequences used to generate the awfmIndex.
 * When this index is built from a fasta, this will be equal to the number
 *  of sequences in that fasta. When built from a text string, this will be one.
 *  Inputs:
 *    index: AwFmIndex to query.
 *
 *  Returns:
 *    Number of sequences contained inside the index.
 */
uint32_t awFmGetNumSequences(const struct AwFmIndex *_RESTRICT_ const index);

#endif /* end of include guard: AW_FM_INDEX_STRUCTS_H */
