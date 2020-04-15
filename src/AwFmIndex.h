#ifndef AW_FM_INDEX_STRUCTS_H
#define AW_FM_INDEX_STRUCTS_H

#include <stdint.h>
#include <stdbool.h>
#include <immintrin.h>
#include <stdio.h>

struct AwFmBacktrace{
  uint64_t position;
  uint64_t _offset;
};

struct AwFmKmerSearchData{
  char                        *kmerString;
  uint64_t                    kmerLength;
  struct AwFmBacktrace        *positionBacktraceList;
  uint32_t                    count;
  uint32_t                    capacity;
};

struct AwFmKmerSearchList{
  size_t                      capacity;
  size_t                      count;
  struct AwFmKmerSearchData   *kmerSearchData;
};

struct AwFmSearchRange{
  uint64_t startPtr;
  uint64_t endPtr;
};

#ifndef AW_FM_PREFETCH_STRATEGY
  #define AW_FM_PREFETCH_STRATEGY _MM_HINT_T2
  //#define AW_FM_PREFETCH_STRATEGY _MM_HINT_T1
  // #define AW_FM_PREFETCH_STRATEGY _MM_HINT_T0
  // #define AW_FM_PREFETCH_STRATEGY _MM_HINT_NTA
#endif

#ifndef AW_FM_NUM_CONCURRENT_QUERIES
  #define AW_FM_NUM_CONCURRENT_QUERIES 8
#endif

#define AW_FM_POSITIONS_PER_FM_BLOCK        256
#define AW_FM_CACHE_LINE_SIZE_IN_BYTES      64

#define AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW 2
#define AW_FM_NUCLEOTIDE_CARDINALITY        4

#define AW_FM_AMINO_VECTORS_PER_WINDOW      5
#define AW_FM_AMINO_CARDINALITY             20



enum AwFmAlphabetType{
  AwFmAlphabetAmino = 1, AwFmAlphabetNucleotide = 2};

enum AwFmBwtType{
  AwFmBwtTypeBackwardOnly = 1, AwFmBwtTypeBiDirectional = 2};

/*Structs*/
struct AwFmAminoBlock{
  __m256i   letterBitVectors[AW_FM_AMINO_VECTORS_PER_WINDOW];
  uint64_t  baseOccurrences[AW_FM_AMINO_CARDINALITY];
};

struct AwFmNucleotideBlock{
  __m256i   letterBitVectors[AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW];
  uint64_t  baseOccurrences[AW_FM_NUCLEOTIDE_CARDINALITY];
};

union AwFmBwtBlockList{
  struct AwFmNucleotideBlock  *asNucleotide;
  struct AwFmAminoBlock       *asAmino;
};

/*Struct for the metadata in the AwFmIndex struct.
* This contains data that helps to build the index, or to determine how the index
* will function.*/
struct AwFmIndexMetadata{
  uint16_t              versionNumber;
  uint8_t               suffixArrayCompressionRatio;
  uint8_t               kmerLengthInSeedTable;
  enum AwFmAlphabetType alphabetType;
  bool                  keepSuffixArrayInMemory;
};

struct AwFmKmerSeedTable{
  struct  AwFmSearchRange *table;
          uint64_t        *sequenceEndingKmerEncodings;
};

struct AwFmIndex{
          uint64_t          bwtLength;
          uint64_t          sentinelCharacterPosition;
  union   AwFmBwtBlockList  bwtBlockList;
          uint64_t          *prefixSums;
  struct  AwFmKmerSeedTable kmerSeedTable;
          uint64_t          *inMemorySuffixArray;
          FILE              *fileHandle;
  struct  AwFmIndexMetadata metadata;
          int               fileDescriptor;
          size_t            suffixArrayFileOffset;
          size_t            sequenceFileOffset;
};


//todo: remove unused return codes
enum AwFmReturnCode{
  AwFmSuccess             = 1,    AwFmFileReadOkay                = 2,    AwFmFileWriteOkay         = 3,
  AwFmGeneralFailure      = -1,   AwFmUnsupportedVersionError     = -2,   AwFmAllocationFailure     = -3,
  AwFmNullPtrError        = -4,   AwFmSuffixArrayCreationFailure  = -5,   AwFmIllegalPositionError  = -6,
  AwFmNoFileSrcGiven      = -7,   AwFmNoDatabaseSequenceGiven     = -8,   AwFmFileFormatError       = -9,
  AwFmFileOpenFail        = -10,  AwFmFileReadFail                = -11,  AwFmFileWriteFail         = -12,
  AwFmErrorDbSequenceNull = -13,  AwFmErrorSuffixArrayNull        = -14,  AwFmFileAlreadyExists     = -15};


/*
 * Function:  awFmCreateIndex
 * --------------------
 * Allocates a new AwFmIndex from the sequence using the given metadata configuration.
 *
 *  Inputs:
 *    index:          Double pointer to a AwFmIndex struct to be allocated and constructed.
 *    metadata:       Fully initialized metadata to construct the index with.
 *      This metadata will be memcpy'd into the created index.
 *    sequence:       Database sequence that the AwFmIndex is built from.
 *    sequenceLength: Length of the sequence.
 *    fileSrc:        File path to write the Index file to.
 *    allowOverwrite: If set, will allow overwriting the file at the given fileSrc.
 *      If allowOverwite is false, will return error code AwFmFileAlreadyExists.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the write. Possible returns are:
 *      AwFmFileWriteOkay on success.
 *      AwFmAllocationFailure if memory could not be allocated during the creation process.
 *      AwFmFileAlreadyExists if a file exists at the given fileSrc, but allowOverwite was false.
 *      AwFmSuffixArrayCreationFailure if an error was caused by divsufsort64 in suffix array creation.
 *      AwFmFileWriteFail if a file write failed.
 */
enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex *restrict *index,
  const struct AwFmIndexMetadata *restrict const metadata, const uint8_t *restrict const sequence, const size_t sequenceLength,
  const char *restrict const fileSrc, const bool allowFileOverwrite);



/*
 * Function:  awFmDeallocIndex
 * --------------------
 * Deallocates the given AwFmIndex struct, as well as all internal arrays.
 *   Performing this deallocation will also close the internally stored file handle.
 *
 *  Inputs:
 *    index:  Pointer to the AwFmIndex struct to deallocate
 */
void          awFmDeallocIndex(struct AwFmIndex *index);


/*
 * Function:  awFmWriteIndexToFile
 * --------------------
 * With a given AwFmIndex and associated data, writes the index to the file.
 *  If you, the user, want to create a new index file, use the awFmCreateIndex
 *  function in AwFmCreate.h
 *
 *  Inputs:
 *    index:          AwFmIndex struct to be written.
 *    suffixArray:    Full (uncompressed) suffix array that the AwFmIndex is built from.
 *    sequence:       Database sequence that the AwFmIndex is built from.
 *    sequenceLength: Length of the sequence.
 *    fileSrc:        File path to write the Index file to.
 *    allowOverwrite: If set, will allow overwriting the file at the given fileSrc.
 *      If allowOverwite is false, will return error code AwFmFileAlreadyExists.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the write. Possible returns are:
 *      AwFmFileWriteOkay on success.
 *      AwFmFileAlreadyExists if a file exists at the given fileSrc, but allowOverwite was false.
 *      AwFmFileWriteFail if a file write failed.
 */
enum AwFmReturnCode awFmWriteIndexToFile(struct AwFmIndex *restrict const index,
  const uint64_t *restrict const suffixArray, const uint8_t *restrict const sequence,
  const uint64_t sequenceLength, const char *restrict const fileSrc, const bool allowOverwrite);


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
 *      AwFmFileFormatError if the header was not correct. This suggests that the file at this location
 *        is not the correct format.
 *      AwFmAllocationFailure on failure to allocated the necessary memory for the index.
 */
enum AwFmReturnCode awFmReadIndexFromFile(struct AwFmIndex *restrict *restrict index,
  const char *fileSrc, const bool keepSuffixArrayInMemory);



/*
 * Function:  awFmCreateKmerSearchList
 * --------------------
 *  Allocates and initializes an AwFmKmerSearchList struct to be used to search
 *  for groups of kmers in a thread-parallel manner that also hides memory read latency
 *  through multiple concurrent queries per thread.
 *  This struct can be used for both count and locate functions.
 *
 *  Note that the kmers inside the searchList are not allocated, and only contain
 *  char pointers that can be set to the kmers you want to query for.
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
 *  Deallocates the given search list struct, and both the internal searchData list
 *    and positionLists inside the searchData.
 *
 *  Note that, since the searchData doesn't own the kmer char strings, it will not try to
 *  deallocate them. If kmers were dynamically allocated externally, it is the caller's responsibility
 *  to deallocate them. This means that if these pointers are the only pointer to the data,
 *  and if they were dynamically allocated, forgetting to deallocate them before calling this function
 *  will leak the data.
 *
 *  Inputs:
 *    searchData:   pointer to the searchData struct to deallocate
 */
void awFmDeallocKmerSearchList(struct AwFmKmerSearchList *restrict const searchList);


/*
 * Function:  awFmParallelSearchLocate
 * --------------------
 *  Using the given index and a searchData struct preloaded with kmers, query the kmers
 *    in a concurrent, thread-parallel manner to find the positions of each occurrence of each
 *    respective kmer. The suggested use case for this function is as follows:
 *
 *    1. Allocate a searchList struct with awFmCreateKmerSearchList().
 *    2. For searching for n kmers, set count to n (must be smaller than capacity!), and
 *      set the first n kmers with the correct char pointer and length.
 *    3. Call this function using the index to search. The positions of each kmer will be loaded
 *      into the positionList member variable of the corresponding AwFmKmerSearchData struct.
 *    4. To query for additional kmers, reuse the searchList struct, starting with step (2).
 *    5. Deallocate the searchData with awFmDeallocKmerSearchList when finished.
 *
 *  Inputs:
 *    index:        pointer to the index to search.
 *    searchList:   pointer to the searchList struct loaded with kmers to search for.
 *    numThreads:   How many threads to direct OpenMP to use. The best value for this argument
 *                    will likely vary from system to system. Suggested default value is 4
 */
void awFmParallelSearchLocate(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, uint8_t numThreads);

/*
 * Function:  awFmParallelSearchCount
 * --------------------
 *  Using the given index and a searchData struct preloaded with kmers, query the kmers
 *    in a concurrent, thread-parallel manner to find the count of occurrences of each
 *    respective kmer. The suggested use case for this function is as follows:
 *
 *    1. Allocate a searchList struct with awFmCreateKmerSearchList().
 *    2. For searching for n kmers, set count to n (must be smaller than capacity!), and
 *      set the first n kmers with the correct char pointer and length.
 *    3. Call this function using the index to search. The count of matching kmers
 *      will be loaded into the AwFmKmerSearchData's count member variable.
 *    4. To query for additional kmers, reuse the searchList struct, starting with step (2).
 *    5. Deallocate the searchData with awFmDeallocKmerSearchList when finished.
 *
 *  Inputs:
 *    index:        pointer to the index to search.
 *    searchList:   pointer to the searchList struct loaded with kmers to search for.
 *    numThreads:   How many threads to direct OpenMP to use. The best value for this argument
 *                    will likely vary from system to system. Suggested default value is 4
 */
void awFmParallelSearchCount(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, uint8_t numThreads);



// //TODO: dreprecated
// /*
//  * Function:  awFmCreateParallelSearchData
//  * --------------------
//  *  Allocates and initializes an AwFmParallelSearchData struct to be used to search
//  *  for groups of kmers in a thread-parallel manner that also hides memory read latency
//  *  through multiple concurrent queries per thread.
//  *
//  *  Note that the kmers inside the searchData are not allocated, and only contain
//  *  char pointers that can be set to the kmers you want to query for.
//  *
//  *  Inputs:
//  *    capacity:     How many kmers the searchData struct can hold.
//  *    numThreads:   number of threads to use to query the searchData.
//  *
//  *  Returns:
//  *    Pointer to the allocated searchData struct, or null on failure.
//  */
// struct AwFmParallelSearchData *awFmCreateParallelSearchData(const size_t capacity,
//   const uint_fast8_t numThreads);


// //TODO: deprecated
// /*
//  * Function:  awFmDeallocParallelSearchData
//  * --------------------
//  *  Deallocates the given search data, and all internally stored dynamically allocated memory.
//  *
//  *  Note that, since the searchData doesn't own the kmer char strings, it will not try to
//  *  deallocate them. If kmers were dynamically allocated externally, it is the caller's responsibility
//  *  to deallocate them. This means that if these pointers are the only pointer to the data,
//  *  and if they were dynamically allocated, forgetting to deallocate them before calling this function
//  *  will leak the data.
//  *
//  *  Inputs:
//  *    searchData:   pointer to the searchData struct to deallocate
//  */
// void awFmDeallocParallelSearchData(struct AwFmParallelSearchData *restrict const searchData);

// //TODO: dreprecated
// /*
//  * Function:  awFmParallelSearch
//  * --------------------
//  *  Using the given index and a searchData struct preloaded with kmers, query the kmers
//  *  in a concurrent, thread-parallel manner. The suggested use case for this function is as follows:
//  *
//  *    1. allocate a searchData struct with awFmCreateParallelSearchData().
//  *    2. for searching for n kmers, set count to n (must be smaller than capacity!), and
//  *      set the first n kmers with the correct char pointer and length.
//  *    3. call this function using the index to search.
//  *    4. use the database sequence positions returned in the sequencePositionLists member.
//  *    5. to query for additional kmers, reuse the searchData struct, starting with step (2).
//  *    6. deallocate the searchData with awFmDeallocParallelSearchData when finished.
//  *
//  *
//  *  Inputs:
//  *    index:        pointer to the index to search.
//  *    searchData:   pointer to the searchData loaded with kmers to search for.
//  */
// void awFmParallelSearch(const struct AwFmIndex *restrict const index,
//   struct AwFmParallelSearchData *restrict const searchData);

/*
 * Function:  awFmReadSequenceFromFile
 * --------------------
 * Given a sequence position, reads a section of sequence surrounding that position from the
 *  corresponding index file.
 *
 *  Inputs:
 *    index:          Pointer to the AwFmIndex struct that contains the file handle to read.
 *      index file that stores the compressed suffix array.
 *    sequencePosition:        Position in the sequence to read.
 *    priorFlankLength: How many character before the given sequence position to include
 *      in the sequence segment.
 *    postFlankLength:  How many characters after the given sequence position to include
 *      in the sequence segement.
 *    sequenceBuffer: Pointer to the buffer to read the sequence segment into. This
 *      buffer  must be large enough to hold (postFlankLength + priorFlankLength +1) characters.
 *      the + 1 on this length is for the null terminator added to the end of the string.
 *
 *  Returns:
 *    AwFmReturnCode represnting the result of the read. Possible returns are:
 *      AwFmFileReadOkay on success.
 *      AwFmFileReadFail if the file could not be read sucessfully.
 */
enum AwFmReturnCode awFmReadSequenceFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankLength, const size_t postFlankLength,
  char *const sequenceBuffer);


#endif /* end of include guard: AW_FM_INDEX_STRUCTS_H */
