# AvxWindowFmIndex (BETA)
A fast, AVX2 accelerated FM index library that utilizes windows of AVX2 vectors to quickly locate exact match Kmers in genetic data. This FM-index is highly optimized for both nucleotide and amino sequences.

This library is currently in beta release, so some bugs are expected. Please report any issues on the git project page.


## Prerequisites
The following is required to build and use this software
* GCC C compiler
* Make
* Cmake (for libdivsufsort submodule build process)


## Building on Linux
To build and the AwFmIndex shared library into the default install location, clone the repo, and cd into the project directory.
```shell
$ make
$ sudo make install
```

### Building to non-default location
Sometimes you may need to install the library to a non-default location, for example, if you do not have sudo privileges to write to /usr/local. To install the shared library into a non-default location, such as ~/usr/local:
```shell
$ destdir="~/usr/local"
$ make
$ make install
```

If AwFmIndex is installed to a non-default location, you may need to set environmental variables to allow your software to find the project .so and .h files at runtime.
```shell
$ export LD_LIBRARY_PATH=~/usr/local/lib
$ export LD_RUN_PATH=~/usr/local/include
```


The first time the project is built, it will attempt to clone the LibDivSufSort library into the project directory. As such, an internet connection is required during the initial call to make. This can be performed manually, if needed, but this is almost always unnecessary.
```shell
$ git submodule init
$ git submodule update
```


## AwFmIndex Quick Start Guide
The full public API can be found in the header src/AwFmIndex.h. Full examples can be found in the tuning/build and tuning/search directories.

### Making an .awfmi index file
To create an .awfmi index, use the awFmCreateIndex() function. This function will take a given sequence, and create an index file that conforms to the given metadata.

``` c
enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex *restrict *index,
  const struct AwFmIndexMetadata *restrict const metadata, const uint8_t *restrict const sequence,
  const size_t sequenceLength, const char *restrict const fileSrc, const bool allowFileOverwrite);
```

The metadata struct is as follows:
``` c
struct AwFmIndexMetadata{
  uint16_t              versionNumber;
  uint8_t               suffixArrayCompressionRatio;
  uint8_t               kmerLengthInSeedTable;
  enum AwFmAlphabetType alphabetType;
  bool                  keepSuffixArrayInMemory;
};
```

* version number represents the version of .awfmi file to create. At the moment, there is only one version, so a versionNumber of 1 should always be given.
* suffixArrayCompressionRatio represents how much to compress the suffix array to reduce the size of the .awfmi file on drive. As an example, a value of 8 will tell AwFmIndex to build a suffix array where 1/8th of the suffix array is sampled, and on average, each hit will take 8 additional backtrace operations to find the actual database sequence position. As the suffix array is never kept in memory queries, it will have no affect on memory usage during index searches.
* kmerLengthInSeedTable represents how long of kmers to memoize in a lookup table to speed up queries. Higher values will speed up searches, but will take exponentially more memory. A value of 12 (268MB lookup table) is recommended for nucleotide indices, and a value of 5 (51MB) is recommended for protein indices. increasing this value by one will result in 4x table size for nucleotide indices, and a 20x table size for protein indices.
* alphabetType allows the user to set the type of index to make. Options are AwFmAlphabetNucleotide and AwFmAlphabetAmino
* keepSuffixArrayInMemory determines if the compressed suffix array is loaded into memory, or left on disk. keeping the suffix array will consume a lot of memory (8 bytes per position in the database sequence), but will speed up searches by not having to go to disk for the final position lookup of each hit. An index made from an average mammalian nucleotide genome with this flag set to true will consume around 28GB of additional memory.


To use awFmCreateIndex, pass a pointer to an uninitialized AwFmIndex struct pointer. The function will allocate memory for the index, build it in memory, and write it to the given fileSrc. The AwFmIndex struct is usable immediately after calling this function, and must be manually deallocated with awFmDeallocIndex().


### Loading an existing Index
To load an existing .awfmi file, use the function
``` c
enum AwFmReturnCode awFmReadIndexFromFile(struct AwFmIndex *restrict *restrict index, const char *fileSrc,
  const bool keepSuffixArrayInMemory);
```

The index argument should be a pointer to an unallocated index pointer, just as with awFmCreateIndex(). Setting the keepSuffixArrayInMemory flag to true will load the compressed suffix array into memory along with the rest of the index.


### Querying batches of kmers in parallel
To query for batches of kmers, create a AwFmKmerSearchList struct with the function
``` c
struct AwFmKmerSearchList *awFmCreateKmerSearchList(const size_t capacity);
```

* capacity is the number of kmers the search data struct can hold.

Once the AwFmKmerSearchList struct had been allocated and initialized with the above function, load the kmer strings you wish to query for.
``` c
struct AwFmKmerSearchList *loadKmers(char **kmerStrings, uint32_t *kmerStringLengths, uint32_t numKmers){
  struct AwFmKmerSearchList * searchList = awFmCreateKmerSearchList(numKmers);
  for(size_t i = 0; i < numKmers; i++){
    searchList->kmerSearchData[i].kmerString = kmerStrings[i];
    searchList->kmerSearchData[i].kmerLength = kmerStringLengths[i];
  }
  searchList->count = numKmers;

  return searchData;
}
```

### Locate() function
To locate all instances of the kmers against an AwFmIndex, use the awFmParallelSearchLocate() function.
``` c
void awFmParallelSearchLocate(const struct AwFmIndex *restrict const index,
 struct AwFmKmerSearchList *restrict const searchList, uint8_t numThreads);
```

When the function returns, each kmer will have been queried, and a corresponding list of locations can be found in the corresponding entry in searchList->kmerSearchData[x].positionBacktraceList.position . The numThreads argument tells the algorithm how many threads to use. The optimal setting will vary from system to system, but 4 or 8 seem to be good choices on many systems.

e.g., to print the positions in the database sequence where a kmer at a given index was found:
``` c
void printKmerHitPositions(struct AwFmKmerSearchList *searchList, size_t kmerIndex)
  const uint32_t numHitPositions = searchList->kmerSearchData[kmerIndex].count;
  struct AwFmBacktrace *positionBacktraceList =  
    searchList->kmerSearchData[kmerIndex].positionBacktraceList;

  for(size_t i = 0; i < numHitPositions; i++){
    printf("kmer at index %zu found at database position %zu."\n,
    kmerIndex, positionBacktraceList[i].position);
  }
}
```

### Count() function
The count() function is used similary to the locate() function.

``` c
void awFmParallelSearchCount(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, uint8_t numThreads);
```

After querying for the count, the number of occurences of each kmer is located in the 'count' member variable of the corresponding AwFmKmerSearchData.

``` c
void printKmerHitCount(struct AwFmKmerSearchList *searchList, size_t kmerIndex)
  const uint32_t kmerCount = searchList->kmerSearchData[kmerIndex].count;
  printf("kmer at position %zu has %u occurrences.\n", kmerIndex, kmerCount);
}
```


### Deallocating the AwFmKmerSearchList
When finished using the AwFmKmerSearchList struct, deallocate it with the awFmDeallocKmerSearchList() function.
``` c
void awFmDeallocKmerSearchList(struct AwFmKmerSearchList *restrict const searchList);
```

### Reading back sections of the database sequence
If you would like to read sections of the database sequence around a given position, use the function:
``` c
enum AwFmReturnCode awFmReadSequenceFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequencePosition, const size_t priorFlankLength, const size_t postFlankLength,
  char *const sequenceBuffer);
```
where
* index is the AwFmIndex to query
* sequencePosition is the position to return a window around.
* priorFlankLength and postFlankLength are the number of characters to include before and after the sequencePosition, respectively
* sequenceBuffer is a preallocated buffer large enough to fit the window described by the flank lengths.

If the flank lengths would force the window to start before the start of the database sequence, or end after, the window is trimmed to stay within the bounds of the database sequence. Regardless, the sequenceBuffer will be null-terminated.
