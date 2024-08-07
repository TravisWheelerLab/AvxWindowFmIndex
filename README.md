# AvxWindowFmIndex

A fast, SIMD accelerated FM-index library that utilizes windows of SIMD
registers to quickly locate exact match kmers in genetic data. This FM-index is
highly optimized for both nucleotide and amino sequences, but is unsuitable for
general text. Despite then name, AvxWindowFmIndex supports SIMD operation on
both x86_64 architectures via AVX2, and ARM_64 architectures via Arm Neon.

Please report any issues or bugs to the "Issues" tab on the project's github
(https://github.com/TravisWheelerLab/AvxWindowFmIndex), and view our published
article here (https://doi.org/10.1186/s13015-021-00204-6).

## Prerequisites

The following is required to build and use this software

* GCC C compiler
* CMake, or optionally,
* Make for the legacy build system


Note that MacOS users will need to install GCC through Homebrew or some other
means, because the built-in Clang compiler, which is aliased as `gcc`, does not
support OpenMP.

```
brew install gcc
```


After cloning the repo, you should initialize and update the submodules. 
```
git clone https://github.com/TravisWheelerLab/AvxWindowFmIndex.git
git submodule update --init --recursive --remote
```

Now the project is ready to build.

## CMake Build

The default way to build AwFmIndex is to use CMake. The usual CMake incantations
will work. This will result in both shared libraries written to the
`build/` directory.

```
cmake .
make
```

to instead make static libraries, set the BUILD_SHARED_LIBS variable as follows:
```
cmake -DBUILD_SHARED_LIBS=OFF .
make
```


If you want both, it's eas easy as running both back-to-back
```
cmake -DBUILD_SHARED_LIBS=OFF .
make
cmake -DBUILD_SHARED_LIBS=ON .
make
```

To install on MacOS, you will need to specify the install location of your
GCC compiler. An Example is given below:

```
cmake -DCMAKE_C_COMPILER=/path/to/gcc .
make
```

The library (both shared and static) is installed to the default install path.
If you'd like to change the installation location, change the PREFIX as follows:
`-DCMAKE_INSTALL_PREFIX=/prefix/path` to the `cmake` commands above.

```
make install
```

## Makefile Build (Legacy)

There is also a custom Makefile included, included to support use-cases where
CMake is disallowed.

### Shared Library

To build and install the AwFmIndex shared library into the default install
location:

```
make -f Makefile_legacy
```

To point the build at a specific version of GCC, which is necessary on a Mac,
use the following:

```
make -f Makefile_legacy CC=/path/to/gcc
```

### Static Library

To build a static library, just use the `static` target:

```
make -f Makefile_legacy static
```

This will generate two static libraries, `libawfmindex.a` and
`libdivsufsort64.a`, plus the associated header files in the build/ directory.

To point the build at a specific version of GCC, which is necessary on a Mac,
use the following:

```
make -f Makefile_legacy static CC=/path/to/gcc
```

### Install

After the library has been built, you can install both shared and static files
so that your compiler can find them. To install into the default location
(`/usr/local/`):

```
sudo -f Makefile_legacy make install
```

To specify a custom install location:

```
make -f Makefile_legacy install PREFIX=~/usr/local
```

If AwFmIndex is installed to a non-default location, you may need to set
environmental variables to allow your software to find the project `.so` and
`.h` files at runtime.

```
export LD_LIBRARY_PATH=~/.local/lib
export LD_RUN_PATH=~/.local/include
```

## AwFmIndex Quick Start Guide

The full public API can be found in the header src/AwFmIndex.h. Full examples
can be found in the tuning/build and tuning/search directories.

### Making an .awfmi index file

To create an .awfmi index, use the awFmCreateIndex() or awFmCreateIndexFromFasta()
functions. These functions will take a given sequence (or fasta containing one or more sequences),
and create an index file that conforms to the given configuration. These functions 
will overwrite the given fileSrc.

``` c
enum AwFmReturnCode awFmCreateIndex(struct AwFmIndex *restrict *index,
  struct AwFmIndexConfiguration *restrict const config, const uint8_t *restrict const sequence,
  const size_t sequenceLength, const char *restrict const fileSrc);
```

or, to generate an index from all sequences in a well-formed fasta file,

``` c
enum AwFmReturnCode awFmCreateIndexFromFasta(struct AwFmIndex *restrict *index,
  struct AwFmIndexConfiguration *restrict const config, const char *fastaSrc,
  const char *restrict const indexFileSrc);
```

Like all functions that return an 'enum AwFmReturnCode', make sure to check the
return code to determine if the result was successful or not. The function prototypes
in the header files contain robust documentation about functions, and their possible
return codes.

The configuration struct is as follows, the fields are described below:

``` c
struct AwFmIndexConfiguration{
  uint8_t               suffixArrayCompressionRatio;
  uint8_t               kmerLengthInSeedTable;
  enum AwFmAlphabetType alphabetType;
  bool                  keepSuffixArrayInMemory;
  bool                  storeOriginalSequence;
};
```

**`suffixArrayCompressionRatio`** represents how much to compress the suffix
*array to reduce the size of the
.awfmi file on drive. As an example, a value of 8 will tell AwFmIndex to build a
suffix array where 1/8th of the suffix array is sampled, and on average, each
hit will take 8 additional backtrace operations to find the actual database
sequence position. As the suffix array is never kept in memory queries, it will
have no affect on memory usage during index searches.

**`kmerLengthInSeedTable`** represents how long of kmers to memoize in a lookup
*table to speed up queries.
Higher values will speed up searches, but will take exponentially more memory. A
value of 12 (268MB lookup table) is recommended for nucleotide indices, and a
value of 5 (51MB) is recommended for protein indices. increasing this value by
one will result in 4x table size for nucleotide indices, and a 20x table size
for protein indices.

**`alphabetType`** allows the user to set the type of index to make. Options are
AwFmAlphabetDna, AwFmAlphabetRna, and AwFmAlphabetAmino

**`keepSuffixArrayInMemory`** determines if the compressed suffix array is
*loaded into memory, or left on
disk. keeping the suffix array will consume a lot of memory (8 bytes per
position in the database sequence), but will speed up searches by not having to
go to disk for the final position lookup of each hit. An index made from an
average mammalian nucleotide genome with this flag set to true will consume
around 28GB of additional memory.

**`storeOriginalSequence`** determines if the original sequence data will be
*saved inside the index. If this
is false, the sequence is omitted, which will generate a smaller index file. If
true, sections of the original sequence can be recalled with the
awFmReadSequenceFromFile() function.

To use `awFmCreateIndex` or `awFmCreateIndexFromFasta`, pass a pointer to an
uninitialized `AwFmIndex` struct. The function will allocate memory for the
index, build it in memory, and write it to the given `fileSrc`. The `AwFmIndex`
struct is usable immediately after calling this function, and must be manually
deallocated with `awFmDeallocIndex`.

### Loading an existing Index

To load an existing .awfmi file, use the function

``` c
enum AwFmReturnCode awFmReadIndexFromFile(struct AwFmIndex *restrict *restrict index, const char *fileSrc,
  const bool keepSuffixArrayInMemory);
```

The index argument should be a pointer to an unallocated index pointer, just as
with awFmCreateIndex(). Setting the keepSuffixArrayInMemory flag to true will
load the compressed suffix array into memory along with the rest of the index.


### Querying batches of kmers in parallel

To query for batches of kmers, create a `AwFmKmerSearchList` struct with the
function:

``` c
struct AwFmKmerSearchList *awFmCreateKmerSearchList(const size_t capacity);
```

**`capacity`** is the number of kmers the search data struct can hold.

Once the `AwFmKmerSearchList` struct had been allocated and initialized with the
above function, load the kmer strings you wish to query for. Example code to do
this is given below.

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

### Locate function

To locate all instances of the kmers against an AwFmIndex, use the
`awFmParallelSearchLocate` function.

``` c
enum AwFmReturnCode awFmParallelSearchLocate(const struct AwFmIndex *restrict const index,
 struct AwFmKmerSearchList *restrict const searchList, uint8_t numThreads);
```

When the function returns, each kmer will have been queried, and a corresponding
list of locations can be found in the corresponding entry in
`searchList->kmerSearchData[x].positionBacktraceList.position`. The `numThreads`
argument tells the algorithm how many threads to use. The optimal setting will
vary from system to system, but 4 or 8 seem to be good choices on many systems.
In the unlikely event of a failure to read from disk, this function will return AwFmFileReadFail.
Otherwise, it will return AwFmSuccess.

To print the positions in the database sequence where a kmer at a given index
was found:

``` c
void printKmerHitPositions(struct AwFmKmerSearchList *searchList, size_t kmerIndex)
  const uint32_t numHitPositions = searchList->kmerSearchData[kmerIndex].count;
  uint64_t *positionList =  
    searchList->kmerSearchData[kmerIndex].positionList;

  for(size_t i = 0; i < numHitPositions; i++){
    printf("kmer at index %zu found at database position %zu."\n,
    kmerIndex, positionList[i]);
  }
}
```

### Count function

The count function is used similary to the locate function.

``` c
void awFmParallelSearchCount(const struct AwFmIndex *restrict const index,
  struct AwFmKmerSearchList *restrict const searchList, uint8_t numThreads);
```

After querying for the count, the number of occurences of each kmer is located
in the `count` member of the corresponding `AwFmKmerSearchData` struct. Example
code is found below.

``` c
void printAllKmerCounts(struct AwFmKmerSearchList *searchList){
  const size_t numKmers = searchList->count;
  for(size_t i = 0; i < numKmers; i++){
    uint32_t thisLmerCount = searchList->kmerSearchData[i].count;
    printf("kmer at index %zu has %u counts\n");
  }
}
```

### Deallocating the AwFmKmerSearchList

When finished using the `AwFmKmerSearchList` struct, deallocate it with the
`awFmDeallocKmerSearchList` function.

``` c
void awFmDeallocKmerSearchList(struct AwFmKmerSearchList *restrict const searchList);
```

### Reading back sections of the database sequence

If you would like to read sections of the database sequence around a given
position, use the function:

``` c
enum AwFmReturnCode awFmReadSequenceFromFile(const struct AwFmIndex *restrict const index,
  const size_t sequenceStartPosition, const size_t sequenceSegmentLength,
  char *const sequenceBuffer);
```

**`index`** is the `AwFmIndex` to query.

**`sequenceStartPosition`** is the position of the first character to be
included in the window.

**`sequenceSegmentLength`** is the length of the sequence segment to read.

**`sequenceBuffer`** is a preallocated buffer large enough to fit sequence
segment. Once populated, the buffer is null terminated.

The total number of characters read from the file equals `sequenceEndPosition -
sequenceStartPosition`. Giving a `sequenceEndPosition` greater than the length
of the sequence can result in undefined behavior. This function will return the
error code `AwFmUnsupportedVersionError` if the index did not store the sequence
data (i.e., the configuration's `storeOriginalSequence` member variable was set
to false when the index was generated.)


### Optional functionality

Depending on the user-set configuration parameters, these functions may apply.

If a `AwFmIndex` is built from a fasta file, it will keep track of the sequence
lengths, and is able to determine which sequence a hit is located in, and the
position inside that sequence.

``` c
enum AwFmReturnCode awFmGetLocalSequencePositionFromIndexPosition(const struct AwFmIndex *restrict const index,
  size_t globalPosition, size_t *sequenceNumber, size_t *localSequencePosition);
```

**`index`** is the AwFmIndex used.

**`globalPosition`** is the position returned as the result of a locate() search.

**`sequenceNumber`** is an out-argument where the index of the sequence the hit
falls inside will be written.

**`localSequencePosition`** is an out-argument where the position in that
sequence that corresponds to the globalPosition will be written.

If this function is called on an index that wasn't built from a fasta file, it
will return the error code AwFmUnsupportedVersionError.

If an index is built from a fasta, the sequence headers can also be retrieved.

``` c
enum AwFmReturnCode awFmGetHeaderStringFromSequenceNumber(const struct AwFmIndex *restrict const index,
  size_t sequenceNumber, char **headerBuffer, size_t *headerLength);
```

**`index`** is the AwFmIndex used,

**`globalPosition`** is the position returned as the result of a locate() search.

**`sequenceNumber`** is the index of the relevant sequence.

**`headerBuffer`** is a pointer to the char* variable that this function will
set to the beginning of the header.

**`headerLength`** is an out-argument where the function will write the length
of the header.

This function also returns `AwFmUnsupportedVersionError` if the index was not
built from a fasta file.




### Usage Example

Here is an example of creating an AwFmIndex and using it to perform a search on queries given
via the argument list parameter.
If you already have an index you'd like to use instead of generating a new one, use awFmReadIndexFromFile().
``` c
int main(int argc, char **argv){

  char *indexFileSrc = "indexFiles/index.awfmi";
  char *fastaInputFileSrc = "fastas/dnaSequence.fasta"
	struct AwFmIndex *index;
	struct AwFmIndexConfiguration config = {.suffixArrayCompressionRatio = 8,
			.kmerLengthInSeedTable																					 = 12,
			.alphabetType																										 = AwFmAlphabetNucleotide,
			.keepSuffixArrayInMemory																				 = false,
			.storeOriginalSequence																					 = true};
	enum AwFmReturnCode returnCode = awFmCreateIndexFromFasta(&index, &config, sequence, sequenceLength, indexFileSrc, true);
  if(awFmReturnCodeIsFailure(returnCode)){
    printf("create index failed with return code %u\n", returnCode);
    exit(1);
  }

  uint64_t numQueries = argc-1; //first argument is the name of the executable, so skip it.
	struct AwFmKmerSearchList *searchList = awFmCreateKmerSearchList(numQueries);
  if(searchList == NULL){
    printf("could not allocate memory for search list\n");
    exit(2);
  }

  searchList->count = numQueries;
  //initialize the queries in the search list
  for(size_t i = 0; i < numQueries; i++){
    searchList->kmerSearchData[i].kmerLength = strlen(argv[i+1]);
    searchList->kmerSearchData[i].kmerString = &argv[i+1];
  }

  //search for the queries

  returnCode awFmParallelSearchLocate(index, searchList, numThreads);
  if(awFmReturnCodeIsFailure(returnCode)){
    printf("parallel search failed with return code %u\n", returnCode);
    exit(3);
  }

  //print all the hits.
  for(size_t kmerIndex = 0; kmerIndex < searchList->count; kmerIndex++){
    const uint32_t numHitPositions = searchList->kmerSearchData[kmerIndex].count;
    uint64_t *positionList =  
      searchList->kmerSearchData[kmerIndex].positionList;

    for(size_t i = 0; i < numHitPositions; i++){
      printf("kmer at index %zu found at database position %zu."\n,
      kmerIndex, positionList[i]);
    }
  }

  //code cleanup
  awFmDeallocKmerSearchList(searchList);
  awFmDeallocIndex(index);

  exit(0);
}
```
