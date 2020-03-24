#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include "AwFmIndex.h"

char *const fullGenomeFile = "chr.all.fa";

int64_t getFullGenomeFromFasta(char *fastaFileSrc, char *buffer);
int64_t getChromosomeFromFullFasta(char *fastaFileSrc, int chromosomeNumber, char *buffer);
void buildChromosomeIndex(int chromosomeNumber, uint8_t suffixArrayCompressionRatio, uint8_t kmerLengthInSeedTable, char *indexFilename);
void buildFullGenomeIndex(uint8_t suffixArrayCompressionRatio, uint8_t kmerLengthInSeedTable, char *indexFilename);
void parseArgs(int argc, char **argv);

//parameters

bool isAminoSequence;
int chromosomeNumber;
bool chromosomeNumberSetInArgs;
bool fullGenomeIndexRequested;
int suffixArrayCompressionRatio;
int kmerLengthInSeedTable;
bool kmerLengthSetInArgs;
char indexFilenameBuffer[256];


int main(int argc, char **argv){
  //init parameters
  strcpy(indexFilenameBuffer, "index.awfmi");
  bool isAminoSequence            = false;
  int chromosomeNumber            = -1;
  bool fullGenomeIndexRequested   = false;
  int suffixArrayCompressionRatio = 8;
  int kmerLengthInSeedTable       = 12;
  bool kmerLengthSetInArgs        = false;
  parseArgs(argc, argv);

  //set the kmer length in seed table based on the alphabet type, if not specified.
  if(!kmerLengthSetInArgs){
    kmerLengthInSeedTable = isAminoSequence? 12: 5;
  }


  if(chromosomeNumber > 22){
    printf("Error: chromosome number was larget than 21. Things might get weird.\n");
  }
  if(chromosomeNumber == -1 && !fullGenomeIndexRequested){
    printf("no chromosome or full genome specified, defaulting to chromosome 21.\n");
    chromosomeNumber = 21;
  }else if(chromosomeNumber != -1 && fullGenomeIndexRequested){
    printf("the full genome flag was set, but a specific chromosome was also requested. Ignoring full genome and making index on chromosome %i\n", chromosomeNumber);
    fullGenomeIndexRequested = false;
  }


  //build the index using the parsed settings
  if(fullGenomeIndexRequested){
    buildFullGenomeIndex(suffixArrayCompressionRatio, kmerLengthInSeedTable, indexFilenameBuffer);
  }
  else{
    buildChromosomeIndex(chromosomeNumber, suffixArrayCompressionRatio,kmerLengthInSeedTable, indexFilenameBuffer);
  }


}


void parseArgs(int argc, char **argv){
  int option = 0;
  while((option = getopt(argc, argv, "ac:gs:k:f:")) != -1){
    switch(option){
      case 'a':
        isAminoSequence = true;
        break;
      case 'c':
        sscanf(optarg, "%i", &chromosomeNumber);
        chromosomeNumberSetInArgs = true;
        break;
      case 'g':
        fullGenomeIndexRequested = true;
        break;
      case 's':
        sscanf(optarg, "%i", &suffixArrayCompressionRatio);
        break;
      case 'k':
        sscanf(optarg, "%i", &kmerLengthInSeedTable);
        kmerLengthSetInArgs = true;
      case 'f':
        strcpy(indexFilenameBuffer, optarg);
        break;
    }
  }
}



void buildChromosomeIndex(int chromosomeNumber, uint8_t suffixArrayCompressionRatio, uint8_t kmerLengthInSeedTable, char *indexFilename){
  char *sequenceBuffer = malloc(1<< 30* sizeof(char));
  int64_t sequenceLength = getChromosomeFromFullFasta(fullGenomeFile, chromosomeNumber, sequenceBuffer);
  if(sequenceLength == 0){
    printf("ERROR: could not get chromosome %i\n", chromosomeNumber);
    exit(-10);
  }

  struct AwFmIndex *index;
  struct AwFmIndexMetadata metadata = {.versionNumber=1, .suffixArrayCompressionRatio=suffixArrayCompressionRatio,
    .kmerLengthInSeedTable=kmerLengthInSeedTable, .alphabetType=AwFmAlphabetNucleotide};
  enum AwFmReturnCode returnCode = awFmCreateIndex(&index, &metadata, (uint8_t*)sequenceBuffer, sequenceLength,
    indexFilename, true);

  if(returnCode < 0){
    printf("Error: awFmCreateIndex returned error code %i", returnCode);
    exit(-1);
  }

  awFmDeallocIndex(index);
  free(sequenceBuffer);
}


void buildFullGenomeIndex(uint8_t suffixArrayCompressionRatio, uint8_t kmerLengthInSeedTable, char *indexFilename){
  char *sequenceBuffer = malloc(1ULL<< 32ULL* sizeof(char));
  int64_t sequenceLength = getFullGenomeFromFasta(fullGenomeFile, sequenceBuffer);
  if(sequenceLength == 0){
    printf("ERROR: could not load full genome\n");
    exit(-999);
  }

  struct AwFmIndex *index;
  struct AwFmIndexMetadata metadata = {.versionNumber=1, .suffixArrayCompressionRatio=suffixArrayCompressionRatio,
    .kmerLengthInSeedTable=kmerLengthInSeedTable, .alphabetType=AwFmAlphabetNucleotide};
  enum AwFmReturnCode returnCode = awFmCreateIndex(&index, &metadata, (uint8_t*)sequenceBuffer, sequenceLength,
    indexFilename, true);

  if(returnCode < 0){
    printf("Error: awFmCreateIndex returned error code %i", returnCode);
    exit(-1);
  }

  awFmDeallocIndex(index);
  free(sequenceBuffer);
}

int64_t getChromosomeFromFullFasta(char *fastaFileSrc, int chromosomeNumber, char *buffer){

  const uint8_t chromosomeNumberOnesPlace = chromosomeNumber % 10;
  const uint8_t chromosomeNumberTensPlace = chromosomeNumber / 10;

  FILE *sequenceFile = fopen(fastaFileSrc, "r");
  if(sequenceFile == NULL){
    printf("ERROR: could not read sequence fasta file at location %s\n", fastaFileSrc);
    exit(-1000);
  }
  size_t sequenceLength = 0;
  bool readingSelectedChromosome = false;

  while(!feof(sequenceFile)){

    //read a line into the buffer.
    char *fgetsResult = fgets(buffer+sequenceLength, 1024, sequenceFile);

    if(fgetsResult == NULL){
      printf("CRITICAL ERROR: fgets returned null. terminating.");
      fclose(sequenceFile);
      exit(-1);
    }
    // check to see if it's a header
    bool isHeader = (buffer[sequenceLength]) == '>';

    //if it was a header, check to see if it's the one we care about
    if(isHeader){
      bool thisHeaderIsCorrectChromosomeHeader =
        (chromosomeNumber <= 9 && buffer[sequenceLength+4] == chromosomeNumberOnesPlace) ||
        (chromosomeNumber > 9 && buffer[sequenceLength+4] == chromosomeNumberTensPlace && buffer[sequenceLength+5] == chromosomeNumberOnesPlace);
        //if we just were reading the correct chromosome, and just hit a header for a new chromosome, we're done.
      if(readingSelectedChromosome && !thisHeaderIsCorrectChromosomeHeader){
        fclose(sequenceFile);
        return sequenceLength;

      }
      readingSelectedChromosome = thisHeaderIsCorrectChromosomeHeader;
    }

    //if we're reading from the correct chromosome, move the buffer offset to preserve the read data.
    if(!isHeader && readingSelectedChromosome){
      sequenceLength += strlen(buffer+sequenceLength);
    }

  }
  fclose(sequenceFile);
  return sequenceLength;
}


int64_t getFullGenomeFromFasta(char *fastaFileSrc, char *buffer){
  FILE *sequenceFile = fopen(fastaFileSrc, "r");
  if(sequenceFile == NULL){
    printf("ERROR: could not read sequence fasta file at location %s\n", fastaFileSrc);
    exit(-2000);
  }
  size_t sequenceLength = 0;
  while(!feof(sequenceFile)){

    //read a line into the buffer
    fgets(buffer+sequenceLength, 1024, sequenceFile);

    // if the line wasn't a header, include it in the buffer.
    //otherwise, don't change bufferOffset and allow it to be overwritten.
    if(buffer[sequenceLength] != '<'){
      sequenceLength += strlen(buffer+sequenceLength);
    }
  }

  fclose(sequenceFile);
  return sequenceLength;
}
