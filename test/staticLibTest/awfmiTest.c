#include "AwFmIndex.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {

  struct AwFmIndex *index;

  struct AwFmIndexConfiguration config;
  config.suffixArrayCompressionRatio = 2;
  config.kmerLengthInSeedTable = 2;
  config.alphabetType = AwFmAlphabetDna;
  config.keepSuffixArrayInMemory = true;
  config.storeOriginalSequence = false;

  enum AwFmReturnCode awfmrc =
      awFmCreateIndexFromFasta(&index, &config, "test.fa", "output.awfmi");

  if (awFmReturnCodeIsFailure(awfmrc)) {
    printf("ERROR: createIndex returned error code %u\n", awfmrc);
    exit(-1);
  } else {
    printf("createIndex successful\n");
  }

  awFmDeallocIndex(index);
}