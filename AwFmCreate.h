#ifndef AW_FM_INDEX_CREATE_H
#define AW_FM_INDEX_CREATE_H

#include "AwFmIndex.h"
#include <stdint.h>


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
enum AwFmReturnCode awFmCreateIndex(const struct AwFmIndex *restrict *index,
  const struct AwFmIndexMetadata *restrict const metadata, const uint8_t *restrict const sequence, const size_t sequenceLength,
  const char *restrict const fileSrc, const bool allowFileOverwrite);

#endif /* end of include guard: AW_FM_INDEX_CREATE_H */
