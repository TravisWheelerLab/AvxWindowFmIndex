#ifndef AW_FM_VECTOR_H
#define AW_FM_VECTOR_H

#include <stdlib.h>
#include <stdbool.h>
#include "AwFmIndex.h"

bool awFmBacktraceVectorCreate(struct AwFmBacktraceVector *restrict const vector);

void awFmBacktraceVectorDealloc(struct AwFmBacktraceVector *const vector);

bool awFmBacktraceVectorSetCount(struct AwFmBacktraceVector *const vector, const size_t newCount);

struct AwFmBacktrace *AwFmBacktraceVectorBacktraceAtIndex(const struct AwFmBacktraceVector *restrict const vector,
  const size_t index);

#endif /* end of include guard: AW_FM_VECTOR_H */
