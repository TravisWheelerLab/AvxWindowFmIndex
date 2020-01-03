#ifndef AW_FM_VECTOR_H
#define AW_FM_VECTOR_H

#include <stdlib.h>
#include <stdbool.h>

struct AwFmBacktrace{
  size_t position;
  size_t offset;
};

struct AwFmBacktraceVector{
  size_t capacity;
  size_t count;
  struct AwFmBacktrace *backtraceArray;
};

bool awFmBacktraceVectorCreate(struct AwFmBacktraceVector *restrict const vector);

void awFmBacktraceVectorDealloc(struct AwFmBacktraceVector *const vector);

bool awFmBacktraceVectorSetCount(struct AwFmBacktraceVector *const vector, const size_t newCount);

struct AwFmBacktrace *AwFmBacktraceVectorBacktraceAtIndex(const struct AwFmBacktraceVector *restrict const vector,
  const size_t index);

#endif /* end of include guard: AW_FM_VECTOR_H */
