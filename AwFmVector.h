#ifndef AW_FM_VECTOR_H
#define AW_FM_VECTOR_H

#include <stdlib.h>
#include <stdbool.h>

struct AwFmVector{
  size_t capacity;
  size_t count;
  size_t elementSize;
  void *arrayPtr;
};

bool awFmVectorCreate(const size_t initialCapacity, const size_t elementSize, struct AwFmVector *restrict const vector);

void awFmVectorDealloc(struct AwFmVector *const vector);

bool awFmVectorAppend(struct AwFmVector *const vector, const void *const data);

#endif /* end of include guard: AW_FM_VECTOR_H */
