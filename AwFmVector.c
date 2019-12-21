#include "AwFmVector.h"
#include <string.h>

#define CACHE_LINE_WIDTH 64


bool awFmVectorCreate(const size_t initialCapacity, const size_t elementSize, struct AwFmVector *restrict const vector){
  vector->capacity    = initialCapacity;
  vector->count       = 0;
  vector->elementSize = elementSize;
  vector->arrayPtr = aligned_alloc(CACHE_LINE_WIDTH, initialCapacity * elementSize);

  return vector->arrayPtr != NULL;
}


void awFmVectorDealloc(struct AwFmVector *const vector){
  free(vector->arrayPtr);
  vector->capacity  = 0;
  vector->count     = 0;
}


bool awFmVectorAppend(struct AwFmVector *const vector, const void *const data){
  if(vector->capacity == vector->count){
    const size_t oldLength = vector->capacity * vector->elementSize;
    vector->capacity = vector->capacity * 2;
    const size_t newLength = vector->capacity * vector->elementSize;

    //perform "aligned realloc"
    void *tmpPtr = aligned_alloc(CACHE_LINE_WIDTH, newLength);
    if(__builtin_expect(tmpPtr == NULL, 0)){
      return false;
    }

    memcpy(tmpPtr, vector->arrayPtr, oldLength);
    free(vector->arrayPtr);
    vector->arrayPtr = tmpPtr;
  }

  const size_t offset = vector->count * vector->elementSize;
  memcpy(vector->arrayPtr + offset, data, vector->elementSize);

  vector->count++;
  return true;
}
