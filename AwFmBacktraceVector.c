#include "AwFmBacktraceVector.h"
#include <string.h>
#include <stdint.h>

#define CACHE_LINE_WIDTH 64
#define DEFAULT_POSITION_LIST_CAPACITY  256

bool awFmBacktraceVectorCreate(struct AwFmBacktraceVector *restrict const vector){
  vector->capacity    = DEFAULT_POSITION_LIST_CAPACITY;
  vector->count       = 0;
  vector->backtraceArray = aligned_alloc(CACHE_LINE_WIDTH, DEFAULT_POSITION_LIST_CAPACITY * sizeof(struct AwFmBacktrace));

  return vector->backtraceArray != NULL;
}


void awFmBacktraceVectorDealloc(struct AwFmBacktraceVector *const vector){
  free(vector->backtraceArray);
  vector->capacity  = 0;
  vector->count     = 0;
}


bool awFmBacktraceVectorSetCount(struct AwFmBacktraceVector *const vector, const size_t newCount){
  if(__builtin_expect(vector->capacity >= newCount, 1)){
    vector->count = newCount;
    return true;

  }
  else{

    const size_t oldLength = vector->capacity * sizeof(struct AwFmBacktrace);
    const size_t newCapacity = newCount * 2;
    const size_t newLength = newCapacity * sizeof(struct AwFmBacktrace);
    //perform "aligned realloc"
    void *tmpPtr = aligned_alloc(CACHE_LINE_WIDTH, newLength);
    if(__builtin_expect(tmpPtr == NULL, 0)){
      return false;
    }
    memcpy(tmpPtr, vector->backtraceArray, oldLength);
    free(vector->backtraceArray);

    vector->count = newCount;
    vector->capacity = newCapacity;
    vector->backtraceArray = tmpPtr;
    return true;
  }
}


struct AwFmBacktrace *AwFmBacktraceVectorBacktraceAtIndex(const struct AwFmBacktraceVector *restrict const vector,
  const size_t index){
  return &(vector->backtraceArray[index]);
}
