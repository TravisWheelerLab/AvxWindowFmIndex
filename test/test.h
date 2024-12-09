#ifndef AW_FM_INDEX_TEST_H
#define AW_FM_INDEX_TEST_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define testAssert(assertion)                                                  \
  (awFmAssertTest(assertion, __FILE__, __func__, __LINE__))
#define testAssertString(assertion, message)                                   \
  (awFmAssertTestString(assertion, message, __FILE__, __func__, __LINE__))
inline void setPrintfColorRed() { printf("\033[1;31m"); }
inline void setPrintfColorDefault() { printf("\033[0m"); }

size_t assertionNumber = 0;
inline void awFmAssertTest(bool assertion, const char *file, const char *func,
                    const int line) {
  if (!assertion) {
    setPrintfColorRed();
    printf("assertion failure at test %zu, file %s, func %s, line %d\n",
           assertionNumber, file, func, line);
    setPrintfColorDefault();
  }
  assertionNumber++;
}

inline void awFmAssertTestString(bool assertion, char *message, const char *file,
                          const char *func, const int line) {
  if (!assertion) {
    setPrintfColorRed();
    printf("assertion failure at test %zu, file %s, func %s, line %d\t %s\n",
           assertionNumber, file, func, line, message);
    setPrintfColorDefault();
  }
  assertionNumber++;
}

inline void resetAssertionNumber() { assertionNumber = 0; }

#endif /* end of include guard: AW_FM_INDEX_TEST_H */
