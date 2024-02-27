TEST_SRC = AwFmLetterTest.c
SRC = $(wildcard ../../src/*.c)

CFLAGS = -std=c11 -Wall -mtune=native -fopenmp -mavx2 -O3
LDLIBS = ../../lib/FastaVector/build/libfastavector_static.a ../../lib/libdivsufsort/build/lib/libdivsufsort64.a -I../../build/

EXE = letterTest.out

backtraceTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(EXE)