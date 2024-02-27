TEST_SRC = AwFmMultiSequenceTest.c
SRC = $(wildcard ../../src/*.c) $(wildcard ../../lib/FastaVector/src/*.c)

CFLAGS = -std=c11  -Wall -mtune=native -fopenmp -mavx2 -O0 -g
LDLIBS = ../../build/libdivsufsort64.a -I../../build/

EXE = multiSequenceTest.out

backtraceTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(EXE)