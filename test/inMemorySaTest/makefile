TEST_SRC = inMemorySaTest.c
SRC = $(wildcard ../../src/*.c)

CFLAGS = -std=c11  -Wall -mtune=native -fopenmp -mavx2 -O3
LDLIBS = ../../build/libfastavector_static.a ../../build/libdivsufsort64.a -I../../build/

EXE = saInMemTest.out

backtraceTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(EXE)