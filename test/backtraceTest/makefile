TEST_SRC = backtraceTest.c
SRC = $(wildcard ../../src/*.c)


CFLAGS = -std=c11 -Wall -mtune=native -fopenmp -mavx2 -O3
LDFLAGS = -L../../build
LDLIBS = -lfastavector_static -ldivsufsort64 -I../../build/

EXE = backtraceTest.out

backtraceTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS) $(LDFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(EXE)