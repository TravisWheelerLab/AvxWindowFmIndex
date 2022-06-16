TEST_SRC	= backtraceTest.c
SRC 			= $(wildcard ../../src/*.c)

CFLAGS 	= -std=c11 -fsanitize=address -Wall -mtune=native -g -fopenmp -ldivsufsort64 -lfastavector -mavx2 -O3 -L../../lib/libdivsufsort/build/lib
EXE 		= backtraceTest.out
backtraceTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
