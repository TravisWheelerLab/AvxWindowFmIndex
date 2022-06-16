TEST_SRC	= parallelSearchTest.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	= -std=c11 -fsanitize=address -Wall -mtune=native -g -fopenmp -ldivsufsort64 -lfastavector -mavx2 -L../../lib/libdivsufsort/build/lib
EXE 		= parallelSearchTest.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
