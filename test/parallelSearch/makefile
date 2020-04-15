TEST_SRC	= parallelSearchTest.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	= -std=c11 -Wall -mtune=native -g -fopenmp -ldivsufsort64 -mavx2
EXE 		= parallelSearchTest.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
