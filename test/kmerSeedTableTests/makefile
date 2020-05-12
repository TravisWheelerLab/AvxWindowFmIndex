TEST_SRC	= kmerSeedTableTests.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	= -std=c11 -Wall -mtune=native -g -fopenmp -ldivsufsort64 -mavx2 -I../../libdivsufsort/build/include
EXE 		= kmerTableTests.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
