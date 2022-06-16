TEST_SRC	= AwFmMultiSequenceTest.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	= -std=c11 -fsanitize=address -Wall -mtune=native -g -lfastavector -ldivsufsort64 -fopenmp -mavx2 -L../../lib/libdivsufsort/build/lib
EXE 		= multiSequenceTest.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
