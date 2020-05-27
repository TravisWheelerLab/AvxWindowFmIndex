TEST_SRC	= occurrenceTests.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	= -std=c11 -fsanitize=address -Wall -mtune=native -g -fopenmp -ldivsufsort64 -mavx2 -I../../libdivsufsort/build/include
EXE 		= occurrenceTests.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
