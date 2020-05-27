TEST_SRC	= bwtTest.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	= -std=c11 -Wall -fsanitize=address -mtune=native -g -fopenmp -ldivsufsort64 -mavx2
EXE 		= bwtTest.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
