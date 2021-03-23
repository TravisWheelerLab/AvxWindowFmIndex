TEST_SRC	= searchTest.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	= -std=c11 -fsanitize=address -Wall -mtune=native -g -fopenmp -ldivsufsort64 -lfastavector -mavx2 -I../../libdivsufsort/build/include
EXE 		= searchTest.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
