TEST_SRC	= saTest.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	=  -O0 -std=c11 -fsanitize=address -Wall -mtune=native -g -fopenmp -ldivsufsort64 -lfastavector -mavx2 -I../../libdivsufsort/build/include
EXE 		= saTest.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
