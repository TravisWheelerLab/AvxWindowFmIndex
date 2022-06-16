TEST_SRC	= saTest.c
SRC 			= $(wildcard ../../src/*.c)
CFLAGS 	=  -O0 -std=c11 -fsanitize=address -Wall -mtune=native -g -fopenmp -ldivsufsort64 -lfastavector -mavx2 -L../../lib/libdivsufsort/build/lib
bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS)
