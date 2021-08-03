TEST_SRC	= saTest.c
SRC 			= $(wildcard ../../src/AwFmSuffixArray.c)
CFLAGS 	= -std=c11 -fsanitize=address -Wall -mtune=native -g -fopenmp  -mavx2 -I../../libdivsufsort/build/include
EXE 		= saTest.out
bwtTest: $(SRC)
	gcc $(TEST_SRC) -O0 $(SRC) -o $(EXE) $(CFLAGS)
