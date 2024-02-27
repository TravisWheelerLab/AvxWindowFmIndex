TEST_SRC	= searchTest.c
SRC 			= $(wildcard ../../src/*.c)

CFLAGS = -std=c11 -Wall -mtune=native -fopenmp -mavx2 -O0 -g
LDLIBS = ../../build/libfastavector_static.a ../../build/libdivsufsort64.a -I../../build/

EXE 		= searchTest.out

bwtTest: $(SRC)
	gcc $(TEST_SRC) $(SRC) -o $(EXE) $(CFLAGS) $(LDLIBS)
