TEST_SRC	= timeSearch.c
INCLUDE_DIR = ~/usr/local/lib
RPATH_DIR		= ~/usr/local/include
CFLAGS 	= -std=c11 -Wall -mtune=native -O3 -fopenmp -ldivsufsort64 -mavx2 -I $(INCLUDE_DIR) -Wl,--rpath=$(RPATH_DIR)
EXE 		= timeSearch.out
bwtTest: $(TEST_SRC)
	gcc $(TEST_SRC) -o $(EXE) $(CFLAGS)
