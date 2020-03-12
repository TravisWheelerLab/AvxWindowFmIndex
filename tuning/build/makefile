TEST_SRC	= buildIndex.c
INCLUDE_DIR = ~/usr/local/lib
RPATH_DIR		= ~/usr/local/include
LIB_FLAGS = -ldivsufsort64 -lawfmindex
CFLAGS 	= -std=c11 -Wall -mtune=native -O3 -fopenmp  -mavx2 $(LIB_FLAGS) -I $(INCLUDE_DIR) -Wl,--rpath=$(RPATH_DIR)
EXE 		= buildIndex.out
bwtTest: $(TEST_SRC)
	gcc $(TEST_SRC) -o $(EXE) $(CFLAGS)
