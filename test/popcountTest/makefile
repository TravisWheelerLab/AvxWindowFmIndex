SRC    	= popcountTest.c
CFLAGS 	= -std=c11 -Wall -mtune=native -O3 -fopenmp -mavx2
EXE 		= popcountSpeedTest.out
popcountTest: $(SRC)
	gcc $(SRC) -o $(EXE) $(CFLAGS)
