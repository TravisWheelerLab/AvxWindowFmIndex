SRC    	= occurrenceTests.c ../../AwFmBacktraceVector.c ../../AwFmCreate.c ../../AwFmFile.c ../../AwFmIndex.c ../../AwFmKmerTable.c ../../AwFmLetter.c ../../AwFmOccurrence.c ../../AwFmParallelSearch.c ../../AwFmSearch.c
CFLAGS 	= -std=c11 -Wall -mtune=native -O0 -fopenmp -ldivsufsort64 -mavx2
EXE 		= occurrenceTest.out
letterTest: $(SRC)
	gcc $(SRC) -o $(EXE) $(CFLAGS)
