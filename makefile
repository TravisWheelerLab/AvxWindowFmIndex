occSpeedTest: AwFmApi.c AwFmApi.h AwFmCreate.c AwFmCreate.h AwFmFile.c AwFmFile.h AwFmGlobals.h AwFmIndex.c AwFmIndex.h AwFmLetter.c AwFmLetter.h AwFmOccupancy.c AwFmOccupancy.h AwFmSearch.c AwFmSearch.h test/occSpeedTest.c
	gcc AwFmIndex.c AwFmLetter.c AwFmOccupancy.c test/occSpeedTest.c -o occSpeedTest -std=c11 -Wall -march=native -O3

