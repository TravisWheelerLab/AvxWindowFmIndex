TEST_SRC = awfmiTest.c

CFLAGS = -std=c11 -Wall -mtune=native -fopenmp -mavx2 -O3
LDFLAGS = -L../../build -I../../build/
LDLIBS = -lawfmindex

EXE = staticLibTest.out

staticLibTest: $(TEST_SRC)
	gcc $(TEST_SRC) $(LDFLAGS)  -o $(EXE) $(CFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(EXE)