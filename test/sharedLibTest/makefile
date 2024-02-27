TEST_SRC = awfmiTest.c

CFLAGS = -std=c11 -Wall -mtune=native -fopenmp -mavx2 -O3
LDLIBS = -lawfmindex -I../../build/

EXE = sharedLibTest.out

sharedLibTest: $(TEST_SRC)
	if [ -z "$$LD_LIBRARY_PATH" ]; then export LD_LIBRARY_PATH=/usr/local/lib; fi
	gcc $(TEST_SRC)  -o $(EXE) $(CFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(EXE)