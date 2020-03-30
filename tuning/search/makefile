TEST_SRC	= timeSearch.c
INCLUDE_DIR 	= $(DESTDIR)/usr/local/include
RPATH_DIR	= $(DESTDIR)/usr/local/lib
LIB_FLAGS 	= -ldivsufsort64 -lawfmindex
CFLAGS 		= -std=c11 -Wall -mtune=native -O3 -fopenmp  -mavx2 
EXE 		= timeSearch.out
TEST_OBJ	= timeSearch.o

.PHONY: all
all: $(TEST_OBJ)
	 gcc $(TEST_OBJ) -o $(EXE)  $(LIB_FLAGS) $(CFLAGS) -L $(RPATH_DIR) 

$(TEST_OBJ): $(TEST_SRC)
	 gcc $(TEST_SRC) -c $(CFLAGS) -I $(INCLUDE_DIR)

.PHONY: clean
clean:
	rm -f $(TEST_OBJ)
	rm -f $(EXE)
