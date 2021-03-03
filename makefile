#for libdivsufsort,  cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="../../lib" .. -DBUILD_DIVSUFSORT64:BOOL=ON

#when calling, allow for setting the DESTDIR env variable
#when calling, allow for the #DONT_INSTALL_DIVSUFSORT env variable
#todo: documentation for calling git clone with --recurse-submodules or use git submodule init & git submodule update
SHELL = /bin/sh

AWFMINDEX_PROJECT_NAME 			= awfmindex
AWFMINDEX_MAJOR_VERSION_NUM = 0
AWFMINDEX_MINOR_VERSION_NUM = 1
AWFMINDEX_VERSION 					= $(MAJOR).$(MINOR)

#filenames
AWFMINDEX_HEADER_FILENAME					= AwFmIndex.h
AWFMINDEX_SHARED_LIB_FILENAME			= libawfmindex.so
AWFMINDEX_STATIC_LIB_FILENAME			=	libawfmindex.a
LIBDIVSUFSORT_HEADER_FILENAME			=	divsufsort64.h
LIBDIVSUFSORT_STATIC_LIB_FILENAME = libdivsufsort64.a
FASTA_VECTOR_HEADER_FILENAME 			= fastaVector.h
FASTA_VECTOR_STATIC_LIB_FILENAME 	= fastaVector.a


#determine the current operating system and architecture.
OS_NAME 	= $(shell uname -s)
ARCH_NAME = $(shell uname -p)	#currently unused.
ifeq ($(OS_NAME), Darwin)
AWFMINDEX_SHARED_LIB_FILENAME = libawfmindex.dylib
endif

#directories
AWFMINDEX_PROJECT_DIR						:= $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
AWFMINDEX_BUILD_DIR 						= $(AWFMINDEX_PROJECT_DIR)/build
AWFMINDEX_SRC_DIR 							= $(AWFMINDEX_PROJECT_DIR)/src
AWFMINDEX_BUILD_INCLUDE_DIR			= $(AWFMINDEX_BUILD_DIR)/include
AWFMINDEX_BUILD_LIBRARY_DIR			= $(AWFMINDEX_BUILD_DIR)/lib
AWFMINDEX_INSTALL_INCLUDE_DIR		= $(DESTDIR)$(PREFIX)/include
AWFMINDEX_INSTALL_LIBRARY_DIR		= $(DESTDIR)$(PREFIX)/lib
LIBDIVSUFSORT_PROJECT_DIR				= $(AWFMINDEX_PROJECT_DIR)/libdivsufsort
LIBDIVSUFSORT_BUILD_DIR					= $(LIBDIVSUFSORT_PROJECT_DIR)/build
LIBDIVSUFSORT_BUILD_INCLUDE_DIR	= $(LIBDIVSUFSORT_BUILD_DIR)/include
LIBDIVSUFSORT_BUILD_LIBRARY_DIR = $(LIBDIVSUFSORT_BUILD_DIR)/lib
FASTA_VECTOR_PROJECT_DIR				= $(AWFMINDEX_PROJECT_DIR)/FastaVector
FASTA_VECTOR_SRC_DIR						= $(FASTA_VECTOR_PROJECT_DIR)/src
FASTA_VECTOR_BUILD_DIR					= $(FASTA_VECTOR_PROJECT_DIR)/build
FASTA_VECTOR_BUILD_INCLUDE_DIR	=	$(FASTA_VECTOR_BUILD_DIR)/include
FASTA_VECTOR_BUILD_LIB_DIR			= $(FASTA_VECTOR_BUILD_DIR)/lib



#file locations
AWFMINDEX_SRC_HEADER_FILE 								= $(AWFMINDEX_SRC_DIR)/$(AWFMINDEX_HEADER_FILENAME)
AWFMINDEX_BUILD_HEADER_FILE								= $(AWFMINDEX_BUILD_INCLUDE_DIR)/$(AWFMINDEX_HEADER_FILENAME)
AWFMINDEX_BUILD_SHARED_LIB_FILE						= $(AWFMINDEX_BUILD_LIBRARY_DIR)/$(AWFMINDEX_SHARED_LIB_FILENAME)
AWFMINDEX_BUILD_STATIC_LIB_FILE						= $(AWFMINDEX_BUILD_LIBRARY_DIR)/$(AWFMINDEX_STATIC_LIB_FILENAME)
AWFMINDEX_INSTALL_HEADER_FILE							= $(AWFMINDEX_INSTALL_INCLUDE_DIR)/$(AWFMINDEX_HEADER_FILENAME)
AWFMINDEX_INSTALL_STATIC_LIB_FILE					= $(AWFMINDEX_INSTALL_LIBRARY_DIR)/$(AWFMINDEX_STATIC_LIB_FILENAME)
AWFMINDEX_INSTALL_SHARED_LIB_FILE					= $(AWFMINDEX_INSTALL_LIBRARY_DIR)/$(AWFMINDEX_SHARED_LIB_FILENAME)
LIBDIVSUFSORT_BUILD_HEADER_FILE						= $(LIBDIVSUFSORT_BUILD_INCLUDE_DIR)/$(LIBDIVSUFSORT_HEADER_FILENAME)
LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE		= $(LIBDIVSUFSORT_BUILD_LIBRARY_DIR)/$(LIBDIVSUFSORT_STATIC_LIB_FILENAME)
FASTA_VECTOR_BUILD_HEADER_FILE						=	$(FASTA_VECTOR_SRC_DIR)/$(FASTA_VECTOR_HEADER_FILENAME)
FASTA_VECTOR_BUILD_STATIC_LIBRARY_FILE		= $(FASTA_VECTOR_BUILD_DIR)/$(FASTA_VECTOR_STATIC_LIB_FILENAME)

#if on a Mac system, CC may need to be overwritten with a makefile argument to compile with
#an actual GCC compiler instead of the clang that ships native with mac.
CC 			= gcc
CFLAGS 	= -std=gnu11 -fpic -O3 -mtune=native -Wall -Wextra -fopenmp



LDFLAGS 	= -shared -L$(LIBDIVSUFSORT_BUILD_LIBRARY_DIR) -I$(LIBDIVSUFSORT_BUILD_INCLUDE_DIR) -ldivsufsort64 # linking flags

SOURCE_FILES 	:= $(wildcard $(AWFMINDEX_SRC_DIR)/*.c)
OBJECT_FILES 	:= $(patsubst $(AWFMINDEX_SRC_DIR)/%, $(AWFMINDEX_BUILD_DIR)/%, $(SOURCE_FILES:.c=.o))

#rules
.PHONY: all
all: $(LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE) $(FASTA_VECTOR_BUILD_STATIC_LIBRARY_FILE) $(AWFMINDEX_BUILD_LIBRARY_DIR) $(OBJECT_FILES) $(AWFMINDEX_BUILD_HEADER_FILE)
	$(CC) $(LDFLAGS) -o $(AWFMINDEX_BUILD_SHARED_LIB_FILE) $(OBJECT_FILES)

.PHONY: install
install:
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	mkdir -p $(DESTDIR)$(PREFIX)/include
	cp $(AWFMINDEX_SRC_HEADER_FILE) $(AWFMINDEX_INSTALL_HEADER_FILE)
	#if either the shared lib or static lib are found, install them
ifneq ("$(wildcard $(AWFMINDEX_BUILD_SHARED_LIB_FILE))","")
	cp $(AWFMINDEX_BUILD_SHARED_LIB_FILE) $(AWFMINDEX_INSTALL_SHARED_LIB_FILE)
endif
ifneq ("$(wildcard $(AWFMINDEX_BUILD_STATIC_LIB_FILE))","")
	cp $(AWFMINDEX_BUILD_STATIC_LIB_FILE) $(AWFMINDEX_INSTALL_STATIC_LIB_FILE)
endif


#make the static libs
.PHONY: static
static: $(LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE) $(FASTA_VECTOR_BUILD_STATIC_LIBRARY_FILE) $(AWFMINDEX_BUILD_LIBRARY_DIR) $(OBJECT_FILES) $(AWFMINDEX_BUILD_HEADER_FILE)
	ar rcs $(AWFMINDEX_BUILD_STATIC_LIB_FILE) $(OBJECT_FILES)

clean:
	rm -rf $(AWFMINDEX_BUILD_DIR)
	cd $(LIBDIVSUFSORT_BUILD_DIR) && make clean

.PHONY: uninstall
uninstall:
	rm -rf  $(AWFMINDEX_INSTALL_SHARED_LIB_FILE)
	rm -rf $(AWFMINDEX_INSTALL_STATIC_LIB_FILE)
	rm -rf $(AWFMINDEX_INSTALL_HEADER_FILE)
	cd $(LIBDIVSUFSORT_BUILD_DIR) && make uninstall


# builds libdivsufsort into a static library
$(LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE): $(LIBDIVSUFSORT_PROJECT_DIR)/.git $(LIBDIVSUFSORT_BUILD_DIR) $(LIBDIVSUFSORT_PROJECT_DIR) $(AWFMINDEX_BUILD_DIR)
	cd $(LIBDIVSUFSORT_BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE="Release" -DUSE_OPENMP="ON" -DBUILD_SHARED_LIBS="OFF" -DBUILD_DIVSUFSORT64:BOOL=ON  -DCMAKE_INSTALL_PREFIX="$(AWFMINDEX_BUILD_DIR)" .. && make

$(FASTA_VECTOR_BUILD_STATIC_LIBRARY_FILE): $(LIBDIVSUFSORT_PROJECT_DIR)/.git
	cd $(FASTA_VECTOR_PROJECT_DIR) $$ make static

#create the object files from each c file in the src directory
$(AWFMINDEX_BUILD_DIR)/%.o: $(AWFMINDEX_SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@ -I $(LIBDIVSUFSORT_BUILD_INCLUDE_DIR) -I $(FASTA_VECTOR_BUILD_INCLUDE_DIR)

#initialize the libdivsufsort project submodule if not already done
$(LIBDIVSUFSORT_PROJECT_DIR)/.git:
	git submodule init
	git submodule update

#copy the AwFmIndex.h header to the build directory
$(AWFMINDEX_BUILD_HEADER_FILE): $(AWFMINDEX_BUILD_INCLUDE_DIR)
	cp $(AWFMINDEX_SRC_HEADER_FILE) $(AWFMINDEX_BUILD_HEADER_FILE)

#make the libdivsufsort build directory if not already done
$(LIBDIVSUFSORT_BUILD_DIR): $(LIBDIVSUFSORT_PROJECT_DIR)/.git
	mkdir -p $(LIBDIVSUFSORT_BUILD_DIR)

#make the AwFmIndex build lib directory
$(AWFMINDEX_BUILD_LIBRARY_DIR):
	mkdir -p $(AWFMINDEX_BUILD_LIBRARY_DIR)

#make the AwFmIndex build lib directory
$(AWFMINDEX_BUILD_DIR):
	mkdir -p $(AWFMINDEX_BUILD_DIR)

#make the AwFmIndex build include directory
$(AWFMINDEX_BUILD_INCLUDE_DIR):
	mkdir -p $(AWFMINDEX_BUILD_INCLUDE_DIR)
