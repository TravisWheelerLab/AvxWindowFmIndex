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



#configurations for different architectures
AWFMINDEX_SIMD_CONFIG_FLAG_ARM64 = AW_FM_SIMD_CONFIG_ARM_NEON
ifeq ($(ARCH), ARM)
AWFMINDEX_SHARED_LIBRARY_FILENAME = libawfmindex.dylib
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


#file locations
AWFMINDEX_SRC_HEADER_FILE 								= $(AWFMINDEX_SRC_DIR)/$(AWFMINDEX_HEADER_FILENAME)
AWFMINDEX_BUILD_HEADER_FILE								= $(AWFMINDEX_BUILD_INCLUDE_DIR)/$(AWFMINDEX_HEADER_FILENAME)
AWFMINDEX_BUILD_SHARED_LIB_FILE						= $(AWFMINDEX_BUILD_LIBRARY_DIR)/$(AWFMINDEX_SHARED_LIB_FILENAME)
AWFMINDEX_BUILD_STATIC_LIB_FILE						= $(AWFMINDEX_BUILD_LIBRARY_DIR)/$(AWFMINDEX_STATIC_LIB_FILENAME)
AWFMINDEX_INSTALL_HEADER_FILE							= $(AWFMINDEX_INSTALL_INCLUDE_DIR)/$(AWFMINDEX_HEADER_FILENAME)
AWFMINDEX_INSTALL_STATIC_LIB_FILE					= $(AWFMINDEX_INSTALL_LIBRARY_DIR)/$(AWFMINDEX_STATIC_LIB_FILENAME)
AWFMINDEX_INSTALL_SHARED_LIB_FILE					= $(AWFMINDEX_INSTALL_LIBRARY_DIR)/$(AWFMINDEX_SHARED_LIBRARY_FILENAME)
LIBDIVSUFSORT_BUILD_HEADER_FILE						= $(LIBDIVSUFSORT_BUILD_INCLUDE_DIR)/$(LIBDIVSUFSORT_HEADER_FILENAME)
# LIBDIVSUFSORT_BUILD_SHARED_LIBRARY_FILE		= $(LIBDIVSUFSORT_BUILD_LIBRARY_DIR)/$(LIBDIVSUFSORT_LIBRARY_FILENAME)
# LIBDIVSUFSORT_INSTALL_HEADER_FILE 				= $(PREFIX)/include/$(LIBDIVSUFSORT_HEADER_FILENAME)
# LIBDIVSUFSORT_INSTALL_SO_FILE		  				= $(PREFIX)/lib/$(LIBDIVSUFSORT_LIBRARY_FILENAME)

#static libs
LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE	= $(LIBDIVSUFSORT_BUILD_LIBRARY_DIR)/$(LIBDIVSUFSORT_STATIC_LIB_FILENAME)
AWFMINDEX_STATIC_LIB_BUILD_SRC					= $(AWFMINDEX_BUILD_DIR)/$(AWFMINDEX_STATIC_LIB_FILENAME)
# AWFMINDEX_STATIC_LIB_DEST_SRC						= $(AWFMINDEX_INSTALL_LIBRARY_DIR)/$(AWFMINDEX_STATIC_LIB_FILENAME)
# LIBDIVSUFSORT_STATIC_LIB_DEST_SRC 			= $(AWFMINDEX_BUILD_DIR)/$(LIBDIVSUFSORT_STATIC_LIB_FILENAME)
# LIBDIVSUFSORT_STATIC_HEADER_DEST_SRC		= $(AWFMINDEX_BUILD_DIR)/$(LIBDIVSUFSORT_HEADER_FILENAME)


CC 				= gcc
CFLAGS 		= -std=gnu11 -fpic -O3 -mtune=native -mavx2 -Wall -Werror -Wextra -fopenmp -ldivsufsort64
CFLAGS_ARM64	=	-std=gnu11 -fpic -O3 -mtune=native -Wall -Werror -Wextra -D $(AWFMINDEX_SIMD_CONFIG_FLAG_ARM64)

ifeq ($(ARCH), ARM)
CFLAGS = $(CFLAGS_ARM64)
endif

LDFLAGS 	= -shared -L$(LIBDIVSUFSORT_BUILD_LIBRARY_DIR) -I$(LIBDIVSUFSORT_BUILD_INCLUDE_DIR) -ldivsufsort64 # linking flags




SOURCE_FILES 	:= $(wildcard $(AWFMINDEX_SRC_DIR)/*.c)
OBJECT_FILES 	:= $(patsubst $(AWFMINDEX_SRC_DIR)/%, $(AWFMINDEX_BUILD_DIR)/%, $(SOURCE_FILES:.c=.o))


#rules
.PHONY: all
all: $(AWFMINDEX_BUILD_DIR) $(LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE) $(OBJECT_FILES)
	$(CC) $(LDFLAGS) -o $(AWFMINDEX_BUILD_SHARED_LIB_FILE) $(OBJECT_FILES)


.PHONY: install
install: $(AWFMINDEX_BUILD_SHARED_LIB_FILE) $(AWFMINDEX_SRC_HEADER_FILE)
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	mkdir -p $(DESTDIR)$(PREFIX)/include
	cp $(AWFMINDEX_SRC_HEADER_FILE) $(AWFMINDEX_INSTALL_HEADER_FILE)
	#if either the shared lib or static lib are found, install them
ifneq ("$(wildcard $(AWFMINDEX_BUILD_SHARED_LIB_FILE))","")
	cp $(AWFMINDEX_BUILD_SHARED_LIB_FILE) $(AWFMINDEX_INSTALL_SHARED_LIB_FILE)
endif
ifneq ("$(wildcard $(AWFMINDEX_STATIC_LIB_BUILD_SRC))","")
	cp $(AWFMINDEX_STATIC_LIB_BUILD_SRC) $(AWFMINDEX_INSTALL_STATIC_LIB_FILE)
endif

	# cp $(AWFMINDEX_BUILD_SHARED_LIB_FILE) $(AWFMINDEX_INSTALL_SHARED_LIB_FILE)
	# cp $(AWFMINDEX_STATIC_LIB_BUILD_SRC) $(AWFMINDEX_INSTALL_STATIC_LIB_FILE)
	#cd into libdivsufsort to make install
	# cd $(LIBDIVSUFSORT_BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE="Release" -DBUILD_DIVSUFSORT64:BOOL=ON  -DCMAKE_INSTALL_PREFIX="$(PREFIX)" .. && make install
	# cp $(LIBDIVSUFSORT_BUILD_HEADER_FILE)  $(LIBDIVSUFSORT_INSTALL_HEADER_FILE)
	# cp $(LIBDIVSUFSORT_BUILD_LIBRARY_FILE) $(LIBDIVSUFSORT_INSTALL_SO_FILE)

#make the static libs
.PHONY: static
static: $(LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE) $(AWFMINDEX_BUILD_LIBRARY_DIR) $(OBJECT_FILES)
	ar rcs $(AWFMINDEX_STATIC_LIB_BUILD_SRC) $(OBJECT_FILES)
	# cp $(LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE) $(LIBDIVSUFSORT_STATIC_LIB_DEST_SRC)
	# cp $(LIBDIVSUFSORT_BUILD_HEADER_FILE) $(LIBDIVSUFSORT_STATIC_HEADER_DEST_SRC)

clean:
	rm -rf $(AWFMINDEX_BUILD_DIR)
	cd $(LIBDIVSUFSORT_BUILD_DIR) && make clean

.PHONY: uninstall
uninstall:
	rm -rf  $(AWFMINDEX_INSTALL_SHARED_LIB_FILE)
	rm -rf $(AWFMINDEX_INSTALL_STATIC_LIB_FILE)
	rm -rf $(AWFMINDEX_INSTALL_HEADER_FILE)
	cd $(LIBDIVSUFSORT_BUILD_DIR) && make uninstall
	# rm -rf $(AWFMINDEX_STATIC_LIB_DEST_SRC)
	# rm -rf $(LIBDIVSUFSORT_STATIC_LIB_DEST_SRC)


# builds libdivsufsort into a static library
$(LIBDIVSUFSORT_BUILD_STATIC_LIBRARY_FILE): $(LIBDIVSUFSORT_PROJECT_DIR)/.git $(LIBDIVSUFSORT_PROJECT_DIR) $(AWFMINDEX_BUILD_DIR)
	cd $(LIBDIVSUFSORT_BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE="Release" -DUSE_OPENMP="ON" -DBUILD_SHARED_LIBS="OFF" -DBUILD_DIVSUFSORT64:BOOL=ON  -DCMAKE_INSTALL_PREFIX="$(AWFMINDEX_BUILD_DIR)" .. && make


#create the object files from each c file in the src directory
$(AWFMINDEX_BUILD_DIR)/%.o: $(AWFMINDEX_SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@ -I $(LIBDIVSUFSORT_BUILD_INCLUDE_DIR)
# ifeq ($(ARCH), M1)
# 	$(CC) $(CFLAGS_M1) -D $(SIMD_M1_MACRO_DEFINE_FLAG) -c $< -o $@ -I $(LIBDIVSUFSORT_BUILD_INCLUDE_DIR)
# else
# endif

# #copies the AwFmIndex.h header to the build directory
# $(AWFMINDEX_BUILD_HEADER_FILE): $(AWFMINDEX_SRC_HEADER_FILE) $(AWFMINDEX_BUILD_INCLUDE_DIR)
# 	cp $(AWFMINDEX_SRC_HEADER_FILE) $(AWFMINDEX_BUILD_HEADER_FILE)

# #building libdivsufsort
# $(LIBDIVSUFSORT_BUILD_HEADER_FILE): $(LIBDIVSUFSORT_BUILD_DIR)
# 	cd $(LIBDIVSUFSORT_BUILD_DIR) &&	cmake -DCMAKE_BUILD_TYPE="Release" -DUSE_OPENMP="ON" -DBUILD_SHARED_LIBS="ON" -DBUILD_DIVSUFSORT64:BOOL=ON  -DCMAKE_INSTALL_PREFIX="$(PREFIX)" .. && make

#initialize the libdivsufsort project submodule if not already done
$(LIBDIVSUFSORT_PROJECT_DIR)/.git:
	git submodule init
	git submodule update

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
