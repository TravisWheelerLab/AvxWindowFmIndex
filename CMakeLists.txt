cmake_minimum_required(VERSION 3.21)
project(
        AWFMIndex
        HOMEPAGE_URL https://github.com/TravisWheelerLab/AvxWindowFmIndex
        LANGUAGES C  # Specify C language explicitly
)

set(CMAKE_C_STANDARD 17)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY build)

set(
        PUBLIC_H_FILES

        src/AwFmIndex.h
)

set(
        H_FILES

        src/AwFmCreate.h
        src/AwFmFile.h
        src/AwFmIndex.h
        src/AwFmIndexStruct.h
        src/AwFmKmerTable.h
        src/AwFmLetter.h
        src/AwFmOccurrence.h
        src/AwFmParallelSearch.h
        src/AwFmSearch.h
        src/AwFmSimdConfig.h
        src/AwFmSuffixArray.h
)
set(
        C_FILES

        src/AwFmCreate.c
        src/AwFmFile.c
        src/AwFmIndexStruct.c
        src/AwFmKmerTable.c
        src/AwFmLetter.c
        src/AwFmOccurrence.c
        src/AwFmParallelSearch.c
        src/AwFmSearch.c
        src/AwFmSimdConfig.c
        src/AwFmSuffixArray.c
)

add_library(
        awfmindex_static STATIC

        ${C_FILES}
)

add_library(
        awfmindex SHARED

        ${C_FILES}
)

set_target_properties(
        awfmindex_static awfmindex

        PROPERTIES
        PUBLIC_HEADER ${PUBLIC_H_FILES}
)

target_compile_options(
        awfmindex_static
        PRIVATE
        -mtune=native
        -march=native
        -mavx2
        -Wall
        -Wextra
        ${OpenMP_C_FLAGS}
)

target_compile_options(
        awfmindex
        PRIVATE
        -mtune=native
        -march=native
        -mavx2 
        -Wall
        -Wextra
        ${OpenMP_C_FLAGS}
)

install(
        TARGETS awfmindex awfmindex_static
)

# ------------
# Dependencies
# ------------

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/)

# FastaVector

target_include_directories(
        awfmindex_static PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/FastaVector/src/
)
target_include_directories(
        awfmindex PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/FastaVector/src/
)

target_link_libraries(awfmindex_static PRIVATE fastavector_static)
target_link_libraries(awfmindex PRIVATE fastavector)

# libdivsufsort

target_include_directories(
        awfmindex PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/libdivsufsort/include/
)
target_include_directories(
        awfmindex_static PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/libdivsufsort/include/
)

target_link_libraries(awfmindex_static PRIVATE divsufsort64)
target_link_libraries(awfmindex PRIVATE divsufsort64)
