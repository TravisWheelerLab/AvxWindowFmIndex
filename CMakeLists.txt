cmake_minimum_required(VERSION 3.21)
project(
        AWFMIndex
        VERSION 1.0.0
        HOMEPAGE_URL https://github.com/TravisWheelerLab/AvxWindowFmIndex
)

set(CMAKE_C_STANDARD 17)
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

add_subdirectory(lib/FastaVector)

find_package(OpenMP REQUIRED)

add_library(
        awfmindex STATIC

        ${C_FILES}
)

target_link_libraries(awfmindex PRIVATE OpenMP)

target_link_libraries(awfmindex PRIVATE fastavector)

#target_link_directories(awfmindex PRIVATE lib/FastaVector/src/)
#target_include_directories(awfmindex PRIVATE lib/FastaVector/src/)

target_link_directories(awfmindex PRIVATE lib/libdivsufsort/lib/)
target_include_directories(awfmindex PRIVATE lib/libdivsufsort/lib/)

install(TARGETS awfmindex ARCHIVE)

# Runs clang-format on the code
add_custom_target(
        format
        COMMAND clang-format -i ${C_FILES} ${H_FILES}
        VERBATIM
)

# Checks whether clang-format needs to be run and sets a non-zero
# exit code if so
add_custom_target(
        check-format
        COMMAND clang-format --dry-run --Werror ${C_FILES} ${H_FILES}
        VERBATIM
)
