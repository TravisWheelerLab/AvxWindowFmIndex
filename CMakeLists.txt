cmake_minimum_required(VERSION 3.21)
project(
        AWFMIndex
        HOMEPAGE_URL https://github.com/TravisWheelerLab/AvxWindowFmIndex
        LANGUAGES C  # Specify C language explicitly
)

set(CMAKE_C_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#add support for OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY build)


#PUBLIC_H_FILES is missing divsufsort64.h, because it doesn't exist until after divsufsort is built
set(
        PUBLIC_H_FILES

        src/AwFmIndex.h
        lib/FastaVector/src/FastaVector.h
        lib/FastaVector/src/FastaVectorMetadataVector.h
        lib/FastaVector/src/FastaVectorString.h
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

# Set public headers for awfmindex_static
set_target_properties(
    awfmindex_static
    PROPERTIES
    PUBLIC_HEADER "${PUBLIC_H_FILES}"
)

# Set public headers for awfmindex
set_target_properties(
    awfmindex
    PROPERTIES
    PUBLIC_HEADER "${PUBLIC_H_FILES}"
)

add_custom_target(
    copy_headers
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/include/AWFMIndex
)


#copy the libdivsufsort and fastavector static libs to the build director
add_custom_command(
    TARGET awfmindex_static
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:divsufsort64>
            $<TARGET_FILE_DIR:awfmindex_static>
)
add_custom_command(
    TARGET awfmindex_static
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:divsufsort64>
            $<TARGET_FILE:fastavector_static>
            $<TARGET_FILE_DIR:awfmindex_static>
)

# Set common compiler options
target_compile_options(
    awfmindex_static
    PRIVATE
    -mtune=native
    -march=native
    -Wall
    -Wextra
    -O3
    -fPIC
    -lgomp
    ${OpenMP_C_FLAGS}
)

target_compile_options(
    awfmindex
    PRIVATE
    -mtune=native
    -march=native
    -Wall
    -Wextra
    -O3
    -fPIC
    -lgomp
    ${OpenMP_C_FLAGS}
)

# Check if the target architecture is x86_64 (Intel) or aarch64 (ARM)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    # Check if the compiler supports AVX2
    include(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-mavx2" COMPILER_SUPPORTS_AVX2)

    if(COMPILER_SUPPORTS_AVX2)
        # Set AVX2 flags for x86_64 (Intel)
        target_compile_options(awfmindex_static PRIVATE -mavx2)
        target_compile_options(awfmindex PRIVATE -mavx2)
    endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
    # Set NEON flags for ARM
    target_compile_options(awfmindex_static PRIVATE -march=armv8-a+simd)
    target_compile_options(awfmindex PRIVATE -march=armv8-a+simd)
endif()

install(
        TARGETS awfmindex awfmindex_static
)

# ------------
# Dependencies
# ------------

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
set(CMAKE_BUILD_TYPE "Release")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/)

# Custom target for building submodules
add_custom_target(build_submodule
    COMMENT "Building submodules"
)


# Custom command to copy public headers after building submodules
foreach(HEADER ${PUBLIC_H_FILES})
    add_custom_command(
        TARGET build_submodule  # Add dependency on build_submodule
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}
                ${CMAKE_CURRENT_BINARY_DIR}/build
    )
endforeach()

add_dependencies(awfmindex_static build_submodule)
add_dependencies(awfmindex_static copy_headers)

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



# Custom target for building submodules
add_custom_target(build_divsufsort
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/lib/
    COMMENT "Building submodule"
)

set_target_properties(build_divsufsort PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Custom command to move the divsufsort64.h file after the end of the build
add_custom_command(
        TARGET awfmindex_static POST_BUILD # Add dependency on build_submodule
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort/include/divsufsort64.h
        ${CMAKE_CURRENT_BINARY_DIR}/build/
)

add_custom_command(
        TARGET awfmindex_static POST_BUILD # Add dependency on build_submodule
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort/lib/build/libdivsufsort64.a
        ${CMAKE_CURRENT_BINARY_DIR}/build/
)

# Add a dependency on build_submodule
add_dependencies(awfmindex_static build_divsufsort)