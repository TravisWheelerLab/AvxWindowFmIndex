# FastaVector

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/FastaVector)

# libdivsufsort

set(BUILD_EXAMPLES OFF CACHE BOOL "Build examples")
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort64")
set(USE_OPENMP ON CACHE BOOL "Use OpenMP for parallelization")
set(WITH_LFS ON CACHE BOOL "Enable Large File Support")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/libdivsufsort)
