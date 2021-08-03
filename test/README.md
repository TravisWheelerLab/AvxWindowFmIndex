These tests are mostly stand-alone, and do not require the AwFmIndex library to be installed. However, some of these tests require libdivsufsort64 to be installed in a find-able location. Before running these tests, install libdivsufsort64 as a shared library in a manner consistent with the library's readme. Note that the 64-bit version is required, and so the cmake call requires the DBUILD_DIVSUFSORT64 flag to be set:

'''
cd build
cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_BUILD_TYPE="Release" -DUSE_OPENMP="ON" -DBUILD_SHARED_LIBS="ON" -DBUILD_DIVSUFSORT64:BOOL=ON  .. && make
'''
