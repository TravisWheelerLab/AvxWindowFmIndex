These tests require that the project be built for static libraries. Nothing needs to be installed for these to function. To build the static libraries, cd into the project's home directory and run
'''
cmake -DBUILD_SHARED_LIBS="OFF" . && make
'''
