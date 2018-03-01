#!/usr/bin/env bash

# Add user project to Chaste/projects/ by
# $ ln -s /path/to/user/project
# Then run this bash script in the build directory (assuming Chaste is located
# at the same parent directory as the build directory) or follow the
# instruction below.

# CMake: to construct a `Makefile` for `Make`. This constains instructions for
# the compiler.
# Use a specific version of cmake 3.6.1 (3.10.1 does not work!)
# Specify the PETSC_DIR and PETSC_ARCH, as it usually cannot find it (and I
# have too many version in it)...
PETSC_DIR=/home/scratch/chaste-libs/petsc-3.6.2 PETSC_ARCH=linux-gnu  ../cmake-3.6.1-Linux-x86_64/bin/cmake ../Chaste

# Make: executes what is in the `Makefile` for a given target
make TestHello

# CTest: run the test!
# -V: verbose
# -R: one single test
# ^,$: specify the start and end of the test name; else it might match extra 
#      stuffs (using RegExp!)
ctest -V -R ^TestHello$

