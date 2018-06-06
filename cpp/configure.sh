#!/bin/bash
CMAKE_BUILD_TYPE=Release
CMAKE_C_COMPILER=gcc-8
#CMAKE_C_COMPILER=gcc
#CMAKE_C_COMPILER=/opt/intel/compilers_and_libraries/linux/bin/intel64/icc

CMAKE_CXX_COMPILER=g++-8
#CMAKE_CXX_COMPILER=g++
#CMAKE_CXX_COMPILER=/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc


cmake \
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} \
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} \
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
    ..
