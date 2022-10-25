cmake_minimum_required(VERSION 3.7)
project(pathassociation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fopenmp -std=c++0x -O3")


set(SOURCE_FILES
        associationrulemining.cpp
        associationrulemining.h
        graph.cpp
        graph.h
        main.cpp
        utility.h)

add_executable(pathassociation ${SOURCE_FILES})
