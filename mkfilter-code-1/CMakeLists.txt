cmake_minimum_required(VERSION 3.3)
project(FilterSuite)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES mkfilter.cpp complex.cpp)
add_executable(mkfilter ${SOURCE_FILES})


