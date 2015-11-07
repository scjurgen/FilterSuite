cmake_minimum_required(VERSION 3.3)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/")

project(FilterSuite)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../../cmake/Modules/")

option(TEST_WITH_FFT "Run tests with qualitative FFT analysis" ON)

add_subdirectory(gmock-1.7.0)
include_directories(SYSTEM gmock-1.7.0/gtest/include)
include_directories(SYSTEM gmock-1.7.0/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra ")
include_directories (${CMAKE_SOURCE_DIR}/..)

add_executable(GenerateFilterSuite
		generator/GenerateFilterSuite.cpp
		ComplexAdditions.cpp
		ComplexAdditions.h
		FilterCoefficients.cpp
		FilterCoefficients.h)

add_executable(FilterCoefficients_unittest
		${GTestLib}
		tests/FilterCoefficients_unittest.cpp
		ComplexAdditions.cpp
		ComplexAdditions.h
		FilterCoefficients.cpp
		FilterCoefficients.h
		)

target_link_libraries(FilterCoefficients_unittest gmock gmock_main)

add_executable(AllFilters_unittest
        ${GTestLib}
        tests/AllFilters_unittest.cpp
		tests/FftBinCompare.cpp
		tests/FftBinCompare.h
		tests/CalcBins.h
        FilterSuite.cpp
        FilterSuite.h
		ComplexAdditions.cpp
		ComplexAdditions.h
        FilterCoefficients.cpp
        FilterCoefficients.h
        )
target_link_libraries(AllFilters_unittest gmock gmock_main)
target_link_libraries (AllFilters_unittest fftw3)



add_executable(FilterSweepSample_Test
		tests/FilterSweepSample_Test.cpp
		tests/FftBinCompare.cpp
		tests/FftBinCompare.h
		tests/CalcBins.h
		ComplexAdditions.cpp
		ComplexAdditions.h
		FilterCoefficients.cpp
		FilterCoefficients.h)
target_link_libraries (FilterSweepSample_Test fftw3)
