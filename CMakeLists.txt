cmake_minimum_required(VERSION 2.8)

# https://cmake.org/cmake/help/v3.0/policy/CMP0048.html#policy:CMP0048
# The purpose of this command is to remove the error message about empty project version string related to this new cmake policy.
if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

set(PROJECT_NAME mc-mini)
set(PROJECT_VERSION 1.1.0)
project(${PROJECT_NAME} VERSION ${PROJECT_VERSION})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Enable/disable testing
option(TESTS_ENABLED "Enable automatic tests" ON)
# Enable/disable test coverage
option(COVERAGE_ENABLED "Enable test coverage" OFF)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories(${CMAKE_SOURCE_DIR}/include)

if(COVERAGE_ENABLED)
  set(CMAKE_CXX_FLAGS_DEBUG "-Wall -fPIC -O0 -ggdb -fprofile-arcs -ftest-coverage")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-Wall -O0 -ggdb")
endif()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(FATAL_ERROR "Compiler does not support C++11")
endif()

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_SOURCE_DIR}/docs/Doxyfile @ONLY)
  add_custom_target(docs
    ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/docs/Doxyfile
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, defaulting to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# Boost, a collection of useful C++ utilities and idioms
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
set(LIBRARIES ${LIBRARIES} ${Boost_LIBRARIES})

# Eigen, a powerful linear algebra header-only library
# Set EIGEN3_INCLUDE_DIR if it's set in ENV. This way TravisCI can use a
# local copy of Eigen instead of the outdated version on apt.
if(DEFINED ENV{EIGEN3_INCLUDE_DIR})
  set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
  message(STATUS "Using local Eigen3 headers")
endif()
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# HDF5, an output library
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})
set(LIBRARIES ${LIBRARIES} ${HDF5_LIBRARIES})

if(TESTS_ENABLED)
  add_subdirectory(test)
endif()

add_subdirectory(source)

# vim:ft=cmake
