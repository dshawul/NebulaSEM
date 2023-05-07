cmake_minimum_required (VERSION 3.18)
set (CMAKE_CXX_STANDARD 11) 
project (solver)

# find required packages
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(MPI COMPONENTS CXX REQUIRED)
find_package(METIS COMPONENTS CXX REQUIRED)

# set release compiler flags
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops -march=native")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -unroll -xHost")
endif()

# add subdirectories
add_subdirectory(src)
add_subdirectory(apps)
