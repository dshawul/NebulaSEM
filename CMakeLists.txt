cmake_minimum_required (VERSION 3.16)
set (CMAKE_CXX_STANDARD 17) 
project (nebulasem VERSION 1.0)

# options
option( USE_OMP "Build with OpenMP enabled" OFF )
option( USE_ACC "Build with OpenACC enabled" OFF )
option( USE_SINGLE "Use single precision" OFF )
option( USE_EXPR_TMPL "Use expression templates" ON )

message(STATUS "USE_OMP ............. ${USE_OMP}")
message(STATUS "USE_ACC ............. ${USE_ACC}")
message(STATUS "USE_SINGLE .......... ${USE_SINGLE}")
message(STATUS "USE_EXPR_TMPL ....... ${USE_EXPR_TMPL}")

# find required packages
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

if(USE_OMP)
    find_package(OpenMP COMPONENTS CXX)
endif()

if(USE_ACC)
    set( OpenACC_ACCEL_TARGET tesla:managed )
    find_package(OpenACC COMPONENTS CXX REQUIRED)
endif()

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
