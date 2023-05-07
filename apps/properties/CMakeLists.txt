target_include_directories(solvers PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(solvers PRIVATE
    properties.cpp
    properties.h
)
