target_include_directories(solver PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(solver PRIVATE
    solve.cpp
    solve.h 
    iteration.h
)
