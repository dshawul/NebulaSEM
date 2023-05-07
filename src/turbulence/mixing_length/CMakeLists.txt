target_include_directories(solver PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(solver
    PRIVATE
    mixing_length.cpp
    mixing_length.h 
)
