target_include_directories(nebulasem PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(nebulasem PRIVATE
    tensor.cpp
    tensor.h 
    types.h
)
