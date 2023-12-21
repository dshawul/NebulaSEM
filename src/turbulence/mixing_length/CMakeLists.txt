target_include_directories(nebulasem PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(nebulasem
    PRIVATE
    mixing_length.cpp
    mixing_length.h 
)
