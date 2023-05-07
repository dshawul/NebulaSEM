target_include_directories(solver PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(solver PRIVATE
    turbulence.cpp
    turbulence.h 
)
add_subdirectory(ke)
add_subdirectory(kw)
add_subdirectory(les)
add_subdirectory(mixing_length)
add_subdirectory(realizableke)
add_subdirectory(rngke)
