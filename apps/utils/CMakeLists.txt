target_include_directories(nebulasem PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(nebulasem PRIVATE
    properties.cpp
    properties.h
    wrapper.cpp
    wrapper.h
    calc_walldist.cpp
    calc_walldist.h
)
