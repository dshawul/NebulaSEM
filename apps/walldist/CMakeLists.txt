add_executable(walldist walldist.cpp)
target_link_libraries(walldist PRIVATE nebulasem)
install(TARGETS walldist)
