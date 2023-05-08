add_executable(walldist walldist.cpp)
target_link_libraries(walldist PRIVATE solver)
install(TARGETS walldist)
