add_executable(euler euler.cpp)
target_link_libraries(euler PRIVATE solver)
install(TARGETS euler)
