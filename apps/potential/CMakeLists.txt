add_executable(potential potential.cpp)
target_link_libraries(potential PRIVATE solver)
install(TARGETS potential)
