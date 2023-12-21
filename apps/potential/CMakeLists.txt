add_executable(potential potential.cpp)
target_link_libraries(potential PRIVATE nebulasem)
install(TARGETS potential)
