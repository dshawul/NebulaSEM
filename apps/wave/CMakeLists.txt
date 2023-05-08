add_executable(wave wave.cpp)
target_link_libraries(wave PRIVATE solver)
install(TARGETS wave)
