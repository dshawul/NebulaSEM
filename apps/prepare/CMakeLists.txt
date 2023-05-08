add_executable(prepare prepareApp.cpp)
target_link_libraries(prepare PRIVATE solver)
install(TARGETS prepare)
