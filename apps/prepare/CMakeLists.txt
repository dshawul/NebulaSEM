add_executable(prepare prepareApp.cpp)
target_link_libraries(prepare PRIVATE nebulasem)
install(TARGETS prepare)
