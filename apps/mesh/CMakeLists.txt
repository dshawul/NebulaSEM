add_executable(mesh meshApp.cpp)
target_link_libraries(mesh PRIVATE nebulasem)
install(TARGETS mesh)
