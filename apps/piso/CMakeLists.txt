add_executable(piso piso.cpp)
target_link_libraries(piso PRIVATE nebulasem)
install(TARGETS piso)
