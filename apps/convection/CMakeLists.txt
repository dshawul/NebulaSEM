add_executable(convection convection.cpp)
target_link_libraries(convection PRIVATE nebulasem)
install(TARGETS convection)
