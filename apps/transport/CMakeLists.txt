add_executable(transport transport.cpp)
target_link_libraries(transport PRIVATE solver)
install(TARGETS transport)
