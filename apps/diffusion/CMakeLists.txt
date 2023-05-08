add_executable(diffusion diffusion.cpp)
target_link_libraries(diffusion PRIVATE solver)
install(TARGETS diffusion)
