add_executable(diffusion diffusion.cpp)
target_link_libraries(diffusion PRIVATE nebulasem)
install(TARGETS diffusion)
