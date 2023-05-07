# mesh generation and manipulation
add_executable(mesh meshApp.cpp)
target_link_libraries(mesh PRIVATE solver)

# grid decomposition
add_executable(prepare prepareApp.cpp)
target_link_libraries(prepare PRIVATE solver)

# different types of solvers
add_executable(solvers solversApp.cpp)
target_link_libraries(solvers PRIVATE solver)

# include directory
target_include_directories(solvers PUBLIC ${CMAKE_SOURCE_DIR})

# subdirectories
add_subdirectory(properties)
add_subdirectory(walldist)
add_subdirectory(potential)
add_subdirectory(convection)
add_subdirectory(diffusion)
add_subdirectory(transport)
add_subdirectory(wave)
add_subdirectory(piso)
add_subdirectory(euler)
add_subdirectory(hydro_balance)

#install
install(TARGETS mesh)
install(TARGETS prepare)
install(TARGETS solvers)

