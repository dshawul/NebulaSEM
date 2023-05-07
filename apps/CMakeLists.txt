# mesh generation and manipulation
add_executable(mesh meshApp.cpp)
target_link_libraries(mesh PRIVATE solver)

# grid decomposition
add_executable(prepare prepareApp.cpp)
target_link_libraries(prepare PRIVATE solver)

# different types of solvers
add_executable(solvers solversApp.cpp)
target_link_libraries(solvers PRIVATE solver)

#install
install(TARGETS mesh)
install(TARGETS prepare)
install(TARGETS solvers)

