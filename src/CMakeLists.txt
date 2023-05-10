add_library(solver STATIC)

add_subdirectory(field)
add_subdirectory(mesh)
add_subdirectory(mp)
add_subdirectory(prepare)
add_subdirectory(solvers)
add_subdirectory(tensor)
add_subdirectory(util)
add_subdirectory(vtk)
add_subdirectory(turbulence)

if(OpenMP_FOUND)
    target_link_libraries(solver PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MPI_FOUND)
    target_link_libraries(solver PUBLIC MPI::MPI_CXX)
endif()

if(METIS_FOUND)
    target_link_libraries(solver PUBLIC ${METIS_LIBRARIES})
endif()

install(TARGETS solver)

