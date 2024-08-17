add_library(nebulasem STATIC)

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
    target_link_libraries(nebulasem PUBLIC OpenMP::OpenMP_CXX)
endif()

if(OpenACC_CXX_FOUND)
    target_link_libraries(nebulasem PUBLIC OpenACC::OpenACC_CXX)
    target_compile_options(nebulasem PUBLIC ${OpenACC_CXX_OPTIONS} "-Minfo=accel")
    target_link_options(nebulasem PUBLIC ${OpenACC_CXX_OPTIONS})
endif()

if(MPI_FOUND)
    target_link_libraries(nebulasem PUBLIC MPI::MPI_CXX)
endif()

if(METIS_FOUND)
    target_include_directories(nebulasem PUBLIC ${METIS_INCLUDE_DIRS})
    target_link_libraries(nebulasem PUBLIC ${METIS_LIBRARIES})
endif()

if(NOT USE_SINGLE)
    target_compile_definitions(nebulasem PUBLIC -DUSE_DOUBLE)
endif()

if(USE_EXPR_TMPL)
    target_compile_definitions(nebulasem PUBLIC -DUSE_EXPR_TMPL)
endif()

if(USE_HEX_REFINEMENT)
    target_compile_definitions(nebulasem PUBLIC -DUSE_HEX_REFINEMENT)
endif()

if(USE_STRONG_FORM)
    target_compile_definitions(nebulasem PUBLIC -DUSE_STRONG_FORM)
endif()

install(TARGETS nebulasem)

