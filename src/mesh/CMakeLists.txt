target_include_directories(nebulasem PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(nebulasem PRIVATE
    hexMesh.cpp
    hexMesh.h 
    mesh.cpp
    mesh.h
    mshMesh.cpp
)
