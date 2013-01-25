#ifndef __MSH_MESH_H
#define __MSH_MESH_H

#include "mesh.h"

void readMshMesh(std::istream& is,Mesh::MeshObject& mo);
void writeMshMesh(std::ostream& os,Mesh::MeshObject& mo);

#endif
