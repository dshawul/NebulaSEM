#ifndef __HEX_MESH_H
#define __HEX_MESH_H

#include "mesh.h"

enum {
    LINEAR, GEOMETRIC, WALL, MIXED
};
enum {
    NONE = 0,ARC,COSINE,QUAD
};

struct Edge {
    int type;
    Scalar theta;
    Scalar L;
    Vector N;
    Vertex v[8];
    Edge() {
        type = NONE;
    }
};

struct MergeObject {
    Vertices vb;
    Facets   fb;
};

void hexMesh(Int* n,Scalar* s,Int* type,Vector* vp,Edge* edges,Mesh::MeshObject& mo);
void merge(Mesh::MeshObject&,MergeObject&,Mesh::MeshObject&);
void remove_duplicate(Mesh::MeshObject&);
void merge(Mesh::MeshObject&,MergeObject&);

#endif
