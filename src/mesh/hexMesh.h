#ifndef __HEX_MESH_H
#define __HEX_MESH_H

#include "mesh.h"

/** Element edges refinement options */
enum {
    LINEAR,     /**< Uniform spacing */
    GEOMETRIC,  /**< Geometrically varying spacing */
    WALL,       /**< Larger spacing at the center */
    MIXED       /**< Mixed mode for different edges of the element */
};

/** Edge shapes */
enum {
    NONE = 0,   /**< Straight edge */
    ARC,        /**< Circular arc */
    COSINE,     /**< Cosine shaped edge */
    QUAD        /**< Quadratic shaped edge */
};

/** Edge of an element */
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
