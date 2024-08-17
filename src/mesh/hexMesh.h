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
    GAUSSIAN,   /**< Gaussian shaped edge */
    RIDGE,      /**< Linear edge */
    SCHAR,      /**< Shcar mountain edge */
};

/** Edge of an element */
struct Edge {
    int type;     /**< Edge shape type */
    union {
        Scalar theta; /**< Subtended angle of arc */
        Scalar scale; /**< Scale of feature */
    };
    Vector center;/**< Center of the arc */
    Vector N;     /**< Unit normal vector of edge */
    Vertex v[4];  /**< Points defining shape of edge */
    Edge() {
        type = NONE;
    }
};

/** Boundary faces and vertices for merging hex grids*/
struct MergeObject {
    Vertices vb;  /**< List of boundary vertices */
    Facets   fb;  /**< List of boundary faces */
};

void hexMesh(Int* n,Scalar* s,Int* type,Vector* vp,Edge* edges,
             Mesh::MeshObject& mo);
void merge(Mesh::MeshObject&,MergeObject&,Mesh::MeshObject&);
void remove_duplicate(Mesh::MeshObject&);
void merge(Mesh::MeshObject&,MergeObject&);

#endif
