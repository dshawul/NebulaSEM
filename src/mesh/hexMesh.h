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

inline Vector center(const Vector& v1,const Vector& v2,const Vector& v3) {
	Vector v12 = v1 - v2;
	Vector v13 = v1 - v3;
	Vector v23 = v2 - v3;
	Scalar d = 2 * magSq(v12 ^ v23);
	Scalar a = magSq(v23) * (v12 & v13) / d;
	Scalar b = magSq(v13) * (-v12 & v23) / d;
	Scalar c = magSq(v12) * (v13 & v23) / d;
	return a * v1 + b * v2 + c * v3;
}

#endif
