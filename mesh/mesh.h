#ifndef __MESH_H
#define __MESH_H

#include <string>
#include <vector>
#include <map>
#include "tensor.h"
#include "util.h"

/*Index by ID instead of pointers */
typedef std::vector<Int>      IntVector;

/*our basic building blocks */
enum ENTITY {
	CELL, FACET, VERTEX
};

/*typdefs*/
typedef Vector           Vertex;
typedef IntVector        Facet;  
typedef IntVector        Cell; 

typedef std::vector<Vertex>   Vertices;
typedef std::vector<Facet>    Facets;
typedef std::vector<Cell>     Cells; 
typedef std::map<std::string,IntVector> Boundaries;

/*global mesh*/
namespace Mesh {
	struct interBoundary {
		IntVector* f;
		int from;
		int to;
		int buffer_index;
	};
	struct MeshObject {
		/*vertices , facets and cells */
		Vertices v;
		Facets   f;
		Cells    c;
		/*other info*/
		std::string name;
		Boundaries  bdry;
		IntVector   fo;
		IntVector   fn;
		std::vector<interBoundary> interMesh;
		/*start of boundary cells,facets & vertices*/
		Int      nv;
		Int      nf;
		Int      nc;
	};
	extern  MeshObject       gMesh;
	extern  std::string&     gMeshName;
	extern  Vertices&		 gVertices;
	extern  Facets&			 gFacets;
	extern  Cells&			 gCells;
	extern  Boundaries&      gBoundaries;
	extern  IntVector&       gFO;
	extern  IntVector&       gFN;
	extern  Int&             gBCellsStart;
	extern  std::vector<interBoundary>& gInterMesh;
	extern  Vertices         probePoints;
	
	bool faceInBoundary(Int);
	void addBoundaryCells();
	void readMesh();
	void enroll();
	int  findNearest(const Vector& v);
	IntVector owner(const IntVector&);
	IntVector neighbor(const IntVector&);
}

/*Boundary condition types*/
namespace Mesh {
	const Int DIRICHLET    = Util::hash_function("DIRICHLET");
	const Int NEUMANN      = Util::hash_function("NEUMANN");
	const Int ROBIN        = Util::hash_function("ROBIN");
	const Int SYMMETRY     = Util::hash_function("SYMMETRY");
    const Int CYCLIC       = Util::hash_function("CYCLIC");
	const Int GHOST        = Util::hash_function("GHOST");
	const Int POWER        = Util::hash_function("POWER");
	const Int LOG          = Util::hash_function("LOG");
    const Int PARABOLIC    = Util::hash_function("PARABOLIC");
	const Int INVERSE      = Util::hash_function("INVERSE");
}
struct BasicBCondition {
	IntVector* bdry;
	Int     fIndex;
	Int     cIndex;
	std::string cname;
	std::string bname;
	std::string fname;
};
template <class type>
struct BCondition : public BasicBCondition {
	type   value;
	Scalar shape;
	type   tI;
	Scalar tIshape;
	Scalar zG;
	std::vector<type> fixed;
	BCondition(std::string tfname) {
		fname = tfname;
		zG = 0;
	}
	void init_indices() {
		bdry = &Mesh::gBoundaries[bname];
		fIndex = Util::hash_function(fname);
		cIndex = Util::hash_function(cname);
	}
};
/*IO*/
template <class type> 
std::ostream& operator << (std::ostream& os, const BCondition<type>& p) {
	os << p.bname << " " << p.cname << " ";
	if(p.cname != "SYMMETRY" && p.cname != "CYCLIC" && p.cname != "GHOST") {
		bool has_shape = (p.cname == "PARABOLIC" || p.cname == "LOG" || p.cname == "POWER" || p.cname == "INVERSE");
		os << p.value << " ";
		if(has_shape)
			os << p.shape << " ";
		if(p.cname == "POWER" || p.cname == "LOG")
			os << p.zG << " ";
		if(p.cname != "NEUMANN") {
			os << p.tI << " ";
			if(has_shape)
				os << p.tIshape << " ";
		}
	}
	return os;
}
template <class type> 
std::istream& operator >> (std::istream& is, BCondition<type>& p) {
	char c;
	p.value = p.tI = type(0);
	p.shape = p.tIshape = Scalar(0);
	is >> p.bname >> p.cname;
	if(p.cname != "SYMMETRY" && p.cname != "CYCLIC" && p.cname != "GHOST") {
		bool has_shape = (p.cname == "PARABOLIC" || p.cname == "LOG" || p.cname == "POWER" || p.cname == "INVERSE");
		is >> p.value;
		if(has_shape)
			is >> p.shape;
		if(p.cname == "POWER" || p.cname == "LOG")
			is >> p.zG;
		if((c = Util::nextc(is)) && isdigit(c)) {
			is >> p.tI;
			if((c = Util::nextc(is)) && isdigit(c))
				is >> p.tIshape;
		}
	}
	p.init_indices();
	return is;
}

/*typedefs*/
typedef BCondition<Scalar>   ScalarBCondition;
typedef BCondition<Vector>   VectorBCondition;

/*list of all BCS*/
namespace Mesh {
	extern  std::vector<BasicBCondition*> AllBConditions;
}
#endif
