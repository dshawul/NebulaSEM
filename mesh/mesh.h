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
	void enroll(Util::ParamList& params);
	int  findNearest(const Vector& v);
	IntVector owner(const IntVector&);
	IntVector neighbor(const IntVector&);
}
/*
 * Model for flow close to the wall (Law of the wall).
 *   1 -> Viscous layer
 *   2 -> Buffer layer
 *   3 -> Log-law layer
 * The wall function model is modified for rough surfaces 
 * using Cebecci and Bradshaw formulae.
 */
struct LawOfWall {
	Scalar E;
	Scalar kappa;
	Scalar ks;
	Scalar cks;

	Scalar yLog;

	LawOfWall() : 
		E(9.793),
		kappa(0.4187),
		ks(0.48),
		cks(0.5)
	{
		init();
	}
	void init() {
		yLog = 11.3f;
		for(Int i = 0;i < 20;i++)
			yLog = log(E * yLog) / kappa;
	}
	Scalar getDB(Scalar ustar,Scalar nu) {
		Scalar dB;
		Scalar ksPlus = (ustar * ks) / nu;
		if(ksPlus < 2.25) {
			dB = 0;
		} else if(ksPlus < 90) {
			dB = (1 / kappa) * log((ksPlus - 2.25) / 87.75 + cks * ksPlus)
				             * sin(0.4258 * (log(ksPlus) - 0.811));
		} else {
			dB = (1 / kappa) * log(1 + cks * ksPlus);
		}
		return dB;
	}
	void write(std::ostream& os) const {
		os << "\tE " << E << std::endl;
		os << "\tkappa " << kappa << std::endl;
		os << "\tks " << ks << std::endl;
		os << "\tcks " << cks << std::endl;
	}
	bool read(std::istream& is,std::string str) {
		using namespace Util;
		if(!compare(str,"E")) {
			is >> E;
		} else if(!compare(str,"kappa")) {
			is >> kappa;
		} else if(!compare(str,"ks")) {
			is >> ks;
		} else if(!compare(str,"cks")) {
			is >> cks;
		} else
			return false;
		return true;
	}
};
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
	const Int NOBC         = Util::hash_function("NONE");
}
struct BasicBCondition {
	IntVector* bdry;
	Int     fIndex;
	Int     cIndex;
	std::string cname;
	std::string bname;
	std::string fname;
	bool isWall;
	LawOfWall low;
};
template <class type>
struct BCondition : public BasicBCondition {
	type   value;
	Scalar shape;
	type   tvalue;
	Scalar tshape;
	Scalar zG;
	Vector dir;
	std::vector<type> fixed;

	BCondition(std::string tfname){
		fname = tfname;
		reset();
	}
	void reset() {
		value = tvalue = type(0);
		shape = tshape = zG = Scalar(0);
		dir = Vector(0,0,1);
	}
	void init_indices() {
		isWall = (bname.find("WALL") != std::string::npos);
		bdry = &Mesh::gBoundaries[bname];
		fIndex = Util::hash_function(fname);
		cIndex = Util::hash_function(cname);
	}
};
/*IO*/
template <class type> 
std::ostream& operator << (std::ostream& os, const BCondition<type>& p) {
	os << p.bname << "\n{\n";
	os << "\ttype " << p.cname << std::endl;
	os << "\tvalue " << p.value << std::endl;
	os << "\tshape " << p.shape << std::endl;
	os << "\ttvalue " << p.tvalue << std::endl;
	os << "\ttshape " << p.tshape << std::endl;	
	os << "\tdir " << p.dir << std::endl;
	os << "\tzG " << p.zG << std::endl;
	if(p.isWall)
		p.low.write(os);
	os << "}\n";
	return os;
}
template <class type> 
std::istream& operator >> (std::istream& is, BCondition<type>& p) {
	using namespace Util;
	std::string str;
	char c;

	p.reset();
	is >> p.bname >> c;

	while(c = Util::nextc(is)) {
		if(c == '}') {
			is >> c;
			break;
		} 
		is >> str;
		if(!compare(str,"type")) {
			is >> p.cname;
		} else if(!compare(str,"value")) {
			is >> p.value;
		} else if(!compare(str,"shape")) {
			is >> p.shape;
		} else if(!compare(str,"tvalue")) {
			is >> p.tvalue;
		} else if(!compare(str,"tshape")) {
			is >> p.tshape;
		} else if(!compare(str,"dir")) {
			is >> p.dir;
		} else if(!compare(str,"zG")) {
			is >> p.zG;
		} else if(p.low.read(is,str)) {
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
