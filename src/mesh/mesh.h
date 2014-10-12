#ifndef __MESH_H
#define __MESH_H

#include <cstring>
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
		Int from;
		Int to;
		Int buffer_index;
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
		/*funcs*/
		void write(std::ostream& os);
		void clear() {
			v.clear();
			f.clear();
			c.clear();
			bdry.clear();
			fo.clear();
			fn.clear();
			interMesh.clear();
		}
	};

	extern std::vector<Vector> _fC;
	extern std::vector<Vector> _cC;
	extern std::vector<Vector> _fN;
	extern std::vector<Scalar> _cV;
	extern std::vector<bool>   _reversed;

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
	extern  IntVector        probeCells;
	
	void clear();
	void addBoundaryCells();
	void calcGeometry();
	void removeBoundary(IntVector&);
	bool readMesh(Int = 0,bool = true);
	void enroll(Util::ParamList& params);
	Int  findNearestCell(const Vector& v);
	Int  findNearestFace(const Vector& v);
	void getProbeCells(IntVector&);
	void getProbeFaces(IntVector&);
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
		E(9.8),
		kappa(0.41),
		ks(0),
		cks(0.5)
	{
		init();
	}
	void init() {
		yLog = 11.3f;
		for(Int i = 0;i < 20;i++)
			yLog = log(E * yLog) / kappa;
	}
	Scalar getUstar(Scalar nu,Scalar U,Scalar y) {
		Scalar a = kappa * U * y / nu;
		Scalar yp = a;
		for(Int i = 0;i < 10;i++)
			yp = (a + yp) / (1 + log(E * yp));
		Scalar ustar = yp * nu / y;
		return ustar;
	}
	Scalar getUp(Scalar ustar,Scalar nu,Scalar yp) {
		Scalar up,dB;
		Scalar ksPlus = (ustar * ks) / nu;
		if(ksPlus < 2.25) {
			dB = 0;
		} else if(ksPlus < 90) {
			dB = (1 / kappa) * log((ksPlus - 2.25) / 87.75 + cks * ksPlus)
				             * sin(0.4258 * (log(ksPlus) - 0.811));
		} else {
			dB = (1 / kappa) * log(1 + cks * ksPlus);
		}
		if(yp > yLog)  up = log(E * yp) / kappa - dB;  
		else           up = yp;  
		return up;
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
	const Int ROUGHWALL    = Util::hash_function("ROUGHWALL");
}
struct BasicBCondition {
	IntVector* bdry;
	Int     fIndex;
	Int     cIndex;
	std::string cname;
	std::string bname;
	std::string fname;
	LawOfWall low;
};
template <class type>
struct BCondition : public BasicBCondition {
	type   value;
	Scalar shape;
	type   tvalue;
	Scalar tshape;
	Scalar zMin;
	Scalar zMax;
	Vector dir;
	bool   first;
	bool   read;
	std::vector<type> fixed;

	BCondition(std::string tfname) {
		fname = tfname;
		reset();
	}
	void reset() {
		value = tvalue = type(0);
		shape = tshape = zMin = zMax = Scalar(0);
		dir = Vector(0,0,1);
	}
	void init_indices() {
		bdry = &Mesh::gBoundaries[bname];
		fixed.resize(bdry->size());
		first = true;
		read = false;
		fIndex = Util::hash_function(fname);
		cIndex = Util::hash_function(cname);
	}
};
/*IO*/
template <class type> 
std::ostream& operator << (std::ostream& os, const BCondition<type>& p) {
	os << p.bname << "\n{\n";
	os << "\ttype " << p.cname << std::endl;
	if(!equal(mag(p.value),Scalar(0)))
		os << "\tvalue " << p.value << std::endl;
	if(!equal(p.shape,Scalar(0)))
		os << "\tshape " << p.shape << std::endl;
	if(!equal(mag(p.tvalue),Scalar(0)))
		os << "\ttvalue " << p.tvalue << std::endl;
	if(!equal(p.tshape,Scalar(0)))
		os << "\ttshape " << p.tshape << std::endl;	
	if(!equal(p.dir,Vector(0,0,1)))
		os << "\tdir " << p.dir << std::endl;
	if(p.zMax > 0) {
		os << "\tzMin " << p.zMin << std::endl;
		os << "\tzMax " << p.zMax << std::endl;
	}
	if(p.read) {
		os << "\tfixed " << p.fixed << std::endl;
	}
	if(p.cIndex == Mesh::ROUGHWALL)
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

	while((c = Util::nextc(is))) {
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
		} else if(!compare(str,"zMin")) {
			is >> p.zMin;
		} else if(!compare(str,"zMax")) {
			is >> p.zMax;
		} else if(!compare(str,"fixed")) {
			is >> p.fixed;
			p.read = true;
		} else if(p.low.read(is,str)) {
		}
	}

	p.init_indices();
	p.low.init();
	return is;
}
/*list of all BCS*/
namespace Mesh {
	extern  std::vector<BasicBCondition*> AllBConditions;
	inline void clearBC() {
		forEach(AllBConditions,i)
			delete AllBConditions[i];
		AllBConditions.clear();
	}
}
#endif
