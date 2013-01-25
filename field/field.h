#ifndef __FIELD_H
#define __FIELD_H

#include <list>
#include <sstream>
#include "mesh.h"
#include "mp.h"

/*******************************************************************************
 *                              Control parameters
 *******************************************************************************/
namespace Controls {

	enum Scheme{
		CDS,UDS,HYBRID,BLENDED,LUD,CDSS,MUSCL,QUICK,
		VANLEER,VANALBADA,MINMOD,SUPERBEE,SWEBY,QUICKL,UMIST,
		DDS,FROMM
	};
	enum NonOrthoScheme {
		NO_CORRECTION,MINIMUM, ORTHOGONAL, OVER_RELAXED
	};
	enum Solvers {
		JACOBI, SOR, PCG
	};
	enum CommMethod {
		BLOCKED, ASYNCHRONOUS
	};
	enum State {
		STEADY, TRANSIENT
	};

	extern Scheme convection_scheme;
	extern Int TVDbruner;
	extern Scheme interpolation_scheme;
	extern NonOrthoScheme nonortho_scheme;
	extern Solvers Solver; 
	extern CommMethod ghost_exchange;
	extern CommMethod parallel_method;
	extern State state;

	extern Scalar SOR_omega;
	extern Scalar tolerance;
	extern Scalar blend_factor;
	extern Scalar time_scheme_factor;
	extern Scalar dt;
	
	extern Int max_iterations;
	extern Int write_interval;
	extern Int start_step;
	extern Int end_step;
}

/* *****************************************************************************
 *                    Field variables defined on mesh                          
 * *****************************************************************************/
enum ACCESS {
	NONE = 0, READ = 1, WRITE = 2,READWRITE = 3
};
template <class type,ENTITY entity> 
class MeshField {
private:
	type*        P;
	int          allocated;
public:
	ACCESS       access;
	Int          fIndex;
	std::string  fName;

	/*common*/
	static Int SIZE;
	static const Int TYPE_SIZE = sizeof(type) / sizeof(Scalar);
	static std::list<MeshField*> fields_;
	static std::list<type*> mem_;

    /*constructors*/
	MeshField(const char* str = "", ACCESS a = NONE) : 
				P(0),allocated(0),access(a),fName(str) {
		construct(str,a);
	}
	MeshField(const MeshField& p) : allocated(0) {
		allocate(); 
		for(Int i = 0;i < SIZE;i++) 
			P[i] = p[i];
	}
	MeshField(const type& p) : allocated(0) {
		allocate(); 
		for(Int i = 0;i < SIZE;i++) 
			P[i] = p;
	}
	explicit MeshField(const bool) : allocated(0) {
	}
	/*allocators*/
	void allocate() {
		if(mem_.empty()) {
			switch(entity) {
				case CELL:   SIZE = Mesh::gCells.size();	 break;
				case FACET:  SIZE = Mesh::gFacets.size();    break;
				case VERTEX: SIZE = Mesh::gVertices.size();  break;
			}
			P = new type[SIZE];
		} else {
			P = mem_.front();
			mem_.pop_front();
		}
		allocated = 1;
	}
	void allocate(std::vector<type>& q) {
		switch(entity) {
				case CELL:   SIZE = Mesh::gCells.size();	 break;
				case FACET:  SIZE = Mesh::gFacets.size();    break;
				case VERTEX: SIZE = Mesh::gVertices.size();  break;
		}
		P = &q[0];
		allocated = 0;
	}
	void construct(const char* str = "", ACCESS a = NONE) {
		access = a;
		fName = str;
		if(Mesh::gCells.size())
			allocate();
		fIndex = Util::hash_function(str);
		if(fIndex)
			fields_.push_back(this);
	}
	/*d'tor re-cycles memory */
	~MeshField() {
		if(allocated && !Util::Terminated) {
			mem_.push_front(P);
			if(fIndex)
				fields_.remove(this);
		}
	}

	/*static functions*/
	void readInternal(std::istream&);
	void read(Int step);
	void write(Int step);
	static void readAll(Int step);
	static void writeAll(Int step);
	static void write_vtk(std::ostream&,bool);
	static int count_writable();

	/*accessors*/
	Int size() const {
		return SIZE;
	}
	type& operator [] (Int i) const {
		return P[i];
	}
	/*unary ops*/
	MeshField operator - () {
		MeshField r;
		for(Int i = 0;i < SIZE;i++) 
			r[i] = -P[i];
		return r;
	}
	friend MeshField<Scalar,entity> operator & (const MeshField& p,const MeshField& q) {
		MeshField<Scalar,entity> r;
		for(Int i = 0;i < SIZE;i++)
			r[i] = p[i] & q[i];
		return r;
	}
	/*unrolled operations*/
#define Op($)													    \
	MeshField& operator $(const MeshField& q) {						\
		for(Int i = 0;i < SIZE;i++)									\
			P[i] $ q[i];											\
		return *this;												\
	}
#define SOp($)														\
	MeshField& operator $(const Scalar& q) {						\
		for(Int i = 0;i < SIZE;i++)									\
			P[i] $ q;												\
		return *this;												\
	}
#define Fp(name)													\
	friend MeshField name(const MeshField& p,const MeshField& s) {	\
		MeshField r;												\
		for(Int i = 0;i < SIZE;i++)								    \
			r[i] = name(p[i],s[i]);									\
		return r;													\
	}
#define Fp1(name)													\
	friend MeshField name(const MeshField& p,const Scalar& s) {		\
		MeshField r;												\
		for(Int i = 0;i < SIZE;i++)								    \
			r[i] = name(p[i],s);									\
		return r;													\
	}
#define Fp2(name)													\
	friend MeshField name(const MeshField& p) {						\
		MeshField r;												\
		for(Int i = 0;i < SIZE;i++)								    \
			r[i] = name(p[i]);										\
		return r;													\
	}
    /*define ops*/
	Op(=);
	Op(+=);
	Op(-=);
	Op(*=);
	Op(/=);
	SOp(=);
	SOp(+=);
	SOp(-=);
	SOp(*=);
	SOp(/=);
	Fp(sdiv);
	/*from math.h*/
	Fp2(acos);
	Fp2(asin);
	Fp2(atan);
	Fp(atan2);
	Fp2(ceil);
	Fp2(cos);
	Fp2(cosh);
	Fp2(exp);
	Fp2(fabs);
	Fp2(floor);
	Fp2(log);
	Fp2(log10);
	Fp1(pow);
	Fp2(sin);
	Fp2(sinh);
	Fp2(sqrt);
	Fp2(tan);
	Fp2(tanh);
	Fp(min);
	Fp(max);
	/*additional*/
	Fp2(unit);
#undef Op
#undef SOp
#undef Fp
#undef Fp1
#undef Fp2
	Operator(MeshField,+);
	Operator(MeshField,-);
	Operator(MeshField,*);
	Operator(MeshField,/);
	/*friend ops*/
	friend MeshField<Scalar,entity> mag(const MeshField& p) {
		MeshField<Scalar,entity> r;
		for(Int i = 0;i < p.size();i++) 
			r[i] = mag(p[i]);
		return r;
	}
	friend MeshField dev(const MeshField& p,const Scalar factor = 1.) {
		MeshField r;
		for(Int i = 0;i < p.size();i++) 
			r[i] = dev(p[i],factor);
		return r;
	}
	friend MeshField hyd(const MeshField& p,const Scalar factor = 1.) {
		MeshField r;
		for(Int i = 0;i < p.size();i++) 
			r[i] = hyd(p[i],factor);
		return r;
	}
	/*relax*/
	void Relax(const MeshField& po,Scalar UR) {
		for(Int i = 0;i < SIZE;i++) 
			P[i] = po[i] + (P[i] - po[i]) * UR;
	}
	/*end*/
};
/***********************************
 *  Specific tensor operations
 ***********************************/

/* Default operator overload for scalar fields*/
#define Op(name,F,S)																			\
	template<class T,ENTITY E>																	\
	MeshField<T,E> name(const MeshField<F,E>& p,const MeshField<S,E>& q) {						\
		MeshField<T,E> r;																		\
		for(Int i = 0;i < p.size();i++)															\
			r[i] = name(p[i],q[i]);																\
		return r;																				\
	}
Op(operator *,Scalar,T);
Op(operator /,Scalar,T);
Op(operator *,T,Scalar);
Op(operator /,T,Scalar);
#undef Op

/*multiply*/
template <ENTITY E>
MeshField<Tensor,E> mul(const MeshField<Vector,E>& p,const MeshField<Vector,E>& q) {
	MeshField<Tensor,E> r;
	for(Int i = 0;i < p.size();i++) 
		r[i] = mul(p[i],q[i]);
	return r;
}
template <ENTITY E> 
inline MeshField<Vector,E> mul(const MeshField<Vector,E>& p,const MeshField<Scalar,E>& q) { 
	return p * q; 
}
template <class T,ENTITY E>
MeshField<T,E> mul(const MeshField<T,E>& p,const MeshField<T,E>& q) {
	MeshField<T,E> r;
	for(Int i = 0;i < p.size();i++) 
		r[i] = mul(p[i],q[i]);
	return r;
}
/*dot*/
template <ENTITY E,Int SIZE> 
MeshField<Vector,E> dot(const MeshField<TTensor<SIZE>,E>& p,const MeshField<Vector,E>& q) {
	MeshField<Vector,E> r;
	for(Int i = 0;i < p.size();i++)
		r[i] = dot(q[i],p[i]);
	return r;
}
template <ENTITY E> 
inline MeshField<Scalar,E> dot(const MeshField<Vector,E>& p,const MeshField<Vector,E>& q) { 
	return p & q; 
}
/*symmetric & skew-symmetric*/
template <ENTITY E>
MeshField<STensor,E> sym(const MeshField<Tensor,E>& p) {
	MeshField<STensor,E> r;
	for(Int i = 0;i < p.size();i++) 
		r[i] = sym(p[i]);
	return r;
}
template <ENTITY E>
MeshField<Tensor,E> skw(const MeshField<Tensor,E>& p) {
	MeshField<Tensor,E> r;
	for(Int i = 0;i < p.size();i++) 
		r[i] = skw(p[i]);
	return r;
}
/*transpose*/
template <ENTITY E>
MeshField<Tensor,E> trn(const MeshField<Tensor,E>& p) {
	MeshField<Tensor,E> r;
	for(Int i = 0;i < p.size();i++) 
		r[i] = trn(p[i]);
	return r;
}
/*static variables*/
template <class T,ENTITY E> 
std::list<MeshField<T,E>*> MeshField<T,E>::fields_;

template <class T,ENTITY E> 
std::list<T*> MeshField<T,E>::mem_;

template <class T,ENTITY E> 
Int MeshField<T,E>::SIZE;

/* typedefs */
typedef MeshField<Scalar,CELL>    ScalarCellField;
typedef MeshField<Scalar,FACET>   ScalarFacetField;
typedef MeshField<Scalar,VERTEX>  ScalarVertexField;
typedef MeshField<Vector,CELL>    VectorCellField;
typedef MeshField<Vector,FACET>   VectorFacetField;
typedef MeshField<Vector,VERTEX>  VectorVertexField;
typedef MeshField<Tensor,CELL>    TensorCellField;
typedef MeshField<Tensor,FACET>   TensorFacetField;
typedef MeshField<Tensor,VERTEX>  TensorVertexField;
typedef MeshField<STensor,CELL>   STensorCellField;
typedef MeshField<STensor,FACET>  STensorFacetField;
typedef MeshField<STensor,VERTEX> STensorVertexField;

/* ***************************************
 * global mesh fields
 * ***************************************/
namespace Mesh {
	extern VectorVertexField vC;
	extern VectorFacetField  fC;
	extern VectorCellField   cC;
	extern VectorFacetField  fN;
	extern ScalarCellField   cV;
	extern ScalarFacetField  fI;
	extern ScalarCellField   yWall;

	void   initGeomMeshFields(bool = true);
	void   write_fields(Int);
	void   read_fields(Int);
	void   calc_walldist(Int);
}
/* **********************************************
 *  Input - output operations
 * **********************************************/
template <class T,ENTITY E> 
void MeshField<T,E>::readInternal(std::istream& is) {
	using namespace Mesh;
	/*size*/
	char c;
	int size;
	std::string str;
	is >> str >> size;

	/*internal field*/
	if((c = Util::nextc(is)) && isalpha(c)) {
		T value = T(0);
		is >> str;
		if(str == "uniform")
			is >> value;
		*this = value;
	} else {
		char symbol;
		is >> size >> symbol;
		for(int i = 0;i < size;i++) {
			is >> (*this)[i];
		}
		is >> symbol;
	}
}
template <class T,ENTITY E> 
void MeshField<T,E>::read(Int step) {
	using namespace Mesh;

	/*open*/
	std::stringstream path;
	path << fName << step;
	std::ifstream is(path.str().c_str());

	if(is.fail()) {
		*this = T(0);
	} else {
		std::cout << "Reading " << fName << step  << std::endl;
		std::cout.flush();

		/*internal*/
		readInternal(is);

		/*boundary*/
		char c;
		BCondition<T>* bc;
		while((c = Util::nextc(is)) && isalpha(c)) {
			bc = new BCondition<T>(this->fName);
			is >> *bc;
			AllBConditions.push_back(bc);
		}
		
		/*update BCs*/
		updateExplicitBCs(*this,true,true);
	}
	/*close*/
	is.close();
}

template <class T,ENTITY E> 
void MeshField<T,E>::readAll(Int step) {
	using namespace Mesh;
	MeshField<T,E>* pf;
	for(typename std::list<MeshField<T,E>*>::iterator it = fields_.begin();it != fields_.end();++it) {
		pf = *it;
		if(pf->access & READ)
			pf->read(step);
	}
}
/*write fields*/
template <class T,ENTITY E> 
void MeshField<T,E>::write(Int step) {
	using namespace Mesh;

	/*open*/
	std::stringstream path;
	path << fName << step;
	std::ofstream of(path.str().c_str());

	/*size*/
	of << "size " << sizeof(T) / sizeof(Scalar) << std::endl;

	/*internal field*/
	of << gBCellsStart << std::endl;
	of << "{" << std::endl;
	for(Int i = 0;i < gBCellsStart;i++)
		of << (*this)[i] << std::endl;
	of << "}" << std::endl;

	/*boundary field*/
	BasicBCondition* bbc;
	BCondition<T>* bc;
	for(Int i = 0;i < AllBConditions.size();i++) {
		bbc = AllBConditions[i];
		if(bbc->fIndex == this->fIndex) {
			bc = static_cast<BCondition<T>*> (bbc);
			of << *bc << std::endl;
		}
	}

	/*close*/
	of.close();
}
template <class T,ENTITY E> 
void MeshField<T,E>::writeAll(Int step) {
	using namespace Mesh;
	MeshField<T,E>* pf;
	for(typename std::list<MeshField<T,E>*>::iterator it = fields_.begin();it != fields_.end();++it) {
		pf = *it;
		if(pf->access & WRITE)
			pf->write(step);
	}
}
/*write in vtk format*/
template <class T,ENTITY E> 
int MeshField<T,E>::count_writable() {
	int count = 0;
	for(typename std::list<MeshField<T,CELL>*>::iterator it = fields_.begin();it != fields_.end();++it) {
		if((*it)->access & WRITE)
			count++;
	}
	return count;
}

template <class T,ENTITY E> 
void MeshField<T,E>::write_vtk(std::ostream& os,bool vertex) {
	MeshField<T,CELL>* pf;
	if(!vertex) {
		for(typename std::list<MeshField<T,CELL>*>::iterator it = fields_.begin();it != fields_.end();++it) {
			pf = *it;
			if(pf->access & WRITE) {
				os << pf->fName <<" "<< TYPE_SIZE <<" "<< Mesh::gBCellsStart << " float" << std::endl;
				for(Int i = 0;i < Mesh::gBCellsStart;i++)
					os << (*pf)[i] << std::endl;
				os << std::endl;
			}
		}
	} else {
		MeshField<T,VERTEX> vf;
		for(typename std::list<MeshField<T,CELL>*>::iterator it = fields_.begin();it != fields_.end();++it) {
			pf = *it;
			if(pf->access & WRITE) {
				vf = cds(cds(*pf));
				os << pf->fName <<" "<< TYPE_SIZE <<" "<< vf.size() << " float" << std::endl;
				for(Int i = 0;i < vf.size();i++)
					os << vf[i] << std::endl;
				os << std::endl;
			}
		}
	}
}
/*IO*/
template <class type,ENTITY entity> 
std::ostream& operator << (std::ostream& os, const MeshField<type,entity>& p) {
	for(Int i = 0;i < p.size();i++)
		os << p[i] << std::endl;
	os << std::endl;
	return os;
}
template <class type,ENTITY entity> 
std::istream& operator >> (std::istream& is, MeshField<type,entity>& p) {
	for(Int i = 0;i < p.size();i++)
		is >> p[i];
	return is;
}
/*********************************************************************************
 *                      matrix class defined on mesh                             
 *********************************************************************************/
template <class type> 
struct MeshMatrix {
	MeshField<type,CELL>* cF;
	MeshField<type,CELL> Su;
	ScalarCellField ap;
	ScalarFacetField an[2];
	Int flags;
	enum FLAG {
		SYMMETRIC = 1
	};
	/*c'tors*/
	MeshMatrix() {
		cF = 0;
		flags = 0;
	}
	MeshMatrix(const MeshMatrix& p) {
		cF = p.cF;
		flags = p.flags;
		ap = p.ap;
		an[0] = p.an[0];
		an[1] = p.an[1];
		Su = p.Su;
	}
	MeshMatrix(const MeshField<type,CELL>& p) {
		cF = 0;
		flags = SYMMETRIC;
		ap = Scalar(0);
		an[0] = Scalar(0);
		an[1] = Scalar(0);
		Su = p;
	}
    /*operators*/
	MeshMatrix operator - () {
		MeshMatrix r;
		r.cF = cF;
		r.flags = flags;
		r.ap = -ap;
		r.an[0] = -an[0];
		r.an[1] = -an[1];
		r.Su = -Su;
		return r;
	}
	MeshMatrix& operator = (const MeshMatrix& q) {
		cF = q.cF;
		flags = q.flags;
		ap = q.ap;
		an[0] = q.an[0];
		an[1] = q.an[1];
		Su = q.Su;
		return *this;
	}
	MeshMatrix& operator += (const MeshMatrix& q) {
		flags &= q.flags;
		ap += q.ap;
		an[0] += q.an[0];
		an[1] += q.an[1];
		Su += q.Su;
		return *this;
	}
	MeshMatrix& operator -= (const MeshMatrix& q) {
		flags &= q.flags;
		ap -= q.ap;
		an[0] -= q.an[0];
		an[1] -= q.an[1];
		Su -= q.Su;
		return *this;
	}
	MeshMatrix& operator *= (const Scalar& q) {
		ap *= q;
		an[0] *= q;
		an[1] *= q;
		Su *= q;
		return *this;
	}
	MeshMatrix& operator /= (const Scalar& q) {
		ap /= q;
		an[0] /= q;
		an[1] /= q;
		Su /= q;
		return *this;
	}
	/*binary ops*/
	Operator(MeshMatrix,+);
	Operator(MeshMatrix,-);
	/*is equal to*/
	friend MeshMatrix operator == (const MeshMatrix& p,const MeshMatrix& q) {
		MeshMatrix r = p;
		r -= q;
		return r;
	}
	/*relax*/
	void Relax(Scalar UR) {
		ap /= UR;
		Su += (*cF) * ap * (1 - UR);
	}
	/*Fix*/
	void Fix(Int c,type value) {
		/*diagonal fix*/
		ap[c] = 10e30;
		Su[c] = value * 10e30;
		/*asymmetric fix*
		using namespace Mesh;
		Cell& mc = gCells[c];
		Int f;
		for(Int i = 0;i < mc.size();i++) {
			f = mc[i];
			if(gFO[f] == c) an[1][f] = 0;
			else an[0][f] = 0;
		}
		Su[c] = ap[c] * value;
		/*symmetric fix*
		using namespace Mesh;
		Cell& mc = gCells[c];
		Int f;
		for(Int i = 0;i < mc.size();i++) {
			f = mc[i];
			if(gFO[f] == c)
				Su[gFN[f]] += an[0][f] * value;
			else
				Su[gFO[f]] += an[1][f] * value;
			an[1][f] = 0;
			an[0][f] = 0;
		}
		Su[c] = ap[c] * value;
		/*end*/
	}
	/*Fix near wall cell values*/
	void FixNearWallValues() {
		using namespace Mesh;
		BasicBCondition* bbc;
		for(Int d = 0;d < AllBConditions.size();d++) {
			bbc = AllBConditions[d];
			if(bbc->isWall && (bbc->fIndex == cF->fIndex)) {
				IntVector& wall_faces = *bbc->bdry;
				if(wall_faces.size()) {
					Int f,c1;
					for(Int i = 0;i < wall_faces.size();i++) {
						f = wall_faces[i];
						c1 = gFO[f];
						Fix(c1,(*cF)[c1]);
					}
				}
			}
		}
	}
};
/* ***************************************
 * Implicit boundary conditions
 * ***************************************/
 template <class T> 
 void applyImplicitBCs(const MeshMatrix<T>& M) {
	 using namespace Mesh;
	 MeshField<T,CELL>& cF = *M.cF;
	 BasicBCondition* bbc;
	 BCondition<T>* bc;

	 /*boundary conditions*/
	 for(Int i = 0;i < AllBConditions.size();i++) {
		 bbc = AllBConditions[i];
		 if(bbc->fIndex == cF.fIndex) {
			 if(bbc->cIndex == NEUMANN ||
				 bbc->cIndex == SYMMETRY)
				 ;
			 else continue;

			 bc = static_cast<BCondition<T>*> (bbc);
			 Int sz = bc->bdry->size();
			 if(sz == 0) continue;

			 for(Int j = 0;j < sz;j++) {
				 Int k = (*bc->bdry)[j];
				 Int c1 = gFO[k];
				 Int c2 = gFN[k];
				 if(bc->cIndex == NEUMANN) {
					 Vector dv = cC[c2] - cC[c1];
					 M.ap[c1] -= M.an[1][k];
					 M.Su[c1] += M.an[1][k] * (bc->value * mag(dv));
					 cF[c2] = 0;
				 } else if(bc->cIndex == ROBIN) {
					 Vector dv = cC[c2] - cC[c1];
					 M.ap[c1] -= (1 - bc->shape) * M.an[1][k];
					 M.Su[c1] += M.an[1][k] * (bc->shape * bc->value + 
						 (1 - bc->shape) * bc->tvalue * mag(dv));
					 cF[c2] = 0;
				 } else if(bc->cIndex == SYMMETRY) {
					 M.ap[c1] -= M.an[1][k];
					 M.Su[c1] += M.an[1][k] * (sym(cF[c1],fN[k]) - cF[c1]);
					 cF[c2] = 0;
				 }
			 }
		 }
	 }
 }
/* ***************************************
 * Explicit boundary conditions
 * **************************************/
template<class T,ENTITY E>
void updateExplicitBCs(const MeshField<T,E>& cF,
							  bool update_ghost = false,
							  bool update_fixed = false
							  ) {
	using namespace Mesh;
	BasicBCondition* bbc;
	BCondition<T>* bc;
	Scalar z = Scalar(0),
		   zmin = Scalar(0),
		   zmax = Scalar(0),
		   zR = Scalar(0);
	Vector C(0);
	bool first = false;

	/*boundary conditions*/
	for(Int i = 0;i < AllBConditions.size();i++) {
		bbc = AllBConditions[i];
		if(bbc->fIndex == cF.fIndex) {
			if(bbc->cIndex == GHOST) 
				continue;

			bc = static_cast<BCondition<T>*> (bbc);
			Int sz = bc->bdry->size();
			if(sz == 0) continue;

			if(update_fixed) {
				if(bc->cIndex == DIRICHLET || 
					bc->cIndex == POWER || 
					bc->cIndex == LOG || 
					bc->cIndex == PARABOLIC ||
					bc->cIndex == INVERSE
					) {
						Int ci,j;
						Scalar r;
						if(bc->zMax > 0) {
							zmin = bc->zMin;
							zmax = bc->zMax;
							zR = zmax - zmin;
						} else {
							zmin = Scalar(10e30);
							zmax = -Scalar(10e30);
							C = Vector(0);
							for(j = 0;j < sz;j++) {
								Facet& f = gFacets[j];
								for(Int k = 0;k < f.size();k++) {
									z = (vC[f[k]] & bc->dir);
									if(z < zmin) 
										zmin = z;
									if(z > zmax) 
										zmax = z;
								}
								C += fC[j];
							}
							C /= Scalar(sz);
							zR = zmax - zmin;

							if(bc->cIndex == PARABOLIC) {
								ci = gFN[(*bc->bdry)[0]];
								zR = magSq(cC[ci] - C);
								for(j = 1;j < sz;j++) {
									ci = gFN[(*bc->bdry)[0]];
									r = magSq(cC[ci] - C);
									if(r < zR) zR = r;
								}
							}
						}
				}
				first = (bc->fixed.size() == 0);
				if(first)
					bc->fixed.resize(sz);
			}
			for(Int j = 0;j < sz;j++) {
				Int k = (*bc->bdry)[j];
				Int c1 = gFO[k];
				Int c2 = gFN[k];
				if(bc->cIndex == NEUMANN) {
					Vector dv = cC[c2] - cC[c1];
					cF[c2] = cF[c1] + bc->value * mag(dv);
				} else if(bc->cIndex == ROBIN) {
					Vector dv = cC[c2] - cC[c1];
					cF[c2] = bc->shape * bc->value + 
						(1 - bc->shape) * (cF[c1] + bc->tvalue * mag(dv));
				} else if(bc->cIndex == SYMMETRY) {
					cF[c2] = sym(cF[c1],fN[k]);
				} else if(bc->cIndex == CYCLIC) {
					Int c22;
					if(j < sz / 2) 
						c22 = gFO[(*bc->bdry)[j + sz/2]];
					else
						c22 = gFO[(*bc->bdry)[j - sz/2]];
					cF[c2] = cF[c22];
				} else { 
					if(update_fixed) {
						T v(0);
						z = (cC[c2] & bc->dir) - zmin;
						if(bc->cIndex == DIRICHLET) {
							v = bc->value;
						} else if(bc->cIndex == POWER) {
							if(z < 0) z = 0;
							if(z > zR) v = bc->value;
							else v = bc->value * pow(z / zR,bc->shape);
						} else if(bc->cIndex == LOG) {
							if(z < 0) z = 0;
							if(z > zR) v = bc->value;
							else v = bc->value * (log(1 + z / bc->shape) / log(1 + zR / bc->shape));
						} else if(bc->cIndex == PARABOLIC) {
							z = magSq(cC[c2] - C);
							v = bc->value * (z / zR);
						} else if(bc->cIndex == INVERSE) {
							v = bc->value / (z + bc->shape);
						}
						if(!first && !equal(mag(bc->tvalue),0)) { 
							T meanTI = v * (bc->tvalue * pow (z / zR,-bc->tshape));
							Scalar rFactor = 4 * ((rand() / Scalar(RAND_MAX)) - 0.5);
							v += ((cF[c2] - v) * 0.9 + (meanTI * rFactor) * 0.1);
						} 
						bc->fixed[j] = cF[c2] = v;
					} else {
						cF[c2] = bc->fixed[j];
					}
				}
			}
		}
	}
	/*ghost cells*/
	if(update_ghost && gInterMesh.size()) {
		exchange_ghost(&cF[0]);
	}
}
/* ***************************************
 * Fill boundary from internal values
 * **************************************/
template<class T,ENTITY E>
void fillBCs(const MeshField<T,E>& cF,
			 bool update_ghost = false) {
	/*neumann update*/
	using namespace Mesh;
	for(Int i = gBCellsStart;i < cF.size();i++)
		cF[i] = cF[gFO[gCells[i][0]]];
	/*ghost cells*/
	if(update_ghost && gInterMesh.size()) {
		exchange_ghost(&cF[0]);
	}
}
/*************************************
 * Exchange ghost cell information
 *************************************/
template <class T> 
void exchange_ghost(T* P) {
	using namespace Mesh;
	Int i,j;
	/*blocked exchange*/
	if(Controls::ghost_exchange == Controls::BLOCKED) {
		MeshField<T,CELL> buffer;
		for(i = 0;i < gInterMesh.size();i++) {
			interBoundary& b = gInterMesh[i];
			IntVector& f = *(b.f);
			if(b.from < b.to) {
				//send
				for(j = 0;j < f.size();j++)
					buffer[j] = P[gFO[f[j]]];
				MP::send(&buffer[0],f.size(),b.to,MP::FIELD);
				//recieve
				MP::recieve(&buffer[0],f.size(),b.to,MP::FIELD);
				for(j = 0;j < f.size();j++)
					P[gFN[f[j]]] = buffer[j];
			} else {
				//recieve
				MP::recieve(&buffer[0],f.size(),b.to,MP::FIELD);
				for(j = 0;j < f.size();j++)
					P[gFN[f[j]]] = buffer[j];
				//send 
				for(j = 0;j < f.size();j++)
					buffer[j] = P[gFO[f[j]]];
				MP::send(&buffer[0],f.size(),b.to,MP::FIELD);
			}
		}
    /*Asynchronous exchange*/
	} else {
		MeshField<T,CELL> sendbuf,recvbuf;
		std::vector<MP::REQUEST> request(2 * gInterMesh.size(),0);
		Int rcount = 0;
		//fill send buffer
		for(i = 0;i < gInterMesh.size();i++) {
			interBoundary& b = gInterMesh[i];
			IntVector& f = *(b.f);
			for(j = 0;j < f.size();j++) 
				sendbuf[b.buffer_index + j] = P[gFO[f[j]]];
		}

		for(i = 0;i < gInterMesh.size();i++) {
			interBoundary& b = gInterMesh[i];
			//non-blocking send/recive
			MP::isend(&sendbuf[b.buffer_index],b.f->size(),
				b.to,MP::FIELD,&request[rcount]);
			rcount++;
			MP::irecieve(&recvbuf[b.buffer_index],b.f->size(),
				b.to,MP::FIELD,&request[rcount]);
			rcount++;
		}
		//wait
		MP::waitall(rcount,&request[0]);
		//recieve buffer
		for(i = 0;i < gInterMesh.size();i++) {
			interBoundary& b = gInterMesh[i];
			IntVector& f = *(b.f);
			for(j = 0;j < f.size();j++)
				P[gFN[f[j]]] = recvbuf[b.buffer_index + j];
		}
	}
	/*end*/
}

/*IO*/
template <class type> 
std::ostream& operator << (std::ostream& os, const MeshMatrix<type>& p) {
	os << p.ap << std::endl << std::endl;
	os << p.an[0] << std::endl << std::endl;
	os << p.an[1] << std::endl << std::endl;
	os << p.Su << std::endl << std::endl;
	return os;
}
template <class type> 
std::istream& operator >> (std::istream& is, MeshMatrix<type>& p) {
	is >> p.ap;
	is >> p.an[0];
	is >> p.an[1];
	is >> p.Su;
	return is;
}
/*typedefs*/
typedef MeshMatrix<Scalar>  ScalarMeshMatrix;
typedef MeshMatrix<Vector>  VectorMeshMatrix;
typedef MeshMatrix<Tensor>  TensorMeshMatrix;
typedef MeshMatrix<STensor> STensorMeshMatrix;

/* *******************************
 * matrix - vector product p * q
 * *******************************/
template <class T> 
MeshField<T,CELL> operator * (const MeshMatrix<T>& p,const MeshField<T,CELL>& q) {
	using namespace Mesh;
	MeshField<T,CELL> r;
	Int c1,c2;
	r = q * p.ap;
	for(Int f = 0;f < gFacets.size();f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] -= q[c2] * p.an[1][f];
		r[c2] -= q[c1] * p.an[0][f];
	}
	return r;
}
/*matrix transopose - vector product pT * q */
template <class T> 
MeshField<T,CELL> operator ^ (const MeshMatrix<T>& p,const MeshField<T,CELL>& q) {
	using namespace Mesh;
	MeshField<T,CELL> r;
	Int c1,c2;
	r = q * p.ap;
	for(Int f = 0;f < gFacets.size();f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c2] -= q[c1] * p.an[1][f];
		r[c1] -= q[c2] * p.an[0][f];
	}
	return r;
}
/* calculate RHS sum */
template <class T> 
MeshField<T,CELL> getRHS(const MeshMatrix<T>& p) {
	using namespace Mesh;
	MeshField<T,CELL> r;
	Int c1,c2;
	r = p.Su;
	for(Int f = 0;f < gFacets.size();f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] += (*p.cF)[c2] * p.an[1][f];
		r[c2] += (*p.cF)[c1] * p.an[0][f];
	}
	return r;
}

/* ********************************
 * Interpolate field operations
 * *******************************/
/*central difference*/
template<class type>
MeshField<type,FACET> cds(const MeshField<type,CELL>& cF) {
	using namespace Mesh;
	MeshField<type,FACET> fF;
	for(Int i = 0;i < fF.size();i++) {
		fF[i] =  (cF[gFO[i]] * (fI[i])) + (cF[gFN[i]] * (1 - fI[i]));
	}
	return fF;
}
/*upwind*/
template<class type>
MeshField<type,FACET> uds(const MeshField<type,CELL>& cF,const ScalarFacetField& flux) {
	using namespace Mesh;
	MeshField<type,FACET> fF;
	for(Int i = 0;i < fF.size();i++) {
		if(flux[i] >= 0) fF[i] = cF[gFO[i]];
		else fF[i] = cF[gFN[i]];
	}
	return fF;
}
/*facet data to vertex data */
template<class type>
MeshField<type,VERTEX> cds(const MeshField<type,FACET>& fF) {
	using namespace Mesh;
	std::vector<Scalar> cnt;
	MeshField<type,VERTEX> vF;
	cnt.assign(vF.size(),0);
	
	vF = type(0);
	for(Int i = 0;i < fF.size();i++) {
		Facet& f = gFacets[i];
		if(gFN[i] < gBCellsStart) {
			for(Int j = 0;j < f.size();j++) {
				vF[f[j]] += fF[i];
				cnt[f[j]]++;
			}
		} else {
			for(Int j = 0;j < f.size();j++) {
				vF[f[j]] += Scalar(10e30) * fF[i];
				cnt[f[j]] += Scalar(10e30);
			}
		}
	}
	for(Int i = 0;i < vF.size();i++) {
		vF[i] /= Scalar(cnt[i]);
		if(mag(vF[i]) < Constants::MachineEpsilon) 
			vF[i] = type(0);
	}
	return vF;
}
/* **************************
 * Integrate field operation
 * **************************/
template<class type>
MeshField<type,CELL> sum(const MeshField<type,FACET>& fF) {
	using namespace Mesh;
	MeshField<type,CELL> cF;
	cF = type(0);
	for(Int i = 0;i < fF.size();i++) {
		cF[gFO[i]] += fF[i];
		cF[gFN[i]] -= fF[i];
	}
	return cF;
}
/**********************************************************************
 * Gradient field operation.
 *   gradV(p) = Sum_f ( fN * p)
 *   grad(p) = gradV(p) / V
 * gradV(p) is integrated over the volume so it can be used directly in 
 * finite volume equations just like div,lap,ddt,src etc...
 * grad(p) returns per-unit volume gradient at the centre.
 **********************************************************************/

/*Explicit*/
inline VectorCellField gradV(const ScalarFacetField& p) {
	return sum(mul(Mesh::fN,p));
}
inline VectorCellField gradV(const ScalarCellField& p) {
	return gradV(cds(p));
}
inline TensorCellField gradV(const VectorFacetField& p) {
	return sum(mul(Mesh::fN,p));
}
inline TensorCellField gradV(const VectorCellField& p) {
	return gradV(cds(p));
}

/*Explicit*/
inline VectorCellField grad(const ScalarFacetField& p) {
	VectorCellField f = gradV(p) / Mesh::cV;
	fillBCs(f,true);
	return f;
}
inline VectorCellField grad(const ScalarCellField& p) {
	return grad(cds(p));
}
inline TensorCellField grad(const VectorFacetField& p) {
	TensorCellField f = gradV(p) / Mesh::cV;
	fillBCs(f,true);
	return f;
}
inline TensorCellField grad(const VectorCellField& p) {
	return grad(cds(p));
}

/* *********************************************
 * Laplacian field operation
 * ********************************************/

/*Implicit*/
template<class type>
MeshMatrix<type> lap(MeshField<type,CELL>& cF,const ScalarFacetField& mu) {
	using namespace Controls;
	using namespace Mesh;
	MeshMatrix<type> m;
	VectorFacetField K;
	Vector dv;
	Int c1,c2;
	Scalar D = 0;
	/*clear*/
	m.cF = &cF;
	m.flags |= m.SYMMETRIC;
	m.Su = type(0);
	m.ap = Scalar(0);
	for(Int i = 0;i < mu.size();i++) {
		c1 = gFO[i];
		c2 = gFN[i];
		dv = cC[c2] - cC[c1];
		/*diffusivity coefficient*/
		if(nonortho_scheme == NO_CORRECTION) {
			D = mag(fN[i]) / mag(dv);
		} else {
			if(nonortho_scheme == OVER_RELAXED) {
				D = ((fN[i] & fN[i]) / (fN[i] & dv));
			} else if(nonortho_scheme == MINIMUM) {
				D = ((fN[i] & dv) / (dv & dv));
			} else if(nonortho_scheme == ORTHOGONAL) {
				D = sqrt((fN[i] & fN[i]) / (dv & dv));
			}
			K[i] = fN[i] - D * dv;
		}
		/*coefficients*/
		m.an[0][i] = D * mu[i];
		m.an[1][i] = D * mu[i];
		m.ap[c1]  += m.an[0][i];
		m.ap[c2]  += m.an[1][i];
	}
	/*non-orthogonality handled through deferred correction*/
	if(nonortho_scheme != NO_CORRECTION) {
		MeshField<type,FACET> r = dot(cds(grad(cF)),K);
		type res;
		for(Int i = 0;i < mu.size();i++) {
			c1 = gFO[i];
			c2 = gFN[i];
			res = m.an[0][i] * (cF[c2] - cF[c1]);
			if(mag(r[i]) > Scalar(0.5) * mag(res)) 
				r[i] = Scalar(0.5) * res;
		}
		m.Su = sum(r);
	}
	/*end*/
	return m;
}

template<class type>
inline MeshMatrix<type> lap(MeshField<type,CELL>& cF,const ScalarCellField& mu) {
	return lap(cF,cds(mu));
}

/* ***************************************************
 * Divergence field operation
 * ***************************************************/ 
/*face flux*/
inline ScalarFacetField flx(const VectorFacetField& p) {
	return dot(p,Mesh::fN);
}
inline ScalarFacetField flx(const VectorCellField& p) {
	return flx(cds(p));
}
inline VectorFacetField flx(const TensorFacetField& p) {
	return dot(p,Mesh::fN);
}
inline VectorFacetField flx(const TensorCellField& p) {
	return flx(cds(p));
}
/* Explicit */
inline ScalarCellField div(const VectorFacetField& p) {
	return sum(flx(p));
}
inline ScalarCellField div(const VectorCellField& p) {
	return sum(flx(p));
}
inline VectorCellField div(const TensorFacetField& p) {
	return sum(flx(p));
}
inline VectorCellField div(const TensorCellField& p) {
	return sum(flx(p));
}
/* Implicit */
template<class type>
MeshMatrix<type> div(MeshField<type,CELL>& cF,const ScalarFacetField& flux,const ScalarFacetField& mu) {
	using namespace Controls;
	using namespace Mesh;
	MeshMatrix<type> m;
	Scalar F,G;
	m.cF = &cF;
	m.flags = 0;
	m.Su = type(0);
	m.ap = Scalar(0);

	/*Implicit convection schemes*/
	bool isImplicit = (
		convection_scheme == CDS ||
		convection_scheme == UDS ||
		convection_scheme == BLENDED ||
		convection_scheme == HYBRID );

	if(isImplicit) {
		ScalarFacetField gamma;
		if(convection_scheme == CDS) 
			gamma = Scalar(1);
		else if(convection_scheme == UDS) 
			gamma = Scalar(0);
		else if(convection_scheme == BLENDED) 
			gamma = Scalar(blend_factor);
		else if(convection_scheme == HYBRID) {
			Scalar D;
			Vector dv;
			for(Int j = 0;j < gFacets.size();j++) {
				/*calc D - uncorrected */
				dv = cC[gFN[j]] - cC[gFO[j]];
				D = (mag(fN[j]) / mag(dv)) * mu[j];
				/*compare F and D */
				F = flux[j];
				if(F < 0) {
					if(-F * fI[j] > D) gamma[j] = 0;
					else gamma[j] = 1;
				} else {
					if(F * (1 - fI[j]) > D) gamma[j] = 0;
					else gamma[j] = 1;
				}
			}
		}
		for(Int i = 0;i < flux.size();i++) {
			F = flux[i];
			G = gamma[i];
			m.an[0][i] = ((G) * (-F * (  fI[i]  )) + (1 - G) * (-Util::max( F,0)));
			m.an[1][i] = ((G) * ( F * (1 - fI[i])) + (1 - G) * (-Util::max(-F,0)));
			m.ap[gFO[i]] += m.an[0][i];
			m.ap[gFN[i]] += m.an[1][i];
		}
	/*deferred correction*/
	} else {
		for(Int i = 0;i < flux.size();i++) {
			F = flux[i];
			m.an[0][i] = -Util::max( F,0);
			m.an[1][i] = -Util::max(-F,0);
			m.ap[gFO[i]] += m.an[0][i];
			m.ap[gFN[i]] += m.an[1][i];
		}

		MeshField<type,FACET> corr;
		if(convection_scheme == CDSS) {
			corr = cds(cF) - uds(cF,flux);
		} else if(convection_scheme == LUD) {
			VectorFacetField R = fC - uds(cC,flux);
			corr = dot(uds(grad(cF),flux),R);
		} else if(convection_scheme == MUSCL) {
			VectorFacetField R = fC - uds(cC,flux);
			corr  = (  blend_factor  ) * (cds(cF) - uds(cF,flux));
			corr += (1 - blend_factor) * (dot(uds(grad(cF),flux),R));
		} else {
			/*
			TVD schemes
			~~~~~~~~~~~
			Reference:
				M.S Darwish and F Moukalled "TVD schemes for unstructured grids"
				Versteeg and Malaskara
			Description:
				phi = phiU + psi(r) * [(phiD - phiC) * (1 - fi)]
			Schemes
				psi(r) = 0 =>UDS
				psi(r) = 1 =>CDS
			R is calculated as ratio of upwind and downwind gradient
				r = phiDC / phiCU
			Further modification to unstructured grid to better fit LUD scheme
				r = (phiDC / phiCU) * (fi / (1 - fi))
			*/
			/*calculate r*/
			MeshField<type,FACET> q,r,phiDC,phiCU;
			ScalarFacetField uFI;
			{
				ScalarFacetField nflux = Scalar(0)-flux;
				phiDC = uds(cF,nflux) - uds(cF,flux);
				for(Int i = 0;i < phiDC.size();i++) {
					if(flux[i] >= 0) G = fI[i];
					else G = 1 - fI[i];
					uFI[i] = G;
				}
				/*Bruner's or Darwish way of calculating r*/
				if(TVDbruner) {
					VectorFacetField R = fC - uds(cC,flux);
					phiCU = 2 * (dot(uds(grad(cF),flux),R));
				} else {
					VectorFacetField R = uds(cC,nflux) - uds(cC,flux);
					phiCU = 2 * (dot(uds(grad(cF),flux),R)) - phiDC;
				}
				/*end*/
			}
			r = (phiCU / phiDC) * (uFI / (1 - uFI));
			for(Int i = 0;i < phiDC.size();i++) {
				if(equal(phiDC[i] * (1 - uFI[i]),type(0)))
					r[i] = type(0);
			}
			/*TVD schemes*/
			if(convection_scheme == VANLEER) {
				q = (r+fabs(r)) / (1+r);
			} else if(convection_scheme == VANALBADA) {
				q = (r+r*r) / (1+r*r);
			} else if(convection_scheme == MINMOD) {
				q = max(type(0),min(r,type(1)));
			} else if(convection_scheme == SUPERBEE) {
				q = max(min(r,type(2)),min(2*r,type(1)));
				q = max(q,type(0));
			} else if(convection_scheme == SWEBY) {
				Scalar beta = 2;
				q = max(min(r,type(beta)),min(beta*r,type(1)));
				q = max(q,type(0));
			} else if(convection_scheme == QUICKL) {
				q = min(2*r,(3+r)/4);
				q = min(q,type(2));
				q = max(q,type(0));
			} else if(convection_scheme == UMIST) {
				q = min(2*r,(3+r)/4);
				q = min(q,(1+3*r)/4);
				q = min(q,type(2));
				q = max(q,type(0));
			} else if(convection_scheme == QUICK) {
				q = (3+r)/4;
			} else if(convection_scheme == DDS) {
				q = 2;
			} else if(convection_scheme == FROMM) {
				q = (1+r)/2;
			}
			corr = q * phiDC * (1 - uFI);
			/*end*/
		}
		m.Su = sum(flux * corr);
	}
	return m;
}

template<class type,ENTITY E>
inline MeshMatrix<type> div(MeshField<type,CELL>& cF,const MeshField<Vector,E>& rhoU,const ScalarFacetField& mu) {
	return div(cF,div(rhoU),mu);
}

/* *******************************
 * Temporal derivative
 * *******************************/
template<class type>
MeshMatrix<type> ddt(MeshField<type,CELL>& cF,const ScalarCellField& rho) {
	MeshMatrix<type> m;
	m.cF = &cF;
	m.flags |= m.SYMMETRIC;
	m.ap = (Mesh::cV * rho) / -Controls::dt;
	m.an[0] = Scalar(0);
	m.an[1] = Scalar(0);
	m.Su = cF * m.ap;
	return m;
}
/* *******************************
 * Linearized source term
 * *******************************/
template<class type>
MeshMatrix<type> src(MeshField<type,CELL>& cF,const ScalarCellField& Sc,const ScalarCellField Sp) {
	MeshMatrix<type> m;
	m.cF = &cF;
	m.flags |= m.SYMMETRIC;
	m.ap = -(Sp * Mesh::cV);
	m.an[0] = Scalar(0);
	m.an[1] = Scalar(0);
	m.Su = (Sc * Mesh::cV);
	return m;
}

/* **************************************
 *   CSR - compressed sparse row format
 *       * Used for on GPU computation
 *       * Propably for AMG too
 * **************************************/
template <class T>
class CSRMatrix {
public:
	std::vector<Int>  rows;
	std::vector<Int>  cols;
	std::vector<Scalar> an;
	std::vector<Scalar> anT;
	std::vector<T> cF;
	std::vector<T> Su;
public:
	template <class T1>
	CSRMatrix(const MeshMatrix<T1>& A) {
		using namespace Mesh;
		const Int N  = A.ap.size();
		const Int NN = A.ap.size() + 
			           A.an[0].size() + 
					   A.an[1].size(); 
		register Int i,j,f;

		//resize
		cF.resize(N);
		Su.resize(N);
		rows.reserve(N + 1);
		cols.reserve(NN);
		an.reserve(NN);
		anT.reserve(NN);

        //source term
		for(i = 0;i < N;i++) {
			Su[i] = A.Su[i];
			cF[i] = (*A.cF)[i];
		}

		//fill matrix in CSR format
		//diagonal element is always at the start of a row
		Int cn = 0;
		for(i = 0;i < N;i++) {
			Cell& c = gCells[i];

			rows.push_back(cn);

			an.push_back(A.ap[i]);
			anT.push_back(A.ap[i]);
			cols.push_back(i);
			cn++;

			for(j = 0;j < c.size();j++) {
				f = c[j];
				if(i == gFO[f]) {
					an.push_back(A.an[1][f]);
					anT.push_back(A.an[0][f]);
					cols.push_back(gFN[f]);
					cn++;
				} else {
					an.push_back(A.an[0][f]);
					anT.push_back(A.an[1][f]);
					cols.push_back(gFO[f]);
					cn++;
				}
			}
		}
		//extra row
		rows.push_back(cn);
	}
	/*IO*/
	friend std::ostream& operator << (std::ostream& os, const CSRMatrix& p) {
		os << p.rows << std::endl;
		os << p.cols << std::endl;
		os << p.an << std::endl;
		os << p.Su << std::endl;
		return os;
	}
	friend std::istream& operator >> (std::istream& is, CSRMatrix& p) {
		is >> p.rows;
		is >> p.cols;
		is >> p.an;
		is >> p.Su;
		return is;
	}
	/*end*/
};
/* ********************
 *        End
 * ********************/
#endif
