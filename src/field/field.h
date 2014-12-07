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
		NONE,MINIMUM, ORTHOGONAL, OVER_RELAXED
	};
	enum TimeScheme {
		EULER, SECOND_ORDER
	};
	enum Solvers {
		JACOBI, SOR, PCG
	};
	enum Preconditioners {
		NOPR,DIAG,SSOR,DILU
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
	extern TimeScheme time_scheme;
	extern Solvers Solver; 
	extern Preconditioners Preconditioner;
	extern CommMethod parallel_method;
	extern State state;

	extern Scalar SOR_omega;
	extern Scalar tolerance;
	extern Scalar blend_factor;
	extern Scalar implicit_factor;
	extern Scalar dt;
	extern Int runge_kutta;
	
	extern Int max_iterations;
	extern Int write_interval;
	extern Int start_step;
	extern Int end_step;
	extern Int n_deferred;
	extern Int save_average;

	extern Vector gravity;
}

namespace DG {
    extern Int Nop[3];
    extern Int NP,NPMAT,NPF;
};
namespace Mesh {
	extern IntVector  probeCells;
	extern Int  gBCSfield; 
	extern Int  gBCSIfield;
	extern Int  gBFSfield;
};
enum ACCESS {
	NO = 0, READ = 1, WRITE = 2,READWRITE = 3,STOREPREV = 4
};

/* *****************************************************************************
 *                    Field variables defined on mesh                          
 * *****************************************************************************/
template <class type,ENTITY entity> 
class MeshField {
private:
	type*        P;
	int          allocated;
	static Int   SIZE;
public:
	ACCESS       access;
	Int          fIndex;
	std::string  fName;

	/*common*/
	static const Int TYPE_SIZE = sizeof(type) / sizeof(Scalar);
	static std::list<MeshField*> fields_;
	static std::list<type*> mem_;

    /*constructors*/
	MeshField(const char* str = "", ACCESS a = NO,bool recycle = true) : 
				P(0),allocated(0),access(a),fName(str) {
		construct(str,a,recycle);
	}
	MeshField(const MeshField& p) : allocated(0) {
		allocate(); 
		forEach(*this,i)
			P[i] = p[i];
	}
	MeshField(const type& p) : allocated(0) {
		allocate(); 
		forEach(*this,i)
			P[i] = p;
	}
	explicit MeshField(const bool) : allocated(0) {
	}
	/*allocators*/
	void allocate(bool recycle = true) {
		if(!recycle || mem_.empty()) {
			switch(entity) {
				case CELL:   SIZE = Mesh::gCells.size() * DG::NP;	 break;
				case FACET:  SIZE = Mesh::gFacets.size() * DG::NPF;  break;
				case VERTEX: SIZE = Mesh::gVertices.size();          break;
				case CELLMAT:SIZE = Mesh::gCells.size() * DG::NPMAT; break;
			}
			Int sz = SIZE;
			if(entity == CELL) 
				sz += 1;
			else if(entity == CELLMAT) 
				sz += DG::NP;
			P = new type[sz];
		} else {
			P = mem_.front();
			mem_.pop_front();
		}
		allocated = 1;
	}
	void allocate(std::vector<type>& q) {
		SIZE = q.size();
		P = &q[0];
		allocated = 0;
	}
	void deallocate(bool recycle = true) {
		if(allocated) {
			allocated = 0;
			if(recycle) {
				mem_.push_front(P);
			} else {
				delete[] P;
				P = 0;
			}
			if(fIndex)
				fields_.remove(this);
		}
	}
	void construct(const char* str = "", ACCESS a = NO, bool recycle = true) {
		access = a;
		fName = str;
		if(Mesh::gCells.size())
			allocate(recycle);
		fIndex = Util::hash_function(str);
		if(fIndex)
			fields_.push_back(this);
	}
	/*d'tor re-cycles memory */
	~MeshField() {
		if(!MP::Terminated)
			deallocate();
	}

	/*static functions*/
	void readInternal(std::istream&);
	void read(Int step);
	void write(Int step);

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
		forEach(*this,i)
			r[i] = -P[i];
		return r;
	}
	friend MeshField<Scalar,entity> operator & (const MeshField& p,const MeshField& q) {
		MeshField<Scalar,entity> r;
		forEach(r,i)
			r[i] = p[i] & q[i];
		return r;
	}
	/*unrolled operations*/
#define Op($)													    \
	MeshField& operator $(const MeshField& q) {						\
		forEach(*this,i)											\
			P[i] $ q[i];											\
		return *this;												\
	}
#define SOp($)														\
	MeshField& operator $(const Scalar& q) {						\
		forEach(*this,i)											\
			P[i] $ q;												\
		return *this;												\
	}
#define Fp(name)													\
	friend MeshField name(const MeshField& p,const MeshField& s) {	\
		MeshField r;												\
		forEach(r,i)												\
			r[i] = name(p[i],s[i]);									\
		return r;													\
	}
#define Fp1(name)													\
	friend MeshField name(const MeshField& p,const Scalar& s) {		\
		MeshField r;												\
		forEach(r,i)												\
			r[i] = name(p[i],s);									\
		return r;													\
	}
#define Fp2(name)													\
	friend MeshField name(const MeshField& p) {						\
		MeshField r;												\
		forEach(r,i)												\
			r[i] = name(p[i]);										\
		return r;													\
	}
    /*define ops*/
	Op(=)
	Op(+=)
	Op(-=)
	Op(*=)
	Op(/=)
	SOp(=)
	SOp(+=)
	SOp(-=)
	SOp(*=)
	SOp(/=)
	Fp(sdiv)
	/*from math.h*/
	Fp2(acos)
	Fp2(asin)
	Fp2(atan)
	Fp(atan2)
	Fp2(ceil)
	Fp2(cos)
	Fp2(cosh)
	Fp2(exp)
	Fp2(fabs)
	Fp2(floor)
	Fp2(log)
	Fp2(log10)
	Fp1(pow)
	Fp2(sin)
	Fp2(sinh)
	Fp2(sqrt)
	Fp2(tan)
	Fp2(tanh)
	Fp(min)
	Fp(max)
	/*additional*/
	Fp2(unit)
#undef Op
#undef SOp
#undef Fp
#undef Fp1
#undef Fp2
	AddOperators(MeshField)
	AddScalarOperators(MeshField)
	/*friend ops*/
	friend MeshField<Scalar,entity> mag(const MeshField& p) {
		MeshField<Scalar,entity> r;
		forEach(r,i)
			r[i] = mag(p[i]);
		return r;
	}
	friend MeshField dev(const MeshField& p,const Scalar factor = 1.) {
		MeshField r;
		forEach(r,i)
			r[i] = dev(p[i],factor);
		return r;
	}
	friend MeshField hyd(const MeshField& p,const Scalar factor = 1.) {
		MeshField r;
		forEach(r,i)
			r[i] = hyd(p[i],factor);
		return r;
	}
	/*relax*/
	void Relax(const MeshField& po,Scalar UR) {
		forEach(*this,i) 
			P[i] = po[i] + (P[i] - po[i]) * UR;
	}
	/*read/write all fields*/
	static void readAll(Int step) {
		forEachIt(typename std::list<MeshField*>, fields_, it) {
			if((*it)->access & READ)
				(*it)->read(step);
		}
	}
	static void writeAll(Int step) {
		forEachIt(typename std::list<MeshField*>, fields_, it) {
			if((*it)->access & WRITE)
				(*it)->write(step);
		}
	} 
	static void removeAll() {
		fields_.clear();
		forEachIt(typename std::list<type*>,mem_,it)
			delete (*it);
		mem_.clear();
	} 
	static int count_writable() {
		int count = 0;
		forEachIt(typename std::list<MeshField*>, fields_, it) {
			if((*it)->access & WRITE)
				count++;
		}
		return count;
	}
	static void writeVtkCellAll(std::ostream& os) {
		MeshField<type,CELL>* pf;
		forEachIt(typename std::list<MeshField*>, fields_, it) {
			pf = *it;
			if(pf->access & WRITE) {
				os << pf->fName <<" "<< TYPE_SIZE <<" "
					<< Mesh::gBCSfield << " float" << std::endl;
				for(Int i = 0;i < Mesh::gBCSfield;i++)
					os << (*pf)[i] << std::endl;
				os << std::endl;
			}
		}
	}
	static void writeVtkVertexAll(std::ostream& os) {
		MeshField<type,VERTEX> vf;
		forEachIt(typename std::list<MeshField*>, fields_, it) {
			if((*it)->access & WRITE) {
				vf = cds(cds(*(*it)));
				os << (*it)->fName <<" "<< TYPE_SIZE <<" "
					<< vf.size() << " float" << std::endl;
				forEach(vf,i)
					os << vf[i] << std::endl;
				os << std::endl;
			}
		}
	}
	/*interpolation*/
	typedef std::list< MeshField<type,VERTEX> > vertexFieldsType;
	static vertexFieldsType* vf_fields_;
	static void interpolateVertexAll() {
		vf_fields_ = new vertexFieldsType;
		vf_fields_->clear();
		MeshField<type,VERTEX> vf;
		forEachIt(typename std::list<MeshField*>, fields_, it) {
			if((*it)->access & WRITE) {
				vf = cds(cds(*(*it)));
				vf_fields_->push_back(vf);
			}
		}
	}
	/*Store previous values*/
	MeshField* tstore;
	void initStore() {
		tstore = new MeshField[2];
		access = ACCESS(int(access) | STOREPREV);
		updateStore();
	}
	void updateStore() {
		tstore[1] = tstore[0];
		tstore[0] = *this;
	}
	/*Time history*/
	static std::vector<std::ofstream*> tseries;
	static std::vector<MeshField*> tavgs;
	static std::vector<MeshField*> tstds;
	
	static void initTimeSeries() {
		MeshField<type,CELL>* pf;
		forEachIt(typename std::list<MeshField*>, fields_, it) {
			pf = *it;
			if(pf->access & WRITE) {
				if(Mesh::probeCells.size()) {
					std::string name = pf->fName + "i";
					std::ofstream* of = new std::ofstream(name.c_str());
					tseries.push_back(of);
				}
				if(Controls::save_average) {
					std::string name;
					name = pf->fName + "avg";
					MeshField* avg = new MeshField(name.c_str(),READWRITE);
					tavgs.push_back(avg);
					name = pf->fName + "std";
					MeshField* std = new MeshField(name.c_str(),READWRITE);
					tstds.push_back(std);
				}
			}
		}
	}
	static void updateTimeSeries(int i) {
		int count = 0;
		MeshField<type,CELL>* pf;
		forEachIt(typename std::list<MeshField*>, fields_, it) {
			pf = *it;
			if(pf->access & WRITE) {
				if(Mesh::probeCells.size()) {
					std::ofstream& of = *tseries[count];
					of << i << " ";
					forEach(Mesh::probeCells,j) 
						of << (*pf)[Mesh::probeCells[j]] << " ";
					of << std::endl;
				}
				if(Controls::save_average) {
					MeshField& avg = *tavgs[count];
					avg += (*pf);
					MeshField& std = *tstds[count];
					std += (*pf) * (*pf);
					count++;
				}
			}
			if(pf->access & STOREPREV) {
				pf->updateStore();
			}
		}
	}
	/*IO*/
	friend std::ostream& operator << (std::ostream& os, const MeshField& p) {
		forEach(p,i)
			os << p[i] << std::endl;
		os << std::endl;
		return os;
	}
	friend std::istream& operator >> (std::istream& is, MeshField& p) {
		forEach(p,i)
			is >> p[i];
		return is;
	}
};
#define forEachCellField(X)	 { 		\
	ScalarCellField::X;    			\
    VectorCellField::X;    			\
	STensorCellField::X;   			\
	TensorCellField::X;    			\
}
#define forEachFacetField(X)	 { 	\
	ScalarFacetField::X;    		\
    VectorFacetField::X;    		\
	STensorFacetField::X;   		\
	TensorFacetField::X;    		\
}
#define forEachVertexField(X)	 { 	\
	ScalarVertexField::X;    		\
    VectorVertexField::X;    		\
	STensorVertexField::X;   		\
	TensorVertexField::X;    		\
}
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

/***********************************
 *  Specific tensor operations
 ***********************************/
/* Default operator overload for scalar fields*/
#define Op(name,F,S)																			\
	template<class T,ENTITY E>																	\
	MeshField<T,E> name(const MeshField<F,E>& p,const MeshField<S,E>& q) {						\
		MeshField<T,E> r;																		\
		forEach(r,i)																			\
			r[i] = name(p[i],q[i]);																\
		return r;																				\
	}
Op(operator *,Scalar,T)
Op(operator /,Scalar,T)
Op(operator *,T,Scalar)
Op(operator /,T,Scalar)
#undef Op
/*multiply*/
template <ENTITY E>
MeshField<Tensor,E> mul(const MeshField<Vector,E>& p,const MeshField<Vector,E>& q) {
	MeshField<Tensor,E> r;
	forEach(r,i)
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
	forEach(r,i) 
		r[i] = mul(p[i],q[i]);
	return r;
}
/*dot*/
template <ENTITY E,Int SIZE> 
MeshField<Vector,E> dot(const MeshField<TTensor<SIZE>,E>& p,const MeshField<Vector,E>& q) {
	MeshField<Vector,E> r;
	forEach(r,i)
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
	forEach(r,i) 
		r[i] = sym(p[i]);
	return r;
}
template <ENTITY E>
MeshField<Tensor,E> skw(const MeshField<Tensor,E>& p) {
	MeshField<Tensor,E> r;
	forEach(r,i) 
		r[i] = skw(p[i]);
	return r;
}
/*transpose*/
template <ENTITY E>
MeshField<Tensor,E> trn(const MeshField<Tensor,E>& p) {
	MeshField<Tensor,E> r;
	forEach(r,i) 
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

template <class T,ENTITY E>
std::vector<std::ofstream*> MeshField<T,E>::tseries;

template <class T,ENTITY E>
std::vector<MeshField<T,E>*> MeshField<T,E>::tavgs;

template <class T,ENTITY E>
std::vector<MeshField<T,E>*> MeshField<T,E>::tstds;

template <class T,ENTITY E> 
typename MeshField<T,E>::vertexFieldsType* MeshField<T,E>::vf_fields_;

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
	extern IntVector         gFO;
	extern IntVector         gFN; 
	
	bool   LoadMesh(Int = 0,bool = true, bool = true);
	void   initGeomMeshFields();
	void   calc_walldist(Int,Int = 1);
	void   write_fields(Int);
	void   read_fields(Int);
	void   remove_fields();
	Int    findNearestCell(const Vector& v);
	Int    findNearestFace(const Vector& v);
	void   getProbeCells(IntVector&);
	void   getProbeFaces(IntVector&);
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
		if(str == "uniform") {
			is >> value;
			*this = value;
		} else if(str == "elliptic") {
			Vector center,radius;
			T perterb;
			is >> value >> perterb >> center >> radius;
			VectorCellField cB = center;
			ScalarCellField R = mag((cC - cB) / radius);
			forEach(R,i) R[i] = min(1.0,R[i]);
			MeshField<T,E> val = value;
			val += MeshField<T,E>(perterb / 2) * 
					(MeshField<Scalar,E>(1.0) + cos(R * Constants::PI));
			*this = val;
		} else if(str == "gaussian") {
			Vector center,radius;
			T perterb;
			is >> value >> perterb >> center >> radius;
			VectorCellField cB = center;
			ScalarCellField R = mag((cC - cB) / radius);
			MeshField<T,E> val = value;
			val += MeshField<T,E>(perterb / 2) *  exp(-R*R);
			*this = val;
		} else if(str == "hydrostatic") {
			T p0;
			Scalar scale, expon;
			is >> p0 >> scale >> expon;
			Vector gu = unit(Controls::gravity);
			ScalarCellField gz = -dot(cC,VectorCellField(gu));
			*this = MeshField<T,E>(p0) * 
				pow(MeshField<Scalar,E>(1.0) - scale * gz, expon);
		} else
			*this = value;
	} else {
		char symbol;
		is >> size >> symbol;
		for(int i = 0;i < size;i++)
			is >> (*this)[i];
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
	if(is.fail())
		return;

	/*start reading*/
	std::cout << "Reading " << fName 
		 << step  << std::endl;
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
	applyExplicitBCs(*this,true,true);
}
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
	Int size = (SIZE == gCells.size() * DG::NP) ? 
	             gBCSfield : 
	             SIZE;
	of << size << std::endl;
	of << "{" << std::endl;
	for(Int i = 0;i < size;i++)
		of << (*this)[i] << std::endl;
	of << "}" << std::endl;

	/*boundary field*/
	BasicBCondition* bbc;
	BCondition<T>* bc;
	forEach(AllBConditions,i) {
		bbc = AllBConditions[i];
		if(bbc->fIndex == this->fIndex) {
			bc = static_cast<BCondition<T>*> (bbc);
			of << *bc << std::endl;
		}
	}
}

/* ********************
 *   DG
 * ********************/
 #include "dg.h"
 
/*********************************************************************************
 *                      matrix class defined on mesh                             
 *********************************************************************************/
template <class type> 
struct MeshMatrix {
	MeshField<type,CELL>*   cF;
	MeshField<type,CELL>    Su;
	MeshField<Scalar,CELL>  ap;
	MeshField<Scalar,FACET> an[2];
	MeshField<Scalar,CELLMAT> adg;
	Int flags;
	enum FLAG {
		SYMMETRIC = 1, DIAGONAL = 2
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
		adg = p.adg;
	}
	MeshMatrix(MeshField<type,CELL>* pcF) {
		cF = pcF;
		flags = (SYMMETRIC | DIAGONAL);
		ap = Scalar(0);
		an[0] = Scalar(0);
		an[1] = Scalar(0);
		Su = type(0);
		adg = Scalar(0);
	}
	MeshMatrix(const MeshField<type,CELL>& p) {
		cF = 0;
		flags = (SYMMETRIC | DIAGONAL);
		ap = Scalar(0);
		an[0] = Scalar(0);
		an[1] = Scalar(0);
		Su = p;
		adg = Scalar(0);
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
		r.adg = -adg;
		return r;
	}
	MeshMatrix& operator = (const MeshMatrix& q) {
		cF = q.cF;
		flags = q.flags;
		ap = q.ap;
		an[0] = q.an[0];
		an[1] = q.an[1];
		Su = q.Su;
		adg = q.adg;
		return *this;
	}
	MeshMatrix& operator += (const MeshMatrix& q) {
		flags &= q.flags;
		ap += q.ap;
		an[0] += q.an[0];
		an[1] += q.an[1];
		Su += q.Su;
		adg += q.adg;
		return *this;
	}
	MeshMatrix& operator -= (const MeshMatrix& q) {
		flags &= q.flags;
		ap -= q.ap;
		an[0] -= q.an[0];
		an[1] -= q.an[1];
		Su -= q.Su;
		adg -= q.adg;
		return *this;
	}
	MeshMatrix& operator *= (const Scalar& q) {
		ap *= q;
		an[0] *= q;
		an[1] *= q;
		Su *= q;
		adg *= q;
		return *this;
	}
	MeshMatrix& operator /= (const Scalar& q) {
		ap /= q;
		an[0] /= q;
		an[1] /= q;
		Su /= q;
		adg /= q;
		return *this;
	}
	/*binary ops*/
	AddOperator(MeshMatrix,+);
	AddOperator(MeshMatrix,-);
	AddScalarOperators(MeshMatrix)
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
	}
	/*IO*/
	friend std::ostream& operator << (std::ostream& os, const MeshMatrix& p) {
		os << p.ap << std::endl << std::endl;
		os << p.an[0] << std::endl << std::endl;
		os << p.an[1] << std::endl << std::endl;
		os << p.Su << std::endl << std::endl;
		os << p.adg << std::endl << std::endl;
		return os;
	}
	friend std::istream& operator >> (std::istream& is, MeshMatrix& p) {
		is >> p.ap;
		is >> p.an[0];
		is >> p.an[1];
		is >> p.Su;
		is >> p.adg;
		return is;
	}
};

/*typedefs*/
typedef MeshMatrix<Scalar>  ScalarCellMatrix;
typedef MeshMatrix<Vector>  VectorCellMatrix;
typedef MeshMatrix<Tensor>  TensorCellMatrix;
typedef MeshMatrix<STensor> STensorCellMatrix;

/*************************************
 *  Asynchronous communication
 *************************************/
template <class T> 
class ASYNC_COMM {
private:
	T* P;
	Int rcount;
	std::vector<MP::REQUEST> request;
public:
	ASYNC_COMM(T* p) : P(p)
	{
	}
	void send() {
		using namespace Mesh;
		using namespace DG;
		
		//---fill send buffer and send
		MeshField<T,CELL> sendbuf;
		request.assign(2 * Mesh::gInterMesh.size(),0);
		rcount = 0;
		forEach(gInterMesh,i) {
			interBoundary& b = gInterMesh[i];
			IntVector& f = *(b.f);
			Int buf_size = b.f->size() * NPF;
			
			//--fill send buffer
			forEach(f,j) {
				Int faceid = f[j];
				for(Int n = 0; n < NPF;n++) {
					Int k = faceid * NPF + n;
					sendbuf[(b.buffer_index + j) * NPF + n] = P[gFO[k]];	
				}															
			}	

			//--non-blocking send/recive
			MP::isend(&sendbuf[b.buffer_index * NPF],buf_size,
				b.to,MP::FIELD,&request[rcount]);
			rcount++;
			MP::irecieve(&P[gFN[f[0]]],buf_size,
				b.to,MP::FIELD,&request[rcount]);
			rcount++;
		}
	}
	void recv() {
 		MP::waitall(rcount,&request[0]);
	}
};

/*************************************
 * Exchange ghost cell information
 *************************************/
template <class T> 
void exchange_ghost(T* P) {
	ASYNC_COMM<T> comm(P);
	comm.send();
	comm.recv();
}
 
/* *******************************
 * matrix - vector product p * q
 * *******************************/
template <class T> 
MeshField<T,CELL> mul (const MeshMatrix<T>& p,const MeshField<T,CELL>& q, const bool sync = false) {
	using namespace Mesh;
	using namespace DG;
	MeshField<T,CELL> r;
	Int c1,c2;
	ASYNC_COMM<T> comm(&q[0]);
	
	if(sync) comm.send();
	
	r = q * p.ap;
	if(NPMAT) {
		for(Int i = 0;i < gBCS;i++) {
			for(Int j = 0;j < NP;j++) {
				T val(Scalar(0));
				for(Int k = 0;k < NP;k++)
					val += q[i * NP + k] * 
						p.adg[i * NPMAT + j * NP + k];
				r[i * NP + j] -= val;
			}
		}
	}
	
	for(Int f = 0; f < gBFS; f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] -= q[c2] * p.an[1][f];
		r[c2] -= q[c1] * p.an[0][f];
	}
	
	if(sync) comm.recv();
	
	for(Int f = gBFS; f < fI.size(); f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] -= q[c2] * p.an[1][f];
	}
	
	return r;
}
/*matrix transopose - vector product pT * q */
template <class T> 
MeshField<T,CELL> mult (const MeshMatrix<T>& p,const MeshField<T,CELL>& q, const bool sync = false) {
	using namespace Mesh;
	using namespace DG;
	MeshField<T,CELL> r;
	Int c1,c2;
	ASYNC_COMM<T> comm(&q[0]);
	
	if(sync) comm.send();
	
	r = q * p.ap;
	if(NPMAT) {
		for(Int i = 0;i < gBCS;i++) {
			for(Int j = 0;j < NP;j++) {
				T val(Scalar(0));
				for(Int k = 0;k < NP;k++)
					val += q[i * NP + k] * 
						p.adg[i * NPMAT + k * NP + j];
				r[i * NP + j] -= val;
			}
		}
	}
	
	for(Int f = 0; f < gBFS; f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c2] -= q[c1] * p.an[1][f];
		r[c1] -= q[c2] * p.an[0][f];
	}
	
	if(sync) comm.recv();
	
	for(Int f = gBFS; f < fI.size(); f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] -= q[c2] * p.an[0][f];
	}
	
	return r;
}
/* calculate RHS sum */
template <class T> 
MeshField<T,CELL> getRHS(const MeshMatrix<T>& p, const bool sync = false) {
	using namespace Mesh;
	using namespace DG;
	MeshField<T,CELL> r;
	Int c1,c2;
	ASYNC_COMM<T> comm(&(*p.cF)[0]);
	
	if(sync) comm.send();
	
	r = p.Su;
	if(NPMAT) {
		for(Int i = 0;i < gBCS;i++) {
			for(Int j = 0;j < NP;j++) {
				T val(Scalar(0));
				for(Int k = 0;k < NP;k++)
					val += (*p.cF)[i * NP + k] * 
						p.adg[i * NPMAT + j * NP + k];
				r[i * NP + j] += val;
			}
		}
	}
	for(Int f = 0; f < gBFS; f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] += (*p.cF)[c2] * p.an[1][f];
		r[c2] += (*p.cF)[c1] * p.an[0][f];
	}
	
	if(sync) comm.recv();
	
	for(Int f = gBFS; f < fI.size(); f++) {
		c1 = gFO[f];
		c2 = gFN[f];
		r[c1] += (*p.cF)[c2] * p.an[1][f];
	}
	
	return r;
}

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
	 forEach(AllBConditions,i) {
		 bbc = AllBConditions[i];
		 if(bbc->fIndex == cF.fIndex) {
			 if(bbc->cIndex == NEUMANN ||
				 bbc->cIndex == SYMMETRY ||
				 bbc->cIndex == ROBIN)
				 ;
			 else continue;

			 bc = static_cast<BCondition<T>*> (bbc);
			 Int sz = bc->bdry->size();
			 if(sz == 0) continue;

			 for(Int j = 0;j < sz;j++) {
				 Int faceid = (*bc->bdry)[j];
				 for(Int n = 0; n < DG::NPF;n++) {
				 	 Int k = faceid * DG::NPF + n;
				 	 
					 Int c1 = gFO[k];
					 Int c2 = gFN[k];
					 if(bc->cIndex == NEUMANN) {
						 Vector dv = cC[c2] - cC[c1];
						 M.ap[c1] -= M.an[1][k];
						 M.Su[c1] += M.an[1][k] * (bc->value * mag(dv));
						 M.an[1][k] = 0;
					 } else if(bc->cIndex == ROBIN) {
						 Vector dv = cC[c2] - cC[c1];
						 M.ap[c1] -= (1 - bc->shape) * M.an[1][k];
						 M.Su[c1] += M.an[1][k] * (bc->shape * bc->value + 
							 (1 - bc->shape) * bc->tvalue * mag(dv));
						 M.an[1][k] = 0;
					 } else if(bc->cIndex == SYMMETRY) {
						 M.ap[c1] -= M.an[1][k];
						 M.Su[c1] += M.an[1][k] * (sym(cF[c1],fN[k]) - cF[c1]);
						 M.an[1][k] = 0;
					 }
				 }
			 }
		 }
	 }
 }
 

 /* ***************************************
 * Explicit boundary conditions
 * **************************************/
template<class T,ENTITY E>
void applyExplicitBCs(const MeshField<T,E>& cF,
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

	/*boundary conditions*/
	forEach(AllBConditions,i) {
		bbc = AllBConditions[i];
		if(bbc->fIndex == cF.fIndex) {
			if(bbc->cIndex == GHOST) 
				continue;

			bc = static_cast<BCondition<T>*> (bbc);
			Int sz = bc->bdry->size();
			if(sz == 0) continue;
			if(!bc->fixed.size())
				bc->fixed.resize(sz * DG::NPF);

			if(update_fixed) {
				if(bc->cIndex == DIRICHLET || 
					bc->cIndex == POWER || 
					bc->cIndex == LOG || 
					bc->cIndex == PARABOLIC ||
					bc->cIndex == INVERSE
					) {
						if(bc->zMax > 0) {
							zmin = bc->zMin;
							zmax = bc->zMax;
							zR = zmax - zmin;
						} else {
							zmin = Scalar(10e30);
							zmax = -Scalar(10e30);
							C = Vector(0);
							for(Int j = 0;j < sz;j++) {
								Facet& f = gFacets[j];
								Vector fc(Scalar(0));
								forEach(f,k) {
									fc += vC[f[k]];
									z = (vC[f[k]] & bc->dir);
									if(z < zmin) 
										zmin = z;
									if(z > zmax) 
										zmax = z;
								}
								C += (fc / f.size());
							}
							C /= Scalar(sz);
							zR = zmax - zmin;

							if(bc->cIndex == PARABOLIC) {
								Int vi = gFacets[(*bc->bdry)[0]][0];
								zR = magSq(vC[vi] - C);
								for(Int j = 1;j < sz;j++) {
									vi = gFacets[(*bc->bdry)[j]][0];
									Scalar r = magSq(vC[vi] - C);
									if(r < zR) zR = r;
								}
							}
						}
				}
			}
			for(Int j = 0;j < sz;j++) {
				Int faceid = (*bc->bdry)[j];
				for(Int n = 0; n < DG::NPF;n++) {
				 	Int k = faceid * DG::NPF + n;

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
					} else if(bc->cIndex == CYCLIC ||
							  bc->cIndex == RECYCLE) {
						Int fi;
						if(j < sz / 2) 
							fi = (*bc->bdry)[j + sz/2];
						else
							fi = (*bc->bdry)[j - sz/2];
						Int k1 = fi * DG::NPF + n;
						Int c11 = gFO[k1];
				 		Int c22 = gFN[k1];
						if(bc->cIndex == CYCLIC) {
							cF[c2] = cF[c11];
							cF[c22] = cF[c1];
						} else {
							if(c2 >= gBCSfield)
								cF[c2] = cF[c11];
							else
								cF[c22] = cF[c1];
						}
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
							if(!bc->first && !equal(mag(bc->tvalue),0)) { 
								T meanTI = v * (bc->tvalue * pow (z / zR,-bc->tshape));
								Scalar rFactor = 4 * ((rand() / Scalar(RAND_MAX)) - 0.5);
								v += ((cF[c2] - v) * 0.9 + (meanTI * rFactor) * 0.1);
							}
							bc->fixed[j * DG::NPF + n] = cF[c2] = v;
						} else {
							cF[c2] = bc->fixed[j * DG::NPF + n];
						}
					}
				}
			}
			bc->first = false;
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
const MeshField<T,E>& fillBCs(const MeshField<T,E>& cF,
			 bool update_ghost = false) {
	/*neumann update*/
	using namespace Mesh;
	forEachS(gCells,i,gBCS) {
		Int faceid = gCells[i][0];
		for(Int n = 0; n < DG::NPF;n++) {
			Int k = faceid * DG::NPF + n;
			cF[gFN[k]] = cF[gFO[k]];
		}
	}
	/*ghost cells*/
	if(update_ghost && gInterMesh.size()) {
		exchange_ghost(&cF[0]);
	}
	return cF;
}

/*******************************************************************************
 * Non-integrated field operations. The field operations to be defined later
 * such as div,lap,grad,ddt,ddt2,src etc.. are values integrated over a volume. 
 * The following macro versions give the corresponding non-integrated cell
 * center values by dividing with the cell volumes
 *******************************************************************************/
#define divi(x)  fillBCs((div(x)   / Mesh::cV),true)
#define lapi(x)  fillBCs((lap(x)   / Mesh::cV),true)
#define gradi(x) fillBCs((grad(x)  / Mesh::cV),true)
#define ddti(x)  fillBCs((ddt(x)   / Mesh::cV),true)
#define ddt2i(x) fillBCs((ddt2(x)  / Mesh::cV),true)
#define srci(x)  fillBCs((src(x)   / Mesh::cV),true)

/* **************************
 * Integrate field operation
 * **************************/
template<class type>
MeshField<type,CELL> sum(const MeshField<type,FACET>& fF) {
	using namespace Mesh;
	MeshField<type,CELL> cF;
	cF = type(0);
	forEach(fF,i) {
		cF[gFO[i]] += fF[i];
		cF[gFN[i]] -= fF[i];
	}
	return cF;
}
/* ********************************
 * Interpolate field operations
 * *******************************/
/*central difference*/
template<class type>
MeshField<type,FACET> cds(const MeshField<type,CELL>& cF) {
	using namespace Mesh;
	MeshField<type,FACET> fF;
	forEach(fF,i) {
		fF[i] =  (cF[gFO[i]] * (fI[i])) + (cF[gFN[i]] * (1 - fI[i]));
	}
	return fF;
}

/*upwind*/
template<class type>
MeshField<type,FACET> uds(const MeshField<type,CELL>& cF,const ScalarFacetField& flux) {
	using namespace Mesh;
	MeshField<type,FACET> fF;
	forEach(fF,i) {
		if(flux[i] >= 0) fF[i] = cF[gFO[i]];
		else fF[i] = cF[gFN[i]];
	}
	return fF;
}

/*facet data to vertex data */
template<class type>
MeshField<type,VERTEX> cds(const MeshField<type,FACET>& fF) {
	using namespace Mesh;
	MeshField<type,VERTEX> vF;
	std::vector<Scalar> cnt;
	
	cnt.assign(vF.size(),Scalar(0));
	vF = type(0);
	
	forEach(fF,i) {
		Facet& f = gFacets[i];
		if(gFN[i] < gBCS) {
			forEach(f,j) {
				Scalar dist = 1.f / magSq(gVertices[f[j]] - fC[i]);
				vF[f[j]] += (fF[i] * dist);
				cnt[f[j]] += dist;
			}
		} else {
			forEach(f,j) {
				vF[f[j]] += Scalar(10e30) * fF[i];
				cnt[f[j]] += Scalar(10e30);
			}
		}
	}
	
	forEach(vF,i) {
		vF[i] /= cnt[i];
		if(mag(vF[i]) < Constants::MachineEpsilon) 
			vF[i] = type(0);
	}
	
	return vF;
}
/* *******************************
 * Linearized source term
 * *******************************/
template<class type>
MeshMatrix<type> src(MeshField<type,CELL>& cF,const MeshField<type,CELL>& Su,const ScalarCellField& Sp) {
	MeshMatrix<type> m;
	m.cF = &cF;
	m.flags |= (m.SYMMETRIC | m.DIAGONAL);
	m.ap = -(Sp * Mesh::cV);
	m.Su = (Su * Mesh::cV);
	//others
	m.adg = Scalar(0);
	m.an[0] = Scalar(0);
	m.an[1] = Scalar(0);
	return m;
}
/*Explicit*/
template<class type>
MeshField<type,CELL> src(const MeshField<type,CELL>& Su) {
	return (Su * Mesh::cV);
}
 
 /***********************************************
 * Gradient field operation.
 ***********************************************/

/*Explicit*/
inline VectorCellField grad(const ScalarFacetField& p) {
	return sum(mul(Mesh::fN,p));
}
inline VectorCellField grad(const ScalarCellField& p) {
	return grad(cds(p));
}
inline TensorCellField grad(const VectorFacetField& p) {
	return sum(mul(Mesh::fN,p));
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
	m.adg = Scalar(0);
	forEach(mu,i) {
		c1 = gFO[i];
		c2 = gFN[i];
		dv = cC[c2] - cC[c1];
		/*diffusivity coefficient*/
		if(nonortho_scheme == NONE) {
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
	if(nonortho_scheme != NONE) {
		MeshField<type,FACET> r = dot(cds(gradi(cF)),K);
		type res;
		forEach(mu,i) {
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
MeshMatrix<type> div(MeshField<type,CELL>& cF,const VectorCellField& flux_cell,
					const ScalarFacetField& flux,const ScalarFacetField& mu) {
	using namespace Controls;
	using namespace Mesh;
	using namespace DG;
	MeshMatrix<type> m;
	Scalar F,G;
	m.cF = &cF;
	m.flags = 0;
	m.Su = type(0);
	m.ap = Scalar(0);

	/*compute differentiation matrix - volume integral*/
	if(NPMAT) {
		for(Int ci = 0; ci < gBCS;ci++) {
			forEachLgl(ii,jj,kk) {
				Int ind1 = INDEX3(ii,jj,kk);
				forEachLgl(i,j,k) {
					Int ind2 = INDEX3(i,j,k);
				
					Int index = INDEX4(ci,ii,jj,kk);
					Tensor& Jin = Jinv[index];
					Scalar dpsi_j_0 = dpsi[0][ii][i] *   psi[1][jj][j] *   psi[2][kk][k];
					Scalar dpsi_j_1 =  psi[0][ii][i] *  dpsi[1][jj][j] *   psi[2][kk][k];
					Scalar dpsi_j_2 =  psi[0][ii][i] *   psi[1][jj][j] *  dpsi[2][kk][k];
					Vector dpsi_j = Vector(dpsi_j_0,dpsi_j_1,dpsi_j_2);
					Scalar dpsiU = flux_cell[index] & dot(dpsi_j,Jin);
					Scalar val = dpsiU * cV[index];
					index = ci * NPMAT + ind2 * NP + ind1;
					if(ind1 == ind2) {
						m.ap[ci * NP + ind1] = val;
						m.adg[index] = 0;
					} else {
						m.adg[index] = -val;
					}
				}
			}
		}
	}
	
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
			forEach(gFacets,faceid) {
				for(Int n = 0; n < NPF;n++) {
					Int k = faceid * NPF + n;
					/*calc D - uncorrected */
					dv = cC[gFN[k]] - cC[gFO[k]];
					D = (mag(fN[k]) / mag(dv)) * mu[k];
					/*compare F and D */
					F = flux[k];
					if(F < 0) {
						if(-F * fI[k] > D) gamma[k] = 0;
						else gamma[k] = 1;
					} else {
						if(F * (1 - fI[k]) > D) gamma[k] = 0;
						else gamma[k] = 1;
					}
				}
			}
		}
		forEach(flux,i) {
			F = flux[i];
			G = gamma[i];
			m.an[0][i] = ((G) * (-F * (  fI[i]  )) + (1 - G) * (-max( F,0)));
			m.an[1][i] = ((G) * ( F * (1 - fI[i])) + (1 - G) * (-max(-F,0)));
			m.ap[gFO[i]] += m.an[0][i];
			m.ap[gFN[i]] += m.an[1][i];
		}
	/*deferred correction*/
	} else {
		forEach(flux,i) {
			F = flux[i];
			m.an[0][i] = -max( F,0);
			m.an[1][i] = -max(-F,0);
			m.ap[gFO[i]] += m.an[0][i];
			m.ap[gFN[i]] += m.an[1][i];
		}

		MeshField<type,FACET> corr;
		if(convection_scheme == CDSS) {
			corr = cds(cF) - uds(cF,flux);
		} else if(convection_scheme == LUD) {
			VectorFacetField R = fC - uds(cC,flux);
			corr = dot(uds(gradi(cF),flux),R);
		} else if(convection_scheme == MUSCL) {
			VectorFacetField R = fC - uds(cC,flux);
			corr  = (  blend_factor  ) * (cds(cF) - uds(cF,flux));
			corr += (1 - blend_factor) * (dot(uds(gradi(cF),flux),R));
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
				forEach(phiDC,i) {
					if(flux[i] >= 0) G = fI[i];
					else G = 1 - fI[i];
					uFI[i] = G;
				}
				/*Bruner's or Darwish way of calculating r*/
				if(TVDbruner) {
					VectorFacetField R = fC - uds(cC,flux);
					phiCU = 2 * (dot(uds(gradi(cF),flux),R));
				} else {
					VectorFacetField R = uds(cC,nflux) - uds(cC,flux);
					phiCU = 2 * (dot(uds(gradi(cF),flux),R)) - phiDC;
				}
				/*end*/
			}
			r = (phiCU / phiDC) * (uFI / (1 - uFI));
			forEach(phiDC,i) {
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
	return div(cF,rhoU,div(rhoU),mu);
}

/* *******************************
 * Temporal derivative
 * *******************************/
template<class type>
MeshMatrix<type> ddt(MeshField<type,CELL>& cF,const ScalarCellField& rho) {
	MeshMatrix<type> m;
	m.cF = &cF;
	m.flags |= (m.SYMMETRIC | m.DIAGONAL);
	//ddt schemes
	if(Controls::time_scheme == Controls::EULER || !(cF.access & STOREPREV)) {
		if(Controls::time_scheme != Controls::EULER) cF.initStore();
		m.ap = (Mesh::cV * rho) / -Controls::dt;
		m.Su = cF * m.ap;
	} else if(Controls::time_scheme == Controls::SECOND_ORDER) {
		m.ap = (1.5 * Mesh::cV * rho) / -Controls::dt;
		m.Su = ((4.0 * cF - cF.tstore[1]) / 3.0) * m.ap;
	}
	//others
	m.adg = Scalar(0);
	m.an[0] = Scalar(0);
	m.an[1] = Scalar(0);
	return m;
}
template<class type>
MeshMatrix<type> ddt2(MeshField<type,CELL>& cF,const ScalarCellField& rho) {
	MeshMatrix<type> m;
	m.cF = &cF;
	m.flags |= (m.SYMMETRIC | m.DIAGONAL);
	if(!(cF.access & STOREPREV)) cF.initStore();
	//diagonal mass matrix
	m.ap = (Mesh::cV * rho) / -(Controls::dt * Controls::dt);
	m.Su = (2.0 * cF - cF.tstore[1]) * m.ap;
	//others
	m.adg = Scalar(0);
	m.an[0] = Scalar(0);
	m.an[1] = Scalar(0);
	return m;
}	
template<class type>
void addTemporal(MeshMatrix<type>& M,const ScalarCellField& rho,Scalar cF_UR) {
	using namespace Controls;
	if(state == STEADY)
		M.Relax(cF_UR);
	else {
		if(!equal(implicit_factor,1)) {
			MeshField<type,CELL> k1 = M.Su - mul(M, *M.cF);
			if(runge_kutta == 1) {
				M = M * (implicit_factor) + k1  * (1 - implicit_factor);
			} else {
				ScalarCellField mdt = Controls::dt / (rho * Mesh::cV);
				if (runge_kutta == 2) {
					MeshField<type, CELL> k2 = M.Su - mul(M, *M.cF + k1 * mdt);
					M = M * (implicit_factor) + ((k2 + k1) / 2) * (1 - implicit_factor);
				} else if (runge_kutta == 3) {
					MeshField<type, CELL> k2 = M.Su - mul(M, *M.cF + k1 * mdt / 2);
					MeshField<type, CELL> k3 = M.Su - mul(M, *M.cF + (2 * k2 - k1) * mdt);
					M = M * (implicit_factor) + ((k3 + 4 * k2 + k1) / 6) * (1 - implicit_factor);
				} else if (runge_kutta == 4) {
					MeshField<type, CELL> k2 = M.Su - mul(M, *M.cF + k1 * mdt / 2);
					MeshField<type, CELL> k3 = M.Su - mul(M, *M.cF + k2 * mdt / 2);
					MeshField<type, CELL> k4 = M.Su - mul(M, *M.cF + k3 * mdt);
					M = M * (implicit_factor) + ((k4 + 2 * k3 + 2 * k2 + k1) / 6) * (1 - implicit_factor);
				}
			}
			if(equal(implicit_factor,0))
				M.flags = (M.SYMMETRIC | M.DIAGONAL);
		}
		M += ddt(*M.cF,rho);
	}
}
/* ************************************************
 * Form transport equation
 *************************************************/
template<class type>
MeshMatrix<type> diffusion(MeshField<type,CELL>& cF,
		const ScalarFacetField& mu,const ScalarCellField& rho,Scalar cF_UR) {
	MeshMatrix<type> M = -lap(cF,mu);
	addTemporal(M,rho,cF_UR);
	return M;
}
template<class type>
MeshMatrix<type> diffusion(MeshField<type,CELL>& cF,
		const ScalarFacetField& mu,const ScalarCellField& rho,Scalar cF_UR, 
		const MeshField<type,CELL>& Su,const ScalarCellField& Sp) {
	MeshMatrix<type> M = -lap(cF,mu) - src(cF,Su,Sp);
	addTemporal(M,rho,cF_UR);
	return M;
}
template<class type>
MeshMatrix<type> convection(MeshField<type,CELL>& cF,const VectorCellField& Fc,
        const ScalarFacetField& F,const ScalarCellField& rho,Scalar cF_UR) {
	MeshMatrix<type> M = div(cF,Fc,F,Scalar(0));
	addTemporal(M,rho,cF_UR);
	return M;
}
template<class type>
MeshMatrix<type> convection(MeshField<type,CELL>& cF,const VectorCellField& Fc,
        const ScalarFacetField& F, const ScalarCellField& rho,Scalar cF_UR, 
		const MeshField<type,CELL>& Su,const ScalarCellField& Sp) {
	MeshMatrix<type> M = div(cF,Fc,F,Scalar(0)) - src(cF,Su,Sp);
	addTemporal(M,rho,cF_UR);
	return M;
}
template<class type>
MeshMatrix<type> transport(MeshField<type,CELL>& cF,const VectorCellField& Fc,const ScalarFacetField& F,
		const ScalarFacetField& mu,const ScalarCellField& rho,Scalar cF_UR) {
	MeshMatrix<type> M = div(cF,Fc,F,mu) - lap(cF,mu);
	addTemporal(M,rho,cF_UR);
	return M;
}
template<class type>
MeshMatrix<type> transport(MeshField<type,CELL>& cF,const VectorCellField& Fc,const ScalarFacetField& F,
		const ScalarFacetField& mu,const ScalarCellField& rho,Scalar cF_UR, 
		const MeshField<type,CELL>& Su,const ScalarCellField& Sp) {
	MeshMatrix<type> M = div(cF,Fc,F,mu) - lap(cF,mu) - src(cF,Su,Sp);
	addTemporal(M,rho,cF_UR);
	return M;
}
/* ********************
 *        End
 * ********************/
#endif
