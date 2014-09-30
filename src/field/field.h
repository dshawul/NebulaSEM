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
		NOP,DIAG,SORP,DILU
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
	extern CommMethod ghost_exchange;
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
    extern Int NP,NPM1;
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
	MeshField(const char* str = "", ACCESS a = NO) : 
				P(0),allocated(0),access(a),fName(str) {
		construct(str,a);
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
	void allocate() {
		if(mem_.empty()) {
			switch(entity) {
				case CELL:   SIZE = Mesh::gCells.size() * DG::NP;	 break;
				case FACET:  SIZE = Mesh::gFacets.size();            break;
				case VERTEX: SIZE = Mesh::gVertices.size();          break;
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
				case CELL:   SIZE = Mesh::gCells.size() * DG::NP;	 break;
				case FACET:  SIZE = Mesh::gFacets.size();            break;
				case VERTEX: SIZE = Mesh::gVertices.size();          break;
		}
		P = &q[0];
		allocated = 0;
	}
	void construct(const char* str = "", ACCESS a = NO) {
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
					<< Mesh::gBCellsStart << " float" << std::endl;
				for(Int i = 0;i < Mesh::gBCellsStart;i++)
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
#define forEachField(X)	 { \
	ScalarCellField::X;    \
    VectorCellField::X;    \
	STensorCellField::X;   \
	TensorCellField::X;    \
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

	void   initGeomMeshFields(bool = true);
	void   write_fields(Int);
	void   read_fields(Int);
	void   calc_walldist(Int,Int = 1);
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
			ScalarCellField R = mag((Mesh::cC - cB) / radius);
			forEach(R,i) R[i] = min(1.0,R[i]);
			MeshField<T,E> val = value;
			val += MeshField<T,E>(perterb / 2) * 
					(MeshField<Scalar,E>(1.0) + cos(R * Constants::PI));
			*this = val;
		} else if(str == "hydrostatic") {
			T p0;
			Scalar scale, expon;
			is >> p0 >> scale >> expon;
			Vector gu = unit(Controls::gravity);
			ScalarCellField gz = -dot(Mesh::cC,VectorCellField(gu));
			*this = MeshField<T,E>(p0) * pow(MeshField<Scalar,E>(1.0) - scale * gz, expon);
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
	updateExplicitBCs(*this,true,true);
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
	Int size = (SIZE == gCells.size()) ? gBCellsStart * DG::NP : SIZE;
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

/*********************************************************************************
 *                      matrix class defined on mesh                             
 *********************************************************************************/
template <class type> 
struct MeshMatrix {
	MeshField<type,CELL>*   cF;
	MeshField<type,CELL>    Su;
	MeshField<Scalar,CELL>  ap;
	MeshField<Scalar,CELL>* adg;
	MeshField<Scalar,FACET>   an[2];
	Int flags;
	enum FLAG {
		SYMMETRIC = 1
	};
	/*c'tors*/
	MeshMatrix() {
		cF = 0;
		flags = 0;
		adg = new MeshField<Scalar,CELL>[DG::NPM1];
	}
	MeshMatrix(const MeshMatrix& p) {
		cF = p.cF;
		flags = p.flags;
		ap = p.ap;
		an[0] = p.an[0];
		an[1] = p.an[1];
		Su = p.Su;
		adg = new MeshField<Scalar,CELL>[DG::NPM1];
		for(Int i = 0; i < DG::NPM1;i++)
			adg[i] = p.adg[i];
	}
	MeshMatrix(const MeshField<type,CELL>& p) {
		cF = 0;
		flags = SYMMETRIC;
		ap = Scalar(0);
		an[0] = Scalar(0);
		an[1] = Scalar(0);
		Su = p;
		adg = new MeshField<Scalar,CELL>[DG::NPM1];
		for(Int i = 0; i < DG::NPM1;i++)
			adg[i] = Scalar(0);
	}
	~MeshMatrix() {
		if(DG::NPM1)
			delete adg;
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
		for(Int i = 0; i < DG::NPM1;i++)
			r.adg[i] = -adg[i];
		return r;
	}
	MeshMatrix& operator = (const MeshMatrix& q) {
		cF = q.cF;
		flags = q.flags;
		ap = q.ap;
		an[0] = q.an[0];
		an[1] = q.an[1];
		Su = q.Su;
		for(Int i = 0; i < DG::NPM1;i++)
			adg[i] = q.adg[i];
		return *this;
	}
	MeshMatrix& operator += (const MeshMatrix& q) {
		flags &= q.flags;
		ap += q.ap;
		an[0] += q.an[0];
		an[1] += q.an[1];
		Su += q.Su;
		for(Int i = 0; i < DG::NPM1;i++)
			adg[i] += q.adg[i];
		return *this;
	}
	MeshMatrix& operator -= (const MeshMatrix& q) {
		flags &= q.flags;
		ap -= q.ap;
		an[0] -= q.an[0];
		an[1] -= q.an[1];
		Su -= q.Su;
		for(Int i = 0; i < DG::NPM1;i++)
			adg[i] -= q.adg[i];
		return *this;
	}
	MeshMatrix& operator *= (const Scalar& q) {
		ap *= q;
		an[0] *= q;
		an[1] *= q;
		Su *= q;
		for(Int i = 0; i < DG::NPM1;i++)
			adg[i] *= q;
		return *this;
	}
	MeshMatrix& operator /= (const Scalar& q) {
		ap /= q;
		an[0] /= q;
		an[1] /= q;
		Su /= q;
		for(Int i = 0; i < DG::NPM1;i++)
			adg[i] /= q;
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
		for(Int i = 0; i < DG::NPM1;i++)
			os << p.adg[i] << " ";
		os << std::endl << std::endl;
		return os;
	}
	friend std::istream& operator >> (std::istream& is, MeshMatrix& p) {
		is >> p.ap;
		is >> p.an[0];
		is >> p.an[1];
		is >> p.Su;
		for(Int i = 0; i < DG::NPM1;i++)
			is >> p.adg[i];
		return is;
	}
};

/*typedefs*/
typedef MeshMatrix<Scalar>  ScalarCellMatrix;
typedef MeshMatrix<Vector>  VectorCellMatrix;
typedef MeshMatrix<Tensor>  TensorCellMatrix;
typedef MeshMatrix<STensor> STensorCellMatrix;

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
		register Int i,f;

		/*resize*/
		cF.resize(N);
		Su.resize(N);
		rows.reserve(N + 1);
		cols.reserve(NN);
		an.reserve(NN);
		anT.reserve(NN);

        /*source term*/
		for(i = 0;i < N;i++) {
			Su[i] = A.Su[i];
			cF[i] = (*A.cF)[i];
		}

		/*fill matrix in CSR format.Diagonal element 
		  is always at the start of a row */
		Int cn = 0;
		for(i = 0;i < N;i++) {
			Cell& c = gCells[i];

			rows.push_back(cn);

			an.push_back(A.ap[i]);
			anT.push_back(A.ap[i]);
			cols.push_back(i);
			cn++;

			forEach(c,j) {
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
		/*push extra row*/
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
 *   FVM
 * ********************/
 #include "fvm.h"
 #include "dg.h"
 
/* ********************
 *        End
 * ********************/
#endif
