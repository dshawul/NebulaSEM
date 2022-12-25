#ifndef __FIELD_H
#define __FIELD_H

#include "mesh.h"
#include "mp.h"

/** Basic building blocks (entities) over which fields are defined */
enum ENTITY {
    CELL,    /**< A cell (volume), each LGL node gets one cell in DG. */
    FACET,   /**< A polygonal face. */
    VERTEX,  /**< A node. */
    CELLMAT  /**< An NXN matrix for each element. */
};

/** AMR parameters */
struct RefineParams {
    Vector dir;
    std::string field;
    Scalar field_max;
    Scalar field_min;
    Int limit;
    RefineParams() {
        dir = Scalar(0);
        field = "U";
        field_max = 0.6;
        field_min = 0.2;
        limit = 100000;
    }
};

/** Decomposition parameters */
struct DecomposeParams {
    Int type;
    IntVector n;
    ScalarVector axis;
    DecomposeParams() {
        type = 2;
        axis.assign(4,0);
        axis[0] = 1;
        n.assign(3,1);
    }
};

namespace Mesh {
    void enroll(Util::ParamList& params);
}

namespace Controls {
    extern RefineParams refine_params;
    extern DecomposeParams decompose_params;
    void enrollRefine(Util::ParamList& params);
    void enrollDecompose(Util::ParamList& params);
}

/*******************************************************************************
 *                              Boundary conditions
 *******************************************************************************/
/**
  \verbatim
  Model for flow close to the wall (Law of the wall).
  1 -> Viscous layer
  2 -> Buffer layer
  3 -> Log-law layer
  The wall function model is modified for rough surfaces 
  using Cebecci and Bradshaw formulae.
  \endverbatim
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
/** 
  Boundary condition types
 */
namespace Mesh {
    const Int DIRICHLET    = Util::hash_function("DIRICHLET");
    const Int NEUMANN      = Util::hash_function("NEUMANN");
    const Int ROBIN        = Util::hash_function("ROBIN");
    const Int SYMMETRY     = Util::hash_function("SYMMETRY");
    const Int CYCLIC       = Util::hash_function("CYCLIC");
    const Int RECYCLE      = Util::hash_function("RECYCLE");
    const Int GHOST        = Util::hash_function("GHOST");
    const Int POWER        = Util::hash_function("POWER");
    const Int LOG          = Util::hash_function("LOG");
    const Int PARABOLIC    = Util::hash_function("PARABOLIC");
    const Int INVERSE      = Util::hash_function("INVERSE");
    const Int ROUGHWALL    = Util::hash_function("ROUGHWALL");
    const Int CALC_NEUMANN = Util::hash_function("CALC_NEUMANN");
}
/** Basic boundary condition */
struct BasicBCondition {
    IntVector* bdry;
    Int     fIndex;
    Int     cIndex;
    std::string cname;
    std::string bname;
    std::string fname;
    LawOfWall low;
};
/** Template boundary condition's class for different tensors */
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
        fixed.clear();
    }
    void init_indices() {
        bdry = &Mesh::gBoundaries[bname];
        first = true;
        read = false;
        fIndex = Util::hash_function(fname);
        cIndex = Util::hash_function(cname);
    }
};
/** Write boundary conditions */
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
/** Read boundary conditions */
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

namespace Mesh {
    /** List of probing points */
    extern  Vertices         probePoints;
    /** List of all boundary conditions*/
    extern  std::vector<BasicBCondition*> AllBConditions;
    /** Clear list of  BCs */
    inline void clearBC() {
        forEach(AllBConditions,i)
            delete AllBConditions[i];
        AllBConditions.clear();
    }
}

/*******************************************************************************
 *                              Control parameters
 *******************************************************************************/
namespace Controls {
    /** Convection differencing schemes */
    enum Scheme{
        CDS,    /**< Central difference */
        UDS,    /**< Upwind difference */
        HYBRID, /**< Hybrid scheme that switches b/n CDS/UDS */
        BLENDED,/**< Blended CDS/UDS difference */
        LUD,    /**< Linear upwind */
        CDSS,   /**< Stabilized central difference */
        MUSCL,  /**< Monotonic upstream centered */
        QUICK,  /**< Upstream weighted quadratic  */
        VANLEER,/**< Van Leer's scheme */
        VANALBADA,  /**< Van Albada's scheme */
        MINMOD,     /**< Minmod scheme */
        SUPERBEE,   /**< Superbee scheme */
        SWEBY,      /**< Sweby scheme */
        QUICKL,     /**< Limited QUICK scheme */
        UMIST,      /**< Lien & Leschziner, 1994 */
        DDS,        /**< DDS scheme */
        FROMM       /**< FROMM scheme */
    };
    /** Non-orthogonality correction schemes (Jasak) */
    enum NonOrthoScheme {
        NONE,           /**< No correction */
        MINIMUM,        /**< Minimum correction */
        ORTHOGONAL,     /**< Orthogonal correction */
        OVER_RELAXED    /**< Over-relaxed correction (default) */
    };
    /** Backward differencing formula */
    enum TimeScheme {
        BDF1,   /**< First order */
        BDF2,   /**< Second order */
        BDF3,   /**< Third order */
        BDF4,   /**< Fourth order */
        BDF5,   /**< Fifth order */
        BDF6    /**< Sixth order */
    };
    /** Iterative solvers */
    enum Solvers {
        JACOBI, /**< Jacobi */
        SOR,    /**< Successive over-relaxation */
        PCG     /**< Pre-conditioned conjugate gradient */
    };
    /** Preconditioners for conjugate gradient */
    enum Preconditioners {
        NOPR,   /**< No preconditioner */
        DIAG,   /**< Diagonal (Jacobi) preconditioner */
        SSOR,   /**< Symmetric SOR preconditioner */
        DILU    /**< Diagonal incomplete LU factorization */
    };
    /** Communication methods */
    enum CommMethod {
        BLOCKED,        /**< Blocked send/recv */
        ASYNCHRONOUS    /**< Asycnhronous communication */
    };
    /** Steady/transient state */
    enum State {
        STEADY,     /**< Steady state solution sought */
        TRANSIENT   /**< Transient solution sought */
    };

    extern Scheme convection_scheme;
    extern Int TVDbruner;
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
    extern Int amr_step;
    extern Int n_deferred;
    extern Int save_average;
    extern Int print_time;

    extern Vector gravity;
}

/** Read/Write access for field */
enum ACCESS {
    NO = 0,         /**< No access allowed */
    READ = 1,       /**< Read only access */
    WRITE = 2,      /**< Write only access */
    READWRITE = 3,  /**< Read write access */
    STOREPREV = 4   /**< Store previous iteration value */
};

namespace DG {
    extern Int Nop[3];
    extern Int NP, NPI, NPMAT, NPF;
};

namespace Mesh {
    extern IntVector  probeCells;
    extern Int  gBCSfield; 
    extern Int  gBCSIfield;
};

/* *****************************************************************************
 *                    Field variables defined on mesh                          
 * *****************************************************************************/

/** Base field class */
class BaseField {   
    public:
        std::string  fName;
    public:
        virtual void deallocate(bool) = 0;
        virtual void refineField(Int,IntVector&,IntVector&) = 0;
        virtual void writeInternal(std::ostream&,IntVector*) = 0;
        virtual void readInternal(std::istream&,IntVector*) = 0;
        virtual void writeBoundary(std::ostream&) = 0;
        virtual void readBoundary(std::istream&) = 0;
        virtual void read(Int step) = 0;
        virtual void write(Int, IntVector* = 0) = 0;
        virtual void norm(BaseField*) = 0;
        virtual ~BaseField() {};

        static std::list<BaseField*> allFields;
        static std::vector<std::string> fieldNames;
        static void destroyFields() {
            std::list<BaseField*> save;
            copyColl(allFields,save);
            forEachIt(save, it)
                (*it)->deallocate(false);
        }
        static BaseField* findField(const std::string& name) {
            forEachIt(allFields, it) {
                if(!Util::compare((*it)->fName,name)) 
                    return (*it);
            }
            return 0;
        }
};

/**
  Template field class for field of type (scalar,vector,tensor)
  defined on entity (vertex,face or cell)
 */
template <class type,ENTITY entity> 
class MeshField : public BaseField, public DVExpr<type,type*>  
{
    private:
        using DVExpr<type,type*>::P;
        static Int   SIZE;
        int          allocated;
    public:
        ACCESS       access;
        Int          fIndex;

        /*common*/
        static const Int TYPE_SIZE = sizeof(type) / sizeof(Scalar);
        static std::list<MeshField*> fields_;
        static std::list<type*> mem_pool;
        static Int n_alloc, n_alloc_max;

        /*constructors*/
        MeshField(const char* str = "", ACCESS a = NO,bool recycle = true) : 
            allocated(0),access(a) {
                P = 0;
                fName = str;
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
            if(!recycle || mem_pool.empty()) {
                switch(entity) {
                    case CELL:   SIZE = Mesh::gCells.size() * DG::NP;    break;
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

                n_alloc++;
                if(n_alloc > n_alloc_max)
                    n_alloc_max = n_alloc;
            } else {
                P = mem_pool.front();
                mem_pool.pop_front();
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
                    mem_pool.push_front(P);
                } else {
                    delete[] P;
                    P = 0;
                    n_alloc--;
                }
                if(fIndex) {
                    fields_.remove(this);
                    allFields.remove(this);
                }

            }
        }
        void construct(const char* str = "", ACCESS a = NO, bool recycle = true) {
            access = a;
            fName = str;
            if(Mesh::gCells.size())
                allocate(recycle);
            fIndex = Util::hash_function(str);
            if(fIndex) {
                fields_.push_back(this);
                allFields.push_back(this);
            }
        }
        /*d'tor re-cycles memory */
        ~MeshField() {
            if(!MP::Terminated)
                deallocate();
        }
        /*accessors*/
        Int size() const {
            return SIZE;
        }
        type& operator [] (Int i) const {
            return P[i];
        }

        /*Assignment from Meshfield and Scalar*/
#define Op($)                                                           \
        MeshField& operator $(const MeshField& q) {                     \
            forEach(*this,i)                                            \
                P[i] $ q[i];                                            \
            return *this;                                               \
        }
#define SOp($)                                                          \
        MeshField& operator $(const Scalar& q) {                        \
            forEach(*this,i)                                            \
                P[i] $ q;                                               \
            return *this;                                               \
        }
        Op(=)
        SOp(=)
        SOp(+=)
        SOp(-=)
        SOp(*=)
        SOp(/=)
#undef Op
#undef SOp

        /*Assignment from expressions*/
        template <class A>
        MeshField(const DVExpr<type,A>& p) {
            allocate();
            forEach(*this,i)
                P[i] = p[i];
        }

#define Op($)                                                           \
        template <class A>                                              \
        MeshField& operator $(const DVExpr<type,A>& q) {                \
            forEach(*this,i)                                            \
                P[i] $ q[i];                                            \
            return *this;                                               \
        }
        Op(=)
        Op(+=)
        Op(-=)
        Op(*=)
        Op(/=)  
#undef Op

        /*other member functions*/
        void calc_neumann(BCondition<type>*);
        void readInternal(std::istream&,IntVector*);
        void readBoundary(std::istream&);
        void writeInternal(std::ostream&,IntVector*);
        void writeBoundary(std::ostream&);
        void read(Int step);
        void write(Int step, IntVector* = 0);

        void norm(BaseField* pnorm) {
            *((MeshField<Scalar,entity>*)pnorm) = mag(*this);
        }
        /*read/write all fields*/
        static void readAll(Int step) {
            forEachIt(fields_, it) {
                if((*it)->access & READ)
                    (*it)->read(step);
            }
        }
        static void writeAll(Int step) {
            forEachIt(fields_, it) {
                if((*it)->access & WRITE)
                    (*it)->write(step);
            }
        } 
        static void removeAll() {
            fields_.clear();
            forEachIt(mem_pool,it)
                delete[] (*it);
            mem_pool.clear();
            n_alloc = 0;
        } 
        static int count_writable() {
            int count = 0;
            forEachIt(fields_, it) {
                if((*it)->access & WRITE)
                    count++;
            }
            return count;
        }
        static void writeVtkCellAll(std::ostream& os) {
            MeshField<type,CELL>* pf;
            forEachIt(fields_, it) {
                pf = *it;
                if(pf->access & WRITE) {
                    os << pf->fName <<" "<< TYPE_SIZE <<" "
                        << Mesh::gBCSfield << " double" << std::endl;
                    for(Int i = 0;i < Mesh::gBCSfield;i++)
                        os << (*pf)[i] << std::endl;
                    os << std::endl;
                }
            }
        }
        static void writeVtkVertexAll(std::ostream& os) {
            MeshField<type,VERTEX> vf;
            forEachIt(fields_, it) {
                if((*it)->access & WRITE) {
                    vf = cds(cds(*(*it)));
                    os << (*it)->fName <<" "<< TYPE_SIZE <<" "
                        << vf.size() << " double" << std::endl;
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
            forEachIt(fields_, it) {
                if((*it)->access & WRITE) {
                    vf = cds(cds(*(*it)));
                    vf_fields_->push_back(vf);
                }
            }
        }
        /*Store previous values*/
        MeshField* tstore;
        Int nstore;
        Int nstored;
        void initStore() {
            nstore = Controls::time_scheme - Controls::BDF1 + 1;
            tstore = new MeshField[nstore];
            access = ACCESS(int(access) | STOREPREV);
            for(Int i = 0;i < nstore;i++)
                tstore[i] = *this;
            nstored = 0;
        }
        void updateStore() {
            for(Int i = 0;i < nstore - 1;i++)
                tstore[nstore - i - 1] = tstore[nstore - i - 2];
            tstore[0] = *this;
            nstored++;
        }
        /*Time history*/
        static std::vector<std::ofstream*> tseries;
        static std::vector<MeshField*> tavgs;
        static std::vector<MeshField*> tstds;

        static void initTimeSeries() {
            MeshField<type,CELL>* pf;
            forEachIt(fields_, it) {
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
            forEachIt(fields_, it) {
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
                //update store
                if(pf->access & STOREPREV) {
                    pf->updateStore();
                }
            }
        }
        /*refine/unrefine field*/
        void refineField(Int step,IntVector& refineMap,IntVector& coarseMap) {
            forEach(coarseMap,i) {
                Int nchildren = coarseMap[i];
                Int id = coarseMap[i + 1];
                for(Int j = 1;j < nchildren;j++)
                    P[id] += P[coarseMap[(i + 1) + j]];
                P[id] /= nchildren;
                i += nchildren;
            }
            if(access & WRITE)
                write(step,&refineMap);
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
        /*Memory usage*/
        static void printUsage() {
            if(n_alloc_max) {
                std::cout << n_alloc_max
                    << " fields of vector of size " 
                    << TYPE_SIZE 
                    << " at " << entity
                    << std::endl;
            }
        }
};

#define forEachCellField(X)  {      \
    ScalarCellField::X;             \
    VectorCellField::X;             \
    STensorCellField::X;            \
    TensorCellField::X;             \
}
#define forEachFacetField(X)     {  \
    ScalarFacetField::X;            \
    VectorFacetField::X;            \
    STensorFacetField::X;           \
    TensorFacetField::X;            \
}
#define forEachVertexField(X)    {  \
    ScalarVertexField::X;           \
    VectorVertexField::X;           \
    STensorVertexField::X;          \
    TensorVertexField::X;           \
}
/** \name Typedef common mesh field types */
//@{
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
//@}

/** \name Static variables of MeshField */
//@{
template <class T,ENTITY E> 
std::list<MeshField<T,E>*> MeshField<T,E>::fields_;

template <class T,ENTITY E> 
std::list<T*> MeshField<T,E>::mem_pool;

template <class T,ENTITY E> 
Int MeshField<T,E>::n_alloc;

template <class T,ENTITY E> 
Int MeshField<T,E>::n_alloc_max;

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
//@}

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
    extern ScalarFacetField  fD;
    extern ScalarCellField   yWall;
    extern IntVector         FO;
    extern IntVector         FN; 

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
    void   calc_courant(const VectorCellField& U, Scalar dt);
    template <class type>
        void   scaleBCs(const MeshField<type,CELL>&, MeshField<type,CELL>&, Scalar);
}

namespace Prepare {
    void createFields(std::vector<std::string>& fields,Int step);
    Int  readFields(std::vector<std::string>& fields,Int step);
    void refineMesh(Int step, bool = false);
    void calcQOI(VectorCellField&);
    void initRefineThreshold();
    int decomposeMesh(Int);
    int mergeFields(Int);
}

/** central difference scheme */
template<class type>
MeshField<type,FACET> cds(const MeshField<type,CELL>& cF) {
    using namespace Mesh;
    MeshField<type,FACET> fF;
    forEach(fF,i) {
        fF[i] =  (cF[FO[i]] * (fI[i])) + (cF[FN[i]] * (1 - fI[i]));
    }
    return fF;
}

template<class type, class A>
MeshField<type,FACET> cds(const DVExpr<type,A>& expr) {
    return cds(MeshField<type,CELL>(expr));
}

/** upwind differencing scheme */
template<class type,class T3>
MeshField<type,FACET> uds(const MeshField<type,CELL>& cF,const MeshField<T3,FACET>& flux) {
    using namespace Mesh;
    MeshField<type,FACET> fF;
    forEach(fF,i) {
        if(dot(flux[i],T3(1)) >= 0) fF[i] = cF[FO[i]];
        else fF[i] = cF[FN[i]];
    }
    return fF;
}

template<class type, class T3, class A>
MeshField<type,FACET> uds(const DVExpr<type,A>& expr,const MeshField<T3,FACET>& flux) {
    return uds(MeshField<type,CELL>(expr),flux);
}

/** interpolate facet data to vertex data */
template<class type>
MeshField<type,VERTEX> cds(const MeshField<type,FACET>& fF) {
    using namespace Mesh;
    MeshField<type,VERTEX> vF;
    std::vector<Scalar> cnt;

    cnt.assign(vF.size(),Scalar(0));
    vF = type(0);

    forEach(fF,i) {
        Facet& f = gFacets[i];
        if(FN[i] < gBCSfield) {
            forEach(f,j) {
                Scalar dist = Scalar(1.0) / magSq(gVertices[f[j]] - fC[i]);
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

/** Integrate field operation */
template<class type>
MeshField<type,CELL> sum(const MeshField<type,FACET>& fF) {
    using namespace Mesh;
    MeshField<type,CELL> cF;
    cF = type(0);
    forEach(fF,i) {
        cF[FO[i]] += fF[i];
        cF[FN[i]] -= fF[i];
    }
    return cF;
}

template<class type, class A>
MeshField<type,CELL> sum(const DVExpr<type,A>& expr) {
    return sum(MeshField<type,FACET>(expr));
}

/** Scale boundary condition by a constant */
template <class type>
void Mesh::scaleBCs(const MeshField<type,CELL>& src, MeshField<type,CELL>& dest, Scalar psi) {
    using namespace Mesh;
    forEach(AllBConditions,i) {
        BCondition<type> *bc, *bc1;
        bc = static_cast<BCondition<type>*> (AllBConditions[i]);
        if(bc->fIndex == src.fIndex) {
            bc1 = new BCondition<type>(dest.fName);
            *bc1 = *bc;
            bc1->fIndex = dest.fIndex;

            if(bc1->cIndex == NEUMANN) {
                bc1->value *= psi;
            } else if(bc1->cIndex == ROBIN) {
                bc1->value *= psi;
                bc1->tvalue *= psi;
            } else if(bc1->cIndex == SYMMETRY ||
                    bc1->cIndex == CYCLIC ||
                    bc1->cIndex == RECYCLE) {
            } else {
                bc1->cIndex = GHOST;
            }

            dest.calc_neumann(bc1);
            AllBConditions.push_back(bc1);
        }
    }
}

/** Calculate neumann boundary condition from initial condition */
template <class T,ENTITY E> 
void MeshField<T,E>::calc_neumann(BCondition<T>* bc) {
    using namespace Mesh;

    Int h = Util::hash_function(bc->cname);
    if(h == CALC_NEUMANN) {
        /*calculate slope*/
        T slope = T(0);
        Int sz = bc->bdry->size();
        for(Int j = 0;j < sz;j++) {
            Int faceid = (*bc->bdry)[j];
            for(Int n = 0; n < DG::NPF;n++) {
                Int k = faceid * DG::NPF + n;
                Int c1 = FO[k];
                Int c2 = FN[k];
                if(!equal(cC[c1],cC[c2])) {
                    slope += ((*this)[c2] - (*this)[c1]) / 
                        mag(cC[c2] - cC[c1]);
                }
            }
        }
        slope /= sz;
        /*set to neumann*/
        bc->cname = "NEUMANN";
        bc->cIndex = NEUMANN;
        bc->value = slope;
    }
}
/* **********************************************
 *  Input - output operations
 * **********************************************/

/** Read internal field */
template <class T,ENTITY E> 
void MeshField<T,E>::readInternal(std::istream& is, IntVector* cMap) {
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
            for(Int i = 0;i < gBCSfield;i++) {
                Scalar R = mag((cC[i] - center) / radius);
                R = min(1.0,R);
                T val = value;
                val += (perterb / 2) * (Scalar(1.0) + cos(R * Constants::PI));
                (*this)[i] = val;
            }
        } else if(str == "gaussian") {
            Vector center,radius;
            T perterb;
            is >> value >> perterb >> center >> radius;
            for(Int i = 0;i < gBCSfield;i++) {
                Scalar R = mag((cC[i] - center) / radius);
                Scalar v = exp(-R*R);
                if(equal(v,Scalar(0))) v = 0;
                T val = value;
                val += perterb *  v;
                (*this)[i] = val;
            }
        } else if(str == "hydrostatic") {
            T p0;
            Scalar scale, expon;
            is >> p0 >> scale >> expon;
            for(Int i = 0;i < gBCSfield;i++) {
                Scalar gz = dot(cC[i],Controls::gravity);
                (*this)[i] = (p0) * pow(Scalar(1.0) + scale * gz, expon);
            }
        } else
            *this = value;
    } else {
        T temp;
        char symbol;
        is >> size >> symbol;
        for(int i = 0;i < size;i++) {
            is >> temp;
            if(cMap) 
                (*this)[(*cMap)[i]] = temp;
            else 
                (*this)[i] = temp;
        }
        is >> symbol;
    }
}

/** Read boundary field */
template <class T,ENTITY E> 
void MeshField<T,E>::readBoundary(std::istream& is) {
    using namespace Mesh;

    /*boundary field*/
    char c;
    BCondition<T>* bc;
    while((c = Util::nextc(is)) && isalpha(c)) {
        bc = new BCondition<T>(this->fName);
        is >> *bc;
        this->calc_neumann(bc);
        AllBConditions.push_back(bc);
    }
}

/** Read field */
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
    if(MP::printOn) {
        std::cout << "Reading " << fName 
            << step  << std::endl;
        std::cout.flush();
    }
    /*internal*/
    readInternal(is,0);
    readBoundary(is);

    /*update BCs*/
    applyExplicitBCs(*this,true,true);
}

/** Write internal field */
template <class T,ENTITY E> 
void MeshField<T,E>::writeInternal(std::ostream& os, IntVector* cMap) {
    using namespace Mesh;

    /*size*/
    os << "size " << sizeof(T) / sizeof(Scalar) << std::endl;

    /*internal field*/
    Int size;
    if(cMap)
        size = cMap->size();
    else
        size = (SIZE == gCells.size() * DG::NP) ? gBCSfield : SIZE;
    os << size << std::endl;
    os << "{" << std::endl;
    for(Int i = 0;i < size;i++) {
        if(cMap)
            os << (*this)[(*cMap)[i]] << std::endl;
        else
            os << (*this)[i] << std::endl;
    }
    os << "}" << std::endl;
}

/** Write boundary field */
template <class T,ENTITY E> 
void MeshField<T,E>::writeBoundary(std::ostream& os) {
    using namespace Mesh;

    /*boundary field*/
    BasicBCondition* bbc;
    BCondition<T>* bc;
    forEach(AllBConditions,i) {
        bbc = AllBConditions[i];
        if(bbc->fIndex == this->fIndex) {
            bc = static_cast<BCondition<T>*> (bbc);
            os << *bc << std::endl;
        }
    }
}

/** Write field */
template <class T,ENTITY E> 
void MeshField<T,E>::write(Int step, IntVector* cMap) {
    /*open*/
    std::stringstream path;
    path << fName << step;
    std::ofstream of(path.str().c_str());

    /*write*/
    writeInternal(of,cMap);
    writeBoundary(of);
}

/* ********************
 *   DG
 * ********************/
#include "dg.h"

/*********************************************************************************
 *                      matrix class defined on mesh                             
 *********************************************************************************/

/**
  Matrix defined on Mesh
 */
template <class T1, class T2 = Scalar, class T3 = T1> 
struct MeshMatrix {
    MeshField<T1,CELL>*   cF;    /**< Solution field X */
    MeshField<T2,CELL>    ap;    /**< Diagonal of the matrix */
    MeshField<T2,FACET>   an[2]; /**< Off-diagonals defined by owner/neighbor of face */
    MeshField<T2,CELLMAT> adg;   /**< In DG this is an NxN matrix tying the nodes in an element */
    MeshField<T3,CELL>    Su;    /**< Source field B */
    Int flags;                   /**< Flags for special matrix */

    /** Special matrix flag */
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
    MeshMatrix(MeshField<T1,CELL>* pcF) {
        cF = pcF;
        flags = (SYMMETRIC | DIAGONAL);
        ap = T2(0);
        an[0] = T2(0);
        an[1] = T2(0);
        Su = T3(0);
        adg = T2(0);
    }
    template<class A>
    MeshMatrix(const DVExpr<T1,A>& p) {
        cF = 0;
        flags = (SYMMETRIC | DIAGONAL);
        ap = T2(0);
        an[0] = T2(0);
        an[1] = T2(0);
        Su = p;
        adg = T2(0);
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
    void Fix(Int c,T1 value) {
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

/** \name Typedef common matrix types*/
//@{
typedef MeshMatrix<Scalar>  ScalarCellMatrix;
typedef MeshMatrix<Vector>  VectorCellMatrix;
typedef MeshMatrix<Tensor>  TensorCellMatrix;
typedef MeshMatrix<STensor> STensorCellMatrix;
//@}

/*************************************
 *  Asynchronous communication
 *************************************/

/** Class for asynchronous communication using MPI */
template <class T> 
class ASYNC_COMM {
    private:
        T* P;
        Int rcount;
        std::vector<MP::REQUEST> request;
        MeshField<T,CELL> recvbuf;
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
                Int buf_size = f.size() * NPF;

                //--fill send buffer
                forEach(f,j) {
                    Int faceid = f[j];
                    for(Int n = 0; n < NPF;n++) {
                        Int k = faceid * NPF + n;
                        sendbuf[(b.buffer_index + j) * NPF + n] = P[FO[k]]; 
                    }                                                           
                }   

                //--non-blocking send/recive
                MP::isend(&sendbuf[b.buffer_index * NPF],buf_size,
                        b.to,MP::FIELD_BLK,&request[rcount]);
                rcount++;
                MP::irecieve(&recvbuf[b.buffer_index * NPF],buf_size,
                        b.to,MP::FIELD_BLK,&request[rcount]);
                rcount++;
            }
        }
        void recv() {
            using namespace Mesh;
            using namespace DG;

            MP::waitall(rcount,&request[0]);

            //--copy from buffer to ghost cells
            forEach(gInterMesh,i) {
                interBoundary& b = gInterMesh[i];
                IntVector& f = *(b.f);

                forEach(f,j) {
                    Int faceid = f[j];
                    for(Int n = 0; n < NPF;n++) {
                        Int k = faceid * NPF + n;
                        P[FN[k]] = recvbuf[(b.buffer_index + j) * NPF + n]; 
                    }                                                           
                }
            }
        }
};
/* ********************************
 *  Tenosor-Product approach
 * ********************************/
#define TensorProduct_(Q,P,tr,$) {                                          \
    for(Int ci = 0;ci < gBCS;ci++) {                                        \
        forEachLgl(ii,jj,kk) {                                              \
            Int index1 = INDEX4(ci,ii,jj,kk);                               \
            T3 val(Scalar(0));                                              \
                                                                            \
            forEachLglX(i) {                                                \
                Int index2 = INDEX4(ci,i,jj,kk);                            \
                Int indexm = ci * NPMAT +                                   \
                    (tr ? INDEX_TX(ii,jj,kk,i) : INDEX_X(ii,jj,kk,i));      \
                val += Q[index2] * P.adg[indexm];                           \
            }                                                               \
            forEachLglY(j) if(j != jj) {                                    \
                Int index2 = INDEX4(ci,ii,j,kk);                            \
                Int indexm = ci * NPMAT +                                   \
                (tr ? INDEX_TY(ii,jj,kk,j) : INDEX_Y(ii,jj,kk,j));          \
                val += Q[index2] * P.adg[indexm];                           \
            }                                                               \
            forEachLglZ(k) if(k != kk) {                                    \
                Int index2 = INDEX4(ci,ii,jj,k);                            \
                Int indexm = ci * NPMAT +                                   \
                (tr ? INDEX_TZ(ii,jj,kk,k) : INDEX_Z(ii,jj,kk,k));          \
                val += Q[index2] * P.adg[indexm];                           \
            }                                                               \
                                                                            \
            r[index1] $ val;                                                \
        }                                                                   \
    }                                                                       \
}

#define TensorProduct(Q,P)  TensorProduct_(Q,P,false,-=)
#define TensorProductT(Q,P) TensorProduct_(Q,P,true,-=)
#define TensorProductM(Q,P) TensorProduct_(Q,P,false,+=)

/* *******************************
 * matrix - vector products
 * *******************************/

/** matrix - vector product = A * x */
template <class T1, class T2, class T3> 
MeshField<T1,CELL> mul (const MeshMatrix<T1,T2,T3>& p,const MeshField<T1,CELL>& q, const bool sync = false) {
    using namespace Mesh;
    using namespace DG;
    MeshField<T3,CELL> r;
    Int c1,c2;
    ASYNC_COMM<T1> comm(&q[0]);

    if(sync) comm.send();

    r = q * p.ap;

    if(NPMAT) {
        TensorProduct(q,p);
    }

    forEach(FN,f) {
        c2 = FN[f];
        if(c2 >= gBCSfield) continue;
        c1 = FO[f];
        r[c1] -= q[c2] * p.an[1][f];
        r[c2] -= q[c1] * p.an[0][f];
    }

    if(sync) comm.recv();

    forEach(FN,f) {
        c2 = FN[f];
        if(c2 < gBCSfield) continue;
        c1 = FO[f];
        r[c1] -= q[c2] * p.an[1][f];
    }

    return r;
}

template <class T1, class T2, class T3, class A> 
MeshField<T1,CELL> mul (const MeshMatrix<T1,T2,T3>& p,const DVExpr<T1,A>& q, const bool sync = false) {
    return mul(p,MeshField<T1,CELL>(q),sync);
}

/** matrix transopose - vector product = A^T * x */
template <class T1, class T2, class T3> 
MeshField<T1,CELL> mult (const MeshMatrix<T1,T2,T3>& p,const MeshField<T1,CELL>& q, const bool sync = false) {
    using namespace Mesh;
    using namespace DG;
    MeshField<T3,CELL> r;
    Int c1,c2;
    ASYNC_COMM<T1> comm(&q[0]);

    if(sync) comm.send();

    r = q * p.ap;

    if(NPMAT) {
        TensorProductT(q,p);
    }

    forEach(FN,f) {
        c2 = FN[f];
        if(c2 >= gBCSfield) continue;
        c1 = FO[f];
        r[c2] -= q[c1] * p.an[1][f];
        r[c1] -= q[c2] * p.an[0][f];
    }

    if(sync) comm.recv();

    forEach(FN,f) {
        c2 = FN[f];
        if(c2 < gBCSfield) continue;
        c1 = FO[f];
        r[c1] -= q[c2] * p.an[0][f];
    }

    return r;
}
/** calculate right-hand-side sum = b - (L + U) * x */
template <class T1, class T2, class T3> 
MeshField<T1,CELL> getRHS(const MeshMatrix<T1,T2,T3>& p, const bool sync = false) {
    using namespace Mesh;
    using namespace DG;
    MeshField<T3,CELL> r;
    MeshField<T1,CELL>& q = (*p.cF);
    Int c1,c2;
    ASYNC_COMM<T1> comm(&q[0]);

    if(sync) comm.send();

    r = p.Su;

    if(NPMAT) {
        TensorProductM(q,p);
    }

    forEach(FN,f) {
        c2 = FN[f];
        if(c2 >= gBCSfield) continue;
        c1 = FO[f];
        r[c1] += q[c2] * p.an[1][f];
        r[c2] += q[c1] * p.an[0][f];
    }

    if(sync) comm.recv();

    forEach(FN,f) {
        c2 = FN[f];
        if(c2 < gBCSfield) continue;
        c1 = FO[f];
        r[c1] += q[c2] * p.an[1][f];
    }

    return r;
}

/* ***************************
 * Apply boundary conditions
 * ***************************/

/** Apply implicit boundary conditions */
template <class T1, class T2, class T3> 
void applyImplicitBCs(const MeshMatrix<T1,T2,T3>& M) {
    using namespace Mesh;
    MeshField<T1,CELL>& cF = *M.cF;
    BasicBCondition* bbc;
    BCondition<T1>* bc;

    /*boundary conditions*/
    forEach(AllBConditions,i) {
        bbc = AllBConditions[i];
        if(bbc->fIndex == cF.fIndex) {
            if(bbc->cIndex == GHOST)
                continue;

            bc = static_cast<BCondition<T1>*> (bbc);
            Int sz = bc->bdry->size();
            if(sz == 0) continue;

            for(Int j = 0;j < sz;j++) {
                Int faceid = (*bc->bdry)[j];
                for(Int n = 0; n < DG::NPF;n++) {
                    Int k = faceid * DG::NPF + n;

                    Int c1 = FO[k];
                    Int c2 = FN[k];
                    /*break connection with boundary cells*/
                    if(bc->cIndex == NEUMANN || bc->cIndex == SYMMETRY ||
                            bc->cIndex == CYCLIC || bc->cIndex == RECYCLE) {
                        M.ap[c1] -= M.an[1][k];
                        M.Su[c1] += M.an[1][k] * (cF[c2] - cF[c1]);
                        M.an[1][k] = 0;
                    } else if(bc->cIndex == ROBIN) {
                        Vector dv = cC[c2] - cC[c1];
                        M.ap[c1] -= (1 - bc->shape) * M.an[1][k];
                        M.Su[c1] += M.an[1][k] * (bc->shape * bc->value + 
                                (1 - bc->shape) * bc->tvalue * mag(dv));
                        M.an[1][k] = 0;
                    } else {
                        M.Su[c1] += M.an[1][k] * cF[c2];
                        M.an[1][k] = 0;
                    }
                }
            }
        }
    }
}

/** Apply explicit boundary conditions */
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
    /*update ghost cells*/
    bool sync = (update_ghost && gInterMesh.size());
    ASYNC_COMM<T> comm(&cF[0]);

    if(sync) comm.send();

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
                    Int c1 = FO[k];
                    Int c2 = FN[k];
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
                        Int c11 = FO[k1];
                        Int c22 = FN[k1];
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

    if(sync) comm.recv();
}

/** Fill boundary values from internals */
template<class T>
const MeshField<T,CELL>& fillBCs(const MeshField<T,CELL>& cF, const bool sync = false, const Int bind = 0) {
    using namespace Mesh;
    forEachS(gCells,i,gBCS) {
        Int faceid = gCells[i][0];
        for(Int n = 0; n < DG::NPF;n++) {
            Int k = faceid * DG::NPF + n;
            cF[FN[k]] = cF[FO[k]];
        }
    }

    //special: set neumann bcs to dirchlet for grad(i)
    if(bind) {
        forEach(AllBConditions,i) {
            BasicBCondition* bbc = AllBConditions[i];
            if(bbc->fIndex == bind) {
                BCondition<T>* bc = static_cast<BCondition<T>*> (bbc);
                Int sz = bc->bdry->size();
                for(Int j = 0;j < sz;j++) {
                    Int faceid = (*bc->bdry)[j];
                    for(Int n = 0; n < DG::NPF;n++) {
                        Int k = faceid * DG::NPF + n;
                        Int c2 = FN[k];
                        if(bc->cIndex == NEUMANN) {
                            cF[c2] = bc->value;
                        } else if(bc->cIndex == SYMMETRY) {
                            cF[c2] = T(0);
                        }
                    }
                }
            }
        }
    }

    if(gInterMesh.size()) {
        ASYNC_COMM<T> comm(&cF[0]);
        comm.send();
        comm.recv();
    }
    return cF;
}
/* *******************************
 * Linearized source term
 * *******************************/

/** Add implicit source terms*/
template <class T1, class T2, class T3> 
MeshMatrix<T1,T2,T3> src(MeshField<T1,CELL>& cF,const MeshField<T3,CELL>& Su, const MeshField<T2,CELL>& Sp) {
    MeshMatrix<T1,T2,T3> m;
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

/** Add explicit source term */
template<class type>
MeshField<type,CELL> srcf(const MeshField<type,CELL>& Su) {
    return (Su * Mesh::cV);
}

#define srci(x) (x)

/***********************************************
 * Gradient field operation.
 ***********************************************/

/*Explicit*/
#define GRADD(im,jm,km) {                           \
    Int index1 = INDEX4(ci,im,jm,km);               \
    Vector dpsi_ij;                                 \
    DPSI(dpsi_ij,im,jm,km);                         \
    dpsi_ij = dot(Jin,dpsi_ij);                     \
    r[index1] -= mul(dpsi_ij,p[index]);             \
}

#define GRAD(T1,T2)                                                                             \
    inline MeshField<T1,CELL> gradf(const MeshField<T2,CELL>& p) {                              \
        using namespace Mesh;                                                                   \
        using namespace DG;                                                                     \
        MeshField<T1,CELL> r;                                                                   \
                                                                                                \
        r = sum(mul(fN,cds(p)));                                                                \
                                                                                                \
        if(NPMAT) {                                                                             \
            for(Int ci = 0; ci < gBCS;ci++) {                                                   \
                forEachLgl(ii,jj,kk) {                                                          \
                    Int index = INDEX4(ci,ii,jj,kk);                                            \
                    Tensor Jin = Jinv[index] * cV[index];                                       \
                    forEachLglX(i) GRADD(i,jj,kk);                                              \
                    forEachLglY(j) if(j != jj) GRADD(ii,j,kk);                                  \
                    forEachLglZ(k) if(k != kk) GRADD(ii,jj,k);                                  \
                }                                                                               \
            }                                                                                   \
        }                                                                                       \
                                                                                                \
        fillBCs(r,false,p.fIndex);                                                              \
                                                                                                \
        return r;                                                                               \
    }

GRAD(Vector,Scalar);
GRAD(Tensor,Vector);
#undef GRAD
#undef GRADD

#define gradi(x) (gradf(x)  / Mesh::cV)

/**
  Compute numerical flux
 */
template<class T1, class T2, class T3, class T4>
void numericalFlux(MeshMatrix<T1,T2,T3>& m, const MeshField<T4,FACET>& flux, const MeshField<T4,CELL>* muc = 0) {
    using namespace Controls;
    using namespace Mesh;
    using namespace DG;

    /*compute surface integral*/
    MeshField<T1,CELL>& cF = *m.cF;
    T4 F;
    Scalar G;

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
        else if(!muc)
            gamma = Scalar(1);
        else if(convection_scheme == HYBRID) {
            MeshField<T4,FACET> mu = cds(*muc);
            forEach(gFacets,faceid) {
                for(Int n = 0; n < NPF;n++) {
                    Int k = faceid * NPF + n;
                    //compare F and D
                    T4 D = fD[k] * mu[k];
                    F = flux[k];
                    if(dot(F,T4(1)) < 0) {
                        if(dot(F * fI[k] + D,T4(1)) >= 0) gamma[k] = 1;
                        else gamma[k] = 0;
                    } else {
                        if(dot(F * (1 - fI[k]) - D,T4(1)) > 0) gamma[k] = 0;
                        else gamma[k] = 1;
                    }
                }
            }
        }
        forEach(flux,i) {
            F = flux[i];
            G = gamma[i];
            m.an[0][i] = ((G) * (-F * (  fI[i]  )) + (1 - G) * (-max( F,T4(0))));
            m.an[1][i] = ((G) * ( F * (1 - fI[i])) + (1 - G) * (-max(-F,T4(0))));
            m.ap[FO[i]] += m.an[0][i];
            m.ap[FN[i]] += m.an[1][i];
        }
        /*deferred correction*/
    } else {
        forEach(flux,i) {
            F = flux[i];
            m.an[0][i] = -max( F,T4(0));
            m.an[1][i] = -max(-F,T4(0));
            m.ap[FO[i]] += m.an[0][i];
            m.ap[FN[i]] += m.an[1][i];
        }

        MeshField<T1,FACET> corr;
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
            /**
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
            MeshField<T1,FACET> q,r,phiDC,phiCU;
            ScalarFacetField uFI;
            {
                MeshField<T4,FACET> nflux = T4(0)-flux;
                phiDC = uds(cF,nflux) - uds(cF,flux);
                forEach(phiDC,i) {
                    if(dot(flux[i],T4(1)) >= 0) G = fI[i];
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
                if(equal(phiDC[i] * (1 - uFI[i]),T1(0)))
                    r[i] = T1(0);
            }
            /*TVD schemes*/
            if(convection_scheme == VANLEER) {
                q = (r+fabs(r)) / (1+r);
            } else if(convection_scheme == VANALBADA) {
                q = (r+r*r) / (1+r*r);
            } else if(convection_scheme == MINMOD) {
                q = max(T1(0),min(r,T1(1)));
            } else if(convection_scheme == SUPERBEE) {
                q = max(min(r,T1(2)),min(2*r,T1(1)));
                q = max(q,T1(0));
            } else if(convection_scheme == SWEBY) {
                Scalar beta = 2;
                q = max(min(r,T1(beta)),min(beta*r,T1(1)));
                q = max(q,T1(0));
            } else if(convection_scheme == QUICKL) {
                q = min(2*r,(3+r)/4);
                q = min(q,T1(2));
                q = max(q,T1(0));
            } else if(convection_scheme == UMIST) {
                q = min(2*r,(3+r)/4);
                q = min(q,(1+3*r)/4);
                q = min(q,T1(2));
                q = max(q,T1(0));
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
}

/**
  Implicit gradient operator
 */
template<class T1, class T2>
MeshMatrix<T1,T2,T2> grad(MeshField<T1,CELL>& cF,const MeshField<T1,CELL>& flux_cell,
        const MeshField<T2,FACET>& flux) {
    using namespace Mesh;
    using namespace DG;
    MeshMatrix<T1,T2,T2> m;
    m.cF = &cF;
    m.flags = 0;
    m.Su = T1(0);
    m.ap = Scalar(0);
    m.adg = Scalar(0);

    /*compute volume integral*/
    if(NPMAT) {
        for(Int ci = 0; ci < gBCS;ci++) {
            forEachLgl(ii,jj,kk) {
                Int index = INDEX4(ci,ii,jj,kk);
                Tensor Jr = Jinv[index];
                T1 F = flux_cell[index] * cV[index];

#define GRADD(im,jm,km) {                                   \
    Int index2 = INDEX4(ci,im,jm,km);                       \
    Vector dpsi_ij;                                         \
    DPSI(dpsi_ij,im,jm,km);                                 \
    T2 val = -dot(Jr,dpsi_ij) * F;                          \
    if(index == index2) {                                   \
        m.ap[index] = -val;                                 \
        m.adg[indexm] = 0;                                  \
    } else {                                                \
        m.adg[indexm] = val;                                \
    }                                                       \
}
                forEachLglX(i) {
                    Int indexm = ci * NPMAT + INDEX_TX(ii,jj,kk,i);
                    GRADD(i,jj,kk);
                }
                forEachLglY(j) if(j != jj) {
                    Int indexm = ci * NPMAT + INDEX_TY(ii,jj,kk,j);
                    GRADD(ii,j,kk);
                }
                forEachLglZ(k) if(k != kk) {
                    Int indexm = ci * NPMAT + INDEX_TZ(ii,jj,kk,k);
                    GRADD(ii,jj,k);
                }
#undef GRADD

            }
        }
    }
    /*compute surface integral*/
    numericalFlux(m,flux);
    
    return m;
}
/* ***************************************************
 * Divergence field operation
 * ***************************************************/ 

/* Explicit */
#define DIVD(im,jm,km) {                            \
    Int index1 = INDEX4(ci,im,jm,km);               \
    Vector dpsi_ij;                                 \
    DPSI(dpsi_ij,im,jm,km);                         \
    dpsi_ij = dot(Jin,dpsi_ij);                     \
    r[index1] -= dot(p[index],dpsi_ij);             \
}

#define DIV(T1,T2)                                                                              \
    inline MeshField<T1,FACET> flx(const MeshField<T2,CELL>& p) {                               \
        return dot(cds(p),Mesh::fN);                                                            \
    }                                                                                           \
    inline MeshField<T1,CELL> divf(const MeshField<T2,CELL>& p) {                               \
        using namespace Mesh;                                                                   \
        using namespace DG;                                                                     \
        MeshField<T1,CELL> r;                                                                   \
                                                                                                \
        r = sum(dot(cds(p),fN));                                                                \
                                                                                                \
        if(NPMAT) {                                                                             \
            for(Int ci = 0; ci < gBCS;ci++) {                                                   \
                forEachLgl(ii,jj,kk) {                                                          \
                    Int index = INDEX4(ci,ii,jj,kk);                                            \
                    Tensor Jin = Jinv[index] * cV[index];                                       \
                    forEachLglX(i) DIVD(i,jj,kk);                                               \
                    forEachLglY(j) if(j != jj) DIVD(ii,j,kk);                                   \
                    forEachLglZ(k) if(k != kk) DIVD(ii,jj,k);                                   \
                }                                                                               \
            }                                                                                   \
        }                                                                                       \
                                                                                                \
        fillBCs(r);                                                                             \
                                                                                                \
        return r;                                                                               \
    }

DIV(Scalar,Vector);
DIV(Vector,Tensor);
#undef DIV
#undef DIVD

#define divi(x)  (divf(x)   / Mesh::cV)

/** 
  Implicit divergence operator
 */
template<class type>
MeshMatrix<type> div(MeshField<type,CELL>& cF,const VectorCellField& flux_cell,
        const ScalarFacetField& flux,const ScalarCellField* muc = 0) {

    using namespace Mesh;
    using namespace DG;
    MeshMatrix<type> m;
    m.cF = &cF;
    m.flags = 0;
    m.Su = type(0);
    m.ap = Scalar(0);
    m.adg = Scalar(0);

    /*compute volume integral*/
    if(NPMAT) {
        for(Int ci = 0; ci < gBCS;ci++) {
            forEachLgl(ii,jj,kk) {
                Int index = INDEX4(ci,ii,jj,kk);
                Vector Jr = dot(Jinv[index],flux_cell[index]) * cV[index];

#define DIVD(im,jm,km) {                                    \
    Int index2 = INDEX4(ci,im,jm,km);                       \
    Vector dpsi_ij;                                         \
    DPSI(dpsi_ij,im,jm,km);                                 \
    Scalar val = -dot(Jr,dpsi_ij);                          \
    if(index == index2) {                                   \
        m.ap[index] = -val;                                 \
        m.adg[indexm] = 0;                                  \
    } else {                                                \
        m.adg[indexm] = val;                                \
    }                                                       \
}
                forEachLglX(i) {
                    Int indexm = ci * NPMAT + INDEX_TX(ii,jj,kk,i);
                    DIVD(i,jj,kk);
                }
                forEachLglY(j) if(j != jj) {
                    Int indexm = ci * NPMAT + INDEX_TY(ii,jj,kk,j);
                    DIVD(ii,j,kk);
                }
                forEachLglZ(k) if(k != kk) {
                    Int indexm = ci * NPMAT + INDEX_TZ(ii,jj,kk,k);
                    DIVD(ii,jj,k);
                }
#undef DIVD

            }
        }
    }
    /*compute surface integral*/
    numericalFlux(m,flux,muc);
    
    return m;
}

/* ***********************************************************
 * Laplacian field operation 
 *    It computes div( mu * grad(p) ), hence it is not really
 *    the laplacian operation for a nonconstant mu
 * ***********************************************************/

/*Explicit*/
template<class type, ENTITY entity>
MeshField<type,entity> lapf(MeshField<type,entity>& cF,const MeshField<Scalar,entity>& mu) {
    return divf(mu * gradi(cF));
}

#define lapi(x,y) (lapf(x,y)  / Mesh::cV)

/**
  Implicit laplacian operator
 */
template<class type>
MeshMatrix<type> lap(MeshField<type,CELL>& cF,const ScalarCellField& muc, const bool penalty = false) {

    using namespace Controls;
    using namespace Mesh;
    using namespace DG;

    MeshMatrix<type> m;
    m.cF = &cF;
    m.flags = m.SYMMETRIC;
    m.Su = type(0);
    m.ap = Scalar(0);
    m.adg = Scalar(0);

    /* diffusion or penalty term */
    {
        ScalarFacetField mu = cds(muc);
        forEach(fN,i) {
            Int c1 = FO[i];
            Int c2 = FN[i];
            /*coefficients*/
            if(penalty || !NPMAT) {
                m.an[0][i] = fD[i] * mu[i];
                m.an[1][i] = fD[i] * mu[i];
            } else{
                m.an[0][i] = mu[i];
                m.an[1][i] = mu[i];
            }
            m.ap[c1]  += m.an[0][i];
            m.ap[c2]  += m.an[1][i];
        }
    }

    if(NPMAT) {

        /*compute volume integral*/
        for(Int ci = 0; ci < gBCS;ci++) {
            forEachLgl(ii,jj,kk) {
                Int index = INDEX4(ci,ii,jj,kk);
                Tensor Jin = Jinv[index];

#define H(in,jn,kn) {                               \
    Int index2 = INDEX4(ci,in,jn,kn);               \
    Vector dpsi_jk;                                 \
    DPSI(dpsi_jk,in,jn,kn);                         \
    dpsi_jk = dot(Jin,dpsi_jk);                     \
    Scalar val = -dot(dpsi_ik,dpsi_jk);             \
    if(index1 == index2) {                          \
        m.ap[index1] += -val;                       \
    } else {                                        \
        m.adg[indexm]  += val;                      \
        m.adg[indexmt] += val;                      \
    }                                               \
}

#define LAPD(im,jm,km) {                                    \
    Int index1 = INDEX4(ci,im,jm,km);                       \
    Vector dpsi_ik,dp1;                                     \
    DPSI(dpsi_ik,im,jm,km);                                 \
    dpsi_ik = dot(Jin,dpsi_ik) * cV[index] * muc[index];    \
    forEachLglX(i2) if(i2 >= im) {                          \
        Int indexm  = ci * NPMAT + INDEX_TX(im,jm,km,i2);   \
        Int indexmt = ci * NPMAT + INDEX_X(im,jm,km,i2);    \
        H(i2,jm,km);                                        \
    }                                                       \
    forEachLglY(j2) if(j2 > jm) {                           \
        Int indexm  = ci * NPMAT + INDEX_TY(im,jm,km,j2);   \
        Int indexmt = ci * NPMAT + INDEX_Y(im,jm,km,j2);    \
        H(im,j2,km);                                        \
    }                                                       \
    forEachLglZ(k2) if(k2 > km) {                           \
        Int indexm  = ci * NPMAT + INDEX_TZ(im,jm,km,k2);   \
        Int indexmt = ci * NPMAT + INDEX_Z(im,jm,km,k2);    \
        H(im,jm,k2);                                        \
    }                                                       \
}
                forEachLglX(i) LAPD(i,jj,kk);
                forEachLglY(j) if(j != jj) LAPD(ii,j,kk);
                forEachLglZ(k) if(k != kk) LAPD(ii,jj,k);
#undef H
#undef LAPD

            }
        }

        /* compute explicit term */
        {
            m.Su += sum(dot(cds(muc * gradi(cF)),fN));
        }
    
    } else {
        /*non-orthogonality*/
        if(nonortho_scheme != NONE) {
            VectorFacetField K;
            forEach(fN,i) {
                Int c1 = FO[i];
                Int c2 = FN[i];
                Vector dv = cC[c2] - cC[c1];
                K[i] = fN[i] - fD[i] * dv;
            }
    
            MeshField<type,FACET> r = dot(cds(muc * gradi(cF)),K);
            forEach(r,i) {
                Int c1 = FO[i];
                Int c2 = FN[i];
                type res = m.an[0][i] * (cF[c2] - cF[c1]);
                if(mag(r[i]) > Scalar(0.5) * mag(res)) 
                    r[i] = Scalar(0.5) * res;
            }
            m.Su = sum(r);
        }
    }
    
    /*end*/
    return m;
}

/* *******************************
 * Temporal derivative
 * *******************************/
#define PREV(k) (rho ? (cF.tstore[k] * rho->tstore[k]) : cF.tstore[k])

/** First derivative with respect to time */
template<class type>
MeshMatrix<type> ddt(MeshField<type,CELL>& cF,ScalarCellField* rho) {
    MeshMatrix<type> m;
    m.cF = &cF;
    m.flags |= (m.SYMMETRIC | m.DIAGONAL);

    //BDF methods
    if(Controls::time_scheme == Controls::BDF1) {
        m.ap = (-1.0 / Controls::dt) * Mesh::cV;
        m.Su = (PREV(0)) * m.ap;
    } else if(Controls::time_scheme == Controls::BDF2) {
        m.ap = (-3.0 / (2.0 * Controls::dt)) * Mesh::cV;
        m.Su = ((4.0 * PREV(0) - PREV(1)) / 3.0) * m.ap;
    } else if(Controls::time_scheme == Controls::BDF3) {
        m.ap = (-11.0 / (6.0 * Controls::dt)) * Mesh::cV;
        m.Su = ((18.0 * PREV(0) - 9.0 * PREV(1) + 2.0 * PREV(2)) / 11.0) * m.ap;
    } else if(Controls::time_scheme == Controls::BDF4) {
        m.ap = (-25.0 / (12.0 * Controls::dt)) * Mesh::cV;
        m.Su = ((48.0 * PREV(0) - 36.0 * PREV(1) + 16.0 * PREV(2) - 3.0 * PREV(3)) / 25.0) * m.ap;
    } else if(Controls::time_scheme == Controls::BDF5) {
        m.ap = (-137.0 / (60.0 * Controls::dt)) * Mesh::cV;
        m.Su = ((300.0 * PREV(0) - 300.0 * PREV(1) + 200.0 * PREV(2) - 75.0 * PREV(3) + 12.0 * PREV(4)) / 137.0) * m.ap;
    } else if(Controls::time_scheme == Controls::BDF6) {
        m.ap = (-147.0 / (60.0 * Controls::dt)) * Mesh::cV;
        m.Su = ((360.0 * PREV(0) - 450.0 * PREV(1) + 400.0 * PREV(2) - 225.0 * PREV(3) + 72.0 * PREV(4) - 10.0 * PREV(5)) / 147.0) * m.ap;
    }
    if(rho) m.ap *= (*rho);

    //others
    m.adg = Scalar(0);
    m.an[0] = Scalar(0);
    m.an[1] = Scalar(0);

    return m;
}

/** Second derivative with respect to time */
template<class type>
MeshMatrix<type> ddt2(MeshField<type,CELL>& cF, ScalarCellField* rho) {
    MeshMatrix<type> m;
    m.cF = &cF;
    m.flags |= (m.SYMMETRIC | m.DIAGONAL);

    //BDF method
    if(Controls::time_scheme == Controls::BDF2) {
        m.ap = (-1.0 / (Controls::dt * Controls::dt)) * Mesh::cV;
        m.Su = (2.0 * PREV(0) - PREV(1)) * m.ap;
    } else if(Controls::time_scheme == Controls::BDF3) {
        m.ap = (-2.0 / (Controls::dt * Controls::dt)) * Mesh::cV;
        m.Su = ((5.0 * PREV(0) - 4.0 * PREV(1) + PREV(2)) / 2.0) * m.ap;
    }
    if(rho) m.ap *= (*rho);

    //others
    m.adg = Scalar(0);
    m.an[0] = Scalar(0);
    m.an[1] = Scalar(0);

    return m;
}   

#undef PREV

/** Time stepper */
template<int order, class type>
void addTemporal(MeshMatrix<type>& M,Scalar cF_UR,ScalarCellField* rho = 0) {
    using namespace Controls;
    if(state == STEADY)
        M.Relax(cF_UR);
    else {
        //store previous values
        if(!(M.cF->access & STOREPREV)) 
            M.cF->initStore();

        //Multistage Runge-Kutta for linearized (constant jacobian)
        if(!equal(implicit_factor,1)) {
            MeshField<type, CELL> k1 = M.Su - mul(M,  M.cF->tstore[0]);
            if(runge_kutta == 1) {
                MeshField<type, CELL> val = k1  * (1 - implicit_factor);
                M = M * (implicit_factor) + val;
            } else {
                ScalarCellField mdt = Controls::dt / (Mesh::cV);
                if (runge_kutta == 2) {
                    MeshField<type, CELL> k2 = k1 - mul(M, k1 * mdt);
                    MeshField<type, CELL> val = ((k2 + k1) / 2) * (1 - implicit_factor);
                    M = M * (implicit_factor) + val;
                } else if (runge_kutta == 3) {
                    MeshField<type, CELL> k2 = k1 - mul(M, k1 * mdt / 2);
                    MeshField<type, CELL> k3 = k1 - mul(M, (2 * k2 - k1) * mdt);
                    MeshField<type, CELL> val = ((k3 + 4 * k2 + k1) / 6) * (1 - implicit_factor);
                    M = M * (implicit_factor) + val;
                } else if (runge_kutta == 4) {
                    MeshField<type, CELL> k2 = k1 - mul(M, k1 * mdt / 2);
                    MeshField<type, CELL> k3 = k1 - mul(M, k2 * mdt / 2);
                    MeshField<type, CELL> k4 = k1 - mul(M, k3 * mdt);
                    MeshField<type, CELL> val = ((k4 + 2 * k3 + 2 * k2 + k1) / 6) * (1 - implicit_factor);
                    M = M * (implicit_factor) + val;
                }
            }
            if(equal(implicit_factor,0))
                M.flags |= (M.SYMMETRIC | M.DIAGONAL);
        }

        //first or second derivative
        if(order == 1)
            M += ddt(*M.cF,rho);
        else
            M += ddt2(*M.cF,rho);
    }
}
/* ************************************************
 * Form transport equation
 *************************************************/
template<class type>
MeshMatrix<type> diffusion(MeshField<type,CELL>& cF,
        const ScalarCellField& mu, Scalar cF_UR, ScalarCellField* rho = 0) {
    MeshMatrix<type> M = -lap(cF,mu,true);
    M.cF = &cF;
    addTemporal<1>(M,cF_UR,rho);
    return M;
}
template<class type>
MeshMatrix<type> diffusion(MeshField<type,CELL>& cF,
        const ScalarCellField& mu, Scalar cF_UR, 
        const MeshField<type,CELL>& Su,const ScalarCellField& Sp, ScalarCellField* rho = 0) {
    MeshMatrix<type> M = -lap(cF,mu,true) - src(cF,Su,Sp);
    M.cF = &cF;
    addTemporal<1>(M,cF_UR,rho);
    return M;
}
template<class type>
MeshMatrix<type> convection(MeshField<type,CELL>& cF,const VectorCellField& Fc,
        const ScalarFacetField& F, Scalar cF_UR, ScalarCellField* rho = 0) {
    MeshMatrix<type> M = div(cF,Fc,F);
    addTemporal<1>(M,cF_UR,rho);
    return M;
}
template<class type>
MeshMatrix<type> convection(MeshField<type,CELL>& cF,const VectorCellField& Fc,
        const ScalarFacetField& F, Scalar cF_UR, 
        const MeshField<type,CELL>& Su,const ScalarCellField& Sp, ScalarCellField* rho = 0) {
    MeshMatrix<type> M = div(cF,Fc,F) - src(cF,Su,Sp);
    addTemporal<1>(M,cF_UR,rho);
    return M;
}
template<class type>
MeshMatrix<type> transport(MeshField<type,CELL>& cF,const VectorCellField& Fc,const ScalarFacetField& F,
        const ScalarCellField& mu, Scalar cF_UR, ScalarCellField* rho = 0) {
    MeshMatrix<type> M = div(cF,Fc,F,&mu) - lap(cF,mu);
    addTemporal<1>(M,cF_UR,rho);
    return M;
}
template<class type>
MeshMatrix<type> transport(MeshField<type,CELL>& cF,const VectorCellField& Fc,const ScalarFacetField& F,
        const ScalarCellField& mu, Scalar cF_UR, 
        const MeshField<type,CELL>& Su,const ScalarCellField& Sp, ScalarCellField* rho = 0) {
    MeshMatrix<type> M = div(cF,Fc,F,&mu) - lap(cF,mu) - src(cF,Su,Sp);
    addTemporal<1>(M,cF_UR,rho);
    return M;
}
/* ********************
 *        End
 * ********************/
#endif
