#ifndef __MESH_H
#define __MESH_H

#include "tensor.h"
#include "util.h"

/** Basic building blocks (entities) */
enum ENTITY {
    CELL,    /**< A cell (volume), each LGL node gets one cell in DG. */
    FACET,   /**< A polygonal face. */
    VERTEX,  /**< A node. */
    CELLMAT  /**< An NXN matrix for each element. */
};

/** \name Indexing of entities is done via integer indices (not pointers). */
//@{
typedef std::vector<Int>      IntVector;
typedef std::vector<Scalar>   ScalarVector;
typedef std::vector<Vector>   VectorVector;
typedef std::vector<bool>     BoolVector;

typedef Vector           Vertex;
typedef IntVector        Facet;  
typedef IntVector        Cell; 

typedef std::vector<Vertex>   Vertices;
typedef std::vector<Facet>    Facets;
typedef std::vector<Cell>     Cells; 
typedef std::map<std::string,IntVector> Boundaries;
//@}

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

/**
Mesh data structure
*/
namespace Mesh {
    /** Inter-processor boundary */
    struct interBoundary {
        IntVector* f;
        Int from;
        Int to;
        Int buffer_index;
    };
    
    /** Boundary patch */
    struct Patch {
        Int from;
        Int to;
        Vector N;
        Vector C;
        Patch() {
            from = 0;
            to = 0;
        }
    };
    
    /** AMR tree node */
    struct Node {
        Int id;
        Int cid;
        Int nchildren;
        Node() {
            id = 0;
            cid = 0;
            nchildren = 0;
        }
        friend std::ostream& operator << (std::ostream& os, const Node& p) {
            os << p.id << " " << p.nchildren << " " << p.cid;
            return os;
        }
        friend std::istream& operator >> (std::istream& is, Node& p) {
            is >> p.id >> p.nchildren >> p.cid;
            return is;
        }
    };
    
    typedef std::vector<interBoundary> InterBoundVector;
    typedef std::vector<Patch> PatchVector;
    typedef std::vector<Node> NodeVector;
    
    /** Mesh object */
    struct MeshObject {
        
        Vertices mVertices; /**< vertices */
        Facets   mFacets;   /**< facets */
        Cells    mCells;    /**< Cells */
        
        std::string name;           /**< File name */
        Boundaries  mBoundaries;    /**< List of boundary patches */
        IntVector   mFOC;           /**< Facet owners of elements */
        IntVector   mFNC;           /**< Facet neighbors of elements */
        
        Int      mNV;   /**< Number of vertices */
        Int      mNF;   /**< Number of faces */
        Int      mBCS;  /**< Number of internal cells */
        Int      mBCSI; /**< Number of internal cells one layer away from boundary */

        PatchVector      mPatches;      /**< List of patches */
        InterBoundVector mInterMesh;    /**< List of inter-processor boundaries */
        
        Cells    mFaceID;   /**< Original face orientation*/

        VectorVector mFC;   /**< Facet centers */
        VectorVector mCC;   /**< Cell centers */
        VectorVector mFN;   /**< Facet normals */
        ScalarVector mCV;   /**< Cell volumes */
        BoolVector   mReversed; /**< Is normal reversed? */
        
        NodeVector   mAmrTree;  /**< AMR tree */
        
        /*functions*/
        void clear();
        void writeMesh(std::ostream&);
        bool readMesh(Int = 0,bool = true);
        void writeMshMesh(std::ostream&);
        void readMshMesh(std::istream&);
        void addBoundaryCells();
        void calcGeometry();
        void removeBoundary(const IntVector&);
        Int  removeUnusedVertices(Int = 0);
        void breakEdges(Int);
        
        void straightenEdges(const Facet&, Facet&, Facet&);
        bool coplanarFaces(const Facet&,const Facet&);
        bool mergeFacets(const Facet&,const Facet&, Facet&);
        void mergeFacetsCell(const Cell&,const IntVector&,Facet&);
        void mergeCells(Cell&,const Cell&,IntVector&);
        void addVerticesToEdge(const int, Facet&, const Facet&);
        void calcFaceCenter(const Facet&,Vector&);
        void calcCellCenter(const Cell&, Vector&);
        void calcUnitNormal(const Facet&,Vector&);
        void initFaceInfo(IntVector&,Cells&,const IntVector&,const Cells&,Int);
        void refineFacet(const Facet&, Facets&, Int, Int); 
        void refineFacets(const IntVector&, IntVector&, const IntVector&, IntVector&, IntVector&,Int);
        void refineCell(const Cell&, IntVector&, Int, IntVector&,
                        IntVector&, IntVector&, Cells&, IntVector&, Int);
        void refineMesh(const IntVector&, const IntVector&, const IntVector&, 
                        const IntVector&, IntVector&, IntVector&);
    };
    
    /** \name Global mesh object with its members */
    //@{
    extern  MeshObject        gMesh;
    extern  std::string&      gMeshName;
    extern  Vertices&         gVertices;
    extern  Facets&           gFacets;
    extern  Cells&            gCells;
    extern  Boundaries&       gBoundaries;
    extern  IntVector&        gFOC;
    extern  IntVector&        gFNC;
    extern  Int&              gBCS;
    extern  Int&              gBCSI;
    extern  Cells&            gFaceID;
    extern  InterBoundVector& gInterMesh;
    extern  NodeVector&       gAmrTree;
    extern  VectorVector&     gFC;
    extern  VectorVector&     gCC;
    extern  VectorVector&     gFN;
    extern  ScalarVector&     gCV;
    //@}
    
    /** probe points */
    extern  Vertices         probePoints;
    
    void enroll(Util::ParamList& params);
}

namespace Controls {
    extern RefineParams refine_params;
    extern DecomposeParams decompose_params;
    void enrollRefine(Util::ParamList& params);
    void enrollDecompose(Util::ParamList& params);
}
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
    
    extern  std::vector<BasicBCondition*> AllBConditions;
    /** Clear list of  BCs */
    inline void clearBC() {
        forEach(AllBConditions,i)
            delete AllBConditions[i];
        AllBConditions.clear();
    }
    /* point in line/polygon*/
    bool pointInLine(const Vector&,const Vector&,const Vector&);
    bool pointInPolygon(const VectorVector&,const IntVector&,const Vector&);
}
#endif
