#ifndef __MESH_H
#define __MESH_H

#include "util.h"

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
        Int id;           /**< Cell id in mCells */ 
        Int cid;          /**< First child id in mAmrTree */
        Int nchildren;    /**< Number of children */
        Int level;        /**< Refinement level */
        Node() {
            id = 0;
            cid = 0;
            nchildren = 0;
            level = 0;
        }
        template<typename Ts>
        friend Ts& operator << (Ts& os, const Node& p) {
            os << p.id << " " << p.nchildren << " " << p.cid << " " << p.level;
            return os;
        }
        template<typename Ts>
        friend Ts& operator >> (Ts& is, Node& p) {
            is >> p.id >> p.nchildren >> p.cid >> p.level;
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
        IntVector   mFOC;           /**< Facet owner cells */
        IntVector   mFNC;           /**< Facet neighbor cells */
        IntVector   mFMC;           /**< Facet mortar flags
                                        0=conforming,
                                        1=facet owner is mortar,
                                        2=facet neighbor is mortar*/
        BoolVector  mHasBoundary;   /**< True if a cell has boundary face */

        Int      mNV;   /**< Number of internal vertices */
        Int      mNF;   /**< Number of internal faces */
        Int      mBCS;  /**< Number of internal cells */

        PatchVector      mPatches;      /**< List of patches */
        InterBoundVector mInterMesh;    /**< List of inter-processor boundaries */

        Cells    mFaceID;   /**< Original face orientation*/

        VectorVector mFC;   /**< Facet centers */
        VectorVector mCC;   /**< Cell centers */
        VectorVector mFN;   /**< Facet normals */
        ScalarVector mCV;   /**< Cell volumes */
        BoolVector   mReversed; /**< Is normal reversed? */

        NodeVector   mAmrTree;  /**< AMR tree */

        Vertices mVerticesNoExtrude; /**< Vertices before extruded on sphere */

        /*functions*/
        void clear();
        void writeMshMesh(std::ostream&);
        void readMshMesh(std::istream&);
        void addBoundaryCells();
        void getHexCorners(const Facet&, const Facet&, IntVector& vp);
        void fixHexCells();
        void calcGeometry();
        void removeBoundary(const IntVector&);
        Int  removeUnusedVertices(Int = 0);
        void breakEdges(Int);
        void ExtrudeMesh();

        void straightenEdges(const Facet&, Facet&, Facet&);
        bool coplanarFaces(const Facet&,const Facet&);
        bool mergeFacets(const Facet&,const Facet&, Facet&);
        void mergeFacetsGroup(const IntVector&,Facet&,const Cell* = 0);
        void mergeCells(Cell&,const Cell&,IntVector&);
        void addVerticesToEdge(const int, Facet&, const Facet&);
        void calcFaceCenter(const Facet&,Vector&);
        void calcCellCenter(const Cell&, Vector&);
        void calcUnitNormal(const Facet&,Vector&);
        void initFaceInfo(IntVector&,Cells&,const IntVector&,const Cells&);
        void refineFacet(const Facet&, Facets&, Int, Int); 
        void refineFacets(const IntVector&, IntVector&, const IntVector&, IntVector&, IntVector&,Int);
        void refineCell(const Cell&, IntVector&, Int, IntVector&,
                IntVector&, IntVector&, Cells&, IntVector&, Int);
        void refineMesh(const IntVector&, const IntVector&, const IntVector&, 
                const IntVector&, IntVector&, IntVector&, IntVector&);

        template<typename Ts>
        friend Ts& operator << (Ts& os, const MeshObject& p) {
            p.writeTextMesh(os);
            return os;
        }
        template<typename Ts>
        friend Ts& operator >> (Ts& is, MeshObject& p) {
            p.readTextMesh(is);
            return is;
        }
        /**
          Read mesh from file in text format
         */
        template<typename Ts>
        bool readTextMesh(Ts& is) {
            using Util::operator>>;
            Int size;
            char symbol;
            /*clear*/
            clear();
            /*read*/
            is >> mVertices;
            is >> mFacets;
            is >> mCells;
            is >> size >> symbol;
            for(Int i = 0; i < size; i++) {
                IntVector index;
                std::string str;
                is >> str;
                is >> index;

                IntVector& gB = mBoundaries[str];
                gB.insert(gB.begin(),index.begin(),index.end());

                /*internal mesh boundaries*/
                if(str.find("interMesh") != std::string::npos) {
                    interBoundary b;
                    b.f = &mBoundaries[str];
                    std::replace(str.begin(), str.end(), '_', ' ');
                    std::stringstream ss(str);
                    ss >> str >> b.from >> b.to;
                    mInterMesh.push_back(b);
                }
            }
            is >> symbol;
            /*start of buffer*/ 
            Int buffer_index = 0;
            forEach(mInterMesh,i) {
                interBoundary& b = mInterMesh[i];
                b.buffer_index = buffer_index;
                buffer_index += b.f->size();
            }
            return true;
        }
        /**
          Write mesh to file in text format
         */
        template<typename Ts>
        void writeTextMesh(Ts& os) const {
            using Util::operator<<;
            os.precision(12);
            if(mVerticesNoExtrude.size())
                os << mVerticesNoExtrude;
            else
                os << mVertices;
            os.precision(6);
            os << mFacets;
            /*cells*/
            Int size = mBCS;
            os << size << "\n{\n";
            for(Int i = 0; i < size; i++)
                os << mCells[i] << "\n";
            os << "}\n";
            /*boundaries*/
            os << mBoundaries;
        }
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
    extern  IntVector&        gFMC;
    extern  BoolVector&       gHasBoundary;
    extern  Int&              gBCS;
    extern  Cells&            gFaceID;
    extern  InterBoundVector& gInterMesh;
    extern  NodeVector&       gAmrTree;
    extern  VectorVector&     gFC;
    extern  VectorVector&     gCC;
    extern  VectorVector&     gFN;
    extern  ScalarVector&     gCV;
    extern  Int               gNCells;
    extern  Int               gNFacets;
    extern  Int               gNVertices;
    //@}

    extern  Vector amr_direction; /**< Direction of AMR */
    extern  Int is_spherical; /**< Mesh is spherical */
    extern Scalar sphere_radius; /**< Radius of the sphere */
    extern Scalar sphere_height; /**< Height of the sphere */

    bool pointInLine(const Vector&,const Vector&,const Vector&);
    bool pointInPolygon(const VectorVector&,const IntVector&,const Vector&);
}

#endif
