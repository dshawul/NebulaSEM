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

        Int      mNV;   /**< Number of internal vertices */
        Int      mNF;   /**< Number of internal faces */
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

    extern  Vector amr_direction; /**< Direction of AMR */

    bool pointInLine(const Vector&,const Vector&,const Vector&);
    bool pointInPolygon(const VectorVector&,const IntVector&,const Vector&);
}

#endif
