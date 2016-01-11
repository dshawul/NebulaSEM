#include "mesh.h"

using namespace std;

// #define RDEBUG

/**
References to global mesh
*/
namespace Mesh {
    MeshObject        gMesh;
    std::string&      gMeshName = gMesh.name;
    Vertices&         gVertices = gMesh.mVertices;
    Facets&           gFacets   = gMesh.mFacets;
    Cells&            gCells    = gMesh.mCells;
    Boundaries&       gBoundaries = gMesh.mBoundaries;
    IntVector&        gFOC = gMesh.mFOC;
    IntVector&        gFNC = gMesh.mFNC;
    Int&              gBCS = gMesh.mBCS;
    Int&              gBCSI = gMesh.mBCSI;
    Cells&            gFaceID = gMesh.mFaceID;
    InterBoundVector& gInterMesh = gMesh.mInterMesh;
    NodeVector&       gAmrTree = gMesh.mAmrTree;
    VectorVector&     gFC = gMesh.mFC;
    VectorVector&     gCC = gMesh.mCC;
    VectorVector&     gFN = gMesh.mFN;
    ScalarVector&     gCV = gMesh.mCV;
    
    vector<BasicBCondition*> AllBConditions;
    Vertices         probePoints;
}
/**
Refinement and domain decomposition parameters
*/
namespace Controls {
    RefineParams refine_params;
    DecomposeParams decompose_params;
}
/**
Enroll refine parameters
*/
void Controls::enrollRefine(Util::ParamList& params) {
    params.enroll("direction",&refine_params.dir);
    params.enroll("field",&refine_params.field);
    params.enroll("field_max",&refine_params.field_max);
    params.enroll("field_min",&refine_params.field_min);
    params.enroll("limit",&refine_params.limit);
}
/**
Enroll domain decomposition parameters
*/
void Controls::enrollDecompose(Util::ParamList& params) {
    params.enroll("n",&decompose_params.n);
    params.enroll("axis",&decompose_params.axis);
    Util::Option* op = new Util::Option(&decompose_params.type, 4, 
            "XYZ","CELLID","METIS","NONE");
    params.enroll("type",op);
}
/**
Clear mesh
*/
void Mesh::clear() {
    gMesh.clear();
    Mesh::clearBC();
    Mesh::probePoints.clear();
}
/**
Clear mesh object
*/
void Mesh::MeshObject::clear() {
    mVertices.clear();
    mFacets.clear();
    mCells.clear();
    mBoundaries.clear();
    mFOC.clear();
    mFNC.clear();
    mInterMesh.clear();
    mPatches.clear();
    mAmrTree.clear();
}
/**
Read mesh at given time step
*/
bool Mesh::MeshObject::readMesh(Int step,bool first) {
    /*open file*/
    stringstream path;
    path << name << "_" << step;
    string str = path.str();
    ifstream is(str.c_str());
    if(is.fail()) {
        if(first) {
            str = name;
            is.open(str.c_str());
        } else
            return false;
    }
    /*read*/
    is >> hex;
    is >> mVertices;
    is >> mFacets;
    is >> mCells;
    while(Util::nextc(is)) {
        IntVector index;
        string str;
        is >> str;
        is >> index;

        IntVector& gB = mBoundaries[str];
        gB.insert(gB.begin(),index.begin(),index.end());

        /*internal mesh boundaries*/
        if(str.find("interMesh") != std::string::npos) {
            interBoundary b;
            sscanf(str.c_str(), "interMesh_%x_%x", &b.from,&b.to);
            b.f    = &mBoundaries[str];
            mInterMesh.push_back(b);
        }
    }
    /*start of buffer*/ 
    Int buffer_index = 0;
    forEach(mInterMesh,i) {
        interBoundary& b = mInterMesh[i];
        b.buffer_index = buffer_index;
        buffer_index += b.f->size();
    }
    is >> dec;
    return true;
}
/**
Write mesh to file
*/
void Mesh::MeshObject::writeMesh(ostream& os) {
    os << hex;
    os.precision(12);
    os << mVertices;
    os.precision(6);
    os << mFacets;
    os << mCells;
    forEachIt(Boundaries,mBoundaries,it)
        os << it->first << " " << it->second << endl;
    os << dec;
}
/**
Add boundary cells around mesh
*/
void Mesh::MeshObject::addBoundaryCells() {
    using namespace Constants;
    
    /*neighbor and owner cells of face*/
    mBCS = mCells.size();
    mFOC.assign(mFacets.size(),MAX_INT);
    mFNC.assign(mFacets.size(),MAX_INT);
    forEach(mCells,i) {
        Cell& c = mCells[i];
        forEach(c,j) {
            Int fi = c[j];
            if(mFOC[fi] == MAX_INT) 
                mFOC[fi] = i;
            else 
                mFNC[fi] = i;
        }
    }
    /*Flag boundary faces not in mBoundaries for auto deletion*/
    {
        IntVector faceInB;
        faceInB.assign(mFacets.size(),0);
        forEachIt(Boundaries,mBoundaries,it) {
            IntVector& mB = it->second; 
            forEach(mB,j)
                faceInB[mB[j]] = 1;
        }
    
        IntVector& mDelete = mBoundaries["delete"];
        forEach(mFNC,i) {
            if(mFNC[i] == MAX_INT) {
                if(!faceInB[i])
                    mDelete.push_back(i);
            }
        }
    }
    /*reorder cells*/
    {   
        Cells bcs;
        IntVector allbs;
        Int count = 0;
        Int bdry_size = 0;
        
        forEachIt(Boundaries,mBoundaries,it)
            bdry_size += it->second.size();
        bcs.resize(bdry_size);
        allbs.resize(bdry_size);
        forEach(mCells,i) {
            Cell& c = mCells[i];
            forEach(c,j) {
                Int fi = c[j];
                if(mFNC[fi] == MAX_INT) {
                    allbs[count] = i;
                    bcs[count] = c;
                    count++;
                    break;
                }
            }
        }
        mBCSI = mBCS - count;
        allbs.resize(count);
        bcs.resize(count);
        erase_indices(mCells,allbs);
        mCells.insert(mCells.end(),bcs.begin(),bcs.end());
    }
    /*add boundary cells*/
    forEachIt(Boundaries,mBoundaries,it) {
        IntVector& mB = it->second;
        forEach(mB,j) {
            Int fi = mB[j];
            if(mFNC[fi] == MAX_INT) {
                Cell c;
                c.push_back(fi);
                mCells.push_back(c);
                mFNC[fi] = mCells.size() - 1;
            }
        }
    }
}
/**
Calculate geometric information
*/
void Mesh::MeshObject::calcGeometry() {
    Int i;
    /*allocate*/
    mFC.assign(mFacets.size(),Vector(0));
    mCC.assign(mCells.size(),Vector(0));
    mFN.assign(mFacets.size(),Vector(0));
    mCV.assign(mCells.size(),Scalar(0));
    mReversed.assign(mFacets.size(),false);

    /* face centre*/
    forEach(mFacets,i) {
        Facet& f = mFacets[i];
        Vector C(0);
        forEach(f,j)
            C += mVertices[f[j]];
        mFC[i] = C / Scalar(f.size());
    }
    /* cell centre */
    forEach(mCells,i) {
        Cell& c = mCells[i];
        Vector C(0);
        forEach(c,j)
            C += mFC[c[j]];
        mCC[i] = C / Scalar(c.size());
    }
    /* face normal */
    forEach(mFacets,i) {
        Facet& f = mFacets[i];
        Vector N(0),C(0),Ci,Ni;
        Scalar Ntot = Scalar(0);
        const Vector& v1 = mFC[i];
        forEach(f,j) {
            const Vector& v2 = mVertices[f[j]];
            const Vector& v3 = mVertices[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]];

            Ni = ((v2 - v1) ^ (v3 - v1));
            Scalar magN = mag(Ni);
            Ci = magN * ((v1 + v2 + v3) / 3);
            
            
            C += Ci;
            Ntot += magN;
            N += Ni;
        }
        mFC[i] = C / Ntot;    /*corrected face centre*/
        Vector v = mFC[i] - mCC[mFOC[i]];
        if((v & N) < 0) {
            N = -N;
            mReversed[i] = true;
        }
        mFN[i] = N / Scalar(2);
    }
    /* cell volumes */
    for(i = 0;i < mBCS;i++) {
        Cell& c = mCells[i];
        Scalar V(0);
        Vector C(0);
        forEach(c,j) {
            Vector v = mCC[i] - mFC[c[j]];
            Scalar Vi = mag(v & mFN[c[j]]);
            C += Vi * (3 * mFC[c[j]] + mCC[i]) / 4;
            V += Vi;
        }
        mCC[i] = C / V;         /*corrected cell centre */
        mCV[i] = V / Scalar(3);
    }
    /*boundary cell centre and volume*/
    forEachS(mCells,i,mBCS) {
        Int fi = mCells[i][0];
        mCV[i] = mCV[mFOC[fi]];
        mCC[i] = mFC[fi];
    }
    /*facet ids*/
    mFaceID.clear();
    forEach(mCells,i) {
        IntVector b;
        forEach(mCells[i],j)
            b.push_back(j);
        mFaceID.push_back(b);
    }
}
/** 
Remove empty boundary (for 2D domain)
*/
void Mesh::MeshObject::removeBoundary(const IntVector& fs) {
    Int count;
    IntVector Idf(mFacets.size(),0);
    IntVector Idc(mCells.size(),0);
    
    /*erase facet reference*/
    forEach(fs,i) {
        Int f = fs[i];
        Cell& co = mCells[mFOC[f]];
        Cell& coid = mFaceID[mFOC[f]];
        forEach(co,j) {
            if(co[j] == f) {
                co.erase(co.begin() + j); 
                coid.erase(coid.begin() + j);
                break; 
            }
        }
        Cell& cn = mCells[mFNC[f]];
        Cell& cnid = mFaceID[mFNC[f]];
        forEach(cn,j) {
            if(cn[j] == f) { 
                cn.erase(cn.begin() + j); 
                cnid.erase(cnid.begin() + j); 
                break; 
            }
        }
    }
    
    /*updated facet id*/
    forEach(fs,i)
        Idf[fs[i]] = Constants::MAX_INT;
    count = 0;
    forEach(mFacets,i) {
        if(Idf[i] != Constants::MAX_INT) 
            Idf[i] = count++;
        else
            mFacets[i].clear();
    }

    /*erase facets*/
    IntVector fzeroIndices;
    forEach(mFacets,i) {
        if(mFacets[i].size() == 0)
            fzeroIndices.push_back(i);
    }
    erase_indices(mFacets,fzeroIndices);
    erase_indices(mFOC,fzeroIndices);
    erase_indices(mFNC,fzeroIndices);
    erase_indices(mFC,fzeroIndices);
    erase_indices(mFN,fzeroIndices);
    /*updated cell id*/
    count = 0;
    forEach(mCells,i) {
        if(mCells[i].size() != 0) 
            Idc[i] = count++;
        else
            Idc[i] = Constants::MAX_INT;
    }
    /*erase cells*/
    IntVector czeroIndices;
    forEach(mCells,i) {
        if(mCells[i].size() == 0)
            czeroIndices.push_back(i);
    }
    erase_indices(mCells,czeroIndices);
    erase_indices(mFaceID,czeroIndices);
    erase_indices(mCC,czeroIndices);
    erase_indices(mCV,czeroIndices);
    
    /*updated facet id*/
    forEach(mCells,i) {
        forEach(mCells[i],j) {
            mCells[i][j] = Idf[mCells[i][j]];
        }
    }
    /*facet owner and neighbor*/
    forEach(mFacets,i) {
        mFOC[i] = Idc[mFOC[i]];
        mFNC[i] = Idc[mFNC[i]];
    }
    /*patches*/
    forEachIt(Boundaries,mBoundaries,it) {
        IntVector& mB = it->second;
        forEach(mB,i)
            mB[i] = Idf[mB[i]];
    }
}
/** 
Remove unused vertices 
*/
Int Mesh::MeshObject::removeUnusedVertices(Int ivBegin) {
    if(!ivBegin) 
        ivBegin = mVertices.size();
    
    IntVector isUsed(mVertices.size(),0);
    forEach(mFacets,i) {
        Facet& f = mFacets[i];
        forEach(f,j)
            isUsed[f[j]] = 1;
    }
    IntVector rVertices;
    Int cnt = 0;
    Int ivCount = 0;
    forEach(isUsed,i) {
        if(!isUsed[i]) {
            rVertices.push_back(i);
            if(i < ivBegin) ivCount++;
        } else
            isUsed[i] = cnt++;
    }
    if(rVertices.size()) {
        forEach(mFacets,i) {
            Facet& f = mFacets[i];
            forEach(f,j)
                f[j] = isUsed[f[j]];
        }
        erase_indices(mVertices,rVertices);
    }
    return ivCount;
}
/** 
Break edges of faces that are not set for refinement 
but has a neighboring refined face
*/
void Mesh::MeshObject::breakEdges(Int ivBegin) {
    forEach(mFacets,i) {
        Facet& f = mFacets[i];
        Facet nf;
        forEach(f,j) {
            nf.push_back(f[j]);
            
            const Vector& v1 = mVertices[f[j]];
            const Vector& v2 = mVertices[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]];

            vector< pair<Int,Scalar> > order;
            for(Int k = ivBegin; k < mVertices.size();k++) {
                const Vector& v = mVertices[k];
                if(pointInLine(v,v1,v2)) {
                    Scalar val = dot((v - v1),(v2 - v1));
                    pair<Int,Scalar> p;
                    p.first = k;
                    p.second = val;
                    order.push_back(p);
                }
            }
            
            if(order.size()) {
                sort(order.begin(), order.end(), compare_pair_second<>());
                
                forEach(order,i)
                    nf.push_back(order[i].first);
            }
        }
        f = nf;
    }
}
/**
Is point inside line segment ?
*/
bool Mesh::pointInLine(const Vector& v,const Vector& v1,const Vector& v2) {
    using namespace Constants;
    
    Vector p = v - v1;
    Vector q = v - v2;
    Scalar e;
    e = p[YY] * q[ZZ] - p[ZZ] * q[YY];
    if(!equal(e,Scalar(0))) return false;
    e = p[ZZ] * q[XX] - p[XX] * q[ZZ];
    if(!equal(e,Scalar(0))) return false;
    e = p[XX] * q[YY] - p[YY] * q[XX];
    if(!equal(e,Scalar(0))) return false;
    
    e = dot((v - v2),(v1 - v2));
    if(e > Scalar(0)) {
        Scalar e1 = dot((v1 - v2),(v1 - v2));
        if(e < e1)
            return true;
    }

    return false;
}
/**
Is point inside convex polygon (face) ?
*/
bool Mesh::pointInPolygon(const VectorVector& points,const IntVector& f,const Vector& C) {
   Scalar sum = 0;

   forEach(f,j) {
      Vector p1 = points[f[j]] - C;
      Vector p2 = points[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]] - C;

      Scalar mg = mag(p1) * mag(p2);
      if (equal(mg,Scalar(0)))
         return true;
      sum += acos(dot(p1,p2) / mg);
   }
   
   if(equal(sum, 2 * Constants::PI))
       return true;
   
   return false;
}
/**
Calculate unit-normal of face
*/
void Mesh::MeshObject::calcUnitNormal(const Facet& f,Vector& N) {    
    const Vector& v1 = mVertices[f[0]];
    for(Int j = 1;j < f.size();j++) {
        const Vector& v2 = mVertices[f[j]];
        const Vector& v3 = mVertices[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]];
        
        if(!Mesh::pointInLine(v2,v1,v3)) {
            N = (v2 - v1) ^ (v3 - v1);
            N = unit(N);
            return;
        }
    }
}
/**
Straighten edges of faces
*/
void Mesh::MeshObject::straightenEdges(const Facet& f, Facet& r, Facet& removed) {    
    Vector v1 = mVertices[f[f.size() - 1]];
    for(Int j = 0;j < f.size();j++) {
        const Vector& v2 = mVertices[f[j]];
        const Vector& v3 = mVertices[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]];
        
        if(!Mesh::pointInLine(v2,v1,v3)) {
            r.push_back(f[j]);
            v1 = v2;
        } else
            removed.push_back(f[j]);
    }
}
/**
Check if faces are co-planar
*/
bool Mesh::MeshObject::coplanarFaces(const Facet& f1,const Facet& f2) {
    Vector N1,N2;
    calcUnitNormal(f1,N1);
    calcUnitNormal(f2,N2);
    Scalar e = magSq(N1 ^ N2);
    if(equal(e,Scalar(0))) {
        Vector v = mVertices[f2[0]] - mVertices[f1[0]];
        e = dot(N1,v);
        if(equal(e,Scalar(0)))
            return true;
    }
    return false;
}
/**
Union of non-overlaping polygons
*/
bool Mesh::MeshObject::mergeFacets(const Facet& f1_,const Facet& f2_, Facet& f) {
    Vector N1,N2;
    calcUnitNormal(f1_,N1);
    calcUnitNormal(f2_,N2);
    
    //make them same direction
    Facet f1 = f1_;
    Facet f2 = f2_;
    if(dot(N1,N2) < 0) {
        forEach(f2_,j)
            f2[f2.size() - j - 1] = f2_[j];
    }
    
    //rotate nodes of faces to first non-shared node
    Int contained = false;
    {
        Int v1 = -1,v2 = -1;
        forEach(f1,i) {
            Int j;
            for(j = 0;j < f2.size();j++) {
                if(f1[i] == f2[j]) 
                    break;
            }
            if(j == f2.size()) {
                v1 = i;
                break;
            }
        }
        forEach(f2,i) {
            Int j;
            for(j = 0;j < f1.size();j++) {
                if(f2[i] == f1[j]) 
                    break;
            }
            if(j == f1.size()) {
                v2 = i;
                break;
            }
        }
        std::rotate(f1.begin(),f1.begin() + v1,f1.end());
        if(v2 != -1)
            std::rotate(f2.begin(),f2.begin() + v2,f2.end());
        else
            contained = true;
    }   
    
    //find start and end of shared edges
    int a[2];
    int b[2];
    Int count;
    
    count = 0;
    forEach(f1,i) {
        forEach(f2,j) {
            if(f1[i] == f2[j]) {
                if(!count) {
                    a[0] = i;
                    b[0] = j;
                } else {
                    a[1] = i;
                    b[1] = j;
                }
                count++;
            }
        }
    }
    
    //should share atleast one  edge
    if(count < 2) 
        return false;
    
    //merge
    f.clear();
    for(Int i = 0;i <= a[0];i++)
        f.push_back(f1[i]);
    if(contained) {
        for(Int i = b[0] + 1;i < b[1];i++)
            f.push_back(f2[i]);
    } else {
        for(Int i = b[0] + 1;i < f2.size();i++)
            f.push_back(f2[i]);
        for(Int i = 0;i < b[1];i++)
            f.push_back(f2[i]);
    }
    for(Int i = a[1];i < f1.size();i++)
        f.push_back(f1[i]);
    
    return true;
}
/**
Merge facets of cell
*/
void Mesh::MeshObject::mergeFacetsCell(const Cell& c1,const IntVector& shared1,Facet& f) {
    IntVector flag;
    flag.assign(shared1.size(),0);
    bool has,mhas;
    f = mFacets[c1[shared1[0]]];
    do {
        has = false;
        mhas = false;
        forEachS(shared1,p,1) {
            if(flag[p]) continue;
            const Facet& f1 = mFacets[c1[shared1[p]]];
            Facet f2 = f;
            if(mergeFacets(f2,f1,f)) {
                flag[p] = 1;
                mhas = true;
            } else {
                has = true;
            }
        }
    } while(has && mhas);

#ifdef RDEBUG
    if(f.size() < 4)  {
        std::cout << "Merge failed.\n";
        std::cout << "=======================\n";
        std::cout << c1 << std::endl;
        std::cout << shared1 << std::endl;
        forEach(shared1,i)
            std::cout << gFacets[c1[shared1[i]]] << std::endl;
        std::cout << "-----------------------\n";
        {
            Facet r,rr;
            straightenEdges(f,r,rr);
            std::cout << "Merged face  : " << f << std::endl;
            std::cout << "Straight Edge: " << r << std::endl;
            std::cout << "Removed      : " << rr << std::endl;
        }
        forEach(f,i)
            std::cout << mVertices[f[i]] << " : ";
        std::cout << std::endl;
        exit(0);
    }
#endif
}
/**
Merge two cells
*/
void Mesh::MeshObject::mergeCells(Cell& c1, const Cell& c2, IntVector& delFacets) {
    forEach(c2,m) {
        Int c2m = c2[m];
        Int n;
        for(n = 0;n < c1.size();n++) {
            if(c2m == c1[n])
                break;
        }
        if(n == c1.size())
            c1.push_back(c2m);
        else {
            delFacets.push_back(c1[n]);
            c1.erase(c1.begin() + n);
        }
    }
}
/**
Add vertices to edges of a face
*/
void Mesh::MeshObject::addVerticesToEdge(const int va, Facet& f, const Facet& fp) {
    const Vector& v1 = mVertices[f[f.size() - 1]];
    const Vector& v2 = mVertices[va];
    forEach(fp,i) {
        const Vector& v = mVertices[fp[i]];
        if(pointInLine(v,v1,v2))
            f.push_back(fp[i]);
    }
}
/**
Calculate face center
*/
void Mesh::MeshObject::calcFaceCenter(const Facet& f,Vector& fCj) {
    Vector C(0);
    forEach(f,j)
        C += mVertices[f[j]];
    fCj = C / Scalar(f.size());

    Vector Ctot(0);
    Scalar Ntot = Scalar(0);
    
    const Vector& v1 = fCj;
    forEach(f,j) {
        const Vector& v2 = mVertices[f[j]];
        const Vector& v3 = mVertices[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]];
        
        Vector Ni = ((v2 - v1) ^ (v3 - v1));
        Scalar magN = mag(Ni);
        Vector Ci = magN * ((v1 + v2 + v3) / 3);

        Ctot += Ci;
        Ntot += magN;
    }
    fCj = Ctot / Ntot;    /*corrected face centre*/
}
/**
Calculate cell center
*/
void Mesh::MeshObject::calcCellCenter(const Cell& c, Vector& cCj) {
    
    //approximate cell centre
    Int cnt = 0;
    cCj = Vector(0);
    forEach(c,i) {
        const Facet& f = mFacets[c[i]];
        forEach(f,j) {
            cCj += mVertices[f[j]];
            cnt++;
        }
    }
    cCj /= cnt;
    
    //exact cell centre
    Scalar Vt(0);
    Vector Ct(0);
    forEach(c,i) {
        const Facet& f = mFacets[c[i]];
    
        Vector fCj,fNj,C;
    
        C = Vector(0);
        forEach(f,j)
            C += mVertices[f[j]];
        fCj = C / Scalar(f.size());
    
        Vector N(0);
        Scalar Ntot = Scalar(0);
    
        C = Vector(0);
        const Vector& v1 = fCj;
        forEach(f,j) {
            const Vector& v2 = mVertices[f[j]];
            const Vector& v3 = mVertices[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]];
        
            Vector Ni = ((v2 - v1) ^ (v3 - v1));
            Scalar magN = mag(Ni);
            Vector Ci = magN * ((v1 + v2 + v3) / 3);
        
            C += Ci;
            Ntot += magN;
            N += Ni;
        }
        fCj = C / Ntot;
        fNj = N / Scalar(2);
    
        Scalar Vi = mag((cCj - fCj) & fNj);
        Ct += Vi * (3 * fCj + cCj) / 4;
        Vt += Vi;
    }

    cCj = Ct / Vt;
}
/**
Refine facet
*/
void Mesh::MeshObject::refineFacet(const Facet& f_, Facets& newf, Int dir, Int ivBegin) {
    using namespace Controls;
    
    const bool refine3D = equal(Scalar(0.0),mag(refine_params.dir));
    
    /*add face center*/
    Vector C;
    calcFaceCenter(f_,C);

    mVertices.push_back(C);
    Int fci = mVertices.size() - 1;

    /*straighten edges of face*/
    Facet f,fr;
    straightenEdges(f_,f,fr);

    IntVector midpts;
    midpts.resize(f.size(),Constants::MAX_INT);

    Int fmid = Constants::MAX_INT;

    forEach(f,j) {
        const Vector& v1 = mVertices[f[j]];
        const Vector& v2 = mVertices[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]];

        /*anisotropic refinement*/
#define RDIR(mDir) {                        \
    Scalar angle = acos(dot(mDir,uLine));   \
    if(angle < Constants::PI / 4 ||         \
       angle >= 3 * Constants::PI / 4)      \
            continue;                       \
}
       
        Vector uLine = unit(v2 - v1);
        
        if(!refine3D) {
            const Vector uDir = unit(refine_params.dir);
            RDIR(uDir);
        }
        
        if(!(dir & 1)) {
            const Vector xDir(1,0,0);
            RDIR(xDir);
        }
        if(!(dir & 2)) {
            const Vector yDir(0,1,0);
            RDIR(yDir);
        }
        if(!(dir & 4)) {
            const Vector zDir(0,0,1);
            RDIR(zDir);
        }
            
#undef RDIR
            
        /*save first mid-point location*/
        if(j < fmid) fmid = j;
        Vector Ce = (v1 + v2) / 2.0;

        /*add midpoint*/
        {
            Int k = Constants::MAX_INT;

            if(fr.size()) {
                forEach(fr,i) {
                    Int m = fr[i];
                    if(equal(Ce,mVertices[m])) {
                        midpts[j] = m;
                        k = m;
                        break;
                    }
                }
            }

            if(k == Constants::MAX_INT) {
                k = ivBegin;
                for(; k < mVertices.size();k++) {
                    if(equal(Ce,mVertices[k])) {
                        midpts[j] = k;
                        break;
                    }
                }
            }

            if(k == mVertices.size()) {
                mVertices.push_back(Ce);
                midpts[j] = k;
            }
        }
    }

    if(fmid != Constants::MAX_INT) {

#define ADD(x) {                    \
    if(fr.size())                   \
        addVerticesToEdge(x,fn,fr); \
    fn.push_back(x);                \
}

#define TEST() {                    \
    ADD(f[k]);                      \
    Int k1 = midpts[k];             \
    if(k1 == Constants::MAX_INT)    \
        skipped = true;             \
    else {                          \
        ADD(k1);                    \
        break;                      \
    }                               \
}
        /*Add faces*/
        forEachS(f,j,fmid) {
            Facet fn;
            fn.push_back(midpts[j]);

            bool skipped = false;

            Int k = j + 1;
            for(;k < f.size();k++)
                TEST();
            if(k == f.size()) {
                for(k = 0;k <= fmid;k++)
                    TEST();
                j = f.size();
            } else {
                j = k - 1;
            }

            if(!skipped)
                fn.push_back(fci);

            newf.push_back(fn);
        }
        
#undef TEST
#undef ADD
        
    } else {
        newf.push_back(f_);
    }

}
/**
Refine list of facets
*/
void Mesh::MeshObject::refineFacets(const IntVector& rFacets,IntVector& refineF, const IntVector& rfDirs,
                             IntVector& startF, IntVector& endF, Int ivBegin) {

    Facets newf;
    forEach(rFacets,i) {
        Int fi = rFacets[i];
        Int dir = rfDirs[i];
        const Facet& f = mFacets[fi];
        
        startF[fi] = mFacets.size() + newf.size();
        
        refineFacet(f,newf,dir,ivBegin);
        
        endF[fi] = mFacets.size() + newf.size();
    }
    
    mFacets.insert(mFacets.end(),newf.begin(),newf.end());
    startF.resize(mFacets.size(),0);
    endF.resize(mFacets.size(),0);
    refineF.resize(mFacets.size(),0);
}
/**
Refine a given cell
*/
void Mesh::MeshObject::refineCell(const Cell& c,IntVector& cr, Int rDir,
                IntVector& refineF,IntVector& startF,IntVector& endF,
                Cells& newc,IntVector& delFacets, Int ivBegin
                ) {
   /* *************************
    *
    * Form new cells (pyramids) 
    *
    * *************************/
    //co-planar faces
    typedef std::map<Int, std::pair<IntVector,IntVector> > typeMap;
    typeMap sharedMap;

    IntVector c1;
    IntVector c1r;
    
    c1.reserve(c.size());
    c1r.reserve(c.size());

    forEach(cr,j) {
        Int crj = cr[j];
        if(crj < Constants::MAX_INT / 2) {
            if(crj != rDir) {
                IntVector& shared = sharedMap[j].first;
                shared.push_back(j);
            } else {
                c1.push_back(c[j]);
                c1r.push_back(0);
            }
            continue;
        }
        crj -= Constants::MAX_INT / 2;

        IntVector& shared = sharedMap[crj].first;
        shared.push_back(j);
    }

    forEachIt(typeMap,sharedMap,it) {
        IntVector& shared = (it->second).first;
        
        Facet f;
        if(shared.size() > 1) {
            mergeFacetsCell(c,shared,f);
        } else {
            f = mFacets[c[shared[0]]];
        }
        
        Facets newf;
        refineFacet(f,newf,rDir,ivBegin);

        forEach(newf,j) {
            Int fi = mFacets.size();
            {
                mFacets.push_back(newf[j]);
                startF.push_back(0);
                endF.push_back(0);
                refineF.push_back(0);
            }
            (it->second).second.push_back(fi);
            c1.push_back(fi);
            c1r.push_back(shared[0] + Constants::MAX_INT / 2);
            delFacets.push_back(fi);
        }
        
#ifdef RDEBUG
        std::cout << it->first << " = " << it->second.first << " = " << it->second.second <<  std::endl;
        std::cout << "---\n";
#endif
        
    }

    //refine cell
    Vector C;
    calcCellCenter(c,C);
    mVertices.push_back(C);
    Int cci = mVertices.size() - 1;
    Int ifBegin = mFacets.size();

#ifdef RDEBUG
    std::cout << "Cell center: " << C << std::endl;
#endif

    IntVector ownerf;
    forEach(c1,j) {
        Int fi = c1[j];
        Int crj = c1r[j];

        /*cell is refined but face is not?*/
        IntVector list;
        if(crj) {
            list.push_back(fi);
            ownerf.push_back(crj);
        } else {
            for(Int k = startF[fi]; k < endF[fi];k++) {
                list.push_back(k);
                ownerf.push_back(fi);
            }
        }

#ifdef RDEBUG
        if(j == 0) std::cout << "====================================\n";
        std::cout << crj << " " << rDir << " ( " << startF[fi] << "  " << endF[fi] << " ) " << std::endl;
#endif
            
        /*refine cells*/
        forEach(list,k) {
            Int fni = list[k];
            Facet f = mFacets[fni];
            
            //straighten edges
            if(f.size() > 4)
            {
                Facet r,rr;
                straightenEdges(f,r,rr);
                f = r;
            }
            
#ifdef RDEBUG
            Vector fc;
            calcFaceCenter(f,fc);
            std::cout << f << std::endl;
            std::cout << fc << std::endl;
            std::cout << "----\n";
#endif
            
            Cell cn;
            cn.push_back(fni);

            forEach(f,l) {
                Int v1i = f[l];
                Int v2i = f[(l + 1 == f.size()) ? 0 : (l + 1)];

                //triangular face
                Facet fn;
                fn.push_back(v1i);
                fn.push_back(v2i);
                fn.push_back(cci);
                //duplicate
                Int k = ifBegin;
                for(; k < mFacets.size();k++) {
                    if(equal(fn,mFacets[k])) {
                        cn.push_back(k);
                        break;
                    }
                }
                if(k == mFacets.size()) {
                    mFacets.push_back(fn);
                    startF.push_back(0);
                    endF.push_back(0);
                    refineF.push_back(0);
                    cn.push_back(k);
                }
            }
            
            newc.push_back(cn);
        }
    }
    
    //Replace co-planar faces
    forEach(newc,i) {
        if(ownerf[i] < Constants::MAX_INT / 2)
            continue;
        
        Cell& nc = newc[i];
        Int fi = nc[0];
        Facet& nf = mFacets[fi];
        
        IntVector forg;
        std::pair<IntVector,IntVector> pair;
        pair = sharedMap[ownerf[i] - Constants::MAX_INT / 2];
        const IntVector& faces1 = pair.first;
        const IntVector& faces2 = pair.second;
        forEach(faces2,j) {
            if(fi == faces2[j]) {
                IntVector org;
                if(faces1.size() > 1) {
                    forEach(faces1, i)
                        org.push_back(c[faces1[i]]);
                } else {
                    Int m = c[faces1[0]];
                    for(Int k = startF[m]; k < endF[m];k++)
                        org.push_back(k);
                }
                
                forEach(org,k) {
                    Int m = org[k];
                    Vector fc;
                    calcFaceCenter(mFacets[m],fc);
                    if(Mesh::pointInPolygon(mVertices,nf,fc))
                        forg.push_back(m);
                }
                
                break;
            }
        }
            
        if(forg.size()) {
            nc.insert(nc.begin() + 1,forg.size() - 1,0);
            Int cnt = 0;
            forEach(forg,j)
                nc[cnt++] = forg[j];
        }
    }
    
    /* *********************
     *
     * Merge new cells 
     *
     * *********************/
    {
        /* Find among new cells on different owner faces,
         * and merge if they have a shared face
         */
#ifdef RDEBUG
        std::cout << "============\n";
        std::cout << " New cells before merge " << newc.size() << std::endl;
        std::cout << newc << std::endl;
        std::cout << "============\n";
#endif
        
        Cells mergec;
        forEach(newc,j) {
            Int o1 = ownerf[j];
            Cell& c1 = newc[j];
            IntVector mg;
            forEachS(newc,k,j+1) {
                Int o2 = ownerf[k];
                Cell& c2 = newc[k];

                if(o1 == o2) continue;
            
                /*Compare faces*/
                forEach(c1,m) {
                    Int f1 = c1[m];
                    forEach(c2,n) {
                        Int f2 = c2[n];
                        if(f1 == f2) {
                            mg.push_back(k);
                            goto END;
                        }
                    }
                }
                END:;
            }
            mergec.push_back(mg);
        }

        bool inserted;
        do {
            inserted = false;
            forEach(mergec,j) {
                IntVector& cm = mergec[j];
                forEach(cm,k) {
                    forEachS(mergec,m,j+1) {
                        IntVector& cn = mergec[m];
                        bool has = false;
                        if(m == cm[k] && cn.size())
                            has = true;
                        else {
                            forEach(cn,z) {
                                if(cn[z] == cm[k]) {
                                    has = true;
                                    break;
                                }
                            }
                        }
                        if(has) {
                            inserted = true;
                            cm.insert(cm.end(),cn.begin(),cn.end());
                            cm.push_back(m);
                            cn.clear();
                        }
                    }
                }
            }
        } while(inserted);
        
        /*merge cells*/
        {
            IntVector erasei;
            erasei.assign(newc.size(),0);
            forEachR(mergec,j) {
                IntVector& cm = mergec[j];
                Cell& c1 = newc[j];
                forEach(cm,k) {
                    if(erasei[cm[k]]) continue;
                    erasei[cm[k]] = 1;
                    Cell& c2 = newc[cm[k]];
        
                    mergeCells(c1,c2,delFacets);
                }
            }
            /*erase cells*/
            IntVector delCells;
            forEach(erasei,j) {
                if(erasei[j])
                    delCells.push_back(j);
            }
            erase_indices(newc,delCells);
        }
        
#ifdef RDEBUG
        std::cout << " New cells " << newc.size() << std::endl;
        std::cout << newc << std::endl;
        std::cout << "============\n";
#endif
    
        /*two cells should share only one face*/
        forEach(newc,j) {
            Cell& c1 = newc[j];
            forEachS(newc,k,j+1) {
                Cell& c2 = newc[k];
                //find shared faces
                IntVector shared1,shared2;
                forEach(c1,m) {
                    Int f1 = c1[m];
                    forEach(c2,n) {
                        Int f2 = c2[n];
                        if(f1 == f2) {
                            shared1.push_back(m);
                            shared2.push_back(n);
                        }
                    }
                }
                //more than one face shared between two cells
                if(shared1.size() > 1) {
                    Facet nf;
                    mergeFacetsCell(c1,shared1,nf);
                    
                    //straighten edges
                    {
                        Facet r,rr;
                        straightenEdges(nf,r,rr);
                        nf = r;
                    }
                    //merge into first face
                    Int fi = c1[shared1[0]];
                    mFacets[fi] = nf;

                    //remove merged facets
                    forEachS(shared1,p,1)
                        delFacets.push_back(c1[shared1[p]]);
                    
                    shared1.erase(shared1.begin());
                    shared2.erase(shared2.begin());
                    sort(shared2.begin(),shared2.end());
                    erase_indices(c1,shared1);
                    erase_indices(c2,shared2);
                }
            }
        }
    }
#ifdef RDEBUG
    std::cout << " After merge new cells " << newc.size() << std::endl;
    std::cout << newc << std::endl;
    forEach(newc,i) {
        Cell& c = newc[i];
        forEach(c,j) {
            if(mFacets[c[j]].size() <= 3) {
                std::cout << "Cell still has triangular face.\n";
                exit(0);
            }
        }
    }
#endif

}
/**
Initialize facet refinement information
*/
void Mesh::MeshObject::initFaceInfo(IntVector& refineF,Cells& crefineF,
                            const IntVector& rCells,const Cells& newCells, Int ivBegin) {

    forEach(rCells,i) {
        const Cell& c = newCells[rCells[i]];
        Cell& cr = crefineF[i];
        
        IntVector flag;
        flag.assign(c.size(),0);

        forEach(c,j) {
            if(flag[j]) continue;

            IntVector shared;
            shared.push_back(j);

            bool straight = false;
            forEachS(c,k,j+1) {
                if(flag[k]) continue;

                if(coplanarFaces(mFacets[c[j]],mFacets[c[k]])) {
                    straight = true;
                    shared.push_back(k);
                }
            }
            if(!straight) {
                refineF[c[j]] = 1;
            } else {
                
                const Int id = Constants::MAX_INT / 2 + shared[0];
                forEach(shared,k) {
                    Int v = shared[k];
                    cr[v] = id;
                    flag[v] = 1;
                } 
            }
        }
    }
}
/**
Refine mesh
*/
void Mesh::MeshObject::refineMesh(const IntVector& rCells,const IntVector& cCells,const IntVector& rLevel,
                                  const IntVector& rDirs, IntVector& cellMap,IntVector& coarseMap) {
    using namespace Controls;
    
    Int ivBegin = mVertices.size();
    IntVector delFacets;
    IntVector delCells;
    
    /*************************
    *
    *        COARSENING
    *
    *************************/
    {
        /*coarsen mesh*/
        IntVector delTree;
        forEach(mAmrTree,i) {
            Node& n = mAmrTree[i];
            if(n.nchildren) {
                bool coarsen = true;
                for(Int j = 0;j < n.nchildren;j++) {
                    Node& cn = mAmrTree[n.cid + j];
                    if(cn.nchildren || !cCells[cn.id]) {
                        coarsen = false;
                        break;
                    }
                }
                if(coarsen) {
                    Int cnid = mAmrTree[n.cid].id;
                    Cell& c1 = mCells[cnid];
                    coarseMap.push_back(n.nchildren);
                    coarseMap.push_back(cnid);
                    for(Int j = 1;j < n.nchildren;j++) {
                        Node& cn = mAmrTree[n.cid + j];
                        Cell& c2 = mCells[cn.id];
                        mergeCells(c1,c2,delFacets);
                        delCells.push_back(cn.id);
                        coarseMap.push_back(cn.id);
                        /*update owner/neighbor info*/
                        forEach(c2,k) {
                            Int fi = c2[k];
                            if(mFOC[fi] == cn.id)
                                mFOC[fi] = cnid;
                            else
                                mFNC[fi] = cnid;
                        }
                    }
                    for(Int j = 0;j < n.nchildren;j++)
                        delTree.push_back(n.cid + j);
                    n.nchildren = 0;
                    n.id = cnid;
                    n.cid = 0;
                }
            }
        }
        /*merge facets*/
        Int myfBegin = mFacets.size();
        forEach(coarseMap,i) {
            Int nchildren = coarseMap[i];
            Int ci = coarseMap[i + 1];
            Cell& c1 = mCells[ci];

#ifdef RDEBUG
            std::cout << "Coarsening faces of cell " << ci << " cC " << mCC[ci]  << std::endl;
#endif

            typedef std::map<Int,IntVector> Pair;
            Pair sharedMap;

            forEach(c1,j) {
                Int f1 = c1[j];
                if(f1 >= myfBegin) continue;
                Int o1 = (mFOC[f1] == ci) ? mFNC[f1] : mFOC[f1];
                if(o1 >= mCells.size()) {
                    Vector u = unit(mFN[f1]);
                    o1 = Int(0.75 * Constants::MAX_INT) +
                         int(Scalar(1000.0) * u[2] +
                             Scalar(100.0) * u[1] +
                             Scalar(10.0) * u[0]);
                }
                IntVector& val = sharedMap[o1];
                val.push_back(j);
            }

#ifdef RDEBUG
            forEachIt(Pair,sharedMap,it)
                std::cout << it->first << " " << it->second << std::endl;
            std::cout << "------\n";
#endif

            IntVector allDel;
            forEachIt(Pair,sharedMap,it) {
                IntVector& faces = it->second;
                if(faces.size() > 1) {
                    //merge facets of cell
                    Facet nf;
                    mergeFacetsCell(c1,faces,nf);

                    //merge into first face
                    Int fi = c1[faces[0]];
                    mFacets[fi] = nf;

                    //erase old faces
                    forEachS(faces,j,1)
                        delFacets.push_back(c1[faces[j]]);
                    allDel.insert(allDel.end(),faces.begin() + 1,faces.end());

                    if(it->first < Constants::MAX_INT / 2) {
                        Cell& c2 = mCells[it->first];
                        forEachS(faces,j,1)
                            eraseValue(c2,c1[faces[j]]);
                    } else {
                        //Assume all faces are on same patch [Fix me later]
                        forEachIt(Boundaries,mBoundaries,it) {
                            IntVector& c2 = it->second;
                            if(std::find(c2.begin(),c2.end(),fi) != c2.end()) {
                                forEachS(faces,j,1)
                                    eraseValue(c2,c1[faces[j]]);
                                break;
                            }
                        }
                    }
                }
            }

            sort(allDel.begin(),allDel.end());
            erase_indices(c1,allDel);

            i += nchildren;
        }
        /*renumber children ids*/
        {
            IntVector cidMap;
            cidMap.assign(mAmrTree.size(),0);
            forEach(delTree,i)
                cidMap[delTree[i]] = 1;
            Int cnt = 0;
            forEach(cidMap,i) {
                if(cidMap[i] == 0)
                    cidMap[i] = cnt++;
            }
            forEach(mAmrTree,i) {
                Node& n = mAmrTree[i];
                if(n.nchildren)
                    n.cid = cidMap[n.cid];
            }
        }
#ifdef RDEBUG
        std::cout << "Coarsening " << delTree << std::endl;
#endif
        /*delete coarsened cells*/
        sort(delTree.begin(),delTree.end());
        erase_indices(mAmrTree,delTree);
    }
    
    /*************************
    *
    *        REFINEMENT
    *
    *************************/
#ifdef RDEBUG
    std::cout << "=======================================\n";
    std::cout << "Refining " << rCells << std::endl;
#endif
        
    /***************************
     * Refine facets
     ***************************/
    IntVector startF;
    IntVector endF;
    IntVector refineF;
    Cells crefineF;
    IntVector rFacets;
    
    /*Initialize face refinement info*/
    startF.assign(mFacets.size(),0);
    endF.assign(mFacets.size(),0);
    refineF.assign(mFacets.size(),0);
    forEach(rCells,i) {
        IntVector v = mCells[rCells[i]];
        forEach(v,k) v[k] = 0;
        crefineF.push_back(v);
    }
    
    initFaceInfo(refineF,crefineF,rCells,mCells,ivBegin);
    
    /*fill refinement face set*/
    forEach(refineF,i) {
        if(refineF[i])
            rFacets.push_back(i);
    }
    
    /*set refine directions and then refine*/
    {
        IntVector rfDirs;
        rfDirs.assign(rFacets.size(),0);
        {
            IntVector rft;
            rft.assign(mFacets.size(),0);
            forEach(rCells,i) {
                Cell& c = mCells[rCells[i]];
                forEach(c,j) {
                    Int fi = c[j];
                    rft[fi] |= rDirs[i];
                }
            }
            Int cnt = 0;
            forEach(refineF,i) {
                if(refineF[i])
                    rfDirs[cnt++] = rft[i];
            }
            
            forEach(rCells,i) {
                Cell& c = mCells[rCells[i]];
                Cell& cr = crefineF[i];
                forEach(cr,j) {
                    if(cr[j] == 0)
                        cr[j] = rft[c[j]];
                }
            }
        }
        
        refineFacets(rFacets,refineF,rfDirs,startF,endF,ivBegin);
    }
    
    /*************************
     * Refine cells
     *************************/
    cellMap.assign(mBCS,0);
    forEach(cellMap,i)
        cellMap[i] = i;
    
    forEach(rCells,i) {
        Cells newc;
        Cells ncrefineF;

        Int ci = rCells[i];
        Int ndir = rDirs[i];
        newc.push_back(mCells[ci]);
        ncrefineF.push_back(crefineF[i]);
        
        Int node_index = -1;
        forEach(mAmrTree,j) {
            Node& n = mAmrTree[j];
            if((n.id == ci) && (n.nchildren == 0)) {
                node_index = j;
                break;
            }
        }

        Int start_index = mCells.size();
        Int NR = rLevel[i];
        
        for(Int m = 0;m < NR;m++) {
            
            if(m) {
                //face refinement
                {
                    IntVector rfFacets;
                    forEach(newc,j)
                        rfFacets.insert(rfFacets.end(),newc[j].begin(),newc[j].end());
                    std::sort(rfFacets.begin(), rfFacets.end());
                    rfFacets.erase(std::unique(rfFacets.begin(),rfFacets.end()), rfFacets.end());
                
                    //handle already refined faces
                    {
                        IntVector alRef;
                        forEach(rfFacets,j) {
                            Int fi = rfFacets[j];
                            if(refineF[fi])
                                alRef.push_back(j);
                        }
                        erase_indices(rfFacets,alRef);
                    }
                    
                    //init face refinement info
                    {
                        IntVector nrCells;
                        ncrefineF.clear();
                        forEach(newc,j) {
                            IntVector v = newc[j];
                            forEach(v,k) v[k] = 0;
                            ncrefineF.push_back(v);
                            nrCells.push_back(j);
                        }
                        
                        initFaceInfo(refineF,ncrefineF,nrCells,newc,ivBegin);
                    }
                    
                    //remove unrefined faces
                    {
                        IntVector alRef;
                        forEach(rfFacets,j) {
                            Int fi = rfFacets[j];
                            if(refineF[fi])
                                rFacets.push_back(fi);
                            else
                                alRef.push_back(j);
                        }
                        erase_indices(rfFacets,alRef);
                    }
                    
                    //set face refinement directions
                    IntVector rfDirs;
                    rfDirs.assign(rfFacets.size(),ndir);
                    
                    forEach(ncrefineF,i) {
                        Cell& cr = ncrefineF[i];
                        forEach(cr,j)
                            if(cr[j] == 0)
                                cr[j] = ndir;
                    }
                    //refine face
                    refineFacets(rfFacets,refineF,rfDirs,startF,endF,ivBegin);
                }
            }

            Cells newcm;
            forEach(newc,j) {
                Cell& c = newc[j];
                Cell& cr = ncrefineF[j];
                Cells newcn;
                
#ifdef RDEBUG
                std::cout << "========================================"
                             "========================================\n";
                std::cout << "Refinement level " << m + 1 << std::endl;
                std::cout << "==================\n";
                std::cout << "Cell " << c << std::endl;
                std::cout << "-------------\n";
                std::cout << "Facets of cell:\n";
                forEach(c,j)
                    std::cout << gFacets[c[j]] << std::endl;
                std::cout << "-------------\n";
                std::cout << "cr " << cr << std::endl;
#endif
                refineCell(c,cr,ndir,refineF,startF,endF,newcn,delFacets,ivBegin);
                
                /*add children to the mAmrTree*/
                {
                    Node* mnd = &mAmrTree[node_index + j];
                    mnd->nchildren = newcn.size();
                    mnd->cid = mAmrTree.size();
                    mnd->id = 0;
                    forEach(newcn,j) {
                        Node n;
                        n.id = start_index + j;
                        mAmrTree.push_back(n);
                    }
                }
                
                newcm.insert(newcm.end(),newcn.begin(),newcn.end());
                if(m == NR - 1)
                    start_index += newcn.size();
            }
            
            /*next level refinement*/
            newc = newcm;
            node_index = mAmrTree.size() - newcm.size();
    
#ifdef RDEBUG
            if(NR > 1) {
                std::cout << "==================\n";
                std::cout << newc << std::endl;
            }
#endif
        }
        
        /*add cells*/
        mCells.insert(mCells.end(),newc.begin(),newc.end());
        forEach(newc,j)
            cellMap.push_back(ci);
    }
    
    /*************************************
     *  Remove refined facets and cells
     *************************************/

    /*adjust facet indexes in cells and boundaries*/
    Int count = 0;
    forEach(delFacets,i)
        refineF[delFacets[i]] = Constants::MAX_INT - 1;
    forEach(refineF,i) {
        if(refineF[i] == 0) refineF[i] = count++;
        else if(refineF[i] == Constants::MAX_INT - 1);
        else refineF[i] = Constants::MAX_INT;
    }

    /*erase old facets and add the new ones*/
#define ERASEADD() {                                        \
    bool has;                                               \
    do {                                                    \
        has = false;                                        \
        IntVector newF,eraseF;                              \
        forEach(c,j) {                                      \
            Int fi = c[j];                                  \
            if(refineF[fi] == Constants::MAX_INT) {         \
                eraseF.push_back(j);                        \
                for(Int k = startF[fi]; k < endF[fi];k++) { \
                    newF.push_back(k);                      \
                    has = true;                             \
                }                                           \
            }                                               \
        }                                                   \
        erase_indices(c,eraseF);                            \
        c.insert(c.end(),newF.begin(),newF.end());          \
    } while (has);                                          \
                                                            \
    forEach(c,j)    {                                       \
        c[j] = refineF[c[j]];                               \
    }                                                       \
}
    forEach(mCells,i) {
        Cell& c = mCells[i];
        ERASEADD();
    }
    forEachIt(Boundaries,mBoundaries,it) {
        IntVector& c = it->second;
        ERASEADD();
    }
#undef ERASEADD
    
    /*construct faces/cells to be removed*/
    IntVector rmFacets = rFacets;
    IntVector rmCells = rCells;
    rmFacets.insert(rmFacets.end(),delFacets.begin(),delFacets.end());
    rmCells.insert(rmCells.end(),delCells.begin(),delCells.end());
    sort(rmFacets.begin(),rmFacets.end());
    sort(rmCells.begin(),rmCells.end());
    
    /*update mAmrTree IDs*/
    {
        IntVector newIDs;
        newIDs.assign(mCells.size(),0);
        forEach(rmCells,i)
            newIDs[rmCells[i]] = 1;
        Int cnt = 0;
        forEach(newIDs,i) {
            if(newIDs[i] == 0)
                newIDs[i] = cnt++;
        }
        forEach(mAmrTree,i) {
            Node& n = mAmrTree[i];
            if(!n.nchildren)
                n.id = newIDs[n.id];
        }
    }
    
    /*erase facets*/
    erase_indices(mFacets,rmFacets);
    
    /*erase cells*/
    erase_indices(mCells,rmCells);
    erase_indices(cellMap,rmCells);
    mBCS = mCells.size();
    
    /*remove unused vertices*/
    Int rem = removeUnusedVertices(ivBegin);
    ivBegin -= rem;
    
    /*break edges of faces*/
    breakEdges(ivBegin);
}

