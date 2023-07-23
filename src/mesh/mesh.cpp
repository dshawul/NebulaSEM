#include "mesh.h"

using namespace std;

//#define RDEBUG

/**
  References to global mesh
 */
namespace Mesh {
    MeshObject        gMesh;
    string&           gMeshName = gMesh.name;
    Vertices&         gVertices = gMesh.mVertices;
    Facets&           gFacets   = gMesh.mFacets;
    Cells&            gCells    = gMesh.mCells;
    Boundaries&       gBoundaries = gMesh.mBoundaries;
    IntVector&        gFOC = gMesh.mFOC;
    IntVector&        gFNC = gMesh.mFNC;
    IntVector&        gFMC = gMesh.mFMC;
    Int&              gBCS = gMesh.mBCS;
    Int&              gBCSI = gMesh.mBCSI;
    Cells&            gFaceID = gMesh.mFaceID;
    InterBoundVector& gInterMesh = gMesh.mInterMesh;
    NodeVector&       gAmrTree = gMesh.mAmrTree;
    VectorVector&     gFC = gMesh.mFC;
    VectorVector&     gCC = gMesh.mCC;
    VectorVector&     gFN = gMesh.mFN;
    ScalarVector&     gCV = gMesh.mCV;

    Vector            amr_direction(0,0,0);
    Int               is_spherical = 0;
    Scalar            sphere_radius = 6371220.0;
}
/**
  Clear mesh object
 */
void Mesh::MeshObject::clear() {
    mVertices.clear();
    mVerticesNoExtrude.clear();
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
        forEachIt(mBoundaries,it) {
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

        forEachIt(mBoundaries,it)
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
    /*recompute mFOC/mFNC*/
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
    /*add boundary cells*/
    forEachIt(mBoundaries,it) {
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
  Find the eight vertices of a hexahedron
 */
void Mesh::MeshObject::getHexCorners(const Facet& f1, const Facet& f2, IntVector& vp) {
    //pick first 4 vertices from face 0 skipping midpoints
    Int i0 = f1.size() - 1;
    Int idx = 0;
    for(Int i = 0; i < f1.size() && idx < 4; i++) {
        Int i1 = (i == f1.size() - 1) ? 0 : i + 1;
        Vector v0 = unit(mVertices[f1[i]] - mVertices[f1[i0]]);
        Vector v1 = unit(mVertices[f1[i1]] - mVertices[f1[i]]);
        Scalar dt = dot(v0,v1);
        if(dt > 1) dt = 1;
        else if(dt < -1) dt = -1;
        Scalar angle = acos(dt);
        const Scalar tol = Constants::PI / 32;
        if(!(angle < tol || angle >= Constants::PI - tol)) {
            vp.push_back(f1[i]);
            i0 = i;
        }
    }
    //pick the other vertices by distance from the first 4
    for(Int i = 0; i < 4; i++) {
        Scalar mind = 1e20;
        Int minj = 0;
        forEach(f2,j) {
            Scalar dist = mag(mVertices[f2[j]] - mVertices[vp[i]]);
            if(dist < mind) {
                mind = dist;
                minj = j;
            }
        }
        vp.push_back(f2[minj]);
    }
}
/**
  Fix hexahedral cells for the sake of DG
 */
void Mesh::MeshObject::fixHexCells() {
    mFMC.assign(mFacets.size(),0);
    mFaceID.clear();

    forEach(mCells,i) {
        Cell& c = mCells[i];
        if(c.size() == 1) {
            /*boundary cell*/
            IntVector b;
            b.push_back(0);
            mFaceID.push_back(b);
        } else {
            /*hex cells with >=6 faces*/
#ifdef RDEBUG
            using namespace Util;
            Cell co;
            forEach(c,i)
                co.push_back(c[i]);
            bool print = true;
            if(print) {
                std::cout << "======================" << std::endl;
                std::cout << co << std::endl;
                forEach(co,i) {
                    Facet& f = mFacets[co[i]];
                    std::cout << "Facet " << co[i] << ": " << f << std::endl;
                    forEach(f,j)
                        std::cout << mVertices[f[j]] << std::endl;
                }
                std::cout << std::endl;
            }
#endif
            /*Group coplanar faces*/
            std::vector<int> groupId(6);
            Cells cng;
            Facets fng;
            for(Int i = 0; i < 6; i++) {
                IntVector indices;
                Cell ci;
            
                ci.push_back(c[0]);
                indices.push_back(0);
                forEachS(c,j,1) {
                    bool coplanar = coplanarFaces(mFacets[c[0]],mFacets[c[j]]);
                    if(coplanar) {
                        ci.push_back(c[j]);
                        indices.push_back(j);
                    }
                }
                erase_indices(c,indices);
                cng.push_back(ci);
            
#ifdef RDEBUG
                using namespace Util;
                std::cout << "===============\n";
                std::cout << "ci: " << ci << std::endl;
                std::cout << "ind: " << indices << std::endl;
                forEach(ci,j) {
                    Facet& f = mFacets[ci[j]];
                    std::cout << f << std::endl;
                    forEach(f,i)
                        std::cout << i << ". " << mVertices[f[i]] << std::endl;
                }
#endif

                Facet fn;
                mergeFacetsGroup(ci,fn);

                fng.push_back(fn);
                groupId[i] = -1;
            }
#ifdef RDEBUG
            if(print) {
                std::cout << "cng: " << cng << std::endl;
                std::cout << "fng: " << fng << std::endl;
            }
#endif
            
            /*order faces so that 0,1 are opposite one another
            same for 2,3 and 4,5*/
            Facet& f0 = fng[0];
            forEach(fng,j) {
                if(groupId[j] >= 0) continue;
            
                Facet& fj = fng[j];
                Int id = j;
                if(j >= 1) {
                    //set face ids
                    int local_id = -1;
                    do {
                        //xz-0 face 2
                        Int count = 0;
                        forEach(fj, k) {
                            if(fj[k] == f0[0] ||
                               fj[k] == f0[1])
                                count++;
                        }
                        if(count >= 2) {
                            local_id = 2;
                            break;
                        }
                        //yz-0 face 4
                        count = 0;
                        forEach(fj, k) {
                            if(fj[k] == f0[0] ||
                               fj[k] == f0[f0.size() - 1])
                                count++;
                        }
                        if(count >= 2) {
                            local_id = 4;
                            break;
                        }
                    } while(0);
                    if(local_id < 0) continue;
                    id = local_id;
                }
            
                groupId[j] = id;
                forEach(fng,k) {
                    if(groupId[k] >= 0) continue;
            
                    Facet& fk = fng[k];
                    bool has_shared = false;
                    forEach(fj, j1) {
                        forEach(fk, k1) {
                            if(fj[j1] == fk[k1]) {
                                has_shared = true;
                                break;
                            }
                        }
                    }
                    if(!has_shared) {
                        groupId[k] = (groupId[j] ^ 1);
                        break;
                    }
                }
            }
            
#ifdef RDEBUG
            if(print) {
                std::cout << "groupId: " << groupId << std::endl;
            }
#endif

            /*fix orientation of vertices */
            IntVector vp;
            getHexCorners(fng[0], fng[1], vp);
            IntVector rots  = {vp[0], vp[4], vp[0], vp[3], vp[0], vp[1]};
            IntVector rote  = {vp[1], vp[5], vp[1], vp[2], vp[3], vp[2]};
            for(Int i = 0; i < 6; i++) {
                auto it = std::find(fng[i].begin(), fng[i].end(), rots[i]);
                std::rotate(fng[i].begin(), it, fng[i].end());

                Scalar dir = dot(unit(mVertices[fng[i][1]] - mVertices[fng[i][0]]),
                                 unit(mVertices[rote[i]] - mVertices[rots[i]]));
                if(dir < 0.99)
                    std::reverse(fng[i].begin()+1, fng[i].end());
            }

            for(Int i = 0; i < 6; i++) {
                Facet& fngo = fng[i];
                Cell& cn = cng[i];

                Facets newf;
                newf.resize(cn.size());
                forEach(fngo,k) {
                    Int vi = fngo[k];
                    forEach(cn,j) {
                        Facet& f = mFacets[cn[j]];
                        if(std::find(f.begin(),f.end(),vi) != f.end())
                            newf[j].push_back(vi);
                    }
                }

                forEach(cn,j) {
                    Facet& f = mFacets[cn[j]];
                    f = newf[j];
                }
            }

            /*for new cell and face IDs*/
            IntVector b;
            for(Int idx = 0; idx < 6; idx++) {
                for(Int i = 0; i < 6; i++) {
                    Int id = groupId[i];
                    if(id != idx) continue;
            
                    Cell& cn = cng[i];
                    forEach(cn,j) {
                        c.push_back(cn[j]);
                        b.push_back(id);
                    }
                }
            }
            mFaceID.push_back(b);
            
#ifdef RDEBUG
            if(print) {
                std::cout << "faceId: " << b << std::endl;
                std::cout << "cell: " << c << std::endl;
                std::cout << "======================" << std::endl;
                if(!equalSet(co,c)) {
                    std::cout << "fixHexCell failed" << std::endl;
                    exit(0);
                }
            }
#endif
       }

       /*for non-conformal*/
       if(c.size() > 6) {
           IntVector& fid = mFaceID[i];
           forEach(c,j) {
               Int f = c[j];
               Int cnt = std::count(fid.begin(), fid.end(), fid[j]);
               if(cnt > 1) {
                   if(mFOC[f] == i)
                       mFMC[f] = 2;
                   else
                       mFMC[f] = 1;
               }
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

    /* face centre: centroid of face vertices*/
    forEach(mFacets,i) {
        Facet& f = mFacets[i];
        Vector C(0);
        forEach(f,j)
            C += mVertices[f[j]];
        mFC[i] = C / Scalar(f.size());
    }
    /* cell centre: centroid of cell vertices */
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
        /*corrected face center: which is the centroid of the solid polygon
          instead of the vertices alone*/
        mFC[i] = C / Ntot;
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
        /*corrected cell center: which is the centroid of the solid control volume
          instead of the vertices alone*/
        mCC[i] = C / V;
        mCV[i] = V / Scalar(3);
        /*correct for sphere*/
        if(is_spherical) {
            Scalar radiusb = mag(mVertices[mFacets[c[0]][0]]);
            Scalar radiust = mag(mVertices[mFacets[c[1]][0]]);
            mCC[i] = ((radiusb + radiust) / (2 * mag(mCC[i]))) * mCC[i];
            mFC[c[0]] = (radiusb / mag(mFC[c[0]])) * mFC[c[0]];
            mFC[c[1]] = (radiust / mag(mFC[c[1]])) * mFC[c[1]];
            for(Int j = 2; j < 6; j++)
                mFC[c[j]] = ((radiusb + radiust) / (2 * mag(mFC[c[j]]))) * mFC[c[j]];
        }
    }
    /*boundary cell centre and volume*/
    forEachS(mCells,i,mBCS) {
        Int fi = mCells[i][0];
        mCV[i] = mCV[mFOC[fi]];
        mCC[i] = mFC[fi];
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
    forEachIt(mBoundaries,it) {
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
        Vector v1 = mVertices[f[f.size() - 1]];
        forEach(f,j) {
            Vector v2 = mVertices[f[j]];
            Vector v3 = mVertices[f[ (j + 1 == f.size()) ? 0 : (j + 1) ]];
            if(!pointInLine(v2,v1,v3))
                isUsed[f[j]] = 1;
            v1 = v2;
        }
    }
    forEach(mFacets,i) {
        Facet& f = mFacets[i];
        IntVector removed;
        forEach(f,j) {
            if(!isUsed[f[j]])
                removed.push_back(j);
        }
        if(removed.size())
            erase_indices(f,removed);
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
 Extrude mesh
 */
void Mesh::MeshObject::ExtrudeMesh() {
    if(is_spherical) {
        if(mVerticesNoExtrude.size()) {
            forEach(mVerticesNoExtrude,i)
                mVertices[i] = mVerticesNoExtrude[i];
        } else {
            mVerticesNoExtrude.resize(mVertices.size());
            forEach(mVerticesNoExtrude,i)
                mVerticesNoExtrude[i] = mVertices[i];
        }

        Scalar minh = 1e30, maxh = 0;
        forEach(mVertices,i) {
           Vertex& v = mVertices[i];
           Scalar h = max(max(fabs(v[0]),fabs(v[1])),fabs(v[2]));
           if(h > maxh) maxh = h;
           if(h < minh) minh = h;
        }

        forEach(mVertices,i) {
           Vertex& v = mVertices[i];
           Scalar h = max(max(fabs(v[0]),fabs(v[1])),fabs(v[2]));
           Scalar f = (h - minh) / (maxh - minh);
           constexpr Scalar radiusi = 6371220;
           constexpr Scalar radiuso = 6381220;
           v = unit(v) * ((f) * radiusi + (1-f) * radiuso);
        }
    }
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
                sort(order.begin(), order.end(), Util::compare_pair_second<>());

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
        Scalar dt = dot(p1,p2) / mg;
        if(dt > 1) dt = 1;
        else if(dt < -1) dt = -1;
        sum += acos(dt);
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
    const Vector& v2 = mVertices[f[1]];
    for(Int j = 1;j < f.size();j++) {
        const Vector& v3 = mVertices[f[f.size() - j]];

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
        Vector v = mVertices[f2[1]] - mVertices[f1[0]];
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
    Facet f1 = f1_;
    Facet f2 = f2_;
    //make them same direction
    {
        Vector N1,N2;
        calcUnitNormal(f1_,N1);
        calcUnitNormal(f2_,N2);
        if(dot(N1,N2) < 0)
            std::reverse(f2.begin()+1,f2.end());
    }
    //rotate nodes of faces to first non-shared node
    Int contained = false;
    {
        int v1 = -1,v2 = -1;
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
        rotate(f1.begin(),f1.begin() + v1,f1.end());
        if(v2 != -1)
            rotate(f2.begin(),f2.begin() + v2,f2.end());
        else
            contained = true;
    }   

    //find start and end of shared edges
    Int a[2];
    Int b[2];
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

    //make sure first vertex is a corner
    while(pointInLine(mVertices[f[0]],mVertices[f[f.size() - 1]], mVertices[f[1]]))
        std::rotate(f.begin(),f.begin() + 1, f.end());

    return true;
}
/**
  Merge facets of cell
 */
void Mesh::MeshObject::mergeFacetsGroup(const IntVector& shared1,Facet& fn, const Cell* c1) {
    IntVector ci;
    if(c1) {
        forEach(shared1,i)
            ci.push_back((*c1)[shared1[i]]);
    } else
        ci = shared1;

    fn = mFacets[ci[0]];
    ci.erase(ci.begin());
    bool repeat = false;
    do {
        IntVector merged;
        repeat = false;
        forEach(ci,j) {
            Facet& f = mFacets[ci[j]];
            Facet fm;
            if(mergeFacets(fn,f,fm)) {
                fn = fm;
                merged.push_back(j);
            } else
                repeat = true;
        }
        erase_indices(ci,merged);
    } while(repeat);

#ifdef RDEBUG
    if(ci.size())  {
        using namespace Util;
        cout << "Merge failed.\n";
        cout << "=======================\n";
        cout << c1 << endl;
        cout << shared1 << endl;
        forEach(shared1,i)
            cout << mFacets[ci[i]] << endl;
        cout << "-----------------------\n";
        {
            Facet r,rr;
            straightenEdges(fn,r,rr);
            cout << "Merged face  : " << fn << endl;
            cout << "Straight Edge: " << r << endl;
            cout << "Removed      : " << rr << endl;
        }
        forEach(fn,i)
            cout << mVertices[fn[i]] << " : ";
        cout << endl;
        exit(0);
    }
#endif
}
/**
  Merge two cells
 */
void Mesh::MeshObject::mergeCells(Cell& c1, const Cell& c2, IntVector& delFacets) {
    forEach(c2,m) {
        auto it = std::find(c1.begin(), c1.end(), c2[m]);
        if(it == c1.end())
            c1.push_back(c2[m]);
        else {
            delFacets.push_back(c2[m]);
            c1.erase(it);
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

    const bool refine3D = equal(Scalar(0.0),mag(amr_direction));

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
    Scalar dt = dot(mDir,uLine);            \
    if(dt > 1) dt = 1;                      \
    else if(dt < -1) dt = -1;               \
    Scalar angle = acos(dt);                \
    if(angle < Constants::PI / 4 ||         \
       angle >= 3 * Constants::PI / 4)      \
    continue;                               \
}

        Vector Ce = (v1 + v2) / 2.0;
        Vector uLine = unit(v2 - v1);

        if(!refine3D) {
            const Vector uDir = is_spherical ? unit(Ce) : unit(amr_direction);
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
        Facets mf;
        forEachS(f,j,fmid) {
            Facet fn;
            Int rot = j;
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

            /*rotate*/
            if(rot == f.size() - 1) {
                auto it = std::find(fn.rbegin(), fn.rend(), f[0]);
                std::rotate(fn.begin(), it.base() - 1, fn.end());
            } else
                std::rotate(fn.rbegin(), fn.rbegin() + rot, fn.rend());
        
            mf.push_back(fn);
        }
        newf.push_back(mf[mf.size() - 1]);
        for(Int j = 0; j < mf.size() - 1; j++)
            newf.push_back(mf[j]);

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

    using namespace Util;

    /* ***************************************
     *
     * Merge and then refine co-planar faces
     *
     * ***************************************/
    map<Int, pair<IntVector,IntVector> > sharedMap;

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

    forEachIt(sharedMap,it) {
        IntVector& shared = (it->second).first;

        Facet f;
        if(shared.size() > 1) {
            mergeFacetsGroup(shared,f,&c);
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
        cout << it->first << " = " << it->second.first << " = " << it->second.second <<  endl;
        cout << "----" << endl;
#endif

    }

    /* *************************
     *
     * Form new cells (pyramids) 
     *
     * *************************/
    Vector C;
    calcCellCenter(c,C);
    mVertices.push_back(C);
    Int cci = mVertices.size() - 1;
    Int ifBegin = mFacets.size();

#ifdef RDEBUG
    cout << "Cell center: " << C << endl;
#endif

    IntVector ownerf;
    forEach(c1,j) {
        Int fi = c1[j];
        Int crj = c1r[j];

        //cell is refined but face is not?
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
        if(j == 0) cout << "====================================\n";
        cout << "------------------------" << endl;
        cout << "Face split id:  " << fi << " rdir " << rDir << " cr " << crj 
             << " subfaces ( " << startF[fi] << " to " << endF[fi] << " ) " << endl;
#endif

        //refine cells
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
            cout << "subface " << fni << " ";
            cout << f << " fC " << fc << endl;
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
                    if(equalSet(fn,mFacets[k])) {
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

    /* *************************
     *
     * Replace co-planar faces
     *
     * *************************/
    forEach(newc,i) {
        if(ownerf[i] < Constants::MAX_INT / 2)
            continue;

        Cell& nc = newc[i];
        Int fi = nc[0];
        Facet& nf = mFacets[fi];

        IntVector forg;
        pair<IntVector,IntVector> pair;
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
#ifdef RDEBUG
        cout << "============\n";
        cout << " New cells before merge " << newc.size() << endl;
        cout << newc << endl;
        cout << "============\n";
#endif

        // Find among new cells on different owner faces,
        // and merge if they have a shared face
        Cells mergec;
        forEach(newc,j) {
            Int o1 = ownerf[j];
            Cell& c1 = newc[j];
            IntVector mg;
            forEachS(newc,k,j+1) {
                Int o2 = ownerf[k];
                Cell& c2 = newc[k];

                if(o1 == o2) continue;

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

        //merge cells
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
        cout << " New cells " << newc.size() << endl;
        cout << newc << endl;
        cout << "============\n";
#endif

        //two cells should share only one face
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
                    mergeFacetsGroup(shared1,nf,&c1);

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
    cout << " New cells after merge " << newc.size() << endl;
    cout << newc << endl;
    forEach(newc, i) {
        Cell& c = newc[i];
        std::cout << "Cell " << c << std::endl;
        forEach(c,j) {
            Facet& f = mFacets[c[j]];
            std::cout << j << ". " << f << std::endl;
            forEach(f,k)
                std::cout << "    " << k << ". " << mVertices[f[k]] << std::endl;
            if(f.size() <= 3) {
                cout << "Cell still has triangular face." << endl;
                exit(0);
            }
        }
        std::cout << "-----" << std::endl;
    }
#endif

}
/**
  Initialize facet refinement information
 */
void Mesh::MeshObject::initFaceInfo(IntVector& refineF,Cells& crefineF,
        const IntVector& rCells,const Cells& newCells) {

    forEach(rCells,i) {
        const Cell& c = newCells[rCells[i]];
        Cell& cr = crefineF[i];

        IntVector flag;
        flag.assign(c.size(),0);

        forEach(c,j) {
            if(flag[j]) continue;

            IntVector shared;
            shared.push_back(j);

            bool coplanar = false;
            forEachS(c,k,j+1) {
                if(flag[k]) continue;

                if(coplanarFaces(mFacets[c[j]],mFacets[c[k]])) {
                    coplanar = true;
                    shared.push_back(k);
                }
            }
            if(!coplanar) {
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
void Mesh::MeshObject::refineMesh(const IntVector& cCells,const IntVector& rCells,const IntVector& rLevel,
        const IntVector& rDirs,IntVector& refineMap,IntVector& coarseMap,IntVector& cellMap) {

    using namespace Util;
    Int ivBegin = mVertices.size();
    IntVector delFacets;
    IntVector delCells;

    /*************************
     *
     *        COARSENING
     *
     *************************/
#ifdef RDEBUG
    cout << "=======================================\n";
    IntVector cCellsL;
    forEach(cCells,i)
        if(cCells[i]) cCellsL.push_back(i);
    cout << "Coarsening " << cCellsL << endl;
#endif

    {
        Int oldgBCS = mCells.size();

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
                    Cell c1;
                    Int pid = mCells.size();
                    coarseMap.push_back(n.nchildren);
                    coarseMap.push_back(pid);
                    for(Int j = 0;j < n.nchildren;j++) {
                        Node& cn = mAmrTree[n.cid + j];
                        Cell& c2 = mCells[cn.id];
                        mergeCells(c1,c2,delFacets);
                        delCells.push_back(cn.id);
                        coarseMap.push_back(cn.id);
                        delTree.push_back(n.cid + j);
                        /*update owner/neighbor info*/
                        forEach(c2,k) {
                            Int fi = c2[k];
                            if(mFOC[fi] == cn.id)
                                mFOC[fi] = pid;
                            else
                                mFNC[fi] = pid;
                        }
                    }
                    mCells.push_back(c1);
                    n.nchildren = 0;
                    n.id = pid;
                    n.cid = 0;
                }
            }
        }

        /*update boundary neighbors*/
        {
            Int inc = mCells.size() - oldgBCS;
            forEachIt(mBoundaries,it) {
                IntVector& mB = it->second;
                forEach(mB,i)
                    mFNC[mB[i]] += inc;
            }
        }

        /*merge facets*/
        forEach(coarseMap,i) {
            Int nchildren = coarseMap[i];
            Int ci = coarseMap[i + 1];
            Cell& c1 = mCells[ci];

#ifdef RDEBUG
            cout << "-----------------------------" << endl;
            cout << "Coarsening faces of cell " << ci << " cC " << mCC[coarseMap[i + 2]]  << endl;
            cout << "children are: ";
            for(Int j = 0; j < nchildren; j++)
                cout << coarseMap[i + 2 + j] << " ";
            cout << std::endl;
            cout << c1 << std::endl;
#endif

            map<Int,IntVector> sharedMap;

            forEach(c1,j) {
                Int f1 = c1[j];
                Int o1 = (mFOC[f1] == ci) ? mFNC[f1] : mFOC[f1];
                if(o1 >= mCells.size()) {
                    Vector u = unit(mFN[f1]);
                    Vector v0 = mVertices[mFacets[f1][1]];
                    Scalar d = dot(v0,u);
                    o1 = round(Scalar(10000.0) * d +
                               Scalar(1000.0) * fabs(u[2]) +
                               Scalar(100.0) * fabs(u[1]) +
                               Scalar(10.0) * fabs(u[0]));
                    o1 += Int(3 * (Constants::MAX_INT >> 2));
                }
                IntVector& val = sharedMap[o1];
                val.push_back(j);
            }

#ifdef RDEBUG
            cout << "Face sharing map with neighbor cell " << std::endl;
            cout << sharedMap;
#endif

            IntVector allDel;
            forEachIt(sharedMap,it) {
                IntVector& faces = it->second;
                if(faces.size() > 1) {
                    //merge facets of cell
                    Facet nf;
                    mergeFacetsGroup(faces,nf,&c1);

                    //merge into first face
                    Int fi = c1[faces[0]];
                    mFacets[fi] = nf;

                    //mark old faces for deletion
                    forEachS(faces,j,1)
                        delFacets.push_back(c1[faces[j]]);
                    allDel.insert(allDel.end(),faces.begin() + 1,faces.end());

                    if(it->first < Constants::MAX_INT / 2) {
                        //erase old faces from neighbor cell
                        Cell& c2 = mCells[it->first];
                        forEachS(faces,j,1)
                            eraseValue(c2,c1[faces[j]]);
                    } else {
                        //erase old faces from neighbor boundary cell
                        forEachIt(mBoundaries,it) {
                            IntVector& c2 = it->second;
                            //Assume all faces are on same patch [Fix me later]
                            if(find(c2.begin(),c2.end(),fi) != c2.end()) {
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

#ifdef RDEBUG
            cout << "-----------------------------" << endl;
            cout << "Cell after coarsening  " << endl;
            cout << c1 << std::endl;
#endif
            i += nchildren + 1;
        }
        /*renumber children ids*/
        {
            IntVector cidMap;
            cidMap.assign(mAmrTree.size(),0);
            forEach(delTree,i)
                cidMap[delTree[i]] = Constants::MAX_INT;
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
        cout << "Coarsening tree to delete " << delTree << endl;
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
    cout << "=======================================\n";
    cout << "Refining " << rCells << endl;
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

    initFaceInfo(refineF,crefineF,rCells,mCells);

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
                forEach(c,j)
                    rft[c[j]] |= rDirs[i];
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

        /*Note: co-planar faces are not refined here*/
        refineFacets(rFacets,refineF,rfDirs,startF,endF,ivBegin);
    }

    /*************************
     * Refine cells
     *************************/
    forEach(rCells,i) {
        Int rci = rCells[i];
        Int rdir = rDirs[i];
        Int rlevel = rLevel[i];

        //find tree node for this cell
        Int start_index = mCells.size();
        Int node_index = Constants::MAX_INT - 1;
        forEach(mAmrTree,j) {
            Node& n = mAmrTree[j];
            if((n.id == rci) && (n.nchildren == 0)) {
                node_index = j;
                break;
            }
        }

        //initialize recursive refinemet of the cell
        Cells newc;
        Cells ncrefineF;
        newc.push_back(mCells[rci]);
        ncrefineF.push_back(crefineF[i]);

        //conduct multiple refinement levels
        for(Int m = 0; m < rlevel; m++) {

            //face refinement for multiple levels only
            if(m) {

                //form unique set of faces to refine
                IntVector rfFacets;
                forEach(newc,j)
                    rfFacets.insert(rfFacets.end(),newc[j].begin(),newc[j].end());
                sort(rfFacets.begin(), rfFacets.end());
                rfFacets.erase(unique(rfFacets.begin(),rfFacets.end()), rfFacets.end());

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

                    initFaceInfo(refineF,ncrefineF,nrCells,newc);
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
                rfDirs.assign(rfFacets.size(),rdir);

                forEach(ncrefineF,i) {
                    Cell& cr = ncrefineF[i];
                    forEach(cr,j)
                        if(cr[j] == 0)
                            cr[j] = rdir;
                }

                //refine the new set of facets
                refineFacets(rfFacets,refineF,rfDirs,startF,endF,ivBegin);
            }

            //cell refinement
            Cells newcm;
            forEach(newc,j) {
                Cell& c = newc[j];
                Cell& cr = ncrefineF[j];

#ifdef RDEBUG
                cout << "========================================"
                        "========================================\n";
                cout << "Refinement level " << m + 1 << endl;
                cout << "==================\n";
                cout << "Cell " << c << endl;
                cout << "-------------\n";
                cout << "Facets of cell:\n";
                forEach(c,j) {
                    Facet& f = gFacets[c[j]];
                    cout << "face " << c[j] << " " << f << endl;
                    forEach(f,k)
                        cout << gVertices[f[k]] << endl;
                }
                cout << "-------------\n";
                cout << "cr " << cr << endl;
#endif

                //refine cell
                Cells newcn;
                refineCell(c,cr,rdir,refineF,startF,endF,newcn,delFacets,ivBegin);

                //add children to the mAmrTree
                {
                    Node* mnd = &mAmrTree[node_index + j];
                    mnd->nchildren = newcn.size();
                    mnd->cid = mAmrTree.size();
                    mnd->id = 0;
                    forEach(newcn,k) {
                        Node n;
                        n.id = start_index + k;
                        n.level = mnd->level + 1;
                        mAmrTree.push_back(n);
                    }
                }

                newcm.insert(newcm.end(),newcn.begin(),newcn.end());
                if(m == rlevel - 1)
                    start_index += newcn.size();
            }

            //next level refinement
            newc = newcm;
            node_index = mAmrTree.size() - newcm.size();

#ifdef RDEBUG
            if(rlevel > 1) {
                cout << "==================\n";
                cout << newc << endl;
            }
#endif
        }

        /*add cells*/
        Int cid = mCells.size();
        mCells.insert(mCells.end(),newc.begin(),newc.end());
        refineMap.push_back(newc.size());
        refineMap.push_back(rci);
        forEach(newc,j)
            refineMap.push_back(cid + j);
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
                                                            \
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
    IntVector eraseF;                                       \
    forEach(c,j)    {                                       \
        Int fi = c[j];                                      \
        if(refineF[fi] != Constants::MAX_INT - 1)           \
            c[j] = refineF[fi];                             \
        else                                                \
            eraseF.push_back(j);                            \
    }                                                       \
    erase_indices(c,eraseF);                                \
}
    forEach(mCells,i) {
        Cell& c = mCells[i];
        ERASEADD();
    }
    forEachIt(mBoundaries,it) {
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
        cellMap.assign(mCells.size(),0);
        forEach(rmCells,i)
            cellMap[rmCells[i]] = Constants::MAX_INT;
        Int cnt = 0;
        forEach(cellMap,i) {
            if(cellMap[i] == 0)
                cellMap[i] = cnt++;
        }
        forEach(mAmrTree,i) {
            Node& n = mAmrTree[i];
            if(!n.nchildren)
                n.id = cellMap[n.id];
        }
    }
    
    /*erase facets*/
    erase_indices(mFacets,rmFacets);
    
    /*erase cells*/
    erase_indices(mCells,rmCells);
    mBCS = mCells.size();
    
    /*remove unused vertices*/
    Int rem = removeUnusedVertices(ivBegin);
    ivBegin -= rem;
    
    /*break edges of faces*/
    breakEdges(ivBegin);
}

