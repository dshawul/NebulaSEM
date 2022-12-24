#include "hexMesh.h"

using namespace Mesh;

/**
  Calculate the circum-circumcenter of triangel (circumcenter of circumscribed circle)
  or crossing point of two perpendicular bisectors. Used for ARC edges.
 */
Vector circumcenter(const Vector& v1,const Vector& v2,const Vector& v3) {
    Vector v12 = v1 - v2;
    Vector v13 = v1 - v3;
    Vector v23 = v2 - v3;
    Scalar d = 2 * magSq(v12 ^ v23);
    Scalar a = magSq(v23) * (v12 & v13) / d;
    Scalar b = magSq(v13) * (-v12 & v23) / d;
    Scalar c = magSq(v12) * (v13 & v23) / d;
    return a * v1 + b * v2 + c * v3;
}
/**
  Add a vertex following shape of a given edge
 */
void ADDV(int w,Scalar m,Edge* edges,Vector* vd) {
    Edge& e = edges[w];
    if(e.type == NONE) {
        vd[w] = (1 - m) * e.v[0] + (m) * e.v[1];
    } else if(e.type == ARC) {
        vd[w] = rotate(e.v[0] - e.v[3],e.N,e.theta * m) + e.v[3];
    } else if(e.type == COSINE) {
        vd[w] = (1 - m) * e.v[0] + (m) * e.v[1] + 
            pow(cos(Constants::PI * (m - 0.5)),2) * e.N;
    } else if(e.type == QUAD) {
        vd[w] = (1 - m) * e.v[0] + (m) * e.v[1] + 
            (4 * m * (1 - m)) * e.N;
    }
}
/**
  Generate hexahedral mesh
 */
void hexMesh(Int* n,Scalar* s,Int* type,Vector* vp,Edge* edges,MeshObject& mo) {
    Int i,j,k,m;

    /*for wall division set twice 
      number of divisions requested*/
    for(j = 0;j < 3;j++) {
        bool found = false;
        for(i = j;i < 12;i+=3) {
            if(type[i] == WALL) {
                if((n[j] % 2) && (n[j] != 1)) {
                    found = true;
                    break;
                }
            }
        }
        if(found) {
            n[j]++;
            for(i = j;i < 12;i+=3) {
                s[i] = 1 / s[i];
            }
        }
    }

    /*calculate scale*/
    Scalar* sc[12];
    for(i = 0;i < 12;i++) {
        Int nt = n[i / 4];
        sc[i] = new Scalar[nt + 1];
        if(type[i] == WALL)
            s[i] = pow(s[i],Scalar(1./(nt / 2.)));
        else
            s[i] = pow(s[i],Scalar(1./nt));
    }
    for(i = 0;i < 12;i++) {
        Int nt = n[i / 4];
        Scalar r = s[i];
        if(nt == 1) {
            sc[i][0] = 0;
            sc[i][1] = 1;
        } else {
            if(type[i] == WALL)
                nt /= 2;
            for(j = 0;j <= nt;j++) {
                if(equal(r,Scalar(1))) 
                    sc[i][j] = Scalar(j) / (nt);
                else
                    sc[i][j] = (1 - pow(r,Scalar(j))) / (1 - pow(s[i],Scalar(nt)));
            }
            if(type[i] == WALL) {
                for(j = 0;j <= nt;j++)
                    sc[i][j] /= 2;
                for(j = 0;j <= nt;j++)
                    sc[i][j + nt] = Scalar(1.0) - sc[i][nt - j];
            }
        }
    }
    for(i = 0;i < 12;i++) {
        Edge& e = edges[i];
        if(e.type == ARC) {
            Vector C = circumcenter(e.v[0],e.v[1],e.v[2]);
            Vector r1 = e.v[0] - C;
            Vector r2 = e.v[1] - C;
            e.theta = acos((r1 & r2) / (mag(r1) * mag(r2)));
            e.v[3] = C;
            e.N = (e.v[2] - e.v[0]) ^ (e.v[1] - e.v[0]);
            e.N = unit(e.N);
        } else if(e.type == COSINE || e.type == QUAD) {
            Vector mid = (e.v[1] + e.v[0]) / 2;
            e.N = e.v[2] - mid;
        }
    }
    /*variables*/
    Int nx = n[0] + 1 , ny = n[1] + 1 , nz = n[2] + 1;
    const Int B1 = (nx - 0) * (ny - 1) * (nz - 1);
    const Int B2 = (nx - 1) * (ny - 0) * (nz - 1);
    const Int B3 = (nx - 1) * (ny - 1) * (nz - 0);
    IntVector VI(nx * ny * nz,0);
    IntVector FI(B1 + B2 + B3,0);

    /*vertices*/
    Vertex v,v1,v2,vd[12],vf[6];
    Scalar rx,ry,rz;

#define I0(i,j,k)  (i * ny * nz + j * nz + k)

#define ADDF(w,rr,rs,i00,i01,i10,i11,ir0,ir1,i0s,i1s) { \
    vf[w] = Interpolate_face(                           \
            rr,rs,                                      \
            vp[i00],vp[i01],vp[i10],vp[i11],            \
            vd[ir0],vd[ir1],vd[i0s],vd[i1s]);           \
}

#define ADDC() {                                        \
    v = Interpolate_cell(                               \
            rx,ry,rz,                                   \
            vp[0],vp[4],vp[3],vp[7],                    \
            vp[1],vp[5],vp[2],vp[6],                    \
            vd[0],vd[3],vd[1],vd[2],                    \
            vd[4],vd[7],vd[5],vd[6],                    \
            vd[8],vd[11],vd[9],vd[10],                  \
            vf[4],vf[5],vf[2],vf[3],vf[0],vf[1]);       \
}

#define ADD() {                                     \
    ADDV(0,sc[0][i],edges,vd);                      \
    ADDV(1,sc[1][i],edges,vd);                      \
    ADDV(2,sc[2][i],edges,vd);                      \
    ADDV(3,sc[3][i],edges,vd);                      \
    ADDV(4,sc[4][j],edges,vd);                      \
    ADDV(5,sc[5][j],edges,vd);                      \
    ADDV(6,sc[6][j],edges,vd);                      \
    ADDV(7,sc[7][j],edges,vd);                      \
    ADDV(8,sc[8][k],edges,vd);                      \
    ADDV(9,sc[9][k],edges,vd);                      \
    ADDV(10,sc[10][k],edges,vd);                    \
    ADDV(11,sc[11][k],edges,vd);                    \
    rx = i / Scalar(nx - 1);                        \
    ry = j / Scalar(ny - 1);                        \
    rz = k / Scalar(nz - 1);                        \
    ADDF(0, rx,ry, 0,3,1,2, 0,1,4,5);               \
    ADDF(1, rx,ry, 4,7,5,6, 3,2,7,6);               \
    ADDF(2, rx,rz, 0,4,1,5, 0,3,8,9);               \
    ADDF(3, rx,rz, 3,7,2,6, 1,2,11,10);             \
    ADDF(4, ry,rz, 0,4,3,7, 4,7,8,11);              \
    ADDF(5, ry,rz, 1,5,2,6, 5,6,9,10);              \
    ADDC();                                         \
};

    /*interior*/
    for(j = 1;j < ny - 1;j++) {
        for(i = 1;i < nx - 1;i++) {
            for(k = 1;k < nz - 1;k++) {
                ADD();
                mo.mVertices.push_back(v);
                VI[I0(i,j,k)] = mo.mVertices.size() - 1;
            }
        }
    }
    mo.mNV = mo.mVertices.size();
    
    /*boundaries*/
    for(i = 0;i < nx; i += (nx - 1)) {
        for(j = 0;j < ny;j++) {
            for(k = 0;k < nz;k++) {
                ADD();
                mo.mVertices.push_back(v);
                VI[I0(i,j,k)] = mo.mVertices.size() - 1;
            }
        }
    }
    for(j = 0;j < ny; j += (ny - 1)) {
        for(i = 1;i < nx - 1;i++) {
            for(k = 0;k < nz;k++) {
                ADD();
                mo.mVertices.push_back(v);
                VI[I0(i,j,k)] = mo.mVertices.size() - 1;
            }
        }
    }
    for(k = 0;k < nz; k += (nz - 1)) {
        for(i = 1;i < nx - 1;i++) {
            for(j = 1;j < ny - 1;j++) {
                ADD();
                mo.mVertices.push_back(v);
                VI[I0(i,j,k)] = mo.mVertices.size() - 1;
            }
        }
    }
/*end*/
#undef ADD
#undef ADDF
#undef ADDE

    delete[] sc[0];
    delete[] sc[1];
    delete[] sc[2];

/*faces*/
#define I1(i,j,k)   (i * (ny - 1) * (nz - 1) + j * (nz - 1) + k)
#define I2(i,j,k)   (i * (ny - 0) * (nz - 1) + j * (nz - 1) + k + B1)
#define I3(i,j,k)   (i * (ny - 1) * (nz - 0) + j * (nz - 0) + k + B1 + B2)

#define ADD(a1,a2,a3,a4) {                          \
    Facet f;                                        \
    m = I0(i,j,k);                                  \
    f.push_back(VI[a1]);                            \
    f.push_back(VI[a2]);                            \
    f.push_back(VI[a3]);                            \
    f.push_back(VI[a4]);                            \
    mo.mFacets.push_back(f);                        \
};

    /*interior*/
    for(i = 0;i < nx - 1;i++) {
        for(j = 0;j < ny - 1;j++) {
            for(k = 1;k < nz - 1;k++) {
                ADD(m, m + ny * nz,m + ny * nz + nz, m + nz);
                FI[I3(i,j,k)] = mo.mFacets.size() - 1;
            }
        }
    }
    for(i = 0;i < nx - 1;i++) {
        for(j = 1;j < ny - 1;j++) {
            for(k = 0;k < nz - 1;k++) {
                ADD(m,m + ny * nz,m + ny * nz + 1,m + 1);
                FI[I2(i,j,k)] = mo.mFacets.size() - 1;
            }
        }
    }
    for(i = 1;i < nx - 1;i++) {
        for(j = 0;j < ny - 1;j++) {
            for(k = 0;k < nz - 1;k++) {
                ADD(m,m + nz,m + nz + 1,m + 1);
                FI[I1(i,j,k)] = mo.mFacets.size() - 1;
            }
        }
    }
    mo.mNF = mo.mFacets.size();
    /*boundaries*/
    for(k = 0;k < nz; k += (nz - 1)) {
        Patch p;
        p.from = mo.mFacets.size();
        for(i = 0;i < nx - 1;i++) {
            for(j = 0;j < ny - 1;j++) {
                ADD(m, m + ny * nz,m + ny * nz + nz, m + nz);
                FI[I3(i,j,k)] = mo.mFacets.size() - 1;
            }
        }
        p.to = mo.mFacets.size();
        mo.mPatches.push_back(p);
    }
    for(j = 0;j < ny;j += (ny - 1)) {
        Patch p;
        p.from = mo.mFacets.size();
        for(i = 0;i < nx - 1;i++) {
            for(k = 0;k < nz - 1;k++) {
                ADD(m,m + ny * nz,m + ny * nz + 1,m + 1);
                FI[I2(i,j,k)] = mo.mFacets.size() - 1;
            }
        }
        p.to = mo.mFacets.size();
        mo.mPatches.push_back(p);
    }
    for(i = 0;i < nx; i += (nx - 1)) {
        Patch p;
        p.from = mo.mFacets.size();
        for(j = 0;j < ny - 1;j++) {
            for(k = 0;k < nz - 1;k++) {
                ADD(m,m + nz,m + nz + 1,m + 1);
                FI[I1(i,j,k)] = mo.mFacets.size() - 1;
            }
        }
        p.to = mo.mFacets.size();
        mo.mPatches.push_back(p);
    }
/*compute normals of mPatches*/
#define NORMAL(i,j,k,l,p) {                     \
    p.N = ((vp[j] - vp[i]) ^ (vp[k] - vp[i]));  \
    p.N /= mag(p.N);                            \
    p.C = (vp[i] + vp[j] + vp[k] + vp[l]) / 4;  \
}
    NORMAL(0,1,2,3,mo.mPatches[0]);
    NORMAL(4,5,6,7,mo.mPatches[1]);
    NORMAL(0,1,5,4,mo.mPatches[2]);
    NORMAL(3,2,6,7,mo.mPatches[3]);
    NORMAL(0,3,7,4,mo.mPatches[4]);
    NORMAL(1,2,6,5,mo.mPatches[5]);
#undef NORMAL

/*end*/
#undef ADD

    /*cells*/
    for(i = 0;i < nx - 1;i++) {
        for(j = 0;j < ny - 1;j++) {
            for(k = 0;k < nz - 1;k++) {
                Cell c;
                m = I3(i,j,k);
                c.push_back(FI[m]);
                c.push_back(FI[m + 1]);
    
                m = I2(i,j,k);
                c.push_back(FI[m]);
                c.push_back(FI[m + (nz - 1)]);
    
                m = I1(i,j,k);
                c.push_back(FI[m]);
                c.push_back(FI[m + (ny - 1) * (nz - 1)]);
    
                mo.mCells.push_back(c);
            }
        }
    }
    mo.mBCS = mo.mCells.size();
#undef I0
#undef I1
#undef I2
#undef I3
    /*remove duplicates*/
    int deformed = 0;
    for(i = 0;i < 8;i++) {
        for(j = i + 1;j < 8;j++) {
            if(equal(vp[i],vp[j])) {
                deformed = 1;
                break;
            }
        }
    }
    if(deformed)
        remove_duplicate(mo);
    /*end*/
}

/**
  Remove duplicate vertices,faces and cells
 */
void remove_duplicate(Mesh::MeshObject& mo) {
    Int i,j,sz,corr;
    int count;
    /*vertices*/
    sz = mo.mVertices.size();
    corr = 0;
    std::vector<int> dup(sz,0);
    for(i = 0;i < sz;i++) {
        for(j = sz - 1;j >= i + 1;j--) {
            if(equal(mo.mVertices[i],mo.mVertices[j])) {
                dup[i] = -int(j);
                if(i < mo.mNV) corr++;
                break;
            }
        }
    }
    mo.mNV -= corr;
    //remove duplicate vertices
    {
        Vertices vt(mo.mVertices.begin(), mo.mVertices.end());
        mo.mVertices.clear();
        count = 0;
        for(i = 0;i < sz;i++) {
            if(!dup[i]) {
                mo.mVertices.push_back(vt[i]);
                dup[i] = count++;
            }
        }
        for(i = 0;i < sz;i++) {
            if(dup[i] < 0) 
                dup[i] = dup[-dup[i]];
        }
    }
    /*faces*/
    sz = mo.mFacets.size();
    for(i = 0;i < sz;i++) {
        Facet& f = mo.mFacets[i];
        forEach(f,j)
            f[j] = dup[f[j]];
    }
    dup.clear();
    dup.assign(sz,0);
    count = 0;
    corr = 0;
    for(i = 0;i < sz;i++) {
        Facet& f = mo.mFacets[i];
        forEach(f,j) {
            forEachS(f,k,j+1) {
                if(f[j] == f[k]) {
                    f.erase(f.begin() + k);
                    k--;
                }
            }
        }
        if(f.size() < 3) {
            dup[i] = -1;
            if(i < mo.mNF) corr++;
        } else {
            dup[i] = count;
            count++;
        }
    }
    mo.mNF -= corr;
    //remove deformed faces
    {
        Facets ft(mo.mFacets.begin(), mo.mFacets.end());
        mo.mFacets.clear();
        for(i = 0;i < sz;i++) {
            if(dup[i] >= 0) mo.mFacets.push_back(ft[i]);
        }
    }
    //adjust bstart
    forEach(mo.mPatches,i) {
        Patch& pi = mo.mPatches[i];
        pi.from = dup[pi.from];
        pi.to = dup[pi.to];
    }
    /*cells*/
    sz = mo.mCells.size();
    for(i = 0;i < sz;i++) {
        Cell& c = mo.mCells[i];
        forEach(c,j) {
            if(dup[c[j]] < 0) {
                c.erase(c.begin() + j);
                j--;
            } else
                c[j] = dup[c[j]];
        }
    }
}

#define MAXNUM 1073741824

/**
  Merge mesh m2 onto m1 (internal) and b (boundary) meshes
 */
void merge(MeshObject& m1,MergeObject& b,MeshObject& m2) {
    Int found,s0,s1,s2,s3;

    //vertices
    {
        s0 = m1.mVertices.size();
        s1 = m2.mNV;
        s2 = m2.mVertices.size();
        s3 = b.vb.size();
        m1.mVertices.insert(m1.mVertices.end(),m2.mVertices.begin(),m2.mVertices.begin() + s1);

        IntVector locv(s2 - s1,MAXNUM);
        for(Int i = s1;i < s2;i++) {
            found = 0;
            for(Int j = 0;j < s3;j++) {
                if(equal(m2.mVertices[i],b.vb[j])) {
                    locv[i - s1] += j;
                    found = 1;
                    break;
                }
            }
            if(!found) {
                b.vb.push_back(m2.mVertices[i]);
                locv[i - s1] += b.vb.size() - 1;
            }
        }
        forEach(m2.mFacets,i) {
            Facet& ft = m2.mFacets[i];
            forEach(ft,j) {
                if(ft[j] >= s1) {
                    ft[j] = locv[ft[j] - s1];
                } else {
                    ft[j] += s0;
                }
            }
        }
    }
    //faces
    {
        s0 = m1.mFacets.size();
        s1 = m2.mNF;
        s2 = m2.mFacets.size();
        s3 = b.fb.size();
        m1.mFacets.insert(m1.mFacets.end(),m2.mFacets.begin(),m2.mFacets.begin() + s1);

        //insert faces
        IntVector index0(s3,0),index1(s2 - s1,0);
        Int count = 0;
        b.fb.reserve(s3 + s2 - s1);
        for(Int j = 0;j < s3;j++) {
            found = 0;
            for(Int i = s1;i < s2;i++) {
                if(!index1[i - s1] && equalSet(m2.mFacets[i],b.fb[j])) {

                    m1.mFacets.push_back(b.fb[j]);
                    index0[j]      = m1.mFacets.size() - 1;
                    index1[i - s1] = m1.mFacets.size() - 1;

                    found = 1;
                    break;
                }
            }
            if(!found) {
                index0[j] = MAXNUM + count;
                b.fb[count] = b.fb[j];
                count++;
            }
        }
        for(Int i = s1;i < s2;i++) {
            if(!index1[i - s1]) {
                index1[i - s1] = MAXNUM + count;

                if(count >= s3) b.fb.push_back(m2.mFacets[i]);
                else b.fb[count] = m2.mFacets[i];
                count++;
            }
        }
        b.fb.resize(count);

        //insert patch
        forEach(m2.mPatches,i) {
            m2.mPatches[i].from += s0 + s3;
            m2.mPatches[i].to += s0 + s3;
        }
        forEach(m1.mPatches,i) {
            m1.mPatches[i].from += s1;
            m1.mPatches[i].to += s1;
        }
        Int npatch = m1.mPatches.size();
        forEach(m2.mPatches,i) {
            Patch& p = m2.mPatches[i];
            bool skip = false;
            Int j = 0;
            for(;j < npatch;j++) {
                if(equal(p.C,m1.mPatches[j].C)) {
                    skip = true;
                    break;
                }
            }
            if(!skip)
                m1.mPatches.push_back(p);
            else {
                Int s4 = m1.mPatches[j].to - m1.mPatches[j].from;
                for(Int k = 0;k < j;k++) {
                    m1.mPatches[k].from += s4;
                    m1.mPatches[k].to += s4;
                }
                for(Int k = i;k < m2.mPatches.size();k++) {
                    m2.mPatches[k].from -= s4;
                    m2.mPatches[k].to -= s4;
                }
                m1.mPatches.erase(m1.mPatches.begin() + j);
            }
        }
        //adjust face ids in cells
        forEach(m1.mCells,i) {
            Cell& ct = m1.mCells[i];
            forEach(ct,j) {
                if(ct[j] >= MAXNUM) {
                    ct[j] = index0[ct[j] - MAXNUM];
                }
            }
        }
        forEach(m2.mCells,i) {
            Cell& ct = m2.mCells[i];
            forEach(ct,j) {
                if(ct[j] >= s1) {
                    ct[j] = index1[ct[j] - s1];
                } else {
                    ct[j] += s0;
                }
            }
        }
    }
    //cells
    {
        m1.mCells.insert(m1.mCells.end(),m2.mCells.begin(),m2.mCells.end());
    }
}
/**
  Merge boundary and internals
 */
void merge(Mesh::MeshObject& m,MergeObject& b) {
    m.mNV = m.mVertices.size();
    m.mNF = m.mFacets.size();
    m.mBCS = m.mCells.size();

    m.mVertices.insert(m.mVertices.end(),b.vb.begin(),b.vb.end());
    m.mFacets.insert(m.mFacets.end(),b.fb.begin(),b.fb.end());
    forEach(m.mFacets,i) {
        Facet& ft = m.mFacets[i];
        forEach(ft,j) {
            if(ft[j] >= MAXNUM) {
                ft[j] -= MAXNUM;
                ft[j] += m.mNV;
            }
        }
    }
    forEach(m.mCells,i) {
        Cell& ct = m.mCells[i];
        forEach(ct,j) {
            if(ct[j] >= MAXNUM) {
                ct[j] -= MAXNUM;
                ct[j] += m.mNF;
            }
        }
    }
}

#undef MAXNUM
