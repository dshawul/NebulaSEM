#include "field.h"

/**
  Nodal Discontinuos Galerkin variables
 */
namespace DG {
    Int Nop[3] = {0, 0, 0};
    Int NPX, NPY, NPZ;
    Int NP, NPI, NPMAT, NPF;
    Scalar Penalty;

    Scalar **psi;
    Scalar **dpsi;
    Scalar **xgl;
    Scalar **wgl;
    TensorCellField Jinv(false);
}
/** 
  Compute the orthogonal legendre polynomial and its first & second derivatives, given p and x
  They are computed using the recurrence relation
    phi_0(x) = 1 -> 0th degree legendre polynomial
    phi_1(x) = x -> 1st degree ...
    phi_n(x) = a * x * phi_{n-1}(x)  - b * phi_{n-2}(x)
  where
    a = (2n - 1) / n
    b = (n - 1) /n
  The first and second derivative can be obtained similarily applying product rule
  of derivaties on the recurrence relation.
 */
void DG::legendre(int p, Scalar x, Scalar& L0, Scalar& L0_1, Scalar& L0_2) {
    Scalar a,b;                
    Scalar L2, L2_1, L2_2;
    Scalar L1, L1_1, L1_2;
    L1 = 0; L1_1 = 0; L1_2 = 0;
    L0 = 1; L0_1 = 0; L0_2 = 0;

    for(int i = 1; i <= p; i++) {
        L2 = L1; L2_1 = L1_1; L2_2 = L1_2;
        L1 = L0; L1_1 = L0_1; L1_2 = L0_2;
        a = (2 * i - 1.0) / i;
        b = (i - 1.0) / i;
        L0 = a * x * L1 - b * L2;
        L0_1 = a * (L1 + x * L1_1) - b * L2_1;
        L0_2 = a * (2 * L1_1 + x * L1_2) - b * L2_2;
    }
}
/**
  Compute Legendre-Gauss-Lobato interpolation points and weights
 */
void DG::legendre_gauss_lobatto(int N, Scalar* xgl, Scalar* wgl) {
    Scalar L0,L0_1,L0_2;
    int p = N - 1; 
    int ph = N / 2;
    Scalar x,dx;

    //quick exit for 0th order polynomial
    if(N == 1) {
        xgl[0] = 0;
        wgl[0] = 2;
        return;
    }

    //compute first half of roots
    for(int i = 0; i < ph; i++) {

        //use roots of Chebyshev polynomial as initial guess for LGL points
        x = cos((2 * i + 1) * Constants::PI / (2 * N));

        //use Newton's method to compute LGL points by computing the roots of (1 - x^2)P'(x)
        //where P'(x) is the derivative of the legendre polynomial.
        for(int k = 1; k <= 20; k++) {
            legendre(p,x,L0,L0_1,L0_2);
            dx = -(1 - x * x) * L0_1 / (-2 * x * L0_1 + (1 - x * x) * L0_2);
            x += dx;
            if(fabs(dx) < 1.0e-20) 
                break;
        }

        //assign interpolation/integration point and associated weight
        xgl[p - i] = x;
        wgl[p - i] = 2 / (p * (p + 1) * L0 * L0);
    }

    //Add 0th point for odd N
    if (N != 2*ph) {
        x = 0;
        legendre(p,x,L0,L0_1,L0_2);
        xgl[ph] = x;
        wgl[ph] = 2 / (p * (p + 1) * L0 * L0);
    }

    //mirror second half of the roots
    for(int i = 0; i < ph; i++) {
        xgl[i] = -xgl[p - i];
        wgl[i] = +wgl[p - i];
    }
}
/**
  Compute lagrange basis functions
 */
void DG::lagrange_basis(int v, int N, Scalar* xgl, Scalar* psi) {
    for(int i = 0;i < N;i++) {
        if(i != v) psi[i] = 0;
        else psi[i] = 1;
    }
}
/**
  Compute lagrange basis derivatives at given point x, given interpolation points
 */
void DG::lagrange_basis_derivative(int v, int N, Scalar* xgl, Scalar* dpsi) {
    Scalar xi,xj,xk,prod;
    Scalar x = xgl[v];
    for(int i = 0;i < N;i++) {
        dpsi[i] = 0;
        xi = xgl[i];
        for(int j = 0;j < N;j++) {
            if(i != j) {
                xj = xgl[j];
                prod = 1;
                for(int k = 0;k < N;k++) {
                    xk = xgl[k];
                    if(k != i && k != j)
                        prod *= ((x - xk) / (xi - xk));
                }
                dpsi[i] += prod / (xi - xj);
            }
        }
    }
}
/**
  Initialize polynomial order
 */
void DG::init_poly() {
    NPX = Nop[0] + 1;
    NPY = Nop[1] + 1;
    NPZ = Nop[2] + 1;
    NP  = NPX * NPY * NPZ;
    NPI = NPX + NPY + NPZ;
    if(NP > 1) 
        NPMAT = NP * NPI;
    else 
        NPMAT = 0;
    if(NPX <= NPY && NPX <= NPZ)
        NPF = NPY * NPZ;
    else if(NPY <= NPX && NPY <= NPZ)
        NPF = NPX * NPZ;
    else
        NPF = NPX * NPY;
}
Int find_surface_node(Int ci, Int cj, Int index0) {
    using namespace Mesh;
    using namespace DG;
    Int face, index1 = Int(-1), minindex;
    Scalar mindist = 1e30;

    for(face = 0; face < 2; face++) {
        Int kk = (face == 0) ? 0 : (NPZ - 1);
        forEachLglXY(ii,jj) {
            index1 = INDEX4(cj,ii,jj,kk);
            if(equal(cC[index0],cC[index1]))
                return index1;
            else {
                Scalar dist = mag(cC[index0] - cC[index1]);
                if(dist < mindist) {
                    mindist = dist;
                    minindex = index1;
                }
            }
        }
    }
    for(face = 2; face < 4; face++) {
        Int jj = (face == 2) ? 0 : (NPY - 1);
        forEachLglXZ(ii,kk) {
            index1 = INDEX4(cj,ii,jj,kk);
            if(equal(cC[index0],cC[index1]))
                return index1;
            else {
                Scalar dist = mag(cC[index0] - cC[index1]);
                if(dist < mindist) {
                    mindist = dist;
                    minindex = index1;
                }
            }
        }
    }
    for(face = 4; face < 6; face++) {
        Int ii = (face == 4) ? 0 : (NPX - 1);
        forEachLglYZ(jj,kk) {
            index1 = INDEX4(cj,ii,jj,kk);
            if(equal(cC[index0],cC[index1]))
                return index1;
            else {
                Scalar dist = mag(cC[index0] - cC[index1]);
                if(dist < mindist) {
                    mindist = dist;
                    minindex = index1;
                }
            }
        }
    }

#if 0
    using namespace Util;
    std::cout << "Not Found: " << cC[index0] << std::endl;
    std::cout << "cell " << ci << " = " << gCells[ci] << std::endl;
    std::cout << "cell " << cj << " = " << gCells[cj] << std::endl;
    forEachLgl(i,j,k) {
        Int index = INDEX4(ci,i,j,k);
        std::cout << ci << ". " << index << " " << cC[index]
                  << " mag " << mag(cC[index]) << std::endl;
    }
    forEachLgl(i,j,k) {
        Int index = INDEX4(cj,i,j,k);
        std::cout << cj << ". " << index << " " << cC[index]
                  << " mag " << mag(cC[index]) << " equal "
                  << " equal " << equal(cC[index],cC[index0])
                  << " mag(diff) " << mag(cC[index]-cC[index0])/mag(cC[index])
                  << std::endl;
    }
    std::cout << "Returning minindex " << minindex << std::endl;
#endif

    return minindex;
}
/**
  Initialize geometry
 */
void DG::init_geom() {
    using namespace Mesh;
    using namespace Constants;

    //compute coordinates of nodes via transfinite interpolation
    FO = gALLfield;
    FN = gALLfield;
    for(Int ci = 0; ci < gBCS;ci++) {
        Cell& c = gCells[ci];
        Facet& f1 = gFacets[c[0]];
        Facet& f2 = gFacets[c[1]];
        //vertices
        Vertex vp[8];
        {
            Vertex vp1[8];
            forEach(f1,i)
                vp1[i + 0] = gVertices[f1[i]];
            forEach(f2,i)
                vp1[i + 4] = gVertices[f2[i]];  
            Int id = gFaceID[ci][0];
            if(id == 2) {
                Int order[8] = {0,1,5,4,3,2,6,7};
                for(Int i = 0;i < 8;i++) 
                    vp[order[i]] = vp1[i];
            } else if(id == 4) {
                Int order[8] = {0,3,7,4,1,2,6,5};
                for(Int i = 0;i < 8;i++) 
                    vp[order[i]] = vp1[i];
            } else {
                for(Int i = 0;i < 8;i++) 
                    vp[i] = vp1[i];
            }
        }   
        //edges
        static const int sides[12][2] = {
            {0,1}, {3,2}, {7,6}, {4,5},
            {0,3}, {1,2}, {5,6}, {4,7},
            {0,4}, {1,5}, {2,6}, {3,7}
        };
        Vertex ev[12][2];
        for(Int i = 0;i < 12;i++) {
            ev[i][0] = vp[sides[i][0]];
            ev[i][1] = vp[sides[i][1]];
        }

        //coordinates
        Vertex v,vd[12],vf[6];
        Scalar rx,ry,rz;

#define ADDV(w,m,ev,vd) {                               \
    vd[w] = (1 - m) * ev[w][0] + (m) * ev[w][1];        \
    if(Controls::is_spherical)                          \
        vd[w] = vd[w] * (((1 -m) * mag(ev[w][0]) +      \
                            (m)  * mag(ev[w][1]))       \
                        / mag(vd[w]));                  \
}

#define ADDF(w,rr,rs,i00,i01,i10,i11,ir0,ir1,i0s,i1s) { \
    vf[w] = Interpolate_face(                           \
            rr,rs,                                      \
            vp[i00],vp[i01],vp[i10],vp[i11],            \
            vd[ir0],vd[ir1],vd[i0s],vd[i1s]);           \
    if(Controls::is_spherical)                          \
        vf[w] = (mag(vd[ir0])/mag(vf[w])) * vf[w];      \
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
    if(Controls::is_spherical)                          \
        v = (mag(vd[8])/mag(v)) * v;                    \
}

#define ADD() {                                         \
    rx = (xgl[0][i] + 1) / 2;                           \
    ry = (xgl[1][j] + 1) / 2;                           \
    rz = (xgl[2][k] + 1) / 2;                           \
    ADDV(0 ,rx,ev,vd);                                  \
    ADDV(1 ,rx,ev,vd);                                  \
    ADDV(2 ,rx,ev,vd);                                  \
    ADDV(3 ,rx,ev,vd);                                  \
    ADDV(4 ,ry,ev,vd);                                  \
    ADDV(5 ,ry,ev,vd);                                  \
    ADDV(6 ,ry,ev,vd);                                  \
    ADDV(7 ,ry,ev,vd);                                  \
    ADDV(8 ,rz,ev,vd);                                  \
    ADDV(9 ,rz,ev,vd);                                  \
    ADDV(10,rz,ev,vd);                                  \
    ADDV(11,rz,ev,vd);                                  \
    ADDF(0, rx,ry, 0,3,1,2, 0,1,4,5);                   \
    ADDF(1, rx,ry, 4,7,5,6, 3,2,7,6);                   \
    ADDF(2, rx,rz, 0,4,1,5, 0,3,8,9);                   \
    ADDF(3, rx,rz, 3,7,2,6, 1,2,11,10);                 \
    ADDF(4, ry,rz, 0,4,3,7, 4,7,8,11);                  \
    ADDF(5, ry,rz, 1,5,2,6, 5,6,9,10);                  \
    ADDC();                                             \
};

        forEachLgl(i,j,k) {
            ADD();

            Scalar wgt = wgl[0][i] * wgl[1][j] * wgl[2][k] / 8;
            Int index = INDEX4(ci,i,j,k);
            cC[index] = v;
            cV[index] *= wgt;
        }
    }   
#undef ADDV
#undef ADDF
#undef ADDC
#undef ADD

    //compute face properties of control volumes formed by LGL nodes
    for(Int ci = 0; ci < gBCS;ci++) {
        Cell& c = gCells[ci];
        forEach(c,mm) {
            Int face = gFaceID[ci][mm];
            Int fi = c[mm];
            Int cj = gFNC[fi];
            if(cj == ci) 
                continue;

#define ADD() {                                             \
    FO[indf] = index0;                                      \
    FN[indf] = index1;                                      \
    if(index1 >= gBCSfield) {                               \
        Vector dir = unit(fN[indf]);                        \
        Scalar d = dot(fC[indf] - cC[index0],dir);          \
        cC[index1] = cC[index0] + d * dir;                  \
        cV[index1] = cV[index0];                            \
    }                                                       \
    Scalar d;                                               \
    if(equal(cC[index0],cC[index1])) d = 0.5;               \
    else d = dot(fC[indf] - cC[index0],fN[indf]) /          \
             dot(cC[index1] - cC[index0],fN[indf]);         \
    fC[indf] = cC[index0] + d * (cC[index1] - cC[index0]);  \
    fN[indf] *= wgt;                                        \
}

            if(face == 0 || face == 1) {
                Int k = (face == 0) ? 0 : (NPZ - 1);
                forEachLglXY(i,j) {
                    Scalar wgt = wgl[0][i] * wgl[1][j] / 4;
                    Int indf = fi * NPF + i * NPY + j;
                    Int index0 = INDEX4(ci,i,j,k);
                    Int index1 = INDEX4(cj,i,j,(NPZ - 1) - k);
                    if(index1 < gBCSfield) index1 = find_surface_node(ci, cj, index0);
                    ADD();
                }
            } else if(face == 2 || face == 3) {
                Int j = (face == 2) ? 0 : (NPY - 1);
                forEachLglXZ(i,k) {
                    Scalar wgt = wgl[0][i] * wgl[2][k] / 4;
                    Int indf = fi * NPF + i * NPZ + k;
                    Int index0 = INDEX4(ci,i,j,k);
                    Int index1 = INDEX4(cj,i,(NPY - 1) - j,k);
                    if(index1 < gBCSfield) index1 = find_surface_node(ci, cj, index0);
                    ADD();
                }
            } else {
                Int i = (face == 4) ? 0 : (NPX - 1);
                forEachLglYZ(j,k) {
                    Scalar wgt = wgl[1][j] * wgl[2][k] / 4;
                    Int indf = fi * NPF + j * NPZ + k;
                    Int index0 = INDEX4(ci,i,j,k);
                    Int index1 = INDEX4(cj,(NPX - 1) - i,j,k);
                    if(index1 < gBCSfield) index1 = find_surface_node(ci, cj, index0);
                    ADD();
                }
            }

#undef ADD
        }
    }
}
/**
  Initialize basis functions
 */
void DG::init_basis() {
    using namespace Mesh;
    using namespace Constants;

    //directional LGL
    xgl = new Scalar*[3];
    wgl = new Scalar*[3];
    psi = new Scalar*[3];
    dpsi = new Scalar*[3];
    for(Int i = 0;i < 3;i++) {
        Int ngl = Nop[i] + 1;
        xgl[i] = new Scalar[ngl];
        wgl[i] = new Scalar[ngl];
        psi[i] = new Scalar[ngl*ngl];
        dpsi[i] = new Scalar[ngl*ngl];
        legendre_gauss_lobatto(ngl,xgl[i],wgl[i]);
        for(Int j = 0;j < ngl;j++) {
            lagrange_basis(j,ngl,xgl[i],&psi[i][j*ngl]);
            lagrange_basis_derivative(j,ngl,xgl[i],&dpsi[i][j*ngl]);
        }
    }

    //init geometry
    init_geom();

    //Compute Jacobian matrix
    Jinv.deallocate(false);
    Jinv.construct();
    #pragma omp parallel for
    #pragma acc parallel loop copyin(gBCS)
    for(Int ci = 0; ci < gBCS;ci++) {

        forEachLgl(ii,jj,kk) {
            Int index = INDEX4(ci,ii,jj,kk);
            Tensor Ji = Tensor(0);
#define JACD(im,jm,km) {                                    \
    Int index1 = INDEX4(ci,im,jm,km);                       \
    Vector dpsi_ij;                                         \
    DPSI(dpsi_ij,im,jm,km);                                 \
    Ji += mul(cC[index1], dpsi_ij);                         \
}
            forEachLglX(i) JACD(i,jj,kk);
            forEachLglY(j) if(j != jj) JACD(ii,j,kk);
            forEachLglZ(k) if(k != kk) JACD(ii,jj,k);
#undef JACD
            //
            // Part I
            // ======
            // Note that the Jacobian is stored inverted and transposed so that it
            // can be used directly for gradient/divergence calculations.
            //
            //   {x,y,z} = J * {a,b,c} where J is
            //
            //       [dx/da dx/db dx/dc]
            //   J = [dy/da dy/db dy/dc]
            //       [dz/da dz/db dz/dc]
            //
            //   Also, {a,b,c} = J' * {x,y,z} where J' is the inverse of J
            //
            //       [da/dx da/dy da/dz]
            //  J' = [db/dx db/dy db/dz]
            //       [dc/dx dc/dy dc/dz]
            //
            //   The transpose of J' is used to compute derivatives of q
            //
            //   {dq/dx, dq/dy, dq/dz} = trn(J') * {dq/da, dq/db, dq/dc}
            //
            //   While the transpose of J is used for the inverse
            //
            //   {dq/da, dq/db, dq/dc} = trn(J) * {dq/dx, dq/dy, dq/dz}
            //
            // Part II
            // =======
            // For 2D problems aligned with one of the axes or is inclined,
            // the Jacobian is singular. We will use the pseudo-inverse (Moore-Penrose)
            // inverse to solve a 2D problem in a 3D space
            //
            Tensor JiT = trn(Ji);
            Ji = mul(JiT,Ji);
            if(NPX == 1) Ji[XX] = 1;
            if(NPY == 1) Ji[YY] = 1;
            if(NPZ == 1) Ji[ZZ] = 1;
            Ji = inv(Ji); /*invert*/
            if(NPX == 1) Ji[XX] = 0;
            if(NPY == 1) Ji[YY] = 0;
            if(NPZ == 1) Ji[ZZ] = 0;
            Ji = mul(Ji,JiT);
            Ji = trn(Ji); /*transpose*/
            Jinv[index] = Ji;
        }
    }
}

