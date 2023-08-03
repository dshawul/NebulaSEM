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
    Scalar **psiRef;
    Scalar **psiCor;
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
  Compute Legendre-Gauss-Lobato interpolation points and associated
  quadrature weights
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
  Compute lagrange basis values for a set of points
 */
void DG::lagrange_basis(int N, Scalar* xgl, int Ns, Scalar* xs, Scalar* psi) {
    Scalar x,xj,xk,prod;
    for(int s = 0;s < Ns;s++) {
        x = xs[s];
        for(int j = 0;j < N;j++) {
            xj = xgl[j];
            prod = 1;
            for(int k = 0;k < N;k++) {
                xk = xgl[k];
                if(k != j)
                    prod *= ((x - xk) / (xj - xk));
            }
            psi[j*Ns+s] = prod;
        }
    }
}
/**
  Compute lagrange basis derivatives for a set of points
 */
void DG::lagrange_basis_derivative(int N, Scalar* xgl, int Ns, Scalar* xs, Scalar* dpsi) {
    Scalar xi,xj,xk,prod;
    for(int s = 0; s < Ns; s++) {
        Scalar x = xs[s];
        for(int i = 0;i < N;i++) {
            dpsi[s*N+i] = 0;
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
                    dpsi[s*N+i] += prod / (xi - xj);
                }
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
/**
  Initialize geometry
 */
void DG::init_geom() {
    using namespace Mesh;
    using namespace Constants;

    static const Int face_map[] = {0, NPZ - 1, 0, NPY - 1, 0, NPX - 1};

    //compute coordinates of nodes via transfinite interpolation
    FO = gALLfield;
    FN = gALLfield;
    for(Int ci = 0; ci < gBCS;ci++) {
        Cell& c = gCells[ci];
        Cell& fids = gFaceID[ci];

        //find first and second face
        Facet f1, f2;
        Int id0 = fids[0], id1 = (id0 ^ 1);
        forEach(c,i) {
            Facet& f = gFacets[c[i]];
            Int id = fids[i];
            Facet fn;
            if(id == id0) {
                if(!f1.size())
                    f1 = f;
                else {
                    gMesh.mergeFacets(f1,f,fn);
                    f1 = fn;
                }
            } else if(id == id1) {
                if(!f2.size())
                    f2 = f;
                else {
                    gMesh.mergeFacets(f2,f,fn);
                    f2 = fn;
                }
            }
        }
#if 0
        {
            using namespace Util;
            std::cout << "=========== cell " << ci << " " << c << " ===========" << std::endl; 
            std::cout << f1 << std::endl;
            forEach(f1,i)
                std::cout << i << ". " << gVertices[f1[i]] << std::endl;
            std::cout << f2 << std::endl;
            forEach(f2,i)
                std::cout << i << ". " << gVertices[f2[i]] << std::endl;
            std::cout << gFaceID[ci] << std::endl;
        }
#endif
        //vertices
        Vertex vp[8];
        {
            //find the 8 corners of hex cells
            Vertex vp1[8];
            IntVector vpi;
            gMesh.getHexCorners(f1,f2,vpi);
            forEach(vpi,i)
                vp1[i] = gVertices[vpi[i]];

            //reorder corners
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
    if(is_spherical)                                    \
        vd[w] = vd[w] * (((1 -m) * mag(ev[w][0]) +      \
                            (m)  * mag(ev[w][1]))       \
                        / mag(vd[w]));                  \
}

#define ADDF(w,rr,rs,i00,i01,i10,i11,ir0,ir1,i0s,i1s) { \
    vf[w] = Interpolate_face(                           \
            rr,rs,                                      \
            vp[i00],vp[i01],vp[i10],vp[i11],            \
            vd[ir0],vd[ir1],vd[i0s],vd[i1s]);           \
    if(is_spherical)                                    \
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
    if(is_spherical)                                    \
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

#undef ADDV
#undef ADDF
#undef ADDC
#undef ADD
    }

    //compute face properties of control volumes formed by LGL nodes
    for(Int ci = 0; ci < gBCS;ci++) {
        Cell& c = gCells[ci];
        IntVector& face_ids = gFaceID[ci];

        forEach(c,mm) {
            Int face_o = face_ids[mm];
            Int fi = c[mm];
            Int cj = gFNC[fi];
            Int fm = gFMC[fi];
            if(cj == ci) 
                continue;

            Int face_n = (face_o ^ 1);
            if(cj < gBCS) {
                Cell& cn = gCells[cj];
                forEach(cn,i) {
                    if(cn[i] == fi) {
                        face_n = gFaceID[cj][i];
                        break;
                    }
                }
            }

#define ADD() {                                             \
    FO[indf] = index0;                                      \
    FN[indf] = index1;                                      \
    if(index1 >= gBCSfield) {                               \
        cC[index1] = cC[index0];                            \
        cV[index1] = cV[index0];                            \
    }                                                       \
    fC[indf] = (fm <= 1) ? cC[index0] : cC[index1];         \
    fN[indf] *= wgt;                                        \
}

            Int vo = face_map[face_o];
            Int vn = face_map[face_n];
            if(face_o == 0 || face_o == 1) {
                forEachLglXY(i,j) {
                    Scalar wgt = wgl[0][i] * wgl[1][j] / 4;
                    Int indf = fi * NPF + i * NPY + j;
                    Int index0 = INDEX4(ci,i,j,vo);
                    Int index1;
                    if(face_n == 0 || face_n == 1)
                        index1 = INDEX4(cj,i,j,vn);
                    else if(face_n == 2 || face_n == 3)
                        index1 = INDEX4(cj,i,vn,j);
                    else
                        index1 = INDEX4(cj,vn,i,j);
                    ADD();
                }
            } else if(face_o == 2 || face_o == 3) {
                forEachLglXZ(i,k) {
                    Scalar wgt = wgl[0][i] * wgl[2][k] / 4;
                    Int indf = fi * NPF + i * NPZ + k;
                    Int index0 = INDEX4(ci,i,vo,k);
                    Int index1;
                    if(face_n == 0 || face_n == 1)
                        index1 = INDEX4(cj,i,k,vn);
                    else if(face_n == 2 || face_n == 3)
                        index1 = INDEX4(cj,i,vn,k);
                    else
                        index1 = INDEX4(cj,vn,i,k);
                    ADD();
                }
            } else {
                forEachLglYZ(j,k) {
                    Scalar wgt = wgl[1][j] * wgl[2][k] / 4;
                    Int indf = fi * NPF + j * NPZ + k;
                    Int index0 = INDEX4(ci,vo,j,k);
                    Int index1;
                    if(face_n == 0 || face_n == 1)
                        index1 = INDEX4(cj,j,k,vn);
                    else if(face_n == 2 || face_n == 3)
                        index1 = INDEX4(cj,j,vn,k);
                    else
                        index1 = INDEX4(cj,vn,j,k);
                    ADD();
                }
            }

#undef ADD
        }
    }

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
/**
  Initialize basis functions
 */
void DG::init_basis() {
    //directional LGL
    xgl = new Scalar*[3];
    wgl = new Scalar*[3];
    psi = new Scalar*[3];
    dpsi = new Scalar*[3];
    psiRef = new Scalar*[3*2];
    psiCor = new Scalar*[3*2];
    for(Int i = 0;i < 3;i++) {
        Int ngl = Nop[i] + 1;
        xgl[i] = new Scalar[ngl];
        wgl[i] = new Scalar[ngl];
        psi[i] = new Scalar[ngl*ngl];
        dpsi[i] = new Scalar[ngl*ngl];
        for(Int j = 0; j < 2; j++) {
            psiRef[i*2+j] = new Scalar[ngl*ngl];
            psiCor[i*2+j] = new Scalar[ngl*ngl];
        }
        //initialize interpolation at LGL nodes
        legendre_gauss_lobatto(ngl,xgl[i],wgl[i]);

        lagrange_basis(ngl,xgl[i],ngl,xgl[i],psi[i]);
        lagrange_basis_derivative(ngl,xgl[i],ngl,xgl[i],dpsi[i]);

        //initialize interpolation for 2:1 amr
        Scalar* xglRef = new Scalar[2*ngl];
        if(ngl == 1) {
            xglRef[0] = 0;
            xglRef[1] = 0;
        } else {
            for(Int j = 0;j < ngl;j++) {
                xglRef[j] = -0.5 + xgl[i][j]/2;
                xglRef[j+ngl] = 0.5 + xgl[i][j]/2;
            }
        }
        for(Int j = 0; j < 2; j++) {
            lagrange_basis(ngl,xgl[i],ngl,&xglRef[j*ngl],psiRef[i*2+j]);
            inv(psiRef[i*2+j],psiCor[i*2+j],ngl);
            //
            // After inverting refine matrix to get coarsening matrix,
            // zero out coefficients that do extrapolation instead of interpolation
            // This is a compact interpolation scheme where each child interpolates
            // points in its local region.
            //
            Int nglh = ngl / 2;
            if(ngl & 1) {
                if(j == 0) {
                    for(Int l = 0; l < ngl; l++) {
                        for(Int k = nglh + 1; k < ngl; k++)
                            psiCor[i*2+j][l * ngl + k] = 0;
                    }
                } else {
                    for(Int l = 0; l < ngl; l++) {
                        for(Int k = 0; k < nglh; k++)
                            psiCor[i*2+j][l * ngl + k] = 0;
                    }
                }
                //Double values except at midpoint
                if(ngl > 1) {
                    for(Int l = 0; l < ngl; l++) {
                        for(Int k = 0; k < ngl; k++)
                            if(k != nglh)
                                psiCor[i*2+j][l * ngl + k] *= 2;
                    }
                }
            } else {
                if(j == 0) {
                    for(Int l = 0; l < ngl; l++) {
                        for(Int k = nglh; k < ngl; k++)
                            psiCor[i*2+j][l * ngl + k] = 0;
                    }
                } else {
                    for(Int l = 0; l < ngl; l++) {
                        for(Int k = 0; k < nglh; k++)
                            psiCor[i*2+j][l * ngl + k] = 0;
                    }
                }
                //Doulbe all values
                for(Int l = 0; l < ngl; l++) {
                    for(Int k = 0; k < ngl; k++)
                        psiCor[i*2+j][l * ngl + k] *= 2;
                }
            }
#if 0
            std::cout << "========" << j << "==========" << std::endl;
            for(Int k = 0; k < ngl; k++) {
                for(Int l = 0; l < ngl; l++)
                    std::cout << psiRef[i*2+j][k*ngl+l] << " ";
                std::cout << std::endl;
            }
            std::cout << "------------" << std::endl;
            for(Int k = 0; k < ngl; k++) {
                for(Int l = 0; l < ngl; l++)
                    std::cout << psiCor[i*2+j][k*ngl+l] << " ";
                std::cout << std::endl;
            }
#endif
        }
    }
}

