#include "field.h"

namespace DG {
	Int Nop[3] = {0, 0, 0};
	Int NPX, NPY, NPZ, NP, NPMAT, NPF;
	
	Scalar **psi[3];
	Scalar **dpsi[3];
	Scalar *xgl[3];
	Scalar *wgl[3];
	TensorCellField J(false);
	TensorCellField Jinv(false);
}

/** 
Compute legendre polynomial and its first & second derivatives
*/
void DG::legendre(int p, Scalar x,
                   Scalar& L0,Scalar& L0_1,Scalar& L0_2) {
    Scalar a,b;                
    Scalar L2,L2_1,L2_2;
    Scalar L1,L1_1,L1_2;
    L1 = 0,L1_1 = 0,L1_2 = 0;
    L0 = 1,L0_1 = 0,L0_2 = 0;

    for(int i = 1;i <= p;i++) {
        L2=L1;L2_1=L1_1;L2_2=L1_2;
        L1=L0;L1_1=L0_1;L1_2=L0_2;
        a=(2*i-1.0)/i;
        b=(i-1.0)/i;
        L0=a*x*L1 - b*L2;
        L0_1=a*(L1 + x*L1_1) - b*L2_1;
        L0_2=a*(2*L1_1 + x*L1_2) - b*L2_2;
    }
}

/**
Compute Legendre-Gauss-Lobato interpolation points and weights
*/
void DG::legendre_gauss_lobatto(int N, Scalar* xgl, Scalar* wgl) {
    Scalar L0,L0_1,L0_2;
	int p = N - 1; 
    int ph = floor( (p+1)/2 );
    Scalar x,dx;

	if(N == 1) {
		xgl[0] = 0;
		wgl[0] = 2;
		return;
	}

    for(int i = 1; i <= ph; i++) {
   		x=cos((2*i-1)*Constants::PI/(2*p+1));
   		for(int k = 1; k <= 20; k++) {
            legendre(p,x,L0,L0_1,L0_2);
      		dx=-(1-x*x)*L0_1/(-2*x*L0_1 + (1-x*x)*L0_2);
      		x=x+dx;
      		if(fabs(dx) < 1.0e-20) 
         		break;
        }
   	    xgl[p+1-i]=x;
        wgl[p+1-i]=2/(p*(p+1)*L0*L0);
    }

    if (p+1 != 2*ph) {
   		x=0;
   		legendre(p,x,L0,L0_1,L0_2);
   		xgl[ph]=x;
   		wgl[ph]=2/(p*(p+1)*L0*L0);
	}
   
    for(int i = 1; i <= ph; i++) {
   		xgl[i-1]=-xgl[p+1-i];
   		wgl[i-1]=+wgl[p+1-i];
	}
}
/**
Compute lagrange basis function at given point x, given LGL points
*/
void DG::lagrange(Scalar x, int N, Scalar* xgl, Scalar* psi) {
	Scalar xi,xj;
	for(int i = 0;i < N;i++) {
		psi[i] = 1;
		xi = xgl[i];
		for(int j = 0;j < N;j++) {
			if(i != j) {
				xj = xgl[j];
				psi[i] *= ((x - xj) / (xi - xj));
			}
		}
	}
}
/**
Compute lagrange basis function derivatives at given point x, given LGL points
*/
void DG::lagrange_der(Scalar x, int N, Scalar* xgl, Scalar* dpsi) {
	Scalar xi,xj,xk,prod;
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
	NP = NPX * NPY * NPZ;
	if(NP > 1) NPMAT = NP * NP;
	else NPMAT = 0;
	if(NPX <= NPY && NPX <= NPZ)
		NPF = NPY * NPZ;
	else if(NPY <= NPX && NPY <= NPZ)
		NPF = NPX * NPZ;
	else
	    NPF = NPX * NPY;
}
/**
Initialize basis functions
*/
void DG::init_basis() {
	using namespace Mesh;
	using namespace Constants;

	//directional LGL
	for(int i = 0;i < 3;i++) {
		Int ngl = Nop[i] + 1;
		
		xgl[i] = new Scalar[ngl];
		wgl[i] = new Scalar[ngl];
		legendre_gauss_lobatto(ngl,xgl[i],wgl[i]);
		psi[i] = new Scalar*[ngl];
		dpsi[i] = new Scalar*[ngl];
		for(int j = 0;j < ngl;j++) {
			psi[i][j] = new Scalar[ngl];
			dpsi[i][j] = new Scalar[ngl];
			lagrange(xgl[i][j],ngl,xgl[i],psi[i][j]);
			lagrange_der(xgl[i][j],ngl,xgl[i],dpsi[i][j]);
		}
	}

	//compute coordinates of nodes via transfinite interpolation
	gFO.assign(gFacets.size() * NPF,0);
	gFN.assign(gFacets.size() * NPF,0);
	for(Int ci = 0; ci < gBCS;ci++) {
		Cell& c = gCells[ci];
		Facet& f1 = gFacets[c[0]];
		Facet& f2 = gFacets[c[1]];
		
		//vertices
		Vertex vp[8];
		forEach(f1,i)
			vp[i + 0] = gVertices[f1[i]];
		forEach(f2,i)
			vp[i + 4] = gVertices[f2[i]];	
			
		//edges
		static const int sides[12][2] = {
			{0,1}, {3,2}, {7,6}, {4,5},
			{0,3}, {1,2}, {5,6}, {4,7},
			{0,4}, {1,5}, {2,6}, {3,7}
		};
		Vertex ev[12][3];
		for(Int i = 0;i < 12;i++) {
			ev[i][0] = vp[sides[i][0]];
			ev[i][1] = vp[sides[i][1]];
		}
		//coordinates
		Vertex v,vd[12],vf[6];
		Scalar rx,ry,rz;
		
#define ADDV(w,mm,ev,vd) {								\
	Scalar m = (mm + 1) / 2;							\
	vd[w] = (1 - m) * ev[w][0] + (m) * ev[w][1];		\
}

#define ADDF(w,rr,rs,i00,i01,i10,i11,ir0,ir1,i0s,i1s) {	\
	vf[w] = Interpolate_face(							\
			rr,rs,										\
			vp[i00],vp[i01],vp[i10],vp[i11],			\
			vd[ir0],vd[ir1],vd[i0s],vd[i1s]);			\
}

#define ADDC() {										\
	v = Interpolate_cell(								\
			rx,ry,rz,									\
			vp[0],vp[4],vp[3],vp[7],					\
			vp[1],vp[5],vp[2],vp[6],					\
			vd[0],vd[3],vd[1],vd[2],					\
			vd[4],vd[7],vd[5],vd[6],					\
			vd[8],vd[11],vd[9],vd[10],					\
			vf[4],vf[5],vf[2],vf[3],vf[0],vf[1]);		\
}

#define ADD() {                                     	\
	rx = xgl[0][i];										\
	ry = xgl[1][j];										\
	rz = xgl[2][k];										\
	ADDV(0 ,rx,ev,vd);									\
	ADDV(1 ,rx,ev,vd);									\
	ADDV(2 ,rx,ev,vd);									\
	ADDV(3 ,rx,ev,vd);									\
	ADDV(4 ,ry,ev,vd);									\
	ADDV(5 ,ry,ev,vd);									\
	ADDV(6 ,ry,ev,vd);									\
	ADDV(7 ,ry,ev,vd);									\
	ADDV(8 ,rz,ev,vd);									\
	ADDV(9 ,rz,ev,vd);									\
	ADDV(10,rz,ev,vd);									\
	ADDV(11,rz,ev,vd);									\
	ADDF(0, rx,ry, 0,3,1,2, 0,1,4,5);					\
	ADDF(1, rx,ry, 4,7,5,6, 3,2,7,6);					\
	ADDF(2, rx,rz, 0,4,1,5, 0,3,8,9);					\
	ADDF(3, rx,rz, 3,7,2,6, 1,2,11,10);					\
	ADDF(4, ry,rz, 0,4,3,7, 4,7,8,11);					\
	ADDF(5, ry,rz, 1,5,2,6, 5,6,9,10);					\
	ADDC();												\
};

		forEachLgl(i,j,k) {
			ADD();
			
			Scalar wgt = wgl[0][i] * wgl[1][j] * wgl[2][k] / 8;
			Int index = INDEX4(ci,i,j,k);
			cC[index] = v;
			cV[index] *= wgt;
		}
		
		forEach(c,face) {
			Int fi = c[face];
			Int cj = gFNC[fi];
			if(cj == ci) 
				continue;
			
			if(face == 0 || face == 1) {
				Int ff = (face == 0) ? 0 : (NPZ - 1);
				forEachLglXY(i,j) {
					Scalar wgt = wgl[0][i] * wgl[1][j] / 4;
					Int indf = fi * NPF + i * NPY + j;
					Int index0 = INDEX4(ci,i,j,ff);
					Int index1 = INDEX4(cj,i,j,NPZ - 1 - ff);
					gFO[indf] = index0;
					gFN[indf] = index1;
					fC[indf] = cC[index0];
					fN[indf] *= wgt;
				}
			} else if(face == 2 || face == 3) {
				Int ff = (face == 2) ? 0 : (NPY - 1);
				forEachLglXZ(i,k) {
					Scalar wgt = wgl[0][i] * wgl[2][k] / 4;
					Int indf = fi * NPF + i * NPZ + k;
					Int index0 = INDEX4(ci,i,ff,k);
					Int index1 = INDEX4(cj,i,NPY - 1 - ff,k);
					gFO[indf] = index0;
					gFN[indf] = index1;
					fC[indf] = cC[index0];
					fN[indf] *= wgt;
				}
			} else {
				Int ff = (face == 4) ? 0 : (NPX - 1);
				forEachLglYZ(j,k) {
					Scalar wgt = wgl[1][j] * wgl[2][k] / 4;
					Int indf = fi * NPF + j * NPZ + k;
					Int index0 = INDEX4(ci,ff,j,k);
					Int index1 = INDEX4(cj,NPX - 1 - ff,j,k);
					gFO[indf] = index0;
					gFN[indf] = index1;
					fC[indf] = cC[index0];
					fN[indf] *= wgt;
				}
			}
		}

#undef ADDV
#undef ADDF
#undef ADDC
#undef ADD
	}
	//boundary cells
	forEachS(gCells,i,gBCS) {
		Int faceid = gCells[i][0];
		Scalar offset = mag(gVertices[gFacets[faceid][0]] - 
			                gVertices[gFacets[faceid][1]]);
		for(Int n = 0; n < NPF;n++) {
			Int k = faceid * NPF + n;
			cV[gFN[k]] = cV[gFO[k]];
			cC[gFN[k]] = cC[gFO[k]] + offset * unit(fN[k]);
		}
	}
	//Compute Jacobian matrix
	J.deallocate(false);
	J.construct();
	Jinv.deallocate(false);
	Jinv.construct();
	for(Int ci = 0; ci < gBCS;ci++) {
		forEachLgl(ii,jj,kk) {
			Tensor Ji(0.0);

			forEachLgl(i,j,k) {
				Int index = INDEX4(ci,i,j,k);
				Vector& C = cC[index];
				Scalar dpsi_0 = dpsi[0][ii][i] *  psi[1][jj][j] *  psi[2][kk][k];
				Scalar dpsi_1 =  psi[0][ii][i] * dpsi[1][jj][j] *  psi[2][kk][k];
				Scalar dpsi_2 =  psi[0][ii][i] *  psi[1][jj][j] * dpsi[2][kk][k];
				Ji[XX] += C[0] * dpsi_0;
				Ji[YX] += C[0] * dpsi_1;
				Ji[ZX] += C[0] * dpsi_2;
				Ji[XY] += C[1] * dpsi_0;
				Ji[YY] += C[1] * dpsi_1;
				Ji[ZY] += C[1] * dpsi_2;
				Ji[XZ] += C[2] * dpsi_0;
				Ji[YZ] += C[2] * dpsi_1;
				Ji[ZZ] += C[2] * dpsi_2;
			}

			Int index = INDEX4(ci,ii,jj,kk);
			J[index] = Ji;
			Jinv[index] = inv(Ji);
		}
	}
			
}

