#include "field.h"

namespace DG {
	Int Nop[3] = {0, 0, 0};
	Int NPX, NPY, NPZ;
	Int NP, NPI, NPMAT, NPF;
	Scalar Penalty;
	
	Scalar **psi[3];
	Scalar **dpsi[3];
	Scalar *xgl[3];
	Scalar *wgl[3];
	TensorCellField Jinv(false);
}
/** 
Compute legendre polynomial and its first & second derivatives, given p and x
*/
void DG::legendre(int p, Scalar x, Scalar& L0, Scalar& L0_1, Scalar& L0_2) {
    Scalar a,b;                
    Scalar L2, L2_1, L2_2;
    Scalar L1, L1_1, L1_2;
    L1 = 0; L1_1 = 0; L1_2 = 0;
    L0 = 1; L0_1 = 0; L0_2 = 0;

    for(int i = 1;i <= p;i++) {
        L2 = L1; L2_1 = L1_1; L2_2 = L1_2;
        L1 = L0; L1_1 = L0_1; L1_2 = L0_2;
        a = (2 * i - 1.0)/i;
        b = (i - 1.0)/i;
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
    int ph = floor( (p+1)/2 );
    Scalar x,dx;

	if(N == 1) {
		xgl[0] = 0;
		wgl[0] = 2;
		return;
	}

    for(int i = 0; i < ph; i++) {
   		x = cos((2 * i + 1) * Constants::PI / (2 * p + 1));
   		for(int k = 1; k <= 20; k++) {
            legendre(p,x,L0,L0_1,L0_2);
      		dx = -(1 - x * x) * L0_1 / (-2 * x * L0_1 + (1 - x * x) * L0_2);
      		x += dx;
      		if(fabs(dx) < 1.0e-20) 
         		break;
        }
   	    xgl[p - i] = x;
        wgl[p - i] = 2 / (p * (p + 1) * L0 * L0);
    }

    if (p+1 != 2*ph) {
   		x = 0;
   		legendre(p,x,L0,L0_1,L0_2);
   		xgl[ph] = x;
   		wgl[ph] = 2 / (p * (p + 1) * L0 * L0);
	}
   
    for(int i = 0; i < ph; i++) {
   		xgl[i] = -xgl[p - i];
   		wgl[i] = +wgl[p - i];
	}
}
/**
Compute Legendre-Gauss interpolation points and weights
*/
void DG::legendre_gauss(int N, Scalar* xgl, Scalar* wgl) {
    Scalar L0,L0_1,L0_2;
	int p = N - 1; 
    int ph = floor( (p+1)/2 );
    Scalar x,dx;

	if(N == 1) {
		xgl[0] = 0;
		wgl[0] = 2;
		return;
	}

    for(int i = 0; i < ph; i++) {
   		x = cos((2 * i + 1) * Constants::PI / (2 * p + 1));
   		for(int k = 1; k <= 20; k++) {
            legendre(p + 1,x,L0,L0_1,L0_2);
      		dx = -L0/L0_1;
      		x += dx;
      		if(fabs(dx) < 1.0e-20) 
         		break;
        }
   	    xgl[p - i] = x;
        wgl[p - i] = 2 / ((1 - x * x) * L0_1 * L0_1);
    }

    if (p+1 != 2*ph) {
   		x = 0;
   		legendre(p + 1,x,L0,L0_1,L0_2);
   		xgl[ph] = x;
   		wgl[ph] = 2 / ((1 - x * x) * L0_1 * L0_1);
	}
   
    for(int i = 0; i < ph; i++) {
   		xgl[i] = -xgl[p - i];
   		wgl[i] = +wgl[p - i];
	}
}
/**
Compute equispaced newton-cotes interpolation points and weights
*/
void DG::newton_cotes(int N, Scalar* xgl, Scalar* wgl) {
	int p = N - 1; 
	int ph = floor( (p+1)/2 );
	
	if(N == 1) {
		xgl[0] = 0;
		wgl[0] = 2;
		return;
	}
	
    for(int i = 0; i < N; i++)
		xgl[i] = -1 + (2.0 * i) / p;
	
	switch(N) {
		case 2:
		wgl[0] = 1;
		break;
		case 3:
		wgl[0] = 1.0/3; wgl[1] = 4.0/3;
		break;
		case 4:
		wgl[0] = 1.0/4; wgl[1] = 3.0/4;
		break;
		case 5:
		wgl[0] = 7.0/45; wgl[1] = 32.0/45; wgl[2] = 12.0/45;
		break;
		case 6:
		wgl[0] = 19.0/144; wgl[1] = 75.0/144; wgl[2] = 50.0/144;
		break;
		case 7:
		wgl[0] = 41.0/420; wgl[1] = 216.0/420; wgl[2] = 27.0/420; wgl[3]=272.0/420;
		break;
		case 8:
		wgl[0] = 751.0/8640; wgl[1] = 3577.0/8640; wgl[2] = 1323.0/8640; wgl[3]=2989.0/8640;
		break;
		case 9:
		wgl[0] = 989.0/14175; wgl[1] = 5888.0/14175; wgl[2] = -928.0/14175; wgl[3]=10496.0/14175; 
		wgl[4] = -4540.0/14175;
		break;
		case 10:
		wgl[0] = 2857.0/44800; wgl[1] = 15741.0/44800; wgl[2] = 1080.0/44800; wgl[3]=19344.0/44800; 
		wgl[4] = 5778.0/44800;
		break;
		case 11:
		wgl[0] = 16067.0/299376; wgl[1] = 106300.0/299376; wgl[2] = -48525.0/299376; wgl[3]=272400.0/299376; 
		wgl[4] = -260550.0/299376; wgl[5] = 427368.0/299376;
		break;
	}
	
    for(int i = 0; i < ph; i++)
   		wgl[p - i] = +wgl[i];
}
/**
Compute cardinal basis functions
*/
void DG::cardinal_basis(int v, int N, Scalar* xgl, Scalar* psi) {
	for(int i = 0;i < N;i++) {
		if(i != v) psi[i] = 0;
		else psi[i] = 1;
	}
}
/**
Compute lagrange basis function derivatives at given point x, given interpolation points
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
Compute legendre basis function derivatives at given point x, given LGL points
*/
void DG::legendre_basis_derivative(int v, int N, Scalar* xgl, Scalar* dpsi) {
	Scalar xi,xj;
	Scalar x = xgl[v];
	
	Scalar* bb = new Scalar[N];
	for(int j = 0;j < N;j++) {
		bb[j] = 0;
		xj = xgl[j];
		for(int i = 0;i < N;i++) {
			if(i != j) {
				xi = xgl[i];
				bb[j] += log(fabs(xj - xi));
			}
		}
	}
	
	Scalar cc = 0;
	for(int i = 0;i < N;i++) {
		if(i != v) {
			xi = xgl[i];
			Scalar val = exp(bb[v] - bb[i]) / (x - xi);
			if((v + i + 2) & 1) val = -val;
			dpsi[i] = val;
			cc += val;
		}
	}
	dpsi[v] = -cc;
	
	delete bb;
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
	
	//compute coordinates of nodes via transfinite interpolation
	gFO.assign(gFacets.size() * NPF,gCells.size() * NP);
	gFN.assign(gFacets.size() * NPF,gCells.size() * NP);
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
		Vertex ev[12][3];
		for(Int i = 0;i < 12;i++) {
			ev[i][0] = vp[sides[i][0]];
			ev[i][1] = vp[sides[i][1]];
		}
		
		//coordinates
		Vertex v,vd[12],vf[6];
		Scalar rx,ry,rz;
		
#define ADDV(w,m,ev,vd) {								\
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
	rx = (xgl[0][i] + 1) / 2;							\
	ry = (xgl[1][j] + 1) / 2;							\
	rz = (xgl[2][k] + 1) / 2;							\
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
	}	
#undef ADDV
#undef ADDF
#undef ADDC
#undef ADD
	
	for(Int ci = 0; ci < gBCS;ci++) {
		Cell& c = gCells[ci];
		forEach(c,mm) {
			Int face = gFaceID[ci][mm];
			Int fi = c[mm];
			Int cj = gFNC[fi];
			if(cj == ci) 
				continue;
			
#define ADD() {												\
	gFO[indf] = index0;										\
	gFN[indf] = index1;										\
	if(index1 >= gBCSfield) {								\
		Vector dir = unit(fN[indf]);						\
		Scalar d = dot(fC[indf] - cC[index0],dir);			\
		cC[index1] = cC[index0] + d * dir;					\
	}														\
	Scalar d;												\
	if(equal(cC[index0],cC[index1])) d = 0.5;				\
	else d = dot(fC[indf] - cC[index0],fN[indf]) / 			\
			 dot(cC[index1] - cC[index0],fN[indf]);			\
	fC[indf] = cC[index0] + d * (cC[index1] - cC[index0]);	\
	fN[indf] *= wgt;										\
}

			if(face == 0 || face == 1) {
				Int ff = (face == 0) ? 0 : (NPZ - 1);
				forEachLglXY(i,j) {
					Scalar wgt = wgl[0][i] * wgl[1][j] / 4;
					Int indf = fi * NPF + i * NPY + j;
					Int index0 = INDEX4(ci,i,j,ff);
					Int index1 = INDEX4(cj,i,j,(NPZ - 1) - ff);
					ADD();
				}
			} else if(face == 2 || face == 3) {
				Int ff = (face == 2) ? 0 : (NPY - 1);
				forEachLglXZ(i,k) {
					Scalar wgt = wgl[0][i] * wgl[2][k] / 4;
					Int indf = fi * NPF + i * NPZ + k;
					Int index0 = INDEX4(ci,i,ff,k);
					Int index1 = INDEX4(cj,i,(NPY - 1) - ff,k);
					ADD();
				}
			} else {
				Int ff = (face == 4) ? 0 : (NPX - 1);
				forEachLglYZ(j,k) {
					Scalar wgt = wgl[1][j] * wgl[2][k] / 4;
					Int indf = fi * NPF + j * NPZ + k;
					Int index0 = INDEX4(ci,ff,j,k);
					Int index1 = INDEX4(cj,(NPX - 1) - ff,j,k);
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
	for(Int i = 0;i < 3;i++) {
		Int ngl = Nop[i] + 1;
		xgl[i] = new Scalar[ngl];
		wgl[i] = new Scalar[ngl];
		legendre_gauss_lobatto(ngl,xgl[i],wgl[i]);
		psi[i] = new Scalar*[ngl];
		dpsi[i] = new Scalar*[ngl];
		for(Int j = 0;j < ngl;j++) {
			psi[i][j] = new Scalar[ngl];
			dpsi[i][j] = new Scalar[ngl];
			cardinal_basis(j,ngl,xgl[i],psi[i][j]);
			lagrange_basis_derivative(j,ngl,xgl[i],dpsi[i][j]);
		}
	}
	
	//init geometry
	init_geom();
	
	//Compute Jacobian matrix
	Jinv.deallocate(false);
	Jinv.construct();
	for(Int ci = 0; ci < gBCS;ci++) {
		forEachLgl(ii,jj,kk) {
			Tensor Ji(Scalar(0));
			
#define JACD(im,jm,km) {									\
	Int index = INDEX4(ci,im,jm,km);						\
	Vector& C = cC[index];									\
	Vector dpsi_ij;											\
	DPSI(dpsi_ij,im,jm,km);									\
	Ji += mul(dpsi_ij,C);									\
}
			forEachLglX(i) JACD(i,jj,kk);
			forEachLglY(j) if(j != jj) JACD(ii,j,kk);
			forEachLglZ(k) if(k != kk) JACD(ii,jj,k);
#undef JACD
			
			if(NPX == 1) {Ji[XX] = 1; Ji[YX] = 0; Ji[ZX] = 0;}
			if(NPY == 1) {Ji[YY] = 1; Ji[XY] = 0; Ji[ZY] = 0;}
			if(NPZ == 1) {Ji[ZZ] = 1; Ji[XZ] = 0; Ji[YZ] = 0;}
			Ji = inv(Ji);
			if(NPX == 1) Ji[XX] = 0;
			if(NPY == 1) Ji[YY] = 0;
			if(NPZ == 1) Ji[ZZ] = 0;
			
			Int index = INDEX4(ci,ii,jj,kk);
			Jinv[index] = Ji;
		}
	}
}

