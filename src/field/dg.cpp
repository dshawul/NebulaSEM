#include "dg.h"

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
void DG::LGL(int N, Scalar* xgl, Scalar* wgl) {
    Scalar L0,L0_1,L0_2;
	int p = N - 1; 
    int ph = floor( (p+1)/2 );
    Scalar x,dx;

    //find roots of legendre polys
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

    //Check for Zero Root
    if (p+1 != 2*ph) {
   		x=0;
   		legendre(p,x,L0,L0_1,L0_2);
   		xgl[ph]=x;
   		wgl[ph]=2/(p*(p+1)*L0*L0);
	}
   
    //Find remainder of roots via symmetry
    for(int i = 1; i <= ph; i++) {
   		xgl[i-1]=-xgl[p+1-i];
   		wgl[i-1]=+wgl[p+1-i];
	}
}