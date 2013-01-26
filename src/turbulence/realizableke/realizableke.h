#ifndef __REALIZABLEKE_H
#define __REALIZABLEKE_H

#include "turbulence.h"

struct REALIZABLE_KE_Model : public KX_Model {
	/*model coefficients*/
	ScalarCellField CmuF;
	ScalarCellField C1;
	ScalarCellField magS;
	Scalar A0;

	/*constructor*/
	REALIZABLE_KE_Model(VectorCellField&,ScalarFacetField&,Scalar&,Scalar&,bool&);

	/*others*/
	virtual void enroll();
	virtual void solve();
	virtual void calcEddyMu() {
		eddy_mu = (rho * CmuF * k * k) / x;
	};
	virtual Scalar calcX(Scalar ustar,Scalar kappa,Scalar y) {
		return pow(ustar,Scalar(3)) / (kappa * y);
	}
	virtual Scalar getCmu(Int i) { 
		return CmuF[i]; 
	}
	virtual void calcEddyViscosity(const TensorCellField& gradU);
};

#endif
