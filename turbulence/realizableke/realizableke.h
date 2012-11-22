#ifndef __REALIZABLEKE_H
#define __REALIZABLEKE_H

#include "turbulence.h"

struct REALIZABLE_KE_Model : public KX_Model {
	/*model coefficients*/
	ScalarCellField CmuF;
	Scalar A0;

	/*constructor*/
	REALIZABLE_KE_Model(VectorCellField&,ScalarFacetField&,Scalar&,Scalar&,bool&);

	void enroll();
	void solve();
	void calcEddyMu() {
		eddy_mu = (rho * CmuF * k * k) / x;
	};
	Scalar calcX(Scalar ustar,Scalar kappa,Scalar y) {
		return pow(ustar,3.0) / (kappa * y);
	}
};

#endif
