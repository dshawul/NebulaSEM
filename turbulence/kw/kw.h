#ifndef __KW_H
#define __KW_H

#include "turbulence.h"

struct KW_Model : public KX_Model {
	/*constructor*/
	KW_Model(VectorCellField&,ScalarFacetField&,Scalar&,Scalar&,bool&);

	/*others*/
	virtual void enroll();
	virtual void solve();
	void calcEddyMu() {
		eddy_mu = (rho * k) / x;
	};
	Scalar calcX(Scalar ustar,Scalar kappa,Scalar y) {
		return ustar / (kappa * y * sqrt(Cmu));
	}
};

#endif
