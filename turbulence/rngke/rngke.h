#ifndef __RNG_KE_H
#define __RNG_KE_H

#include "ke.h"

struct RNG_KE_Model : public KE_Model {
	/*model coefficients*/
	Scalar eta0;
	Scalar beta;

	/*constructor*/
	RNG_KE_Model(VectorCellField&,ScalarFacetField&,Scalar&,Scalar&,bool&);

	void enroll();
	void solve();
};

#endif
