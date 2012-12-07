#ifndef __LES_H
#define __LES_H

#include "mixing_length.h"

struct LES_Model : public MixingLength_Model {
	/*model coefficients*/
	Scalar Cs;

	/*constructor*/
	LES_Model(VectorCellField&,ScalarFacetField&,Scalar&,Scalar&,bool&);

	/*others*/
	virtual void enroll();
	virtual void calcLengthScale();
};

#endif
