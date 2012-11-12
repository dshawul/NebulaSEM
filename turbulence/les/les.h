#ifndef __LES_H
#define __LES_H

#include "turbulence.h"

struct LES_Model : public Turbulence_Model {
	/*model coefficients*/
	Scalar Cs;

	/*turbulence fields*/
	ScalarCellField eddy_mu;  

	/*constructor*/
	LES_Model(VectorCellField&,ScalarFacetField&,Scalar&,Scalar&,bool&);

	/*others*/
	virtual void enroll();
	virtual void calcEddyMu(const TensorCellField& gradU);
	virtual void addTurbulentStress(VectorMeshMatrix&);
};

#endif
