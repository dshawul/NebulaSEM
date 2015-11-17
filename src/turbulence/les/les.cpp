#include "les.h"
/*
References:
    http://www.cfd-online.com/Wiki/Smagorinsky-Lilly_model
*/
LES_Model::LES_Model(VectorCellField& tU,ScalarFacetField& tF,ScalarCellField& trho,Scalar& tnu) :
	MixingLength_Model(tU,tF,trho,tnu),
	Cs(0.11)
{
}
void LES_Model::enroll() {
	params.enroll("Cs",&Cs);
	MixingLength_Model::enroll();
}
void LES_Model::calcLengthScale() {
	ScalarCellField delta = pow(Mesh::cV,Scalar(1./3));
	lm = Cs * delta;
}


