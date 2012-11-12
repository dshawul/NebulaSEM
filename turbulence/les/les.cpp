#include "les.h"
/*
References:
    http://www.cfd-online.com/Wiki/Smagorinsky-Lilly_model
*/
LES_Model::LES_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
	Turbulence_Model(tU,tF,trho,tnu,tSteady),
	Cs(0.11)
{
}
void LES_Model::enroll() {
	Util::ScalarParams::enroll("Cs",&Cs);
	Turbulence_Model::enroll();
}
void LES_Model::calcEddyMu(const TensorCellField& gradU) {
	using namespace Mesh;
	STensorCellField S = sym(gradU);
	ScalarCellField delta = pow(cV,Scalar(1./3));
	eddy_mu = rho * pow(Cs * delta,Scalar(2)) * sqrt(2 * (S & S));
	eddy_mu.FillBoundaryValues();
}
void LES_Model::addTurbulentStress(VectorMeshMatrix& M) {
	TensorCellField gradU = grad(U);
	calcEddyMu(gradU);

	/*sub-grid scale stresses*/
	ScalarCellField eff_mu = eddy_mu + rho * nu;
	M -= lap(U,eff_mu);
	M -= sum(div(eddy_mu * dev(trn(gradU),2)));
}
