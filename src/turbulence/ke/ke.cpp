#include "ke.h"
/*
References:
    http://www.cfd-online.com/Wiki/Standard_k-epsilon_model
*/
KE_Model::KE_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu) :
	KX_Model(tU,tF,trho,tnu,"e")
{
	Cmu = 0.09;
	SigmaK = 1;
	SigmaX = 1.314;
	C1x = 1.44;
	C2x = 1.92;
}
void KE_Model::enroll() {
	using namespace Util;
	KX_Model::enroll();
	params.enroll("Cmu",&Cmu);
	params.enroll("SigmaK",&SigmaK);
	params.enroll("SigmaE",&SigmaX);
	params.enroll("C1e",&C1x);
	params.enroll("C2e",&C2x);
}
void KE_Model::solve() {
	ScalarCellMatrix M;
	ScalarCellField mu;

	/*turbulent dissipation*/
	mu = eddy_mu / SigmaX + rho * nu;
	M = transport(x, rho * U, F, mu, rho, x_UR,
			(C1x * Pk * x / k),
			-(C2x * rho * x / k));
	FixNearWallValues(M);
	Solve(M);
	x = max(x,Constants::MachineEpsilon);

	/*turbulent kinetic energy*/
	mu = eddy_mu / SigmaK + rho * nu;
	M = transport(k, rho * U, F, mu, rho, k_UR,
				Pk,
				-(rho * x / k));
	if(wallModel == STANDARD)
		FixNearWallValues(M);
	Solve(M);
	k = max(k,Constants::MachineEpsilon);
}
