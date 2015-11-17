#include "realizableke.h"

/*
References:
    http://www.cfd-online.com/Wiki/Realisable_k-epsilon_model
	http://www.laturbolenza.com/?p=92
*/
REALIZABLE_KE_Model::REALIZABLE_KE_Model(VectorCellField& tU,ScalarFacetField& tF,ScalarCellField& trho,Scalar& tnu) :
	KX_Model(tU,tF,trho,tnu,"e"),
	CmuF(Scalar(0.09)),
	A0(4.04)
{
	SigmaK = 1.0;
	SigmaX = 1.2;
	C2x = 1.9;
}
void REALIZABLE_KE_Model::enroll() {
	using namespace Util;
	KX_Model::enroll();
	params.enroll("SigmaK",&SigmaK);
	params.enroll("SigmaE",&SigmaX);
	params.enroll("C2e",&C2x);
}
void REALIZABLE_KE_Model::calcEddyViscosity(const TensorCellField& gradU) {
	/*calculate CmuF*/
	STensorCellField S = sym(gradU);
	{
		TensorCellField O = skw(gradU);
		ScalarCellField Ustar = sqrt((S & S) + (O & O));
		ScalarCellField Sbar = sqrt(S & S);
		ScalarCellField W = ((mul(S,S) & S) / pow(Sbar,3.0)) * sqrt(6.0);
		W = min(max(W,-1.0),1.0);
		ScalarCellField As = sqrt(6.0) * cos(acos(W) / 3.0);
		CmuF = 1.0 / (A0 + As * Ustar * k / x);
		CmuF = min(CmuF,0.09);
	}
	/*calculate C1*/
	magS = sqrt((S & S) * 2.0);
	{
		ScalarCellField eta = magS * (k / x);
		C1 = max(eta/(eta + 5.0),0.43);
	}
	/*calculate viscosity*/
	KX_Model::calcEddyViscosity(gradU);
}
void REALIZABLE_KE_Model::solve() {
	ScalarCellMatrix M;
	ScalarCellField mu;

	/*turbulent dissipation*/
	mu = eddy_mu / SigmaX + rho * nu;
	M = transport(x, U, F, mu, x_UR,
				(C1 * rho * magS * x),
				-(C2x * rho * x / (k + sqrt(nu * x))), &rho);
	FixNearWallValues(M);
	Solve(M);
	x = max(x,Constants::MachineEpsilon);

	/*turbulent kinetic energy*/
	mu = eddy_mu / SigmaK + rho * nu;
	M = transport(k, U, F, mu, k_UR,
					Pk,
					-(rho * x / k), &rho);
	if(wallModel == STANDARD)
		FixNearWallValues(M);
	Solve(M);
	k = max(k,Constants::MachineEpsilon);
}
