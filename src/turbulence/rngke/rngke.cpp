#include "rngke.h"
/*
References:
    http://www.cfd-online.com/Wiki/RNG_k-epsilon_model
*/
RNG_KE_Model::RNG_KE_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
	KE_Model(tU,tF,trho,tnu,tSteady),
	eta0(4.38),
	beta(0.012)
{
	Cmu = 0.0845;
	SigmaK = 0.7194;
	SigmaX = 0.7194;
	C1x = 1.42;
	C2x = 1.68;
}
void RNG_KE_Model::enroll() {
	using namespace Util;
	KE_Model::enroll();
	params.enroll("eta0",&eta0);
	params.enroll("beta",&beta);
}
void RNG_KE_Model::calcEddyViscosity(const TensorCellField& gradU) {
	/*calculate C2eStar*/
	{
		ScalarCellField eta = sqrt(getS2(gradU)) * (k / x);
		C2eStar = C2x + Cmu * pow(eta,3.0) * (1 - eta / eta0) / 
			(1 + beta * pow(eta,3.0));
		C2eStar = max(C2eStar,0.0);
	}
	/*calculate viscosity*/
	KE_Model::calcEddyViscosity(gradU);
}
void RNG_KE_Model::solve() {
	ScalarMeshMatrix M;
	ScalarFacetField mu;

	/*turbulent dissipation*/
	mu = cds(eddy_mu) / SigmaX + rho * nu;;
	M = div(x,F,mu) 
		- lap(x,mu);
	M -= src(x,
		(C1x * Pk * x / k),                        //Su
		-(C2eStar * rho * x / k)                   //Sp  
		);
	if(Steady)
		M.Relax(x_UR);
	else
		M += ddt(x,rho);
	M.FixNearWallValues();
	Solve(M);
	x = max(x,Constants::MachineEpsilon);

	/*turbulent kinetic energy*/
	mu = cds(eddy_mu) / SigmaK + rho * nu;
	M = div(k,F,mu) 
		- lap(k,mu);
	M -= src(k,
		Pk,                                        //Su
		-(rho * x / k)                             //Sp
		);
	if(Steady)
		M.Relax(k_UR);
	else
		M += ddt(k,rho);
	if(wallModel == STANDARD)
		M.FixNearWallValues();
	Solve(M);
	k = max(k,Constants::MachineEpsilon);
}
