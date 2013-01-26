#include "ke.h"
/*
References:
    http://www.cfd-online.com/Wiki/Standard_k-epsilon_model
*/
KE_Model::KE_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
	KX_Model(tU,tF,trho,tnu,tSteady,"e")
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
	ScalarMeshMatrix M;
	ScalarFacetField mu;

	/*turbulent dissipation*/
	mu = cds(eddy_mu) / SigmaX + rho * nu;
	M = div(x,F,mu) 
		- lap(x,mu);
	M -= src(x,
		(C1x * Pk * x / k),                        //Su
		-(C2x * rho * x / k)                       //Sp  
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
