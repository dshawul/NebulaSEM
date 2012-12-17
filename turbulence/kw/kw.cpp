#include "kw.h"
/*
References:
    http://www.cfd-online.com/Wiki/Wilcox%27s_k-omega_model
*/
KW_Model::KW_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
	KX_Model(tU,tF,trho,tnu,tSteady,"w")
{
	Cmu = 0.09;
	SigmaK = 2;
	SigmaX = 2;
	C1x = 5./9;
	C2x = 3./40;
}
void KW_Model::enroll() {
	using namespace Util;
	KX_Model::enroll();
	params.enroll("Cmu",&Cmu);
	params.enroll("SigmaK",&SigmaK);
	params.enroll("SigmaW",&SigmaX);
	params.enroll("C1w",&C1x);
	params.enroll("C2w",&C2x);
}
void KW_Model::solve() {
	ScalarMeshMatrix M;
	ScalarFacetField mu;

	/*turbulent dissipation*/
	mu = cds(eddy_mu) / SigmaX + rho * nu;;
	M = div(x,F,mu) 
		- lap(x,mu);
	M -= src(x,
		(C1x * G * x / k),                         //Su
		-(C2x * x * rho)                           //Sp  
		);
	if(Steady)
		M.Relax(x_UR);
	else
		M += ddt(x,rho);
	Solve(M);
	x = max(x,Constants::MachineEpsilon);

	/*turbulent kinetic energy*/
	mu = cds(eddy_mu) / SigmaK + rho * nu;;
	M = div(k,F,mu) 
		- lap(k,mu);
	M -= src(k,
		G,                                         //Su
		-(Cmu * x * rho)                           //Sp
		);
	if(Steady)
		M.Relax(k_UR);
	else
		M += ddt(k,rho);
	Solve(M);
	k = max(k,Constants::MachineEpsilon);
}
