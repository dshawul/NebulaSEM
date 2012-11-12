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
	ScalarParams::enroll("Cmu",&Cmu);
	ScalarParams::enroll("SigmaK",&SigmaK);
	ScalarParams::enroll("SigmaE",&SigmaX);
	ScalarParams::enroll("C1e",&C1x);
	ScalarParams::enroll("C2e",&C2x);
}
void KE_Model::solve() {
	ScalarMeshMatrix M;
	ScalarFacetField mu;
	ScalarCellField w = (x / k);

	/*turbulent dissipation*/
	mu = cds(eddy_mu) / SigmaX;
	M = div(x,F,mu) 
		- lap(x,mu);
	M -= src(x,
		(C1x * w * G),                             //Su
		-(C2x * rho * w)                           //Sp  
		);
	if(Steady)
		M.Relax(UR);
	else
		M += ddt(x,rho);
	Solve(M);

	/*turbulent kinetic energy*/
	mu = cds(eddy_mu) / SigmaK;
	M = div(k,F,mu) 
		- lap(k,mu);
	M -= src(k,
		G,                                         //Su
		-(rho * w)                                 //Sp
		);
	if(Steady)
		M.Relax(UR);
	else
		M += ddt(k,rho);
	Solve(M);
}
