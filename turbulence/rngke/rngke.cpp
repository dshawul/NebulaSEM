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
	ScalarParams::enroll("eta0",&eta0);
	ScalarParams::enroll("beta",&beta);
}
void RNG_KE_Model::solve() {
	/*calculate C2eStar*/
	ScalarCellField C2eStar;
	{
		ScalarCellField eta;
		{
			STensorCellField S = sym(grad(U));
			ScalarCellField magS = sqrt((S & S) * 2);
			eta = magS * (k / x);
		}
		Scalar c;
		for(Int i = 0;i < eta.size();i++) {
			c = C2x + Cmu * pow(eta[i],3) * (1 - eta[i] / eta0) / 
			                       (1 + beta * pow(eta[i],3));
			if(c < 0) c = 0;
			C2eStar[i] = c;
		}
	}

	/*solve*/
	ScalarMeshMatrix M;
	ScalarFacetField mu;
	ScalarCellField w = (x / k);

	/*turbulent dissipation*/
	mu = cds(eddy_mu) / SigmaX;
	M = div(x,F,mu) 
		- lap(x,mu);
	M -= src(x,
		(C1x * w * G),                             //Su
		-(C2eStar * rho * w)                       //Sp  
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
