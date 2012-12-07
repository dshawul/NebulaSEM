#include "realizableke.h"

/*
References:
    http://www.cfd-online.com/Wiki/Realisable_k-epsilon_model
	http://www.laturbolenza.com/?p=92
*/
REALIZABLE_KE_Model::REALIZABLE_KE_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
	KX_Model(tU,tF,trho,tnu,tSteady,"e"),
	CmuF(0.09),
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
void REALIZABLE_KE_Model::solve() {
	ScalarCellField C1;
	ScalarCellField magS;
	/*calculate model constants*/
	{
		STensorCellField S = sym(grad(U));
		magS = sqrt((S & S) * 2.0);
		ScalarCellField eta = magS * (k / x);
		C1 = max(eta/(eta + 5.0),0.43);
	}

	/*solve*/
	ScalarMeshMatrix M;
	ScalarFacetField mu;

	/*turbulent dissipation*/
	mu = cds(eddy_mu) / SigmaX + rho * nu;
	M = div(x,F,mu) 
		- lap(x,mu);
	M -= src(x,
		(C1 * rho * magS * x),                     //Su
		-(C2x * rho * x / (k + sqrt(nu * x)))      //Sp  
		);
	if(Steady)
		M.Relax(x_UR);
	else
		M += ddt(x,rho);
	Solve(M);

	/*turbulent kinetic energy*/
	mu = cds(eddy_mu) / SigmaK + rho * nu;
	M = div(k,F,mu) 
		- lap(k,mu);
	M -= src(k,
		G,                                         //Su
		-(rho * x / k)                             //Sp
		);
	if(Steady)
		M.Relax(k_UR);
	else
		M += ddt(k,rho);
	Solve(M);

	/*calculate CmuF*/
	{
		TensorCellField gradU = grad(U);
		STensorCellField S = sym(gradU);
		TensorCellField O = skw(gradU);
		ScalarCellField Ustar = sqrt((S & S) + (O & O));
		ScalarCellField Sbar = sqrt(S & S);
		ScalarCellField W = ((mul(S,S) & S) / pow(Sbar,3.0)) * sqrt(6.0);
		W = min(max(W,-1.0),1.0);
		ScalarCellField As = sqrt(6.0) * cos(acos(W) / 3.0);
		CmuF = 1.0 / (A0 + As * Ustar * k / x);
		CmuF = min(CmuF,0.09);
	}
	/*end*/
}
