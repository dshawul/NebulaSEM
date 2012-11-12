#include "realizableke.h"

/*
References:
    http://www.cfd-online.com/Wiki/Realisable_k-epsilon_model
	http://www.laturbolenza.com/?p=92
*/
REALIZABLE_KE_Model::REALIZABLE_KE_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
	KX_Model(tU,tF,trho,tnu,tSteady,"e"),
	CmuF(0.09)
{
	SigmaK = 1.0;
	SigmaX = 1.2;
	C2x = 1.9;
}
void REALIZABLE_KE_Model::enroll() {
	using namespace Util;
	KX_Model::enroll();
	ScalarParams::enroll("SigmaK",&SigmaK);
	ScalarParams::enroll("SigmaE",&SigmaX);
	ScalarParams::enroll("C2e",&C2x);
}
void REALIZABLE_KE_Model::solve() {
	ScalarCellField C1;
	ScalarCellField magS;
	/*calculate model constants*/
	{
		STensorCellField S = sym(grad(U));
		magS = sqrt((S & S) * 2);
		ScalarCellField eta = magS * (k / x);
		Scalar r;
		for(Int i = 0;i < C1.size();i++) {
			r = eta[i] / (5 + eta[i]);
			if(r < 0.43) r = 0.43;
			C1[i] = r;
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
		(C1 * rho * magS * x),                     //Su
		-(C2x * rho * x / (k + sqrt(nu * x)))      //Sp  
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

	/*calculate CmuF*/
	{
		const Scalar A0 = 4.04;
		TensorCellField gradU = grad(U);
		STensorCellField S = sym(gradU);
		TensorCellField O = skw(gradU);
		ScalarCellField Ustar = sqrt((S & S) + (O & O));
		ScalarCellField Sbar = sqrt(S & S);
		ScalarCellField W = ((mul(S,S) & S) / pow(Sbar,Scalar(3))) * sqrt(6.);
		for(Int i = 0;i < W.size();i++) {
			if(W[i] < Scalar(-1)) W[i] = -1;
			else if(W[i] > Scalar(1)) W[i] = 1;
		}
		ScalarCellField As = sqrt(6.) * cos(acos(W) / 3);
		CmuF = Scalar(1) / (A0 + As * Ustar / w);
		for(Int i = 0;i < CmuF.size();i++) {
			if(CmuF[i] > 0.09) CmuF[i] = 0.09;
		}
	}
	/*end*/
}
