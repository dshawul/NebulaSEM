#include "mixing_length.h"
/**
\verbatim
References:
    Book by Pope pg. 369
Description:
	Velocity and time scales are modelled as:
		l* = lm
		u* = lm * |S|
		eddy_nu = u*l* 
		        = (lm^2) * |S|
    Generalization of the mixing length model for 3D flows:
	   by Smagorinsky (1963).
		  eddy_nu = (lm^2) * sqrt(2 * (S & S))
	   by Baldwin & Lomax (1978)
		  eddy_nu = (lm^2) * sqrt(2 * (O & O))
	The turbulent kinetic energy k can be approximated by equating turbulent 
	viscosity eddy_nu with the one from Prandtl/Smagorinsky one equation models. 
	    u* = C * k^1/2
		eddy_nu = C * k^1/2 * lm
	Equating with the above eqn yields
	    k = (lm / C)^2 * (2 * (S & S))

	For high-Re flows, the mixing length close to the wall is set :
	    lm = kappa * y_wall
	Thus for Smagornsky LES model 
	    lm = min(Cs * Delta, kappa * y_wall)
\endverbatim
*/
MixingLength_Model::MixingLength_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu) :
	EddyViscosity_Model(tU,tF,trho,tnu),
	mixingLength(0),
	C(0.55),
	wallDamping(1),
	kappa(0.41)
{
}
void MixingLength_Model::enroll() {
	using namespace Util;
	params.enroll("mixing_length",&mixingLength);
	Option* op = new BoolOption(&wallDamping);
	params.enroll("wall_damping",op);
	params.enroll("kappa",&kappa);
	params.enroll("C",&C);
	EddyViscosity_Model::enroll();
}
ScalarCellField MixingLength_Model::getK() {
	return pow(lm / C,2.0) * getS2(gradi(U));
}
void MixingLength_Model::calcEddyViscosity(const TensorCellField& gradU) {
	calcLengthScale();
	if(wallDamping)
		lm = min(kappa * Mesh::yWall,lm);
	eddy_mu = rho * pow(lm,Scalar(2)) * sqrt(getS2(gradU));
}
void MixingLength_Model::applyWallFunction(Int f,LawOfWall& low) {
	using namespace Mesh;
	Int c1 = gFO[f];
	Int c2 = gFN[f];

	/*calc ustar*/
	Scalar ustar = 0.0;
	Scalar y = mag(unit(fN[f]) & (cC[c1] - cC[c2]));
	if(wallModel == STANDARD)
		ustar = low.getUstar(nu,mag(U[c1]),y);

	/* calculate eddy viscosity*/
	Scalar yp = (ustar * y) / nu;
	Scalar up = low.getUp(ustar,nu,yp);                                      	
	eddy_mu[c1] = (rho * nu) * (yp / up - 1);
}


