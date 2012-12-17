#ifndef __TURBULENCE_H
#define __TURBULENCE_H

#include "field.h"
#include "solve.h"
/*
	 Description of RANS turbulence models
	 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 Navier Stokes without source term:
	     d(rho*u)/dt + div(rho*uu) = -grad(p) + div(mu*gu)
	 RANS:
	     d(rho*U)/dt + div(rho*UU) + div(rho*u'u') = -grad(P) + div(mu*gU)
	     d(rho*U)/dt + div(rho*UU) = -grad(P) + div(mu*gU) - div(rho*u'u')
	     d(rho*U)/dt + div(rho*UU) = -grad(P) + div(V + R)
	 where Viscous (V) and Reynolds (R) stress tensors are
		 V =  mu*gU
	     R = -rho*u'u'
	 Boussinesq model for R:
	     Traceless(R) = 2 * emu * Traceless(S) 
		 where S = (gU + gUt) / 2
	     R - R_ii/3 = 2 * emu * (S - S_ii/3)
	     R = 2 * emu * (S - S_ii/3) + R_ii/3
	       = 2 * emu * ((gU + gUt)/2 - gU_ii/3) + R_ii/3
	       = emu * gU + emu * (gUt - 2/3*gUt_ii) + R_ii/3
	       = emu * gU + emu * dev(gUt,2) - 2/3*rho*k*I
	 Viscous and Reynolds stress together:
	     V + R = {mu * gU} + {emu * gU + emu * dev(gUt,2) - 2/3*rho*k*I}
	           = (mu + emu) * gU + emu * dev(gUt,2) - 2/3*rho*k*I
	           = ( eff_mu ) * gU + emu * dev(gUt,2) - 2/3*rho*k*I
	 Volume integrated V+R i.e force:
	     div(V + R) = div(eff_mu*gU) + div(emu * dev(gUt,2)) - div(2/3*rho*k*I)
	                     Implicit           Explicit          Absored in pressure 
	                                                          p_m = p + 2/3*k*rho
	 Final RANS equation after substituting div(V+R):
		 d(rho*U)/dt + div(rho*UU) = -grad(P) + div(V + R)
	     d(rho*U)/dt + div(rho*UU) = -grad(P_m) + div(eff_mu*gU) + div(emu * dev(gUt,2))
	 Since the k term is absorbed into the pressure gradient, we only need models for
	 turbulent diffusivity emu.
*/

/*
     Base turbulence model:
         This default class has no turbulence model so it is a laminar solver. 
         Only the viscous stress V is added to the NS equations. Turbulence models 
         derived from this class add a model for Reynold's stress R usually by solving 
         'turbulence transport' equations.
*/
struct Turbulence_Model {

	VectorCellField& U;
	ScalarFacetField& F;
	Scalar& rho;
	Scalar& nu;
	bool& Steady;

	Util::ParamList params;
	/*constructor*/
	Turbulence_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
		U(tU),
		F(tF),
		rho(trho),
		nu(tnu),
		Steady(tSteady),
		params("turbulence")
	{
	}
	/*overridable functions*/
	virtual void enroll() {};
	virtual void solve() {};
	virtual void addTurbulentStress(VectorMeshMatrix& M) {
		ScalarFacetField mu = rho * nu;
		M -= lap(U,mu);
	};
	/* V */
	STensorCellField getViscousStress() {
		STensorCellField V = 2 * rho * nu * sym(grad(U));
		return V;
	}
	/* R */
	virtual STensorCellField getReynoldsStress() {
		return STensor(0);
	}
	/* TKE */
	virtual ScalarCellField getK() {
		return Scalar(0);
	}
};
/*
 * Eddy viscosity models based on Boussinesq's assumption
 * that the action of Reynolds and Viscous stress are similar.
 */
struct EddyViscosity_Model : public Turbulence_Model {
	ScalarCellField eddy_mu; 
	Int modelType;

	/*constructor*/
	EddyViscosity_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
		Turbulence_Model(tU,tF,trho,tnu,tSteady),
		modelType(0)
	{
	}
	/*Register options*/
	virtual void enroll() {
		using namespace Util;
		params.enroll("modelType",&modelType);
		Turbulence_Model::enroll();
	}
	/*eddy_mu*/
	virtual void calcEddyViscosity(const TensorCellField& gradU) = 0;

	/* V + R */
	virtual void addTurbulentStress(VectorMeshMatrix& M) {
		TensorCellField gradU = grad(U);
		calcEddyViscosity(gradU);

		ScalarCellField eff_mu = eddy_mu + rho * nu;
		M -= lap(U,eff_mu);
		M -= div(eddy_mu * dev(trn(gradU),2));
	};
	/* R */
	virtual STensorCellField getReynoldsStress() {
		STensorCellField R = 2 * eddy_mu * dev(sym(grad(U))) - 
			     STensorCellField(Constants::I_ST) * (2 * rho * getK() / 3);
		return R;
	}
	/* S2 */
	ScalarCellField getS2(const TensorCellField& gradU) {
		ScalarCellField magS;
		if(modelType == 0) {
			STensorCellField S = sym(gradU);
			magS = S & S;
		} else if(modelType == 1) {
			TensorCellField O = skw(gradU);
			magS = O & O;
		} else {
			STensorCellField S = sym(gradU);
			TensorCellField O = skw(gradU);
			magS = sqrt((S & S) * (O & O));
		}
		return (2 * magS);
	}
};
/*
 * Model for flow close to the wall (Law of the wall).
 *   1 -> Viscous layer
 *   2 -> Buffer layer
 *   3 -> Log-law layer
 * The wall function model is modified for rough surfaces 
 * using Cebecci and Bradshaw formulae.
 */
struct LawOfWall {
	Scalar E;
	Scalar kappa;
	Scalar ks;
	Scalar cks;

	Scalar yLog;

	LawOfWall() : 
		E(9.793),
		kappa(0.4187),
		ks(0.48),
		cks(0.5)
	{
		init();
	}
	void enroll(Util::ParamList params) {
		params.enroll("E",&E);
		params.enroll("Kappa",&kappa);
		params.enroll("Ks",&ks);
		params.enroll("Cks",&cks);
	}
	void init() {
		yLog = 11.3f;
		for(Int i = 0;i < 20;i++)
			yLog = log(E * yLog) / kappa;
	}
	Scalar getDB(Scalar ustar,Scalar nu) {
		Scalar dB;
		Scalar ksPlus = (ustar * ks) / nu;
		if(ksPlus < 2.25) {
			dB = 0;
		} else if(ksPlus < 90) {
			dB = (1 / kappa) * log((ksPlus - 2.25) / 87.75 + cks * ksPlus)
				             * sin(0.4258 * (log(ksPlus) - 0.811));
		} else {
			dB = (1 / kappa) * log(1 + cks * ksPlus);
		}
		return dB;
	}
};
/*
 * Base two equation K-X turbulence model
 */ 
struct KX_Model : public EddyViscosity_Model {
	/*model coefficients*/
	Scalar Cmu;
	Scalar SigmaK;
	Scalar SigmaX;
	Scalar C1x;
	Scalar C2x;

	Scalar k_UR;
	Scalar x_UR;

	/*turbulence fields*/
	ScalarCellField k;         
	ScalarCellField x;        
	ScalarCellField G;       

	/*constructor*/
	KX_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady,const char* xname) :
		EddyViscosity_Model(tU,tF,trho,tnu,tSteady),
		k_UR(0.7),
		x_UR(0.7),
		k("k",READWRITE),
		x(xname,READWRITE)
	{
	}
	/*TKE*/
	virtual ScalarCellField getK() { return k; }
	/*Register options*/
	virtual void enroll() {
		using namespace Util;
		params.enroll("k_UR",&k_UR);
		params.enroll("x_UR",&x_UR);
		EddyViscosity_Model::enroll();
	}
	/* eddy viscosity*/
	virtual void calcEddyViscosity(const TensorCellField& gradU) {
		calcEddyMu();
		G = getS2(gradU) * eddy_mu;
		addWallContribution();
	}
	/* k-x model specific over-ridables*/
	virtual void calcEddyMu() = 0;
	virtual Scalar calcX(Scalar ustar,Scalar kappa,Scalar y) = 0;
	virtual Scalar getCmu(Int i) { 
		return Cmu; 
	}
	/* Calculate wall turbulence generation */
	void addWallContribution() {
		using namespace Mesh;
		BasicBCondition* bbc;
		for(Int d = 0;d < AllBConditions.size();d++) {
			bbc = AllBConditions[d];
			if(bbc->isWall && (bbc->fIndex == U.fIndex)) {
				IntVector& wall_faces = *bbc->bdry;
				LawOfWall law;
				Scalar E = law.E;
				Scalar kappa = law.kappa;
				Scalar yLog = law.yLog;

				if(wall_faces.size()) {
					Vector dv;
					Scalar y,ustar,yp,up,gU;
					Int f,c1,c2;

					/*calc eddy viscosity*/
					for(Int i = 0;i < wall_faces.size();i++) {
						f = wall_faces[i];
						c1 = gFO[f];
						c2 = gFN[f];

						/*viscous and log-law layer*/
						y = mag(unit(fN[f]) & (cC[c1] - cC[c2]));
						ustar = pow(getCmu(c1),Scalar(0.25)) * sqrt(k[c1]);

						yp = (ustar * y) / nu;
						if(yp > yLog)  up = log(E * yp) / kappa - law.getDB(ustar,nu);  
						else           up = yp;                                         
						eddy_mu[c1] = rho * nu * (yp / up);

						/*wall dissipation and generation*/
						x[c1] = calcX(ustar,kappa,y);
						if(yp > yLog) {
							gU = mag((U[c2] - U[c1]) / y);
							G[c1] = eddy_mu[c1] * gU * ustar / (kappa * y);
						}
						/*set boundary values*/
						x[c2] = x[c1];
						G[c2] = 0;
					}
				}
			}
		}
	}
	/*end*/
};

#endif
