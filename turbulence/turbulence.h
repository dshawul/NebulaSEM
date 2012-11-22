#ifndef __TURBULENCE_H
#define __TURBULENCE_H

#include "field.h"
#include "solve.h"

/*
 * Base turbulence model
 *  -> This default class has no turbulence model so it is a laminar solver. 
 *     It adds only the viscous stresses. Turbulence models derived from this class 
 *     add a model for Reynold's stress term.
 */
struct Turbulence_Model {

	VectorCellField& U;
	ScalarFacetField& F;
	Scalar& rho;
	Scalar& nu;
	bool& Steady;

	/*constructor*/
	Turbulence_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady) :
		U(tU),
		F(tF),
		rho(trho),
		nu(tnu),
		Steady(tSteady)
	{
	}
	/*others*/
	virtual void enroll() {};
	virtual void solve() {};
	virtual void addTurbulentStress(VectorMeshMatrix& M) {
		/*add viscous stress*/
		ScalarFacetField mu = rho * nu;
		M -= lap(U,mu);
	};
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
	Scalar yLog;

	Scalar ks;
	Scalar cks;

	LawOfWall() {
		init();
	}
	LawOfWall(Scalar Et, Scalar kappat) {
		init(Et,kappat);
	}
	LawOfWall(Scalar Et, Scalar kappat,Scalar kst,Scalar ckst){
		init(Et,kappat,kst,ckst);
	}

	void init(Scalar Et = 9.793, Scalar kappat = 0.4187,
	          Scalar kst = 0.1,Scalar ckst = 0.5) {
		E = Et;
		kappa = kappat;
		ks = kst;
		cks = ckst;

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

struct WallFunction {
	IntVector* faces;
	LawOfWall func;
};
/*
 * Base two equation K-X turbulence model
 */ 
struct KX_Model : public Turbulence_Model {
	/*model coefficients*/
	Scalar Cmu;
	Scalar SigmaK;
	Scalar SigmaX;
	Scalar C1x;
	Scalar C2x;

	Scalar UR;
	Int    katoLaunder;

	/*turbulence fields*/
	ScalarCellField k;         
	ScalarCellField x;        
	ScalarCellField eddy_mu;  
	ScalarCellField G;       

	/*wall functions*/
	std::vector<WallFunction> wall_functions;

	/*constructor*/
	KX_Model(VectorCellField& tU,ScalarFacetField& tF,Scalar& trho,Scalar& tnu,bool& tSteady,const char* xname) :
		Turbulence_Model(tU,tF,trho,tnu,tSteady),
		UR(0.7),
		katoLaunder(0),
		k("k",READWRITE),
		x(xname,READWRITE)
	{
		/*
		 * Any patch with a name that has 'WALL' in it will be
		 * considered for application of the law of the wall
		 */
		using namespace Mesh;
		for(Boundaries::iterator it = gBoundaries.begin();it != gBoundaries.end();++it) {
			if(it->first.find("WALL") != std::string::npos) {
				WallFunction wall_func;
				wall_func.faces = &it->second;
				wall_functions.push_back(wall_func);
			}
		}
	}
	/*Register options*/
	virtual void enroll() {
		using namespace Util;
		ScalarParams::enroll("turbulence_UR",&UR);
		Util::Option* op = new Util::Option(&katoLaunder,2,"NO","YES");
		Util::OptionParams::enroll("katoLaunder",op);
	}
	/*
	 * Kato-Launder modifaction to reduce over-production of turbulence.
	 *   Replace one of the strain rates S with the vorticity O
	 */
	void calcG(const TensorCellField& gradU) {
		if(katoLaunder) {
            STensorCellField S = sym(gradU);
			TensorCellField O = skw(gradU);
			G = sqrt((S & S) * (O & O)) * 2 * eddy_mu;
		} else {
			STensorCellField S = sym(gradU);
			G = (S & S) * 2 * eddy_mu;
		}
	}
	void addTurbulentStress(VectorMeshMatrix& M) {
		using namespace Mesh;

		TensorCellField gradU = grad(U);
		calcEddyMu();
		calcG(gradU);
		addWallContribution();
		/*
		 * Viscous and Reynolds stress
		 *   Reynolds stress = R = -mu * U'U'
		 *   Viscous stress  = V =  mu * gU
		 * Boussinesq approximation:
		 *   Traceless(R) = 2 * emu * Traceless(S)
		 *   R - R_ii/3 = 2 * emu * (S - S_ii/3)
		 *   R = 2 * emu * (S - S_ii/3) + R_ii/3
		 *     = 2 * emu * ((gU + gUt)/2 - gU_ii/3) + R_ii/3
		 *     = emu * gU + emu * (gUt - 2/3*gUt_ii) + R_ii/3
		 *     = emu * gU + emu * dev(gUt,2) - 2/3*rho*k*I
		 * Viscous and Reynolds stress
		 *   V + R = {mu * gU} + {emu * gU + emu * dev(gUt,2) - 2/3*rho*k*I
		 *         = (mu + emu) * gU + emu * dev(gUt,2) - 2/3*rho*k*I
		 * Volume integrated V and R
		 *   div(V + R) = div(eff_mu*gU) + div(emu * dev(gUt,2)) - div(2/3*rho*k*I)
		 *                   Implicit             Explicit       Absored in pressure 
		 *                                                       p_m = p + 2/3*k*rho
		 * Navier Stokes.
		 *     du/dt + div(rho*u*u) = -grad(p) + div(mu*gU)
		 * RANS
		 *     dU/dt + div(rho*U*U) = -grad(P_m) + div(eff_mu*gU) + div(emu * dev(gUt,2))
		 */
		{
			ScalarCellField eff_mu = eddy_mu + rho * nu;
			M -= lap(U,eff_mu);
			M -= sum(div(eddy_mu * dev(trn(gradU),2)));
		}
	}
	/* Reynolds stress tensor */
	STensorCellField getReynoldsStress() {
		STensorCellField R = 2 * eddy_mu * dev(sym(grad(U))) - 
			                STensorCellField(Constants::I_ST) * (2 * rho * k / 3);
		return R;
	}
	/* virtual functions*/
	virtual void calcEddyMu() = 0;
	virtual Scalar calcX(Scalar ustar,Scalar kappa,Scalar y) = 0;
	virtual Scalar getCmu(Int i) { return Cmu; }
	/* Calculate wall turbulence generation*/
	void addWallContribution() {
		using namespace Mesh;
		for(Int d = 0;d < wall_functions.size();d++) {
			WallFunction& wall_func = wall_functions[d];
			IntVector& wall_faces = *wall_func.faces;
			LawOfWall& law = wall_func.func;
			Scalar E = law.E;
			Scalar kappa = law.kappa;
			Scalar yLog = law.yLog;

			if(wall_faces.size()) {
				Vector dv;
				Scalar y,ustar,yp,up,gU;
				Int f,c1,c2;

				/*accumulate turbulence generation*/
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

					/*dissipation*/
					x[c1] = calcX(ustar,kappa,y);
					if(yp > yLog) {
						gU = mag((U[c2] - U[c1]) / y);
						G[c1] = eddy_mu[c1] * gU * ustar / (kappa * y);
					}
					/*assume neumann*/
					x[c2] = x[c1];
					G[c2] = G[c1];
				}
			}
		}
	}
};

#endif
