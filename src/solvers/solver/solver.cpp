#include "field.h"
#include "turbulence.h"
#include "mp.h"
#include "system.h"
#include "solve.h"

using namespace std;

/*general properties*/
namespace GENERAL {
	Scalar density = 1.2041;
	Scalar viscosity = 1.461e-5;
	Scalar Pr = 0.9;
	Scalar Prt = 0.7;
	Scalar beta = 3.33e-3;
	Scalar T0 = 300;
	Scalar P0 = 10000;
	
	void enroll(Util::ParamList& params) {
		params.enroll("rho", &density);
		params.enroll("viscosity", &viscosity);
		params.enroll("Pr", &Pr);
		params.enroll("Prt", &Prt);
		params.enroll("beta", &beta);
		params.enroll("T0", &T0);
		params.enroll("P0", &P0);
	}
}

/*solvers*/
void piso(istream&);
void diffusion(istream&);
void convection(istream&);
void potential(istream&);
void transport(istream&);
void walldist(istream&);
/**
 \verbatim
 Main application entry point for different solvers.
 \endverbatim
*/
int main(int argc, char* argv[]) {

	/*message passing object*/
	MP mp(argc, argv);
	if(!strcmp(argv[1],"-h")) {
		std::cout << "Usage:\n"
				  << "  ./solver <inputfile>\n"
				  << "Options:\n"
				  << "  -h          --  Display this message\n\n";
		return 0;
	} 
	ifstream input(argv[1]);

	/*main options*/
	Util::ParamList params("general");
	string sname;
	params.enroll("solver", &sname);
	params.enroll("mesh", &Mesh::gMeshName);
	Mesh::enroll(params);
	GENERAL::enroll(params);
	params.read(input);

	/*Mesh*/
	if (mp.n_hosts > 1) {
		stringstream s;
		s << Mesh::gMeshName << mp.host_id;
		if (!System::cd(s.str()))
			return 1;
	}
	Mesh::LoadMesh();
	atexit(Util::cleanup);

	/*call solver*/
	if (!Util::compare(sname, "piso")) {
		piso(input);
	} else if (!Util::compare(sname, "diffusion")) {
		diffusion(input);
	} else if (!Util::compare(sname, "convection")) {
		convection(input);
	} else if (!Util::compare(sname, "transport")) {
		transport(input);
	} else if (!Util::compare(sname, "potential")) {
		potential(input);
	} else if (!Util::compare(sname, "walldist")) {
		walldist(input);
	}

	return 0;
}
/**
 Iteration object that does common book keeping stuff
 for all solvers.
*/
class Iteration {
private:
	Int starti;
	Int endi;
	Int i;
	Int n_deferred;
	Int idf;
public:
	Iteration() {
		Int step = Controls::start_step / Controls::write_interval;
		starti = Controls::write_interval * step + 1;
		endi = Controls::end_step;
		n_deferred = Controls::n_deferred;
		i = starti;
		idf = 0;

		Mesh::read_fields(step);
		Mesh::getProbeCells(Mesh::probeCells);
		forEachCellField (initTimeSeries());
	}
	bool start() {
		return (i == starti);
	}
	bool end() {
		if(i > endi)
			return true;
		/*iteration number*/
		if(MP::host_id == 0) {
			if(Controls::state == Controls::STEADY)
				MP::printH("Step %d\n",i);
			else
				MP::printH("Time %f\n",i * Controls::dt);
		}
		return false;
	}
	void next() {
		idf++;
		if(idf <= n_deferred)
			return;

		/*update time series*/
		forEachCellField(updateTimeSeries(i));

		/*write result to file*/
		if((i % Controls::write_interval) == 0) {
			Int step = i / Controls::write_interval;
			Mesh::write_fields(step);
		}

		/*increment*/
		i++;
	}
	~Iteration() {
	}
	static Int get_start() {
		return Controls::start_step / Controls::write_interval;
	}
	static Int get_end() {
		return Controls::end_step / Controls::write_interval;
	}
};
/**
 \verbatim
 Navier stokes solver using PISO algorithm
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 References:
	Hrvoje Jasak, "Error analysis and estimation of FVM with 
	applications to fluid flow".
 Description:
    The PISO algorithm is used to solve NS equations on collocated grids 
	using Rhie-Chow interpolation to avoid wiggles in pressure field.

	Prediction
	~~~~~~~~~~
	Discretize and solve the momenum equation with current values of pressure. 
	The velocities obtained will not satisfy continuity unless exact pressure 
	happened to be specified. 

	Correction
	~~~~~~~~~~
	Step 1) 
	  Determine velocity with all terms included except pressure gradient source contribution.
	      ap * Up = H(U) - grad(p)
		  Up = H(U) / ap - grad(p) / ap
      Droping grad(p) term:
          Ua = H(U) / ap
          Up = Ua - grad(p) / ap
	  One jacobi sweep is done to find Ua.
    Step 2)
      Solve poisson pressure equation to satisfy continuity with fluxes calculated 
	  from interpolated Ua.
	      div(Up) = 0
		  div(1/ap * grad(p)) = div(H(U)/ap)
		  lap(p,1/ap) = div(Ua)
    Step 3)
	  Correct the velocity with gradient of newly found pressure
	      U = Ua - grad(p) / ap
    These steps are repeated two or more times for transient solutions.
	For steady state problems once is enough.
	\endverbatim
*/
void piso(istream& input) {
	/*Solver specific parameters*/
	Scalar& rho = GENERAL::density;
	Scalar& nu = GENERAL::viscosity;
	Scalar velocity_UR = Scalar(0.8);
	Scalar pressure_UR = Scalar(0.5);
	Scalar t_UR = Scalar(0.8);
	Int n_PISO = 1;
	Int n_ORTHO = 0;
	/*Include buoyancy?*/
	enum BOUYANCY {
		NONE, BOUSSINESQ_T1, BOUSSINESQ_T2, 
		BOUSSINESQ_THETA1, BOUSSINESQ_THETA2,
	};
	BOUYANCY buoyancy = NONE;
	/*piso options*/
	Util::ParamList params("piso");
	Util::Option* op = new Util::Option(&buoyancy, 5, 
			"NONE", "BOUSSINESQ_T1","BOUSSINESQ_T2",
			"BOUSSINESQ_THETA1","BOUSSINESQ_THETA2");
	params.enroll("buoyancy", op);
	params.enroll("velocity_UR", &velocity_UR);
	params.enroll("pressure_UR", &pressure_UR);
	params.enroll("t_UR", &t_UR);
	params.enroll("n_PISO", &n_PISO);
	params.enroll("n_ORTHO", &n_ORTHO);

	VectorCellField U("U", READWRITE);
	ScalarCellField p("p", READWRITE);
	ScalarCellField T(false);

	/*turbulence model*/
	ScalarFacetField F;
	Turbulence_Model::RegisterTable(params);
	params.read(input);
	Turbulence_Model* turb = Turbulence_Model::Select(U, F, rho, nu);
	turb->enroll();

	/*read parameters*/
	Util::read_params(input);

	/*temperature*/
	if(buoyancy != NONE)
		T.construct("T",READWRITE);

	/*wall distance*/
	if (turb->needWallDist())
		Mesh::calc_walldist(Iteration::get_start());

	/*Calculate for each time step*/
	Iteration it;
	ScalarCellField po = p;
	VectorCellField gP = -grad(p);
	F = flx(rho * U);

	for (; !it.end(); it.next()) {
		/*
		 * Prediction
		 */
		VectorCellMatrix M;
		{
			VectorCellField Sc = turb->getTurbSource();
			ScalarCellField eddy_mu = turb->getTurbVisc();
			/* Add buoyancy in two ways
			 *  1. pm = p - rho*g*h
			 *  2. pm = p - rho_m*g*h
			 */
			if (buoyancy != NONE) {
				Scalar beta;
				if(buoyancy <= BOUSSINESQ_T2) 
					beta = GENERAL::beta;
				else
					beta = 1 / GENERAL::T0;
				if (buoyancy == BOUSSINESQ_T1 || buoyancy == BOUSSINESQ_THETA1) {  
					ScalarCellField rhok = rho * (0 - beta * (T - GENERAL::T0));
					Sc += (rhok * VectorCellField(Controls::gravity));
				} else if(buoyancy == BOUSSINESQ_T2 || buoyancy == BOUSSINESQ_THETA2) {
					ScalarCellField gz = dot(Mesh::cC,VectorCellField(Controls::gravity));
					Sc += gz * (rho * beta) * gradi(T);
				}
			}
			/*momentum prediction*/
			{
				ScalarFacetField mu = cds(eddy_mu) + rho * nu;
				M = transport(U, U, F, mu, rho, velocity_UR, Sc, Scalar(0));
				Solve(M == gP);
			}
			/*energy predicition*/
			if (buoyancy != NONE) {
				ScalarFacetField mu = cds(eddy_mu) / GENERAL::Prt + (rho * nu) / GENERAL::Pr;
				ScalarCellMatrix Mt = transport(T, U, F, mu, rho, t_UR);
				Solve(Mt);
				T = max(T, Constants::MachineEpsilon);
			}
		}

		/*
		 * Correction
		 */
		ScalarCellField api = fillBCs(1.0 / M.ap, true);
		ScalarCellField rmu = rho * api * Mesh::cV;

		/*PISO loop*/
		for (Int j = 0; j < n_PISO; j++) {
			/* Ua = H(U) / ap*/
			U = getRHS(M) * api;
			updateExplicitBCs(U, true);

			/*solve pressure poisson equation to satisfy continuity*/
			{
				ScalarCellField rhs = div(rho * U);
				for (Int k = 0; k <= n_ORTHO; k++)
					Solve(lap(p, rmu) += rhs);
			}

			/*explicit velocity correction : add pressure contribution*/
			gP = -grad(p);
			U -= gP * api;
			updateExplicitBCs(U, true);
		}

		/*update fluctuations*/
		updateExplicitBCs(U, true, true);
		F = flx(rho * U);

		/*solve turbulence transport equations*/
		turb->solve();

		/*explicitly under relax pressure*/
		if (Controls::state == Controls::STEADY) {
			p.Relax(po, pressure_UR);
			gP = -grad(p);
			po = p;
		}
	}

	/*write calculated turbulence fields*/
	if (turb->writeStress) {
		ScalarCellField K("Ksgs", WRITE);
		STensorCellField R("Rsgs", WRITE);
		STensorCellField V("Vsgs", WRITE);
		K = turb->getK();
		R = turb->getReynoldsStress();
		V = turb->getViscousStress();
		Mesh::write_fields(Iteration::get_end());
	}
}
/**
 \verbatim
 Diffusion solver
 ~~~~~~~~~~~~~~~~
 Solver for pdes of parabolic heat equation type:
       dT/dt = lap(T,DT)
 \endverbatim
*/
void diffusion(istream& input) {
	/*Solver specific parameters*/
	Scalar& rho = GENERAL::density;
	Scalar DT = Scalar(1);
	Scalar t_UR = Scalar(1);

	/*diffusion*/
	Util::ParamList params("diffusion");
	params.enroll("DT", &DT);
	params.enroll("t_UR", &t_UR);

	ScalarCellField T("T", READWRITE);

	/*read parameters*/
	Util::read_params(input);

	/*Calculate for each time step*/
	ScalarFacetField mu = rho * DT;

	for (Iteration it; !it.end(); it.next()) {
		ScalarCellMatrix M;
		M = diffusion(T, mu, rho, t_UR);
		Solve(M);
	}
}
/**
  \verbatim
  Convection solver
  ~~~~~~~~~~~~~~~~~~~~~~~~~~
  Given a flow field (U), the solver determines the distribution of a 
  scalar by convection.
     dT/dt + div(T,F,0) = 0
  \endverbatim
 */
void convection(istream& input) {
	/*Solver specific parameters*/
	Scalar& rho = GENERAL::density;
	Scalar t_UR = Scalar(1);

	/*transport*/
	Util::ParamList params("convection");
	params.enroll("t_UR", &t_UR);

	VectorCellField U("U", READWRITE);
	ScalarCellField T("T", READWRITE);

	/*read parameters*/
	Util::read_params(input);

	/*Calculate for each time step*/
	Iteration it;
	ScalarFacetField F = flx(rho * U), mu = Scalar(0);
	for (; !it.end(); it.next()) {
		ScalarCellMatrix M;
		M = convection(T, U, F, mu, rho, t_UR);
		Solve(M);
	}
}
/**
  \verbatim
  Transport equation solver
  ~~~~~~~~~~~~~~~~~~~~~~~~~
  Given a flow field (U), the solver determines the distribution of a 
  scalar by convection and diffusion.
     dT/dt + div(T,F,DT) = lap(T,DT)
  \endverbatim
 */
void transport(istream& input) {
	/*Solver specific parameters*/
	Scalar& rho = GENERAL::density;
	Scalar DT = Scalar(1.0e-4);
	Scalar t_UR = Scalar(1);

	/*transport*/
	Util::ParamList params("transport");
	params.enroll("t_UR", &t_UR);
	params.enroll("DT", &DT);

	VectorCellField U("U", READWRITE);
	ScalarCellField T("T", READWRITE);

	/*read parameters*/
	Util::read_params(input);

	/*Calculate for each time step*/
	Iteration it;
	ScalarFacetField F = flx(rho * U), mu = rho * DT;
	for (; !it.end(); it.next()) {
		ScalarCellMatrix M;
		M = transport(T, U, F, mu, rho, t_UR);
		Solve(M);
	}
}
/**
  \verbatim
    Potential flow solver
	~~~~~~~~~~~~~~~~~~~~~
	In potential flow the velocity field is irrotational (vorticity = curl(U) = 0).
	This assumption fails for boundary layers and wakes that exhibit strong vorticity,
	but it can still be used to initialize the flow field for further simulations.

	For incompressible flow
	       div(U) = 0
    Velocity is the gradient of velocity potential phi
	       U = grad(phi)
		   div(grad(phi)) = 0
		   lap(phi) = 0
    phi is pressure for this solver. The initial flow field will inevitably not satisfy 
	continuity due to imposed boundary conditons. Therefore we solve a pressure poisson 
	equation and then correct the velocity with the gradient of p.
	       lap(p) = div(U)
		   U -= grad(p)
  \endverbatim
*/
void potential(istream& input) {
	/*Solver specific parameters*/
	Int n_ORTHO = 0;

	/*potential*/
	Util::ParamList params("potential");
	params.enroll("n_ORTHO", &n_ORTHO);

	VectorCellField U("U", READWRITE);
	ScalarCellField p("p", READ);

	/*read parameters*/
	Util::read_params(input);

	/*set internal field to zero*/
	for (Int i = 0; i < Mesh::gBCSfield; i++) {
		U[i] = Vector(0, 0, 0);
		p[i] = Scalar(0);
	}
	updateExplicitBCs(U, true);
	updateExplicitBCs(p, true);

	for (Iteration it; it.start(); it.next()) {
		/*solve potential equation*/
		ScalarCellField divU = div(U);
		ScalarFacetField one = Scalar(1);
		for (Int k = 0; k <= n_ORTHO; k++)
			Solve(lap(p, one) == divU);

		/*correct velocity*/
		U -= gradi(p);
		updateExplicitBCs(U, true);
	}
}
/**
 \verbatim
 Wall distance
 ~~~~~~~~~~~~~
	Reference:
	   D.B.Spalding, Calculation of turbulent heat transfer in cluttered spaces
    Description:
	   Poisson equation is solved to get approximate nearest wall distance.
	         lap(phi,1) = -cV
	   The boundary conditions are phi=0 at walls, and grad(phi) = 0 elsewhere.
 \endverbatim
*/
void walldist(istream& input) {
	/*Solver specific parameters*/
	Int n_ORTHO = 0;

	/*walldist options*/
	Util::ParamList params("walldist");
	params.enroll("n_ORTHO", &n_ORTHO);
	Util::read_params(input);

	/*solve*/
	Mesh::calc_walldist(Iteration::get_start(), n_ORTHO);
}
void Mesh::calc_walldist(Int step, Int n_ORTHO) {
	ScalarCellField& phi = yWall;
	/*poisson equation*/
	{
		ScalarFacetField one = Scalar(1);
		for (Int k = 0; k <= n_ORTHO; k++)
			Solve(lap(phi, one) == -cV);
	}
	/*wall distance*/
	{
		VectorCellField g = gradi(phi);
		yWall = sqrt((g & g) + 2 * phi) - mag(g);
	}
	/*write it*/
	yWall.write(step);
}
