#include "field.h"
#include "turbulence.h"
#include "mp.h"
#include "system.h"
#include "solve.h"

using namespace std;

/*general properties*/
namespace GENERAL {
	Scalar density = 1.177;
	Scalar viscosity = 1.568e-5;
	Scalar Pr = 0.9;
	Scalar Prt = 0.7;
	Scalar beta = 3.33e-3;
	Scalar T0 = 300;
	Scalar P0 = 101325;
	Scalar cp = 1004.67;
	Scalar cv = 715.5;
	
	void enroll(Util::ParamList& params) {
		params.enroll("rho", &density);
		params.enroll("viscosity", &viscosity);
		params.enroll("Pr", &Pr);
		params.enroll("Prt", &Prt);
		params.enroll("beta", &beta);
		params.enroll("T0", &T0);
		params.enroll("P0", &P0);
		params.enroll("cp", &cp);
		params.enroll("cv", &cv);
	}
}

/*solvers*/
void piso(istream&);
void diffusion(istream&);
void convection(istream&);
void potential(istream&);
void transport(istream&);
void walldist(istream&);
void euler(istream&);
void wave(istream&);
/**
 \verbatim
 Main application entry point for different solvers.
 \endverbatim
*/
int main(int argc, char* argv[]) {

	/*message passing object*/
	MP mp(argc, argv);
	MP::printOn = (MP::host_id == 0);
	if(!strcmp(argv[1],"-h")) {
		std::cout << "Usage:\n"
				  << "  ./solver <inputfile>\n"
				  << "Options:\n"
				  << "  -h          --  Display this message\n\n";
		return 0;
	} 
	ifstream input(argv[1]);

	/*General options*/
	string sname;
	{
		Util::ParamList params("general");
		params.enroll("solver", &sname);
		params.enroll("mesh", &Mesh::gMeshName);
		Mesh::enroll(params);
		GENERAL::enroll(params);
		params.read(input);
	}
	/*AMR options*/
	{
		Util::ParamList params("refinement");
		Controls::enrollRefine(params);
		params.read(input);
	}
	/*Decompose options*/
	{
		Util::ParamList params("decomposition");
		Controls::enrollDecompose(params);
		params.read(input); 
	}
	/*Fields*/
	{
		Util::ParamList params("prepare");
		params.enroll("fields",&BaseField::fieldNames);
		params.read(input);
	}
	/*cleanup*/
	atexit(MP::cleanup);

	/*call solver*/
	if (!Util::compare(sname, "piso")) {
		piso(input);
	} else if (!Util::compare(sname, "euler")) {
		euler(input);
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
	} else if (!Util::compare(sname, "wave")) {
		wave(input);
	}
	
#ifdef _DEBUG
	/*print memory usage*/
	std::cout << "====================================" << std::endl;
	std::cout << "Memory Usage:" << std::endl;
	forEachCellField(printUsage());
	forEachFacetField(printUsage());
	forEachVertexField(printUsage());
	std::cout << "====================================" << std::endl;
#endif	
	
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
	Iteration(Int step) {
		starti = Controls::write_interval * step + 1;
		endi = Controls::write_interval * (step + Controls::amr_step);
		if(endi > Controls::end_step) endi = Controls::end_step;
		n_deferred = Controls::n_deferred;
		i = starti;
		idf = 0;
		if(MP::printOn)
			cout << "--------------------------------------------\n";
		Mesh::read_fields(step);
		Mesh::getProbeCells(Mesh::probeCells);
		forEachCellField (initTimeSeries());
		if(MP::printOn) {
			cout << "--------------------------------------------\n";
			MP::printH("Starting iterations.\n");
		}
	}
	bool start() {
		return (i == starti);
	}
	bool end() {
		if(i > endi)
			return true;
		/*iteration number*/
		if(MP::printOn && idf == 0) {
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
		idf = 0;
		
		/*set output printing*/
		MP::printOn = (MP::host_id == 0 && 
			(MP::hasElapsed(Controls::print_time) || i == Controls::end_step - 1)); 
		
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
};
/**
 Iterator for AMR
*/
class AmrIteration {
private:
	Int starti;
	Int endi;
	Int i;
public:
	AmrIteration() {
		starti = Controls::start_step / Controls::write_interval;
		endi = Controls::end_step / Controls::write_interval;
		if(!Controls::amr_step)
			Controls::amr_step = endi;
		i = starti;

		if (MP::n_hosts > 1) {
			/*decompose*/
			if (MP::host_id == 0)
				Prepare::decomposeMesh(i);
			/*wait*/
			MP::barrier();
			/*change directory*/
			stringstream s;
			s << Mesh::gMeshName << MP::host_id;
			System::cd(s.str());
		}
		Mesh::LoadMesh(i);
	}
	bool start() {
		return (i == starti);
	}
	bool end() {
		if(i >= endi)
			return true;
		return false;
	}
	void next() {
		i += Controls::amr_step;
		if(i < endi) {
			if (MP::host_id == 0) {
				if(MP::n_hosts > 1) {
					System::cd(MP::workingDir);
					Prepare::mergeFields(i);
				}
				
				if(i == starti + Controls::amr_step)
					Prepare::initRefineThreshold();
				
				Prepare::refineMesh(i);

				if(MP::n_hosts > 1) {

					/*decompose mesh*/
					Prepare::decomposeMesh(i);

					/*change directory*/
					stringstream s;
					s << Mesh::gMeshName << MP::host_id;
					System::cd(s.str());
				}
			}
			MP::barrier();
			Mesh::LoadMesh(i);	
		}
	}
	~AmrIteration() {
	}
	Int get_step() {
		return i;
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
	
	Turbulence_Model::RegisterTable(params);
	params.read(input);
	
	/*AMR iteration*/
	for (AmrIteration ait; !ait.end(); ait.next()) {
		ScalarCellField rho = GENERAL::density;
		VectorCellField U("U", READWRITE);
		ScalarCellField p("p", READWRITE);
		ScalarCellField T(false);

		/*temperature*/
		if(buoyancy != NONE)
			T.construct("T",READWRITE);
		
		/*turbulence model*/
		ScalarFacetField F;
		Turbulence_Model* turb = Turbulence_Model::Select(U, F, rho, nu);

		/*read parameters*/
		if(ait.start()) {
			turb->enroll();
			Util::read_params(input,MP::printOn);
		}

		/*wall distance*/
		if (turb->needWallDist())
			Mesh::calc_walldist(ait.get_step(), 2);

		/*Time loop*/
		Iteration it(ait.get_step());
		ScalarCellField po = p;
		VectorCellField gP = -gradf(p);
		VectorCellField Fc;
	
		Fc = rho * U;
		F = flx(Fc);

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
					ScalarCellField mu = eddy_mu + rho * nu;
					M = transport(U, Fc, F, mu, velocity_UR, Sc, Scalar(0));
					Solve(M == gP);
				}
				/*energy predicition*/
				if (buoyancy != NONE) {
					ScalarCellField mu = eddy_mu / GENERAL::Prt + (rho * nu) / GENERAL::Pr;
					ScalarCellMatrix Mt = transport(T, Fc, F, mu, t_UR);
					Solve(Mt);
					T = max(T, Constants::MachineEpsilon);
				}
			}

			/*
			 * Correction
			 */
			ScalarCellField api = fillBCs(1.0 / M.ap);
			ScalarCellField rmu = rho * api * Mesh::cV;
		
			/*PISO loop*/
			for (Int j = 0; j < n_PISO; j++) {
				/* Ua = H(U) / ap*/
				U = getRHS(M) * api;
				applyExplicitBCs(U, true);
			
				/*solve pressure poisson equation to satisfy continuity*/
				{
					ScalarCellField rhs = divf(rho * U);
					for (Int k = 0; k <= n_ORTHO; k++)
						Solve(lap(p, rmu) += rhs);
				}
			
				/*explicit velocity correction : add pressure contribution*/
				gP = -gradf(p);
				U -= gP * api;
				applyExplicitBCs(U, true);
			}
			/*update fluctuations*/
			applyExplicitBCs(U, true, true);
			Fc = rho * U;
			F = flx(Fc);

			/*solve turbulence transport equations*/
			turb->solve();

			/*explicitly under relax pressure*/
			if (Controls::state == Controls::STEADY) {
				p.Relax(po, pressure_UR);
				gP = -gradf(p);
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
			Mesh::write_fields(ait.get_step());
		}
		
		delete turb;
	}
}
/**
  \verbatim
  Euler equations solver
  ~~~~~~~~~~~~~~~~~~~~~~
  d(rho*(1))/dt + div(1,F,0) = 0
  d(rho*(U))/dt + div(U,F,0) = -grad(p) + rho * gravity
  d(rho*(T))/dt + div(T,F,0) = 0
  \endverbatim
 */
void euler(istream& input) {
	/*Solver specific parameters*/
	Scalar pressure_UR = Scalar(0.5);
	Scalar velocity_UR = Scalar(0.8);
	Scalar t_UR = Scalar(0.8);
	Int buoyancy = 1;
	Int diffusion = 1;

	/*transport*/
	Util::ParamList params("euler");
	params.enroll("velocity_UR", &velocity_UR);
	params.enroll("pressure_UR", &pressure_UR);
	params.enroll("t_UR", &t_UR);
	Util::Option* op;
	op = new Util::BoolOption(&buoyancy);
	params.enroll("buoyancy",op);
	op = new Util::BoolOption(&diffusion);
	params.enroll("diffusion",op);
	
	/*read parameters*/
	Util::read_params(input,MP::printOn);
	
	/*AMR iteration*/
	for (AmrIteration ait; !ait.end(); ait.next()) {
		
		ScalarCellField p("p",READWRITE);
		VectorCellField U("U",READWRITE);
		ScalarCellField T("T",READWRITE);
		ScalarCellField rho("rho",WRITE);
		
		/*Read fields*/
		Iteration it(ait.get_step());
		ScalarFacetField F;
		ScalarCellField mu;

		/*gas constants*/
		Scalar p_factor, p_gamma, R, psi, iPr;
		using namespace GENERAL;
		R = cp - cv;
		p_gamma = cp / cv;
		p_factor = P0 * pow(R / P0, p_gamma);
		psi = 1 / (R * T0);
		iPr = 1 / Pr;

		/*calculate rho*/
		rho = pow(p / p_factor, 1 / p_gamma) / T;
		Mesh::scaleBCs<Scalar>(p,rho,psi);
		rho.write(0);
	
		/*Time loop*/
		for (; !it.end(); it.next()) {
			/*fluxes*/
			F = flx(rho * U);
			/*rho-equation*/
			{
				ScalarCellMatrix M;
				M = convection(rho, U, flx(U), pressure_UR);
				Solve(M);
			}
			/*artificial viscosity*/
			if(diffusion) mu = rho * viscosity;
			else mu = Scalar(0);
			/*T-equation*/
			{
				ScalarCellMatrix M,Mt;
				M = transport(T, rho * U, F, mu * iPr, t_UR, &rho);
				Solve(M);
			}
			/*calculate p*/
			p = p_factor * pow(rho * T, p_gamma);
			/*U-equation*/
			{
				VectorCellField Sc = Vector(0);
				/*buoyancy*/
				if(buoyancy)
					Sc += rho * VectorCellField(Controls::gravity);
				/*solve*/
				VectorCellMatrix M;
				M = transport(U, rho * U, F, mu, velocity_UR, Sc, Scalar(0), &rho);
				Solve(M == -gradf(p));
			}
			/*courant number*/
			if(Controls::state != Controls::STEADY)
				Mesh::calc_courant(U,Controls::dt);
		}
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
	Scalar t_UR = Scalar(1);

	/*transport*/
	Util::ParamList params("convection");
	params.enroll("t_UR", &t_UR);

	/*read parameters*/
	Util::read_params(input,MP::printOn);
	
	/*AMR iteration*/
	for (AmrIteration ait; !ait.end(); ait.next()) {
		
		VectorCellField U("U", READWRITE);
		ScalarCellField T("T", READWRITE);

		/*Time loop*/
		Iteration it(ait.get_step());
		ScalarFacetField F = flx(U);
		for (; !it.end(); it.next()) {
			ScalarCellMatrix M;
			M = convection(T, U, F, t_UR);
			Solve(M);
		}
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
	Scalar DT = Scalar(1);
	Scalar t_UR = Scalar(1);

	/*diffusion*/
	Util::ParamList params("diffusion");
	params.enroll("DT", &DT);
	params.enroll("t_UR", &t_UR);

	/*read parameters*/
	Util::read_params(input,MP::printOn);
	
	/*AMR iteration*/
	for (AmrIteration ait; !ait.end(); ait.next()) {
		
		ScalarCellField T("T", READWRITE);
		
		/*Time loop*/
		ScalarCellField mu = DT;
		for (Iteration it(ait.get_step()); !it.end(); it.next()) {
			ScalarCellMatrix M;
			M = diffusion(T, mu, t_UR);
			Solve(M);
		}
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
	Scalar DT = Scalar(1.0e-4);
	Scalar t_UR = Scalar(1);

	/*transport*/
	Util::ParamList params("transport");
	params.enroll("t_UR", &t_UR);
	params.enroll("DT", &DT);

	/*read parameters*/
	Util::read_params(input,MP::printOn);
	
	/*AMR iteration*/
	for (AmrIteration ait; !ait.end(); ait.next()) {
		
		VectorCellField U("U", READWRITE);
		ScalarCellField T("T", READWRITE);
		
		/*Time loop*/
		Iteration it(ait.get_step());
		ScalarFacetField F = flx(U);
		ScalarCellField mu = DT;
		for (; !it.end(); it.next()) {
			ScalarCellMatrix M;
			M = transport(T, U, F, mu, t_UR);
			Solve(M);
		}
	}
}
/**
 \verbatim
 Wave equation solver
 ~~~~~~~~~~~~~~~~~~~~
 Solver for pdes of hyperbolic wave equation type:
       d2T/dt2 = c^2 * lap(T)
 \endverbatim
*/
void wave(istream& input) {
	/*Solver specific parameters*/
	Scalar C2 = Scalar(1);
	Scalar t_UR = Scalar(1);

	/*diffusion*/
	Util::ParamList params("wave");
	params.enroll("C2", &C2);
	params.enroll("t_UR", &t_UR);

	/*read parameters*/
	Util::read_params(input,MP::printOn);
	
	/*AMR iteration*/
	for (AmrIteration ait; !ait.end(); ait.next()) {
		
		ScalarCellField T("T", READWRITE);
		
		/*Time loop*/
		ScalarCellField mu = C2;
		for (Iteration it(ait.get_step()); !it.end(); it.next()) {
			ScalarCellMatrix M = -lap(T,mu);
			M.cF = &T;
			addTemporal<2>(M,t_UR);
			Solve(M);
		}
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

	/*read parameters*/
	Util::read_params(input,MP::printOn);

	/*AMR iteration*/
	for (AmrIteration ait; !ait.end(); ait.next()) {
		
		VectorCellField U("U", READWRITE);
		ScalarCellField p("p", READ);
	
		/*set internal field to zero*/
		for (Int i = 0; i < Mesh::gBCSfield; i++) {
			U[i] = Vector(0, 0, 0);
			p[i] = Scalar(0);
		}
		applyExplicitBCs(U, true);
		applyExplicitBCs(p, true);

		/*Time loop*/
		for (Iteration it(ait.get_step()); it.start(); it.next()) {
			/*solve potential equation*/
			ScalarCellField divU = divf(U);
			ScalarCellField one = Scalar(1);
			for (Int k = 0; k <= n_ORTHO; k++)
				Solve(lap(p, one, true) == divU);

			/*correct velocity*/
			U -= gradi(p);
			applyExplicitBCs(U, true);
		}
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
	Util::read_params(input,MP::printOn);

	/*solve*/
	for (AmrIteration ait; !ait.end(); ait.next()) {
		Mesh::calc_walldist(ait.get_step(), n_ORTHO);
	}
}
void Mesh::calc_walldist(Int step, Int n_ORTHO) {
	ScalarCellField& phi = yWall;
	/*poisson equation*/
	{
		ScalarCellField one = Scalar(1);
		for (Int k = 0; k <= n_ORTHO; k++)
			Solve(lap(phi, one, true) == -cV);
	}
	/*wall distance*/
	{
		VectorCellField g = gradi(phi);
		yWall = sqrt((g & g) + 2 * phi) - mag(g);
	}
	/*write it*/
	yWall.write(step);
}
