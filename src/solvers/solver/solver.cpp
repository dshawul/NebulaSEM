#include "field.h"
#include "turbulence.h"
#include "mp.h"
#include "system.h"
#include "solve.h"

using namespace std;

/*general properties*/
namespace GENERAL {
	Scalar density = 1;
	Scalar viscosity = 1e-5;
	Scalar conductivity = 1e-4;
	Vector gravity = Vector(0,0,-9.81);

	void enroll(Util::ParamList& params) {
		params.enroll("rho",&density);
		params.enroll("viscosity",&viscosity);
		params.enroll("conductivity",&conductivity);
		params.enroll("gravity",&gravity);
	}
};

/*solvers*/
void piso(istream&);
void diffusion(istream&);
void potential(istream&);
void transport(istream&);
void walldist(istream&);

/**
\verbatim
 Main application entry point for different solvers.
 \endverbatim
 */
int main(int argc,char* argv[]) {

	/*message passing object*/
	MP mp(argc,argv);
	ifstream input(argv[1]);

	/*main options*/
	Util::ParamList params("general");
	string sname;
	params.enroll("solver",&sname);
	params.enroll("mesh",&Mesh::gMeshName);
	Mesh::enroll(params);
	GENERAL::enroll(params);
	params.read(input);

	/*Mesh*/
	if(mp.n_hosts > 1) {
		stringstream s;
		s << Mesh::gMeshName << mp.host_id;
		if(!System::cd(s.str()))
			return 1;
	}
	Mesh::readMesh();
	Mesh::initGeomMeshFields();
	atexit(Util::cleanup);

	/*call solver*/
	if(!Util::compare(sname,"piso")) {
		piso(input);
	} else if(!Util::compare(sname,"diffusion")) {
		diffusion(input);
	} else if(!Util::compare(sname,"transport")) {
		transport(input);
	} else if(!Util::compare(sname,"potential")) {
		potential(input);
	} else if(!Util::compare(sname,"walldist")) {
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
		forEachField(initTimeSeries());
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
		forEachField(updateTimeSeries(i));

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
	  One jacobi sweep is done to find Ua.
    Step 2)
      Solve poisson pressure equation to satisfy continuity with fluxes calculated 
	  from interpolated Ua.
	      div(Up) = 0
		  div(1/ap * grad(p)) = div(H(U)/ap)
		  lap(p,1/ap) = div(Ua)
    Step 3)
	  Correct the velocity with gradient of newly found pressure
	      U -= grad(p)
    These steps are repeated two or more times for transient solutions.
	For steady state problems once is enough.
	\endverbatim
*/
void piso(istream& input) {
	/*Solver specific parameters*/
	Scalar& rho = GENERAL::density;
	Scalar& viscosity = GENERAL::viscosity;
	Scalar velocity_UR = Scalar(0.8);
	Scalar pressure_UR = Scalar(0.5);
	Int n_PISO = 1;
	Int n_ORTHO = 0;

	/*piso options*/
	Util::ParamList params("piso");
	params.enroll("velocity_UR",&velocity_UR);
	params.enroll("pressure_UR",&pressure_UR);
	params.enroll("n_PISO",&n_PISO);
	params.enroll("n_ORTHO",&n_ORTHO);

	VectorCellField U("U",READWRITE);
	ScalarCellField p("p",READWRITE);

	/*turbulence model*/
	ScalarFacetField F;
	bool Steady;
	Turbulence_Model::RegisterTable(params);
	params.read(input);
	Turbulence_Model* turb = 
		Turbulence_Model::Select(U,F,rho,viscosity,Steady);
	turb->enroll();

	/*read parameters*/
	Util::read_params(input);

	/*wall distance*/
	if(turb->needWallDist())
		Mesh::calc_walldist(Iteration::get_start());

	/*time*/
	Scalar time_factor = Controls::time_scheme_factor;
	Steady = (Controls::state == Controls::STEADY);

	/*Calculate for each time step*/
	Iteration it;
	ScalarCellField po = p;
	VectorCellField gP = -gradV(p);
	F = flx(rho * U); 

	for(;!it.end();it.next()) {
		/*Form Navier-stokes equation*/
		VectorMeshMatrix M;

		/*convection*/
		{
			ScalarFacetField mu = rho * viscosity;
			M = div(U,F,mu);
		}

		/*viscous/turbulent stress*/
		turb->addTurbulentStress(M);

		/*relax if steady state otherwise add time contribution*/
		if(Steady)
			M.Relax(velocity_UR);
		else {
			/*crank nicolson*/
			if(!equal(time_factor,1)) {
				VectorCellField po = M * U;
				M *= time_factor;
				M.Su -= (1 - time_factor) * po;
			}
			/*time derivative*/
			M += ddt(U,rho);
		}

		/*solve momentum equation*/
		Solve(M == gP);

		/*1/ap*/
		ScalarCellField api = (1 / M.ap);
		fillBCs(api,true);
		ScalarCellField rmu = rho * api * Mesh::cV;

		/*PISO loop*/
		for(Int j = 0;j < n_PISO;j++) {
			/* Ua = H(U) / ap*/
			U = getRHS(M) * api;
			updateExplicitBCs(U,true);

			/*solve pressure poisson equation to satisfy continuity*/
			{
				ScalarCellField rhs = div(rho * U);
				for(Int k = 0;k <= n_ORTHO;k++)
					Solve(lap(p,rmu) += rhs);
			}

			/*explicit velocity correction : add pressure contribution*/
			gP = -gradV(p);
			U -= gP * api;
			updateExplicitBCs(U,true);
		}

		/*update fluctuations*/
		updateExplicitBCs(U,true,true);
		F = flx(rho * U);

		/*solve turbulence transport equations*/
		turb->solve();

		/*explicitly under relax pressure*/
		if(Steady) {
			p.Relax(po,pressure_UR);
			gP = -gradV(p);
			po = p;
		}
	}

	/*write calculated turbulence fields*/
	if(turb->writeStress) {
		ScalarCellField K("Ksgs",WRITE);
		STensorCellField R("Rsgs",WRITE);
		STensorCellField V("Vsgs",WRITE);
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
       d(rho*u)/dt = lap(T,rho*DT)
 \endverbatim
*/
void diffusion(istream& input) {
	/*Solver specific parameters*/
	Scalar& rho = GENERAL::density;
	Scalar DT = Scalar(1);
	Scalar t_UR = Scalar(1);

	/*diffusion*/
	Util::ParamList params("diffusion");
	params.enroll("DT",&DT);
	params.enroll("t_UR",&t_UR);

	ScalarCellField T("T",READWRITE);

	/*read parameters*/
	Util::read_params(input);

	/*time*/
	Scalar time_factor = Controls::time_scheme_factor;
	bool Steady = (Controls::state == Controls::STEADY);

	/*Calculate for each time step*/
	ScalarFacetField mu = rho * DT;

	for(Iteration it;!it.end();it.next()) {
		ScalarMeshMatrix M;

		M = -lap(T,mu);

		if(Steady)
			M.Relax(t_UR);
		else {
			if(!equal(time_factor,1)) {
				ScalarCellField po = M * T;
				M *= time_factor;
				M.Su -= (1 - time_factor) * po;
			}
			M += ddt(T,rho);
		}

		Solve(M);
	}
}
/**
  \verbatim
  Transport equation solver
  ~~~~~~~~~~~~~~~~~~~~~~~~~
  Given a flow field (U) and values of a scalar at the boundaries, 
  the solver determines the distribution of the scalar.
     dT/dt + div(T,F,mu) = lap(T,mu)
  \endverbatim
 */
void transport(istream& input) {
	/*Solver specific parameters*/
	Scalar& rho = GENERAL::density;
	Scalar DT = Scalar(4e-2);
	Scalar t_UR = Scalar(1);

	/*transport*/
	Util::ParamList params("transport");
	params.enroll("DT",&DT);
	params.enroll("t_UR",&t_UR);

	VectorCellField U("U",READWRITE);
	ScalarCellField T("T",READWRITE);

	/*read parameters*/
	Util::read_params(input);

	/*time*/
	Scalar time_factor = Controls::time_scheme_factor;
	bool Steady = (Controls::state == Controls::STEADY);

	/*Calculate for each time step*/
	ScalarFacetField F,mu = rho * DT,gamma;

	for(Iteration it;!it.end();it.next()) {
		ScalarMeshMatrix M;

		F = flx(rho * U); 
		M = div(T,F,mu) 
			- lap(T,mu);

		if(Steady)
			M.Relax(t_UR);
		else {
			if(!equal(time_factor,1)) {
				ScalarCellField po = M * T;
				M *= time_factor;
				M.Su -= (1 - time_factor) * po;
			}
			M += ddt(T,rho);
		}
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
	params.enroll("n_ORTHO",&n_ORTHO);

	VectorCellField U("U",READWRITE);
	ScalarCellField p("p",READ);
	
	/*read parameters*/
	Util::read_params(input);

	/*set internal field to zero*/
	for(Int i = 0;i < Mesh::gBCellsStart;i++) {
		U[i] = Vector(0,0,0);
		p[i] = Scalar(0);
	}
	updateExplicitBCs(U,true);
	updateExplicitBCs(p,true);

	for(Iteration it;it.start();it.next()) {
		/*solve potential equation*/
		ScalarCellField divU = div(U);
		ScalarFacetField one = Scalar(1);
		for(Int k = 0;k <= n_ORTHO;k++)
			Solve(lap(p,one) == divU);

		/*correct velocity*/
		U -= grad(p);
		updateExplicitBCs(U,true);
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
	params.enroll("n_ORTHO",&n_ORTHO);
	Util::read_params(input);
	
	/*solve*/
	Mesh::calc_walldist(Iteration::get_start(),n_ORTHO);
}
void Mesh::calc_walldist(Int step,Int n_ORTHO) {
	ScalarCellField& phi = yWall;
    /*poisson equation*/
	{
		ScalarFacetField one = Scalar(1);
		for(Int k = 0;k <= n_ORTHO;k++)
			Solve(lap(phi,one) == -cV);
	}
	/*wall distance*/
	{
		VectorCellField g = grad(phi);
		yWall = sqrt((g & g) + 2 * phi) - mag(g);
	}
	/*write it*/
	yWall.write(step);
}
