#include "field.h"
#include "turbulence.h"
#include "ke.h"
#include "rngke.h"
#include "realizableke.h"
#include "kw.h"
#include "les.h"
#include "mp.h"
#include "system.h"
#include "solve.h"
#include "vtk.h"

using namespace std;

/*solvers*/
void piso(istream&);
void diffusion(istream&);
void potential(istream&);
void transport(istream&);

/********************
 * Main application
 *******************/
int main(int argc,char* argv[]) {

	/*message passing object*/
	MP mp(argc,argv);
	ifstream input(argv[1]);

	string sname;
	Util::StringParams::enroll("solver",&sname);
	Util::read_param(input,"solver");

	/*Mesh*/
	Mesh::enroll();
	Util::read_param(input,"mesh");
	if(mp.n_hosts > 1) {
		stringstream s;
		s << Mesh::gMeshName << mp.host_id;
		Mesh::gMeshName = s.str();
		if(!System::cd(Mesh::gMeshName))
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
	}

	return 0;
}

/********************************************
 * Navier stokes solver using PISO algorithm
 ********************************************/
void piso(istream& input) {
	/*Solver specific parameters*/
	Scalar rho = 1;
	Scalar viscosity = 1e-5;
	Scalar velocity_UR = Scalar(0.8);
	Scalar pressure_UR = Scalar(0.5);
	Int n_PISO = 1;
	Int n_ORTHO = 0;
	Int LESaverage = 0;

	Util::ScalarParams::enroll("rho",&rho);
	Util::ScalarParams::enroll("viscosity",&viscosity);
	Util::ScalarParams::enroll("velocity_UR",&velocity_UR);
	Util::ScalarParams::enroll("pressure_UR",&pressure_UR);
	Util::IntParams::enroll("n_PISO",&n_PISO);
	Util::IntParams::enroll("n_ORTHO",&n_ORTHO);

	VectorCellField U("U",READWRITE);
	ScalarCellField p("p",READWRITE);

	/*turbulence model*/
	enum TurbModel {NONE, KE, RNG_KE,REALIZABLE_KE, KW, LES};
	TurbModel turb_model = KE;
	Util::Option* op;
	op = new Util::Option(&turb_model,6,"NONE","KE","RNG_KE","REALIZABLE_KE","KW","LES");
	Util::OptionParams::enroll("turbulence_model",op);
	op = new Util::Option(&LESaverage,2,"NO","YES");
	Util::OptionParams::enroll("les_average",op);

	/*Select turbulence model*/
	ScalarFacetField F;
	bool Steady;

	Turbulence_Model* turb;
	Util::read_param(input,"turbulence_model");
	switch(turb_model) {
		case KE:   
			turb = new KE_Model(U,F,rho,viscosity,Steady); 
			break;
		case RNG_KE:   
			turb = new RNG_KE_Model(U,F,rho,viscosity,Steady); 
			break;
		case REALIZABLE_KE:   
			turb = new REALIZABLE_KE_Model(U,F,rho,viscosity,Steady); 
			break;
		case KW:   
			turb = new KW_Model(U,F,rho,viscosity,Steady); 
			break;
		case LES:  
			turb = new LES_Model(U,F,rho,viscosity,Steady); 
			break;
		default:
			turb = new Turbulence_Model(U,F,rho,viscosity,Steady); 
			break;
	}
	turb->enroll();

	/*read parameters*/
	std::cout << "**************************\n";
	Util::read_params(input);
	std::cout << "==========================\n";
	
	/*average statistics for LES */
	VectorCellField Uavg(0),Ustd(0);
	ScalarCellField pavg(0),pstd(0);
	if(LESaverage) {
		Uavg.construct("Uavg",READWRITE);
		Ustd.construct("Ustd",READWRITE);
		pavg.construct("pavg",READWRITE);
		pstd.construct("pstd",READWRITE);
	}

	/*instantaneous values*/
	IntVector probe_points;
	for(Int j = 0;j < Mesh::probePoints.size();j++) {
		Vector v = Mesh::probePoints[j];
		Int index = Mesh::findNearest(v);
		probe_points.push_back(index);
	}
	Int probe = probe_points.size();
	ofstream oUi,opi;
	if(probe) {
		oUi.open("Ui");
		opi.open("pi");
	}

	/*Read at selected start time step*/
	Int step,start;
	step = Controls::start_step / Controls::write_interval;
	start = Controls::write_interval * step + 1;
	Mesh::read_fields(step);
	Util::write_vtk(step);

	/*time*/
	Scalar time_factor = Controls::time_scheme_factor;
	Steady = (Controls::state == Controls::STEADY);

	/*Calculate for each time step*/
	VectorCellField gP = -src(grad(p));
	F = div(rho * U); 

	for(Int i = start; i <= Controls::end_step; i++) {
		/*Print step*/
		if(MP::host_id == 0) {
			if(Steady)
				MP::printH("Step %d\n",i);
			else
				MP::printH("Time %f\n",i * Controls::dt);
		}
		/*Momentum and pressure solution*/
		{
			VectorMeshMatrix M;
			{
				/*convection*/
				{
					ScalarFacetField mu = rho * viscosity;
					M = div(U,F,mu);
				}
				/*turbulent stress*/
				turb->addTurbulentStress(M);
				/*end*/
			}
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

			/*PISO loop*/
			for(Int j = 0;j < n_PISO;j++) {
				/* Separate pressure gradient from the rest U = H(U) / ap - grad(p) / ap
				 * Velocity without the effect of pressure is therefore U = H(U) / ap
				 * Later we will correct it by adding the second term
				 * H(U) is calculated with previous estimate of U by doing a jacobi sweep
				 */
				U = getRHS(M) * api;
				updateExplicitBCs(U,true);
				/*solve pressure poisson equation to satisfy continuity*/
				{
					ScalarCellField po;
					if(Steady)
						po = p;
					F = div(rho * U);
					for(Int k = 0;k <= n_ORTHO;k++)
						Solve((lap(p,rho * api * Mesh::cV) += sum(F)));
					if(Steady)
						p.Relax(po,pressure_UR);
				}
				gP = -src(grad(p));
				/*explicit velocity correction : add pressure contribution*/
				U -= gP * api;
				updateExplicitBCs(U,true);
				/*end*/
			}
		}
		/*update fluctuations*/
		updateExplicitBCs(U,true,true);
		F = div(rho * U);
		
		/*solve transport equations*/
		turb->solve();

		/*average*/
		if(LESaverage) {
			Uavg += U;
			pavg += p;
			Ustd += (U * U);
			pstd += (p * p);
		}

		/*store instantaneous values for some locations*/
		if(probe) {
			oUi << i << " ";
			opi << i << " ";
			for(Int j = 0;j < probe_points.size();j++) {
				oUi << U[probe_points[j]] << " ";
				opi << p[probe_points[j]] << " ";
			}
			oUi << endl;
			opi << endl;
		}
		/*write result to file*/
		if((i % Controls::write_interval) == 0) {
			step = i / Controls::write_interval;
	
			if(LESaverage) {
				VectorCellField Ua = Uavg,Us = Ustd;
				ScalarCellField pa = pavg,ps = pstd;
				Scalar n = Scalar(i);
				Uavg /= n;
				pavg /= n;
				Ustd += Uavg * (n * Uavg - 2 * Ua);
                pstd += pavg * (n * pavg - 2 * pa);
				Ustd = sqrt(Ustd / n);
				pstd = sqrt(pstd / n);

				Mesh::write_fields(step);
				Util::write_vtk(step);

				Uavg = Ua;
				Ustd = Us;
				pavg = pa;
				pstd = ps;
			} else {
				Mesh::write_fields(step);
				Util::write_vtk(step);
			}
		}
		/*end*/
	}
}
/**************************
 * potential flow solver
 **************************/
void potential(istream& input) {
	/*Solver specific parameters*/
	Int n_ORTHO = 0;
	Util::IntParams::enroll("n_ORTHO",&n_ORTHO);

	VectorCellField U("U",READWRITE);
	ScalarCellField p("p",READ);
	
	/*read parameters*/
	std::cout << "**************************\n";
	Util::read_params(input);
	std::cout << "==========================\n";

	/*Read at selected start time step*/
	Int step,start;
	step = Controls::start_step / Controls::write_interval;
	start = Controls::write_interval * step + 1;
	Mesh::read_fields(step);

	/*set internal field to zero*/
	for(Int i = 0;i < Mesh::gBCellsStart;i++) {
		U[i] = Vector(0,0,0);
		p[i] = Scalar(0);
	}
	updateExplicitBCs(U,true);
	updateExplicitBCs(p,true);

	/*solve potential equation*/
	ScalarFacetField F = div(U),one = Scalar(1);
	for(Int k = 0;k <= n_ORTHO;k++)
		Solve(lap(p,one) == sum(F));

	/*correct velocity*/
	U -= grad(p);
	updateExplicitBCs(U,true);

	/*write result to file*/
	Mesh::write_fields(step);
	Util::write_vtk(step);
}
/********************************************
 * laplace diffusion equation solver
 ********************************************/
void diffusion(istream& input) {
	/*Solver specific parameters*/
	Scalar DT = Scalar(1);
	Scalar t_UR = Scalar(1);

	Util::ScalarParams::enroll("DT",&DT);
	Util::ScalarParams::enroll("t_UR",&t_UR);

	ScalarCellField T("T",READWRITE);

	/*read parameters*/
	std::cout << "**************************\n";
	Util::read_params(input);
	std::cout << "==========================\n";

	/*Read at selected start time step*/
	Int step,start;
	step = Controls::start_step / Controls::write_interval;
	start = Controls::write_interval * step + 1;
	Mesh::read_fields(step);
	Util::write_vtk(step);

	/*time*/
	Scalar time_factor = Controls::time_scheme_factor;
	bool Steady = (Controls::state == Controls::STEADY);

	/*Calculate for each time step*/
	ScalarFacetField mu = DT;
	ScalarCellField rho = Scalar(1);

	for(Int i = start; i <= Controls::end_step; i++) {
		/*Print step*/
		if(MP::host_id == 0)
			MP::printH("Time %f\n",i * Controls::dt);

		/*solve*/
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
		
		/*write result to file*/
		if((i % Controls::write_interval) == 0) {
			step = i / Controls::write_interval;
			Mesh::write_fields(step);
			Util::write_vtk(step);
		}
	}
}
/***********************************************
 * transport of a scalar with given flow field
 **********************************************/
void transport(istream& input) {
	/*Solver specific parameters*/
	Scalar DT = Scalar(4e-2);
	Scalar t_UR = Scalar(1);

	Util::ScalarParams::enroll("DT",&DT);
	Util::ScalarParams::enroll("t_UR",&t_UR);

	VectorCellField U("U",READWRITE);
	ScalarCellField T("T",READWRITE);

	/*read parameters*/
	std::cout << "**************************\n";
	Util::read_params(input);
	std::cout << "==========================\n";

	/*Read at selected start time step*/
	Int step,start;
	step = Controls::start_step / Controls::write_interval;
	start = Controls::write_interval * step + 1; 
	Mesh::read_fields(step);
	Util::write_vtk(step);

	/*time*/
	Scalar time_factor = Controls::time_scheme_factor;
	bool Steady = (Controls::state == Controls::STEADY);

	/*Calculate for each time step*/
	ScalarFacetField F,mu = DT,gamma;
	ScalarCellField rho = Scalar(1);

	for(Int i = start; i <= Controls::end_step; i++) {
		/*Print step*/
		if(MP::host_id == 0)
			MP::printH("Time %f\n",i * Controls::dt);

		/*solve*/
		ScalarMeshMatrix M;

		F = div(rho * U); 
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
		
		/*write result to file*/
		if((i % Controls::write_interval) == 0) {
			step = i / Controls::write_interval;
			Mesh::write_fields(step);
			Util::write_vtk(step);
		}
	}
}



