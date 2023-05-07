#include "field.h"
#include "turbulence.h"
#include "mp.h"
#include "system.h"
#include "solve.h"
#include "iteration.h"

using namespace std;

/**
  general properties
 */
namespace General {
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
void hydro_balance(istream&);
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
            << "  ./solvers <inputfile>\n"
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
        General::enroll(params);
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
    } else if (!Util::compare(sname, "hydro_balance")) {
        hydro_balance(input);
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
    Scalar velocity_UR = Scalar(0.8);
    Scalar pressure_UR = Scalar(0.5);
    Scalar t_UR = Scalar(0.8);
    Int n_PISO = 1;
    Int n_ORTHO = 0;
    Int momentum_predictor = 0;
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
    op = new Util::BoolOption(&momentum_predictor);
    params.enroll("momentum_predictor",op);

    Turbulence_Model::RegisterTable(params);
    params.read(input);

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {
        ScalarCellField rho = General::density;
        ScalarCellField mu = rho * General::viscosity;
        VectorCellField U("U", READWRITE);
        ScalarCellField p("p", READWRITE);
        ScalarCellField T(false);

        /*temperature*/
        if(buoyancy != NONE)
            T.construct("T",READWRITE);

        /*turbulence model*/
        VectorCellField Fc;
        ScalarFacetField F;
        Turbulence_Model* turb = Turbulence_Model::Select(U, Fc, F, rho, mu);

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

        Fc = flxc(rho * U);
        F = flx(rho * U);

        for (; !it.end(); it.next()) {
            /*
             * Prediction
             */
            VectorCellMatrix M;
            {
                VectorCellField Sc = turb->getExplicitStresses();
                const ScalarCellField eddy_mu = turb->getTurbVisc();
                /* Add buoyancy in two ways
                 *  1. pm = p - rho*g*h
                 *  2. pm = p - rho_m*g*h
                 */
                if (buoyancy != NONE) {
                    Scalar beta;
                    if(buoyancy <= BOUSSINESQ_T2) 
                        beta = General::beta;
                    else
                        beta = 1 / General::T0;
                    if (buoyancy == BOUSSINESQ_T1 || buoyancy == BOUSSINESQ_THETA1) {  
                        ScalarCellField rhok = rho * (0 - beta * (T - General::T0));
                        Sc += (rhok * VectorCellField(Controls::gravity));
                    } else if(buoyancy == BOUSSINESQ_T2 || buoyancy == BOUSSINESQ_THETA2) {
                        ScalarCellField gz = dot(Mesh::cC,VectorCellField(Controls::gravity));
                        Sc += gz * (rho * beta) * gradi(T);
                    }
                }
                /*momentum prediction*/
                {
                    const ScalarCellField eff_mu = eddy_mu + mu;
                    M = transport(U, Fc, F, eff_mu, velocity_UR, Sc, Scalar(0));
                    if(momentum_predictor) {
                        Solve(M == gP);
                    }
                }
            }

            /*
             * Correction
             */
            const ScalarCellField api = fillBCs<Scalar>(1.0 / M.ap);
            const ScalarCellField rmu = rho * api * Mesh::cV;

            /*PISO loop*/
            for (Int j = 0; j < n_PISO; j++) {
                /* Ua = H(U) / ap*/
                U = getRHS(M) * api;
                applyExplicitBCs(U, true);

                /*solve pressure poisson equation to satisfy continuity*/
                {
                    const ScalarCellField rhs = divf(rho * U);
                    for (Int k = 0; k <= n_ORTHO; k++)
                        Solve(lap(p, rmu, true) += rhs);
                }

                /*explicit velocity correction : add pressure contribution*/
                gP = -gradf(p);
                U -= gP * api;
                applyExplicitBCs(U, true);
            }

            /*update fluctuations*/
            applyExplicitBCs(U, true, true);
            Fc = flxc(rho * U);
            F = flx(rho * U);

            /*solve turbulence transport equations*/
            turb->solve();

            /*solve energy transport*/
            if (buoyancy != NONE) {
                const ScalarCellField eddy_mu = turb->getTurbVisc();
                const ScalarCellField eff_mu = eddy_mu / General::Prt + mu / General::Pr;
                ScalarCellMatrix Mt = transport(T, Fc, F, eff_mu, t_UR);
                Solve(Mt);
                T = max(T, Constants::MachineEpsilon);
            }

            /*explicitly under relax pressure*/
            if (Controls::state == Controls::STEADY) {
                p = po + (p - po) * pressure_UR;
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
        VectorCellField Fc;
        ScalarFacetField F;
        ScalarCellField mu;

        /*gas constants*/
        Scalar p_factor, p_gamma, R, psi, iPr;
        using namespace General;
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
            Fc = flxc(rho * U);
            F = flx(rho * U);
            /*rho-equation*/
            {
                ScalarCellMatrix M;
                M = convection(rho, flxc(U), flx(U), pressure_UR);
                Solve(M);
            }
            /*artificial viscosity*/
            if(diffusion) mu = rho * viscosity;
            else mu = Scalar(0);
            /*T-equation*/
            {
                ScalarCellMatrix M,Mt;
                M = transport(T, Fc, F, mu * iPr, t_UR, &rho);
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
                M = transport(U, Fc, F, mu, velocity_UR, Sc, Scalar(0), &rho);
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
        VectorCellField Fc = flxc(U);
        ScalarFacetField F = flx(U);
        for (; !it.end(); it.next()) {
            ScalarCellMatrix M;
            M = convection(T, Fc, F, t_UR);
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
        VectorCellField Fc = flxc(U);
        ScalarFacetField F = flx(U);
        ScalarCellField mu = DT;
        for (; !it.end(); it.next()) {
            ScalarCellMatrix M;
            M = transport(T, Fc, F, mu, t_UR);
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
  but the theory can still be used to initialize flow field for complex simulations.

  The potential flow assumption is that velocity is the gradient of a scalar field,
  which is the velocity potential (phi) 
         U = grad(phi)
         curl(U) = curl(grad(phi))
         curl(U) = 0
  where the last step is possible due to the vector identity curl(grad(phi)) = 0.
  Hence, defining the velocity as a gradient of a scalar ensures vorticity is zero.
  Let us now take the divergence instead as
         U = grad(phi)
         div(U) = div(grad(phi))
         div(U) = lap(phi)
  Since div(U)=0 for incompressible, the poisson equation becomes a laplace equation
         lap(phi) = 0

  What the solver does
  ~~~~~~~~~~~~~~~~~~~~~
  Given an initial velocity field (Ua) that does not satisfy continuity (div(Ua) != 0),
  we can correct Ua with the pressure gradient to get a divergence free velocity field U
         U = Ua - grad(p)
         div(Ua - grad(p)) = div(U) = 0
         div(Ua) = div(grad(p))
  If Ua is irrotational, i.e. curl(Ua) = 0 and Ua = grad(phia), then so is U because
         U = Ua - grad(p)
         U = grad(phia) - grad(p)
         U = grad(phia - p)
  where U = grad(phi) such that phi = phia - p.
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
        ScalarCellField p("p", READWRITE);

        /*Time loop*/
        for (Iteration it(ait.get_step()); it.start(); it.next()) {
            const ScalarCellField divU = divf(U);
            const ScalarCellField one = Scalar(1);

            /*solve pressure poisson equation for correction*/
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
  Hydrostatic
  ~~~~~~~~~~~~~
     Solver for hydrostatic balance
         grad(p) = -rho*g
     Using gravitational potential theory
        div(grad(p)) = div(-rho*g)
  \endverbatim
*/
void hydro_balance(istream& input) {
    /*Solver specific parameters*/
    Int n_ORTHO = 0;

    /*potential*/
    Util::ParamList params("hydro_balance");
    params.enroll("n_ORTHO", &n_ORTHO);

    /*read parameters*/
    Util::read_params(input,MP::printOn);

    /*solve*/

    for (AmrIteration ait; !ait.end(); ait.next()) {

        ScalarCellField p("p", READWRITE);

        for (Iteration it(ait.get_step()); it.start(); it.next()) {

            const ScalarCellField one = Scalar(1);
            const VectorCellField rhog = General::density * Controls::gravity;
            const ScalarCellField ndivRhoG = -divf(rhog);

            /*solve poisson equation*/
            for (Int k = 0; k <= n_ORTHO; k++)
                Solve(lap(p, one, true) == ndivRhoG);
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
/**
  Calculate wall distance at given time step
 */
void Mesh::calc_walldist(Int step, Int n_ORTHO) {
    ScalarCellField& phi = yWall;
    /*poisson equation*/
    {
        const ScalarCellField one = Scalar(1);
        for (Int k = 0; k <= n_ORTHO; k++)
            Solve(lap(phi, one, true) == -cV);
    }
    /*wall distance*/
    {
        const VectorCellField g = gradi(phi);
        yWall = sqrt((g & g) + 2 * phi) - mag(g);
    }
    /*write it*/
    yWall.write(step);
}
