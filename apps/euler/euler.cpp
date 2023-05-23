#include "solve.h"
#include "iteration.h"
#include "wrapper.h"
#include "properties.h"

/**
  \verbatim
  Euler equations solver
  ~~~~~~~~~~~~~~~~~~~~~~
  d(rho*(1))/dt + div(1,F,0) = 0
  d(rho*(U))/dt + div(U,F,0) = -grad(p) + rho * gravity
  d(rho*(T))/dt + div(T,F,0) = 0
  \endverbatim
 */
void euler(std::istream& input) {
    /*Solver specific parameters*/
    Scalar pressure_UR = Scalar(0.5);
    Scalar velocity_UR = Scalar(0.8);
    Scalar t_UR = Scalar(0.8);
    Int buoyancy = 1;
    Int diffusion = 1;

    /*fluid properties*/
    {
        Util::ParamList params("general");
        Fluid::enroll(params);
        params.read(input);
    }

    /*euler options*/
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
        using namespace Fluid;
        Scalar p_gamma, R, psi, iPr;
        R = cp - cv;
        p_gamma = cp / cv;
        psi = 1 / (R * T0);
        iPr = 1 / Pr;

        /*calculate rho*/
        rho = (P0 / (R*T)) * pow(p / P0, 1 / p_gamma);
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
                ScalarCellMatrix M;
                M = transport(T, Fc, F, mu * iPr, t_UR, &rho);
                Solve(M);
            }
            /*calculate p*/
            p = P0 * pow((rho*T*R) / P0, p_gamma);
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
  Main application entry point for euler solver.
  \endverbatim
 */
int main(int argc, char* argv[]) {
   MP mp(argc, argv);
   Solver::Initialize(argc, argv);
   euler(Solver::input);
   Solver::Finalize();
   return 0;
}
