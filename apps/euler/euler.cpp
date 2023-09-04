#include "solve.h"
#include "iteration.h"
#include "wrapper.h"
#include "properties.h"

//#define LINEARIZE

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

    /*total mass and energy of system*/
    Scalar mass0, energy0;

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {

        ScalarCellField p("p",READWRITE);
        VectorCellField U("U",READWRITE);
        ScalarCellField T("T",READWRITE);
        ScalarCellField rho("rho",WRITE);
        VectorCellField g(false);

        /*Read fields*/
        Iteration it(ait.get_step());
        VectorCellField Fc;
        ScalarFacetField F;
        ScalarCellField rhof;
        ScalarCellField mu;
        ScalarFacetField lambdaMax;

        /*gas constants*/
        using namespace Fluid;
        Scalar p_gamma, R, psi, iPr;
        R = cp - cv;
        p_gamma = cp / cv;
        psi = 1 / (R * T0);
        iPr = 1 / Pr;

        /*buoyancy*/
        ScalarCellField p_ref, rho_ref, gh;
        if(buoyancy) {
            /*gravity vector*/
            g.construct("gravity",WRITE);
            if(Mesh::is_spherical) {
                g = -unit(Mesh::cC) * mag(Controls::gravity);
                gh = -(mag(Mesh::cC) - Mesh::sphere_radius) * mag(Controls::gravity);
            } else {
                g = Controls::gravity;
                gh = dot(g,Mesh::cC);
            }
            Mesh::fixedBCs<Vector>(U,g);
            applyExplicitBCs(g,true);
            if(ait.get_step() == 0)
                g.write(0);

            /*compute hydrostatic pressure*/
            p_ref = P0 * pow(1.0 + gh / (cp * T0), cp / R);
            p += P0 * pow(1.0 + gh / (cp * (T + T0)), cp / R);
        } else {
            p_ref = P0;
            p += P0;
        }

        Mesh::scaleBCs<Scalar>(p,p_ref,1.0);
        applyExplicitBCs(p_ref,true);
        applyExplicitBCs(p,true);

        /*calculate rho*/
        rho_ref = (P0 / (R*T0)) * pow(p_ref / P0, 1 / p_gamma);
        Mesh::scaleBCs<Scalar>(p,rho_ref,psi);
        applyExplicitBCs(rho_ref,true);

        rho = (P0 / (R*(T + T0))) * pow(p / P0, 1 / p_gamma);
        Mesh::scaleBCs<Scalar>(p,rho,psi);
        applyExplicitBCs(rho,true);

        /*compute total mass and energy*/
        Scalar mass0, energy0;
        {
            ScalarCellField sf = rho * Mesh::cV;
            mass0 = reduce_sum(sf);
            //potential + kinetic + internal
            sf = gh +
                 0.5 * magSq(U) +
                 pow(p / P0, R / cp) * (T + T0) * cv;
            sf = rho * Mesh::cV * sf;
            energy0 = reduce_sum(sf);
        }

        /*write calculated rho*/
        rho -= rho_ref;
        if(ait.get_step() == 0) {
            rho.write(0);
        }

        /*Time loop*/
        for (; !it.end(); it.next()) {

            /*add reference values*/
            rho += rho_ref;
            T += T0;

            /*fluxes*/
            Fc = flxc(rho * U);
            F = flx(rho * U);
#ifndef LINEARIZE
            lambdaMax = ( cds(mag(U)) + cds(sqrt(p_gamma*R*T)) ) / 2;
#endif

            /*artificial viscosity*/
            if(diffusion) mu = rho * viscosity;
            else mu = Scalar(0);

            /*rho-equation*/
            rhof = rho;
            {
                ScalarCellMatrix M;
#ifdef LINEARIZE
                M = convection(rho, flxc(U), flx(U), pressure_UR);
#else
                {
                    VectorCellField fq = U * rho;
                    ScalarFacetField Fr = flx(U);
                    M = divf(fq,false,&Fr,&rho,&lambdaMax);
                    M.cF = &rho;
                    addTemporal<1>(M,pressure_UR);
                }
#endif
                Solve(M);
            }

            /*calculate p*/
            p = P0 * pow((rho*T*R) / P0, p_gamma);
            applyExplicitBCs(p,true);
            p -= p_ref;

            /*U-equation*/
            {
                /*buoyancy*/
                ScalarCellField Sp = Scalar(0);
#ifdef LINEARIZE
                VectorCellField Sc = -gradf<true>(p,true);
#else
                VectorCellField Sc = Vector(0.0);
#endif
                if(buoyancy)
                    Sc += (rho - rho_ref) * g;
                /*solve*/
                VectorCellMatrix M;
#ifdef LINEARIZE
                M = transport(U, Fc, F, mu, velocity_UR, Sc, Sp, &rho, &rhof);
#else
                {
                    TensorCellField fq = mul(Fc,U) + TensorCellField(Constants::I_T) * p;
                    VectorCellField q = rho * U;
                    M = divf(fq,false,&F,&q,&lambdaMax) - lap(U,mu) - src(U,Sc,Sp);
                    M.cF = &U;
                    addTemporal<1>(M,velocity_UR,&rho,&rhof);
                }
#endif
                Solve(M);
            }

            /*T-equation*/
            {
                ScalarCellMatrix M;
#ifdef LINEARIZE
                M = transport(T, Fc, F, mu * iPr, t_UR, &rho, &rhof);
#else
                {
                    ScalarCellField imu = mu * iPr;
                    VectorCellField fq = Fc * T;
                    ScalarCellField q = rho * T;
                    M = divf(fq,false,&F,&q,&lambdaMax) - lap(T,imu);
                    M.cF = &T;
                    addTemporal<1>(M,t_UR,&rho,&rhof);
                }
#endif
                Solve(M);
            }

            /*courant number*/
            if(Controls::state != Controls::STEADY) {
                if(MP::printOn) {
                    Vector courant = Mesh::calc_courant(U,Controls::dt);
                    MP::printH("Courant number: Max: %g Min: %g Avg: %g\n",
                        courant[0], courant[1], courant[2]);

                    Scalar mass, energy;
                    {
                        ScalarCellField sf = rho * Mesh::cV;
                        mass = reduce_sum(sf);
                        //potential + kinetic + internal
                        sf = gh +
                             0.5 * magSq(U) +
                             pow((p + p_ref) / P0, R / cp) * T * cv;
                        sf = rho * Mesh::cV * sf;
                        energy = reduce_sum(sf);
                    }
                    MP::printH("Mass loss: %g Energy loss %g\n",
                        (mass0 - mass) / mass0,
                        (energy0 - energy) / energy0);
                }
            }

            /*save perturbation values*/
            rho -= rho_ref;
            T -= T0;
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
