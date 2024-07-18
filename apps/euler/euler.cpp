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
    enum PROBLEM_INIT {
        NONE, ISENTROPIC_VORTEX
    };
    PROBLEM_INIT problem_init = NONE;

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
    op = new Util::Option(&problem_init,
            {"NONE", "ISENTROPIC_VORTEX"});
    params.enroll("problem_init", op);

    /*read parameters*/
    Util::read_params(input,MP::printOn);

    /*total mass and energy of system*/
    Scalar mass0, energy0, volume0;

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {
        ScalarCellField p("p",READWRITE);
        VectorCellField U("U",READWRITE);
        ScalarCellField T("T",READWRITE);
        ScalarCellField rho("rho",READWRITE);
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

        /*special initializations*/
        if(ait.start() && problem_init != NONE) {
            /*isentropic vortex*/
            auto isentropic_vortex = [&]() {
                using namespace Constants;
                constexpr Scalar beta = 5;
                ScalarCellField r = mag(Mesh::cC);
                T = (-((p_gamma - 1) * beta * beta) / (8 * p_gamma * PI * PI)) * exp(1 - r*r);
                #pragma omp parallel for
                #pragma acc parallel loop
                for(Int i = 0; i < Mesh::gBCSfield; i++) {
                    U[i][0] += (beta / (2 * PI)) * exp((1 - r[i]*r[i])/2.0) * -Mesh::cC[i][1];
                    U[i][1] += (beta / (2 * PI)) * exp((1 - r[i]*r[i])/2.0) *  Mesh::cC[i][0];
                }
                p = pow(T + T0, p_gamma / (p_gamma - 1)) - P0;
                rho = (P0 / (R*(T + T0))) * pow((p + P0) / P0, 1 / p_gamma) - (P0 / (R*T0));
            };

            /*choose init type*/
            if(problem_init == ISENTROPIC_VORTEX)
                isentropic_vortex();

            /*write initial values*/
            T.write(0);
            U.write(0);
            p.write(0);
            rho.write(0);
        }

        /*buoyancy*/
        ScalarCellField p_ref, rho_ref, gh;
        if(buoyancy) {
            /*gravity*/
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

            /*reference hydrostatic state*/
            p_ref = P0 * pow(1.0 + gh / (cp * T0), cp / R);
            rho_ref = (P0 / (R*T0)) * pow(p_ref / P0, 1 / p_gamma);
        } else {
            /*gravity*/
            gh = Scalar(0.0);

            /*reference state*/
            p_ref = P0;
            rho_ref = (P0 / (R*T0));
        }

#if 0
        if(ait.start()) {
#endif
            /*calculate rho from p*/
            p += p_ref;
            Mesh::scaleBCs<Scalar>(p,p_ref,1.0);
            applyExplicitBCs(p_ref,true);
            applyExplicitBCs(p,true);

            /*calculate rho*/
            Mesh::scaleBCs<Scalar>(p,rho_ref,psi);
            applyExplicitBCs(rho_ref,true);

            rho = (P0 / (R*(T + T0))) * pow(p / P0, 1 / p_gamma);
            Mesh::scaleBCs<Scalar>(p,rho,psi);
            applyExplicitBCs(rho,true);
            rho -= rho_ref;
            rho.write(0);

            p -= p_ref;
#if 0
        } else {
            /*calculate p from rho*/
            rho += rho_ref;
            Mesh::scaleBCs<Scalar>(rho,rho_ref,1.0);
            applyExplicitBCs(rho_ref,true);
            applyExplicitBCs(rho,true);

            /*calculate p*/
            Mesh::scaleBCs<Scalar>(rho,p_ref,1.0/psi);
            applyExplicitBCs(p_ref,true);

            p = P0 * pow((rho*(T+T0)*R) / P0, p_gamma);
            Mesh::scaleBCs<Scalar>(rho,p,1.0/psi);
            applyExplicitBCs(p,true);
            p -= p_ref;

            rho -= rho_ref;
        }
#endif

        /*compute total mass and energy*/
        if(ait.get_step() == 0)
        {
            ScalarCellField sf = (rho + rho_ref) * Mesh::cV;
            mass0 = reduce_sum(sf);
            //potential + kinetic + internal
            sf = gh +
                 0.5 * magSq(U) +
                 pow((p + p_ref) / P0, R / cp) * (T + T0) * cv;
            sf = (rho + rho_ref) * Mesh::cV * sf;
            energy0 = reduce_sum(sf);
            volume0 = reduce_sum(Mesh::cV);
        }

        /*Time loop*/
        for (; !it.end(); it.next()) {
            /*add reference values*/
            rho += rho_ref;
            T += T0;

            /*fluxes*/
            Fc = flxc(rho * U);
            F = flx(rho * U);
            lambdaMax = ( cds(mag(U)) + cds(sqrt(p_gamma*R*T)) ) / 2;

            /*artificial viscosity*/
            if(diffusion) mu = rho * viscosity;
            else mu = Scalar(0);

            /*rho-equation*/
            rhof = rho;
            {
                ScalarCellMatrix M;
#ifdef LINEARIZE
                M = convection(rho, flxc(U), flx(U), pressure_UR, &lambdaMax);
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
                M = transport(U, Fc, F, mu, velocity_UR, Sc, Sp, &lambdaMax, &rho, &rhof);
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
                M = transport(T, Fc, F, mu * iPr, t_UR, &lambdaMax, &rho, &rhof);
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
                Vector courant = Mesh::calc_courant(U,Controls::dt);
                Scalar mass, energy, volume;
                {
                    ScalarCellField sf = rho * Mesh::cV;
                    mass = reduce_sum(sf);
                    //potential + kinetic + internal
                    sf = gh +
                         0.5 * magSq(U) +
                         pow((p + p_ref) / P0, R / cp) * T * cv;
                    sf = rho * Mesh::cV * sf;
                    energy = reduce_sum(sf);
                    volume = reduce_sum(Mesh::cV);
                }
                if(MP::printOn) {
                    MP::printH("Courant number: Max: %g Min: %g Avg: %g\n",
                        courant[0], courant[1], courant[2]);
                    MP::printH("Mass loss: %.12g Energy loss %.12g Volume loss %.12g\n",
                        (mass0 - mass) / mass0,
                        (energy0 - energy) / energy0,
                        (volume0 - volume) / volume0);
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
