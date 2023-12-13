#include "solve.h"
#include "iteration.h"
#include "wrapper.h"

//#define LINEARIZE

/**
  \verbatim
  Convection solver
  ~~~~~~~~~~~~~~~~~~~~~~~~~~
  Given a flow field (U), the solver determines the distribution of a 
  scalar by convection.
     dT/dt + div(T,F,0) = 0
  \endverbatim
 */
void convection(std::istream& input) {
    /*Solver specific parameters*/
    Scalar t_UR = Scalar(1);
    enum PROBLEM_INIT {
        NONE, LEVEQUE, LAURITZEN_0, LAURITZEN_1
    };
    PROBLEM_INIT problem_init = NONE;

    /*advection options*/
    Util::ParamList params("convection");
    params.enroll("t_UR", &t_UR);
    Util::Option* op = new Util::Option(&problem_init,
            {"NONE", "LEVEQUE", "LAURITZEN_0", "LAURITZEN_1"});
    params.enroll("problem_init", op);

    /*read parameters*/
    Util::read_params(input,MP::printOn);

    /*total scalar quanity of system*/
    Scalar scalar0, volume0;

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {

        VectorCellField U("U", READWRITE);
        ScalarCellField T("T", READWRITE);

        /**
        Wind field intialization for Scalar advection
        */
        auto init_wind_field = [&] (Scalar time, Scalar etime) {
            using namespace Mesh;
            using namespace Constants;
        
            Scalar radius = sphere_radius;
            Scalar period = etime;
            Scalar RoT = radius / period; 

            #pragma omp parallel for
            #pragma acc parallel loop
            forEach(cC,i) {
                Scalar x = cC[i][0], y = cC[i][1], z = cC[i][2];
                Scalar u,v;
        
                if(Mesh::is_spherical) {
                    Vector s = cart_to_sphere(cC[i]);
                    Scalar lat = s[1], lon = s[2]; 
                    Scalar lambda = lon - 2.0 * PI * time / period;

                    if(problem_init == LAURITZEN_0) {
                        u = 10.0 * RoT * pow(sin(lambda),2.0) * sin(2.0 * lat) * 
                               cos(PI * time / period) + 2.0 * PI * RoT * cos(lat);
                        v = 10.0 * RoT * sin(2.0 * lambda) * cos(lat) * cos(PI * time / period);
                    } else if (problem_init == LAURITZEN_1) {
                        u = -5.0 * RoT * pow(sin(0.5 * lambda),2.0) * sin(2.0 * lat) * pow(cos(lat),2.0) *
                               cos(PI * time / period) + 2.0 * PI * RoT * cos(lat);
                        v = 2.5 * RoT * sin(lambda) * pow(cos(lat),3.0) * cos(PI * time / period);
                    }

                    U[i] = wind_field(u,v,lat,lon);
                } else {
                    if(problem_init == LEVEQUE) {
                        u = pow(sin(PI * x), 2.0) * sin(2 * PI * y) * cos(PI * time / period);
                        v = -pow(sin(PI * y), 2.0) * sin(2 * PI * x) * cos(PI * time / period);
                    }

                    U[i] = Vector(u,v,0);
                }
            }
        };

        /*special initializations*/
        Iteration it(ait.get_step());
        VectorCellField Fc;
        ScalarFacetField F;
        ScalarFacetField lambdaMax;

        /*special initializations*/
        if(ait.start() && problem_init != NONE) {
            int step = it.get_step();
            if(step < 0) step = 0;
            Scalar ctime = step * Controls::dt;
            Scalar etime = Controls::end_step * Controls::dt;
            init_wind_field(ctime, etime);
            U.write(0);
        }
        Fc = flxc(U);
        F = flx(U);
        lambdaMax = cds(mag(U)) / 2;

        /*Compute total scalar*/
        if(ait.get_step() == 0)
        {
            ScalarCellField sf = T * Mesh::cV;
            scalar0 = reduce_sum(sf);
            volume0 = reduce_sum(Mesh::cV);
        }

        /*Time loop*/
        for (; !it.end(); it.next()) {

            /*overwrite wind field at all time steps*/
            if(problem_init != NONE) {
                Scalar ctime = it.get_step() * Controls::dt;
                Scalar etime = Controls::end_step * Controls::dt;
                init_wind_field(ctime, etime);
                Fc = flxc(U);
                F = flx(U);
                lambdaMax = cds(mag(U)) / 2;
            }

            /*solve*/
            ScalarCellMatrix M;
#ifdef LINEARIZE
            M = convection(T, Fc, F, t_UR, &lambdaMax);
#else
            {
                VectorCellField fq = Fc * T;
                M = divf(fq,false,&F,&T,&lambdaMax);
                M.cF = &T;
                addTemporal<1>(M,t_UR);
            }
#endif
            Solve(M);

            /*compute scalar loss*/
            if(MP::printOn) {
                Scalar scalar, volume;
                ScalarCellField sf = T * Mesh::cV;
                scalar = reduce_sum(sf);
                volume = reduce_sum(Mesh::cV);
                MP::printH("Scalar loss: %.12g Volume loss: %.12g\n",
                        (scalar0 - scalar) / scalar0,
                        (volume0 - volume) / volume0);
            }
        }
    }
}

/**
  \verbatim
  Main application entry point for convection solver.
  \endverbatim
 */
int main(int argc, char* argv[]) {
   MP mp(argc, argv);
   Solver::Initialize(argc, argv);
   convection(Solver::input);
   Solver::Finalize();
   return 0;
}
