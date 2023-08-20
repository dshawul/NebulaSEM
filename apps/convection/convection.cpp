#include "solve.h"
#include "iteration.h"
#include "wrapper.h"

//#define OVERWRITE_WIND

/**
Wind field intialization for Scalar advection
*/
#ifdef OVERWRITE_WIND
void init_wind_field(Scalar time, VectorCellField& U) {
    using namespace Mesh;
    using namespace Constants;

    Scalar radius = 6371220.0;
    Scalar period = 12 * 24 * 3600;
    //Scalar period = 5.0;
    Scalar RoT = radius / period; 
    forEach(cC,i) {
        Scalar x = cC[i][0], y = cC[i][1], z = cC[i][2];
        Vector s = cart_to_sphere(cC[i]);
        Scalar lat = s[1], lon = s[2]; 
        Scalar lambda = lon - 2.0 * PI * time / period;
#if 0
        Scalar u = 40.0, v = 0.0;
#elif 1
        Scalar u = 10.0 * RoT * pow(sin(lambda),2.0) * sin(2.0 * lat) * 
                   cos(PI * time / period) + 2.0 * PI * RoT * cos(lat);
        Scalar v = 10.0 * RoT * sin(2.0 * lambda) * cos(lat) * cos(PI * time / period);
#elif 0
        Scalar u = -5.0 * RoT * pow(sin(0.5 * lambda),2.0) * sin(2.0 * lat) * pow(cos(lat),2.0) *
                   cos(PI * time / period) + 2.0 * PI * RoT * cos(lat);
        Scalar v = 2.5 * RoT * sin(lambda) * pow(cos(lat),3.0) * cos(PI * time / period);
#elif 0
        Scalar u = pow(sin(PI * x), 2.0) * sin(2 * PI * y) * cos(PI * time / period);
        Scalar v = -pow(sin(PI * y), 2.0) * sin(2 * PI * x) * cos(PI * time / period);
#endif
        if(Mesh::is_spherical)
            U[i] = wind_field(u,v,lat,lon);
        else
            U[i] = Vector(u,v,0);
    }
}
#endif

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

    /*transport*/
    Util::ParamList params("convection");
    params.enroll("t_UR", &t_UR);

    /*read parameters*/
    Util::read_params(input,MP::printOn);

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {

        VectorCellField U("U", READWRITE);
        ScalarCellField T("T", READWRITE);
        VectorCellField Fc;
        ScalarFacetField F;

        /*Time loop*/
        Iteration it(ait.get_step());
#ifdef OVERWRITE_WIND
        Scalar time = it.get_step() * Controls::dt;
        init_wind_field(time, U);
#endif
        Fc = flxc(U);
        F = flx(U);
        for (; !it.end(); it.next()) {
#ifdef OVERWRITE_WIND
            Scalar time = it.get_step() * Controls::dt;
            init_wind_field(time, U);
            Fc = flxc(U);
            F = flx(U);
#endif
            ScalarCellMatrix M;
            M = convection(T, Fc, F, t_UR);
            Solve(M);
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
