#include "convection.h"
#include "solve.h"
#include "iteration.h"

using namespace std;

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
