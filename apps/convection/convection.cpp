#include "solve.h"
#include "iteration.h"
#include "wrapper.h"

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
