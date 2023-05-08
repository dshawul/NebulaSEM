#include "solve.h"
#include "iteration.h"
#include "wrapper.h"

/**
  \verbatim
  Wave equation solver
  ~~~~~~~~~~~~~~~~~~~~
  Solver for pdes of hyperbolic wave equation type:
        d2T/dt2 = c^2 * lap(T)
  \endverbatim
 */
void wave(std::istream& input) {
    /*Solver specific parameters*/
    Scalar C2 = Scalar(1);
    Scalar t_UR = Scalar(1);

    /*wave options*/
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
  Main application entry point for wave solver.
  \endverbatim
 */
int main(int argc, char* argv[]) {
   MP mp(argc, argv);
   Solver::Initialize(argc, argv);
   wave(Solver::input);
   Solver::Finalize();
   return 0;
}
