#include "solve.h"
#include "iteration.h"
#include "wrapper.h"

/**
  \verbatim
  Diffusion solver
  ~~~~~~~~~~~~~~~~
  Solver for pdes of parabolic heat equation type:
        dT/dt = lap(T,DT)
  \endverbatim
 */
void diffusion(std::istream& input) {
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
  Main application entry point for diffusion solver.
  \endverbatim
 */
int main(int argc, char* argv[]) {
   MP mp(argc, argv);
   Solver::Initialize(argc, argv);
   diffusion(Solver::input);
   Solver::Finalize();
   return 0;
}
