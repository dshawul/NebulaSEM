#include "solve.h"
#include "iteration.h"
#include "calc_walldist.h"
#include "wrapper.h"

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
void walldist(std::istream& input) {
    /*Solver specific parameters*/
    Int n_ORTHO = 0;

    /*walldist options*/
    Util::ParamList params("walldist");
    params.enroll("n_ORTHO", &n_ORTHO);
    Util::read_params(input,MP::printOn);

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {
        Mesh::calc_walldist(ait.get_step(), n_ORTHO);
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
   walldist(Solver::input);
   Solver::Finalize();
   return 0;
}
