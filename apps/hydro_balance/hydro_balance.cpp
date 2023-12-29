#include "solve.h"
#include "iteration.h"
#include "properties.h"
#include "wrapper.h"

/**
  \verbatim
  Hydrostatic
  ~~~~~~~~~~~~~
     Solver for hydrostatic balance
         grad(p) = -rho*g
     Using gravitational potential theory
        div(grad(p)) = div(-rho*g)
  \endverbatim
*/
void hydro_balance(std::istream& input) {
    /*Solver specific parameters*/
    Int n_ORTHO = 0;

    /*fluid properties*/
    {
        Util::ParamList params("general");
        Fluid::enroll(params);
        params.read(input);
    }

    /*hydro_balance options*/
    Util::ParamList params("hydro_balance");
    params.enroll("n_ORTHO", &n_ORTHO);

    /*read parameters*/
    Util::read_params(input,MP::printOn);

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {
        ScalarCellField p("p", READWRITE);

        for (Iteration it(ait.get_step()); it.start(); it.next()) {
            const ScalarCellField one = Scalar(1);
            const VectorCellField rhog = Fluid::density * Controls::gravity;
            const ScalarCellField ndivRhoG = -divf(rhog);

            /*solve poisson equation*/
            for (Int k = 0; k <= n_ORTHO; k++)
                Solve(lap(p, one, true) == ndivRhoG);
        }
    }
}

/**
  \verbatim
  Main application entry point for hydrostatic balance solver.
  \endverbatim
 */
int main(int argc, char* argv[]) {
   MP mp(argc, argv);
   Solver::Initialize(argc, argv);
   hydro_balance(Solver::input);
   Solver::Finalize();
   return 0;
}
