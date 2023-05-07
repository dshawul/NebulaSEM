#include "walldist.h"
#include "solve.h"
#include "iteration.h"

using namespace std;

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
void walldist(istream& input) {
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
  Calculate wall distance at given time step
 */
void Mesh::calc_walldist(Int step, Int n_ORTHO) {
    ScalarCellField& phi = yWall;
    /*poisson equation*/
    {
        const ScalarCellField one = Scalar(1);
        for (Int k = 0; k <= n_ORTHO; k++)
            Solve(lap(phi, one, true) == -cV);
    }
    /*wall distance*/
    {
        const VectorCellField g = gradi(phi);
        yWall = sqrt((g & g) + 2 * phi) - mag(g);
    }
    /*write it*/
    yWall.write(step);
}
