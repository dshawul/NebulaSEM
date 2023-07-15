#include "calc_walldist.h"
#include "solve.h"

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
        const VectorCellField g = gradf(phi,true);
        yWall = sqrt((g & g) + 2 * phi) - mag(g);
    }
    /*write it*/
    yWall.write(step);
}
