#include "potential.h"
#include "solve.h"
#include "iteration.h"

using namespace std;

/**
  \verbatim
  Potential flow solver
  ~~~~~~~~~~~~~~~~~~~~~
  In potential flow the velocity field is irrotational (vorticity = curl(U) = 0).
  This assumption fails for boundary layers and wakes that exhibit strong vorticity,
  but the theory can still be used to initialize flow field for complex simulations.

  The potential flow assumption is that velocity is the gradient of a scalar field,
  which is the velocity potential (phi) 
         U = grad(phi)
         curl(U) = curl(grad(phi))
         curl(U) = 0
  where the last step is possible due to the vector identity curl(grad(phi)) = 0.
  Hence, defining the velocity as a gradient of a scalar ensures vorticity is zero.
  Let us now take the divergence instead as
         U = grad(phi)
         div(U) = div(grad(phi))
         div(U) = lap(phi)
  Since div(U)=0 for incompressible, the poisson equation becomes a laplace equation
         lap(phi) = 0

  What the solver does
  ~~~~~~~~~~~~~~~~~~~~~
  Given an initial velocity field (Ua) that does not satisfy continuity (div(Ua) != 0),
  we can correct Ua with the pressure gradient to get a divergence free velocity field U
         U = Ua - grad(p)
         div(Ua - grad(p)) = div(U) = 0
         div(Ua) = div(grad(p))
  If Ua is irrotational, i.e. curl(Ua) = 0 and Ua = grad(phia), then so is U because
         U = Ua - grad(p)
         U = grad(phia) - grad(p)
         U = grad(phia - p)
  where U = grad(phi) such that phi = phia - p.
  \endverbatim
 */
void potential(istream& input) {
    /*Solver specific parameters*/
    Int n_ORTHO = 0;

    /*potential options*/
    Util::ParamList params("potential");
    params.enroll("n_ORTHO", &n_ORTHO);

    /*read parameters*/
    Util::read_params(input,MP::printOn);

    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {

        VectorCellField U("U", READWRITE);
        ScalarCellField p("p", READWRITE);

        /*Time loop*/
        for (Iteration it(ait.get_step()); it.start(); it.next()) {
            const ScalarCellField divU = divf(U);
            const ScalarCellField one = Scalar(1);

            /*solve pressure poisson equation for correction*/
            for (Int k = 0; k <= n_ORTHO; k++)
                Solve(lap(p, one, true) == divU);

            /*correct velocity*/
            U -= gradi(p);
            applyExplicitBCs(U, true);
        }
    }
}
