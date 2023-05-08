#include "solve.h"
#include "iteration.h"
#include "wrapper.h"

/**
  \verbatim
  Transport equation solver
  ~~~~~~~~~~~~~~~~~~~~~~~~~
  Given a flow field (U), the solver determines the distribution of a 
  scalar by convection and diffusion.
     dT/dt + div(T,F,DT) = lap(T,DT)
  \endverbatim
 */
void transport(std::istream& input) {
    /*Solver specific parameters*/
    Scalar DT = Scalar(1.0e-4);
    Scalar t_UR = Scalar(1);

    /*transport options*/
    Util::ParamList params("transport");
    params.enroll("t_UR", &t_UR);
    params.enroll("DT", &DT);

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
        ScalarCellField mu = DT;
        for (; !it.end(); it.next()) {
            ScalarCellMatrix M;
            M = transport(T, Fc, F, mu, t_UR);
            Solve(M);
        }
    }
}

/**
  \verbatim
  Main application entry point for transport solver.
  \endverbatim
 */
int main(int argc, char* argv[]) {
   MP mp(argc, argv);
   Solver::Initialize(argc, argv);
   transport(Solver::input);
   Solver::Finalize();
   return 0;
}
