# Writing new solvers

Writing new solvers for NebulaSEM is easy. Convenience wrappers are provided to make
the process comfortable. Lets discusss the advection-diffusion solver called `transport`

```C++
#include "solve.h"
#include "iteration.h"
#include "wrapper.h"
```
Include these header files that provide solvers for linear system of equations, iteration
objects for AMR and time stepping, and other wrappers.
 
```C++
/**
  \verbatim
  Transport equation solver
  ~~~~~~~~~~~~~~~~~~~~~~~~~
  Given a flow field (U), the solver determines the distribution of a 
  scalar by convection and diffusion.
     dT/dt + div(T,F,DT) = lap(T,DT)
  \endverbatim
 */
```
Describe the PDE we want to solve. For this solver we need a temporal discretization scheme,
the divergence and laplacian operators.

```C++
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
```
Declare and register parameters needed for the solver. The `controls` file can be used
to specify values for these parameters
```C++
    /*AMR iteration*/
    for (AmrIteration ait; !ait.end(); ait.next()) {
```
The AMR iteration loop that refines/coarsens the grid. When AMR is not required, this loop
runs once.
```C++
        VectorCellField U("U", READWRITE);
        ScalarCellField T("T", READWRITE);
```
The solver require a velocity vector field, and a scalar field for the transported tracer.
When the grid is modified by AMR, these fields are automatically resized and field values
transferred to the new grid.
```C++
        /*Time loop*/
        Iteration it(ait.get_step());
        VectorCellField Fc = flxc(U);
        ScalarFacetField F = flx(U);
        ScalarCellField mu = DT;
```
Compute additional fields such as the volume and surface flux.
```C++
        for (; !it.end(); it.next()) {
```
The time step iteration loop.
```C++
            ScalarCellMatrix M;
            M = transport(T, Fc, F, mu, t_UR);
            Solve(M);
```
Define the matrix `M` to store the coefficients of the discretized linear system of equations.
Since advection-diffusion equation is a common PDE used for example in several turbulence models etc.
a convenience wrapper `transport` is provided. Apply the `Solve` function on the matrix to get solutions
for the timestep.
```C++
        }
    }
}
```

```C++
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
```
The main function declares and initializes an MPI object to distribute the work
across nodes.
