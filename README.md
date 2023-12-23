[![NebulaSEM](https://github.com/dshawul/NebulaSEM/actions/workflows/NebulaSEM.yml/badge.svg)](https://github.com/dshawul/NebulaSEM/actions/workflows/NebulaSEM.yml)
[![APIdoc](https://github.com/dshawul/NebulaSEM/actions/workflows/APIdoc.yml/badge.svg)](https://github.com/dshawul/NebulaSEM/actions/workflows/APIdoc.yml)
[![JOSS](https://github.com/dshawul/NebulaSEM/actions/workflows/paper.yml/badge.svg)](https://github.com/dshawul/NebulaSEM/actions/workflows/paper.yml)

## NebulaSEM

NebulaSEM is an experimental finite volume (FV) and discontinous galerkin spectral element (dGSEM) code 
for solving partial differential equations (PDEs) of fluid dynamics. It comes with solvers for compressible
and incompressible flow, and also provides infrastructure for writing PDE solvers easily with support 
for adaptive mesh refinement (AMR). The primary focus of the software is research in a high-order 
non-hydrostatic dynamical core for atmospheric simulations. Several examples are provided to demonstrate this
capability.


### Build and install

Clone the repository

    git clone https://github.com/dshawul/NebulaSEM.git

To build and install NebulaSEM

    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=.. ..
    make && make install

Additional options to enable fine-grained parallelization on CPUs and GPUs

    -DUSE_ACC=ON
    -DUSE_OMP=ON

This will install tools for pre-processing, solution and post-processing.
The tool `mesh` generates the grid, `prepare` does various pre- and post-processing,
and several other binaries for solving PDEs e.g. euler, convection etc

### Requirements
- [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for domain decomposition
- [OpenMPI](http://open-mpi.org) or other MPI library compatible with your C++ compiler
- [CMake](https://cmake.org) for makefile generation
- [GCC](https://gcc.gnu.org/) or other C++ compiler with atleast C++17 standard support

Optional packages
- [DOxygen](http://www.doxygen.nl/) for API documentation
- [ParaView](https://www.paraview.org/) for visualization
- [NVHPC](https://developer.nvidia.com/hpc-sdk) or other OpenACC compiler

### Testing

A testing script `test.sh` is provided. By default it runs the lid-driven test case
under `examples/cavity`..

To run a test case, execute the script specifying the number of MPI ranks if greater than 1, the
test case name.

    ./test.sh -n 2 -c examples/atmo/advection-leveque/bubble

The path to the test case should point to the grid file, in this case `bubble` not to the directory itself.

    Usage: ./test.sh [options]
    
       -n,--np       Number of processors to use.
       -c,--case     Path to grid file name that is under a test case directory.
       -b,--bin-path Path to binaries: mesh, prepare and solvers.
       -s,--steps    Number of time steps, which overwrites the one in control file.
       -h,--help     Display this help message.



#### Lid-driven cavity flow

This test case uses the Pressure Implicit Splitting of Operators (PISO) solver for incompressible
flow at low Reynolds number i.e. no turbulence.

    $ ./test.sh -n 1 -c examples/cavity/cavity-amr

This will generate a `run-examples-cavity-amr` directory in which you can find the results including VTK
files for visualization by paraview.

Here are images of the decompostion using METIS with 12 mpi ranks, and the magnitude of
velocity plots.

<p align="center">
  <img width="300px" src="./images/cavity-decomp.png"/>
  <img width="300px" src="./images/cavity-velocity.png"/>
</p>

#### Pitz-Daily test case

Another test case, namely the Pitz and Daily, solved using LES is shown below.
You can see the formation of eddies at the backward facing step and later convection towards
the outlet.

<p align="center">
  <img width="900px" src="./images/pitz-velocity.gif"/>
</p>

The same test case simulated with the ke turbulence model is shown below. It is a Reynolds-average
turbulence scheme so only mean state is displayed.

<p align="center">
  <img width="900px" src="./images/pitz-velocity-ke.png"/>
</p>

#### Rising thermal bubble

This is a popular test case for numerical weather prediction models that solve the Euler equations
using explicit time-stepping unlike other CFD applications that often use implicit solvers.
Moreover this test cases uses Discontinous Galerkin method (spectral-element version) on hexahedral
grids, and adaptive mesh refinement. Thanks to my postdoc supervisor Francis X. Giraldo, from 
whom I learned this stuff!

A thermal bubble of gaussian distribution rises up due to bouyancy, while deforming on the way,
and collides with the top boundary.

<p align="center">
  <img width="500px" src="./images/rtb-temp.gif"/>
</p>

## Contribution

Users are welcome to suggest and implement new features such as new solvers, turbulence models, new
test cases, as well as report bugs or issues they encounter. Please feel free to [open an issue on this repository](https://github.com/dshawul/NebulaSEM/issues/new/choose) 
describing your desired change/bug-fix. Pull requests are also welcome!

