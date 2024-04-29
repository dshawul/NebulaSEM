[![NebulaSEM](https://github.com/dshawul/NebulaSEM/actions/workflows/NebulaSEM.yml/badge.svg)](https://github.com/dshawul/NebulaSEM/actions/workflows/NebulaSEM.yml)
[![APIdoc](https://github.com/dshawul/NebulaSEM/actions/workflows/APIdoc.yml/badge.svg)](https://github.com/dshawul/NebulaSEM/actions/workflows/APIdoc.yml)
[![JOSS](https://github.com/dshawul/NebulaSEM/actions/workflows/paper.yml/badge.svg)](https://github.com/dshawul/NebulaSEM/actions/workflows/paper.yml)
[![DOI](https://zenodo.org/badge/28114781.svg)](https://zenodo.org/doi/10.5281/zenodo.11088209)

## NebulaSEM

NebulaSEM is an experimental finite volume (FV) and discontinuous Galerkin spectral element (dGSEM) code 
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
- [GCC](https://gcc.gnu.org/) or other C++ compiler with at least C++17 standard support

Optional packages
- [DOxygen](http://www.doxygen.nl/) for API documentation
- [ParaView](https://www.paraview.org/) for visualization
- [NVHPC](https://developer.nvidia.com/hpc-sdk) or other OpenACC compiler


### Documentation
The API documentation built with Doxygen can be found [here](https://dshawul.github.io/html/index.html).
Also READMEs are provided for writing new solvers under `apps/` directory, and for setting up
test cases under `examples/` directory.

### Testing

A testing script `test.sh` is provided. By default it runs the lid-driven test case
under `examples/cavity` using the binaries installed with the `make install` command.

To run a test case, execute the `test.sh` script specifying the number of MPI ranks, 
if greater than, 1 and the test case name.

    ./test.sh -n 2 -c examples/atmo/advection-leveque/bubble

If the binaries are installed other than `$PWD/bin`, pass the installation location to it using the `--bin-path` option

    ./test.sh -n 2 --bin-path /usr/local/bin -c examples/atmo/advection-leveque/bubble

The path to the test case should point to the grid file, in this case `bubble`, not to the directory itself.

    Usage: ./test.sh [options]
    
       -n,--np       Number of MPI processes to launch.
       -c,--case     Path to grid file name that is under a test case directory.
       -b,--bin-path Path to binaries: mesh, prepare and solvers.
       -s,--steps    Number of time steps, which overwrites the one in control file.
       -h,--help     Display this help message.


#### Lid-driven cavity flow

This test case uses the Pressure Implicit Splitting of Operators (PISO) solver for incompressible
flow at low Reynolds number i.e. no turbulence.

    $ ./test.sh -n 1 -c examples/cavity-amr/cavity

This will generate a `run-examples-cavity-amr` directory in which you can find the results including VTK
files for visualization by Paraview.

Here are images of the decomposition using METIS with 12 mpi ranks, and the magnitude of
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

The same test case simulated with the k-e turbulence model is shown below. It is a Reynolds-average
turbulence scheme so only mean state is displayed.

<p align="center">
  <img width="900px" src="./images/pitz-velocity-ke.png"/>
</p>

#### Rising thermal bubble

This is a popular test case for numerical weather prediction models that solve the Euler equations
using explicit time-stepping unlike other CFD applications that often use implicit solvers.
Moreover this test cases uses discontinuous Galerkin method (spectral-element version) on hexahedral
grids, and adaptive mesh refinement. Thanks to my postdoc supervisor Francis X. Giraldo, from 
whom I learned this stuff!

A thermal bubble of Gaussian distribution rises up due to buoyancy, while deforming on the way,
and collides with the top boundary.

<p align="center">
  <img width="500px" src="./images/rtb-temp.gif"/>
</p>

### A note about MPI/OpenMP/OpenACC parallelization

NebulaSEM is able to exploit multi-core CPUs either using a pure MPI approach or a hybrid MPI+OpenMP approach.
CFD codes often utilize a pure MPI approach because that scales better than OpenMP. However, when extreme scalability
on supercomputers is required, the hybrid MPI + OpenMP parallelization can become beneficial by reducing inter-process
communication. To compile NebulaSEM for a hybrid MPI+OpenMP

    mkdir build-omp && cd build-omp
    cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=.. -DUSE_OMP=ON ..
    make && make install

The number of threads for OpenMP is controlled by the `OMP_NUM_THREADS` environment variables.
Then we can run the lid-driven cavity flow as usual specifying only 1 MPI rank and 2 OpenMP threads
for a pure OpenMP approach as:

    $ export OMP_NUM_THREADS=2
    $ ./test.sh -n 1 -c examples/cavity/cavity

Or 2 mpi ranks + 2 openmp threads per rank for a hybrid MPI+OpenMP approach

    $ export OMP_NUM_THREADS=2
    $ ./test.sh -n 2 -c examples/cavity/cavity

Or 2 mpi ranks with 1 threads per rank for a pure MPI approach

    $ export OMP_NUM_THREADS=1
    $ ./test.sh -n 2 -c examples/cavity/cavity

Note that to obtain significant speedups from either MPI or OpenMP, the problem size should be large enough.
For the lid-driven cavity flow we can increase the problem size to `100x100` by editing `examples/cavity/vaity`
```diff
-8{0 1 2 3 4 5 6 7} linear 3{20 20 1}
+8{0 1 2 3 4 5 6 7} linear 3{100 100 1}
```
Run single MPI rank run
```
$ export OMP_NUM_THREADS=1
$ ./test.sh -n 1 -c examples/cavity/cavity 
...
9018 [0] Time 5.000000
9025 [0] SYMM-FULL-SSOR-PCG :Iterations 1 Initial Residual 1.47888e-11 Final Residual 1.14971e-11
9027 [0] SYMM-FULL-SSOR-PCG :Iterations 1 Initial Residual 1.44703e-11 Final Residual 1.12991e-11
9073 [0] Exiting application run with 1 processes
```
Takes about 9073 milliseconds.

Run 2-mpi ranks
```
$ ./test.sh -n 2 -c examples/cavity/cavity
....
4623 [0] Time 5.000000
4626 [0] SYMM-FULL-SSOR-PCG :Iterations 1 Initial Residual 3.55206e-11 Final Residual 3.16645e-11
4628 [0] SYMM-FULL-SSOR-PCG :Iterations 1 Initial Residual 3.38580e-11 Final Residual 3.00837e-11
4672 [0] Exiting application run with 2 processes
```
Takes about 4672 milliseconds for a speedup of  `1.94x` out of 2 which is good.

----
Lets do the same with the OpenMP implementation. Make sure to use the `build-omp` binaries by doing:
```
cd build-omp && make install
```
We can now run the test case with 1 mpi rank
```
$ export OMP_NUM_THREADS=2
$ ./test.sh -n 1 -c examples/cavity/cavity
....
6686 [0] Time 5.000000
6691 [0] SYMM-FULL-SSOR-PCG :Iterations 1 Initial Residual 1.47888e-11 Final Residual 1.14971e-11
6693 [0] SYMM-FULL-SSOR-PCG :Iterations 1 Initial Residual 1.44703e-11 Final Residual 1.12991e-11
6739 [0] Exiting application run with 1 processes
```
It took about 6739 milliseconds. It is slower than the MPI implementation but still faster than 
the serial implementation with a speedup of `1.34x` times.

## Contribution

Users are welcome to suggest and implement new features such as new solvers, turbulence models, new
test cases, as well as report bugs or issues they encounter. Please feel free to [open an issue on this repository](https://github.com/dshawul/NebulaSEM/issues/new/choose) 
describing your desired change/bug-fix. Pull requests are also welcome!

