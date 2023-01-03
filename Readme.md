### Solver

Solver is an experimental finite volume and discontinous galerkin code 
for solving partial differential equations.

### Build

To build Solver type

    make

followed by

    make install

This will install three tools for pre-processing, solution and post-processing.
The tool 'mesh' generates the grid, 'solver' does the solution and 'prepare' does
various post-processing.

### Testing

A testing script `test.sh` is provided. By default it runs the lid-driven test case
under `examples/cavity`, but you can modify the script to run any test case.
To run the test case execute is specifying the number of processors if >1 or a 
different test case.

    Usage: ./test.sh [options]
    
       -n,--np       Number of processors to use.
       -c,--case     Path to grid name under a test case directory.
       -h,--help     Display this help message.

    $ ./test.sh -n 1 -c examples/cavity/cavity
    $ ./test.sh

This will generate a `run1` directory in which you can find the results including VTK
files needed for visualization by paraview.
