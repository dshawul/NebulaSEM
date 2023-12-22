# Setting up an experiment 

The contents of an experiment directory should look like this

```
├── atmo
│   ├── acoustic-sphere
│   │   ├── controls
│   │   ├── p0.txt
│   │   ├── sphere
│   │   ├── T0.txt
│   │   └── U0.txt
│   ├── acoustic-sphere-amr
│   │   ├── controls
│   │   ├── p0.txt
│   │   ├── sphere
│   │   ├── T0.txt
│   │   └── U0.txt
```

The `controls` file specifies parameters such as the solver name, start and end of timestep,
turbulence model, refinement criteria etc

```
general
{
######################################
#  mesh file name and solver type
######################################
    solver                       piso
    mesh                         grid

######################################
#  Fluid properties specific to solver
######################################
    rho                          1
    viscosity                    1e-1
#############################
#  Time increment
#############################

    state                        TRANSIENT
    start_step                   0
    end_step                     1000
    write_interval               100
    dt                           0.005
...
```

The file `sphere` for this particular examples specifies how to generate the grid.
For this particular cases, it is a cubed-sphere grid.
The file contains a list of all vertices, and specifies how they connect to form
hexahedral elements, boundary faces and their names etc.

Note that cubed-sphere grids are extruded to the surface of a sphere, unlike other grids.
When the keyword `delete` is specified as a boundary name, the face is deleted. This is useful
for running 2D simulations on a 3D grid.

```
# List of all vertices
# Number of vertices followed by (x y z)
#

16
{
-100.0 -100.0 -100.0
 100.0 -100.0 -100.0
 100.0  100.0 -100.0
-100.0  100.0 -100.0
-100.0 -100.0  100.0
 100.0 -100.0  100.0
 100.0  100.0  100.0
-100.0  100.0  100.0

-200.0 -200.0 -200.0
 200.0 -200.0 -200.0
 200.0  200.0 -200.0
-200.0  200.0 -200.0
-200.0 -200.0  200.0
 200.0 -200.0  200.0
 200.0  200.0  200.0
-200.0  200.0  200.0
}

# meshing
# Three types of divisions :
#     a) linear b) geometric c) wall 
#
#     8{0 1 2 3 4 5 6 7} linear    3{10 10 1}
#     8{0 1 2 3 4 5 6 7} geometric 3{10 10 1} 3{0.25 0.25 1}
#     8{0 1 2 3 4 5 6 7} wall      3{10 10 1} 3{0.25 0.25 1}
#

8{0 1 2 3 8 9 10 11} linear 3{10 10 1}
8{4 5 6 7 12 13 14 15} linear 3{10 10 1}
8{0 3 7 4 8 11 15 12} linear 3{10 10 1}
8{1 2 6 5 9 10 14 13} linear 3{10 10 1}
8{0 1 5 4 8 9 13 12} linear 3{10 10 1}
8{3 2 6 7 11 10 14 15} linear 3{10 10 1}

# boundaries
# All faces contained in the specified plane will be included.
# So it is enough to specify 3 vertices.
#    top   4{3 2 6 7}
#    top   3{3 2 6}
#

delete  4{8 9 10 11} 4{12 13 14 15} 4{8 11 15 12} 4{9 10 14 13} 4{8 9 13 12} 4{11 10 14 15} 
        4{0 1 2 3} 4{4 5 6 7} 4{0 3 7 4} 4{1 2 6 5} 4{0 1 5 4} 4{3 2 6 7}

```

Other files such as `p0.txt`, `U0.txt` and `T0.txt` specify initial conditions for the fields
used by a specific solver, `euler` solver for this particular example. The format of the files
is as follows
```
size 1
internal 1
{
    cosine 0 0.5    500 50 350   250 1000 250
}
boundary 3
{
    top {
        type NEUMANN
    }
    bottom {
        type NEUMANN
    }
    sides {
        type NEUMANN
    }
}

```
The `size` parameter specifies the size of the field, scalars take 1, and vectors take 3.
The `internal` entry specifies how to initialize internal fields. Here, you can list all the
field values taken by all cells in the grid, which is the most common way of initialization.

For ideal test case runs, additional ways of intitialization are provided with special keyworks.
The details of the implementations can be seen in `src/field/filed.h`
```
uniform - for uniform initialization with a given value
cosine - cosine bell distribution as used in the SRTB test case
gaussian- gaussian bell
linear - linear bell
hydrostatic - initialize with hydrostatic equilibrium
```

