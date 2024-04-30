---
title: 'NebulaSEM: A high-order discontinuous Galerkin spectral element code for atmospheric modeling'

tags:
  - discontinuos Galerking spectral element
  - finite volume
  - atmospheric modelling
  - global circulation model
  - computational fluid dynamics
  - adaptive mesh refinement
  - C++ expression templates

authors:
  - name: Daniel S. Abdi
    orcid: 0000-0003-0671-2941
    corresponding: true
    affiliation: 1

affiliations:
 - name: Cooperative Institute for Research in Environmental Sciences
   index: 1

date: 20 December 2023

bibliography: paper.bib
---

# Summary

NebulaSEM is an advanced computational fluid dynamics (CFD) model employing the high-order discontinuous Galerkin 
spectral element (dGSEM) method. Originally designed for generic CFD applications with polyhedral elements on unstructured grids, 
NebulaSEM has recently evolved into an atmospheric simulation code. This paper presents the key features, capabilities, 
and applications of NebulaSEM, highlighting its high-order discretization, oct-tree cell-based adaptive mesh refinement (AMR), 
support for new solver development, turbulence models, and parallelization strategies. The software addresses the need for 
high-fidelity simulations in diverse areas, offering an efficient and scalable solution with a focus on atmospheric modeling.


# NebulaSEM Features

## High-order dGSEM discretization
NebulaSEM supports arbitrarily high-order discontinuous Galerkin spectral element discretization of PDEs besides
finite-volume discretization, which is a subset of dGSEM with the lowest polynomial order of zero. The spectral element discretization
is significantly more efficient than standard discontinuous Galerkin on unstructured grids because the former exploits the tensor-product nature of
computations to reduce computations from O(N^3) to O(3N). The dGSEM possesses several desirable characteristics [@abdi8] 
such as: high-order accuracy, larger geometrical flexibility compared to global spectral methods, high scalability due to 
the high arithmetic intensity per element, suitability for GPU acceleration, and support for both h- and p- refinement.

## Oct-tree cell-based AMR
Modeling of the atmosphere is challenging in that a range of spatial and temporal scales are involved [@stefanick1981].
Adaptive Mesh Refinement (AMR) provides the tool to address multi-scale phenomena efficiently by focusing resources
where they are needed. NebulaSEM implements oct-tree cell-based AMR using the forest-of-octrees approach pioneered in [@p4est].
The `AmrIteration` class provides a high-level interface to enable AMR for any solver written using the library.
A single loop enclosing the timestep iterations and declaration of fields involved in the PDE is enough to provide AMR 
capability for any solver. The details of regridding the domain, memory management, resizing and transferring fields in 
a conservative manner etc. are all taken care of behind the scenes by the library.

## Support for new solver development
NebulaSEM offers a range of operators for spatial and temporal discretization, streamlining the development of 
solvers for Partial Differential Equations (PDEs). The provided code snippet includes an example solver for 
the advection-diffusion equation.

```C++
void transport() {
  /*AMR iteration loop with object (ait)*/
    for (AmrIteration ait; !ait.end(); ait.next()) {
         VectorCellField U("U", READWRITE);      /*Velocity field over the grid*/
         ScalarCellField T("T", READWRITE);      /*Scalar field*/
         ScalarFacetField F = flx(U);            /*Compute flux field*/
         ScalarCellField mu = 1;                 /*Diffusion parameter*/
         /*Time loop with support for deferred correction */
         for (Iteration it(ait.get_step()); !it.end(); it.next()) {
             ScalarCellMatrix M;                 /*Matrix for the PDE discretization*/
             M = div(T,U,F,&mu) - lap(T,mu);     /*Divergence & Laplacian terms*/
             addTemporal<1>(M);                  /*Add temporal derivative*/
             Solve(M);                           /*Solve the matrix */
         }
    }
}
```

Spatial operators within dGSEM encompass divergence, gradient, Laplacian, and more. Temporal discretization is 
accomplished through explicit and implicit schemes, including first-order Euler explicit and implicit schemes, 
linear multi-step methods such as Adams-Moulton and Adams-Bashforth, the Runge-Kutta method up to 4th order, 
and fully-implicit Backward Differencing (BDF) methods.

## Turbulence models
The software includes turbulence models for high-Reynolds CFD simulations including a suite of Reynolds Averaged Navier Stokes (RANS) 
[@henk1972] and Large Eddy Simulation (LES) models [@les]. The list includes a mixing-length model, k-epsilon, k-omega, RNG k-epsilon, RNG k-omega
and the Smagornisky-Lilly LES model. These turbulence models have been utilized to evaluate the aerodynamic roughness of the 
built environment and complex terrain in [@abdi5].

## Parallelization with MPI+OpenMP/OpenACC
NebulaSEM achieves scalability on supercomputers through a combination of coarse- and fine-grained parallelism. 
The Message Passing Interface (MPI) is employed for distributed computing, while directive-based threading libraries 
such as OpenMP for CPUs and OpenACC for GPUs optimize fine-grained parallelism, minimizing communication overhead.
Details of parallelization of the linear system of equations solvers, such as the preconditioned gradient solver, 
can be found in [@abdi9]. In addition, NebulaSEM implements a unique approach of asynchronous parallelization that is
not commonly found in CFD applications.

Efficient GPU implementation of dGSEM is achieved through offloading of all field computations to the GPU [@abdi8], 
using a memory pool to recycle previously allocated memory by fields that went out of scope, utilization of
managed memory to simplify the data transfer logic between CPU and GPU etc.

# Showcases
We showcase the capabilities of NebulaSEM through two example applications:

a) Generic CFD Application:
We employ the Pressure-Implicit Splitting of Operators (PISO) solver for incompressible fluid flow to solve the 
Pitz-Daily problem [@pitz1983], utilizing the Smagornisky-Lilly large eddy simulation. The instantaneous velocity 
profiles are depicted in \autoref{fig:pitz-daily}.

b) Atmospheric Simulation:
NebulaSEM can serve as a non-hydrostatic dynamical core for atmospheric simulations.
To evaluate the dynamical core, we solve the rising thermal bubble problem [@robert1993] using the non-hydrostatic 
Euler equations solver. The test involves adaptive mesh refinement and discontinuous Galerkin spectral element 
discretization with polynomial order of 4. The result is depicted in \autoref{fig:srtb}

![Simulation results for the Pitz-Daily problem [@pitz1983] that aims to evaluate the effect of combustion
on mixing layer growth. A snapshot of large eddy simulation (LES) results using finite-volume method of NebulaSEM is presented.
\label{fig:pitz-daily}](pitz-daily.png)

![Simulation results for the Robert Rising Thermal Bubble (RRTB) problem [@robert1993] using oct-tree cell-based AMR
with two levels of refinement. The discontinuous Galerkin spectral element method is used with polynomial order of 4.
\label{fig:srtb}](srtb-amr.png){width=30%, height=30%}

# Statement of need

Over the years, NebulaSEM has undergone a transformation from being a purely Computational Fluid Dynamics (CFD) application 
to primarily serve as an atmospheric simulation (AtmoSim) code. Consequently, it addresses specific needs inherent to both types of
applications, as outlined below.

Many CFD codes prioritize robustness over other considerations. As a result, they often utilize first or at best second-order accurate 
finite volume methods on unstructured grids. In contrast, NebulaSEM offers high-order discretization characterized 
by high accuracy and minimal dissipation. This feature proves invaluable for tasks such as accurately capturing shocks and
discontinuities, conducting highly accurate large eddy simulations, simulating turbulent flows with precision, and 
facilitating high-fidelity aeroacoustic simulations.

While high-order methods are more commonly associated with atmospheric modeling, such endeavors often rely on
finite-difference discretization, which is less geometrically flexible, or global spectral methods on latitude-longitude grids, 
which lack scalability. NebulaSEM provides geometrical flexibility due to its CFD roots, allowing atmospheric simulations on any type of grid.
The element-based design, as opposed to global spectral-methods, ensures high-scalability for large scale simulations.
The high-order dGSEM discretization it employs delivers the accuracy necessary for achieving high-fidelity atmospheric simulations. 
In addition, NebulaSEM incorporates dynamic adaptive mesh refinement capabilities, a feature not commonly found in traditional atmospheric modeling approaches.


# Acknowledgements

I would like to thank my postdoc supervisor Francis X. Giraldo who introduced me to the high-order discontinuous Galerkin method.

# References
