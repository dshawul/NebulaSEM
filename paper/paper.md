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
spectral element method. Originally designed for generic CFD applications with polyhedral elements on unstructured grids, 
NebulaSEM has evolved into a powerful atmospheric simulation code. This paper presents the key features, capabilities, 
and applications of NebulaSEM, highlighting its high-order discretization, oct-tree cell-based adaptive mesh refinement (AMR), 
new solver development support, turbulence models, and parallelization strategies. The software addresses the need for 
high-fidelity simulations in diverse areas, offering an efficient and scalable solution with a focus on atmospheric modeling.


# NebulaSEM Features

## High-order dGSEM discretization
NebulaSEM supports arbitrarily high-order discontinuous Galerkin spectral element discretization of PDEs besides
finite-volume discretization, which is a subset of dGSEM with the lowest polynomial order of zero. The spectral element discretization
is significantly more efficient than dG on unstructured grid because it exploits the tensor-product nature of
computations to reduce computations from O(N^3) to O(3N). The dGSEM possesses several desirable characteristics 
such as: high-order accuracy, geometrical flexibility compared to global spectral methods, high scalability due to 
the high arithmetic intensity per element, suitability for GPU acceleration, and support for both h- and p- refinement.

## Oct-tree cell-based AMR
Modeling of the atmosphere is challenging in that a range of spatial and temporal scales are involved [@wallace2006].
Adaptive Mesh Refinement (AMR) provides the tool to address multi-scale phenomenon efficiently by focusing resources
where they are needed. NebulaSEM implements oct-tree cell-based AMR using the forest-of-octrees approach pioneered in [@p4est].
The AmrIteration class provides a high-level interface to enable AMR for any solver written using NebulaSEM.
A single loop enclosing the timestep iterations and declaration of fields involved in the PDE is enough to provide AMR 
capability for any solver.

## Support for new solver development
NebulaSEM offers a range of operators for spatial and temporal discretization, streamlining the development of 
solvers for Partial Differential Equations (PDEs). The provided code snippet includes an example solver for 
the advection-diffusion equation.

```C++
void transport() {
    for (AmrIteration ait; !ait.end(); ait.next()) {     /*AMR iteration object (ait) and loop*/
         VectorCellField U("U", READWRITE);               /*Velocity field defined over the grid*/
         ScalarCellField T("T", READWRITE);               /*Scalar field*/
         ScalarFacetField F = flx(U);                     /*Compute flux field*/
         ScalarCellField mu = 1;                          /*Diffusion parameter*/
         for (Iteration it(ait.get_step()); !it.end(); it.next()) { /*Time loop with support for deferred correction */
             ScalarCellMatrix M;                          /*Matrix for the PDE discretization*/
             M = div(T,U,F,&mu) - lap(T,mu);              /*Divergence & Laplacian terms*/
             addTemporal<1>(M);                           /*Add temporal derivative*/
             Solve(M);                                    /*Solve the matrix */
         }
    }
}
```

Spatial operators within dGSEM encompass divergence, gradient, laplacian, and more. Temporal discretization is 
accomplished through explicit and implicit schemes, including first-order Euler explicit and implicit schemes, 
linear multi-step methods such as Adams-Moulton and Adams-Bashforth, the Runge-Kutta method up to 4th order, 
and fully-implicit Backward Differencing (BDF) methods.

## Turbulence models
The software includes turbulence models for high-Reynolds CFD simulations including a suite of Reynolds Averaged Navier Stokes (RANS) 
[@henk1972] and Large Eddy Simulation (LES) models [@les]. The list includes a mixing-length model, k-epsilon, k-omega, RNG k-epsilon, RNG k-omega
and the Smagornisky-Lilly LES model. Theese turbulence models has been utilized to evaluate the aerodynamic roughness of the 
built environment and complex terrain in [@abdi5].

## Parallelization with MPI+OpenMP/OpenACC
NebulaSEM achieves scalability on supercomputers through a combination of coarse- and fine-grained parallelism. 
The Message Passing Interface (MPI) is employed for distributed computing, while directive-based threading libraries 
such as OpenMP for CPUs and OpenACC for GPUs optimize fine-grained parallelism, minimizing communication overhead.


# Showcases
We demonstrate the capabilities of NebulaSEM using two example applications: a) as a generic CFD application
depicted in \autoref{fig:pitz-daily} b) as a non-hydrostatic dynamical core for atmospheric simulation \autoref{fig:srtb}.

![Simulation results for the Pitz-Daily problem [@pitz1983] that aims to evaluate the effect of combustion
on mixing layer growth. A snapshot of large eddy simulation (LES) results using finite-volume method of NebulaSEM is presented.
\label{fig:pitz-daily}](pitz-daily.png)

![Simulation results for the Robert Rising Thermal Bubble (RRTB) problem [@robert1993] using oct-tree cell-based AMR
with two levels of refinement. The discontinuous Galerkin spectral element method is used with polynomial order of 4.
\label{fig:srtb}](srtb-amr.png){width=50%, height=50%}

# Statement of need

Over the years, NebulaSEM has undergone a transformation from being a purely Computational Fluid Dynamics (CFD) application 
to primarily serve as an atmospheric simulation code. Consequently, it addresses specific needs inherent to both types of
applications, as outlined below.

Many CFD codes predominantly employ first or second-order accurate finite volume methods. NebulaSEM, on the other hand,
offers high-order discretization characterized by exceptional accuracy and minimal dissipation. This feature proves 
invaluable for tasks such as accurately capturing shocks and discontinuities, conducting highly accurate large eddy 
simulations, simulating turbulent flows with precision, and facilitating high-fidelity aeroacoustic simulations.

While high-order methods are more commonly associated with atmospheric modeling, such endeavors often rely on either 
finite-difference discretization, which is less geometrically flexible, or global spectral methods on latitude-longitude grids, 
which lack scalability. NebulaSEM, utilizing the discontinuous Galerkin spectral element method (dGSEM) for discretization, 
distinguishes itself by providing both geometric flexibility and high scalability. Moreover, NebulaSEM incorporates 
dynamic adaptive mesh refinement capabilities, a feature not commonly found in traditional atmospheric modeling approaches.


# Acknowledgements

I would like to thank my postdoc supervisor Francis X. Giraldo who introduced me to the high-order discontinuous Galerkin method.

# References
