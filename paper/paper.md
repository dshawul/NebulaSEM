---
title: 'ERF: Energy Research and Forecasting'

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

NebulaSEM is a model that employs the high-order discontinuous Galerkin spectral element method for solving
the governing equations of fluid dynamics. It was originally designed as a generic Computational Fluid Dynamics (CFD)
application using an unstructured grid with polyhedral elements. This was first used for the author's research in evaluation of
aerodynamic roughness of the built environment and complex terrain [@abdi5]. Later, it is extended with
dGSEM spatial discretization that has several desirable characteristics such as: high-order accuracy, geometrical
flexibility compared to global spectral methods, high scalability due to the high arithimetic intensity per element,
suitability for GPU acceleration, and support for both h- and p- refinement.

# NebulaSEM Features

## High-order dGSEM discretization
## Oct-tree cell-based AMR
## Eases writing of new solvers for PDEs
## A suite of RANS and LES turbulence models
## Parallelization with MPI+OpenMP/OpenACC

# Statement of need

NebulaSEM has evolved over the years from a purely CFD application to primarily an atmospheric simulation code.
Many CFD codes are mostly written with first or second-order accurate finite volume method. NebulaSEM provides
high-order discretization that is minimally dissipative. This can be useful for example, capturing shock/discontinuity
accurately, highly accurate large eddy simulations and turbulent flow simulations in general, high-fidelity aeroacoustic
simulations etc.

Although it is more common to see high-order methods in the area of atmospheric modelling, it is often carried out
using either finite-difference discretiation that is less geometrically flexible, or using globabl spectral methods on
latitude-longitude grids that is not scalable. The dGSEM discretization in NebulaSEM provides both geometrical flexibility
and high scalability. Additionally, NebulaSEM provides dynamic adaptive mesh refinement capabilities that is not prevalent
in atmospheric modelling. 

# Acknowledgements

I would like to thank my postdoc supervisor Francis X. Giraldo who introduced me to the high-order discontinuos Galerkin method.

# References
