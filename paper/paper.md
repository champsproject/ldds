---
title: 'LDDS: A python package for computing and visualizing Lagrangian Descriptors for Dynamical Systems'
date: 8 January 2020
bibliography: 
- paper.bib
- chapter1part1.bib
- ham_dyn.bib
- hyper.bib
- KAMChapter.bib
- lag_des.bib
- LDs.bib
output:
  pdf_document:
    fig_caption: yes
    fig_height: 3
  html_document:
    fig_caption: yes
    fig_height: 3
tags:
- Dynamical systems
- Lagrangian descriptors
affiliations:
- index: 1
  name: School of Mathematics, University of Bristol, Fry Building, Woodland Road,
    Bristol BS8 1UG, United Kingdom
---

## Statement of Need

Nonlinear dynamical systems are ubiquitous in natural and engineering sciences, for example, fluid mechanics, theoretical chemistry, ship dynamics, rigid body dynamics, atomic physics, solid mechanics, condensed matter physics, mathematical biology, oceanography, meteorology, celestial mechanics [see @wiggins1994normally for a list of collated references]. There has been many advances in understanding phenomena across these disciplines using the geometric viewpoint of the solutions and the underlying structures in the phase space. Chief among these phase space structures are the invariant manifolds that form a barrier between dynamically distinct solutions. In most nonlinear systems, the invariant manifolds are computed using numerical techniques that rely on some form of linearization around equilibrium points followed by continuation and globalization. However, these methods become computationally expensive and challenging when applied to the high-dimensional phase space of chemical reactions [rephrase the introduction from the CNSNS paper] or vector field defined using numerical simulation or experimental data. This points to the need for techniques that can be paired with trajectory calculations, without excessive computational overhead and at the same time can be visualized along with trajectory data. The Python package, `LDDS`, serves this need for analyzing deterministic, stochastic, high-dimensional nonlinear dynamical systems described either by an analytical vector field or a trajectory data obtained from numerical simulations or experiments.

## Summary

`LDDS` has the following features:




### Example systems {#examples}



## Visualization of Lagrangian Descriptors



## Relation to ongoing research projects

We are developing geometric methods of phase space transport in the context of chemical reaction dynamics that rely on identifying and computing the UPOs.


## Acknowledgements

We acknowledge the support of EPSRC Grant No. EP/P021123/1 and Office of Naval Research (Grant No. N00014-01-1-0769). 


## References

