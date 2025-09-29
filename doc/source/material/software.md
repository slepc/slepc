# Software Packages that Use SLEPc

The list below shows software packages that use SLEPc or expose part of its functionality in a higher-level context.

## Finite elements and high-level computational toolkits

  1. [FEniCS](https://fenicsproject.org): a toolkit for the Automation of Computational Mathematical Modeling (ACMM).
  2. [RBniCS](https://www.rbnicsproject.org/): reduced order modelling in FEniCS.
  3. [dolfin-adjoint](http://www.dolfin-adjoint.org): automatic computation of adjoint linear models using DOLFIN (FEniCS). See also [tlm_adjoint](https://github.com/jrmaddison/tlm_adjoint).
  4. [Firedrake](https://www.firedrakeproject.org/): an automated system for the solution of PDEs with FEM. See [here](https://www.firedrakeproject.org/demos/qgbasinmodes.py) an example that uses SLEPc.
  5. [libMesh](https://libmesh.github.io): a C++ framework for the numerical simulation of partial differential equations.
  6. [deal.II](https://dealii.org/): a finite element Differential Equations Analysis Library. See [here](https://dealii.org/9.0.0/doxygen/deal.II/step_36) an example that uses SLEPc.
  7. [GetDP](https://www.getdp.info): a General Environment for the Treatment of Discrete Problems. See [here](https://gitlab.onelab.info/doc/models/-/wikis/Bloch-modes-in-periodic-waveguides) an example that computes eigenvalues.
  8. [MOOSE](https://mooseframework.inl.gov/): Multiphysics Object Oriented Simulation Environment.
  9. [PHAML](https://www.nist.gov/programs-projects/parallel-hierarchical-adaptive-multilevel-project-phaml): adaptive finite elements for elliptic PDEs.
  10. [FreeFEM](https://freefem.org): A high level multiphysics finite element software.
  11. [MFEM](https://mfem.org): a free, lightweight, scalable C++ library for finite element methods.
  12. [OOFEM](http://www.oofem.org): an Object Oriented Finite Element code.
  13. [PHG](https://lsec.cc.ac.cn/phg/index_en.htm): Parallel Hierarchical Grid, an adaptive mesh refinement FEM framework.
  14. [Feel++](https://docs.feelpp.org/home/index.html): a C++ library for partial differential equation solves using generalized Galerkin methods.
  15. [FEMuS](https://github.com/eaulisa/MyFEMuS): open-source Finite Element C++ library.
  16. [OpenCMISS](https://opencmiss.org/): Open Continuum Mechanics, Imaging, Signal processing and System identification.
  17. [SfePy](https://sfepy.org/doc-devel/): Simple Finite Elements in Python.
  18. [ff-bifbox](https://github.com/cmdoug/ff-bifbox?tab=readme-ov-file): FreeFEM scripts for numerical continuation, bifurcation analysis, resolvent analysis, and time-integration of large-scale time-dependent nonlinear PDEs on adaptively refined meshes.

## Many-body calculations, quantum systems, photonics

  1. [ELSI](https://wordpress.elsi-interchange.org/): ELectronic Structure Infrastructure.
  2. [DFT-FE](https://sites.google.com/umich.edu/dftfe): real-space DFT calculations using Finite Elements.
  3. [TiberCAD](http://www.tibercad.org): multiscale device simulator.
  4. [NEMO5](https://engineering.purdue.edu/gekcogrp/software-projects/nemo5): NanoElectronics MOdeling Tools, which is the basis of other tools such as [Quantum Dot Lab](https://nanohub.org/tools/qdot).
  5. [Femwell](https://helgegehring.github.io/femwell/): simulation tool for integrated circuits, electric and photonic.
  6. [Hammer](http://www.thphys.nuim.ie/hammer): numerical tools for treating systems of strongly interacting quantum many body systems.
  7. [Yambo](https://www.yambo-code.org/): many-body calculations in solid state and molecular physics.
  8. [pyCTQW](https://pyctqw.readthedocs.io): Continuous-Time Quantum Walk simulator.
  9. [PsiQuaSP](https://github.com/modmido/psiquasp): Permutation symmetry for identical Quantum Systems Package.
  10. [dynamite](https://dynamite.readthedocs.io): fast full quantum dynamics.
  11. [quimb](https://quimb.readthedocs.io): python library for quantum information and many-body calculations.
  12. [DanceQ](https://DanceQ.gitlab.io/danceq): Divide-And-conquer Number Conserving Exact diagonalization for Quantum systems.

## Plasma physics, nuclear engineering

  1. [GENE](https://www.genecode.org): Gyrokinetic Electromagnetic Numerical Experiment.
  2. [GYRO](https://gafusion.github.io/doc): the General Atomics TGYRO code suite.
  3. [PB3D](https://pb3d.github.io): Peeling-Ballooning in 3-D.
  4. [VERA](https://vera.ornl.gov): Virtual Environment for Reactor Applications.
  5. [BOUT++](https://boutproject.github.io): Plasma simulation in curvilinear coordinate systems.
  6. [FEMFFUSION](https://femffusion.webs.upv.es/): a finite element method code for nuclear reactor modelling.
  7. [Milonga](https://www.seamplex.com/milonga): a free nuclear reactor core analysis code.

## Other

  1. [SALSA](https://icl.utk.edu/salsa/): Self-Adapting Large-scale Solver Architecture.
  2. [Cubica](http://www.tkim.graphics/cubica/): a toolkit for subspace deformations.
  3. [Dome](http://faraday1.ucd.ie/dome.html): a power system analysis toolbox.
  4. [ncpaprop](https://github.com/chetzer-ncpa/ncpaprop-release): NCPA Infrasound Propagation Modeling Package.
  5. [EasterEig](https://github.com/nennigb/EasterEig): parametric eigenvalue problem depending on a parameter.
  6. [pyGPCCA](https://github.com/msmdev/pyGPCCA): Generalized Perron Cluster Cluster Analysis.
  7. [cmdtools](https://github.com/zib-cmd/cmdtools): a suite of tools used and/or developed in the Computational Molecular Design group of ZIB.
  8. [BROADCAST](https://broadcast.readthedocs.io): a Python software that discretizes the compressible Navier-Stokes equations and extracts the linearized state derivative operators.
  9. [VFI-MEMS](https://gitlab.tuwien.ac.at/andre.gesing/non_linear_eigenvalue): solve the dynamics of a MEMS resonator in a fluid with a nonlinear eigensolver.
  10. [helmholtz-x](https://doi.org/10.17863/CAM.112694): a python library built upon DOLFINx to solve a non-homogeneous Helmholtz equation, specifically thermoacoustic Helmholtz.
  11. [RSVD-{math}`\Delta`t](https://github.com/AliFarghadan/RSVD-Delta-t/tree/Resolvent-analysis): Randomized Singular Value Decomposition with Time-stepping for large-scale resolvent analysis
  12. [OpenParEM](https://openparem.org): Open-source parallel electromagnetic simulators
