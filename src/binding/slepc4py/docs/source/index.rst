================
SLEPc for Python
================

.. only:: html or man

   :Authors:      Lisandro Dalcin, Jose E. Roman
   :Contact:      dalcinl@gmail.com, jroman@dsic.upv.es
   :Web Site:     https://gitlab.com/slepc/slepc
   :Date:         |today|

.. topic:: Abstract

   This document describes slepc4py_, a Python_ wrapper to the SLEPc_
   libraries.

   SLEPc_ is a software package for the parallel solution of
   large-scale eigenvalue problems. It can be used for computing
   eigenvalues and eigenvectors of large, sparse matrices, or matrix
   pairs, and also for computing singular values and vectors of a
   rectangular matrix. It also provides more advanced functionality
   such as solvers for nonlinear eigenvalue problems (either polynomial
   or general) and matrix functions.

   SLEPc_ relies on PETSc_ for basic functionality such as the
   representation of matrices and vectors, and the solution of linear
   systems of equations. Thus, slepc4py_ must be used together with
   its companion petsc4py_.

Contents
========

.. include:: links.txt

.. toctree::
   :maxdepth: 2

   overview
   tutorial
   install
   citing

.. toctree::
   :caption: Python specifics
   :maxdepth: 2

   reference
   demo/demo
