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

Discussion and Support
======================

+ You can write to the `PETSc Users Mailing List
  <https://petsc.org/release/community/mailing/#mailing-lists>`__
+ See also the `contact information for SLEPc
  <https://slepc.upv.es/release/contact>`__
+ Issue Tracker:  https://gitlab.com/slepc/slepc/-/issues
+ Git Repository: https://gitlab.com/slepc/slepc.git
  (the source code is in ``src/binding/slepc4py``)
+ Current release in PyPI: https://pypi.org/project/slepc4py/

Acknowledgments
===============

L. Dalcin was partially supported by the
Extreme Computing Research Center (ECRC),
Division of Computer, Electrical, and
Mathematical Sciences & Engineering (CEMSE),
King Abdullah University of Science and Technology (KAUST).

See additional acknowledgements related to `the SLEPc project
<https://slepc.upv.es/release/contact/acknowledgement.html>`__.

Contents
========

.. include:: links.txt

.. toctree::
   :maxdepth: 2

   overview
   tutorial
   install
   citing
   changes

.. toctree::
   :caption: Python specifics
   :maxdepth: 2

   reference
   demo/demo
