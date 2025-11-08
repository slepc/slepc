SLEPc for Python
================

Python bindings for SLEPc_.

.. _SLEPc:    https://slepc.upv.es

Install
-------

If you have a working MPI implementation and the ``mpicc`` compiler
wrapper is on your search path, it is highly recommended to install
``mpi4py`` first::

  $ pip install mpi4py

Ensure you have NumPy installed::

  $ pip install numpy

and finally::

  $ pip install petsc petsc4py slepc slepc4py


Citations
---------

If SLEPc for Python has been significant to a project that leads to an
academic publication, please acknowledge that fact by citing the
project.

* L. Dalcin, P. Kler, R. Paz, and A. Cosimo,
  *Parallel Distributed Computing using Python*,
  Advances in Water Resources, 34(9):1124-1139, 2011.
  https://doi.org/10.1016/j.advwatres.2011.04.013

* V. Hernandez, J.E. Roman and V. Vidal.
  *SLEPc: A Scalable and Flexible Toolkit for the
  Solution of Eigenvalue Problems*,
  ACM Transactions on Mathematical Software, 31(3):351-362, 2005.
  https://doi.org/10.1145/1089014.1089019
