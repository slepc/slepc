Installation
============

Install from PyPI using **pip**
-------------------------------

You can use :program:`pip` to install slepc4py_ and its
dependencies.

If you have a working MPI implementation and the ``mpicc`` compiler
wrapper is on your search path, it is highly recommended to install
mpi4py_ first::

  $ python -m pip install mpi4py

Ensure you have NumPy_ and petsc4py_ installed::

  $ python -m pip install numpy petsc petsc4py

Finally, install slepc4py_::

  $ python -m pip install slepc slepc4py

If you already have working PETSc and SLEPc installs, set environment
variables :envvar:`SLEPC_DIR` and :envvar:`PETSC_DIR` (and perhaps
:envvar:`PETSC_ARCH` for non-prefix installs) to appropriate values
and next use :program:`pip`::

  $ export SLEPC_DIR=/path/to/slepc
  $ export PETSC_DIR=/path/to/petsc
  $ export PETSC_ARCH=arch-linux2-c-opt
  $ python -m pip install petsc4py slepc4py

Install from the SLEPc source tree
----------------------------------

If you also want to install petsc4py_ from the PETSc source tree, follow
the instructions in the `petsc4py installation`_ page.

Set the :envvar:`PETSC_DIR` and :envvar:`PETSC_ARCH`
environment variables, as well as :envvar:`SLEPC_DIR`.
Follow the instructions to `build SLEPc`_. Then
:file:`cd` to the top of the SLEPc source tree and run::

  $ python -m pip install src/binding/slepc4py

The installation of slepc4py_ supports multiple :envvar:`PETSC_ARCH`
in the form of a colon separated list::

  $ PETSC_ARCH='arch-0:...:arch-N' python -m pip install src/binding/slepc4py

If you are cross-compiling, and the :mod:`numpy` module cannot be loaded on
your build host, then before invoking :program:`pip`, set the
:envvar:`NUMPY_INCLUDE` environment variable to the path that would be returned
by :samp:`import numpy; numpy.get_include()`::

  $ export NUMPY_INCLUDE=/usr/lib/pythonX/site-packages/numpy/core/include

.. _NumPy: https://www.numpy.org
.. _mpi4py: https://mpi4py.readthedocs.io
.. _petsc4py: https://petsc.org/release/petsc4py/reference.html
.. _petsc4py installation: https://petsc.org/release/petsc4py/install.html
.. _slepc4py: https://slepc.upv.es/release/slepc4py/reference.html
.. _build SLEPc: https://slepc.upv.es/release/installation/index.html
