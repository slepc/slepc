Release 3.24
------------

- Update to SLEPc 3.24.

- Support (opt-in via setting the environment variable
  ``SLEPC4PY_BUILD_PYSABI=1``) for building with ``Py_LIMITED_API``
  (Python Stable ABI) under Python 3.10+ (requires Cython 3.1+).

- Add support for standard Python operators for `BV` and `FN` classes.

- Add new `LME` class.


Release 3.23
------------

- Update to SLEPc 3.23.


Release 3.22
------------

- Update to SLEPc 3.22.
- In slepc4py now `EPS.getEigenpair()` and `EPS.getEigenvalue()` will return a real value
  instead of a complex, if the problem is of Hermitian or generalized Hermitian type.


Release 3.21
------------

- Update to SLEPc 3.21.


Release 3.20
------------

- Update to SLEPc 3.20.


Release 3.19
------------

- Update to SLEPc 3.19.


Release 3.18
------------

- Update to SLEPc 3.18.


Release 3.17
------------

- Update to SLEPc 3.17.


Release 3.16
------------

- Update to SLEPc 3.16.


Release 3.15
------------

- Update to SLEPc 3.15.
- Updates in installation scripts.


Release 3.14
------------

- Update to SLEPc 3.14.


Release 3.13
------------

- Update to SLEPc 3.13.


Release 3.12
------------

- Update to SLEPc 3.12.


Release 3.11
------------

- Update to SLEPc 3.11.


Release 3.10
------------

- Update to SLEPc 3.10.


Release 3.9
-----------

- Update to SLEPc 3.9.


Release 3.8
-----------

- Update to SLEPc 3.8.


Release 3.7
-----------

- Update to SLEPc 3.7.


Release 3.6
-----------

- Update to SLEPc 3.6.


Release 3.5
-----------

- Update to SLEPc 3.5.
- Add RG class introduced in SLEPc 3.5 release.
- Add PySlepcXXX_New/Get C API functions.
- Fix compilation problem with complex scalars on OS X.
- Fix outdated SWIG interface file.


Release 3.4
-----------

- Update to SLEPc 3.4.


Release 3.3
-----------

- Update to SLEPc 3.3.
- Regenerate the wrappers using Cython 0.18 and fix binary
  compatibility issues with petsc4py 3.3.1.


Release 1.2
-----------

- Update to SLEPc 3.2.


Release 1.1
-----------

- Support for new QEP quadratic eigenproblem solver in SLEPc.
- Support for ``pip install slepc4py`` to download and install SLEPc.
- Support for PETSc/SLEPc static library builds (Linux-only).
- Preliminary support for Python 3.


Release 1.0.0
-------------

- This is the fist release of the all-new, Cython-based, implementation
  of *SLEPc for Python*.
