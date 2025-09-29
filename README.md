
SLEPc: Scalable Library for Eigenvalue Problem Computations
===========================================================

- Authors: Jose E. Roman, Carmen Campos, Lisandro Dalcin, Eloy Romero, Andres Tomas
- Organization: Universitat Politecnica de Valencia, Spain
- Website: https://slepc.upv.es
- Contact: slepc-maint@upv.es


Overview
--------

SLEPc, the Scalable Library for Eigenvalue Problem Computations, is a software package for the solution of *large sparse eigenvalue problems* on parallel computers. It can be used for the solution of problems formulated in either standard or generalized form, as well as other related problems such as the singular value decomposition and the nonlinear eigenproblem.

The emphasis of the software is on methods and techniques appropriate for problems in which the associated matrices are sparse. Therefore, most of the methods offered by the library are projection methods such as Krylov-Schur or Jacobi-Davidson. It also provides built-in support for spectral transformation such as shift-and-invert.

SLEPc is built on top of [PETSc](https://petsc.org/), the Portable Extensible Toolkit for Scientific Computation. It can be considered an extension of PETSc providing all the functionality necessary for the solution of eigenvalue problems.

### Supported problem classes

| Problem class                 | Model equation                               | Module |
| ----------------------------- | -------------------------------------------- | ------ |
| Linear eigenvalue problem     | $`Ax=\lambda x,\quad Ax=\lambda Bx`$         | `EPS`  |
| Quadratic eigenvalue problem  | $`(K+\lambda C+\lambda^2M)x=0`$              | `PEP`  |
| Polynomial eigenvalue problem | $`(A_0+\lambda A_1+\cdots+\lambda^dA_d)x=0`$ | `PEP`  |
| Nonlinear eigenvalue problem  | $`T(\lambda)x=0`$                            | `NEP`  |
| Singular value decomposition  | $`Av=\sigma u`$                              | `SVD`  |
| Matrix function (action of)   | $`y=f(A)v`$                                  | `MFN`  |
| Linear matrix equation        | $`AXE+DXB=C`$                                | `LME`  |


Documentation
-------------

The Users Manual as well as the HTML man pages for the detailed reference of each individual SLEPc routines are included in the SLEPc distribution and can also be found at the [online documentation](https://slepc.upv.es/release/documentation).

The main reference for SLEPc is the following paper (see [other references at the SLEPc website](https://slepc.upv.es/release/material)):

- V. Hernandez, J. E. Roman, and V. Vidal, *SLEPc: A scalable and flexible toolkit for the solution of eigenvalue problems*, ACM Trans. Math. Software 31: 351-362 (2005). [DOI](https://doi.org/10.1145%2F1089014.1089019)


Installation
------------

The installation procedure of SLEPc is very similar to that of PETSc. Briefly, the environment variables `$SLEPC_DIR` and `$PETSC_DIR` must be set, then the `configure` script is executed and finally the libraries are built with the command `make`. More details can be found in the Users Manual or in the online [installation instructions](https://slepc.upv.es/release/installation).


License
-------

SLEPc is distributed under a 2-clause BSD license (see [LICENSE.md](LICENSE.md)). Note that this license DOES NOT apply to any software that may be obtained via the `--download-package` option of the SLEPc configuration. Each of those packages are covered by their own licenses.

Funding
-------

The development of SLEPc has been partially supported by the following grants:

- Oficina de Ciencia i Tecnologia, Generalitat Valenciana, CTIDB/2002/54.
- D. G. Investigacio i Transf. de Tecnologia, Generalitat Valenciana, GV06/091.
- Ministerio de Ciencia e Innovacion, TIN2009-07519.
- Ministerio de Economia y Competitividad, TIN2013-41049-P.
- Agencia Estatal de Investigacion, TIN2016-75985-P.
- Agencia Estatal de Investigacion, PID2019-107379RB-I00.
- Agencia Estatal de Investigacion, PID2022-139568NB-I00.

