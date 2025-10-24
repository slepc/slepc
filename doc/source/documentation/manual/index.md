# SLEPc Users Manual

```{only} html
**J. E. Roman** {sup}`1`,
**C. Campos** {sup}`1,2`,
**L. Dalcin** {sup}`3`,
**E. Romero** {sup}`1,4`,
**A. Tomas** {sup}`1`

{sup}`1` Universitat Polit&egrave;cnica de Val&egrave;ncia<br>
{sup}`2` Universitat de Val&egrave;ncia<br>
{sup}`3` Extreme Computing Research Center, King Abdullah University of Science and Technology<br>
{sup}`4` Jefferson Lab
```

(abstract)=
**Abstract**

This manual describes SLEPc, the *Scalable Library for Eigenvalue Problem Computations*, a software package for the solution of large sparse eigenproblems on parallel computers. It can be used to solve various types of eigenvalue problems, including linear and nonlinear, as well as other related problems such as the singular value decomposition (see a summary of supported problem classes in table [](#tab:modules)). SLEPc is a general library in the sense that it covers both Hermitian and non-Hermitian problems, with either real or complex arithmetic.

The emphasis of the software is on methods and techniques appropriate for problems in which the associated matrices are large and sparse, for example, those arising after the discretization of partial differential equations. Thus, most of the methods offered by the library are projection methods, including different variants of Krylov and Davidson iterations. In addition to its own solvers, SLEPc provides transparent access to some external software packages such as ARPACK. Apart from the solvers, SLEPc also provides built-in support for some operations commonly used in the context of eigenvalue computations, such as preconditioning or the shift-and-invert spectral transformation.

SLEPc is built on top of PETSc, the Portable, Extensible Toolkit for Scientific Computation. SLEPc extends PETSc with all the functionality necessary for the solution of eigenvalue problems. This means that PETSc must be previously installed in order to use SLEPc. PETSc users will find SLEPc very easy to use, since it enforces the same programming paradigm. Those readers that are not acquainted with PETSc are highly recommended to familiarize with it before proceeding with SLEPc.

{#how-to-get-slepc}
```{only} latex
% Only for latex as it would be duplicated on the web

**How to Get SLEPc**

All the information related to SLEPc can be found at the following web site:

> ::: center
> <https://slepc.upv.es>.
> :::

The distribution file is available for download at this site. Instructions for installing the software can be found in section [](#sec:inst).

PETSc can be downloaded from <https://petsc.org>. PETSc is supported, and information on contacting support can be found at that site.
```

{#additional-documentation .unnumbered}
**Additional Documentation**

```{only} latex
This manual provides a general description of SLEPc. In addition, manual pages for individual routines are available on-line at <https://slepc.upv.es> (both the C/Fortran API and the Python API). These manual pages provide hyperlinked access to the source code and enable easy movement among related topics. Finally, there are also several hands-on exercises available, which are intended for learning the basic concepts easily, see <https://slepc.upv.es/documentation>.
```
```{only} html
This manual provides a general description of SLEPc. In addition, manual pages for individual routines are available through the top bar menu (both the [C/Fortran API](../../manualpages/index) and the [Python API](../../slepc4py/index)). These manual pages provide hyperlinked access to the source code and enable easy movement among related topics. Finally, there are also several [hands-on exercises available](../../documentation/hands-on/index), which are intended for learning the basic concepts easily.
```

{#how-to-read-this-manual .unnumbered}
**How to Read this Manual**

Users that are already familiar with PETSc can read chapter [](#ch:int) very fast. Section [](#sec:eig) provides a brief overview of eigenproblems and the general concepts used by eigensolvers, so it can be skipped by experienced users. Chapters [](#ch:eps) to [](#ch:aux) describe the main SLEPc functionality. Some of them include an advanced usage section that can be skipped at a first reading. Finally, chapter [](#ch:add) contains less important, additional information.

{#slepc-technical-reports .unnumbered}
**SLEPc Technical Reports**

```{only} latex
The information contained in this manual is complemented by a set of Technical Reports available at <https://slepc.upv.es/documentation>. They provide technical details that normal users typically do not need to know but may be useful for experts in order to identify the particular method implemented in SLEPc.
```
```{only} html
The information contained in this manual is complemented by a set of [Technical Reports](../../documentation/index). They provide technical details that normal users typically do not need to know but may be useful for experts in order to identify the particular method implemented in SLEPc.
```

{#acknowledgments .unnumbered}
```{only} latex
**Acknowledgments**

The current version contains code contributed by several people. A list is given at <https://slepc.upv.es/contact>. That page also lists the funding sources that have made the development of SLEPc possible. We are thankful to all of them.
```

{#license-and-copyright}
```{only} latex
**License and Copyright**

Starting from version 3.8, SLEPc is released under a 2-clause BSD license (see `LICENSE` file).

> ::: sffamily
> Copyright 2002--{{release_year}} Universitat Polit&egrave;cnica de Val&egrave;ncia, Spain
> :::
```

{#supported-problem-classes .unnumbered}
**Supported Problem Classes**

The following table provides an overview of the functionality offered by SLEPc, organized by problem classes.

:::{table} SLEPc modules
:name: tab:modules

  Problem class                  |               Model equation               | Module | Chapter
  -------------------------------|--------------------------------------------|--------|-------------------------------------------------------------------
  Linear eigenvalue problem      |     $Ax=\lambda x,\quad Ax=\lambda Bx$     | `EPS`  | [](#ch:eps)
  Polynomial eigenvalue problem  | $(A_0+\lambda A_1+\cdots+\lambda^dA_d)x=0$ | `PEP`  | [](#ch:pep)
  Nonlinear eigenvalue problem   |              $T(\lambda)x=0$               | `NEP`  | [](#ch:nep)
  Singular value decomposition   |               $Av=\sigma u$                | `SVD`  | [](#ch:svd)
  Matrix function (action of)    |                 $y=f(A)v$                  | `MFN`  | [](#ch:mfn)
  Linear matrix equation         |                $AXE+DXB=C$                 | `LME`  | [](#ch:lme)
:::

In order to solve a given problem, one should create a solver object corresponding to the solver class (module) that better fits the problem (the less general one; e.g., we do not recommend using `NEP` to solve a linear eigenproblem).

(notes)=
:::{admonition} Notes

-   Most users are typically interested in linear eigenproblems only.

-   In each problem class there may exist several subclasses (problem types in SLEPc terminology), for instance symmetric-definite generalized eigenproblem in `EPS`.

-   The solver class (module) is named after the problem class. For historical reasons, the one for linear eigenvalue problems is called `EPS` rather than `LEP`.

-   In addition to the SVD shown in the table, the `SVD` module also supports other related problems such as the GSVD and the HSVD.

-   For the action of a matrix function (`MFN`), in SLEPc we focus on methods that are closely related to methods for eigenvalue problems.
:::

```{raw} latex
\sphinxtableofcontents
```

```{toctree}
:maxdepth: 3
:hidden:

intro
eps
st
svd
pep
nep
mfn
lme
aux
extra
```
