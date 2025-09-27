# User manual

(abstract)=
**Abstract**

This document describes SLEPc, the *Scalable Library for Eigenvalue Problem Computations*, a software package for the solution of large sparse eigenproblems on parallel computers. It can be used for the solution of various types of eigenvalue problems, including linear and nonlinear, as well as other related problems such as the singular value decomposition (see a summary of supported problem classes in table [](#tab:modules). SLEPc is a general library in the sense that it covers both Hermitian and non-Hermitian problems, with either real or complex arithmetic.

The emphasis of the software is on methods and techniques appropriate for problems in which the associated matrices are large and sparse, for example, those arising after the discretization of partial differential equations. Thus, most of the methods offered by the library are projection methods, including different variants of Krylov and Davidson iterations. In addition to its own solvers, SLEPc provides transparent access to some external software packages such as ARPACK. These packages are optional and their installation is not required to use SLEPc, see section [](#sec:wrap) for details. Apart from the solvers, SLEPc also provides built-in support for some operations commonly used in the context of eigenvalue computations, such as preconditioning or the shift-and-invert spectral transformation.

SLEPc is built on top of PETSc, the Portable, Extensible Toolkit for Scientific Computation {cite:p}`Balay:PUM`. It can be considered an extension of PETSc providing all the functionality necessary for the solution of eigenvalue problems. This means that PETSc must be previously installed in order to use SLEPc. PETSc users will find SLEPc very easy to use, since it enforces the same programming paradigm. Those readers that are not acquainted with PETSc are highly recommended to familiarize with it before proceeding with SLEPc.

{#how-to-get-slepc}
```{only} latex
% Only for latex as it would be duplicated on the web

**How to Get SLEPc**

All the information related to SLEPc can be found at the following web site:

> ::: center
> <https://slepc.upv.es>.
> :::

The distribution file is available for download at this site. Other information is provided there, such as installation instructions and contact information. Instructions for installing the software can also be found in section [](#sec:inst).

PETSc can be downloaded from <https://petsc.org>. PETSc is supported, and information on contacting support can be found at that site.
```

{#additional-documentation .unnumbered}
**Additional Documentation**

This manual provides a general description of SLEPc. In addition, manual pages for individual routines are included in the distribution file in hypertext format, and are also available on-line at <https://slepc.upv.es/documentation>. These manual pages provide hyperlinked access to the source code and enable easy movement among related topics. Finally, there are also several hands-on exercises available, which are intended for learning the basic concepts easily.

{#how-to-read-this-manual .unnumbered}
**How to Read this Manual**

Users that are already familiar with PETSc can read chapter [](#ch:int) very fast. Section [](#sec:eig) provides a brief overview of eigenproblems and the general concepts used by eigensolvers, so it can be skipped by experienced users. Chapters [](#ch:eps)--[](#ch:mfn) describe the main SLEPc functionality. Some of them include an advanced usage section that can be skipped at a first reading. Finally, chapter [](#ch:add) contains less important, additional information.

{#slepc-technical-reports .unnumbered}
**SLEPc Technical Reports**

The information contained in this manual is complemented by a [set of Technical Reports](#str), which provide technical details that normal users typically do not need to know but may be useful for experts in order to identify the particular method implemented in SLEPc. These reports are not included in the SLEPc distribution file but can be accessed via the SLEPc web site.

{#acknowledgments .unnumbered}
**Acknowledgments**

The current version contains code contributed by: A. Lamas Davi&ntilde;a (CUDA code), F. Alvarruiz (restarted Lanczos for the GSVD, structured BSE solvers), B. Mellado-Pinto (structured BSE solvers), Y. Maeda, T. Sakurai (CISS solvers), M. Moldaschl, W. Gansterer (BDC subroutines), F. Kong (nonlinear inverse iteration), H. Fang, Y. Saad (`FILTLAN` polynomial filter).

Development of SLEPc has been partially funded by the following grants:

-   Innovation Study ISOLV-BSE has received funding through the Inno4scale project, which is funded by the European High-Performance Computing Joint Undertaking (JU) under Grant Agreement No 101118139. The JU receives support from the European Union's Horizon Europe Programme.

-   Agencia Estatal de Investigaci&oacute;n (Spain), grant no. PID2022-139568NB-I00, PI: Jos&eacute; E. Rom&aacute;n.

-   Agencia Estatal de Investigaci&oacute;n (Spain), grant no. PID2019-107379RB-I00, PI: Jos&eacute; E. Rom&aacute;n.

-   Agencia Estatal de Investigaci&oacute;n (Spain), grant no. TIN2016-75985-P, PI: Jos&eacute; E. Rom&aacute;n.

-   Ministerio de Econom&imath;&#769;a y Comp. (Spain), grant no. TIN2013-41049-P, PI: Jos&eacute; E. Rom&aacute;n.

-   Ministerio de Ciencia e Innovaci&oacute;n (Spain), grant no. TIN2009-07519, PI: Jos&eacute; E. Rom&aacute;n.

-   Valencian Regional Government, grant no. GV06/091, PI: Jos&eacute; E. Rom&aacute;n.

-   Valencian Regional Government, grant no. CTIDB/2002/54, PI: Vicente Hern&aacute;ndez.

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
  Quadratic eigenvalue problem   |       $(K+\lambda C+\lambda^2M)x=0$        |   --   | --
  Polynomial eigenvalue problem  | $(A_0+\lambda A_1+\cdots+\lambda^dA_d)x=0$ | `PEP`  | [](#ch:pep)
  Nonlinear eigenvalue problem   |              $T(\lambda)x=0$               | `NEP`  | [](#ch:nep)
  Singular value decomposition   |               $Av=\sigma u$                | `SVD`  | [](#ch:svd)
  Matrix function (action of)    |                 $y=f(A)v$                  | `MFN`  | [](#ch:mfn)
  Linear matrix equation         |                $AXE+DXB=C$                 | `LME`  | See [Notes](#notes)
:::

In order to solve a given problem, one should create a solver object corresponding to the solver class (module) that better fits the problem (the less general one; e.g., we do not recommend using `NEP` to solve a linear eigenproblem).

(notes)=
:::{admonition} Notes

-   Most users are typically interested in linear eigenproblems only.

-   In each problem class there may exist several subclasses (problem types in SLEPc terminology), for instance symmetric-definite generalized eigenproblem in `EPS`.

-   The solver class (module) is named after the problem class. For historical reasons, the one for linear eigenvalue problems is called `EPS` rather than `LEP`.

-   In addition to the SVD shown in the table, the `SVD` module also supports other related problems such as the GSVD and the HSVD.

-   In previous SLEPc versions there was a `QEP` module for quadratic eigenproblems. It has been replaced by `PEP`.

-   For the action of a matrix function (`MFN`), in SLEPc we focus on methods that are closely related to methods for eigenvalue problems.

-   The solver class `LME` is still experimental and it is not covered in this manual yet.
:::

```{rubric} Footnotes
```

```{eval-rst}
.. bibliography::
   :filter: docname in docnames
```

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
extra
```
