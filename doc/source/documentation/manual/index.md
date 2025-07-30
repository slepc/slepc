::: titlepage
:::

# User manual

{#abstract .unnumbered}
## Abstract

This document describes SLEPc, the *Scalable Library for Eigenvalue Problem Computations*, a software package for the solution of large sparse eigenproblems on parallel computers. It can be used for the solution of various types of eigenvalue problems, including linear and nonlinear, as well as other related problems such as the singular value decomposition (see a summary of supported problem classes in table [](#tab:modules). SLEPc is a general library in the sense that it covers both Hermitian and non-Hermitian problems, with either real or complex arithmetic.

The emphasis of the software is on methods and techniques appropriate for problems in which the associated matrices are large and sparse, for example, those arising after the discretization of partial differential equations. Thus, most of the methods offered by the library are projection methods, including different variants of Krylov and Davidson iterations. In addition to its own solvers, SLEPc provides transparent access to some external software packages such as ARPACK. These packages are optional and their installation is not required to use SLEPc, see section [](#sec:wrap) for details. Apart from the solvers, SLEPc also provides built-in support for some operations commonly used in the context of eigenvalue computations, such as preconditioning or the shift-and-invert spectral transformation.

SLEPc is built on top of PETSc, the Portable, Extensible Toolkit for Scientific Computation {cite:p}`Balay:PUM`. It can be considered an extension of PETSc providing all the functionality necessary for the solution of eigenvalue problems. This means that PETSc must be previously installed in order to use SLEPc. PETSc users will find SLEPc very easy to use, since it enforces the same programming paradigm. Those readers that are not acquainted with PETSc are highly recommended to familiarize with it before proceeding with SLEPc.

{#how-to-get-slepc .unnumbered}
## How to Get SLEPc

All the information related to SLEPc can be found at the following web site:

> ::: center
> <https://slepc.upv.es>.
> :::

The distribution file is available for download at this site. Other information is provided there, such as installation instructions and contact information. Instructions for installing the software can also be found in section [](#sec:inst).

PETSc can be downloaded from <https://petsc.org>. PETSc is supported, and information on contacting support can be found at that site.

{#additional-documentation .unnumbered}
## Additional Documentation

This manual provides a general description of SLEPc. In addition, manual pages for individual routines are included in the distribution file in hypertext format, and are also available on-line at <https://slepc.upv.es/documentation>. These manual pages provide hyperlinked access to the source code and enable easy movement among related topics. Finally, there are also several hands-on exercises available, which are intended for learning the basic concepts easily.

{#how-to-read-this-manual .unnumbered}
## How to Read this Manual

Users that are already familiar with PETSc can read chapter [](#cap:int) very fast. Section [](#sec:eig) provides a brief overview of eigenproblems and the general concepts used by eigensolvers, so it can be skipped by experienced users. Chapters [](#cap:eps)--[](#cap:mfn) describe the main SLEPc functionality. Some of them include an advanced usage section that can be skipped at a first reading. Finally, chapter [](#cap:add) contains less important, additional information.

{#slepc-technical-reports .unnumbered}
## SLEPc Technical Reports

The information contained in this manual is complemented by a set of Technical Reports, which provide technical details that normal users typically do not need to know but may be useful for experts in order to identify the particular method implemented in SLEPc. These reports are not included in the SLEPc distribution file but can be accessed via the SLEPc web site. A [list of available reports](#str) is included at the end of the Bibliography.

{#acknowledgments .unnumbered}
## Acknowledgments

The current version contains code contributed by: A. Lamas Davi&ntilde;a (CUDA code), F. Alvarruiz (restarted Lanczos for the GSVD, structured BSE solvers), B. Mellado-Pinto (structured BSE solvers), Y. Maeda, T. Sakurai (CISS solvers), M. Moldaschl, W. Gansterer (BDC subroutines), F. Kong (nonlinear inverse iteration), H. Fang, Y. Saad ([filtlan]{.smallcaps} polynomial filter).

Development of SLEPc has been partially funded by the following grants:

-   Innovation Study ISOLV-BSE has received funding through the Inno4scale project, which is funded by the European High-Performance Computing Joint Undertaking (JU) under Grant Agreement No 101118139. The JU receives support from the European Union's Horizon Europe Programme.

-   Agencia Estatal de Investigaci&oacute;n (Spain), grant no. PID2022-139568NB-I00, PI: Jos&eacute; E. Rom&aacute;n.

-   Agencia Estatal de Investigaci&oacute;n (Spain), grant no. PID2019-107379RB-I00, PI: Jos&eacute; E. Rom&aacute;n.

-   Agencia Estatal de Investigaci&oacute;n (Spain), grant no. TIN2016-75985-P, PI: Jos&eacute; E. Rom&aacute;n.

-   Ministerio de Econom&imath;&#769;a y Comp. (Spain), grant no. TIN2013-41049-P, PI: Jos&eacute; E. Rom&aacute;n.

-   Ministerio de Ciencia e Innovaci&oacute;n (Spain), grant no. TIN2009-07519, PI: Jos&eacute; E. Rom&aacute;n.

-   Valencian Regional Government, grant no. GV06/091, PI: Jos&eacute; E. Rom&aacute;n.

-   Valencian Regional Government, grant no. CTIDB/2002/54, PI: Vicente Hern&aacute;ndez.

{#license-and-copyright .unnumbered}
## License and Copyright

Starting from version 3.8, SLEPc is released under a 2-clause BSD license (see `LICENSE` file).

> ::: sffamily
> Copyright 2002--2025 Universitat Polit&egrave;cnica de Valencia, Spain
> :::

{#supported-problem-classes .unnumbered}
## Supported Problem Classes

The following table provides an overview of the functionality offered by SLEPc, organized by problem classes.

:::{table} SLEPc modules
:name: tab:modules

  Problem class                  |               Model equation               | Module | Chapter
  -------------------------------|--------------------------------------------|--------|-------------------------------------------------------------------
  Linear eigenvalue problem      |     $Ax=\lambda x,\quad Ax=\lambda Bx$     | `EPS`  | [](#cap:eps)
  Quadratic eigenvalue problem   |       $(K+\lambda C+\lambda^2M)x=0$        |   --   | --
  Polynomial eigenvalue problem  | $(A_0+\lambda A_1+\cdots+\lambda^dA_d)x=0$ | `PEP`  | [](#cap:pep)
  Nonlinear eigenvalue problem   |              $T(\lambda)x=0$               | `NEP`  | [](#cap:nep)
  Singular value decomposition   |               $Av=\sigma u$                | `SVD`  | [](#cap:svd)
  Matrix function (action of)    |                 $y=f(A)v$                  | `MFN`  | [](#cap:mfn)
  Linear matrix equation         |                $AXE+DXB=C$                 | `LME`  | See [](#notes)
:::

In order to solve a given problem, one should create a solver object corresponding to the solver class (module) that better fits the problem (the less general one; e.g., we do not recommend using `NEP` to solve a linear eigenproblem).\
Notes:

:::{admonition} Notes
:name: notes

-   Most users are typically interested in linear eigenproblems only.

-   In each problem class there may exist several subclasses (problem types in SLEPc terminology), for instance symmetric-definite generalized eigenproblem in `EPS`.

-   The solver class (module) is named after the problem class. For historical reasons, the one for linear eigenvalue problems is called `EPS` rather than `LEP`.

-   In addition to the SVD shown in the table, the `SVD` module also supports other related problems such as the GSVD and the HSVD.

-   In previous SLEPc versions there was a `QEP` module for quadratic eigenproblems. It has been replaced by `PEP`.

-   For the action of a matrix function (`MFN`), in SLEPc we focus on methods that are closely related to methods for eigenvalue problems.

-   The solver class `LME` is still experimental and it is not covered in this manual yet.
:::

```{toctree}
:maxdepth: 2
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

::: thebibliography
40 natexlab

Alvarruiz, F., C. Campos, and J. E. Roman (2024). Thick-restarted joint Lanczos bidiagonalization for the GSVD. , 440:115506.

Alvarruiz, F., B. Mellado-Pinto, and J. E. Roman (2025). Variants of thick-restart Lanczos for the Bethe-Salpeter eigenvalue problem. .

Anderson, E., Z. Bai, C. Bischof, L. S. Blackford, J. Demmel, J. Dongarra, J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney, and D. Sorensen (1999). . SIAM, Philadelphia, third edition.

Auckenthaler, T., V. Blum, H.-J. Bungartz, T. Huckle, R. Johanni, L. Kr&auml;mer, B. Lang, H. Lederer, and P. R. Willems (2011). Parallel solution of partial symmetric eigenvalue problems from electronic structure calculations. , 37:783--794.

Bai, Z., J. Demmel, J. Dongarra, A. Ruhe, and H. van der Vorst (eds.) (2000). . SIAM, Philadelphia.

Balay, S., S. Abhyankar, M. F. Adams, S. Benson, J. Brown, P. Brune, K. Buschelman, E. Constantinescu, L. Dalcin, A. Dener, V. Eijkhout, J. Faibussowitsch, W. D. Gropp, V. Hapla, T. Isaac, P. Jolivet, D. Karpeyev, D. Kaushik, M. G. Knepley, F. Kong, S. Kruger, D. A. May, L. C. McInnes, R. T. Mills, L. Mitchell, T. Munson, J. E. Roman, K. Rupp, P. Sanan, J. Sarich, B. F. Smith, H. Suh, S. Zampini, H. Zhang, H. Zhang, and J. Zhang (2025). Users Manual Revision 3.23. Argonne Technical Memorandum ANL-21/39.

Balay, S., W. D. Gropp, L. C. McInnes, and B. F. Smith (1997). Efficient management of parallelism in object oriented numerical software libraries. In *Modern Software Tools in Scientific Computing* (edited by E. Arge, A. M. Bruaset, and H. P. Langtangen), pp. 163--202. Birkha&uuml;ser.

Betcke, T. (2008). Optimal scaling of generalized and polynomial eigenvalue problems. , 30(4):1320--1338.

Betcke, T. and D. Kressner (2011). Perturbation, extraction and refinement of invariant pairs for matrix polynomials. , 435(3):514--536.

Bj&ouml;rck, &angst;. (1996). . SIAM, Philadelphia.

Blackford, L. S., J. Choi, A. Cleary, E. D'Azevedo, J. Demmel, I. Dhillon, J. Dongarra, S. Hammarling, G. Henry, A. Petitet, K. Stanley, D. Walker, and R. C. Whaley (1997). . SIAM, Philadelphia.

Campos, C. and J. E. Roman (2012). Strategies for spectrum slicing based on restarted Lanczos methods. , 60(2):279--295.

Campos, C. and J. E. Roman (2016a). Parallel iterative refinement in polynomial eigenvalue problems. , 23(4):730--745.

Campos, C. and J. E. Roman (2016b). Parallel Krylov solvers for the polynomial eigenvalue problem in SLEPc. , 38(5):S385--S411.

Campos, C. and J. E. Roman (2020a). A polynomial Jacobi--Davidson solver with support for non-monomial bases and deflation. , 60:295--318.

Campos, C. and J. E. Roman (2020b). Inertia-based spectrum slicing for symmetric quadratic eigenvalue problems. , 27(4):e2293.

Campos, C. and J. E. Roman (2021). NEP: a module for the parallel solution of nonlinear eigenvalue problems in SLEPc. , 47(3):23:1--23:29.

Canning, A., L. W. Wang, A. Williamson, and A. Zunger (2000). Parallel empirical pseudopotential electronic structure calculations for million atom systems. , 160(1):29--41.

Chen, T.-Y. and J. W. Demmel (2000). Balancing sparse matrices for computing eigenvalues. , 309(1--3):261--287.

Eiermann, M. and O. G. Ernst (2006). A restarted Krylov subspace method for the evaluation of matrix functions. , 44(6):2481--2504.

Ericsson, T. and A. Ruhe (1980). The spectral transformation Lanczos method for the numerical solution of large sparse generalized symmetric eigenvalue problems. , 35(152):1251--1268.

Fang, H. and Y. Saad (2012). A filtered Lanczos procedure for extreme and interior eigenvalue problems. , 34(4):A2220--A2246.

Golub, G. H. and H. A. van der Vorst (2000). Eigenvalue computation in the 20th century. , 123(1-2):35--65.

Golub, G. H. and C. F. van Loan (1996). . The Johns Hopkins University Press, Baltimore, MD, third edition.

Grimes, R. G., J. G. Lewis, and H. D. Simon (1994). A shifted block Lanczos algorithm for solving sparse symmetric generalized eigenproblems. , 15(1):228--272.

G&uuml;ttel, S. and F. Tisseur (2017). The nonlinear eigenvalue problem. , 26:1--94.

Hansen, P. C. (1998). . SIAM, Philadelphia, PA.

Higham, N. J. and A. H. Al-Mohy (2010). Computing matrix functions. , 19:159--208.

Knyazev, A. V., M. E. Argentati, I. Lashuk, and E. E. Ovtchinnikov (2007). lock Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX) in HYPRE and PETSc. , 29(5):2224--2239.

Lehoucq, R. B. and A. G. Salinger (2001). Large-scale eigenvalue calculations for stability analysis of steady flows on massively parallel computers. , 36:309--327.

Lehoucq, R. B., D. C. Sorensen, and C. Yang (1998). . SIAM, Philadelphia.

Li, R., Xi, Y., Erlandson, L., and Y. Saad (2019). The eigenvalues slicing library (EVSL): Algorithms, implementation, and software. , 41(4):C393-C415.

Maschhoff, K. J. and D. C. Sorensen (1996). : An efficient portable large scale eigenvalue package for distributed memory parallel architectures. , 1184:478--486.

Meerbergen, K. and A. Spence (1997). Implicitly restarted Arnoldi with purification for the shift-invert transformation. , 66(218):667--689.

Meerbergen, K., A. Spence, and D. Roose (1994). Shift-invert and Cayley transforms for detection of rightmost eigenvalues of nonsymmetric matrices. , 34(3):409--423.

Mehrmann, V. and H. Voss (2004). Nonlinear eigenvalue problems: a challenge for modern eigenvalue methods. , 27(2):121--152.

Morgan, R. B. and M. Zeng (2006). A harmonic restarted Arnoldi algorithm for calculating eigenvalues and determining multiplicity. , 415(1):96--113.

MPI Forum (1994). : a message-passing interface standard. , 8(3/4):159--416.

Nour-Omid, B., B. N. Parlett, T. Ericsson, and P. S. Jensen (1987). How to implement the spectral transformation. , 48(178):663--673.

Onn, R., A. O. Steinhardt, and A. Bojanczyk (1991). The hyperbolic singular value decomposition and applications. , 39(7):1575--1588.

Parlett, B. N. (1980). . Prentice-Hall, Englewood Cliffs, NJ. Reissued with revisions by SIAM, Philadelphia, 1998.

Polizzi, E. (2009). Density-matrix-based algorithm for solving eigenvalue problems. , 79(11):115112.

Poulson, J., B. Marker, R. A. van de Geijn, J. R. Hammond, and N. A. Romero (2013). Elemental: a new framework for distributed memory dense matrix computations. , 39(2):1--24.

Romero, E. and J. E. Roman (2014). A parallel implementation of Davidson methods for large-scale eigenvalue problems in SLEPc. , 40(2):13:1--13:29.

Saad, Y. (1992). . John Wiley and Sons, New York.

Sangalli, D., et al. (2019). Many-body perturbation theory calculations using the yambo code. , 31(32):325902.

Scott, D. S. (1982). The advantages of inverted operators in Rayleigh-Ritz approximations. , 3(1):68--75.

Sidje, R. B. (1998). Expokit: a software package for computing matrix exponentials. , 24(1):130--156.

Stathopoulos, A. and J. R. McCombs (2010). : methods and software description. , 37(2):21:1--21:30.

Stewart, G. W. (2001). . SIAM, Philadelphia.

Sukkari, D., H. Ltaief, A. Esposito, and D. Keyes (2019). A QDWH-based SVD software framework on distributed-memory manycore systems. , 45(2):18:1--21.

Tisseur, F. (2000). Backward error and condition of polynomial eigenvalue problems. , 309(1--3):339--361.

Tisseur, F. and K. Meerbergen (2001). The quadratic eigenvalue problem. , 43(2):235--286.

Winkelmann, J., P. Springer, and E. Di Napoli (2019). : Chebyshev accelerated subspace iteration eigensolver for sequences of Hermitian eigenvalue problems. , 45(2):1----34.

Wu, K. and H. Simon (2000). Thick-restart Lanczos method for large symmetric eigenvalue problems. , 22(2):602--616.
:::

{#str}
## SLEPc Technical Reports

(Note: these reports are available through the [SLEPc web site](https://slepc.upv.es).)

::: list
V. Hern&aacute;ndez, J. E. Rom&aacute;n, A. Tom&aacute;s, V. Vidal. "Orthogonalization Routines in SLEPc."

V. Hern&aacute;ndez, J. E. Rom&aacute;n, A. Tom&aacute;s, V. Vidal. "Single Vector Iteration Methods in SLEPc."

V. Hern&aacute;ndez, J. E. Rom&aacute;n, A. Tom&aacute;s, V. Vidal. "Subspace Iteration in SLEPc."

V. Hern&aacute;ndez, J. E. Rom&aacute;n, A. Tom&aacute;s, V. Vidal. "Arnoldi Methods in SLEPc."

V. Hern&aacute;ndez, J. E. Rom&aacute;n, A. Tom&aacute;s, V. Vidal. "Lanczos Methods in SLEPc."

V. Hern&aacute;ndez, J. E. Rom&aacute;n, A. Tom&aacute;s, V. Vidal. "A Survey of Software for Sparse Eigenvalue Problems."

V. Hern&aacute;ndez, J. E. Rom&aacute;n, A. Tom&aacute;s, V. Vidal. "Krylov-Schur Methods in SLEPc."

V. Hern&aacute;ndez, J. E. Rom&aacute;n, A. Tom&aacute;s. "Restarted Lanczos Bidiagonalization for the SVD in SLEPc."

J. E. Rom&aacute;n. "Practical Implementation of Harmonic Krylov-Schur."

M. E. Hochstenbach, E. Romero, J. E. Roman. "Davidson Type Subspace Expansions for the Linear Eigenvalue Problem."

Y. Maeda, T. Sakurai, J. E. Roman. "Contour Integral Spectrum Slicing Method in SLEPc."
:::

```{rubric} Footnotes
```

```{eval-rst}
.. bibliography::
   :filter: docname in docnames
```
