# About

%SLEPc,
The Scalable Library for Eigenvalue Problem Computations is a software library for the solution of large sparse eigenproblems on parallel computers. It can be used for the solution of linear eigenvalue problems formulated in either standard or generalized form, as well as other related problems such as the singular value decomposition. It can also be used to solve nonlinear eigenvalue problems, either those formulated as polynomial eigenproblems or more general nonlinear problems. Finally, SLEPc provides solvers for the computation of the action of a matrix function on a vector.

The emphasis of the software is on methods and techniques appropriate for problems in which the associated matrices are sparse, for example, those arising after the discretization of partial differential equations. Therefore, most of the methods offered by the library are projection methods or other methods with similar properties. Examples of these methods are Krylov-Schur, Jacobi-Davidson, and Conjugate Gradient, to name a few. SLEPc provides implementations of these methods. It also provides built-in support for spectral transformations such as the shift-and-invert technique. SLEPc is a general library in the sense that it covers both Hermitian and non-Hermitian problems, with either real or complex arithmetic.

SLEPc is built on top of [PETSc](https://petsc.org), the Portable, Extensible Toolkit for Scientific Computation. It can be considered an extension of PETSc providing all the functionality necessary for the solution of eigenvalue problems. This means that PETSc must be previously installed in order to use SLEPc. PETSc users will find SLEPc very easy to use, since it enforces the same programming paradigm. For those users who are not familiar with PETSc yet, our recommendation is to fully understand its basic concepts before proceeding with SLEPc. Parallelism in SLEPc is obtaned by means of MPI, and in addition there is some support for GPU computing.

SLEPc can be used directly, and also indirectly via many [software packages that interface to it](../material/software). In turn, SLEPc interfaces to other external software packages such as [ARPACK](https://github.com/opencollab/arpack-ng) or [BLOPEX](https://github.com/lobpcg/blopex).

%## Description

%SLEPc, the Scalable Library for Eigenvalue Problem Computations, is a software library for the solution of large, sparse eigenproblems on parallel computers. It can be used for the solution of linear eigenvalue problems formulated in either standard or generalized form, both Hermitian and non-Hermitian, with either real or complex arithmetic, as well as other related problems such as the singular value decomposition or the nonlinear eigenvalue problem.

%SLEPc focuses on sparse problems, for example, those arising after the discretization of partial differential equations. It provides projection methods or other methods with similar properties, such as Krylov-Schur or Jacobi-Davidson. It also provides built-in support for different types of problems and spectral transformations such as the shift-and-invert technique.

%SLEPc is built on top of [PETSc](https://petsc.org) (Portable, Extensible Toolkit for Scientific Computation) and extends it with all the functionality necessary for the solution of eigenvalue problems.

## Linear Eigensolvers

> Parallel Solution of Large-scale Sparse Linear Eigenvalue Problems
>```{math}
>Ax=\lambda x
>```
>```{math}
>Ax=\lambda Bx
>```

The solvers for linear eigenvalue problems currently provided by SLEPc are the following:

  * **Krylov-Schur** , which is a variation of Arnoldi/Lanczos with a very effective restarting technique.
  * Generalized Davidson, a simple iteration based on the subspace expansion by the preconditioned residual.
  * Jacobi-Davidson, a preconditioned eigensolver with an effective correction equation.
  * Locally optimal block preconditioned conjugate gradient (LOBPCG).
  * Rayleigh quotient conjugate gradient minimization with Rayleigh-Ritz acceleration.
  * Contour integral spectrum slicing, based on Sakurai-Sugiura method.
  * Basic solvers:
    * Power Iteration with deflation. When combined with shift-and-invert, it is equivalent to the Inverse Iteration. Also, this solver embeds the Rayleigh Quotient Iteration (RQI) by allowing variable shifts.
    * Subspace Iteration with Rayleigh-Ritz projection and locking.
    * Arnoldi method with explicit restart and deflation.
    * Lanczos method for symmetric/Hermitian problems, with explicit restart and deflation, using different reorthogonalization strategies.

## Nonlinear Eigensolvers

>Polynomial and General Nonlinear Eigenproblems
>```{math}
>P(\lambda)x=0
>```
>```{math}
>T(\lambda)x=0
>```

Apart from the linear eigenvalue solvers, SLEPc also provides specific solvers for nonlinear eigenvalue problems, either polynomial (PEP) or general nonlinear (NEP).

We provide several Krylov subspace methods specific for the polynomial eigenvalue problem, mainly Q-Arnoldi (and derivatives thereof) and a polynomial Jacobi-Davidson. It is also possible to linearize the polynomial problem and use any of the above mentioned linear eigensolvers.

For general nonlinear eigenproblems, the user provides the nonlinear function _T( *)_ either by means of callback routines or in split form (matrices combined with simple nonlinear functions). So far, the following solvers are available:

  * NLEIGS.
  * Contour integral.
  * Polynomial interpolation that uses a PEP solver.
  * Basic solvers:
    * Residual inverse iteration with varying shift.
    * Successive linear problems, a Newton-type iteration based on first order Taylor approximation.
    * Nonlinear Arnoldi.

## Solvers for SVD and Matrix Functions

>Partial SVD and Matrix Function
>```{math}
>Av=\sigma u
>```
>```{math}
>w=f(\alpha A)v
>```

In addition to eigensolvers, SLEPc includes functionality for computing the partial singular value decomposition (SVD) of a large sparse rectangular matrix, and also for computing the action of the function of a large sparse matrix on a vector (MFN).

For the computation of singular values and vectors, any of the above mentioned linear eigensolvers can be used, together with some specific SVD solvers such as Lanczos and thick-restart Lanczos bidiagonalization.

Regarding the matrix function computation, SLEPc currently provides a basic Krylov method for computing the action of a matrix function on a vector, where the function can be the exponential, the square root, a rational function, or combinations thereof.

## Highlights

  * Easy programming with PETSc's object-oriented style. See an example in {{'[C](https://slepc.upv.es/{}/src/eps/tutorials/ex1.c.html)'.format(branch)}} and {{'[Fortran](https://slepc.upv.es/{}/src/eps/tutorials/ex1f.F90.html)'.format(branch)}}. More examples are available in the documentation section.
  * Data-structure neutral implementation. Problems can be solved with matrices stored in parallel and serial, sparse and dense formats, and even without explicit storage.
  * Run-time flexibility, giving full control over the solution process. See some command-line examples below.
  * Portability to a wide range of parallel platforms, including Linux clusters, IBM Bluegene, MacOSX, and many others.
  * Usable from code written in C, C++, Fortran77, and Fortran90, as well as python (note that slepc4py is now part of the SLEPc distribution).
  * Extensive documentation, including a users manual, on-line tutorial exercises, example programs and on-line manual pages for every subroutine.
  * Seamless integration of other eigensolver packages such as [ARPACK](https://github.com/opencollab/arpack-ng) or [BLOPEX](https://github.com/lobpcg/blopex).

## Command-line Examples

The following examples illustrate the use of a SLEPc program via command-line options. All options are also available with an equivalent programmatic interface.

```{code} console
$ ./exeps -eps_view -eps_monitor
$ ./exeps -eps_type krylovschur -eps_smallest_real
$ ./exeps -eps_type lanczos -eps_nev 6 -eps_ncv 24
$ ./exeps -eps_type krylovschur -eps_nev 100 -eps_mpd 40
$ ./exeps -eps_type arnoldi -eps_tol 1e-8 -eps_max_it 2000
$ ./exeps -eps_type gd -eps_gd_plusk 1
$ ./exeps -eps_type jd -eps_jd_blocksize 2 -eps_jd_krylov_start
$ ./exeps -eps_type subspace -eps_hermitian -log_view
$ ./exeps -eps_type arpack -eps_view_values draw -draw_pause -1

$ ./exeps -eps_target 2.1
$ ./exeps -eps_target 2.1 -eps_harmonic
$ ./exeps -eps_target 2.1 -st_type sinvert

$ ./exsvd -svd_type trlanczos -svd_nsv 4
$ ./exsvd -svd_type cross -svd_cross_eps_type krylovschur -svd_smallest -svd_monitor_lg -svd_ncv 30

$ ./expep -pep_nev 6 -pep_tol 1e-5
$ ./expep -pep_type linear -pep_linear_linearization 1,0 -pep_linear_explicitmatrix
$ ./expep -pep_type toar
```

## Scheme of SLEPc Classes

The following scheme represents the functionality provided by SLEPc and how it relates to PETSc.

  * SLEPc provides five main objects: EPS (Eigenvalue Problem Solver), SVD (Singular Value Decomposition), PEP (Polynomial Eigenvalue Problem), NEP (Nonlinear Eigenvalue Problem), and MFN (Matrix Function).
  * These objects occupy a level of abstraction similar to other PETSc solvers such as KSP or SNES and use low-level infrastructure such as Mat and Vec.
  * The shaded blocks represent the generic interface of the object while the white boxes represent different implementations. The programmer usually interacts with the object via its interface and the particular implementation is typically picked at run time.
  * ST (Spectral Transformation) is used in combination of most solvers to compute interior eigenvalues.
  * Additionally, several auxiliary classes are provided: BV (Basis Vectors), DS (Dense System), RG (Region), FN (Mathematical Function).

![petsc-slepc](../_static/images/petsc-slepc-3.7.png)

%```{toctree}
%:maxdepth: 2
%:hidden:
%
%to_do
%```
