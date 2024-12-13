# Changelog of SLEPc versions

## [unreleased]

### Added

- `SVD`: new support for a stopping test based on a threshold, so that it is possible
  to compute singular values based on how they decay until a certain limit, without
  knowing a priori how many will be finally obtained. See `SVDSetThreshold()`.
- `EPS`: similar threshold stopping test for eigenvalues, see `EPSSetThreshold()`.
- New interface to external package ChASE for standard Hermitian eigenproblems.

## [3.22] - 2024-09-29

### Added

- `EPS`: added support for structured eigenvalue problems in the Krylov-Schur solver.
  In particular, the problem type `EPS_BSE` is intended to solve eigenproblems with
  Bethe-Salpeter structure, whose matrix has been created with `MatCreateBSE()`.
- `BV`: new function `BVGetSplitRows()` to split a `BV` row-wise using index sets.
- Added support for AMD GPUs via HIP and ROCm, see section 8.3 of the Users Manual.

### Changed

- Some interface functions now support `PETSC_CURRENT` to preserve the current value of
  the arguments. This implies that `PETSC_DEFAULT` is deprecated for setting a default
  value in those cases, and `PETSC_DETERMINE` or `PETSC_DECIDE` should be used instead.
  This is the case of, e.g., `EPSSetDimensions()`.
- `DS`: various improvements including optimization of memory usage.
- In slepc4py now `EPS.getEigenpair()` and `EPS.getEigenvalue()` will return a real value
  instead of a complex, if the problem is of Hermitian or generalized Hermitian type.

## [3.21] - 2024-03-30

### Added

- `ST`: new function `STShellSetApplyHermitianTranspose()` and other improvements for handling
  the conjugate transpose of the operator.

### Changed

- `SVDGetSignature()` now takes a `Vec` argument instead of `Vec*`. The former interface
  could not be used from Fortran.

## [3.20] - 2023-09-29

### Changed

- `BV`: now most `BV` implementations use the concept of leading dimension for internal
  storage. The implication of this is that codes relying on `BVGetArray()` must be
  adapted to use the value of the leading dimension, which may be different from the
  local number of rows. Use `BVGetLeadingDimension()` for this purpose.
- `BV`: now the default `BV` type is `BVMAT` instead of `BVSVEC`. This change should not
  affect application code.

## [3.19] - 2023-03-31

### Added

- New `DSHSVD` type for the hyperbolic SVD.
- Support for the hyperbolic SVD (HSVD) in the `SVDTRLANCZOS` solver.
- Add interface to dense SVD solvers in the external package KSVD.

## [3.18] - 2022-10-01

### Added

- `SVD`: add scaling parameter in thick-restart Lanczos for the GSVD, and also enable
  one-sided orthogonalization for GSVD.
- `SVD`: support for the hyperbolic SVD (HSVD), with a signature matrix passed with
  `SVDSetSignature()`.
- `FN`: add support for matrix function computation on the GPU, whenever the `Mat`
  arguments of `FNEvaluateFunctionMat()` are of type `MATSEQDENSECUDA` (and similarly
  for `FNEvaluateFunctionMatVec()`).
- `DS`: new functions `DSGetMatAndColumn()`/`DSRestoreMatAndColumn()`.

### Changed

- In configure, now the arguments intended to store link flags (and libraries) instead of
  being represented as comma-separated lists should be specified via a quoted string, e.g.,
  `--with-arpack-lib="-L/opt/soft/arpack -Wl,-rpath,/opt/soft/arpack -lparpack -larpack"`.
- Interface of `BVMatLanczos()` changed to receive a `Mat` instead of pointers, in a
  similar way to `BVMatArnoldi()`.
- `DS`: now every call to `DSGetMat()` must be paired with a call to `DSRestoreMat()`, it
  is not longer possible to just `MatDestroy()` the matrix.

### Removed

- Remove python2 support, python-3.4+ is now required.
- `DS`: remove `DSCopyMat()`.

## [3.17] - 2022-03-31

### Added

- `NEP`: user-defined preconditioner for the split form, see `NEPSetSplitPreconditioner()`.
  Support for preconditioner matrix different from the function matrix has also been
  improved for the callback interface.
- Add support for split preconditioner in `EPS` and `PEP` via `STSetSplitPreconditioner()`.

### Changed

- In C programs, instead of using `CHKERRQ()`, the recommended way of doing error
  checking is wrapping with `PetscCall()`, for instance, `PetscCall(EPSSolve(eps))`.
  `CHKERRQ()` will continue to work as before. In Fortran, error checking has not changed.
- `SVD`: clarify convergence test for GSVD, now the default is relative to the matrix norms,
  and relative to the singular value is not allowed.

## [3.16] - 2021-09-30

### Added

- `SVD`: support for the GSVD in `SVDTRLANCZOS`, `SVDCROSS` and `SVDCYCLIC` solvers.
- `NEP`: the `NEPCISS` solver (that was deprecated in version 3.14) has been recovered
  and enriched with a new Rayleigh-Ritz variant that is used by default now, see
  `NEPCISSSetExtraction()`. It now also supports subcommunicators.
- `PEP`: new solver `PEPCISS` with similar functionality as `NEPCISS`.
- Fortran: add support for ``SlepcInitialize(filename,help,ierr)``, ``SlepcInitialize(ierr)``
  in addition to current ``SlepcInitialize(filename,ierr)``.
- New `DSGSVD` type for the generalized SVD.
- New contour integral method in `DSNEP`.
- New function `BVSetRandomSign()` for random values equal to -1 or 1 with equal probability.
- New functions `RGIsAxisymmetric()`, `RGCanUseConjugates()` and `RGComputeQuadrature()`.
- New function `SlepcHasExternalPackage()` to check whether SLEPc was configured with the
  given external package.

### Changed

- `DSSetDimensions()` now has one argument less. Use `DSSVDSetDimensions()` to specify
  the number of columns in `DSSVD`. And similarly for `DSGSVD`.
- One of the `DS` matrices has been renamed from `DS_MAT_VT` to `DS_MAT_V` since it is
  no longer transposed.
- Interface of `BVMatArnoldi()` changed to receive a `Mat` instead of a pointer.
- `STShellGetContext()` now takes `void*` as return argument, as well as `XXXGetMonitorContext()`.
- Rename `BVMultInPlaceTranspose()` to `BVMultInPlaceHermitianTranspose()` to be
  consistent with PETSc naming convention.

### Removed

- The interface to BLZPACK has been removed since it is no longer available for download.

## [3.15] - 2021-03-31

### Added

- `SVD`: support for the GSVD added. This implies changes in the interface, e.g.,
  `SVDSetOperators()` should be used instead of `SVDSetOperator()`. New functions have
  been added such as `SVDSetProblemType()` and `SVDIsGeneralized()`.
- `SVD`: New solver `SVDRANDOMIZED`.
- `EPS`: New interface to external package EVSL.
- Import SLEPc4py sources into SLEPc source tree.
- New function `EPSKrylovSchurGetKSP()` to extract the `KSP` object used in the case
  of spectrum slicing runs. This allows setting MUMPS options even with partitions.
- New tolerance to check a definite inner product, see `BVSetDefiniteTolerance()`.
- New function `STSetPreconditionerMat()` to pass the matrix for building the
  preconditioner.
- New function `BVSetRandomNormal()` for random initial vectors with normal distribution.
- New function `DSGetTruncateSize()` to be used with `DSTruncate()`.

### Changed

- Configure option `--download-slepc4py` renamed to `--with-slepc4py`.
- `EPSLANCZOS` now requires that `EPSLanczosSetReorthog()` is called explicitly, or the
  corresponding command-line option is given. The intention is that inexperienced users
  realize that they should use `EPSKRYLOVSCHUR` instead of `EPSLANCZOS`.
- `ST`: now the coefficient/preconditioner matrix has the same prefix as the `KSP` object,
  so setting MUMPS options in the command-line requires the appropriate prefix, e.g.,
  `-st_mat_mumps_icntl_13 1`.
- `DSTruncate()` now has an additional argument `trim`.
- Monitors have been reorganized as in PETSc. Now use `-eps_monitor draw::draw_lg`
  instead of `-eps_monitor_lg`.

### Deprecated

- `STPrecondSetMatForPC()` is superseded by `STSetPreconditionerMat()`.

## [3.14] - 2020-09-30

### Added

- Add interfaces to dense eigensolvers and SVD solvers in ScaLAPACK, ELPA and Elemental.
- Reintroduce interface to FEAST external solver via MKL.
- New configure option `--download-slepc4py`.
- New functions `BVNormalize`, `VecCheckOrthonormality`, `STApplyMat`, `STMatMatSolve`.

### Changed

- `EPSLOBPCG` and `EPSCISS` now use `MatProduct()`/`KSPMatSolve()` instead of
  `MatMult()`/`KSPSolve()`.
- Rename `XXXReasonView()` to `XXXConvergedReasonView()`, and `XXXReasonViewFromOptions()`
  to `XXXConvergedReasonViewFromOptions()`.

### Deprecated

- The `NEPCISS` solver has been deprecated and will be completely removed in future versions.

## [3.13] - 2020-03-31

### Added

- `EPS`: new solver `lyapii` (Lyapunov inverse iteration) to compute rightmost eigenvalues.
- `NEP`: add a two-sided version of SLP.
- New functions: `EPSSetLeftInitialSpace()`, `MFNSolveTranspose()`.
- Interface to `PCHPDDM` and `KSPHPDDM` using `--download-hpddm`.

### Changed

- `configure`: the options `--with-xxx-flags` have been renamed to `--with-xxx-lib` (see
  `./configure --help`).
- `STGetOperator()` now has to be used with `STRestoreOperator()`.
- Now the arrays returned by `RGPolygonGetVertices()` must be freed by the user.

## [3.12] - 2019-09-30

- `ST`: now the default is to use Cholesky instead of LU if the problem matrices are `SBAIJ` or the
  user flags them as symmetric.
- Removed interface to external solver FEAST.
- Removed support for CMake and legacy builds.

## [3.11] - 2019-03-29

- `EPS`: reintroduce support for computing left eigenvectors, see `EPSSetTwoSided()`.
- `PEP`: enable JD solver in real arithmetic.
- `PEP`: spectrum slicing can now be used with symmetric non-hyperbolic problems.
- `NEP`: support for automatic determination of singularities in NLEIGS with the AAA algorithm.
- `NEP`: support for left eigenvectors and approximation of the resolvent.
- `FN`: implemented matrix logarithm.
- Enable compilation with non-GNU compilers in Windows platforms.
- Added support for `DESTDIR` in prefix install.

## [3.10] - 2018-09-18

- `NEP`: add support for deflation in SLP, RII and NArnoldi solvers.
- `FN`: improved matrix function evaluation of Phi functions.
- Interface changes: `PEPLinearSetCompanionForm()` has been deprecated, use
  `PEPLinearSetLinearization()` instead.

## [3.9] - 2018-04-12

- `PEP`: add spectrum slicing for Hermitian quadratic eigenvalue problems via STOAR.
- `PEP`: add support for non-monomial bases in JD.
- `BV`: add `BVTENSOR` to represent subspace bases with compact representation $`V = (I \otimes U) S`$.
- `BV`: improved implementation of block orthogonalization, including TSQR and SVQB.
- `BV`: add `BVGetSplit()` operation to split a given `BV` object in two.

## [3.8] - 2017-10-20

- Released under 2-clause BSD license.
- New solver class `LME` for linear matrix equations (such as Lyapunov or Sylvester) whose
  solution has low rank.
- `NEP`: added specific support for rational eigenvalue problems with `NEPSetProblemType()`.
- Added nonlinear inverse iteration as an option of `EPSPOWER`.
- Added a preliminary implementation of polynomial filters in `STFILTER` to compute interior
  eigenvalues of symmetric problems without factorizing a matrix.
- `SVD`: add wrapper to PRIMME SVD solver.
- Improved evaluation of matrix functions, both in `FN` and `MFN`.
- GPU support improved in all solver classes.
- Simplified Fortran usage as in PETSc.
- Interface changes: `DSNormalize()` has been removed; `NEPInterpolSetDegree()` has been
  renamed to `NEPInterpolSetInterpolation()`, and takes an additional argument.

## [3.7] - 2016-05-16

- `NEP`: new solver `nleigs` that implements a (rational) Krylov method operating on a
  companion-type linearization of a rational interpolant of the nonlinear function.
- `PEP`: the `jd` solver can now compute more than one eigenpair.
- `MFN`: added a new solver `krylov` that works for different functions, not only the exponential.
  The solver available in previous versions has been renamed to `expokit`.
- `EPS`: in spectrum slicing in multi-communicator mode now it is possible to update the problem
  matrices directly on the sub-communicators.
- `EPS`: the contour integral solver now provides Chebyshev quadrature rule and Hankel extraction;
  all options are documented in STR-11.
- Now most solvers allow a user-defined criterion to stop iterating based on a callback function.
- Optimized solution of linear systems in Newton refinement for `PEP` and `NEP`.
- Added download option for most external packages.
- GPU support updated to use `VECCUDA` instead of `VECCUSP`, now including complex scalars.
- Interface changes: `EPS_CONV_EIG` has been renamed to `EPS_CONV_REL`;
  `BVAXPY()` has been removed, use `BVMult()` instead;
  `BVSetRandom()` no longer takes a `PetscRandom` argument, use `BVSetRandomContext()` instead.

## [3.6] - 2015-06-12

- New `EPS` solver: locally optimal block preconditioned conjugate gradient (LOBPCG).
- New `PEP` solvers: Jacobi-Davidson (JD), and symmetric TOAR (STOAR).
- New `NEP` solver: contour integral spectrum slice (CISS).
- Improved BLOPEX interface by adding hard locking. Now the user can specify the block size.
- Spectrum slicing in `EPS` can now be run in multi-communicator mode, where each process group
  computes a sub-interval.
- Added functions and command-line options to view the computed solution after solve, e.g.
  `-eps_view_vectors binary:myvecs.bin`
- The `MFN` solver class now takes an `FN` object to define the function. The functionality for
  computing functions of small, dense matrices has been moved from `DS` to `FN`.
  `MFNSetScaleFactor()` has been removed, now this scale factor must be specified in `FN`.
- `FN` now allows the definition of functions by combining two other functions.
- Added two new `RG` regions: `ring` a stripe along an ellipse with optional start and end angles;
  `polygon` an arbitrary polygon made up of a list of vertices.
- User makefiles must now include `${SLEPC_DIR}/lib/slepc/conf/slepc_common`.
- Interface changes: `EPSComputeRelativeError()` and `EPSComputeResidualNorm()` have been
  deprecated (use `EPSComputeError()` instead); the same for `PEP`, `NEP` and `SVD`;
  `PEPSetScale()` now allows passing two `Vec`s for diagonal scaling; `XXXPrintSolution()` has
  been replaced with `XXXErrorView()`; `STGetOperationCounters()` has been removed, its
  functionality is available via the log summary and with `KSPGetTotalIterations()`;
  `BVGetVec()` has been renamed to `BVCreateVec()`; the interface for defining `FN` functions of
  rational type has changed; `BVSetOrthogonalization()` takes one more argument.

## [3.5] - 2014-07-29

- A new solver class `PEP` for polynomial eigenvalue problems has been added. It replaces the
  former `QEP` class, that has been removed. `PEP` contains a new solver TOAR that can handle
  polynomials of arbitrary degree. Q-Lanczos has been removed since it did not have guaranteed
  stability.
- New `NEP` solver: polynomial interpolation using `PEP`.
- Added Newton iterative refinement in both `PEP` and `NEP`.
- A new auxiliary class `RG` allows the user to define a region in the complex plane. This can
  be used for filtering out unwanted eigenvalues in `EPS` and `PEP`.
- The auxiliary object `IP` has been removed and a new object `BV` has been added that subsumes
  its functionality.
- Support for requesting left eigenvectors has been removed, since none of the solvers were
  computing them internally.
- The `STFOLD` spectral transformation has been removed. Example `ex24.c` reproduces this
  functionality.
- Interface changes: `SVDSetTransposeMode()` has been renamed to `SVDSetImplicitTranspose()`;
  in `STSHIFT` the sign of sigma is the opposite to previous versions; `EPSJDSetBOrth()` now takes
  a boolean argument instead of an enum.

## [3.4] - 2013-07-05

- Added new class of solvers `NEP` for the nonlinear eigenvalue problem.
- Added new class of solvers `MFN` for computing the action of a matrix function on a vector.
- New `EPS` solver: Contour integral spectrum slice (CISS). Allows to compute all eigenvalues
  inside a region.
- New `QEP` solver: Q-Lanczos is a specialized variant of Q-Arnoldi for problems with symmetric
  matrices.
- Added support for shift-and-invert in `QEP`.
- Added a new auxiliary class `FN`: Mathematical Function, to be used in the definition of
  nonlinear eigenproblems.
- Added interface to external solver FEAST.
- Changed options `-xxx_monitor_draw` to `-xxx_monitor_lg`, and similarly for `-xxx_monitor_draw_all`.

## [3.3] - 2012-08-03

- New `EPS` solver: Rayleigh quotient conjugate gradient (RQCG). This is the first CG-type
  eigensolver in SLEPc. It can be used for computing smallest eigenvalues of symmetric-definite
  matrix pairs without inverting any matrix (a preconditioner can be used instead).
- Added a customizable parameter to specify how to restart in Krylov-Schur, see
  `EPSKrylovSchurSetRestart`. Tuning this parameter may speed up convergence significantly in
  some cases.
- Added support for generalized symmetric-indefinite eigenproblems in Krylov-Schur and the Davidson
  solvers. To use this, set the problem type to `EPS_GHIEP`.
- New variant of Generalized Davidson for generalized eigenproblems that expands the subspace with
  two vectors (GD2). It can be activated with `-eps_gd_double_expansion`.
- Added experimental support for arbitrary selection of eigenpairs, where the solver chooses the
  most wanted ones based on a user-defined function of the eigenvectors rather than simply sorting
  the eigenvalues.
- Added a new auxiliary class `DS`: Direct Solver (or Dense System), which is intended for
  developers rather than normal users.

## [3.2] - 2011-10-28

- Computational intervals for symmetric eigenproblems, that activate a spectrum slicing mechanism
  to obtain all eigenvalues in a given interval, see `EPSSetInterval()`.
- Partial support (experimental) for GPU computing via PETSc's `VECCUSP` and `MATAIJCUSP`.
- Improved performance and robustness of GD and JD solvers in (generalized) Hermitian problems.
- Performance improvement of solvers with explicit matrix such as `SVDCYCLIC` and `QEPLINEAR`
  (now use matrix preallocation).
- Added Matlab interface.
- Added support for parallel builds with CMake.
- Added support for quad precision (configure PETSc `--with-precision=__float128` with gcc-4.6 or
  later).
- Interface changes: now all `XXXDestroy()` routines take a pointer to the object.

## [3.1] - 2010-08-04

- New `EPS` solvers: Generalized Davidson (GD) and Jacobi-Davidson (JD). These are the first
  eigensolvers in SLEPc that belong to the class of preconditioned eigensolvers.
- Added a new instance of `ST` called `STPRECOND`. This is not really a spectral transformation
  but rather a convenient way of handling the preconditioners in the new eigensolvers.
- Added a new class `QEP` for solving quadratic eigenvalue problems. Currently, it contains two
  solvers: the Q-Arnoldi method and another one that provides a linearization of the problem and
  then invokes an eigensolver from `EPS`.
- Added support for balancing of non-Hermitian problems, see `EPSSetBalance()`.
- Improved sorting of eigenvalues, now with the possibility of sorting with respect to a target
  value. With shift-and-invert, now the ordering of eigenvalues is the expected one, relative to
  the target. Also added support for user-defined orderings, see `EPSSetWhichEigenpairs()`.
- Added support for user-defined convergence tests, see `EPSSetConvergenceTest()`. Several
  predefined convergence criteria are available. Also, there is a new flag for computing the true
  residual for the convergence test, see `EPSSetTrueResidual()`.
- Monitors have been reorganized and now more possibilities are available.
- Changes in user interface: `EPSAttachDeflationSpace()` has been renamed to `EPSSetDeflationSpace()`,
  `EPSSetLeftVectorsWanted()` replaces `EPSSetClass()` for requesting left eigenvectors; Change in
  arguments: `EPSSetDeflationSpace()`; Deprecated function: `EPSSetInitialVector()`, replaced by
  `EPSSetInitialSpace()`; `STSINV` has been renamed to `STSINVERT`.

## [3.0.0] - 2009-02-04

- Released under GNU LGPL license.
- Improved support for the case that many eigenpairs are to be computed. This is especially so in
  the default eigensolver (Krylov-Schur) for symmetric problems, as well as for SVD computations. The
  user can control the behaviour of the solver with a new parameter `mpd`, see `EPSSetDimensions()`.
- Support for harmonic projection in the default eigensolver (Krylov-Schur), see `EPSSetExtraction()`.
  This can be useful for computing interior or rightmost eigenvalues without the need of a spectral
  transformation.
- Memory usage has been optimized in most solvers. In some cases, memory requirements have been
  halved with respect to the previous versions.
- In the spectral transformations (`ST`) the linear solver used internally has been changed to a
  direct solver by default. The user can still employ an iterative linear solver by setting the
  appropriate options.
- Added better support for Fortran 90.
- Improved support for `make install`.

## [2.3.3] - 2007-06-01

- A new solver class `SVD` has been introduced for computing the singular value decomposition of
  a rectangular matrix. The structure of this new type is very similar to that of `EPS`, and it
  simplifies the computation of singular values and vectors.
- Better support for generalized problems. Eigenvector purification has been added to improve
  accuracy in the case of generalized eigenproblems with singular B. Also, a new problem type
  (`EPS_PGNHEP`) has been added for better addressing generalized problems in which A is non-Hermitian
  but B is Hermitian and positive definite.
- Now `make install` is available thus facilitating system-wide installation.

## [2.3.2] - 2006-10-02

- A new `krylovschur` eigensolver has been added, that implements the Krylov-Schur method. This
  method is related to the Arnoldi and Lanczos algorithms, but incorporates a new restarting scheme
  that makes it competitive with respect to implicit restart. This eigensolver is now the default
  for both symmetric and non-symmetric problems.
- A new wrapper has been developed to interface with the PRIMME library. This library provides
  Davidson-type eigensolvers.
- The `lanczos` solver has been improved, in particular, the different reorthogonalization
  strategies are now more robust.
- Now the `arnoldi` eigensolver supports the computation of eigenvalues other than those of largest
  magnitude.
- `EPSGetLinearIterations()` has been replaced with `EPSGetOperationCounters()`, providing more
  statistics about the solution process.
- `EPSGetIterationNumber()` now returns the number corresponding to outer iterations.
- The `lobpcg` wrapper has been renamed to `blopex`.
- The `planso` wrapper has been removed since PLANSO is no longer being distributed.

## [2.3.1] - 2006-03-03

- New variant of the Arnoldi method added to the `arnoldi` eigensolver (with delayed
  reorthogonalization, see `EPSArnoldiSetDelayed()`).
- Several optimizations for improving performance and scalability, in particular the
  orthogonalization steps.

## [2.3.0] - 2005-06-30

- New `lanczos` eigensolver, an explicitly restarted version of the Lanczos method for symmetric
  eigenproblems. Allows the user to choose among 5 reorthogonalization strategies.
- New spectrum folding spectral transformation.
- New configuration system, similar to PETSc's `configure.py`.
- New interface to an external eigensolver: LOBPCG implemented in Hypre.
- Added graphical convergence monitor (with `-eps_xmonitor`).
- Improvement of Arnoldi solver in terms of efficiency and robustness.
- Now the `lapack` solver uses specific Lapack routines for symmetric and generalized problems.
- Bug fix in the ARPACK interface.

## [2.2.1] - 2004-08-21

- The `power` eigensolver has been replaced by a simpler implementation.
- The `rqi` eigensolver has been removed. Now the Rayleigh Quotient Iteration is embedded in the
  `power` method.
- The `subspace` eigensolver has been rewritten. Now it follows the SRRIT implementation, which is
  much faster than the old one.
- The `arnoldi` eigensolver has been re-implemented as well. The new implementation is much more
  robust and efficient.
- A new Spectral Transformation (`ST`) has been added: the generalized Cayley transform.
- Support for user-provided deflation subspaces has been added, see `EPSAttachDeflationSpace()`.
- Support for preservation of symmetry in eigensolvers. For this feature, the user must explicitly
  call `EPSSetProblemType()` in symmetric problems.
- The two types of monitors (error estimates and values) have been merged in a single one.
- New function `EPSGetInvariantSubspace()`.
- Better support for spectrum slicing in `blzpack`.

## [2.2.0] - 2004-04-13

- `EPSSolve()` does not return the number of iterations. Use `EPSGetIterationNumber()` for this
  purpose.
- `EPSGetSolution()` has been replaced by `EPSGetEigenpair()` with a cleaner interface.
- `EPSComputeError()` has been replaced by `EPSComputeRelativeError()` and `EPSComputeResidualNorm()`
  with better error computing for zero eigenvalues. These functions now are oriented to single
  eigenpairs, as well as `EPSGetErrorEstimate()`.
- The possibilities of `EPSSetWhichEigenpairs()` have been reduced and now are more coherent across
  problem types.
- Removed `STNONE` spectral transformation. The default of `STSHIFT` with 0 shift is equivalent.
- Added `STSinvertSetMatStructure()` to optimize performance of `MatAXPY()` in shift-and-invert
  transformation.
- Classical and modified Gram-Schmidt orthogonalization use iterative refinement, with user options
  for parameter adjustment.

[unreleased]: https://gitlab.com/slepc/slepc/compare/v3.22...main
[3.22]: https://gitlab.com/slepc/slepc/compare/v3.21...v3.22
[3.21]: https://gitlab.com/slepc/slepc/compare/v3.20...v3.21
[3.20]: https://gitlab.com/slepc/slepc/compare/v3.19...v3.20
[3.19]: https://gitlab.com/slepc/slepc/compare/v3.18...v3.19
[3.18]: https://gitlab.com/slepc/slepc/compare/v3.17...v3.18
[3.17]: https://gitlab.com/slepc/slepc/compare/v3.16...v3.17
[3.16]: https://gitlab.com/slepc/slepc/compare/v3.15...v3.16
[3.15]: https://gitlab.com/slepc/slepc/compare/v3.14...v3.15
[3.14]: https://gitlab.com/slepc/slepc/compare/v3.13...v3.14
[3.13]: https://gitlab.com/slepc/slepc/compare/v3.12...v3.13
[3.12]: https://gitlab.com/slepc/slepc/compare/v3.11...v3.12
[3.11]: https://gitlab.com/slepc/slepc/compare/v3.10...v3.11
[3.10]: https://gitlab.com/slepc/slepc/compare/v3.9...v3.10
[3.9]: https://gitlab.com/slepc/slepc/compare/v3.8...v3.9
[3.8]: https://gitlab.com/slepc/slepc/compare/v3.7...v3.8
[3.7]: https://gitlab.com/slepc/slepc/compare/v3.6...v3.7
[3.6]: https://gitlab.com/slepc/slepc/compare/v3.5...v3.6
[3.5]: https://gitlab.com/slepc/slepc/compare/v3.4...v3.5
[3.4]: https://gitlab.com/slepc/slepc/compare/v3.3...v3.4
[3.3]: https://gitlab.com/slepc/slepc/compare/v3.2...v3.3
[3.2]: https://gitlab.com/slepc/slepc/compare/v3.1...v3.2
[3.1]: https://gitlab.com/slepc/slepc/compare/v3.0.0...v3.1
[3.0.0]: https://gitlab.com/slepc/slepc/compare/v2.3.3...v3.0.0
[2.3.3]: https://gitlab.com/slepc/slepc/compare/v2.3.2...v2.3.3
[2.3.2]: https://gitlab.com/slepc/slepc/compare/v2.3.1...v2.3.2
[2.3.1]: https://gitlab.com/slepc/slepc/compare/v2.3.0...v2.3.1
[2.3.0]: https://gitlab.com/slepc/slepc/compare/v2.2.1...v2.3.0
[2.2.1]: https://gitlab.com/slepc/slepc/compare/v2.2.0...v2.2.1
[2.2.0]: https://gitlab.com/slepc/slepc/-/tags/v2.2.0
