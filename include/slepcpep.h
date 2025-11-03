/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for SLEPc's polynomial eigenvalue solvers
*/

#pragma once

#include <slepceps.h>

/* SUBMANSEC = PEP */

SLEPC_EXTERN PetscErrorCode PEPInitializePackage(void);
SLEPC_EXTERN PetscErrorCode PEPFinalizePackage(void);

/*S
   PEP - SLEPc object that manages all the polynomial eigenvalue problem solvers.

   Level: beginner

.seealso: [](ch:pep), `PEPCreate()`
S*/
typedef struct _p_PEP* PEP;

/*J
   PEPType - String with the name of a polynomial eigensolver.

   Level: beginner

.seealso: [](ch:pep), `PEPSetType()`, `PEP`
J*/
typedef const char *PEPType;
#define PEPTOAR      "toar"
#define PEPSTOAR     "stoar"
#define PEPQARNOLDI  "qarnoldi"
#define PEPLINEAR    "linear"
#define PEPJD        "jd"
#define PEPCISS      "ciss"

/* Logging support */
SLEPC_EXTERN PetscClassId PEP_CLASSID;

/*E
   PEPProblemType - Determines the type of the polynomial eigenproblem.

   Values:
+  `PEP_GENERAL`    - polynomial eigenproblem with no particular structure
.  `PEP_HERMITIAN`  - polynomial eigenproblem with all coefficient matrices Hermitian
.  `PEP_HYPERBOLIC` - quadratic eigenproblem with hyperbolic structure
-  `PEP_GYROSCOPIC` - quadratic eigenproblem with gyroscopic structure

   Note:
   By default, no particular structure is assumed (`PEP_GENERAL`).

   Level: intermediate

.seealso: [](ch:pep), `PEPSetProblemType()`, `PEPGetProblemType()`
E*/
typedef enum { PEP_GENERAL    = 1,
               PEP_HERMITIAN  = 2,
               PEP_HYPERBOLIC = 3,
               PEP_GYROSCOPIC = 4
             } PEPProblemType;

/*MC
   PEP_GENERAL - A polynomial eigenproblem with no particular structure.

   Note:
   This is the default problem type.

   Level: intermediate

.seealso: [](ch:pep), `PEPProblemType`, `PEPSetProblemType()`, `PEP_HERMITIAN`, `PEP_HYPERBOLIC`, `PEP_GYROSCOPIC`
M*/

/*MC
   PEP_HERMITIAN - A polynomial eigenvalue problem with all coefficient matrices Hermitian.

   Notes:
   This is used when all $A_i$ matrices passed in `PEPSetOperators()`
   are Hermitian.

   Currently there is support only for quadratic eigenvalue problems,
   $(K+\lambda C+\lambda^2M)x=0$ with $K$, $C$, $M$ Hermitian.

   Level: intermediate

.seealso: [](ch:pep), `PEPProblemType`, `PEPSetProblemType()`, `PEP_GENERAL`, `PEP_HYPERBOLIC`, `PEP_GYROSCOPIC`
M*/

/*MC
   PEP_HYPERBOLIC - A quadratic eigenvalue problem with hyperbolic structure.

   Note:
   This is reserved for the case of a quadratic eigenvalue problem
   $(K+\lambda C+\lambda^2M)x=0$ with Hermitian coefficient matrices, and in
   addition $M$ is positive definite and $(x^*Cx)^2>4(x^*Mx)(x^*Kx)$ for all
   nonzero $x\in\mathbb{C}^n$. All eigenvalues are real, and form two separate
   groups of $n$ eigenvalues, each of them having linearly independent eigenvectors.

   Level: intermediate

.seealso: [](ch:pep), `PEPProblemType`, `PEPSetProblemType()`, `PEP_GENERAL`, `PEP_HERMITIAN`, `PEP_GYROSCOPIC`
M*/

/*MC
   PEP_GYROSCOPIC - A quadratic eigenvalue problem with gyroscopic structure.

   Notes:
   This is reserved for the case of a quadratic eigenvalue problem
   $(K+\lambda C+\lambda^2M)x=0$ with $M$, $K$ Hermitian, $M>0$, and
   $C$ skew-Hermitian.

   Currently there is support for this problem type only in `PEPLINEAR`, using
   a general eigensolver, without exploiting the structure.

   Level: intermediate

.seealso: [](ch:pep), `PEPProblemType`, `PEPSetProblemType()`, `PEP_GENERAL`, `PEP_HERMITIAN`, `PEP_HYPERBOLIC`, `PEPLINEAR`
M*/

/*E
   PEPWhich - Determines which part of the spectrum is requested.

   Values:
+  `PEP_LARGEST_MAGNITUDE`  - largest $|\lambda|$
.  `PEP_SMALLEST_MAGNITUDE` - smallest $|\lambda|$
.  `PEP_LARGEST_REAL`       - largest $\mathrm{Re}(\lambda)$
.  `PEP_SMALLEST_REAL`      - smallest $\mathrm{Re}(\lambda)$
.  `PEP_LARGEST_IMAGINARY`  - largest $\mathrm{Im}(\lambda)$
.  `PEP_SMALLEST_IMAGINARY` - smallest $\mathrm{Im}(\lambda)$
.  `PEP_TARGET_MAGNITUDE`   - smallest $|\lambda-\tau|$
.  `PEP_TARGET_REAL`        - smallest $|\mathrm{Re}(\lambda-\tau)|$
.  `PEP_TARGET_IMAGINARY`   - smallest $|\mathrm{Im}(\lambda-\tau)|$
.  `PEP_ALL`                - all $\lambda\in[a,b]$ or $\lambda\in\Omega$
-  `PEP_WHICH_USER`         - user-defined sorting criterion

   Notes:
   If SLEPc is compiled for real scalars `PEP_LARGEST_IMAGINARY` and
   `PEP_SMALLEST_IMAGINARY` use the absolute value of the imaginary part
   for eigenvalue selection.

   The target $\tau$ is a scalar value provided with `PEPSetTarget()`.

   The case `PEP_ALL` needs an interval $[a,b]$ given with `PEPSetInterval()`
   or a region $\Omega$ specified with an `RG` object.

   Level: intermediate

.seealso: [](ch:pep), `PEPSetWhichEigenpairs()`, `PEPSetTarget()`, `PEPSetInterval()`
E*/
typedef enum { PEP_LARGEST_MAGNITUDE  = 1,
               PEP_SMALLEST_MAGNITUDE = 2,
               PEP_LARGEST_REAL       = 3,
               PEP_SMALLEST_REAL      = 4,
               PEP_LARGEST_IMAGINARY  = 5,
               PEP_SMALLEST_IMAGINARY = 6,
               PEP_TARGET_MAGNITUDE   = 7,
               PEP_TARGET_REAL        = 8,
               PEP_TARGET_IMAGINARY   = 9,
               PEP_ALL                = 10,
               PEP_WHICH_USER         = 11 } PEPWhich;

/*E
   PEPBasis - The type of polynomial basis used to represent the polynomial
   eigenproblem.

   Values:
+  `PEP_BASIS_MONOMIAL`   - monomial basis
.  `PEP_BASIS_CHEBYSHEV1` - Chebyshev polynomials of the 1st kind
.  `PEP_BASIS_CHEBYSHEV2` - Chebyshev polynomials of the 2nd kind
.  `PEP_BASIS_LEGENDRE`   - Legendre polynomials
.  `PEP_BASIS_LAGUERRE`   - Laguerre polynomials
-  `PEP_BASIS_HERMITE`    - Hermite polynomials

   Notes:
   The default is to work with the monomial basis to represent the polynomial,
   i.e., $1, x, x^2, \dots, x^d$. For large degree $d$, numerical
   difficulties may arise. In that case, a different basis is recommended.
   The user is responsible for providing the coefficient matrices to
   `PEPSetOperators()` represented in the selected basis.

   Level: intermediate

.seealso: [](ch:pep), `PEPSetBasis()`
E*/
typedef enum { PEP_BASIS_MONOMIAL,
               PEP_BASIS_CHEBYSHEV1,
               PEP_BASIS_CHEBYSHEV2,
               PEP_BASIS_LEGENDRE,
               PEP_BASIS_LAGUERRE,
               PEP_BASIS_HERMITE } PEPBasis;
SLEPC_EXTERN const char *PEPBasisTypes[];

/*E
   PEPScale - The scaling strategy.

   Values:
+  `PEP_SCALE_NONE`     - no scaling
.  `PEP_SCALE_SCALAR`   - multiply by a scalar value
.  `PEP_SCALE_DIAGONAL` - multiply by two diagonal matrices
-  `PEP_SCALE_BOTH`     - both scalar and diagonal scaling

   Note:
   See section [](#sec:scaling) for a discussion of the different scaling strategies.

   Level: intermediate

.seealso: [](ch:pep), [](#sec:scaling), `PEPSetScale()`
E*/
typedef enum { PEP_SCALE_NONE,
               PEP_SCALE_SCALAR,
               PEP_SCALE_DIAGONAL,
               PEP_SCALE_BOTH } PEPScale;
SLEPC_EXTERN const char *PEPScaleTypes[];

/*E
   PEPRefine - The type of Newton iterative refinement.

   Values:
+  `PEP_REFINE_NONE`     - no refinement
.  `PEP_REFINE_SIMPLE`   - refinement of each converged eigenpair individually
-  `PEP_REFINE_MULTIPLE` - refinement of the invariant pair as a whole

   Note:
   See section [](#sec:refine) for a discussion of the different refinement strategies.

   Level: intermediate

.seealso: [](ch:pep), [](#sec:refine), `PEPSetRefine()`
E*/
typedef enum { PEP_REFINE_NONE,
               PEP_REFINE_SIMPLE,
               PEP_REFINE_MULTIPLE } PEPRefine;
SLEPC_EXTERN const char *PEPRefineTypes[];

/*E
   PEPRefineScheme - The scheme used for solving linear systems during iterative refinement.

   Values:
+  `PEP_REFINE_SCHEME_SCHUR`    - use the Schur complement
.  `PEP_REFINE_SCHEME_MBE`      - use the mixed block elimination (MBE) scheme
-  `PEP_REFINE_SCHEME_EXPLICIT` - build the full matrix explicitly

   Note:
   Iterative refinement may be very costly, due to the expensive linear
   solves. These linear systems have a particular structure that can be
   exploited in different ways, as described in {cite:p}`Cam16b`. See
   `PEPSetRefine()` for additional details.

   Level: intermediate

.seealso: [](ch:pep), [](#sec:refine), `PEPSetRefine()`
E*/
typedef enum { PEP_REFINE_SCHEME_SCHUR    = 1,
               PEP_REFINE_SCHEME_MBE      = 2,
               PEP_REFINE_SCHEME_EXPLICIT = 3 } PEPRefineScheme;
SLEPC_EXTERN const char *PEPRefineSchemes[];

/*E
   PEPExtract - The eigenvector extraction strategy.

   Values:
+  `PEP_EXTRACT_NONE`       - trivial extraction
.  `PEP_EXTRACT_NORM`       - extraction based on the norm
.  `PEP_EXTRACT_RESIDUAL`   - extraction based on the residual
-  `PEP_EXTRACT_STRUCTURED` - extraction using a linear combination of all the blocks

   Note:
   This is relevant for solvers based on linearization. Once the solver has
   converged, the polynomial eigenvectors can be extracted from the
   eigenvectors of the linearized problem in different ways. See the
   discussion in section [](#sec:pepextr).

   Level: intermediate

.seealso: [](ch:pep), [](#sec:pepextr), `PEPSetExtract()`
E*/
typedef enum { PEP_EXTRACT_NONE       = 1,
               PEP_EXTRACT_NORM       = 2,
               PEP_EXTRACT_RESIDUAL   = 3,
               PEP_EXTRACT_STRUCTURED = 4 } PEPExtract;
SLEPC_EXTERN const char *PEPExtractTypes[];

/*MC
   PEP_EXTRACT_NONE - Trivial eigenvector extraction.

   Note:
   Given the eigenvector of the linearization, $y$, the eigenvector of the
   polynomial eigenproblem $x$ is taken from the first block of $y$.

   Level: intermediate

.seealso: [](ch:pep), [](#sec:pepextr), `PEPExtract`, `PEPSetExtract()`, `PEP_EXTRACT_NORM`, `PEP_EXTRACT_RESIDUAL`, `PEP_EXTRACT_STRUCTURED`
M*/

/*MC
   PEP_EXTRACT_NORM - Eigenvector extraction based on the norm.

   Notes:
   Given the eigenvector of the linearization, $y$, the eigenvector of the
   polynomial eigenproblem $x$ is obtained from the $i$th block for which
   $|\phi_i(\lambda)|$ is maximum, where $\phi_i$ is the $i$th element
   of the polynomial basis. In the case of the monomial basis, it will
   select the first or last block depending on $|\lambda|\geq 1$ or $|\lambda|<1$.

   This is the default extraction.

   Level: intermediate

.seealso: [](ch:pep), [](#sec:pepextr), `PEPExtract`, `PEPSetExtract()`, `PEP_EXTRACT_NONE`, `PEP_EXTRACT_RESIDUAL`, `PEP_EXTRACT_STRUCTURED`
M*/

/*MC
   PEP_EXTRACT_RESIDUAL - Eigenvector extraction based on the residual.

   Note:
   Given the eigenvector of the linearization, $y$, the eigenvector of the
   polynomial eigenproblem $x$ is the block of $y$ that minimizes the
   residual norm.

   Level: intermediate

.seealso: [](ch:pep), [](#sec:pepextr), `PEPExtract`, `PEPSetExtract()`, `PEP_EXTRACT_NONE`, `PEP_EXTRACT_NORM`, `PEP_EXTRACT_STRUCTURED`
M*/

/*MC
   PEP_EXTRACT_STRUCTURED - Eigenvector extraction using a linear combination of
   all the blocks.

   Note:
   Given the eigenvector of the linearization, $y$, the eigenvector of the
   polynomial eigenproblem $x$ is obtained as a linear combination of all
   the blocks, such that it minimizes a certain norm.

   Level: intermediate

.seealso: [](ch:pep), [](#sec:pepextr), `PEPExtract`, `PEPSetExtract()`, `PEP_EXTRACT_NONE`, `PEP_EXTRACT_NORM`, `PEP_EXTRACT_RESIDUAL`
M*/

/*E
   PEPErrorType - The error type used to assess the accuracy of computed solutions.

   Values:
+  `PEP_ERROR_ABSOLUTE` - compute error bound as $\|r\|$
.  `PEP_ERROR_RELATIVE` - compute error bound as $\|r\|/|\lambda|$
-  `PEP_ERROR_BACKWARD` - compute error bound as $\|r\|/(\sum_j\|A_j\||\lambda_i|^j)$

   Level: intermediate

.seealso: [](ch:pep), `PEPComputeError()`
E*/
typedef enum { PEP_ERROR_ABSOLUTE,
               PEP_ERROR_RELATIVE,
               PEP_ERROR_BACKWARD } PEPErrorType;
SLEPC_EXTERN const char *PEPErrorTypes[];

/*E
   PEPConv - The convergence criterion to be used by the solver.

   Values:
+  `PEP_CONV_ABS`  - absolute convergence criterion, $\|r\|$
.  `PEP_CONV_REL`  - convergence criterion relative to eigenvalue, $\|r\|/|\lambda|$
.  `PEP_CONV_NORM` - convergence criterion relative to matrix norms, $\|r\|/(\sum_j\|A_j\||\lambda|^j)$
-  `PEP_CONV_USER` - convergence dictated by user-provided function

   Level: intermediate

.seealso: [](ch:pep), `PEPSetConvergenceTest()`, `PEPSetConvergenceTestFunction()`
E*/
typedef enum { PEP_CONV_ABS,
               PEP_CONV_REL,
               PEP_CONV_NORM,
               PEP_CONV_USER } PEPConv;

/*E
   PEPStop - The stopping test to decide the termination of the outer loop
   of the eigensolver.

   Values:
+  `PEP_STOP_BASIC` - default stopping test
-  `PEP_STOP_USER`  - user-provided stopping test

   Level: advanced

.seealso: [](ch:pep), `PEPSetStoppingTest()`, `PEPSetStoppingTestFunction()`
E*/
typedef enum { PEP_STOP_BASIC,
               PEP_STOP_USER } PEPStop;

/*MC
   PEP_STOP_BASIC - The default stopping test.

   Note:
   By default, the termination of the outer loop is decided by calling
   `PEPStoppingBasic()`, which will stop if all requested eigenvalues are converged,
   or if the maximum number of iterations has been reached.

   Level: advanced

.seealso: [](ch:pep), `PEPStop`, `PEPSetStoppingTest()`, `PEPStoppingBasic()`
M*/

/*MC
   PEP_STOP_USER - The user-provided stopping test.

   Note:
   Customized stopping test using the user-provided function given with
   `PEPSetStoppingTestFunction()`.

   Level: advanced

.seealso: [](ch:pep), `PEPStop`, `PEPSetStoppingTest()`, `PEPSetStoppingTestFunction()`
M*/

/*E
   PEPConvergedReason - Reason a polynomial eigensolver was determined to have converged
   or diverged.

   Values:
+  `PEP_CONVERGED_TOL`          - converged up to tolerance
.  `PEP_CONVERGED_USER`         - converged due to a user-defined condition
.  `PEP_DIVERGED_ITS`           - exceeded the maximum number of allowed iterations
.  `PEP_DIVERGED_BREAKDOWN`     - generic breakdown in method
.  `PEP_DIVERGED_SYMMETRY_LOST` - pseudo-Lanczos was not able to keep symmetry
-  `PEP_CONVERGED_ITERATING`    - the solver is still running

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPGetConvergedReason()`, `PEPSetTolerances()`
E*/
typedef enum {/* converged */
              PEP_CONVERGED_TOL                =  1,
              PEP_CONVERGED_USER               =  2,
              /* diverged */
              PEP_DIVERGED_ITS                 = -1,
              PEP_DIVERGED_BREAKDOWN           = -2,
              PEP_DIVERGED_SYMMETRY_LOST       = -3,
              PEP_CONVERGED_ITERATING          =  0} PEPConvergedReason;
SLEPC_EXTERN const char *const*PEPConvergedReasons;

/*MC
   PEP_CONVERGED_TOL - The computed error estimates, based on residual norms,
   for all requested eigenvalues are below the tolerance.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPGetConvergedReason()`, `PEPConvergedReason`
M*/

/*MC
   PEP_CONVERGED_USER - The solver was declared converged due to a user-defined condition.

   Note:
   This happens only when a user-defined stopping test has been set with
   `PEPSetStoppingTestFunction()`.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPGetConvergedReason()`, `PEPConvergedReason`, `PEPSetStoppingTestFunction()`
M*/

/*MC
   PEP_DIVERGED_ITS - Exceeded the maximum number of allowed iterations
   before the convergence criterion was satisfied.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPGetConvergedReason()`, `PEPConvergedReason`
M*/

/*MC
   PEP_DIVERGED_BREAKDOWN - A breakdown in the solver was detected so the
   method could not continue.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPGetConvergedReason()`, `PEPConvergedReason`
M*/

/*MC
   PEP_DIVERGED_SYMMETRY_LOST - The selected solver uses a pseudo-Lanczos recurrence,
   which is numerically unstable, and a symmetry test revealed that instability
   had appeared so the solver could not continue.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPGetConvergedReason()`, `PEPConvergedReason`
M*/

/*MC
   PEP_CONVERGED_ITERATING - This value is returned if `PEPGetConvergedReason()` is called
   while `PEPSolve()` is still running.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPGetConvergedReason()`, `PEPConvergedReason`
M*/

SLEPC_EXTERN PetscErrorCode PEPCreate(MPI_Comm,PEP*);
SLEPC_EXTERN PetscErrorCode PEPDestroy(PEP*);
SLEPC_EXTERN PetscErrorCode PEPReset(PEP);
SLEPC_EXTERN PetscErrorCode PEPSetType(PEP,PEPType);
SLEPC_EXTERN PetscErrorCode PEPGetType(PEP,PEPType*);
SLEPC_EXTERN PetscErrorCode PEPSetProblemType(PEP,PEPProblemType);
SLEPC_EXTERN PetscErrorCode PEPGetProblemType(PEP,PEPProblemType*);
SLEPC_EXTERN PetscErrorCode PEPSetOperators(PEP,PetscInt,Mat[]);
SLEPC_EXTERN PetscErrorCode PEPGetOperators(PEP,PetscInt,Mat*);
SLEPC_EXTERN PetscErrorCode PEPGetNumMatrices(PEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode PEPSetTarget(PEP,PetscScalar);
SLEPC_EXTERN PetscErrorCode PEPGetTarget(PEP,PetscScalar*);
SLEPC_EXTERN PetscErrorCode PEPSetInterval(PEP,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPGetInterval(PEP,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPSetFromOptions(PEP);
SLEPC_EXTERN PetscErrorCode PEPSetDSType(PEP);
SLEPC_EXTERN PetscErrorCode PEPSetUp(PEP);
SLEPC_EXTERN PetscErrorCode PEPSolve(PEP);
SLEPC_EXTERN PetscErrorCode PEPView(PEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode PEPViewFromOptions(PEP,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode PEPErrorView(PEP,PEPErrorType,PetscViewer);
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "PEPErrorView()", ) static inline PetscErrorCode PEPPrintSolution(PEP pep,PetscViewer v) {return PEPErrorView(pep,PEP_ERROR_BACKWARD,v);}
SLEPC_EXTERN PetscErrorCode PEPErrorViewFromOptions(PEP);
SLEPC_EXTERN PetscErrorCode PEPConvergedReasonView(PEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode PEPConvergedReasonViewFromOptions(PEP);
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "PEPConvergedReasonView()", ) static inline PetscErrorCode PEPReasonView(PEP pep,PetscViewer v) {return PEPConvergedReasonView(pep,v);}
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "PEPConvergedReasonViewFromOptions()", ) static inline PetscErrorCode PEPReasonViewFromOptions(PEP pep) {return PEPConvergedReasonViewFromOptions(pep);}
SLEPC_EXTERN PetscErrorCode PEPValuesView(PEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode PEPValuesViewFromOptions(PEP);
SLEPC_EXTERN PetscErrorCode PEPVectorsView(PEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode PEPVectorsViewFromOptions(PEP);
SLEPC_EXTERN PetscErrorCode PEPSetBV(PEP,BV);
SLEPC_EXTERN PetscErrorCode PEPGetBV(PEP,BV*);
SLEPC_EXTERN PetscErrorCode PEPSetRG(PEP,RG);
SLEPC_EXTERN PetscErrorCode PEPGetRG(PEP,RG*);
SLEPC_EXTERN PetscErrorCode PEPSetDS(PEP,DS);
SLEPC_EXTERN PetscErrorCode PEPGetDS(PEP,DS*);
SLEPC_EXTERN PetscErrorCode PEPSetST(PEP,ST);
SLEPC_EXTERN PetscErrorCode PEPGetST(PEP,ST*);
SLEPC_EXTERN PetscErrorCode PEPRefineGetKSP(PEP,KSP*);

SLEPC_EXTERN PetscErrorCode PEPSetTolerances(PEP,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode PEPGetTolerances(PEP,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode PEPGetConvergedReason(PEP,PEPConvergedReason*);

SLEPC_EXTERN PetscErrorCode PEPSetDimensions(PEP,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode PEPGetDimensions(PEP,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode PEPSetScale(PEP,PEPScale,PetscReal,Vec,Vec,PetscInt,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPGetScale(PEP,PEPScale*,PetscReal*,Vec*,Vec*,PetscInt*,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPSetRefine(PEP,PEPRefine,PetscInt,PetscReal,PetscInt,PEPRefineScheme);
SLEPC_EXTERN PetscErrorCode PEPGetRefine(PEP,PEPRefine*,PetscInt*,PetscReal*,PetscInt*,PEPRefineScheme*);
SLEPC_EXTERN PetscErrorCode PEPSetExtract(PEP,PEPExtract);
SLEPC_EXTERN PetscErrorCode PEPGetExtract(PEP,PEPExtract*);
SLEPC_EXTERN PetscErrorCode PEPSetBasis(PEP,PEPBasis);
SLEPC_EXTERN PetscErrorCode PEPGetBasis(PEP,PEPBasis*);

SLEPC_EXTERN PetscErrorCode PEPGetConverged(PEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode PEPGetEigenpair(PEP,PetscInt,PetscScalar*,PetscScalar*,Vec,Vec);
SLEPC_EXTERN PetscErrorCode PEPComputeError(PEP,PetscInt,PEPErrorType,PetscReal*);
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "PEPComputeError()", ) static inline PetscErrorCode PEPComputeRelativeError(PEP pep,PetscInt i,PetscReal *r) {return PEPComputeError(pep,i,PEP_ERROR_BACKWARD,r);}
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "PEPComputeError() with PEP_ERROR_ABSOLUTE", ) static inline PetscErrorCode PEPComputeResidualNorm(PEP pep,PetscInt i,PetscReal *r) {return PEPComputeError(pep,i,PEP_ERROR_ABSOLUTE,r);}
SLEPC_EXTERN PetscErrorCode PEPGetErrorEstimate(PEP,PetscInt,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPGetIterationNumber(PEP,PetscInt*);

SLEPC_EXTERN PetscErrorCode PEPSetInitialSpace(PEP,PetscInt,Vec[]);
SLEPC_EXTERN PetscErrorCode PEPSetWhichEigenpairs(PEP,PEPWhich);
SLEPC_EXTERN PetscErrorCode PEPGetWhichEigenpairs(PEP,PEPWhich*);

SLEPC_EXTERN PetscErrorCode PEPSetTrackAll(PEP,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPGetTrackAll(PEP,PetscBool*);

/*S
   PEPMonitorFn - A function prototype for functions provided to `PEPMonitorSet()`.

   Calling Sequence:
+  pep    - the polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - relative error estimates for each eigenpair
.  nest   - number of error estimates
-  ctx    - optional monitoring context, as provided with `PEPMonitorSet()`

   Level: intermediate

.seealso: [](ch:pep), `PEPMonitorSet()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode PEPMonitorFn(PEP pep,PetscInt its,PetscInt nconv,PetscScalar eigr[],PetscScalar eigi[],PetscReal errest[],PetscInt nest,void *ctx);

/*S
   PEPMonitorRegisterFn - A function prototype for functions provided to `PEPMonitorRegister()`.

   Calling Sequence:
+  pep    - the polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - relative error estimates for each eigenpair
.  nest   - number of error estimates
-  ctx    - `PetscViewerAndFormat` object

   Level: advanced

   Note:
   This is a `PEPMonitorFn` specialized for a context of `PetscViewerAndFormat`.

.seealso: [](ch:pep), `PEPMonitorSet()`, `PEPMonitorRegister()`, `PEPMonitorFn`, `PEPMonitorRegisterCreateFn`, `PEPMonitorRegisterDestroyFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode PEPMonitorRegisterFn(PEP pep,PetscInt its,PetscInt nconv,PetscScalar eigr[],PetscScalar eigi[],PetscReal errest[],PetscInt nest,PetscViewerAndFormat *ctx);

/*S
   PEPMonitorRegisterCreateFn - A function prototype for functions that do the
   creation when provided to `PEPMonitorRegister()`.

   Calling Sequence:
+  viewer - the viewer to be used with the `PEPMonitorRegisterFn`
.  format - the format of the viewer
.  ctx    - a context for the monitor
-  result - a `PetscViewerAndFormat` object

   Level: advanced

.seealso: [](ch:pep), `PEPMonitorRegisterFn`, `PEPMonitorSet()`, `PEPMonitorRegister()`, `PEPMonitorFn`, `PEPMonitorRegisterDestroyFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode PEPMonitorRegisterCreateFn(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **result);

/*S
   PEPMonitorRegisterDestroyFn - A function prototype for functions that do the after
   use destruction when provided to `PEPMonitorRegister()`.

   Calling Sequence:
.  vf - a `PetscViewerAndFormat` object to be destroyed, including any context

   Level: advanced

.seealso: [](ch:pep), `PEPMonitorRegisterFn`, `PEPMonitorSet()`, `PEPMonitorRegister()`, `PEPMonitorFn`, `PEPMonitorRegisterCreateFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode PEPMonitorRegisterDestroyFn(PetscViewerAndFormat **result);

SLEPC_EXTERN PetscErrorCode PEPMonitor(PEP,PetscInt,PetscInt,PetscScalar[],PetscScalar[],PetscReal[],PetscInt);
SLEPC_EXTERN PetscErrorCode PEPMonitorSet(PEP,PEPMonitorFn,void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode PEPMonitorCancel(PEP);
SLEPC_EXTERN PetscErrorCode PEPGetMonitorContext(PEP,void*);

SLEPC_EXTERN PetscErrorCode PEPMonitorSetFromOptions(PEP,const char[],const char[],void*,PetscBool);
SLEPC_EXTERN PEPMonitorRegisterFn        PEPMonitorFirst;
SLEPC_EXTERN PEPMonitorRegisterFn        PEPMonitorFirstDrawLG;
SLEPC_EXTERN PEPMonitorRegisterCreateFn  PEPMonitorFirstDrawLGCreate;
SLEPC_EXTERN PEPMonitorRegisterFn        PEPMonitorAll;
SLEPC_EXTERN PEPMonitorRegisterFn        PEPMonitorAllDrawLG;
SLEPC_EXTERN PEPMonitorRegisterCreateFn  PEPMonitorAllDrawLGCreate;
SLEPC_EXTERN PEPMonitorRegisterFn        PEPMonitorConverged;
SLEPC_EXTERN PEPMonitorRegisterCreateFn  PEPMonitorConvergedCreate;
SLEPC_EXTERN PEPMonitorRegisterFn        PEPMonitorConvergedDrawLG;
SLEPC_EXTERN PEPMonitorRegisterCreateFn  PEPMonitorConvergedDrawLGCreate;
SLEPC_EXTERN PEPMonitorRegisterDestroyFn PEPMonitorConvergedDestroy;

SLEPC_EXTERN PetscErrorCode PEPSetOptionsPrefix(PEP,const char[]);
SLEPC_EXTERN PetscErrorCode PEPAppendOptionsPrefix(PEP,const char[]);
SLEPC_EXTERN PetscErrorCode PEPGetOptionsPrefix(PEP,const char*[]);

SLEPC_EXTERN PetscFunctionList PEPList;
SLEPC_EXTERN PetscFunctionList PEPMonitorList;
SLEPC_EXTERN PetscFunctionList PEPMonitorCreateList;
SLEPC_EXTERN PetscFunctionList PEPMonitorDestroyList;
SLEPC_EXTERN PetscErrorCode PEPRegister(const char[],PetscErrorCode(*)(PEP));
SLEPC_EXTERN PetscErrorCode PEPMonitorRegister(const char[],PetscViewerType,PetscViewerFormat,PEPMonitorRegisterFn*,PEPMonitorRegisterCreateFn*,PEPMonitorRegisterDestroyFn*);

SLEPC_EXTERN PetscErrorCode PEPSetWorkVecs(PEP,PetscInt);
SLEPC_EXTERN PetscErrorCode PEPAllocateSolution(PEP,PetscInt);

/*S
   PEPConvergenceTestFn - A prototype of a `PEP` convergence test function that
   would be passed to `PEPSetConvergenceTestFunction()`.

   Calling Sequence:
+  pep    - the polynomial eigensolver context
.  eigr   - real part of the eigenvalue
.  eigi   - imaginary part of the eigenvalue
.  res    - residual norm associated to the eigenpair
.  errest - [output] computed error estimate
-  ctx    - optional convergence context, as set by `PEPSetConvergenceTestFunction()`

   Level: advanced

.seealso: [](ch:pep), `PEPSetConvergenceTest()`, `PEPSetConvergenceTestFunction()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode PEPConvergenceTestFn(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx);

SLEPC_EXTERN PetscErrorCode PEPSetConvergenceTest(PEP,PEPConv);
SLEPC_EXTERN PetscErrorCode PEPGetConvergenceTest(PEP,PEPConv*);
SLEPC_EXTERN PEPConvergenceTestFn PEPConvergedAbsolute;
SLEPC_EXTERN PEPConvergenceTestFn PEPConvergedRelative;
SLEPC_EXTERN PEPConvergenceTestFn PEPConvergedNorm;
SLEPC_EXTERN PetscErrorCode PEPSetConvergenceTestFunction(PEP,PEPConvergenceTestFn*,void*,PetscCtxDestroyFn*);

/*S
   PEPStoppingTestFn - A prototype of a `PEP` stopping test function that would be
   passed to `PEPSetStoppingTestFunction()`.

   Calling Sequence:
+  pep    - the polynomial eigensolver context
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged eigenpairs
.  nev    - number of requested eigenpairs
.  reason - [output] result of the stopping test
-  ctx    - optional stopping context, as set by `PEPSetStoppingTestFunction()`

   Note:
   A positive value of `reason` indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, `reason` must be set to `PEP_CONVERGED_ITERATING`
   (zero).

   Level: advanced

.seealso: [](ch:pep), `PEPSetStoppingTestFunction()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode PEPStoppingTestFn(PEP pep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,PEPConvergedReason *reason,void *ctx);

SLEPC_EXTERN PetscErrorCode PEPSetStoppingTest(PEP,PEPStop);
SLEPC_EXTERN PetscErrorCode PEPGetStoppingTest(PEP,PEPStop*);
SLEPC_EXTERN PEPStoppingTestFn PEPStoppingBasic;
SLEPC_EXTERN PetscErrorCode PEPSetStoppingTestFunction(PEP,PEPStoppingTestFn*,void*,PetscCtxDestroyFn*);

SLEPC_EXTERN PetscErrorCode PEPSetEigenvalueComparison(PEP,SlepcEigenvalueComparisonFn*,void*);

/* --------- options specific to particular eigensolvers -------- */

SLEPC_EXTERN PetscErrorCode PEPLinearSetLinearization(PEP,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPLinearGetLinearization(PEP,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPLinearSetExplicitMatrix(PEP,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPLinearGetExplicitMatrix(PEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode PEPLinearSetEPS(PEP,EPS);
SLEPC_EXTERN PetscErrorCode PEPLinearGetEPS(PEP,EPS*);
PETSC_DEPRECATED_FUNCTION(3, 10, 0, "PEPLinearSetLinearization()", ) static inline PetscErrorCode PEPLinearSetCompanionForm(PEP pep,PetscInt cform) {return (cform==1)?PEPLinearSetLinearization(pep,1.0,0.0):PEPLinearSetLinearization(pep,0.0,1.0);}
PETSC_DEPRECATED_FUNCTION(3, 10, 0, "PEPLinearGetLinearization()", ) static inline PetscErrorCode PEPLinearGetCompanionForm(PEP pep,PetscInt *cform) {(void)pep; if (cform) *cform=1; return PETSC_SUCCESS;}

SLEPC_EXTERN PetscErrorCode PEPQArnoldiSetRestart(PEP,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPQArnoldiGetRestart(PEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPQArnoldiSetLocking(PEP,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPQArnoldiGetLocking(PEP,PetscBool*);

SLEPC_EXTERN PetscErrorCode PEPTOARSetRestart(PEP,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPTOARGetRestart(PEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPTOARSetLocking(PEP,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPTOARGetLocking(PEP,PetscBool*);

SLEPC_EXTERN PetscErrorCode PEPSTOARSetLinearization(PEP,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPSTOARGetLinearization(PEP,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPSTOARSetLocking(PEP,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPSTOARGetLocking(PEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode PEPSTOARSetDetectZeros(PEP,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPSTOARGetDetectZeros(PEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode PEPSTOARGetInertias(PEP,PetscInt*,PetscReal*[],PetscInt*[]);
SLEPC_EXTERN PetscErrorCode PEPSTOARSetDimensions(PEP,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode PEPSTOARGetDimensions(PEP,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode PEPSTOARSetCheckEigenvalueType(PEP,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPSTOARGetCheckEigenvalueType(PEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode PEPCheckDefiniteQEP(PEP,PetscReal*,PetscReal*,PetscInt*,PetscInt*);

/*E
   PEPJDProjection - The type of projection to be used in the Jacobi-Davidson solver.

   Values:
+  `PEP_JD_PROJECTION_HARMONIC`   - oblique projection
-  `PEP_JD_PROJECTION_ORTHOGONAL` - orthogonal projection

   Level: advanced

.seealso: [](ch:pep), `PEPJDSetProjection()`
E*/
typedef enum { PEP_JD_PROJECTION_HARMONIC,
               PEP_JD_PROJECTION_ORTHOGONAL } PEPJDProjection;
SLEPC_EXTERN const char *PEPJDProjectionTypes[];

SLEPC_EXTERN PetscErrorCode PEPJDSetRestart(PEP,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPJDGetRestart(PEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPJDSetFix(PEP,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPJDGetFix(PEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPJDSetReusePreconditioner(PEP,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPJDGetReusePreconditioner(PEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode PEPJDSetMinimalityIndex(PEP,PetscInt);
SLEPC_EXTERN PetscErrorCode PEPJDGetMinimalityIndex(PEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode PEPJDSetProjection(PEP,PEPJDProjection);
SLEPC_EXTERN PetscErrorCode PEPJDGetProjection(PEP,PEPJDProjection*);

/*E
   PEPCISSExtraction - The extraction technique used in the CISS solver.

   Values:
+  `PEP_CISS_EXTRACTION_RITZ`   - Rayleigh-Ritz extraction
.  `PEP_CISS_EXTRACTION_HANKEL` - block Hankel method
-  `PEP_CISS_EXTRACTION_CAA`    - communication-avoiding Arnoldi method

   Level: advanced

.seealso: [](ch:pep), `PEPCISSSetExtraction()`, `PEPCISSGetExtraction()`
E*/
typedef enum { PEP_CISS_EXTRACTION_RITZ,
               PEP_CISS_EXTRACTION_HANKEL,
               PEP_CISS_EXTRACTION_CAA    } PEPCISSExtraction;
SLEPC_EXTERN const char *PEPCISSExtractions[];

#if defined(PETSC_USE_COMPLEX) || defined(PETSC_CLANG_STATIC_ANALYZER)
SLEPC_EXTERN PetscErrorCode PEPCISSSetExtraction(PEP,PEPCISSExtraction);
SLEPC_EXTERN PetscErrorCode PEPCISSGetExtraction(PEP,PEPCISSExtraction*);
SLEPC_EXTERN PetscErrorCode PEPCISSSetSizes(PEP,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPCISSGetSizes(PEP,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*);
SLEPC_EXTERN PetscErrorCode PEPCISSSetThreshold(PEP,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode PEPCISSGetThreshold(PEP,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode PEPCISSSetRefinement(PEP,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode PEPCISSGetRefinement(PEP,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode PEPCISSGetKSPs(PEP,PetscInt*,KSP*[]);
#else
#define SlepcPEPCISSUnavailable(pep) do { \
    PetscFunctionBegin; \
    SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"%s() not available with real scalars",PETSC_FUNCTION_NAME); \
    } while (0)
static inline PetscErrorCode PEPCISSSetExtraction(PEP pep,PETSC_UNUSED PEPCISSExtraction ex) {SlepcPEPCISSUnavailable(pep);}
static inline PetscErrorCode PEPCISSGetExtraction(PEP pep,PETSC_UNUSED PEPCISSExtraction *ex) {SlepcPEPCISSUnavailable(pep);}
static inline PetscErrorCode PEPCISSSetSizes(PEP pep,PETSC_UNUSED PetscInt ip,PETSC_UNUSED PetscInt bs,PETSC_UNUSED PetscInt ms,PETSC_UNUSED PetscInt npart,PETSC_UNUSED PetscInt bsmax,PETSC_UNUSED PetscBool realmats) {SlepcPEPCISSUnavailable(pep);}
static inline PetscErrorCode PEPCISSGetSizes(PEP pep,PETSC_UNUSED PetscInt *ip,PETSC_UNUSED PetscInt *bs,PETSC_UNUSED PetscInt *ms,PETSC_UNUSED PetscInt *npart,PETSC_UNUSED PetscInt *bsmak,PETSC_UNUSED PetscBool *realmats) {SlepcPEPCISSUnavailable(pep);}
static inline PetscErrorCode PEPCISSSetThreshold(PEP pep,PETSC_UNUSED PetscReal delta,PETSC_UNUSED PetscReal spur) {SlepcPEPCISSUnavailable(pep);}
static inline PetscErrorCode PEPCISSGetThreshold(PEP pep,PETSC_UNUSED PetscReal *delta,PETSC_UNUSED PetscReal *spur) {SlepcPEPCISSUnavailable(pep);}
static inline PetscErrorCode PEPCISSSetRefinement(PEP pep,PETSC_UNUSED PetscInt inner,PETSC_UNUSED PetscInt blsize) {SlepcPEPCISSUnavailable(pep);}
static inline PetscErrorCode PEPCISSGetRefinement(PEP pep,PETSC_UNUSED PetscInt *inner,PETSC_UNUSED PetscInt *blsize) {SlepcPEPCISSUnavailable(pep);}
static inline PetscErrorCode PEPCISSGetKSPs(PEP pep,PETSC_UNUSED PetscInt *nsolve,PETSC_UNUSED KSP *ksp[]) {SlepcPEPCISSUnavailable(pep);}
#undef SlepcPEPCISSUnavailable
#endif
