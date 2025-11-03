/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the SLEPc linear eigenvalue solvers
*/

#pragma once

#include <slepcst.h>
#include <slepcbv.h>
#include <slepcds.h>
#include <slepcrg.h>
#include <slepclme.h>
#include <petscsnes.h>

/* SUBMANSEC = EPS */

SLEPC_EXTERN PetscErrorCode EPSInitializePackage(void);
SLEPC_EXTERN PetscErrorCode EPSFinalizePackage(void);

/*S
   EPS - SLEPc object that manages all the linear eigenvalue problem solvers.

   Level: beginner

.seealso: [](ch:eps), `EPSCreate()`, `ST`
S*/
typedef struct _p_EPS* EPS;

/*J
   EPSType - String with the name of a linear eigensolver.

   Level: beginner

.seealso: [](ch:eps), `EPSSetType()`, `EPS`
J*/
typedef const char *EPSType;
#define EPSPOWER       "power"
#define EPSSUBSPACE    "subspace"
#define EPSARNOLDI     "arnoldi"
#define EPSLANCZOS     "lanczos"
#define EPSKRYLOVSCHUR "krylovschur"
#define EPSGD          "gd"
#define EPSJD          "jd"
#define EPSRQCG        "rqcg"
#define EPSLOBPCG      "lobpcg"
#define EPSCISS        "ciss"
#define EPSLYAPII      "lyapii"
#define EPSLAPACK      "lapack"
#define EPSARPACK      "arpack"
#define EPSBLOPEX      "blopex"
#define EPSPRIMME      "primme"
#define EPSFEAST       "feast"
#define EPSSCALAPACK   "scalapack"
#define EPSELPA        "elpa"
#define EPSELEMENTAL   "elemental"
#define EPSEVSL        "evsl"
#define EPSCHASE       "chase"

/* Logging support */
SLEPC_EXTERN PetscClassId EPS_CLASSID;

/*E
   EPSProblemType - Determines the type of eigenvalue problem.

   Values:
+  `EPS_HEP`    - Hermitian
.  `EPS_GHEP`   - generalized Hermitian
.  `EPS_NHEP`   - non-Hermitian
.  `EPS_GNHEP`  - generalized non-Hermitian
.  `EPS_PGNHEP` - generalized non-Hermitian with positive (semi-)definite $B$
.  `EPS_GHIEP`  - generalized Hermitian-indefinite
.  `EPS_BSE`    - structured Bethe-Salpeter
-  `EPS_HAMILT` - structured Hamiltonian

   Note:
   In real scalars, one should read the term Hermitian as symmetric.

   Level: intermediate

.seealso: [](ch:eps), `EPSSetProblemType()`, `EPSGetProblemType()`
E*/
typedef enum { EPS_HEP    = 1,
               EPS_GHEP   = 2,
               EPS_NHEP   = 3,
               EPS_GNHEP  = 4,
               EPS_PGNHEP = 5,
               EPS_GHIEP  = 6,
               EPS_BSE    = 7,
               EPS_HAMILT = 8 } EPSProblemType;

/*MC
   EPS_HEP - A Hermitian eigenvalue problem.

   Note:
   The problem is formulated as $Ax=\lambda x$, where $A$ is real symmetric
   or complex Hermitian.

   Level: intermediate

.seealso: [](ch:eps), `EPSProblemType`, `EPSSetProblemType()`
M*/

/*MC
   EPS_GHEP - A generalized Hermitian eigenvalue problem.

   Note:
   The problem is formulated as $Ax=\lambda Bx$, where $A$ and $B$ are real
   symmetric or complex Hermitian, and $B$ is positive (semi-)definite.

   Level: intermediate

.seealso: [](ch:eps), `EPSProblemType`, `EPSSetProblemType()`
M*/

/*MC
   EPS_NHEP - A non-Hermitian eigenvalue problem.

   Note:
   The problem is formulated as $Ax=\lambda x$, where $A$ is non-symmetric
   (or non-Hermitian).

   Level: intermediate

.seealso: [](ch:eps), `EPSProblemType`, `EPSSetProblemType()`
M*/

/*MC
   EPS_GNHEP - A generalized non-Hermitian eigenvalue problem.

   Note:
   The problem is formulated as $Ax=\lambda Bx$, where $A$ or $B$ are
   non-symmetric (or non-Hermitian).

   Level: intermediate

.seealso: [](ch:eps), `EPSProblemType`, `EPSSetProblemType()`
M*/

/*MC
   EPS_PGNHEP - A generalized non-Hermitian eigenvalue problem with positive
   (semi-)definite $B$.

   Notes:
   The problem is formulated as $Ax=\lambda Bx$, where $A$ is non-symmetric
   (or non-Hermitian), but $B$ is symmetric (or Hermitian) and positive
   (semi-)definite.

   The problem will be solved with a non-Hermitian solver, but using an
   inner product induced by matrix $B$.

   Level: intermediate

.seealso: [](ch:eps), `EPSProblemType`, `EPSSetProblemType()`
M*/

/*MC
   EPS_GHIEP - A generalized Hermitian-indefinite eigenvalue problem.

   Notes:
   The problem is formulated as $Ax=\lambda Bx$, where both $A$ and $B$ are
   real symmetric or complex Hermitian, but $B$ is indefinite.

   The solver will try to exploit the symmetry by using an indefinite
   inner product, which may turn the computation numerically unstable.
   To avoid this, solve the problem as non-Hermitian.

   Level: intermediate

.seealso: [](ch:eps), `EPSProblemType`, `EPSSetProblemType()`
M*/

/*MC
   EPS_BSE - A structured Bethe-Salpeter eigenvalue problem.

   Notes:
   The problem is formulated as $Hx=\lambda x$, where $H$ has a Bethe-Salpeter
   structure,
     $$H = \begin{bmatrix}
        R & C \\
        -C^* & -R^T
        \end{bmatrix},$$
   where $R$ is Hermitian and $C$ is complex symmetric. Can also be used in
   the case of real matrices.

   A description of the properties of this problem can be found in {cite:p}`Alv25`
   and references therein.

   Level: intermediate

.seealso: [](ch:eps), [](sec:structured), `EPSProblemType`, `EPSSetProblemType()`
M*/

/*MC
   EPS_HAMILT - A structured Hamiltonian eigenvalue problem.

   Note:
   The problem is formulated as $Hx=\lambda x$, where $H$ has a Hamiltonian
   structure,
     $$H = \begin{bmatrix}
        A & B \\
        C & -A^*
        \end{bmatrix},$$
   where $A$, $B$ and $C$ are either real with $B=B^T$, $C=C^T$, or complex with
   $B=B^*$, $C=C^*$.

   Level: intermediate

.seealso: [](ch:eps), [](sec:structured), `EPSProblemType`, `EPSSetProblemType()`
M*/

/*E
   EPSExtraction - Determines the type of extraction technique employed
   by the eigensolver.

   Values:
+  `EPS_RITZ`              - Rayleigh-Ritz extraction
.  `EPS_HARMONIC`          - harmonic Ritz extraction
.  `EPS_HARMONIC_RELATIVE` - harmonic Ritz extraction relative to the eigenvalue
.  `EPS_HARMONIC_RIGHT`    - harmonic Ritz extraction for rightmost eigenvalues
.  `EPS_HARMONIC_LARGEST`  - harmonic Ritz extraction for largest magnitude (without target)
.  `EPS_REFINED`           - refined Ritz extraction
-  `EPS_REFINED_HARMONIC`  - refined harmonic Ritz extraction

   Level: advanced

.seealso: [](ch:eps), `EPSSetExtraction()`, `EPSGetExtraction()`
E*/
typedef enum { EPS_RITZ,
               EPS_HARMONIC,
               EPS_HARMONIC_RELATIVE,
               EPS_HARMONIC_RIGHT,
               EPS_HARMONIC_LARGEST,
               EPS_REFINED,
               EPS_REFINED_HARMONIC } EPSExtraction;

/*MC
   EPS_RITZ - The standard Rayleigh-Ritz extraction.

   Note:
   This is the default way of computing eigenpair approximations from a
   given subspace.

   Level: advanced

.seealso: [](ch:eps), `EPSExtraction`, `EPSSetExtraction()`
M*/

/*MC
   EPS_HARMONIC - The harmonic Ritz extraction.

   Notes:
   This extraction method may provide better convergence when computing
   interior eigenvalues close to a given target.

   For the particular case of Krylov-Schur, a detailed description can
   be found in {cite:p}`Rom09`.

   Level: advanced

.seealso: [](ch:eps), `EPSExtraction`, `EPSSetExtraction()`, `EPSSetTarget()`
M*/

/*MC
   EPS_HARMONIC_RELATIVE - The harmonic Ritz extraction relative to the eigenvalue.

   Note:
   This is a variation of `EPS_HARMONIC`, used in Davidson methods only.

   Level: advanced

.seealso: [](ch:eps), `EPSExtraction`, `EPSSetExtraction()`, `EPSSetTarget()`
M*/

/*MC
   EPS_HARMONIC_RIGHT - The harmonic Ritz extraction for rightmost eigenvalues.

   Note:
   This is a variation of `EPS_HARMONIC`, used in Davidson methods only.

   Level: advanced

.seealso: [](ch:eps), `EPSExtraction`, `EPSSetExtraction()`, `EPSSetTarget()`
M*/

/*MC
   EPS_HARMONIC_LARGEST - The harmonic Ritz extraction for largest magnitude
   eigenvalues (without target).

   Note:
   This is a variation of `EPS_HARMONIC`, used in Davidson methods only.

   Level: advanced

.seealso: [](ch:eps), `EPSExtraction`, `EPSSetExtraction()`
M*/

/*MC
   EPS_REFINED - The refined Ritz extraction method {cite:p}`Jia97`.

   Note:
   Currently implemented only in `EPSARNOLDI`.

   Level: advanced

.seealso: [](ch:eps), `EPSExtraction`, `EPSSetExtraction()`
M*/

/*MC
   EPS_REFINED_HARMONIC - The refined harmonic Ritz extraction.

   Note:
   This is a combination of `EPS_HARMONIC` and `EPS_REFINED`.

   Developer Note:
   Currently not implemented, reserved for future use.

   Level: advanced

.seealso: [](ch:eps), `EPSExtraction`, `EPSSetExtraction()`
M*/

/*E
   EPSWhich - Determines which part of the spectrum is requested.

   Values:
+  `EPS_LARGEST_MAGNITUDE`  - largest $|\lambda|$
.  `EPS_SMALLEST_MAGNITUDE` - smallest $|\lambda|$
.  `EPS_LARGEST_REAL`       - largest $\mathrm{Re}(\lambda)$
.  `EPS_SMALLEST_REAL`      - smallest $\mathrm{Re}(\lambda)$
.  `EPS_LARGEST_IMAGINARY`  - largest $\mathrm{Im}(\lambda)$
.  `EPS_SMALLEST_IMAGINARY` - smallest $\mathrm{Im}(\lambda)$
.  `EPS_TARGET_MAGNITUDE`   - smallest $|\lambda-\tau|$
.  `EPS_TARGET_REAL`        - smallest $|\mathrm{Re}(\lambda-\tau)|$
.  `EPS_TARGET_IMAGINARY`   - smallest $|\mathrm{Im}(\lambda-\tau)|$
.  `EPS_ALL`                - all $\lambda\in[a,b]$ or $\lambda\in\Omega$
-  `EPS_WHICH_USER`         - user-defined sorting criterion

   Notes:
   If SLEPc is compiled for real scalars `EPS_LARGEST_IMAGINARY` and
   `EPS_SMALLEST_IMAGINARY` use the absolute value of the imaginary part
   for eigenvalue selection.

   The target $\tau$ is a scalar value provided with `EPSSetTarget()`.

   The case `EPS_ALL` needs an interval $[a,b]$ given with `EPSSetInterval()`
   or a region $\Omega$ specified with an `RG` object.

   Level: intermediate

.seealso: [](ch:eps), `EPSSetWhichEigenpairs()`, `EPSSetTarget()`, `EPSSetInterval()`
E*/
typedef enum { EPS_LARGEST_MAGNITUDE  = 1,
               EPS_SMALLEST_MAGNITUDE = 2,
               EPS_LARGEST_REAL       = 3,
               EPS_SMALLEST_REAL      = 4,
               EPS_LARGEST_IMAGINARY  = 5,
               EPS_SMALLEST_IMAGINARY = 6,
               EPS_TARGET_MAGNITUDE   = 7,
               EPS_TARGET_REAL        = 8,
               EPS_TARGET_IMAGINARY   = 9,
               EPS_ALL                = 10,
               EPS_WHICH_USER         = 11 } EPSWhich;

/*E
   EPSBalance - The type of balancing used for non-Hermitian problems.

   Values:
+  `EPS_BALANCE_NONE`    - no balancing matrix is used
.  `EPS_BALANCE_ONESIDE` - balancing matrix $D$ is computed with a one-sided Krylov method
.  `EPS_BALANCE_TWOSIDE` - balancing matrix $D$ is computed with a two-sided Krylov method
-  `EPS_BALANCE_USER`    - use a balancing matrix $D$ provided by the user

   Level: intermediate

.seealso: [](ch:eps), [](sec:balancing), `EPSSetBalance()`
E*/
typedef enum { EPS_BALANCE_NONE,
               EPS_BALANCE_ONESIDE,
               EPS_BALANCE_TWOSIDE,
               EPS_BALANCE_USER } EPSBalance;
SLEPC_EXTERN const char *EPSBalanceTypes[];

/*E
   EPSErrorType - The error type used to assess the accuracy of computed solutions.

   Values:
+  `EPS_ERROR_ABSOLUTE` - compute error bound as $\|r\|$
.  `EPS_ERROR_RELATIVE` - compute error bound as $\|r\|/|\lambda|$
-  `EPS_ERROR_BACKWARD` - compute error bound as $\|r\|/(\|A\|+|\lambda|\|B\|)$

   Level: intermediate

.seealso: [](ch:eps), `EPSComputeError()`
E*/
typedef enum { EPS_ERROR_ABSOLUTE,
               EPS_ERROR_RELATIVE,
               EPS_ERROR_BACKWARD } EPSErrorType;
SLEPC_EXTERN const char *EPSErrorTypes[];

/*E
   EPSConv - The convergence criterion to be used by the solver.

   Values:
+  `EPS_CONV_ABS`  - absolute convergence criterion, $\|r\|$
.  `EPS_CONV_REL`  - convergence criterion relative to eigenvalue, $\|r\|/|\lambda|$
.  `EPS_CONV_NORM` - convergence criterion relative to matrix norms, $\|r\|/(\|A\|+|\lambda|\|B\|)$
-  `EPS_CONV_USER` - convergence dictated by user-provided function

   Level: intermediate

.seealso: [](ch:eps), `EPSSetConvergenceTest()`, `EPSSetConvergenceTestFunction()`
E*/
typedef enum { EPS_CONV_ABS,
               EPS_CONV_REL,
               EPS_CONV_NORM,
               EPS_CONV_USER } EPSConv;

/*E
   EPSStop - The stopping test to decide the termination of the outer loop
   of the eigensolver.

   Values:
+  `EPS_STOP_BASIC`     - default stopping test
.  `EPS_STOP_USER`      - user-provided stopping test
-  `EPS_STOP_THRESHOLD` - threshold stopping test

   Level: advanced

.seealso: [](ch:eps), `EPSSetStoppingTest()`, `EPSSetStoppingTestFunction()`
E*/
typedef enum { EPS_STOP_BASIC,
               EPS_STOP_USER,
               EPS_STOP_THRESHOLD } EPSStop;

/*MC
   EPS_STOP_BASIC - The default stopping test.

   Note:
   By default, the termination of the outer loop is decided by calling
   `EPSStoppingBasic()`, which will stop if all requested eigenvalues are converged,
   or if the maximum number of iterations has been reached.

   Level: advanced

.seealso: [](ch:eps), `EPSStop`, `EPSSetStoppingTest()`, `EPSStoppingBasic()`
M*/

/*MC
   EPS_STOP_USER - The user-provided stopping test.

   Note:
   Customized stopping test using the user-provided function given with
   `EPSSetStoppingTestFunction()`.

   Level: advanced

.seealso: [](ch:eps), `EPSStop`, `EPSSetStoppingTest()`, `EPSSetStoppingTestFunction()`
M*/

/*MC
   EPS_STOP_THRESHOLD - The threshold stopping test.

   Note:
   When a threshold has been provided with `EPSSetThreshold()`, the termination
   of the outer loop is decided by calling `EPSStoppingThreshold()`, which will
   stop when one of the computed eigenvalues is not above/below the threshold.
   If a number of wanted eigenvalues has been specified via `EPSSetDimensions()`
   then it is also taken into account, and the solver will stop when one of the
   two conditions (threshold or number of converged values) is met.

   Level: advanced

.seealso: [](ch:eps), `EPSStop`, `EPSSetStoppingTest()`, `EPSStoppingThreshold()`, `EPSSetThreshold()`, `EPSSetDimensions()`
M*/

/*E
   EPSConvergedReason - Reason an eigensolver was determined to have converged
   or diverged.

   Values:
+  `EPS_CONVERGED_TOL`          - converged up to tolerance
.  `EPS_CONVERGED_USER`         - converged due to a user-defined condition
.  `EPS_DIVERGED_ITS`           - exceeded the maximum number of allowed iterations
.  `EPS_DIVERGED_BREAKDOWN`     - generic breakdown in method
.  `EPS_DIVERGED_SYMMETRY_LOST` - pseudo-Lanczos was not able to keep symmetry
-  `EPS_CONVERGED_ITERATING`    - the solver is still running

   Level: intermediate

.seealso: [](ch:eps), `EPSSolve()`, `EPSGetConvergedReason()`, `EPSSetTolerances()`
E*/
typedef enum {/* converged */
              EPS_CONVERGED_TOL                =  1,
              EPS_CONVERGED_USER               =  2,
              /* diverged */
              EPS_DIVERGED_ITS                 = -1,
              EPS_DIVERGED_BREAKDOWN           = -2,
              EPS_DIVERGED_SYMMETRY_LOST       = -3,
              EPS_CONVERGED_ITERATING          =  0} EPSConvergedReason;
SLEPC_EXTERN const char *const*EPSConvergedReasons;

/*MC
   EPS_CONVERGED_TOL - The computed error estimates, based on residual norms,
   for all requested eigenvalues are below the tolerance.

   Level: intermediate

.seealso: [](ch:eps), `EPSSolve()`, `EPSGetConvergedReason()`, `EPSConvergedReason`
M*/

/*MC
   EPS_CONVERGED_USER - The solver was declared converged due to a user-defined condition.

   Note:
   This happens only when a user-defined stopping test has been set with
   `EPSSetStoppingTestFunction()`.

   Level: intermediate

.seealso: [](ch:eps), `EPSSolve()`, `EPSGetConvergedReason()`, `EPSConvergedReason`, `EPSSetStoppingTestFunction()`
M*/

/*MC
   EPS_DIVERGED_ITS - Exceeded the maximum number of allowed iterations
   before the convergence criterion was satisfied.

   Level: intermediate

.seealso: [](ch:eps), `EPSSolve()`, `EPSGetConvergedReason()`, `EPSConvergedReason`
M*/

/*MC
   EPS_DIVERGED_BREAKDOWN - A breakdown in the solver was detected so the
   method could not continue.

   Level: intermediate

.seealso: [](ch:eps), `EPSSolve()`, `EPSGetConvergedReason()`, `EPSConvergedReason`
M*/

/*MC
   EPS_DIVERGED_SYMMETRY_LOST - The selected solver uses a pseudo-Lanczos recurrence,
   which is numerically unstable, and a symmetry test revealed that instability
   had appeared so the solver could not continue.

   Level: intermediate

.seealso: [](ch:eps), `EPSSolve()`, `EPSGetConvergedReason()`, `EPSConvergedReason`
M*/

/*MC
   EPS_CONVERGED_ITERATING - This value is returned if `EPSGetConvergedReason()` is called
   while `EPSSolve()` is still running.

   Level: intermediate

.seealso: [](ch:eps), `EPSSolve()`, `EPSGetConvergedReason()`, `EPSConvergedReason`
M*/

/*S
   EPSStoppingCtx - Data structure (C struct) to hold additional information to
   be used in some stopping test functions.

   Level: advanced

.seealso: [](ch:eps), `EPSSetStoppingTestFunction()`
S*/
struct _n_EPSStoppingCtx {
  PetscReal firstev;    /* the (absolute) value of the first converged eigenvalue */
  PetscReal lastev;     /* the (absolute) value of the last converged eigenvalue */
  PetscReal thres;      /* threshold set with EPSSetThreshold() */
  PetscBool threlative; /* threshold is relative */
  EPSWhich  which;      /* which eigenvalues are being computed */
};
typedef struct _n_EPSStoppingCtx* EPSStoppingCtx;

SLEPC_EXTERN PetscErrorCode EPSCreate(MPI_Comm,EPS*);
SLEPC_EXTERN PetscErrorCode EPSDestroy(EPS*);
SLEPC_EXTERN PetscErrorCode EPSReset(EPS);
SLEPC_EXTERN PetscErrorCode EPSSetType(EPS,EPSType);
SLEPC_EXTERN PetscErrorCode EPSGetType(EPS,EPSType*);
SLEPC_EXTERN PetscErrorCode EPSSetProblemType(EPS,EPSProblemType);
SLEPC_EXTERN PetscErrorCode EPSGetProblemType(EPS,EPSProblemType*);
SLEPC_EXTERN PetscErrorCode EPSSetExtraction(EPS,EPSExtraction);
SLEPC_EXTERN PetscErrorCode EPSGetExtraction(EPS,EPSExtraction*);
SLEPC_EXTERN PetscErrorCode EPSSetBalance(EPS,EPSBalance,PetscInt,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSGetBalance(EPS,EPSBalance*,PetscInt*,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSSetOperators(EPS,Mat,Mat);
SLEPC_EXTERN PetscErrorCode EPSGetOperators(EPS,Mat*,Mat*);
SLEPC_EXTERN PetscErrorCode EPSSetFromOptions(EPS);
SLEPC_EXTERN PetscErrorCode EPSSetDSType(EPS);
SLEPC_EXTERN PetscErrorCode EPSSetUp(EPS);
SLEPC_EXTERN PetscErrorCode EPSSolve(EPS);
SLEPC_EXTERN PetscErrorCode EPSView(EPS,PetscViewer);
SLEPC_EXTERN PetscErrorCode EPSViewFromOptions(EPS,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode EPSErrorView(EPS,EPSErrorType,PetscViewer);
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "EPSErrorView()", ) static inline PetscErrorCode EPSPrintSolution(EPS eps,PetscViewer v) {return EPSErrorView(eps,EPS_ERROR_RELATIVE,v);}
SLEPC_EXTERN PetscErrorCode EPSErrorViewFromOptions(EPS);
SLEPC_EXTERN PetscErrorCode EPSConvergedReasonView(EPS,PetscViewer);
SLEPC_EXTERN PetscErrorCode EPSConvergedReasonViewFromOptions(EPS);
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "EPSConvergedReasonView()", ) static inline PetscErrorCode EPSReasonView(EPS eps,PetscViewer v) {return EPSConvergedReasonView(eps,v);}
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "EPSConvergedReasonViewFromOptions()", ) static inline PetscErrorCode EPSReasonViewFromOptions(EPS eps) {return EPSConvergedReasonViewFromOptions(eps);}
SLEPC_EXTERN PetscErrorCode EPSValuesView(EPS,PetscViewer);
SLEPC_EXTERN PetscErrorCode EPSValuesViewFromOptions(EPS);
SLEPC_EXTERN PetscErrorCode EPSVectorsView(EPS,PetscViewer);
SLEPC_EXTERN PetscErrorCode EPSVectorsViewFromOptions(EPS);

SLEPC_EXTERN PetscErrorCode EPSSetTarget(EPS,PetscScalar);
SLEPC_EXTERN PetscErrorCode EPSGetTarget(EPS,PetscScalar*);
SLEPC_EXTERN PetscErrorCode EPSSetInterval(EPS,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSGetInterval(EPS,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSSetST(EPS,ST);
SLEPC_EXTERN PetscErrorCode EPSGetST(EPS,ST*);
SLEPC_EXTERN PetscErrorCode EPSSetBV(EPS,BV);
SLEPC_EXTERN PetscErrorCode EPSGetBV(EPS,BV*);
SLEPC_EXTERN PetscErrorCode EPSSetRG(EPS,RG);
SLEPC_EXTERN PetscErrorCode EPSGetRG(EPS,RG*);
SLEPC_EXTERN PetscErrorCode EPSSetDS(EPS,DS);
SLEPC_EXTERN PetscErrorCode EPSGetDS(EPS,DS*);
SLEPC_EXTERN PetscErrorCode EPSSetTolerances(EPS,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGetTolerances(EPS,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSSetDimensions(EPS,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGetDimensions(EPS,PetscInt*,PetscInt*,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSGetConvergedReason(EPS,EPSConvergedReason*);

SLEPC_EXTERN PetscErrorCode EPSGetConverged(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSGetEigenpair(EPS,PetscInt,PetscScalar*,PetscScalar*,Vec,Vec);
SLEPC_EXTERN PetscErrorCode EPSGetEigenvalue(EPS,PetscInt,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode EPSGetEigenvector(EPS,PetscInt,Vec,Vec);
SLEPC_EXTERN PetscErrorCode EPSGetLeftEigenvector(EPS,PetscInt,Vec,Vec);

SLEPC_EXTERN PetscErrorCode EPSComputeError(EPS,PetscInt,EPSErrorType,PetscReal*);
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "EPSComputeError()", ) static inline PetscErrorCode EPSComputeRelativeError(EPS eps,PetscInt i,PetscReal *r) {return EPSComputeError(eps,i,EPS_ERROR_RELATIVE,r);}
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "EPSComputeError() with EPS_ERROR_ABSOLUTE", ) static inline PetscErrorCode EPSComputeResidualNorm(EPS eps,PetscInt i,PetscReal *r) {return EPSComputeError(eps,i,EPS_ERROR_ABSOLUTE,r);}
SLEPC_EXTERN PetscErrorCode EPSGetInvariantSubspace(EPS,Vec[]);
SLEPC_EXTERN PetscErrorCode EPSGetErrorEstimate(EPS,PetscInt,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSGetIterationNumber(EPS,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSSetWhichEigenpairs(EPS,EPSWhich);
SLEPC_EXTERN PetscErrorCode EPSGetWhichEigenpairs(EPS,EPSWhich*);
SLEPC_EXTERN PetscErrorCode EPSSetThreshold(EPS,PetscReal,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetThreshold(EPS,PetscReal*,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSSetTwoSided(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetTwoSided(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSSetTrueResidual(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetTrueResidual(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSSetPurify(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetPurify(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSIsGeneralized(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSIsHermitian(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSIsPositive(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSIsStructured(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSSetTrackAll(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetTrackAll(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSSetDeflationSpace(EPS,PetscInt,Vec[]);
SLEPC_EXTERN PetscErrorCode EPSSetInitialSpace(EPS,PetscInt,Vec[]);
SLEPC_EXTERN PetscErrorCode EPSSetLeftInitialSpace(EPS,PetscInt,Vec[]);

/*S
   EPSMonitorFn - A function prototype for functions provided to `EPSMonitorSet()`.

   Calling Sequence:
+  eps    - the linear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - relative error estimates for each eigenpair
.  nest   - number of error estimates
-  ctx    - optional monitoring context, as provided with `EPSMonitorSet()`

   Level: intermediate

.seealso: [](ch:eps), `EPSMonitorSet()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode EPSMonitorFn(EPS eps,PetscInt its,PetscInt nconv,PetscScalar eigr[],PetscScalar eigi[],PetscReal errest[],PetscInt nest,void *ctx);

/*S
   EPSMonitorRegisterFn - A function prototype for functions provided to `EPSMonitorRegister()`.

   Calling Sequence:
+  eps    - the linear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - relative error estimates for each eigenpair
.  nest   - number of error estimates
-  ctx    - `PetscViewerAndFormat` object

   Level: advanced

   Note:
   This is an `EPSMonitorFn` specialized for a context of `PetscViewerAndFormat`.

.seealso: [](ch:eps), `EPSMonitorSet()`, `EPSMonitorRegister()`, `EPSMonitorFn`, `EPSMonitorRegisterCreateFn`, `EPSMonitorRegisterDestroyFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode EPSMonitorRegisterFn(EPS eps,PetscInt its,PetscInt nconv,PetscScalar eigr[],PetscScalar eigi[],PetscReal errest[],PetscInt nest,PetscViewerAndFormat *ctx);

/*S
   EPSMonitorRegisterCreateFn - A function prototype for functions that do the
   creation when provided to `EPSMonitorRegister()`.

   Calling Sequence:
+  viewer - the viewer to be used with the `EPSMonitorRegisterFn`
.  format - the format of the viewer
.  ctx    - a context for the monitor
-  result - a `PetscViewerAndFormat` object

   Level: advanced

.seealso: [](ch:eps), `EPSMonitorRegisterFn`, `EPSMonitorSet()`, `EPSMonitorRegister()`, `EPSMonitorFn`, `EPSMonitorRegisterDestroyFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode EPSMonitorRegisterCreateFn(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **result);

/*S
   EPSMonitorRegisterDestroyFn - A function prototype for functions that do the after
   use destruction when provided to `EPSMonitorRegister()`.

   Calling Sequence:
.  vf - a `PetscViewerAndFormat` object to be destroyed, including any context

   Level: advanced

.seealso: [](ch:eps), `EPSMonitorRegisterFn`, `EPSMonitorSet()`, `EPSMonitorRegister()`, `EPSMonitorFn`, `EPSMonitorRegisterCreateFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode EPSMonitorRegisterDestroyFn(PetscViewerAndFormat **result);

SLEPC_EXTERN PetscErrorCode EPSMonitor(EPS,PetscInt,PetscInt,PetscScalar[],PetscScalar[],PetscReal[],PetscInt);
SLEPC_EXTERN PetscErrorCode EPSMonitorSet(EPS,EPSMonitorFn,void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode EPSMonitorCancel(EPS);
SLEPC_EXTERN PetscErrorCode EPSGetMonitorContext(EPS,void*);

SLEPC_EXTERN PetscErrorCode EPSMonitorSetFromOptions(EPS,const char[],const char[],void*,PetscBool);
SLEPC_EXTERN EPSMonitorRegisterFn        EPSMonitorFirst;
SLEPC_EXTERN EPSMonitorRegisterFn        EPSMonitorFirstDrawLG;
SLEPC_EXTERN EPSMonitorRegisterCreateFn  EPSMonitorFirstDrawLGCreate;
SLEPC_EXTERN EPSMonitorRegisterFn        EPSMonitorAll;
SLEPC_EXTERN EPSMonitorRegisterFn        EPSMonitorAllDrawLG;
SLEPC_EXTERN EPSMonitorRegisterCreateFn  EPSMonitorAllDrawLGCreate;
SLEPC_EXTERN EPSMonitorRegisterFn        EPSMonitorConverged;
SLEPC_EXTERN EPSMonitorRegisterCreateFn  EPSMonitorConvergedCreate;
SLEPC_EXTERN EPSMonitorRegisterFn        EPSMonitorConvergedDrawLG;
SLEPC_EXTERN EPSMonitorRegisterCreateFn  EPSMonitorConvergedDrawLGCreate;
SLEPC_EXTERN EPSMonitorRegisterDestroyFn EPSMonitorConvergedDestroy;

SLEPC_EXTERN PetscErrorCode EPSSetOptionsPrefix(EPS,const char[]);
SLEPC_EXTERN PetscErrorCode EPSAppendOptionsPrefix(EPS,const char[]);
SLEPC_EXTERN PetscErrorCode EPSGetOptionsPrefix(EPS,const char*[]);

SLEPC_EXTERN PetscFunctionList EPSList;
SLEPC_EXTERN PetscFunctionList EPSMonitorList;
SLEPC_EXTERN PetscFunctionList EPSMonitorCreateList;
SLEPC_EXTERN PetscFunctionList EPSMonitorDestroyList;
SLEPC_EXTERN PetscErrorCode EPSRegister(const char[],PetscErrorCode(*)(EPS));
SLEPC_EXTERN PetscErrorCode EPSMonitorRegister(const char[],PetscViewerType,PetscViewerFormat,EPSMonitorRegisterFn*,EPSMonitorRegisterCreateFn*,EPSMonitorRegisterDestroyFn*);

SLEPC_EXTERN PetscErrorCode EPSSetWorkVecs(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSAllocateSolution(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSReallocateSolution(EPS,PetscInt);

/*S
   EPSConvergenceTestFn - A prototype of an `EPS` convergence test function that
   would be passed to `EPSSetConvergenceTestFunction()`.

   Calling Sequence:
+  eps    - the linear eigensolver context
.  eigr   - real part of the eigenvalue
.  eigi   - imaginary part of the eigenvalue
.  res    - residual norm associated to the eigenpair
.  errest - [output] computed error estimate
-  ctx    - optional convergence context, as set by `EPSSetConvergenceTestFunction()`

   Level: advanced

.seealso: [](ch:eps), `EPSSetConvergenceTest()`, `EPSSetConvergenceTestFunction()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode EPSConvergenceTestFn(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx);

SLEPC_EXTERN PetscErrorCode EPSSetConvergenceTest(EPS,EPSConv);
SLEPC_EXTERN PetscErrorCode EPSGetConvergenceTest(EPS,EPSConv*);
SLEPC_EXTERN EPSConvergenceTestFn EPSConvergedAbsolute;
SLEPC_EXTERN EPSConvergenceTestFn EPSConvergedRelative;
SLEPC_EXTERN EPSConvergenceTestFn EPSConvergedNorm;
SLEPC_EXTERN PetscErrorCode EPSSetConvergenceTestFunction(EPS,EPSConvergenceTestFn*,void*,PetscCtxDestroyFn*);

/*S
   EPSStoppingTestFn - A prototype of an `EPS` stopping test function that would
   be passed to `EPSSetStoppingTestFunction()`.

   Calling Sequence:
+  eps    - the linear eigensolver context
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged eigenpairs
.  nev    - number of requested eigenpairs
.  reason - [output] result of the stopping test
-  ctx    - optional stopping context, as set by `EPSSetStoppingTestFunction()`

   Note:
   A positive value of `reason` indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, `reason` must be set to `EPS_CONVERGED_ITERATING`
   (zero).

   Level: advanced

.seealso: [](ch:eps), `EPSSetStoppingTest()`, `EPSSetStoppingTestFunction()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode EPSStoppingTestFn(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx);

SLEPC_EXTERN PetscErrorCode EPSSetStoppingTest(EPS,EPSStop);
SLEPC_EXTERN PetscErrorCode EPSGetStoppingTest(EPS,EPSStop*);
SLEPC_EXTERN EPSStoppingTestFn EPSStoppingBasic;
SLEPC_EXTERN EPSStoppingTestFn EPSStoppingThreshold;
SLEPC_EXTERN PetscErrorCode EPSSetStoppingTestFunction(EPS,EPSStoppingTestFn*,void*,PetscCtxDestroyFn*);

SLEPC_EXTERN PetscErrorCode EPSSetEigenvalueComparison(EPS,SlepcEigenvalueComparisonFn*,void*);
SLEPC_EXTERN PetscErrorCode EPSSetArbitrarySelection(EPS,SlepcArbitrarySelectionFn*,void*);

/* --------- options specific to particular eigensolvers -------- */

/*E
   EPSPowerShiftType - The type of shift used in the Power iteration solver.

   Values:
+  `EPS_POWER_SHIFT_CONSTANT`  - constant shift
.  `EPS_POWER_SHIFT_RAYLEIGH`  - variable shift using Rayleigh quotient
-  `EPS_POWER_SHIFT_WILKINSON` - variable shift using Wilkinson's approach

   Note:
   Details of the three variants can be found in {cite:p}`Her05`.

   Level: advanced

.seealso: [](ch:eps), `EPSPowerSetShiftType()`, `EPSPowerGetShiftType()`
E*/
typedef enum { EPS_POWER_SHIFT_CONSTANT,
               EPS_POWER_SHIFT_RAYLEIGH,
               EPS_POWER_SHIFT_WILKINSON } EPSPowerShiftType;
SLEPC_EXTERN const char *EPSPowerShiftTypes[];

/*MC
   EPS_POWER_SHIFT_CONSTANT - The power iteration will use a constant shift.

   Note:
   Together with `STSINVERT`, the `EPSPOWER` solver implements the inverse iteration
   method, i.e., it will apply $(A-\sigma I)^{-1}$ at each iteration, by solving
   a linear system. By default, the shift $\sigma$ is constant and given by the
   user with `EPSSetTarget()`.

   Details of the three variants can be found in {cite:p}`Her05`.

   Level: advanced

.seealso: [](ch:eps), `EPSPowerShiftType`, `EPSPowerSetShiftType()`, `STSetShift()`, `EPSSetTarget()`, `EPS_POWER_SHIFT_RAYLEIGH`, `EPS_POWER_SHIFT_WILKINSON`
M*/

/*MC
   EPS_POWER_SHIFT_RAYLEIGH - The power iteration will use a variable shift
   computed with the Rayleigh quotient.

   Notes:
   Together with `STSINVERT`, the `EPSPOWER` solver implements the inverse iteration
   method, i.e., it will apply $(A-\sigma I)^{-1}$ at each iteration, by solving
   a linear system. With this strategy, the value of the shift will be updated at
   each iteration as $\sigma=\frac{x^*Ax}{x^*x}$, where $x$ is the current eigenvector
   approximation. The resulting iteration is the RQI method.

   Updating the shift may involve a high computational cost if the linear solve
   is done via a factorization.

   Details of the three variants can be found in {cite:p}`Her05`.

   Level: advanced

.seealso: [](ch:eps), `EPSPowerShiftType`, `EPSPowerSetShiftType()`, `STSetShift()`, `EPSSetTarget()`, `EPS_POWER_SHIFT_CONSTANT`, `EPS_POWER_SHIFT_WILKINSON`
M*/

/*MC
   EPS_POWER_SHIFT_WILKINSON - The power iteration will use a variable shift
   computed with Wilkinson's approach.

   Note:
   Together with `STSINVERT`, the `EPSPOWER` solver implements the inverse iteration
   method, i.e., it will apply $(A-\sigma I)^{-1}$ at each iteration, by solving
   a linear system. With this strategy, the value of the shift will be updated at
   each iteration as proposed by Wilkinson, see {cite:p}`Par80{8.10}`.

   Updating the shift may involve a high computational cost if the linear solve
   is done via a factorization.

   Details of the three variants can be found in {cite:p}`Her05`.

   Level: advanced

.seealso: [](ch:eps), `EPSPowerShiftType`, `EPSPowerSetShiftType()`, `STSetShift()`, `EPSSetTarget()`, `EPS_POWER_SHIFT_CONSTANT`, `EPS_POWER_SHIFT_RAYLEIGH`
M*/

SLEPC_EXTERN PetscErrorCode EPSPowerSetShiftType(EPS,EPSPowerShiftType);
SLEPC_EXTERN PetscErrorCode EPSPowerGetShiftType(EPS,EPSPowerShiftType*);
SLEPC_EXTERN PetscErrorCode EPSPowerSetNonlinear(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSPowerGetNonlinear(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSPowerSetUpdate(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSPowerGetUpdate(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSPowerSetSignNormalization(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSPowerGetSignNormalization(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSPowerSetSNES(EPS,SNES);
SLEPC_EXTERN PetscErrorCode EPSPowerGetSNES(EPS,SNES*);

SLEPC_EXTERN PetscErrorCode EPSArnoldiSetDelayed(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSArnoldiGetDelayed(EPS,PetscBool*);

/*E
   EPSKrylovSchurBSEType - The method to be used in the Krylov-Schur solver
   for the case of BSE structured eigenproblems.

   Values:
+  `EPS_KRYLOVSCHUR_BSE_SHAO`         - a Lanczos method proposed by Shao and coauthors
.  `EPS_KRYLOVSCHUR_BSE_GRUNING`      - a Lanczos method proposed by Gruning and coauthors
-  `EPS_KRYLOVSCHUR_BSE_PROJECTEDBSE` - a Lanczos method resulting is a projected problem with BSE structure

   Note:
   All variants are implemented in combination with a thick restart, see
   the details in {cite:p}`Alv25`.

   Level: advanced

.seealso: [](ch:eps), `EPSKrylovSchurSetBSEType()`, `EPSKrylovSchurGetBSEType()`
E*/
typedef enum { EPS_KRYLOVSCHUR_BSE_SHAO,
               EPS_KRYLOVSCHUR_BSE_GRUNING,
               EPS_KRYLOVSCHUR_BSE_PROJECTEDBSE } EPSKrylovSchurBSEType;
SLEPC_EXTERN const char *EPSKrylovSchurBSETypes[];

SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetBSEType(EPS,EPSKrylovSchurBSEType);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetBSEType(EPS,EPSKrylovSchurBSEType*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetRestart(EPS,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetRestart(EPS,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetLocking(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetLocking(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetPartitions(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetPartitions(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetDetectZeros(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetDetectZeros(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetDimensions(EPS,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetDimensions(EPS,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetSubintervals(EPS,PetscReal[]);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetSubintervals(EPS,PetscReal*[]);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetInertias(EPS,PetscInt*,PetscReal*[],PetscInt*[]);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetSubcommInfo(EPS,PetscInt*,PetscInt*,Vec*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetSubcommPairs(EPS,PetscInt,PetscScalar*,Vec);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetSubcommMats(EPS,Mat*,Mat*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurUpdateSubcommMats(EPS,PetscScalar,PetscScalar,Mat,PetscScalar,PetscScalar, Mat,MatStructure,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetKSP(EPS,KSP*);

/*E
   EPSLanczosReorthogType - The type of reorthogonalization used in `EPSLANCZOS`.

   Values:
+  `EPS_LANCZOS_REORTHOG_LOCAL`     - local orthogonalization, only involves previous two vectors
.  `EPS_LANCZOS_REORTHOG_FULL`      - full orthogonalization against all previous vectors
.  `EPS_LANCZOS_REORTHOG_SELECTIVE` - selective reorthogonalization against nearly converged Ritz vectors
.  `EPS_LANCZOS_REORTHOG_PERIODIC`  - periodic reorthogonalization against all previous vectors
.  `EPS_LANCZOS_REORTHOG_PARTIAL`   - partial reorthogonalization against as subset of previous vectors
-  `EPS_LANCZOS_REORTHOG_DELAYED`   - full orthogonalization with delayed reorthogonalization

   Note:
   Details of the different reorthogonalization strategies can be found in
   {cite:p}`Her06`.

   Level: advanced

.seealso: [](ch:eps), `EPSLanczosSetReorthog()`, `EPSLanczosGetReorthog()`
E*/
typedef enum { EPS_LANCZOS_REORTHOG_LOCAL,
               EPS_LANCZOS_REORTHOG_FULL,
               EPS_LANCZOS_REORTHOG_SELECTIVE,
               EPS_LANCZOS_REORTHOG_PERIODIC,
               EPS_LANCZOS_REORTHOG_PARTIAL,
               EPS_LANCZOS_REORTHOG_DELAYED } EPSLanczosReorthogType;
SLEPC_EXTERN const char *EPSLanczosReorthogTypes[];

SLEPC_EXTERN PetscErrorCode EPSLanczosSetReorthog(EPS,EPSLanczosReorthogType);
SLEPC_EXTERN PetscErrorCode EPSLanczosGetReorthog(EPS,EPSLanczosReorthogType*);

/*E
   EPSPRIMMEMethod - The method selected in the PRIMME library.

   Note:
   See the documentation of PRIMME {cite:p}`Sta10` for a description of the methods.

   Level: advanced

.seealso: [](ch:eps), `EPSPRIMMESetMethod()`, `EPSPRIMMEGetMethod()`
E*/
typedef enum { EPS_PRIMME_DYNAMIC             = 1,
               EPS_PRIMME_DEFAULT_MIN_TIME    = 2,
               EPS_PRIMME_DEFAULT_MIN_MATVECS = 3,
               EPS_PRIMME_ARNOLDI             = 4,
               EPS_PRIMME_GD                  = 5,
               EPS_PRIMME_GD_PLUSK            = 6,
               EPS_PRIMME_GD_OLSEN_PLUSK      = 7,
               EPS_PRIMME_JD_OLSEN_PLUSK      = 8,
               EPS_PRIMME_RQI                 = 9,
               EPS_PRIMME_JDQR                = 10,
               EPS_PRIMME_JDQMR               = 11,
               EPS_PRIMME_JDQMR_ETOL          = 12,
               EPS_PRIMME_SUBSPACE_ITERATION  = 13,
               EPS_PRIMME_LOBPCG_ORTHOBASIS   = 14,
               EPS_PRIMME_LOBPCG_ORTHOBASISW  = 15 } EPSPRIMMEMethod;
SLEPC_EXTERN const char *EPSPRIMMEMethods[];

SLEPC_EXTERN PetscErrorCode EPSPRIMMESetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSPRIMMEGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSPRIMMESetMethod(EPS,EPSPRIMMEMethod);
SLEPC_EXTERN PetscErrorCode EPSPRIMMEGetMethod(EPS,EPSPRIMMEMethod*);

SLEPC_EXTERN PetscErrorCode EPSGDSetKrylovStart(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGDGetKrylovStart(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSGDSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGDGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSGDSetRestart(EPS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGDGetRestart(EPS,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSGDSetInitialSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGDGetInitialSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSGDSetBOrth(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGDGetBOrth(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSGDSetDoubleExpansion(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGDGetDoubleExpansion(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSJDSetKrylovStart(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSJDGetKrylovStart(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSJDSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSJDGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSJDSetRestart(EPS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSJDGetRestart(EPS,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSJDSetInitialSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSJDGetInitialSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSJDSetFix(EPS,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSJDGetFix(EPS,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSJDSetConstCorrectionTol(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSJDGetConstCorrectionTol(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSJDSetBOrth(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSJDGetBOrth(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSRQCGSetReset(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSRQCGGetReset(EPS,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSLOBPCGSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGSetRestart(EPS,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGGetRestart(EPS,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGSetLocking(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGGetLocking(EPS,PetscBool*);

/*E
   EPSCISSQuadRule - The quadrature rule used in the `EPSCISS` solver.

   Values:
+  `EPS_CISS_QUADRULE_TRAPEZOIDAL` - trapezoidal rule
-  `EPS_CISS_QUADRULE_CHEBYSHEV`   - Gauss quadrature on Chebyshev points

   Note:
   For a detailed description see {cite:p}`Mae16`.

   Level: advanced

.seealso: [](ch:eps), `EPSCISSSetQuadRule()`, `EPSCISSGetQuadRule()`
E*/
typedef enum { EPS_CISS_QUADRULE_TRAPEZOIDAL = 1,
               EPS_CISS_QUADRULE_CHEBYSHEV   = 2 } EPSCISSQuadRule;
SLEPC_EXTERN const char *EPSCISSQuadRules[];

/*E
   EPSCISSExtraction - The extraction technique used in the `EPSCISS` solver.

   Values:
+  `EPS_CISS_EXTRACTION_RITZ`   - Ritz approximations from Rayleigh-Ritz projection
-  `EPS_CISS_EXTRACTION_HANKEL` - use Hankel pencil as in the original Sakurai-Sugiura method

   Note:
   For a detailed description see {cite:p}`Mae16`.

   Level: advanced

.seealso: [](ch:eps), `EPSCISSSetExtraction()`, `EPSCISSGetExtraction()`
E*/
typedef enum { EPS_CISS_EXTRACTION_RITZ,
               EPS_CISS_EXTRACTION_HANKEL } EPSCISSExtraction;
SLEPC_EXTERN const char *EPSCISSExtractions[];

SLEPC_EXTERN PetscErrorCode EPSCISSSetExtraction(EPS,EPSCISSExtraction);
SLEPC_EXTERN PetscErrorCode EPSCISSGetExtraction(EPS,EPSCISSExtraction*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetQuadRule(EPS,EPSCISSQuadRule);
SLEPC_EXTERN PetscErrorCode EPSCISSGetQuadRule(EPS,EPSCISSQuadRule*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetSizes(EPS,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSCISSGetSizes(EPS,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetThreshold(EPS,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSCISSGetThreshold(EPS,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetRefinement(EPS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSCISSGetRefinement(EPS,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetUseST(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSCISSGetUseST(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSCISSGetKSPs(EPS,PetscInt*,KSP*[]);

SLEPC_EXTERN PetscErrorCode EPSLyapIISetLME(EPS,LME);
SLEPC_EXTERN PetscErrorCode EPSLyapIIGetLME(EPS,LME*);
SLEPC_EXTERN PetscErrorCode EPSLyapIISetRanks(EPS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSLyapIIGetRanks(EPS,PetscInt*,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSBLOPEXSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSBLOPEXGetBlockSize(EPS,PetscInt*);

/*E
   EPSEVSLDOSMethod - The method to approximate the density of states (DOS)
   in the `EPSEVSL` solver.

   Values:
+  `EPS_EVSL_DOS_KPM`     - Kernel Polynomial Method
-  `EPS_EVSL_DOS_LANCZOS` - Lanczos method

   Note:
   See the documentation of EVSL {cite:p}`Li19` for an explanation.

   Level: advanced

.seealso: [](ch:eps), `EPSEVSLSetDOSParameters()`, `EPSEVSLGetDOSParameters()`
E*/
typedef enum { EPS_EVSL_DOS_KPM,
               EPS_EVSL_DOS_LANCZOS } EPSEVSLDOSMethod;
SLEPC_EXTERN const char *EPSEVSLDOSMethods[];

/*E
   EPSEVSLDamping - The damping type used in the `EPSEVSL` solver.

   Values:
+  `EPS_EVSL_DAMPING_NONE`    - no damping
.  `EPS_EVSL_DAMPING_JACKSON` - Jackson damping
-  `EPS_EVSL_DAMPING_SIGMA`   - Lanczos damping

   Note:
   See the documentation of EVSL {cite:p}`Li19` for an explanation.

   Level: advanced

.seealso: [](ch:eps), `EPSEVSLSetDOSParameters()`, `EPSEVSLGetDOSParameters()`
E*/
typedef enum { EPS_EVSL_DAMPING_NONE,
               EPS_EVSL_DAMPING_JACKSON,
               EPS_EVSL_DAMPING_SIGMA } EPSEVSLDamping;
SLEPC_EXTERN const char *EPSEVSLDampings[];

SLEPC_EXTERN PetscErrorCode EPSEVSLSetRange(EPS,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSEVSLGetRange(EPS,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSEVSLSetSlices(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSEVSLGetSlices(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSEVSLSetDOSParameters(EPS,EPSEVSLDOSMethod,PetscInt,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSEVSLGetDOSParameters(EPS,EPSEVSLDOSMethod*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSEVSLSetPolParameters(EPS,PetscInt,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSEVSLGetPolParameters(EPS,PetscInt*,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSEVSLSetDamping(EPS,EPSEVSLDamping);
SLEPC_EXTERN PetscErrorCode EPSEVSLGetDamping(EPS,EPSEVSLDamping*);

SLEPC_EXTERN PetscErrorCode EPSFEASTSetNumPoints(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSFEASTGetNumPoints(EPS,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSCHASESetDegree(EPS,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSCHASEGetDegree(EPS,PetscInt*,PetscBool*);
