/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for SLEPc's singular value solvers
*/

#pragma once

#include <slepceps.h>
#include <slepcbv.h>
#include <slepcds.h>

/* SUBMANSEC = SVD */

SLEPC_EXTERN PetscErrorCode SVDInitializePackage(void);
SLEPC_EXTERN PetscErrorCode SVDFinalizePackage(void);

/*S
   SVD - SLEPc object that manages all the singular value problem solvers.

   Level: beginner

.seealso: [](ch:svd), `SVDCreate()`
S*/
typedef struct _p_SVD* SVD;

/*J
   SVDType - String with the name of a singular value solver.

   Level: beginner

.seealso: [](ch:svd), `SVDSetType()`, `SVD`
J*/
typedef const char *SVDType;
#define SVDCROSS       "cross"
#define SVDCYCLIC      "cyclic"
#define SVDLANCZOS     "lanczos"
#define SVDTRLANCZOS   "trlanczos"
#define SVDRANDOMIZED  "randomized"
#define SVDLAPACK      "lapack"
#define SVDSCALAPACK   "scalapack"
#define SVDKSVD        "ksvd"
#define SVDELEMENTAL   "elemental"
#define SVDPRIMME      "primme"

/* Logging support */
SLEPC_EXTERN PetscClassId SVD_CLASSID;

/*E
   SVDProblemType - Determines the type of the singular value problem.

   Values:
+  `SVD_STANDARD`    - standard SVD (SVD)
.  `SVD_GENERALIZED` - generalized SVD (GSVD)
-  `SVD_HYPERBOLIC`  - hyperbolic SVD (HSVD)

   Level: beginner

.seealso: [](ch:svd), `SVDSetProblemType()`, `SVDGetProblemType()`
E*/
typedef enum { SVD_STANDARD    = 1,
               SVD_GENERALIZED = 2,
               SVD_HYPERBOLIC  = 3
             } SVDProblemType;

/*MC
   SVD_STANDARD - A (standard) singular value problem (SVD).

   Note:
   The problem is formulated as $A=U\Sigma V^*$, where $A$ is a possibly
   non-square matrix.

   Level: beginner

.seealso: [](ch:svd), `SVDProblemType`, `SVDSetProblemType()`, `SVD_GENERALIZED`, `SVD_HYPERBOLIC`
M*/

/*MC
   SVD_GENERALIZED - A generalized singular value problem (GSVD).

   Note:
   The problem is formulated as $U^*AX=C$, $V^*BX=S$, where $A$ and $B$
   have the same number of columns.

   Level: beginner

.seealso: [](ch:svd), `SVDProblemType`, `SVDSetProblemType()`, `SVD_STANDARD`, `SVD_HYPERBOLIC`
M*/

/*MC
   SVD_HYPERBOLIC - A hyperbolic singular value problem (HSVD).

   Note:
   The problem is formulated as $A=U\Sigma V^*$, with $U^*\Omega U=\tilde\Omega$,
   where $A$ is a possibly non-square matrix, and $\Omega$, $\tilde\Omega$
   are signature matrices.

   Level: beginner

.seealso: [](ch:svd), `SVDProblemType`, `SVDSetProblemType()`, `SVD_STANDARD`, `SVD_GENERALIZED`
M*/

/*E
   SVDWhich - Determines whether largest or smallest singular values
   are to be computed.

   Values:
+  `SVD_LARGEST`  - largest singular values
-  `SVD_SMALLEST` - smallest singular values

   Level: intermediate

.seealso: [](ch:svd), `SVDSetWhichSingularTriplets()`, `SVDGetWhichSingularTriplets()`
E*/
typedef enum { SVD_LARGEST,
               SVD_SMALLEST } SVDWhich;

/*MC
   SVD_LARGEST - The solver is configured to compute largest singular values.

   Note:
   This is the default.

   Level: intermediate

.seealso: [](ch:svd), `SVDWhich`, `SVDSetWhichSingularTriplets()`, `SVD_SMALLEST`
M*/

/*MC
   SVD_SMALLEST - The solver is configured to compute smallest singular values.

   Note:
   Computing small singular values is generally more difficult than computing
   largest ones, because in many cases these values are very small and
   tightly clustered together. In the case of rank-deficient matrices, smallest
   singular values are zero, and this may pose difficulties to the solvers.

   Level: intermediate

.seealso: [](ch:svd), `SVDWhich`, `SVDSetWhichSingularTriplets()`, `SVD_LARGEST`
M*/

/*E
   SVDErrorType - The error type used to assess accuracy of computed solutions.

   Values:
+  `SVD_ERROR_ABSOLUTE` - compute error bound as $\|r\|$
.  `SVD_ERROR_RELATIVE` - compute error bound as $\|r\|/\sigma$
-  `SVD_ERROR_NORM`     - compute error bound as $\|r\|/\max\{\|A\|,\|B\|\}$

   Note:
   The residual norm $\|r\|$ is actually computed from two parts, such as
   $\sqrt{\eta_1^2+\eta_2^2}$ with $\eta_1 = \|Av-\sigma u\|_2$ and
   $\eta_2 = \|A^*u-\sigma v\|_2$, see more details at `SVDComputeError()`.
   There is also a normalization factor related to the norm of the vectors,
   which also varies with the problem type.

   Level: intermediate

.seealso: [](ch:svd), `SVDComputeError()`, `SVDProblemType`
E*/
typedef enum { SVD_ERROR_ABSOLUTE,
               SVD_ERROR_RELATIVE,
               SVD_ERROR_NORM } SVDErrorType;
SLEPC_EXTERN const char *SVDErrorTypes[];

/*E
   SVDConv - The convergence criterion to be used by the solver.

   Values:
+  `SVD_CONV_ABS`   - absolute convergence criterion, $\|r\|$
.  `SVD_CONV_REL`   - convergence criterion relative to singular value, $\|r\|/\sigma$
.  `SVD_CONV_NORM`  - convergence criterion relative to matrix norms, $\|r\|/\max\{\|A\|,\|B\|\}$
.  `SVD_CONV_MAXIT` - no convergence until maximum number of iterations has been reached
-  `SVD_CONV_USER`  - convergence dictated by user-provided function

   Note:
   The `SVD_CONV_MAXIT` convergence criterion is used only in `SVDRANDOMIZED`.

   Level: intermediate

.seealso: [](ch:svd), `SVDSetConvergenceTest()`, `SVDSetConvergenceTestFunction()`, `SVDSetTolerances()`
E*/
typedef enum { SVD_CONV_ABS,
               SVD_CONV_REL,
               SVD_CONV_NORM,
               SVD_CONV_MAXIT,
               SVD_CONV_USER } SVDConv;

/*E
   SVDStop - The stopping test to decide the termination of the outer loop
   of the singular value solver.

   Values:
+  `SVD_STOP_BASIC`     - default stopping test
.  `SVD_STOP_USER`      - user-provided stopping test
-  `SVD_STOP_THRESHOLD` - threshold stopping test

   Level: advanced

.seealso: [](ch:svd), `SVDSetStoppingTest()`, `SVDSetStoppingTestFunction()`
E*/
typedef enum { SVD_STOP_BASIC,
               SVD_STOP_USER,
               SVD_STOP_THRESHOLD } SVDStop;

/*MC
   SVD_STOP_BASIC - The default stopping test.

   Note:
   By default, the termination of the outer loop is decided by calling
   `SVDStoppingBasic()`, which will stop if all requested singular values are converged,
   or if the maximum number of iterations has been reached.

   Level: advanced

.seealso: [](ch:svd), `SVDStop`, `SVDSetStoppingTest()`, `SVDStoppingBasic()`
M*/

/*MC
   SVD_STOP_USER - The user-provided stopping test.

   Note:
   Customized stopping test using the user-provided function given with
   `SVDSetStoppingTestFunction()`.

   Level: advanced

.seealso: [](ch:svd), `SVDStop`, `SVDSetStoppingTest()`, `SVDSetStoppingTestFunction()`
M*/

/*MC
   SVD_STOP_THRESHOLD - The threshold stopping test.

   Note:
   When a threshold has been provided with `SVDSetThreshold()`, the termination
   of the outer loop is decided by calling `SVDStoppingThreshold()`, which will
   stop when one of the computed singular values is not above/below the threshold.
   If a number of wanted singular values has been specified via `SVDSetDimensions()`
   then it is also taken into account, and the solver will stop when one of the
   two conditions (threshold or number of converged values) is met.

   Level: advanced

.seealso: [](ch:svd), `SVDStop`, `SVDSetStoppingTest()`, `SVDStoppingThreshold()`, `SVDSetThreshold()`, `SVDSetDimensions()`
M*/

/*E
   SVDConvergedReason - Reason a singular value solver was determined to have
   converged or diverged.

   Values:
+  `SVD_CONVERGED_TOL`          - converged up to tolerance
.  `SVD_CONVERGED_USER`         - converged due to a user-defined condition
.  `SVD_CONVERGED_MAXIT`        - reached maximum number of iterations with `SVD_CONV_MAXIT` criterion
.  `SVD_DIVERGED_ITS`           - exceeded the maximum number of allowed iterations
.  `SVD_DIVERGED_BREAKDOWN`     - generic breakdown in method
.  `SVD_DIVERGED_SYMMETRY_LOST` - underlying indefinite eigensolver was not able to keep symmetry
-  `SVD_CONVERGED_ITERATING`    - the solver is still running

   Level: intermediate

.seealso: [](ch:svd), `SVDSolve()`, `SVDGetConvergedReason()`, `SVDSetTolerances()`
E*/
typedef enum {/* converged */
              SVD_CONVERGED_TOL                =  1,
              SVD_CONVERGED_USER               =  2,
              SVD_CONVERGED_MAXIT              =  3,
              /* diverged */
              SVD_DIVERGED_ITS                 = -1,
              SVD_DIVERGED_BREAKDOWN           = -2,
              SVD_DIVERGED_SYMMETRY_LOST       = -3,
              SVD_CONVERGED_ITERATING          =  0 } SVDConvergedReason;
SLEPC_EXTERN const char *const*SVDConvergedReasons;

/*MC
   SVD_CONVERGED_TOL - The computed error estimates, based on residual norms,
   for all requested singular values are below the tolerance.

   Level: intermediate

.seealso: [](ch:svd), `SVDSolve()`, `SVDGetConvergedReason()`, `SVDConvergedReason`
M*/

/*MC
   SVD_CONVERGED_USER - The solver was declared converged due to a user-defined condition.

   Note:
   This happens only when a user-defined stopping test has been set with
   `SVDSetStoppingTestFunction()`.

   Level: intermediate

.seealso: [](ch:svd), `SVDSolve()`, `SVDGetConvergedReason()`, `SVDConvergedReason`, `SVDSetStoppingTestFunction()`
M*/

/*MC
   SVD_CONVERGED_MAXIT - The solver has reached the maximum number of iterations
   with the `SVD_CONV_MAXIT` criterion.

   Note:
   This is considered a successful exit, because the user wanted to do a fixed
   number of iterations. But be aware that the computed solution may be inaccurate,
   in particular, individual singular vectors will not have good residual. This
   is available in `SVDRANDOMIZED` only.

   Level: intermediate

.seealso: [](ch:svd), `SVDSolve()`, `SVDGetConvergedReason()`, `SVDConvergedReason`, `SVD_CONV_MAXIT`, `SVDRANDOMIZED`
M*/

/*MC
   SVD_DIVERGED_ITS - Exceeded the maximum number of allowed iterations
   before the convergence criterion was satisfied.

   Level: intermediate

.seealso: [](ch:svd), `SVDSolve()`, `SVDGetConvergedReason()`, `SVDConvergedReason`
M*/

/*MC
   SVD_DIVERGED_BREAKDOWN - A breakdown in the solver was detected so the
   method could not continue.

   Level: intermediate

.seealso: [](ch:svd), `SVDSolve()`, `SVDGetConvergedReason()`, `SVDConvergedReason`
M*/

/*MC
   SVD_DIVERGED_SYMMETRY_LOST - The selected solver uses a pseudo-Lanczos recurrence,
   which is numerically unstable, and a symmetry test revealed that instability
   had appeared so the solver could not continue.

   Level: intermediate

.seealso: [](ch:svd), `SVDSolve()`, `SVDGetConvergedReason()`, `SVDConvergedReason`
M*/

/*MC
   SVD_CONVERGED_ITERATING - This value is returned if `SVDGetConvergedReason()` is called
   while `SVDSolve()` is still running.

   Level: intermediate

.seealso: [](ch:svd), `SVDSolve()`, `SVDGetConvergedReason()`, `SVDConvergedReason`
M*/

/*S
   SVDStoppingCtx - Data structure (C struct) to hold additional information to
   be used in some stopping test functions.

   Level: advanced

.seealso: [](ch:svd), `SVDSetStoppingTestFunction()`
S*/
struct _n_SVDStoppingCtx {
  PetscReal firstsv;    /* the value of the first converged singular value */
  PetscReal lastsv;     /* the value of the last converged singular value */
  PetscReal thres;      /* threshold set with SVDSetThreshold() */
  PetscBool threlative; /* threshold is relative */
  SVDWhich  which;      /* which singular values are being computed */
};
typedef struct _n_SVDStoppingCtx* SVDStoppingCtx;

SLEPC_EXTERN PetscErrorCode SVDCreate(MPI_Comm,SVD*);
SLEPC_EXTERN PetscErrorCode SVDSetBV(SVD,BV,BV);
SLEPC_EXTERN PetscErrorCode SVDGetBV(SVD,BV*,BV*);
SLEPC_EXTERN PetscErrorCode SVDSetDS(SVD,DS);
SLEPC_EXTERN PetscErrorCode SVDGetDS(SVD,DS*);
SLEPC_EXTERN PetscErrorCode SVDSetType(SVD,SVDType);
SLEPC_EXTERN PetscErrorCode SVDGetType(SVD,SVDType*);
SLEPC_EXTERN PetscErrorCode SVDSetProblemType(SVD,SVDProblemType);
SLEPC_EXTERN PetscErrorCode SVDGetProblemType(SVD,SVDProblemType*);
SLEPC_EXTERN PetscErrorCode SVDIsGeneralized(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDIsHyperbolic(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDSetOperators(SVD,Mat,Mat);
PETSC_DEPRECATED_FUNCTION(3, 15, 0, "SVDSetOperators()", ) static inline PetscErrorCode SVDSetOperator(SVD svd,Mat A) {return SVDSetOperators(svd,A,PETSC_NULLPTR);}
SLEPC_EXTERN PetscErrorCode SVDGetOperators(SVD,Mat*,Mat*);
PETSC_DEPRECATED_FUNCTION(3, 15, 0, "SVDGetOperators()", ) static inline PetscErrorCode SVDGetOperator(SVD svd,Mat *A) {return SVDGetOperators(svd,A,PETSC_NULLPTR);}
SLEPC_EXTERN PetscErrorCode SVDSetSignature(SVD,Vec);
SLEPC_EXTERN PetscErrorCode SVDGetSignature(SVD,Vec);
SLEPC_EXTERN PetscErrorCode SVDSetInitialSpaces(SVD,PetscInt,Vec[],PetscInt,Vec[]);
PETSC_DEPRECATED_FUNCTION(3, 1, 0, "SVDSetInitialSpaces()", ) static inline PetscErrorCode SVDSetInitialSpace(SVD svd,PetscInt nr,Vec *isr) {return SVDSetInitialSpaces(svd,nr,isr,0,PETSC_NULLPTR);}
PETSC_DEPRECATED_FUNCTION(3, 1, 0, "SVDSetInitialSpaces()", ) static inline PetscErrorCode SVDSetInitialSpaceLeft(SVD svd,PetscInt nl,Vec *isl) {return SVDSetInitialSpaces(svd,0,PETSC_NULLPTR,nl,isl);}
SLEPC_EXTERN PetscErrorCode SVDSetImplicitTranspose(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDGetImplicitTranspose(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDSetDimensions(SVD,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDGetDimensions(SVD,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDSetTolerances(SVD,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDGetTolerances(SVD,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDSetWhichSingularTriplets(SVD,SVDWhich);
SLEPC_EXTERN PetscErrorCode SVDGetWhichSingularTriplets(SVD,SVDWhich*);
SLEPC_EXTERN PetscErrorCode SVDSetThreshold(SVD,PetscReal,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDGetThreshold(SVD,PetscReal*,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDSetFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDSetOptionsPrefix(SVD,const char[]);
SLEPC_EXTERN PetscErrorCode SVDAppendOptionsPrefix(SVD,const char[]);
SLEPC_EXTERN PetscErrorCode SVDGetOptionsPrefix(SVD,const char*[]);
SLEPC_EXTERN PetscErrorCode SVDSetDSType(SVD);
SLEPC_EXTERN PetscErrorCode SVDSetUp(SVD);
SLEPC_EXTERN PetscErrorCode SVDSolve(SVD);
SLEPC_EXTERN PetscErrorCode SVDGetIterationNumber(SVD,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDGetConvergedReason(SVD,SVDConvergedReason*);
SLEPC_EXTERN PetscErrorCode SVDGetConverged(SVD,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDGetSingularTriplet(SVD,PetscInt,PetscReal*,Vec,Vec);
SLEPC_EXTERN PetscErrorCode SVDComputeError(SVD,PetscInt,SVDErrorType,PetscReal*);
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "SVDComputeError()", ) static inline PetscErrorCode SVDComputeRelativeError(SVD svd,PetscInt i,PetscReal *r) {return SVDComputeError(svd,i,SVD_ERROR_RELATIVE,r);}
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "SVDComputeError() with SVD_ERROR_ABSOLUTE", ) static inline PetscErrorCode SVDComputeResidualNorms(SVD svd,PetscInt i,PetscReal *r1,PETSC_UNUSED PetscReal *r2) {return SVDComputeError(svd,i,SVD_ERROR_ABSOLUTE,r1);}
SLEPC_EXTERN PetscErrorCode SVDView(SVD,PetscViewer);
SLEPC_EXTERN PetscErrorCode SVDViewFromOptions(SVD,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode SVDErrorView(SVD,SVDErrorType,PetscViewer);
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "SVDErrorView()", ) static inline PetscErrorCode SVDPrintSolution(SVD svd,PetscViewer v) {return SVDErrorView(svd,SVD_ERROR_RELATIVE,v);}
SLEPC_EXTERN PetscErrorCode SVDErrorViewFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDConvergedReasonView(SVD,PetscViewer);
SLEPC_EXTERN PetscErrorCode SVDConvergedReasonViewFromOptions(SVD);
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "SVDConvergedReasonView()", ) static inline PetscErrorCode SVDReasonView(SVD svd,PetscViewer v) {return SVDConvergedReasonView(svd,v);}
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "SVDConvergedReasonViewFromOptions()", ) static inline PetscErrorCode SVDReasonViewFromOptions(SVD svd) {return SVDConvergedReasonViewFromOptions(svd);}
SLEPC_EXTERN PetscErrorCode SVDValuesView(SVD,PetscViewer);
SLEPC_EXTERN PetscErrorCode SVDValuesViewFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDVectorsView(SVD,PetscViewer);
SLEPC_EXTERN PetscErrorCode SVDVectorsViewFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDDestroy(SVD*);
SLEPC_EXTERN PetscErrorCode SVDReset(SVD);
SLEPC_EXTERN PetscErrorCode SVDSetWorkVecs(SVD,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDSetTrackAll(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDGetTrackAll(SVD,PetscBool*);

/*S
   SVDMonitorFn - A function prototype for functions provided to `SVDMonitorSet()`.

   Calling Sequence:
+  svd    - the singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets
.  sigma  - singular values
.  errest - relative error estimates for each singular triplet
.  nest   - number of error estimates
-  ctx    - optional monitoring context, as provided with `SVDMonitorSet()`

   Level: intermediate

.seealso: [](ch:svd), `SVDMonitorSet()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode SVDMonitorFn(SVD svd,PetscInt its,PetscInt nconv,PetscReal sigma[],PetscReal errest[],PetscInt nest,void *ctx);

/*S
   SVDMonitorRegisterFn - A function prototype for functions provided to `SVDMonitorRegister()`.

   Calling Sequence:
+  svd    - the singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets
.  sigma  - singular values
.  errest - relative error estimates for each singular triplet
.  nest   - number of error estimates
-  ctx    - `PetscViewerAndFormat` object

   Level: advanced

   Note:
   This is an `SVDMonitorFn` specialized for a context of `PetscViewerAndFormat`.

.seealso: [](ch:svd), `SVDMonitorSet()`, `SVDMonitorRegister()`, `SVDMonitorFn`, `SVDMonitorRegisterCreateFn`, `SVDMonitorRegisterDestroyFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode SVDMonitorRegisterFn(SVD svd,PetscInt its,PetscInt nconv,PetscReal sigma[],PetscReal errest[],PetscInt nest,PetscViewerAndFormat *ctx);

/*S
   SVDMonitorRegisterCreateFn - A function prototype for functions that do
   the creation when provided to `SVDMonitorRegister()`.

   Calling Sequence:
+  viewer - the viewer to be used with the `SVDMonitorRegisterFn`
.  format - the format of the viewer
.  ctx    - a context for the monitor
-  result - a `PetscViewerAndFormat` object

   Level: advanced

.seealso: [](ch:svd), `SVDMonitorRegisterFn`, `SVDMonitorSet()`, `SVDMonitorRegister()`, `SVDMonitorFn`, `SVDMonitorRegisterDestroyFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode SVDMonitorRegisterCreateFn(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **result);

/*S
   SVDMonitorRegisterDestroyFn - A function prototype for functions that do the after
   use destruction when provided to `SVDMonitorRegister()`.

   Calling Sequence:
.  vf - a `PetscViewerAndFormat` object to be destroyed, including any context

   Level: advanced

.seealso: [](ch:svd), `SVDMonitorRegisterFn`, `SVDMonitorSet()`, `SVDMonitorRegister()`, `SVDMonitorFn`, `SVDMonitorRegisterCreateFn`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode SVDMonitorRegisterDestroyFn(PetscViewerAndFormat **result);

SLEPC_EXTERN PetscErrorCode SVDMonitor(SVD,PetscInt,PetscInt,PetscReal[],PetscReal[],PetscInt);
SLEPC_EXTERN PetscErrorCode SVDMonitorSet(SVD,SVDMonitorFn,void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode SVDMonitorCancel(SVD);
SLEPC_EXTERN PetscErrorCode SVDGetMonitorContext(SVD,void*);

SLEPC_EXTERN PetscErrorCode SVDMonitorSetFromOptions(SVD,const char[],const char[],void*,PetscBool);
SLEPC_EXTERN SVDMonitorRegisterFn        SVDMonitorFirst;
SLEPC_EXTERN SVDMonitorRegisterFn        SVDMonitorFirstDrawLG;
SLEPC_EXTERN SVDMonitorRegisterCreateFn  SVDMonitorFirstDrawLGCreate;
SLEPC_EXTERN SVDMonitorRegisterFn        SVDMonitorAll;
SLEPC_EXTERN SVDMonitorRegisterFn        SVDMonitorAllDrawLG;
SLEPC_EXTERN SVDMonitorRegisterCreateFn  SVDMonitorAllDrawLGCreate;
SLEPC_EXTERN SVDMonitorRegisterFn        SVDMonitorConverged;
SLEPC_EXTERN SVDMonitorRegisterCreateFn  SVDMonitorConvergedCreate;
SLEPC_EXTERN SVDMonitorRegisterFn        SVDMonitorConvergedDrawLG;
SLEPC_EXTERN SVDMonitorRegisterCreateFn  SVDMonitorConvergedDrawLGCreate;
SLEPC_EXTERN SVDMonitorRegisterDestroyFn SVDMonitorConvergedDestroy;
SLEPC_EXTERN SVDMonitorRegisterFn        SVDMonitorConditioning;

SLEPC_EXTERN PetscFunctionList SVDList;
SLEPC_EXTERN PetscFunctionList SVDMonitorList;
SLEPC_EXTERN PetscFunctionList SVDMonitorCreateList;
SLEPC_EXTERN PetscFunctionList SVDMonitorDestroyList;
SLEPC_EXTERN PetscErrorCode SVDRegister(const char[],PetscErrorCode(*)(SVD));
SLEPC_EXTERN PetscErrorCode SVDMonitorRegister(const char[],PetscViewerType,PetscViewerFormat,SVDMonitorRegisterFn*,SVDMonitorRegisterCreateFn*,SVDMonitorRegisterDestroyFn*);

SLEPC_EXTERN PetscErrorCode SVDAllocateSolution(SVD,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDReallocateSolution(SVD,PetscInt);

/*S
   SVDConvergenceTestFn - A prototype of an `SVD` convergence test function that
   would be passed to `SVDSetConvergenceTestFunction()`.

   Calling Sequence:
+  svd    - the singular value solver context
.  sigma  - computed singular value
.  res    - residual norm associated to the singular triplet
.  errest - [output] computed error estimate
-  ctx    - optional convergence context, as set by `SVDSetConvergenceTestFunction()`

   Level: advanced

.seealso: [](ch:svd), `SVDSetConvergenceTest()`, `SVDSetConvergenceTestFunction()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode SVDConvergenceTestFn(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx);

SLEPC_EXTERN PetscErrorCode SVDSetConvergenceTest(SVD,SVDConv);
SLEPC_EXTERN PetscErrorCode SVDGetConvergenceTest(SVD,SVDConv*);
SLEPC_EXTERN SVDConvergenceTestFn SVDConvergedAbsolute;
SLEPC_EXTERN SVDConvergenceTestFn SVDConvergedRelative;
SLEPC_EXTERN SVDConvergenceTestFn SVDConvergedNorm;
SLEPC_EXTERN SVDConvergenceTestFn SVDConvergedMaxIt;
SLEPC_EXTERN PetscErrorCode SVDSetConvergenceTestFunction(SVD,SVDConvergenceTestFn*,void*,PetscCtxDestroyFn*);

/*S
   SVDStoppingTestFn - A prototype of an `SVD` stopping test function that would
   be passed to `SVDSetStoppingTestFunction()`.

   Calling Sequence:
+  svd    - the singular value solver context
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged singular triplets
.  nsv    - number of requested singular triplets
.  reason - [output] result of the stopping test
-  ctx    - optional stopping context, as set by `SVDSetStoppingTestFunction()`

   Note:
   A positive value of `reason` indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, `reason` must be set to `SVD_CONVERGED_ITERATING`
   (zero).

   Level: advanced

.seealso: [](ch:svd), `SVDSetStoppingTest()`, `SVDSetStoppingTestFunction()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode SVDStoppingTestFn(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nsv,SVDConvergedReason *reason,void *ctx);

SLEPC_EXTERN PetscErrorCode SVDSetStoppingTest(SVD,SVDStop);
SLEPC_EXTERN PetscErrorCode SVDGetStoppingTest(SVD,SVDStop*);
SLEPC_EXTERN SVDStoppingTestFn SVDStoppingBasic;
SLEPC_EXTERN SVDStoppingTestFn SVDStoppingThreshold;
SLEPC_EXTERN PetscErrorCode SVDSetStoppingTestFunction(SVD,SVDStoppingTestFn*,void*,PetscCtxDestroyFn*);

/* --------- options specific to particular solvers -------- */

SLEPC_EXTERN PetscErrorCode SVDCrossSetExplicitMatrix(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDCrossGetExplicitMatrix(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDCrossSetEPS(SVD,EPS);
SLEPC_EXTERN PetscErrorCode SVDCrossGetEPS(SVD,EPS*);

SLEPC_EXTERN PetscErrorCode SVDCyclicSetExplicitMatrix(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDCyclicGetExplicitMatrix(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDCyclicSetEPS(SVD,EPS);
SLEPC_EXTERN PetscErrorCode SVDCyclicGetEPS(SVD,EPS*);

SLEPC_EXTERN PetscErrorCode SVDLanczosSetOneSide(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDLanczosGetOneSide(SVD,PetscBool*);

/*E
   SVDTRLanczosGBidiag - The choice of bidiagonalization for the `SVDTRLANCZOS` GSVD solver.

   Values:
+  `SVD_TRLANCZOS_GBIDIAG_SINGLE` - single bidiagonalization ($Q_A$)
.  `SVD_TRLANCZOS_GBIDIAG_UPPER`  - joint bidiagonalization, both $Q_A$ and $Q_B$ in upper bidiagonal form
-  `SVD_TRLANCZOS_GBIDIAG_LOWER`  - joint bidiagonalization, $Q_A$ lower bidiagonal, $Q_B$ upper bidiagonal

   Note:
   The different variants are described in {cite:p}`Alv24`.

   Level: advanced

.seealso: [](ch:svd), `SVDTRLanczosSetGBidiag()`, `SVDTRLanczosGetGBidiag()`
E*/
typedef enum {
  SVD_TRLANCZOS_GBIDIAG_SINGLE,
  SVD_TRLANCZOS_GBIDIAG_UPPER,
  SVD_TRLANCZOS_GBIDIAG_LOWER
} SVDTRLanczosGBidiag;
SLEPC_EXTERN const char *SVDTRLanczosGBidiags[];

SLEPC_EXTERN PetscErrorCode SVDTRLanczosSetGBidiag(SVD,SVDTRLanczosGBidiag);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosGetGBidiag(SVD,SVDTRLanczosGBidiag*);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosSetOneSide(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosGetOneSide(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosSetKSP(SVD,KSP);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosGetKSP(SVD,KSP*);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosSetRestart(SVD,PetscReal);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosGetRestart(SVD,PetscReal*);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosSetLocking(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosGetLocking(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosSetExplicitMatrix(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosGetExplicitMatrix(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosSetScale(SVD,PetscReal);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosGetScale(SVD,PetscReal*);

/*E
   SVDPRIMMEMethod - The SVD method selected in the PRIMME library.

   Note:
   See the documentation of PRIMME {cite:p}`Sta10` for a description of the methods.

   Level: advanced

.seealso: [](ch:svd), `SVDPRIMMESetMethod()`, `SVDPRIMMEGetMethod()`
E*/
typedef enum { SVD_PRIMME_HYBRID          = 1,
               SVD_PRIMME_NORMALEQUATIONS = 2,
               SVD_PRIMME_AUGMENTED       = 3 } SVDPRIMMEMethod;
SLEPC_EXTERN const char *SVDPRIMMEMethods[];

SLEPC_EXTERN PetscErrorCode SVDPRIMMESetBlockSize(SVD,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDPRIMMEGetBlockSize(SVD,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDPRIMMESetMethod(SVD,SVDPRIMMEMethod);
SLEPC_EXTERN PetscErrorCode SVDPRIMMEGetMethod(SVD,SVDPRIMMEMethod*);

/*E
   SVDKSVDEigenMethod - The method to solve the eigenproblem within the KSVD library.

   Note:
   See the documentation of KSVD {cite:p}`Suk19` for a description of the methods.

   Level: advanced

.seealso: [](ch:svd), `SVDKSVDSetEigenMethod()`, `SVDKSVDGetEigenMethod()`
E*/
typedef enum { SVD_KSVD_EIGEN_MRRR = 1,
               SVD_KSVD_EIGEN_DC   = 2,
               SVD_KSVD_EIGEN_ELPA = 3 } SVDKSVDEigenMethod;
SLEPC_EXTERN const char *SVDKSVDEigenMethods[];

/*E
   SVDKSVDPolarMethod - The method to compute the polar decomposition within the KSVD library.

   Note:
   See the documentation of KSVD {cite:p}`Suk19` for a description of the methods.

   Level: advanced

.seealso: [](ch:svd), `SVDKSVDSetPolarMethod()`, `SVDKSVDGetPolarMethod()`
E*/
typedef enum { SVD_KSVD_POLAR_QDWH   = 1,
               SVD_KSVD_POLAR_ZOLOPD = 2 } SVDKSVDPolarMethod;
SLEPC_EXTERN const char *SVDKSVDPolarMethods[];

SLEPC_EXTERN PetscErrorCode SVDKSVDSetEigenMethod(SVD,SVDKSVDEigenMethod);
SLEPC_EXTERN PetscErrorCode SVDKSVDGetEigenMethod(SVD,SVDKSVDEigenMethod*);
SLEPC_EXTERN PetscErrorCode SVDKSVDSetPolarMethod(SVD,SVDKSVDPolarMethod);
SLEPC_EXTERN PetscErrorCode SVDKSVDGetPolarMethod(SVD,SVDKSVDPolarMethod*);
