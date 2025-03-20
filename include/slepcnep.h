/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for SLEPc's nonlinear eigenvalue solvers
*/

#pragma once

#include <slepceps.h>
#include <slepcpep.h>
#include <slepcfn.h>

/* SUBMANSEC = NEP */

SLEPC_EXTERN PetscErrorCode NEPInitializePackage(void);
SLEPC_EXTERN PetscErrorCode NEPFinalizePackage(void);

/*S
     NEP - Abstract SLEPc object that manages all solvers for
     nonlinear eigenvalue problems.

   Level: beginner

.seealso:  NEPCreate()
S*/
typedef struct _p_NEP* NEP;

/*J
    NEPType - String with the name of a nonlinear eigensolver

   Level: beginner

.seealso: NEPSetType(), NEP
J*/
typedef const char *NEPType;
#define NEPRII       "rii"
#define NEPSLP       "slp"
#define NEPNARNOLDI  "narnoldi"
#define NEPCISS      "ciss"
#define NEPINTERPOL  "interpol"
#define NEPNLEIGS    "nleigs"

/* Logging support */
SLEPC_EXTERN PetscClassId NEP_CLASSID;

/*E
    NEPProblemType - Determines the type of the nonlinear eigenproblem

    Level: intermediate

.seealso: NEPSetProblemType(), NEPGetProblemType()
E*/
typedef enum { NEP_GENERAL  = 1,
               NEP_RATIONAL = 2     /* NEP defined in split form with all f_i rational */
             } NEPProblemType;

/*E
    NEPWhich - Determines which part of the spectrum is requested

    Level: intermediate

.seealso: NEPSetWhichEigenpairs(), NEPGetWhichEigenpairs()
E*/
typedef enum { NEP_LARGEST_MAGNITUDE  = 1,
               NEP_SMALLEST_MAGNITUDE = 2,
               NEP_LARGEST_REAL       = 3,
               NEP_SMALLEST_REAL      = 4,
               NEP_LARGEST_IMAGINARY  = 5,
               NEP_SMALLEST_IMAGINARY = 6,
               NEP_TARGET_MAGNITUDE   = 7,
               NEP_TARGET_REAL        = 8,
               NEP_TARGET_IMAGINARY   = 9,
               NEP_ALL                = 10,
               NEP_WHICH_USER         = 11 } NEPWhich;

/*E
    NEPErrorType - The error type used to assess accuracy of computed solutions

    Level: intermediate

.seealso: NEPComputeError()
E*/
typedef enum { NEP_ERROR_ABSOLUTE,
               NEP_ERROR_RELATIVE,
               NEP_ERROR_BACKWARD } NEPErrorType;
SLEPC_EXTERN const char *NEPErrorTypes[];

/*E
    NEPRefine - The refinement type

    Level: intermediate

.seealso: NEPSetRefine()
E*/
typedef enum { NEP_REFINE_NONE,
               NEP_REFINE_SIMPLE,
               NEP_REFINE_MULTIPLE } NEPRefine;
SLEPC_EXTERN const char *NEPRefineTypes[];

/*E
    NEPRefineScheme - The scheme used for solving linear systems during iterative refinement

    Level: intermediate

.seealso: NEPSetRefine()
E*/
typedef enum { NEP_REFINE_SCHEME_SCHUR    = 1,
               NEP_REFINE_SCHEME_MBE      = 2,
               NEP_REFINE_SCHEME_EXPLICIT = 3 } NEPRefineScheme;
SLEPC_EXTERN const char *NEPRefineSchemes[];

/*E
    NEPConv - Determines the convergence test

    Level: intermediate

.seealso: NEPSetConvergenceTest(), NEPSetConvergenceTestFunction()
E*/
typedef enum { NEP_CONV_ABS,
               NEP_CONV_REL,
               NEP_CONV_NORM,
               NEP_CONV_USER } NEPConv;

/*E
    NEPStop - Determines the stopping test

    Level: advanced

.seealso: NEPSetStoppingTest(), NEPSetStoppingTestFunction()
E*/
typedef enum { NEP_STOP_BASIC,
               NEP_STOP_USER } NEPStop;

/*E
    NEPConvergedReason - Reason a nonlinear eigensolver was said to
         have converged or diverged

    Level: intermediate

.seealso: NEPSolve(), NEPGetConvergedReason(), NEPSetTolerances()
E*/
typedef enum {/* converged */
              NEP_CONVERGED_TOL                =  1,
              NEP_CONVERGED_USER               =  2,
              /* diverged */
              NEP_DIVERGED_ITS                 = -1,
              NEP_DIVERGED_BREAKDOWN           = -2,
                    /* unused                  = -3 */
              NEP_DIVERGED_LINEAR_SOLVE        = -4,
              NEP_DIVERGED_SUBSPACE_EXHAUSTED  = -5,
              NEP_CONVERGED_ITERATING          =  0} NEPConvergedReason;
SLEPC_EXTERN const char *const*NEPConvergedReasons;

SLEPC_EXTERN PetscErrorCode NEPCreate(MPI_Comm,NEP*);
SLEPC_EXTERN PetscErrorCode NEPDestroy(NEP*);
SLEPC_EXTERN PetscErrorCode NEPReset(NEP);
SLEPC_EXTERN PetscErrorCode NEPSetType(NEP,NEPType);
SLEPC_EXTERN PetscErrorCode NEPGetType(NEP,NEPType*);
SLEPC_EXTERN PetscErrorCode NEPSetProblemType(NEP,NEPProblemType);
SLEPC_EXTERN PetscErrorCode NEPGetProblemType(NEP,NEPProblemType*);
SLEPC_EXTERN PetscErrorCode NEPSetTarget(NEP,PetscScalar);
SLEPC_EXTERN PetscErrorCode NEPGetTarget(NEP,PetscScalar*);
SLEPC_EXTERN PetscErrorCode NEPSetFromOptions(NEP);
SLEPC_EXTERN PetscErrorCode NEPSetDSType(NEP);
SLEPC_EXTERN PetscErrorCode NEPSetUp(NEP);
SLEPC_EXTERN PetscErrorCode NEPSolve(NEP);
SLEPC_EXTERN PetscErrorCode NEPView(NEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPViewFromOptions(NEP,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode NEPErrorView(NEP,NEPErrorType,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPErrorViewFromOptions(NEP);
SLEPC_EXTERN PetscErrorCode NEPConvergedReasonView(NEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPConvergedReasonViewFromOptions(NEP);
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "NEPConvergedReasonView()", ) static inline PetscErrorCode NEPReasonView(NEP nep,PetscViewer v) {return NEPConvergedReasonView(nep,v);}
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "NEPConvergedReasonViewFromOptions()", ) static inline PetscErrorCode NEPReasonViewFromOptions(NEP nep) {return NEPConvergedReasonViewFromOptions(nep);}
SLEPC_EXTERN PetscErrorCode NEPValuesView(NEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPValuesViewFromOptions(NEP);
SLEPC_EXTERN PetscErrorCode NEPVectorsView(NEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPVectorsViewFromOptions(NEP);

/*S
  NEPFunctionFn - A prototype of a NEP function evaluation function that would be passed to NEPSetFunction()

  Calling Sequence:
+   nep    - eigensolver context obtained from NEPCreate()
.   lambda - the scalar argument where T(.) must be evaluated
.   T      - matrix that will contain T(lambda)
.   P      - [optional] different matrix to build the preconditioner
-   ctx    - [optional] user-defined context for private data for the
             function evaluation routine (may be NULL)

  Level: beginner

.seealso: NEPSetFunction()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(NEPFunctionFn)(NEP nep,PetscScalar lambda,Mat T,Mat P,void *ctx);

/*S
  NEPJacobianFn - A prototype of a NEP Jacobian evaluation function that would be passed to NEPSetJacobian()

  Calling Sequence:
+   nep    - eigensolver context obtained from NEPCreate()
.   lambda - the scalar argument where T'(.) must be evaluated
.   J      - matrix that will contain T'(lambda)
-   ctx    - [optional] user-defined context for private data for the
             Jacobian evaluation routine (may be NULL)

  Level: beginner

.seealso: NEPSetJacobian()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(NEPJacobianFn)(NEP nep,PetscScalar lambda,Mat J,void *ctx);

SLEPC_EXTERN PetscErrorCode NEPSetFunction(NEP,Mat,Mat,NEPFunctionFn*,void*);
SLEPC_EXTERN PetscErrorCode NEPGetFunction(NEP,Mat*,Mat*,NEPFunctionFn**,void**);
SLEPC_EXTERN PetscErrorCode NEPSetJacobian(NEP,Mat,NEPJacobianFn*,void*);
SLEPC_EXTERN PetscErrorCode NEPGetJacobian(NEP,Mat*,NEPJacobianFn**,void**);
PETSC_DEPRECATED_FUNCTION(3, 12, 0, "NEPSetFunction() and NEPSetJacobian()", ) static inline PetscErrorCode NEPSetDerivatives(NEP nep,Mat A,PetscErrorCode (*fun)(NEP,PetscScalar,PetscInt,Mat,void*),void *ctx) {(void)A;(void)fun;(void)ctx;SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Not implemented in this version");}
PETSC_DEPRECATED_FUNCTION(3, 12, 0, "NEPGetFunction() and NEPGetJacobian()", ) static inline PetscErrorCode NEPGetDerivatives(NEP nep,Mat *A,PetscErrorCode (**fun)(NEP,PetscScalar,PetscInt,Mat,void*),void **ctx) {(void)A;(void)fun;(void)ctx;SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Not implemented in this version");}
SLEPC_EXTERN PetscErrorCode NEPSetSplitOperator(NEP,PetscInt,Mat[],FN[],MatStructure);
SLEPC_EXTERN PetscErrorCode NEPGetSplitOperatorTerm(NEP,PetscInt,Mat*,FN*);
SLEPC_EXTERN PetscErrorCode NEPGetSplitOperatorInfo(NEP,PetscInt*,MatStructure*);
SLEPC_EXTERN PetscErrorCode NEPSetSplitPreconditioner(NEP,PetscInt,Mat[],MatStructure);
SLEPC_EXTERN PetscErrorCode NEPGetSplitPreconditionerTerm(NEP,PetscInt,Mat*);
SLEPC_EXTERN PetscErrorCode NEPGetSplitPreconditionerInfo(NEP,PetscInt*,MatStructure*);

SLEPC_EXTERN PetscErrorCode NEPSetBV(NEP,BV);
SLEPC_EXTERN PetscErrorCode NEPGetBV(NEP,BV*);
SLEPC_EXTERN PetscErrorCode NEPSetRG(NEP,RG);
SLEPC_EXTERN PetscErrorCode NEPGetRG(NEP,RG*);
SLEPC_EXTERN PetscErrorCode NEPSetDS(NEP,DS);
SLEPC_EXTERN PetscErrorCode NEPGetDS(NEP,DS*);
SLEPC_EXTERN PetscErrorCode NEPRefineGetKSP(NEP,KSP*);
SLEPC_EXTERN PetscErrorCode NEPSetTolerances(NEP,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPGetTolerances(NEP,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPSetConvergenceTest(NEP,NEPConv);
SLEPC_EXTERN PetscErrorCode NEPGetConvergenceTest(NEP,NEPConv*);
SLEPC_EXTERN PetscErrorCode NEPConvergedAbsolute(NEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode NEPConvergedRelative(NEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode NEPConvergedNorm(NEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode NEPSetStoppingTest(NEP,NEPStop);
SLEPC_EXTERN PetscErrorCode NEPGetStoppingTest(NEP,NEPStop*);
SLEPC_EXTERN PetscErrorCode NEPStoppingBasic(NEP,PetscInt,PetscInt,PetscInt,PetscInt,NEPConvergedReason*,void*);
SLEPC_EXTERN PetscErrorCode NEPSetDimensions(NEP,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPGetDimensions(NEP,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPSetRefine(NEP,NEPRefine,PetscInt,PetscReal,PetscInt,NEPRefineScheme);
SLEPC_EXTERN PetscErrorCode NEPGetRefine(NEP,NEPRefine*,PetscInt*,PetscReal*,PetscInt*,NEPRefineScheme*);

SLEPC_EXTERN PetscErrorCode NEPGetConverged(NEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPGetEigenpair(NEP,PetscInt,PetscScalar*,PetscScalar*,Vec,Vec);
SLEPC_EXTERN PetscErrorCode NEPGetLeftEigenvector(NEP,PetscInt,Vec,Vec);

SLEPC_EXTERN PetscErrorCode NEPComputeError(NEP,PetscInt,NEPErrorType,PetscReal*);
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "NEPComputeError()", ) static inline PetscErrorCode NEPComputeRelativeError(NEP nep,PetscInt i,PetscReal *r) {return NEPComputeError(nep,i,NEP_ERROR_RELATIVE,r);}
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "NEPComputeError() with NEP_ERROR_ABSOLUTE", ) static inline PetscErrorCode NEPComputeResidualNorm(NEP nep,PetscInt i,PetscReal *r) {return NEPComputeError(nep,i,NEP_ERROR_ABSOLUTE,r);}
SLEPC_EXTERN PetscErrorCode NEPGetErrorEstimate(NEP,PetscInt,PetscReal*);

SLEPC_EXTERN PetscErrorCode NEPComputeFunction(NEP,PetscScalar,Mat,Mat);
SLEPC_EXTERN PetscErrorCode NEPComputeJacobian(NEP,PetscScalar,Mat);
SLEPC_EXTERN PetscErrorCode NEPApplyFunction(NEP,PetscScalar,Vec,Vec,Vec,Mat,Mat);
SLEPC_EXTERN PetscErrorCode NEPApplyAdjoint(NEP,PetscScalar,Vec,Vec,Vec,Mat,Mat);
SLEPC_EXTERN PetscErrorCode NEPApplyJacobian(NEP,PetscScalar,Vec,Vec,Vec,Mat);
SLEPC_EXTERN PetscErrorCode NEPProjectOperator(NEP,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPGetIterationNumber(NEP,PetscInt*);

SLEPC_EXTERN PetscErrorCode NEPSetInitialSpace(NEP,PetscInt,Vec[]);
SLEPC_EXTERN PetscErrorCode NEPSetWhichEigenpairs(NEP,NEPWhich);
SLEPC_EXTERN PetscErrorCode NEPGetWhichEigenpairs(NEP,NEPWhich*);
SLEPC_EXTERN PetscErrorCode NEPSetTwoSided(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPGetTwoSided(NEP,PetscBool*);

SLEPC_EXTERN PetscErrorCode NEPApplyResolvent(NEP,RG,PetscScalar,Vec,Vec);

SLEPC_EXTERN PetscErrorCode NEPSetTrackAll(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPGetTrackAll(NEP,PetscBool*);

SLEPC_EXTERN PetscErrorCode NEPGetConvergedReason(NEP,NEPConvergedReason*);

SLEPC_EXTERN PetscErrorCode NEPMonitor(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPMonitorSet(NEP,PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*),void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode NEPMonitorCancel(NEP);
SLEPC_EXTERN PetscErrorCode NEPGetMonitorContext(NEP,void*);

SLEPC_EXTERN PetscErrorCode NEPMonitorSetFromOptions(NEP,const char[],const char[],void*,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPMonitorFirst(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode NEPMonitorFirstDrawLG(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode NEPMonitorFirstDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode NEPMonitorAll(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode NEPMonitorAllDrawLG(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode NEPMonitorAllDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode NEPMonitorConverged(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode NEPMonitorConvergedCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode NEPMonitorConvergedDrawLG(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode NEPMonitorConvergedDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode NEPMonitorConvergedDestroy(PetscViewerAndFormat**);

SLEPC_EXTERN PetscErrorCode NEPSetOptionsPrefix(NEP,const char*);
SLEPC_EXTERN PetscErrorCode NEPAppendOptionsPrefix(NEP,const char*);
SLEPC_EXTERN PetscErrorCode NEPGetOptionsPrefix(NEP,const char*[]);

SLEPC_EXTERN PetscFunctionList NEPList;
SLEPC_EXTERN PetscFunctionList NEPMonitorList;
SLEPC_EXTERN PetscFunctionList NEPMonitorCreateList;
SLEPC_EXTERN PetscFunctionList NEPMonitorDestroyList;
SLEPC_EXTERN PetscErrorCode NEPRegister(const char[],PetscErrorCode(*)(NEP));
SLEPC_EXTERN PetscErrorCode NEPMonitorRegister(const char[],PetscViewerType,PetscViewerFormat,PetscErrorCode(*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscErrorCode(*)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**),PetscErrorCode(*)(PetscViewerAndFormat**));

SLEPC_EXTERN PetscErrorCode NEPSetWorkVecs(NEP,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPAllocateSolution(NEP,PetscInt);

/*S
  NEPConvergenceTestFn - A prototype of a NEP convergence test function that would be passed to NEPSetConvergenceTestFunction()

  Calling Sequence:
+   nep    - eigensolver context obtained from NEPCreate()
.   eigr   - real part of the eigenvalue
.   eigi   - imaginary part of the eigenvalue
.   res    - residual norm associated to the eigenpair
.   errest - [output] computed error estimate
-   ctx    - [optional] user-defined context for private data for the
             convergence test routine (may be NULL)

  Level: advanced

.seealso: NEPSetConvergenceTestFunction()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(NEPConvergenceTestFn)(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx);

/*S
  NEPStoppingTestFn - A prototype of a NEP stopping test function that would be passed to NEPSetStoppingTestFunction()

  Calling Sequence:
+   nep    - eigensolver context obtained from NEPCreate()
.   its    - current number of iterations
.   max_it - maximum number of iterations
.   nconv  - number of currently converged eigenpairs
.   nev    - number of requested eigenpairs
.   reason - [output] result of the stopping test
-   ctx    - [optional] user-defined context for private data for the
             stopping test routine (may be NULL)

  Level: advanced

.seealso: NEPSetStoppingTestFunction()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(NEPStoppingTestFn)(NEP nep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,NEPConvergedReason *reason,void *ctx);

SLEPC_EXTERN PetscErrorCode NEPSetConvergenceTestFunction(NEP,NEPConvergenceTestFn*,void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode NEPSetStoppingTestFunction(NEP,NEPStoppingTestFn*,void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode NEPSetEigenvalueComparison(NEP,SlepcEigenvalueComparisonFn*,void*);

/* --------- options specific to particular eigensolvers -------- */

SLEPC_EXTERN PetscErrorCode NEPRIISetMaximumIterations(NEP,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPRIIGetMaximumIterations(NEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPRIISetLagPreconditioner(NEP,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPRIIGetLagPreconditioner(NEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPRIISetConstCorrectionTol(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPRIIGetConstCorrectionTol(NEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPRIISetHermitian(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPRIIGetHermitian(NEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPRIISetDeflationThreshold(NEP,PetscReal);
SLEPC_EXTERN PetscErrorCode NEPRIIGetDeflationThreshold(NEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode NEPRIISetKSP(NEP,KSP);
SLEPC_EXTERN PetscErrorCode NEPRIIGetKSP(NEP,KSP*);

SLEPC_EXTERN PetscErrorCode NEPSLPSetDeflationThreshold(NEP,PetscReal);
SLEPC_EXTERN PetscErrorCode NEPSLPGetDeflationThreshold(NEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode NEPSLPSetEPS(NEP,EPS);
SLEPC_EXTERN PetscErrorCode NEPSLPGetEPS(NEP,EPS*);
SLEPC_EXTERN PetscErrorCode NEPSLPSetEPSLeft(NEP,EPS);
SLEPC_EXTERN PetscErrorCode NEPSLPGetEPSLeft(NEP,EPS*);
SLEPC_EXTERN PetscErrorCode NEPSLPSetKSP(NEP,KSP);
SLEPC_EXTERN PetscErrorCode NEPSLPGetKSP(NEP,KSP*);

SLEPC_EXTERN PetscErrorCode NEPNArnoldiSetKSP(NEP,KSP);
SLEPC_EXTERN PetscErrorCode NEPNArnoldiGetKSP(NEP,KSP*);
SLEPC_EXTERN PetscErrorCode NEPNArnoldiSetLagPreconditioner(NEP,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPNArnoldiGetLagPreconditioner(NEP,PetscInt*);

/*E
    NEPCISSExtraction - determines the extraction technique in the CISS solver

    Level: advanced

.seealso: NEPCISSSetExtraction(), NEPCISSGetExtraction()
E*/
typedef enum { NEP_CISS_EXTRACTION_RITZ,
               NEP_CISS_EXTRACTION_HANKEL,
               NEP_CISS_EXTRACTION_CAA    } NEPCISSExtraction;
SLEPC_EXTERN const char *NEPCISSExtractions[];

#if defined(PETSC_USE_COMPLEX) || defined(PETSC_CLANG_STATIC_ANALYZER)
SLEPC_EXTERN PetscErrorCode NEPCISSSetExtraction(NEP,NEPCISSExtraction);
SLEPC_EXTERN PetscErrorCode NEPCISSGetExtraction(NEP,NEPCISSExtraction*);
SLEPC_EXTERN PetscErrorCode NEPCISSSetSizes(NEP,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPCISSGetSizes(NEP,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPCISSSetThreshold(NEP,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode NEPCISSGetThreshold(NEP,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode NEPCISSSetRefinement(NEP,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPCISSGetRefinement(NEP,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPCISSGetKSPs(NEP,PetscInt*,KSP*[]);
#else
#define SlepcNEPCISSUnavailable(nep) do { \
    PetscFunctionBegin; \
    SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"%s() not available with real scalars",PETSC_FUNCTION_NAME); \
    } while (0)
static inline PetscErrorCode NEPCISSSetExtraction(NEP nep,PETSC_UNUSED NEPCISSExtraction ex) {SlepcNEPCISSUnavailable(nep);}
static inline PetscErrorCode NEPCISSGetExtraction(NEP nep,PETSC_UNUSED NEPCISSExtraction *ex) {SlepcNEPCISSUnavailable(nep);}
static inline PetscErrorCode NEPCISSSetSizes(NEP nep,PETSC_UNUSED PetscInt ip,PETSC_UNUSED PetscInt bs,PETSC_UNUSED PetscInt ms,PETSC_UNUSED PetscInt npart,PETSC_UNUSED PetscInt bsmax,PETSC_UNUSED PetscBool realmats) {SlepcNEPCISSUnavailable(nep);}
static inline PetscErrorCode NEPCISSGetSizes(NEP nep,PETSC_UNUSED PetscInt *ip,PETSC_UNUSED PetscInt *bs,PETSC_UNUSED PetscInt *ms,PETSC_UNUSED PetscInt *npart,PETSC_UNUSED PetscInt *bsmak,PETSC_UNUSED PetscBool *realmats) {SlepcNEPCISSUnavailable(nep);}
static inline PetscErrorCode NEPCISSSetThreshold(NEP nep,PETSC_UNUSED PetscReal delta,PETSC_UNUSED PetscReal spur) {SlepcNEPCISSUnavailable(nep);}
static inline PetscErrorCode NEPCISSGetThreshold(NEP nep,PETSC_UNUSED PetscReal *delta,PETSC_UNUSED PetscReal *spur) {SlepcNEPCISSUnavailable(nep);}
static inline PetscErrorCode NEPCISSSetRefinement(NEP nep,PETSC_UNUSED PetscInt inner,PETSC_UNUSED PetscInt blsize) {SlepcNEPCISSUnavailable(nep);}
static inline PetscErrorCode NEPCISSGetRefinement(NEP nep,PETSC_UNUSED PetscInt *inner,PETSC_UNUSED PetscInt *blsize) {SlepcNEPCISSUnavailable(nep);}
static inline PetscErrorCode NEPCISSGetKSPs(NEP nep,PETSC_UNUSED PetscInt *nsolve,PETSC_UNUSED KSP *ksp[]) {SlepcNEPCISSUnavailable(nep);}
#undef SlepcNEPCISSUnavailable
#endif

SLEPC_EXTERN PetscErrorCode NEPInterpolSetPEP(NEP,PEP);
SLEPC_EXTERN PetscErrorCode NEPInterpolGetPEP(NEP,PEP*);
SLEPC_EXTERN PetscErrorCode NEPInterpolSetInterpolation(NEP,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPInterpolGetInterpolation(NEP,PetscReal*,PetscInt*);

SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetSingularitiesFunction(NEP,PetscErrorCode (*)(NEP,PetscInt*,PetscScalar*,void*),void*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetSingularitiesFunction(NEP,PetscErrorCode (**)(NEP,PetscInt*,PetscScalar*,void*),void**);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetRestart(NEP,PetscReal);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetRestart(NEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetLocking(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetLocking(NEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetInterpolation(NEP,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetInterpolation(NEP,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetRKShifts(NEP,PetscInt,PetscScalar[]);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetRKShifts(NEP,PetscInt*,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetKSPs(NEP,PetscInt*,KSP*[]);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetFullBasis(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetFullBasis(NEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetEPS(NEP,EPS);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetEPS(NEP,EPS*);
