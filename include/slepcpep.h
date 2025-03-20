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
     PEP - Abstract SLEPc object that manages all the polynomial eigenvalue
     problem solvers.

   Level: beginner

.seealso:  PEPCreate()
S*/
typedef struct _p_PEP* PEP;

/*J
    PEPType - String with the name of a polynomial eigensolver

   Level: beginner

.seealso: PEPSetType(), PEP
J*/
typedef const char *PEPType;
#define PEPLINEAR    "linear"
#define PEPQARNOLDI  "qarnoldi"
#define PEPTOAR      "toar"
#define PEPSTOAR     "stoar"
#define PEPJD        "jd"
#define PEPCISS      "ciss"

/* Logging support */
SLEPC_EXTERN PetscClassId PEP_CLASSID;

/*E
    PEPProblemType - Determines the type of the polynomial eigenproblem

    PEP_HERMITIAN is used when all A_i are Hermitian,
    PEP_HYPERBOLIC is reserved for a QEP with Hermitian matrices, M>0, (x'Cx)^2 > 4(x'Mx)(x'Kx),
    PEP_GYROSCOPIC is for aQEP with M, K  Hermitian, M>0, C skew-Hermitian.

    Level: intermediate

.seealso: PEPSetProblemType(), PEPGetProblemType()
E*/
typedef enum { PEP_GENERAL    = 1,
               PEP_HERMITIAN  = 2,
               PEP_HYPERBOLIC = 3,
               PEP_GYROSCOPIC = 4
             } PEPProblemType;

/*E
    PEPWhich - Determines which part of the spectrum is requested

    Level: intermediate

.seealso: PEPSetWhichEigenpairs(), PEPGetWhichEigenpairs()
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
    eigenproblem

    Level: intermediate

.seealso: PEPSetBasis()
E*/
typedef enum { PEP_BASIS_MONOMIAL,
               PEP_BASIS_CHEBYSHEV1,
               PEP_BASIS_CHEBYSHEV2,
               PEP_BASIS_LEGENDRE,
               PEP_BASIS_LAGUERRE,
               PEP_BASIS_HERMITE } PEPBasis;
SLEPC_EXTERN const char *PEPBasisTypes[];

/*E
    PEPScale - The scaling strategy

    Level: intermediate

.seealso: PEPSetScale()
E*/
typedef enum { PEP_SCALE_NONE,
               PEP_SCALE_SCALAR,
               PEP_SCALE_DIAGONAL,
               PEP_SCALE_BOTH } PEPScale;
SLEPC_EXTERN const char *PEPScaleTypes[];

/*E
    PEPRefine - The refinement type

    Level: intermediate

.seealso: PEPSetRefine()
E*/
typedef enum { PEP_REFINE_NONE,
               PEP_REFINE_SIMPLE,
               PEP_REFINE_MULTIPLE } PEPRefine;
SLEPC_EXTERN const char *PEPRefineTypes[];

/*E
    PEPRefineScheme - The scheme used for solving linear systems during iterative refinement

    Level: intermediate

.seealso: PEPSetRefine()
E*/
typedef enum { PEP_REFINE_SCHEME_SCHUR    = 1,
               PEP_REFINE_SCHEME_MBE      = 2,
               PEP_REFINE_SCHEME_EXPLICIT = 3 } PEPRefineScheme;
SLEPC_EXTERN const char *PEPRefineSchemes[];

/*E
    PEPExtract - The extraction type

    Level: intermediate

.seealso: PEPSetExtract()
E*/
typedef enum { PEP_EXTRACT_NONE       = 1,
               PEP_EXTRACT_NORM       = 2,
               PEP_EXTRACT_RESIDUAL   = 3,
               PEP_EXTRACT_STRUCTURED = 4 } PEPExtract;
SLEPC_EXTERN const char *PEPExtractTypes[];

/*E
    PEPErrorType - The error type used to assess accuracy of computed solutions

    Level: intermediate

.seealso: PEPComputeError()
E*/
typedef enum { PEP_ERROR_ABSOLUTE,
               PEP_ERROR_RELATIVE,
               PEP_ERROR_BACKWARD } PEPErrorType;
SLEPC_EXTERN const char *PEPErrorTypes[];

/*E
    PEPConv - Determines the convergence test

    Level: intermediate

.seealso: PEPSetConvergenceTest(), PEPSetConvergenceTestFunction()
E*/
typedef enum { PEP_CONV_ABS,
               PEP_CONV_REL,
               PEP_CONV_NORM,
               PEP_CONV_USER } PEPConv;

/*E
    PEPStop - Determines the stopping test

    Level: advanced

.seealso: PEPSetStoppingTest(), PEPSetStoppingTestFunction()
E*/
typedef enum { PEP_STOP_BASIC,
               PEP_STOP_USER } PEPStop;

/*E
    PEPConvergedReason - Reason an eigensolver was said to
         have converged or diverged

    Level: intermediate

.seealso: PEPSolve(), PEPGetConvergedReason(), PEPSetTolerances()
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
SLEPC_EXTERN PetscErrorCode PEPSetConvergenceTest(PEP,PEPConv);
SLEPC_EXTERN PetscErrorCode PEPGetConvergenceTest(PEP,PEPConv*);
SLEPC_EXTERN PetscErrorCode PEPConvergedAbsolute(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode PEPConvergedRelative(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode PEPConvergedNorm(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode PEPSetStoppingTest(PEP,PEPStop);
SLEPC_EXTERN PetscErrorCode PEPGetStoppingTest(PEP,PEPStop*);
SLEPC_EXTERN PetscErrorCode PEPStoppingBasic(PEP,PetscInt,PetscInt,PetscInt,PetscInt,PEPConvergedReason*,void*);
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

SLEPC_EXTERN PetscErrorCode PEPMonitor(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt);
SLEPC_EXTERN PetscErrorCode PEPMonitorSet(PEP,PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*),void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode PEPMonitorCancel(PEP);
SLEPC_EXTERN PetscErrorCode PEPGetMonitorContext(PEP,void*);

SLEPC_EXTERN PetscErrorCode PEPMonitorSetFromOptions(PEP,const char[],const char[],void*,PetscBool);
SLEPC_EXTERN PetscErrorCode PEPMonitorFirst(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode PEPMonitorFirstDrawLG(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode PEPMonitorFirstDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode PEPMonitorAll(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode PEPMonitorAllDrawLG(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode PEPMonitorAllDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode PEPMonitorConverged(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode PEPMonitorConvergedCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode PEPMonitorConvergedDrawLG(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode PEPMonitorConvergedDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode PEPMonitorConvergedDestroy(PetscViewerAndFormat**);

SLEPC_EXTERN PetscErrorCode PEPSetOptionsPrefix(PEP,const char*);
SLEPC_EXTERN PetscErrorCode PEPAppendOptionsPrefix(PEP,const char*);
SLEPC_EXTERN PetscErrorCode PEPGetOptionsPrefix(PEP,const char*[]);

SLEPC_EXTERN PetscFunctionList PEPList;
SLEPC_EXTERN PetscFunctionList PEPMonitorList;
SLEPC_EXTERN PetscFunctionList PEPMonitorCreateList;
SLEPC_EXTERN PetscFunctionList PEPMonitorDestroyList;
SLEPC_EXTERN PetscErrorCode PEPRegister(const char[],PetscErrorCode(*)(PEP));
SLEPC_EXTERN PetscErrorCode PEPMonitorRegister(const char[],PetscViewerType,PetscViewerFormat,PetscErrorCode(*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscErrorCode(*)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**),PetscErrorCode(*)(PetscViewerAndFormat**));

SLEPC_EXTERN PetscErrorCode PEPSetWorkVecs(PEP,PetscInt);
SLEPC_EXTERN PetscErrorCode PEPAllocateSolution(PEP,PetscInt);

/*S
  PEPConvergenceTestFn - A prototype of a PEP convergence test function that would be passed to PEPSetConvergenceTestFunction()

  Calling Sequence:
+   pep    - eigensolver context obtained from PEPCreate()
.   eigr   - real part of the eigenvalue
.   eigi   - imaginary part of the eigenvalue
.   res    - residual norm associated to the eigenpair
.   errest - [output] computed error estimate
-   ctx    - [optional] user-defined context for private data for the
             convergence test routine (may be NULL)

  Level: advanced

.seealso: PEPSetConvergenceTestFunction()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(PEPConvergenceTestFn)(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx);

/*S
  PEPStoppingTestFn - A prototype of a PEP stopping test function that would be passed to PEPSetStoppingTestFunction()

  Calling Sequence:
+   pep    - eigensolver context obtained from PEPCreate()
.   its    - current number of iterations
.   max_it - maximum number of iterations
.   nconv  - number of currently converged eigenpairs
.   nev    - number of requested eigenpairs
.   reason - [output] result of the stopping test
-   ctx    - [optional] user-defined context for private data for the
             stopping test routine (may be NULL)

  Level: advanced

.seealso: PEPSetStoppingTestFunction()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(PEPStoppingTestFn)(PEP pep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,PEPConvergedReason *reason,void *ctx);

SLEPC_EXTERN PetscErrorCode PEPSetConvergenceTestFunction(PEP,PEPConvergenceTestFn*,void*,PetscCtxDestroyFn*);
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
    PEPJDProjection - The type of projection to be used in the Jacobi-Davidson solver

    Level: intermediate

.seealso: PEPJDSetProjection()
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
    PEPCISSExtraction - determines the extraction technique in the CISS solver

    Level: advanced

.seealso: PEPCISSSetExtraction(), PEPCISSGetExtraction()
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
