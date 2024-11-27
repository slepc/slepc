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
     SVD - Abstract SLEPc object that manages all the singular value
     problem solvers.

   Level: beginner

.seealso:  SVDCreate()
S*/
typedef struct _p_SVD* SVD;

/*J
    SVDType - String with the name of a SLEPc singular value solver

   Level: beginner

.seealso: SVDSetType(), SVD
J*/
typedef const char* SVDType;
#define SVDCROSS       "cross"
#define SVDCYCLIC      "cyclic"
#define SVDLAPACK      "lapack"
#define SVDLANCZOS     "lanczos"
#define SVDTRLANCZOS   "trlanczos"
#define SVDRANDOMIZED  "randomized"
#define SVDSCALAPACK   "scalapack"
#define SVDKSVD        "ksvd"
#define SVDELEMENTAL   "elemental"
#define SVDPRIMME      "primme"

/* Logging support */
SLEPC_EXTERN PetscClassId SVD_CLASSID;

/*E
    SVDProblemType - Determines the type of singular value problem

    Level: beginner

.seealso: SVDSetProblemType(), SVDGetProblemType()
E*/
typedef enum { SVD_STANDARD=1,
               SVD_GENERALIZED,   /* GSVD */
               SVD_HYPERBOLIC     /* HSVD */
             } SVDProblemType;

/*E
    SVDWhich - Determines whether largest or smallest singular triplets
    are to be computed

    Level: intermediate

.seealso: SVDSetWhichSingularTriplets(), SVDGetWhichSingularTriplets()
E*/
typedef enum { SVD_LARGEST,
               SVD_SMALLEST } SVDWhich;

/*E
    SVDErrorType - The error type used to assess accuracy of computed solutions

    Level: intermediate

.seealso: SVDComputeError()
E*/
typedef enum { SVD_ERROR_ABSOLUTE,
               SVD_ERROR_RELATIVE,
               SVD_ERROR_NORM } SVDErrorType;
SLEPC_EXTERN const char *SVDErrorTypes[];

/*E
    SVDConv - Determines the convergence test

    Level: intermediate

.seealso: SVDSetConvergenceTest(), SVDSetConvergenceTestFunction()
E*/
typedef enum { SVD_CONV_ABS,
               SVD_CONV_REL,
               SVD_CONV_NORM,
               SVD_CONV_MAXIT,
               SVD_CONV_USER } SVDConv;

/*E
    SVDStop - Determines the stopping test

    Level: advanced

.seealso: SVDSetStoppingTest(), SVDSetStoppingTestFunction()
E*/
typedef enum { SVD_STOP_BASIC,
               SVD_STOP_USER,
               SVD_STOP_THRESHOLD } SVDStop;

/*E
    SVDConvergedReason - Reason a singular value solver was said to
         have converged or diverged

   Level: intermediate

.seealso: SVDSolve(), SVDGetConvergedReason(), SVDSetTolerances()
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

/*S
   SVDStoppingCtx - Data structure (C struct) to hold additional information to
   be used in some stopping test functions.

   Level: advanced

.seealso: SVDSetStoppingTestFunction()
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
SLEPC_EXTERN PetscErrorCode SVDSetOptionsPrefix(SVD,const char*);
SLEPC_EXTERN PetscErrorCode SVDAppendOptionsPrefix(SVD,const char*);
SLEPC_EXTERN PetscErrorCode SVDGetOptionsPrefix(SVD,const char*[]);
SLEPC_EXTERN PetscErrorCode SVDSetDSType(SVD);
SLEPC_EXTERN PetscErrorCode SVDSetUp(SVD);
SLEPC_EXTERN PetscErrorCode SVDSolve(SVD);
SLEPC_EXTERN PetscErrorCode SVDGetIterationNumber(SVD,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDSetConvergenceTest(SVD,SVDConv);
SLEPC_EXTERN PetscErrorCode SVDGetConvergenceTest(SVD,SVDConv*);
SLEPC_EXTERN PetscErrorCode SVDConvergedAbsolute(SVD,PetscReal,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode SVDConvergedRelative(SVD,PetscReal,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode SVDConvergedNorm(SVD,PetscReal,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode SVDConvergedMaxIt(SVD,PetscReal,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode SVDSetStoppingTest(SVD,SVDStop);
SLEPC_EXTERN PetscErrorCode SVDGetStoppingTest(SVD,SVDStop*);
SLEPC_EXTERN PetscErrorCode SVDStoppingBasic(SVD,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*);
SLEPC_EXTERN PetscErrorCode SVDStoppingThreshold(SVD,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*);
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

SLEPC_EXTERN PetscErrorCode SVDMonitor(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDMonitorSet(SVD,PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*),void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode SVDMonitorCancel(SVD);
SLEPC_EXTERN PetscErrorCode SVDGetMonitorContext(SVD,void*);

SLEPC_EXTERN PetscErrorCode SVDMonitorSetFromOptions(SVD,const char[],const char[],void*,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDMonitorLGCreate(MPI_Comm,const char[],const char[],const char[],PetscInt,const char*[],int,int,int,int,PetscDrawLG*);
SLEPC_EXTERN PetscErrorCode SVDMonitorFirst(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode SVDMonitorFirstDrawLG(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode SVDMonitorFirstDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode SVDMonitorAll(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode SVDMonitorAllDrawLG(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode SVDMonitorAllDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode SVDMonitorConverged(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode SVDMonitorConvergedCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode SVDMonitorConvergedDrawLG(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode SVDMonitorConvergedDrawLGCreate(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode SVDMonitorConvergedDestroy(PetscViewerAndFormat**);
SLEPC_EXTERN PetscErrorCode SVDMonitorConditioning(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);

SLEPC_EXTERN PetscFunctionList SVDList;
SLEPC_EXTERN PetscFunctionList SVDMonitorList;
SLEPC_EXTERN PetscFunctionList SVDMonitorCreateList;
SLEPC_EXTERN PetscFunctionList SVDMonitorDestroyList;
SLEPC_EXTERN PetscErrorCode SVDRegister(const char[],PetscErrorCode(*)(SVD));
SLEPC_EXTERN PetscErrorCode SVDMonitorRegister(const char[],PetscViewerType,PetscViewerFormat,PetscErrorCode(*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscErrorCode(*)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**),PetscErrorCode(*)(PetscViewerAndFormat**));

SLEPC_EXTERN PetscErrorCode SVDAllocateSolution(SVD,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDReallocateSolution(SVD,PetscInt);

/*S
  SVDConvergenceTestFn - A prototype of an SVD convergence test function that would be passed to SVDSetConvergenceTestFunction()

  Calling Sequence:
+   svd    - singular value solver context obtained from SVDCreate()
.   sigma  - computed singular value
.   res    - residual norm associated to the singular triplet
.   errest - [output] computed error estimate
-   ctx    - [optional] user-defined context for private data for the
             convergence test routine (may be NULL)

  Level: advanced

.seealso: SVDSetConvergenceTestFunction()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(SVDConvergenceTestFn)(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx);

/*S
  SVDStoppingTestFn - A prototype of an SVD stopping test function that would be passed to SVDSetStoppingTestFunction()

  Calling Sequence:
+   svd    - singular value solver context obtained from SVDCreate()
.   its    - current number of iterations
.   max_it - maximum number of iterations
.   nconv  - number of currently converged singular triplets
.   nsv    - number of requested singular triplets
.   reason - [output] result of the stopping test
-   ctx    - [optional] user-defined context for private data for the
             stopping test routine (may be NULL)

  Level: advanced

.seealso: SVDSetStoppingTestFunction()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(SVDStoppingTestFn)(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nsv,SVDConvergedReason *reason,void *ctx);

SLEPC_EXTERN PetscErrorCode SVDSetConvergenceTestFunction(SVD,SVDConvergenceTestFn*,void*,PetscCtxDestroyFn*);
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
    SVDTRLanczosGBidiag - determines the bidiagonalization choice for the
    TRLanczos GSVD solver

    Level: advanced

.seealso: SVDTRLanczosSetGBidiag(), SVDTRLanczosGetGBidiag()
E*/
typedef enum {
  SVD_TRLANCZOS_GBIDIAG_SINGLE, /* single bidiagonalization (Qa) */
  SVD_TRLANCZOS_GBIDIAG_UPPER,  /* joint bidiagonalization, both Qa and Qb in upper bidiagonal form */
  SVD_TRLANCZOS_GBIDIAG_LOWER   /* joint bidiagonalization, Qa lower bidiagonal, Qb upper bidiagonal */
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
    SVDPRIMMEMethod - determines the SVD method selected in the PRIMME library

    Level: advanced

.seealso: SVDPRIMMESetMethod(), SVDPRIMMEGetMethod()
E*/
typedef enum { SVD_PRIMME_HYBRID=1,
               SVD_PRIMME_NORMALEQUATIONS,
               SVD_PRIMME_AUGMENTED } SVDPRIMMEMethod;
SLEPC_EXTERN const char *SVDPRIMMEMethods[];

SLEPC_EXTERN PetscErrorCode SVDPRIMMESetBlockSize(SVD,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDPRIMMEGetBlockSize(SVD,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDPRIMMESetMethod(SVD,SVDPRIMMEMethod);
SLEPC_EXTERN PetscErrorCode SVDPRIMMEGetMethod(SVD,SVDPRIMMEMethod*);

/*E
    SVDKSVDEigenMethod - determines the method to solve the eigenproblem within KSVD

    Level: advanced

.seealso: SVDKSVDSetEigenMethod(), SVDKSVDGetEigenMethod()
E*/
typedef enum { SVD_KSVD_EIGEN_MRRR=1,
               SVD_KSVD_EIGEN_DC,
               SVD_KSVD_EIGEN_ELPA } SVDKSVDEigenMethod;
SLEPC_EXTERN const char *SVDKSVDEigenMethods[];

/*E
    SVDKSVDPolarMethod - determines the method to compute the polar decomposition within KSVD

    Level: advanced

.seealso: SVDKSVDSetPolarMethod(), SVDKSVDGetPolarMethod()
E*/
typedef enum { SVD_KSVD_POLAR_QDWH=1,
               SVD_KSVD_POLAR_ZOLOPD } SVDKSVDPolarMethod;
SLEPC_EXTERN const char *SVDKSVDPolarMethods[];

SLEPC_EXTERN PetscErrorCode SVDKSVDSetEigenMethod(SVD,SVDKSVDEigenMethod);
SLEPC_EXTERN PetscErrorCode SVDKSVDGetEigenMethod(SVD,SVDKSVDEigenMethod*);
SLEPC_EXTERN PetscErrorCode SVDKSVDSetPolarMethod(SVD,SVDKSVDPolarMethod);
SLEPC_EXTERN PetscErrorCode SVDKSVDGetPolarMethod(SVD,SVDKSVDPolarMethod*);
