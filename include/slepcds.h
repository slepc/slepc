/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the direct solver object in SLEPc
*/

#pragma once

#include <slepcsc.h>
#include <slepcfn.h>
#include <slepcrg.h>

/* SUBMANSEC = DS */

#define DS_MAX_SOLVE 6

SLEPC_EXTERN PetscErrorCode DSInitializePackage(void);
SLEPC_EXTERN PetscErrorCode DSFinalizePackage(void);

/*S
    DS - Direct solver (or dense system), to represent low-dimensional
    eigenproblems that must be solved within iterative solvers. This is an
    auxiliary object and is not normally needed by application programmers.

    Level: beginner

.seealso:  DSCreate()
S*/
typedef struct _p_DS* DS;

/*J
    DSType - String with the name of the type of direct solver. Roughly,
    there are as many types as problem types are available within SLEPc.

    Level: advanced

.seealso: DSSetType(), DS
J*/
typedef const char* DSType;
#define DSHEP    "hep"
#define DSNHEP   "nhep"
#define DSGHEP   "ghep"
#define DSGHIEP  "ghiep"
#define DSGNHEP  "gnhep"
#define DSNHEPTS "nhepts"
#define DSSVD    "svd"
#define DSHSVD   "hsvd"
#define DSGSVD   "gsvd"
#define DSPEP    "pep"
#define DSNEP    "nep"

/* Logging support */
SLEPC_EXTERN PetscClassId DS_CLASSID;

/*E
    DSStateType - Indicates in which state the direct solver is

    Level: advanced

.seealso: DSSetState()
E*/
typedef enum { DS_STATE_RAW,
               DS_STATE_INTERMEDIATE,
               DS_STATE_CONDENSED,
               DS_STATE_TRUNCATED } DSStateType;
SLEPC_EXTERN const char *DSStateTypes[];

/*E
    DSMatType - Used to refer to one of the matrices stored internally in DS

    Notes:
    The matrices preferentially refer to
+   DS_MAT_A  - first matrix of eigenproblem/singular value problem
.   DS_MAT_B  - second matrix of a generalized eigenproblem
.   DS_MAT_C  - third matrix of a quadratic eigenproblem (deprecated)
.   DS_MAT_T  - tridiagonal matrix
.   DS_MAT_D  - diagonal matrix
.   DS_MAT_Q  - orthogonal matrix of (right) Schur vectors
.   DS_MAT_Z  - orthogonal matrix of left Schur vectors
.   DS_MAT_X  - right eigenvectors
.   DS_MAT_Y  - left eigenvectors
.   DS_MAT_U  - left singular vectors
.   DS_MAT_V  - right singular vectors
.   DS_MAT_W  - workspace matrix
-   DS_MAT_Ex - extra matrices (x=0,..,9)

    All matrices can have space to hold ld x ld elements, except for
    DS_MAT_T that has space for 3 x ld elements (ld = leading dimension)
    and DS_MAT_D that has space for just ld elements.

    In DSPEP problems, matrices A, B, W can have space for d*ld x d*ld,
    where d is the polynomial degree, and X can have ld x d*ld.
    Also DSNEP has exceptions. Check the manual page of each DS type
    for details.

    Level: advanced

.seealso: DSAllocate(), DSGetArray(), DSGetArrayReal(), DSVectors()
E*/
typedef enum { DS_MAT_A,
               DS_MAT_B,
               DS_MAT_C,
               DS_MAT_T,
               DS_MAT_D,
               DS_MAT_Q,
               DS_MAT_Z,
               DS_MAT_X,
               DS_MAT_Y,
               DS_MAT_U,
               DS_MAT_V,
               DS_MAT_W,
               DS_MAT_E0,
               DS_MAT_E1,
               DS_MAT_E2,
               DS_MAT_E3,
               DS_MAT_E4,
               DS_MAT_E5,
               DS_MAT_E6,
               DS_MAT_E7,
               DS_MAT_E8,
               DS_MAT_E9,
               DS_NUM_MAT } DSMatType;

/* Convenience for indexing extra matrices */
SLEPC_EXTERN DSMatType DSMatExtra[];
#define DS_NUM_EXTRA  10

/*E
    DSParallelType - Indicates the parallel mode that the direct solver will use

    Level: advanced

.seealso: DSSetParallel()
E*/
typedef enum { DS_PARALLEL_REDUNDANT,
               DS_PARALLEL_SYNCHRONIZED,
               DS_PARALLEL_DISTRIBUTED } DSParallelType;
SLEPC_EXTERN const char *DSParallelTypes[];

SLEPC_EXTERN PetscErrorCode DSCreate(MPI_Comm,DS*);
SLEPC_EXTERN PetscErrorCode DSSetType(DS,DSType);
SLEPC_EXTERN PetscErrorCode DSGetType(DS,DSType*);
SLEPC_EXTERN PetscErrorCode DSSetOptionsPrefix(DS,const char *);
SLEPC_EXTERN PetscErrorCode DSAppendOptionsPrefix(DS,const char *);
SLEPC_EXTERN PetscErrorCode DSGetOptionsPrefix(DS,const char *[]);
SLEPC_EXTERN PetscErrorCode DSSetFromOptions(DS);
SLEPC_EXTERN PetscErrorCode DSView(DS,PetscViewer);
SLEPC_EXTERN PetscErrorCode DSViewFromOptions(DS,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode DSViewMat(DS,PetscViewer,DSMatType);
SLEPC_EXTERN PetscErrorCode DSDestroy(DS*);
SLEPC_EXTERN PetscErrorCode DSReset(DS);
SLEPC_EXTERN PetscErrorCode DSDuplicate(DS,DS*);

SLEPC_EXTERN PetscErrorCode DSAllocate(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSReallocate(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGetLeadingDimension(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSSetState(DS,DSStateType);
SLEPC_EXTERN PetscErrorCode DSGetState(DS,DSStateType*);
SLEPC_EXTERN PetscErrorCode DSSetDimensions(DS,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGetDimensions(DS,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSSetBlockSize(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGetBlockSize(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSGetTruncateSize(DS,PetscInt,PetscInt,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSTruncate(DS,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode DSSetIdentity(DS,DSMatType);
SLEPC_EXTERN PetscErrorCode DSSetMethod(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGetMethod(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSSetParallel(DS,DSParallelType);
SLEPC_EXTERN PetscErrorCode DSGetParallel(DS,DSParallelType*);
SLEPC_EXTERN PetscErrorCode DSSetCompact(DS,PetscBool);
SLEPC_EXTERN PetscErrorCode DSGetCompact(DS,PetscBool*);
SLEPC_EXTERN PetscErrorCode DSSetExtraRow(DS,PetscBool);
SLEPC_EXTERN PetscErrorCode DSGetExtraRow(DS,PetscBool*);
SLEPC_EXTERN PetscErrorCode DSSetRefined(DS,PetscBool);
SLEPC_EXTERN PetscErrorCode DSGetRefined(DS,PetscBool*);
SLEPC_EXTERN PetscErrorCode DSGetMat(DS,DSMatType,Mat*);
SLEPC_EXTERN PetscErrorCode DSRestoreMat(DS,DSMatType,Mat*);
SLEPC_EXTERN PetscErrorCode DSGetMatAndColumn(DS,DSMatType,PetscInt,Mat*,Vec*);
SLEPC_EXTERN PetscErrorCode DSRestoreMatAndColumn(DS,DSMatType,PetscInt,Mat*,Vec*);
SLEPC_EXTERN PetscErrorCode DSGetArray(DS,DSMatType,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode DSRestoreArray(DS,DSMatType,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode DSGetArrayReal(DS,DSMatType,PetscReal*[]);
SLEPC_EXTERN PetscErrorCode DSRestoreArrayReal(DS,DSMatType,PetscReal*[]);
SLEPC_EXTERN PetscErrorCode DSVectors(DS,DSMatType,PetscInt*,PetscReal*);
SLEPC_EXTERN PetscErrorCode DSSolve(DS,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode DSSort(DS,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSSortWithPermutation(DS,PetscInt*,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode DSSynchronize(DS,PetscScalar*,PetscScalar*);
PETSC_DEPRECATED_FUNCTION(3, 18, 0, "DSGetMat()+MatDenseGetSubMatrix()+MatCopy()", ) static inline PetscErrorCode DSCopyMat(DS ds,DSMatType m,PetscInt mr,PetscInt mc,Mat A,PetscInt Ar,PetscInt Ac,PetscInt rows,PetscInt cols,PetscBool out)
{
  Mat M,M0,A0;

  PetscFunctionBegin;
  PetscCall(DSGetMat(ds,m,&M));
  PetscCall(MatDenseGetSubMatrix(M,mr,mr+rows,mc,mc+cols,&M0));
  PetscCall(MatDenseGetSubMatrix(A,Ar,Ar+rows,Ac,Ac+cols,&A0));
  if (out) PetscCall(MatCopy(M0,A0,SAME_NONZERO_PATTERN));
  else PetscCall(MatCopy(A0,M0,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseRestoreSubMatrix(M,&M0));
  PetscCall(MatDenseRestoreSubMatrix(A,&A0));
  PetscCall(DSRestoreMat(ds,m,&M));
  PetscFunctionReturn(PETSC_SUCCESS);
}
SLEPC_EXTERN PetscErrorCode DSMatGetSize(DS,DSMatType,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSMatIsHermitian(DS,DSMatType,PetscBool*);
SLEPC_EXTERN PetscErrorCode DSSetSlepcSC(DS,SlepcSC);
SLEPC_EXTERN PetscErrorCode DSGetSlepcSC(DS,SlepcSC*);
SLEPC_EXTERN PetscErrorCode DSUpdateExtraRow(DS);
SLEPC_EXTERN PetscErrorCode DSCond(DS,PetscReal*);
SLEPC_EXTERN PetscErrorCode DSTranslateHarmonic(DS,PetscScalar,PetscReal,PetscBool,PetscScalar*,PetscReal*);
SLEPC_EXTERN PetscErrorCode DSTranslateRKS(DS,PetscScalar);
SLEPC_EXTERN PetscErrorCode DSOrthogonalize(DS,DSMatType,PetscInt,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSPseudoOrthogonalize(DS,DSMatType,PetscInt,PetscReal*,PetscInt*,PetscReal*);

/* --------- options specific to particular solvers -------- */

SLEPC_EXTERN PetscErrorCode DSSVDSetDimensions(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSSVDGetDimensions(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSGSVDSetDimensions(DS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGSVDGetDimensions(DS,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSHSVDSetDimensions(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSHSVDGetDimensions(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSHSVDSetReorthogonalize(DS,PetscBool);
SLEPC_EXTERN PetscErrorCode DSHSVDGetReorthogonalize(DS,PetscBool*);

SLEPC_EXTERN PetscErrorCode DSPEPSetDegree(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSPEPGetDegree(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSPEPSetCoefficients(DS,PetscReal*);
SLEPC_EXTERN PetscErrorCode DSPEPGetCoefficients(DS,PetscReal**);

SLEPC_EXTERN PetscErrorCode DSNEPSetFN(DS,PetscInt,FN*);
SLEPC_EXTERN PetscErrorCode DSNEPGetFN(DS,PetscInt,FN*);
SLEPC_EXTERN PetscErrorCode DSNEPGetNumFN(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSNEPSetMinimality(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSNEPGetMinimality(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSNEPSetRefine(DS,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode DSNEPGetRefine(DS,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSNEPSetIntegrationPoints(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSNEPGetIntegrationPoints(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSNEPSetSamplingSize(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSNEPGetSamplingSize(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSNEPSetRG(DS,RG);
SLEPC_EXTERN PetscErrorCode DSNEPGetRG(DS,RG*);

/*S
  DSNEPMatrixFunctionFn - A prototype of a DSNEP compute matrix function that would be passed to DSNEPSetComputeMatrixFunction()

  Calling Sequence:
+   ds     - the direct solver object
.   lambda - point where T(lambda) or T'(lambda) must be evaluated
.   deriv  - if true compute T'(lambda), otherwise compute T(lambda)
.   mat    - the DS matrix where the result must be stored
-   ctx    - [optional] user-defined context for private data for the
             matrix evaluation routine (may be `NULL`)

  Level: developer

.seealso: DSNEPSetComputeMatrixFunction()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(DSNEPMatrixFunctionFn)(DS ds,PetscScalar lambda,PetscBool deriv,DSMatType mat,void *ctx);

SLEPC_EXTERN PetscErrorCode DSNEPSetComputeMatrixFunction(DS,DSNEPMatrixFunctionFn*,void*);
SLEPC_EXTERN PetscErrorCode DSNEPGetComputeMatrixFunction(DS,DSNEPMatrixFunctionFn**,void**);

SLEPC_EXTERN PetscFunctionList DSList;
SLEPC_EXTERN PetscErrorCode DSRegister(const char[],PetscErrorCode(*)(DS));
