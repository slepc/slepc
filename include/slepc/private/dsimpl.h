/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcds.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = DS */

SLEPC_EXTERN PetscBool DSRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode DSRegisterAll(void);
SLEPC_EXTERN PetscLogEvent DS_Solve,DS_Vectors,DS_Synchronize,DS_Other;
SLEPC_INTERN const char *DSMatName[];

typedef struct _DSOps *DSOps;

struct _DSOps {
  PetscErrorCode (*allocate)(DS,PetscInt);
  PetscErrorCode (*setfromoptions)(DS,PetscOptionItems*);
  PetscErrorCode (*view)(DS,PetscViewer);
  PetscErrorCode (*vectors)(DS,DSMatType,PetscInt*,PetscReal*);
  PetscErrorCode (*solve[DS_MAX_SOLVE])(DS,PetscScalar*,PetscScalar*);
  PetscErrorCode (*sort)(DS,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*);
  PetscErrorCode (*sortperm)(DS,PetscInt*,PetscScalar*,PetscScalar*);
  PetscErrorCode (*gettruncatesize)(DS,PetscInt,PetscInt,PetscInt*);
  PetscErrorCode (*truncate)(DS,PetscInt,PetscBool);
  PetscErrorCode (*update)(DS);
  PetscErrorCode (*cond)(DS,PetscReal*);
  PetscErrorCode (*transharm)(DS,PetscScalar,PetscReal,PetscBool,PetscScalar*,PetscReal*);
  PetscErrorCode (*transrks)(DS,PetscScalar);
  PetscErrorCode (*destroy)(DS);
  PetscErrorCode (*matgetsize)(DS,DSMatType,PetscInt*,PetscInt*);
  PetscErrorCode (*hermitian)(DS,DSMatType,PetscBool*);
  PetscErrorCode (*synchronize)(DS,PetscScalar*,PetscScalar*);
  PetscErrorCode (*setcompact)(DS,PetscBool);
  PetscErrorCode (*reallocate)(DS,PetscInt);
};

struct _p_DS {
  PETSCHEADER(struct _DSOps);
  /*------------------------- User parameters --------------------------*/
  DSStateType    state;              /* the current state */
  PetscInt       method;             /* identifies the variant to be used */
  PetscBool      compact;            /* whether the matrices are stored in compact form */
  PetscBool      refined;            /* get refined vectors instead of regular vectors */
  PetscBool      extrarow;           /* assume the matrix dimension is (n+1) x n */
  PetscInt       ld;                 /* leading dimension */
  PetscInt       l;                  /* number of locked (inactive) leading columns */
  PetscInt       n;                  /* current dimension */
  PetscInt       k;                  /* intermediate dimension (e.g. position of arrow) */
  PetscInt       t;                  /* length of decomposition when it was truncated */
  PetscInt       bs;                 /* block size */
  SlepcSC        sc;                 /* sorting criterion */
  DSParallelType pmode;              /* parallel mode (redundant, synchronized, distributed) */

  /*----------------- Status variables and working data ----------------*/
  Mat            omat[DS_NUM_MAT];   /* the matrices (PETSc object) */
  PetscInt       *perm;              /* permutation */
  void           *data;              /* placeholder for solver-specific stuff */
  PetscBool      scset;              /* the sc was provided by the user */
  PetscScalar    *work;
  PetscReal      *rwork;
  PetscBLASInt   *iwork;
  PetscInt       lwork,lrwork,liwork;
};

/*
    Macros to test valid DS arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define DSCheckAlloc(h,arg) do {(void)(h);} while (0)
#define DSCheckSolved(h,arg) do {(void)(h);} while (0)
#define DSCheckValidMat(ds,m,arg) do {(void)(ds);} while (0)
#define DSCheckValidMatReal(ds,m,arg) do {(void)(ds);} while (0)

#else

#define DSCheckAlloc(h,arg) \
  do { \
    PetscCheck((h)->ld,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"Must call DSAllocate() first: Parameter #%d",arg); \
  } while (0)

#define DSCheckSolved(h,arg) \
  do { \
    PetscCheck((h)->state>=DS_STATE_CONDENSED,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"Must call DSSolve() first: Parameter #%d",arg); \
  } while (0)

#define DSCheckValidMat(ds,m,arg) \
  do { \
    PetscCheck((m)<DS_NUM_MAT,PetscObjectComm((PetscObject)(ds)),PETSC_ERR_ARG_WRONG,"Invalid matrix: Parameter #%d",arg); \
    PetscCheck((ds)->omat[m],PetscObjectComm((PetscObject)(ds)),PETSC_ERR_ARG_WRONGSTATE,"Requested matrix was not created in this DS: Parameter #%d",arg); \
  } while (0)

#define DSCheckValidMatReal(ds,m,arg) \
  do { \
    PetscCheck((m)==DS_MAT_T || (m)==DS_MAT_D,PetscObjectComm((PetscObject)(ds)),PETSC_ERR_ARG_WRONG,"Invalid matrix, can only be used for T and D: Parameter #%d",arg); \
    PetscCheck((ds)->omat[m],PetscObjectComm((PetscObject)(ds)),PETSC_ERR_ARG_WRONGSTATE,"Requested matrix was not created in this DS: Parameter #%d",arg); \
  } while (0)

#endif

SLEPC_INTERN PetscErrorCode DSAllocateMat_Private(DS,DSMatType);
SLEPC_INTERN PetscErrorCode DSReallocateMat_Private(DS,DSMatType,PetscInt);
SLEPC_INTERN PetscErrorCode DSAllocateWork_Private(DS,PetscInt,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode DSSortEigenvalues_Private(DS,PetscScalar*,PetscScalar*,PetscInt*,PetscBool);
SLEPC_INTERN PetscErrorCode DSSortEigenvaluesReal_Private(DS,PetscReal*,PetscInt*);
SLEPC_INTERN PetscErrorCode DSPermuteColumns_Private(DS,PetscInt,PetscInt,PetscInt,DSMatType,PetscInt*);
SLEPC_INTERN PetscErrorCode DSPermuteColumnsTwo_Private(DS,PetscInt,PetscInt,PetscInt,DSMatType,DSMatType,PetscInt*);
SLEPC_INTERN PetscErrorCode DSPermuteRows_Private(DS,PetscInt,PetscInt,PetscInt,DSMatType,PetscInt*);
SLEPC_INTERN PetscErrorCode DSPermuteBoth_Private(DS,PetscInt,PetscInt,PetscInt,PetscInt,DSMatType,DSMatType,PetscInt*);
SLEPC_INTERN PetscErrorCode DSGetTruncateSize_Default(DS,PetscInt,PetscInt,PetscInt*);

SLEPC_INTERN PetscErrorCode DSGHIEPOrthogEigenv(DS,DSMatType,PetscScalar*,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode DSPseudoOrthog_HR(PetscInt*,PetscScalar*,PetscInt,PetscReal*,PetscScalar*,PetscInt,PetscBLASInt*,PetscBLASInt*,PetscBool*,PetscScalar*);
SLEPC_INTERN PetscErrorCode DSGHIEPComplexEigs(DS,PetscInt,PetscInt,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode DSGHIEPInverseIteration(DS,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode DSIntermediate_GHIEP(DS);
SLEPC_INTERN PetscErrorCode DSSwitchFormat_GHIEP(DS,PetscBool);
SLEPC_INTERN PetscErrorCode DSGHIEPRealBlocks(DS);
SLEPC_INTERN PetscErrorCode DSSolve_GHIEP_HZ(DS,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode DSArrowTridiag(PetscBLASInt,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt);

SLEPC_INTERN PetscErrorCode DSSolve_NHEP_Private(DS,DSMatType,DSMatType,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode DSSort_NHEP_Total(DS,DSMatType,DSMatType,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode DSSortWithPermutation_NHEP_Private(DS,PetscInt*,DSMatType,DSMatType,PetscScalar*,PetscScalar*);

SLEPC_INTERN PetscErrorCode BDC_dibtdc_(const char*,PetscBLASInt,PetscBLASInt,PetscBLASInt*,PetscReal*,PetscBLASInt,PetscBLASInt,PetscReal*,PetscBLASInt*,PetscBLASInt,PetscBLASInt,PetscReal,PetscReal*,PetscReal*,PetscBLASInt,PetscReal*,PetscBLASInt,PetscBLASInt*,PetscBLASInt,PetscBLASInt*,PetscBLASInt);
SLEPC_INTERN PetscErrorCode BDC_dlaed3m_(const char*,const char*,PetscBLASInt,PetscBLASInt,PetscBLASInt,PetscReal*,PetscReal*,PetscBLASInt,PetscReal,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt,PetscBLASInt);
SLEPC_INTERN PetscErrorCode BDC_dmerg2_(const char*,PetscBLASInt,PetscBLASInt,PetscReal*,PetscReal*,PetscBLASInt,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt,PetscReal*,PetscBLASInt,PetscBLASInt,PetscReal*,PetscBLASInt,PetscBLASInt*,PetscReal,PetscBLASInt*,PetscBLASInt);
SLEPC_INTERN PetscErrorCode BDC_dsbtdc_(const char*,const char*,PetscBLASInt,PetscBLASInt,PetscBLASInt*,PetscReal*,PetscBLASInt,PetscBLASInt,PetscReal*,PetscBLASInt,PetscBLASInt,PetscReal,PetscReal,PetscReal,PetscReal*,PetscReal*,PetscBLASInt,PetscReal*,PetscBLASInt,PetscBLASInt*,PetscBLASInt,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt,PetscBLASInt);
SLEPC_INTERN PetscErrorCode BDC_dsrtdf_(PetscBLASInt*,PetscBLASInt,PetscBLASInt,PetscReal*,PetscReal*,PetscBLASInt,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscReal,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
