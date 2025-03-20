/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcfn.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = FN */

SLEPC_EXTERN PetscBool FNRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode FNRegisterAll(void);
SLEPC_EXTERN PetscLogEvent FN_Evaluate;

typedef struct _FNOps *FNOps;

struct _FNOps {
  PetscErrorCode (*evaluatefunction)(FN,PetscScalar,PetscScalar*);
  PetscErrorCode (*evaluatederivative)(FN,PetscScalar,PetscScalar*);
  PetscErrorCode (*evaluatefunctionmat[FN_MAX_SOLVE])(FN,Mat,Mat);
  PetscErrorCode (*evaluatefunctionmatcuda[FN_MAX_SOLVE])(FN,Mat,Mat);
  PetscErrorCode (*evaluatefunctionmatvec[FN_MAX_SOLVE])(FN,Mat,Vec);
  PetscErrorCode (*evaluatefunctionmatveccuda[FN_MAX_SOLVE])(FN,Mat,Vec);
  PetscErrorCode (*setfromoptions)(FN,PetscOptionItems);
  PetscErrorCode (*view)(FN,PetscViewer);
  PetscErrorCode (*duplicate)(FN,MPI_Comm,FN*);
  PetscErrorCode (*destroy)(FN);
};

#define FN_MAX_W 6

struct _p_FN {
  PETSCHEADER(struct _FNOps);
  /*------------------------- User parameters --------------------------*/
  PetscScalar    alpha;          /* inner scaling (argument) */
  PetscScalar    beta;           /* outer scaling (result) */
  PetscInt       method;         /* the method to compute matrix functions */
  FNParallelType pmode;          /* parallel mode (redundant or synchronized) */

  /*---------------------- Cached data and workspace -------------------*/
  Mat            W[FN_MAX_W];    /* workspace matrices */
  PetscInt       nw;             /* number of allocated W matrices */
  PetscInt       cw;             /* current W matrix */
  void           *data;
};

/*
  FN_AllocateWorkMat - Allocate a work Mat of the same dimension of A and copy
  its contents. The work matrix is returned in M and should be freed with
  FN_FreeWorkMat().
*/
static inline PetscErrorCode FN_AllocateWorkMat(FN fn,Mat A,Mat *M)
{
  PetscInt       n,na;
  PetscBool      create=PETSC_FALSE;

  PetscFunctionBegin;
  *M = NULL;
  PetscCheck(fn->cw<FN_MAX_W,PETSC_COMM_SELF,PETSC_ERR_SUP,"Too many requested work matrices %" PetscInt_FMT,fn->cw);
  if (fn->nw<=fn->cw) {
    create=PETSC_TRUE;
    fn->nw++;
  } else {
    PetscCall(MatGetSize(fn->W[fn->cw],&n,NULL));
    PetscCall(MatGetSize(A,&na,NULL));
    if (n!=na) {
      PetscCall(MatDestroy(&fn->W[fn->cw]));
      create=PETSC_TRUE;
    }
  }
  if (create) PetscCall(MatDuplicate(A,MAT_COPY_VALUES,&fn->W[fn->cw]));
  else PetscCall(MatCopy(A,fn->W[fn->cw],SAME_NONZERO_PATTERN));
  *M = fn->W[fn->cw];
  fn->cw++;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  FN_FreeWorkMat - Release a work matrix created with FN_AllocateWorkMat().
*/
static inline PetscErrorCode FN_FreeWorkMat(FN fn,Mat *M)
{
  PetscFunctionBegin;
  PetscCheck(fn->cw,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"There are no work matrices");
  fn->cw--;
  PetscCheck(fn->W[fn->cw]==*M,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Work matrices must be freed in the reverse order of their creation");
  *M = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_INTERN PetscErrorCode FNSqrtmSchur(FN,PetscBLASInt,PetscScalar*,PetscBLASInt,PetscBool);
SLEPC_INTERN PetscErrorCode FNSqrtmDenmanBeavers(FN,PetscBLASInt,PetscScalar*,PetscBLASInt,PetscBool);
SLEPC_INTERN PetscErrorCode FNSqrtmNewtonSchulz(FN,PetscBLASInt,PetscScalar*,PetscBLASInt,PetscBool);
SLEPC_INTERN PetscErrorCode FNSqrtmSadeghi(FN,PetscBLASInt,PetscScalar*,PetscBLASInt);
SLEPC_INTERN PetscErrorCode SlepcNormAm(PetscBLASInt,PetscScalar*,PetscInt,PetscScalar*,PetscRandom,PetscReal*);
SLEPC_INTERN PetscErrorCode FNEvaluateFunctionMat_Private(FN,Mat,Mat,PetscBool);
SLEPC_INTERN PetscErrorCode FNEvaluateFunctionMatVec_Private(FN,Mat,Vec,PetscBool);
SLEPC_INTERN PetscErrorCode FNEvaluateFunctionMat_Exp_Higham(FN,Mat,Mat); /* used in FNPHI */
#if defined(PETSC_HAVE_CUDA)
SLEPC_INTERN PetscErrorCode FNSqrtmDenmanBeavers_CUDAm(FN,PetscBLASInt,PetscScalar*,PetscBLASInt,PetscBool);
SLEPC_INTERN PetscErrorCode FNSqrtmNewtonSchulz_CUDA(FN,PetscBLASInt,PetscScalar*,PetscBLASInt,PetscBool);
SLEPC_INTERN PetscErrorCode FNSqrtmSadeghi_CUDAm(FN,PetscBLASInt,PetscScalar*,PetscBLASInt);
#endif
