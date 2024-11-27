/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the Elemental SVD solver
*/

#include <slepc/private/svdimpl.h>    /*I "slepcsvd.h" I*/
#include <petsc/private/petscelemental.h>

typedef struct {
  Mat Ae;        /* converted matrix */
} SVD_Elemental;

static PetscErrorCode SVDSetUp_Elemental(SVD svd)
{
  SVD_Elemental  *ctx = (SVD_Elemental*)svd->data;
  PetscInt       M,N;

  PetscFunctionBegin;
  SVDCheckStandard(svd);
  SVDCheckDefinite(svd);
  if (svd->nsv==0) svd->nsv = 1;
  PetscCall(MatGetSize(svd->A,&M,&N));
  PetscCheck(M==N,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Not implemented for rectangular matrices");
  svd->ncv = N;
  if (svd->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(svd,"Warning: parameter mpd ignored\n"));
  if (svd->max_it==PETSC_DETERMINE) svd->max_it = 1;
  svd->leftbasis = PETSC_TRUE;
  SVDCheckIgnored(svd,SVD_FEATURE_STOPPING);
  PetscCall(SVDAllocateSolution(svd,0));

  /* convert matrix */
  PetscCall(MatDestroy(&ctx->Ae));
  PetscCall(MatConvert(svd->OP,MATELEMENTAL,MAT_INITIAL_MATRIX,&ctx->Ae));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_Elemental(SVD svd)
{
  SVD_Elemental  *ctx = (SVD_Elemental*)svd->data;
  Mat            A = ctx->Ae,Z,Q,U,V;
  Mat_Elemental  *a = (Mat_Elemental*)A->data,*q,*z;
  PetscInt       i,rrank,ridx,erow;

  PetscFunctionBegin;
  El::DistMatrix<PetscReal,El::STAR,El::VC> sigma(*a->grid);
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&Z));
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&Q));
  z = (Mat_Elemental*)Z->data;
  q = (Mat_Elemental*)Q->data;

  El::SVD(*a->emat,*z->emat,sigma,*q->emat);

  for (i=0;i<svd->ncv;i++) {
    P2RO(A,1,i,&rrank,&ridx);
    RO2E(A,1,rrank,ridx,&erow);
    svd->sigma[i] = sigma.Get(erow,0);
  }
  PetscCall(BVGetMat(svd->U,&U));
  PetscCall(MatConvert(Z,MATDENSE,MAT_REUSE_MATRIX,&U));
  PetscCall(BVRestoreMat(svd->U,&U));
  PetscCall(MatDestroy(&Z));
  PetscCall(BVGetMat(svd->V,&V));
  PetscCall(MatConvert(Q,MATDENSE,MAT_REUSE_MATRIX,&V));
  PetscCall(BVRestoreMat(svd->V,&V));
  PetscCall(MatDestroy(&Q));

  svd->nconv  = svd->ncv;
  svd->its    = 1;
  svd->reason = SVD_CONVERGED_TOL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDDestroy_Elemental(SVD svd)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(svd->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDReset_Elemental(SVD svd)
{
  SVD_Elemental  *ctx = (SVD_Elemental*)svd->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->Ae));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_Elemental(SVD svd)
{
  SVD_Elemental  *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  svd->data = (void*)ctx;

  svd->ops->solve          = SVDSolve_Elemental;
  svd->ops->setup          = SVDSetUp_Elemental;
  svd->ops->destroy        = SVDDestroy_Elemental;
  svd->ops->reset          = SVDReset_Elemental;
  PetscFunctionReturn(PETSC_SUCCESS);
}
