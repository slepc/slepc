/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the ScaLAPACK SVD solver
*/

#include <slepc/private/svdimpl.h>    /*I "slepcsvd.h" I*/
#include <slepc/private/slepcscalapack.h>

typedef struct {
  Mat As;        /* converted matrix */
} SVD_ScaLAPACK;

static PetscErrorCode SVDSetUp_ScaLAPACK(SVD svd)
{
  SVD_ScaLAPACK  *ctx = (SVD_ScaLAPACK*)svd->data;
  PetscInt       M,N;

  PetscFunctionBegin;
  SVDCheckStandard(svd);
  SVDCheckDefinite(svd);
  if (svd->nsv==0) svd->nsv = 1;
  PetscCall(MatGetSize(svd->A,&M,&N));
  svd->ncv = N;
  if (svd->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(svd,"Warning: parameter mpd ignored\n"));
  if (svd->max_it==PETSC_DETERMINE) svd->max_it = 1;
  svd->leftbasis = PETSC_TRUE;
  SVDCheckIgnored(svd,SVD_FEATURE_STOPPING);
  PetscCall(SVDAllocateSolution(svd,0));

  /* convert matrix */
  PetscCall(MatDestroy(&ctx->As));
  PetscCall(MatConvert(svd->OP,MATSCALAPACK,MAT_INITIAL_MATRIX,&ctx->As));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_ScaLAPACK(SVD svd)
{
  SVD_ScaLAPACK  *ctx = (SVD_ScaLAPACK*)svd->data;
  Mat            A = ctx->As,Z,Q,QT,U,V;
  Mat_ScaLAPACK  *a = (Mat_ScaLAPACK*)A->data,*q,*z;
  PetscScalar    *work,minlwork;
  PetscBLASInt   info,lwork=-1,one=1;
  PetscInt       M,N,m,n,mn;
#if defined(PETSC_USE_COMPLEX)
  PetscBLASInt   lrwork;
  PetscReal      *rwork,dummy;
#endif

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));
  mn = (M>=N)? n: m;
  PetscCall(MatCreate(PetscObjectComm((PetscObject)A),&Z));
  PetscCall(MatSetSizes(Z,m,mn,PETSC_DECIDE,PETSC_DECIDE));
  PetscCall(MatSetType(Z,MATSCALAPACK));
  PetscCall(MatAssemblyBegin(Z,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Z,MAT_FINAL_ASSEMBLY));
  z = (Mat_ScaLAPACK*)Z->data;
  PetscCall(MatCreate(PetscObjectComm((PetscObject)A),&QT));
  PetscCall(MatSetSizes(QT,mn,n,PETSC_DECIDE,PETSC_DECIDE));
  PetscCall(MatSetType(QT,MATSCALAPACK));
  PetscCall(MatAssemblyBegin(QT,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(QT,MAT_FINAL_ASSEMBLY));
  q = (Mat_ScaLAPACK*)QT->data;

  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined(PETSC_USE_COMPLEX)
  /* allocate workspace */
  PetscCallBLAS("SCALAPACKgesvd",SCALAPACKgesvd_("V","V",&a->M,&a->N,a->loc,&one,&one,a->desc,svd->sigma,z->loc,&one,&one,z->desc,q->loc,&one,&one,q->desc,&minlwork,&lwork,&info));
  PetscCheckScaLapackInfo("gesvd",info);
  PetscCall(PetscBLASIntCast((PetscInt)minlwork,&lwork));
  PetscCall(PetscMalloc1(lwork,&work));
  /* call computational routine */
  PetscCallBLAS("SCALAPACKgesvd",SCALAPACKgesvd_("V","V",&a->M,&a->N,a->loc,&one,&one,a->desc,svd->sigma,z->loc,&one,&one,z->desc,q->loc,&one,&one,q->desc,work,&lwork,&info));
  PetscCheckScaLapackInfo("gesvd",info);
  PetscCall(PetscFree(work));
#else
  /* allocate workspace */
  PetscCallBLAS("SCALAPACKgesvd",SCALAPACKgesvd_("V","V",&a->M,&a->N,a->loc,&one,&one,a->desc,svd->sigma,z->loc,&one,&one,z->desc,q->loc,&one,&one,q->desc,&minlwork,&lwork,&dummy,&info));
  PetscCheckScaLapackInfo("gesvd",info);
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(minlwork),&lwork));
  lrwork = 1+4*PetscMax(a->M,a->N);
  PetscCall(PetscMalloc2(lwork,&work,lrwork,&rwork));
  /* call computational routine */
  PetscCallBLAS("SCALAPACKgesvd",SCALAPACKgesvd_("V","V",&a->M,&a->N,a->loc,&one,&one,a->desc,svd->sigma,z->loc,&one,&one,z->desc,q->loc,&one,&one,q->desc,work,&lwork,rwork,&info));
  PetscCheckScaLapackInfo("gesvd",info);
  PetscCall(PetscFree2(work,rwork));
#endif
  PetscCall(PetscFPTrapPop());

  PetscCall(MatHermitianTranspose(QT,MAT_INITIAL_MATRIX,&Q));
  PetscCall(MatDestroy(&QT));
  PetscCall(BVGetMat(svd->U,&U));
  PetscCall(BVGetMat(svd->V,&V));
  if (M>=N) {
    PetscCall(MatConvert(Z,MATDENSE,MAT_REUSE_MATRIX,&U));
    PetscCall(MatConvert(Q,MATDENSE,MAT_REUSE_MATRIX,&V));
  } else {
    PetscCall(MatConvert(Q,MATDENSE,MAT_REUSE_MATRIX,&U));
    PetscCall(MatConvert(Z,MATDENSE,MAT_REUSE_MATRIX,&V));
  }
  PetscCall(BVRestoreMat(svd->U,&U));
  PetscCall(BVRestoreMat(svd->V,&V));
  PetscCall(MatDestroy(&Z));
  PetscCall(MatDestroy(&Q));

  svd->nconv  = svd->ncv;
  svd->its    = 1;
  svd->reason = SVD_CONVERGED_TOL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDDestroy_ScaLAPACK(SVD svd)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(svd->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDReset_ScaLAPACK(SVD svd)
{
  SVD_ScaLAPACK  *ctx = (SVD_ScaLAPACK*)svd->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->As));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_ScaLAPACK(SVD svd)
{
  SVD_ScaLAPACK  *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  svd->data = (void*)ctx;

  svd->ops->solve          = SVDSolve_ScaLAPACK;
  svd->ops->setup          = SVDSetUp_ScaLAPACK;
  svd->ops->destroy        = SVDDestroy_ScaLAPACK;
  svd->ops->reset          = SVDReset_ScaLAPACK;
  PetscFunctionReturn(PETSC_SUCCESS);
}
