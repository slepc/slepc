/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the KSVD SVD solver
*/

#include <slepc/private/svdimpl.h>    /*I "slepcsvd.h" I*/
#include <slepc/private/slepcscalapack.h>
#include <ksvd.h>

typedef struct {
  Mat                As;        /* converted matrix */
  SVDKSVDEigenMethod eigen;
  SVDKSVDPolarMethod polar;
} SVD_KSVD;

static PetscErrorCode SVDSetUp_KSVD(SVD svd)
{
  SVD_KSVD  *ctx = (SVD_KSVD*)svd->data;
  PetscInt  M,N;

  PetscFunctionBegin;
  SVDCheckStandard(svd);
  SVDCheckDefinite(svd);
  if (svd->nsv==0) svd->nsv = 1;
  PetscCall(MatGetSize(svd->A,&M,&N));
  PetscCheck(M==N,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"The interface to KSVD does not support rectangular matrices");
  svd->ncv = N;
  if (svd->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(svd,"Warning: parameter mpd ignored\n"));
  if (svd->max_it==PETSC_DETERMINE) svd->max_it = 1;
  svd->leftbasis = PETSC_TRUE;
  SVDCheckIgnored(svd,SVD_FEATURE_STOPPING);
  PetscCall(SVDAllocateSolution(svd,0));

  /* default methods */
  if (!ctx->eigen) ctx->eigen = SVD_KSVD_EIGEN_MRRR;
  if (!ctx->polar) ctx->polar = SVD_KSVD_POLAR_QDWH;

  /* convert matrix */
  PetscCall(MatDestroy(&ctx->As));
  PetscCall(MatConvert(svd->OP,MATSCALAPACK,MAT_INITIAL_MATRIX,&ctx->As));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_KSVD(SVD svd)
{
  SVD_KSVD       *ctx = (SVD_KSVD*)svd->data;
  Mat            A = ctx->As,Z,Q,U,V;
  Mat_ScaLAPACK  *a = (Mat_ScaLAPACK*)A->data,*q,*z;
  PetscScalar    *work,minlwork;
  PetscBLASInt   info,lwork=-1,*iwork,liwork=-1,minliwork,one=1;
  PetscInt       M,N,m,n,mn;
  const char     *eigen;
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
  PetscCall(MatCreate(PetscObjectComm((PetscObject)A),&Q));
  PetscCall(MatSetSizes(Q,mn,n,PETSC_DECIDE,PETSC_DECIDE));
  PetscCall(MatSetType(Q,MATSCALAPACK));
  PetscCall(MatAssemblyBegin(Q,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Q,MAT_FINAL_ASSEMBLY));
  q = (Mat_ScaLAPACK*)Q->data;

  /* configure solver */
  setPolarAlgorithm((int)ctx->polar);
  switch (ctx->eigen) {
    case SVD_KSVD_EIGEN_MRRR: eigen = "r"; break;
    case SVD_KSVD_EIGEN_DC:   eigen = "d"; break;
    case SVD_KSVD_EIGEN_ELPA: eigen = "e"; break;
  }

  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  /* allocate workspace */
  PetscStackCallExternalVoid("pdgeqsvd",pdgeqsvd("V","V",eigen,a->M,a->N,a->loc,one,one,a->desc,svd->sigma,z->loc,one,one,z->desc,q->loc,one,one,q->desc,&minlwork,lwork,&minliwork,liwork,&info));
  PetscCheck(!info,PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"Error in KSVD subroutine pdgeqsvd: info=%d",(int)info);
  PetscCall(PetscBLASIntCast((PetscInt)minlwork,&lwork));
  PetscCall(PetscBLASIntCast(minliwork,&liwork));
  PetscCall(PetscMalloc2(lwork,&work,liwork,&iwork));
  /* call computational routine */
  PetscStackCallExternalVoid("pdgeqsvd",pdgeqsvd("V","V",eigen,a->M,a->N,a->loc,one,one,a->desc,svd->sigma,z->loc,one,one,z->desc,q->loc,one,one,q->desc,work,lwork,iwork,liwork,&info));
  PetscCheck(!info,PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"Error in KSVD subroutine pdgeqsvd: info=%d",(int)info);
  PetscCall(PetscFree2(work,iwork));
  PetscCall(PetscFPTrapPop());

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

static PetscErrorCode SVDDestroy_KSVD(SVD svd)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(svd->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDKSVDSetEigenMethod_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDKSVDGetEigenMethod_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDKSVDSetPolarMethod_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDKSVDGetPolarMethod_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDReset_KSVD(SVD svd)
{
  SVD_KSVD  *ctx = (SVD_KSVD*)svd->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->As));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDView_KSVD(SVD svd,PetscViewer viewer)
{
  PetscBool      isascii;
  SVD_KSVD       *ctx = (SVD_KSVD*)svd->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  eigensolver method: %s\n",SVDKSVDEigenMethods[ctx->eigen]));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  polar decomposition method: %s\n",SVDKSVDPolarMethods[ctx->polar]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSetFromOptions_KSVD(SVD svd,PetscOptionItems *PetscOptionsObject)
{
  SVD_KSVD           *ctx = (SVD_KSVD*)svd->data;
  SVDKSVDEigenMethod eigen;
  SVDKSVDPolarMethod polar;
  PetscBool          flg;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"SVD KSVD Options");

    PetscCall(PetscOptionsEnum("-svd_ksvd_eigen_method","Method for solving the internal eigenvalue problem","SVDKSVDSetEigenMethod",SVDKSVDEigenMethods,(PetscEnum)ctx->eigen,(PetscEnum*)&eigen,&flg));
    if (flg) PetscCall(SVDKSVDSetEigenMethod(svd,eigen));
    PetscCall(PetscOptionsEnum("-svd_ksvd_polar_method","Method for computing the polar decomposition","SVDKSVDSetPolarMethod",SVDKSVDPolarMethods,(PetscEnum)ctx->polar,(PetscEnum*)&polar,&flg));
    if (flg) PetscCall(SVDKSVDSetPolarMethod(svd,polar));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDKSVDSetEigenMethod_KSVD(SVD svd,SVDKSVDEigenMethod eigen)
{
  SVD_KSVD *ctx = (SVD_KSVD*)svd->data;

  PetscFunctionBegin;
  ctx->eigen = eigen;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDKSVDSetEigenMethod - Sets the method to solve the internal eigenproblem within the KSVD solver.

   Logically Collective

   Input Parameters:
+  svd - the singular value solver context
-  eigen - method that will be used by KSVD for the eigenproblem

   Options Database Key:
.  -svd_ksvd_eigen_method - Sets the method for the KSVD eigensolver

   If not set, the method defaults to SVD_KSVD_EIGEN_MRRR.

   Level: advanced

.seealso: SVDKSVDGetEigenMethod(), SVDKSVDEigenMethod
@*/
PetscErrorCode SVDKSVDSetEigenMethod(SVD svd,SVDKSVDEigenMethod eigen)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,eigen,2);
  PetscTryMethod(svd,"SVDKSVDSetEigenMethod_C",(SVD,SVDKSVDEigenMethod),(svd,eigen));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDKSVDGetEigenMethod_KSVD(SVD svd,SVDKSVDEigenMethod *eigen)
{
  SVD_KSVD *ctx = (SVD_KSVD*)svd->data;

  PetscFunctionBegin;
  *eigen = ctx->eigen;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDKSVDGetEigenMethod - Gets the method for the KSVD eigensolver.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  eigen - method that will be used by KSVD for the eigenproblem

   Level: advanced

.seealso: SVDKSVDSetEigenMethod(), SVDKSVDEigenMethod
@*/
PetscErrorCode SVDKSVDGetEigenMethod(SVD svd,SVDKSVDEigenMethod *eigen)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(eigen,2);
  PetscUseMethod(svd,"SVDKSVDGetEigenMethod_C",(SVD,SVDKSVDEigenMethod*),(svd,eigen));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDKSVDSetPolarMethod_KSVD(SVD svd,SVDKSVDPolarMethod polar)
{
  SVD_KSVD *ctx = (SVD_KSVD*)svd->data;

  PetscFunctionBegin;
  ctx->polar = polar;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDKSVDSetPolarMethod - Sets the method to compute the polar decomposition within the KSVD solver.

   Logically Collective

   Input Parameters:
+  svd - the singular value solver context
-  polar - method that will be used by KSVD for the polar decomposition

   Options Database Key:
.  -svd_ksvd_polar_method - Sets the method for the KSVD polar decomposition

   If not set, the method defaults to SVD_KSVD_POLAR_QDWH.

   Level: advanced

.seealso: SVDKSVDGetPolarMethod(), SVDKSVDPolarMethod
@*/
PetscErrorCode SVDKSVDSetPolarMethod(SVD svd,SVDKSVDPolarMethod polar)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,polar,2);
  PetscTryMethod(svd,"SVDKSVDSetPolarMethod_C",(SVD,SVDKSVDPolarMethod),(svd,polar));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDKSVDGetPolarMethod_KSVD(SVD svd,SVDKSVDPolarMethod *polar)
{
  SVD_KSVD *ctx = (SVD_KSVD*)svd->data;

  PetscFunctionBegin;
  *polar = ctx->polar;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDKSVDGetPolarMethod - Gets the method for the KSVD polar decomposition.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  polar - method that will be used by KSVD for the polar decomposition

   Level: advanced

.seealso: SVDKSVDSetPolarMethod(), SVDKSVDPolarMethod
@*/
PetscErrorCode SVDKSVDGetPolarMethod(SVD svd,SVDKSVDPolarMethod *polar)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(polar,2);
  PetscUseMethod(svd,"SVDKSVDGetPolarMethod_C",(SVD,SVDKSVDPolarMethod*),(svd,polar));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_KSVD(SVD svd)
{
  SVD_KSVD  *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  svd->data = (void*)ctx;

  svd->ops->solve          = SVDSolve_KSVD;
  svd->ops->setup          = SVDSetUp_KSVD;
  svd->ops->setfromoptions = SVDSetFromOptions_KSVD;
  svd->ops->destroy        = SVDDestroy_KSVD;
  svd->ops->reset          = SVDReset_KSVD;
  svd->ops->view           = SVDView_KSVD;

  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDKSVDSetEigenMethod_C",SVDKSVDSetEigenMethod_KSVD));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDKSVDGetEigenMethod_C",SVDKSVDGetEigenMethod_KSVD));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDKSVDSetPolarMethod_C",SVDKSVDSetPolarMethod_KSVD));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDKSVDGetPolarMethod_C",SVDKSVDGetPolarMethod_KSVD));
  PetscFunctionReturn(PETSC_SUCCESS);
}
