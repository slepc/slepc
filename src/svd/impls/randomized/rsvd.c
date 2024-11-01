/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "randomized"

   Method: RSVD

   Algorithm:

       Randomized singular value decomposition.

   References:

       [1] N. Halko, P.-G. Martinsson, and J. A. Tropp, "Finding
           structure with randomness: Probabilistic algorithms for
           constructing approximate matrix decompositions", SIAM Rev.,
           53(2):217-288, 2011.
*/

#include <slepc/private/svdimpl.h>                /*I "slepcsvd.h" I*/

static PetscErrorCode SVDSetUp_Randomized(SVD svd)
{
  PetscInt       N;

  PetscFunctionBegin;
  SVDCheckStandard(svd);
  SVDCheckDefinite(svd);
  SVDCheckUnsupported(svd,SVD_FEATURE_THRESHOLD);
  PetscCheck(svd->which==SVD_LARGEST,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"This solver supports only largest singular values");
  PetscCall(MatGetSize(svd->A,NULL,&N));
  PetscCall(SVDSetDimensions_Default(svd));
  PetscCheck(svd->ncv>=svd->nsv,PetscObjectComm((PetscObject)svd),PETSC_ERR_USER_INPUT,"The value of ncv must not be smaller than nsv");
  if (svd->max_it==PETSC_DETERMINE) svd->max_it = PetscMax(N/svd->ncv,100);
  svd->leftbasis = PETSC_TRUE;
  svd->mpd = svd->ncv;
  PetscCall(SVDAllocateSolution(svd,0));
  PetscCall(DSSetType(svd->ds,DSSVD));
  PetscCall(DSAllocate(svd->ds,svd->ncv));
  PetscCall(SVDSetWorkVecs(svd,1,1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDRandomizedResidualNorm(SVD svd,PetscInt i,PetscScalar sigma,PetscReal *res)
{
  PetscReal      norm1,norm2;
  Vec            u,v,wu,wv;

  PetscFunctionBegin;
  *res = 1.0;
  if (svd->conv!=SVD_CONV_MAXIT) {
    wu = svd->swapped? svd->workr[0]: svd->workl[0];
    wv = svd->swapped? svd->workl[0]: svd->workr[0];
    PetscCall(BVGetColumn(svd->V,i,&v));
    PetscCall(BVGetColumn(svd->U,i,&u));
    /* norm1 = ||A*v-sigma*u||_2 */
    PetscCall(MatMult(svd->A,v,wu));
    PetscCall(VecAXPY(wu,-sigma,u));
    PetscCall(VecNorm(wu,NORM_2,&norm1));
    /* norm2 = ||A^T*u-sigma*v||_2 */
    PetscCall(MatMult(svd->AT,u,wv));
    PetscCall(VecAXPY(wv,-sigma,v));
    PetscCall(VecNorm(wv,NORM_2,&norm2));
    PetscCall(BVRestoreColumn(svd->V,i,&v));
    PetscCall(BVRestoreColumn(svd->U,i,&u));
    *res = SlepcAbs(norm1,norm2);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* If A is a virtual Hermitian transpose, then BVMatMult will fail if PRODUCT_AhB is not implemented */
static PetscErrorCode BlockMatMult(BV V,Mat A,BV Y,Mat AT)
{
  PetscFunctionBegin;
  if (!PetscDefined(USE_COMPLEX)) PetscCall(BVMatMult(V,A,Y));
  else {
    PetscBool flg=PETSC_FALSE;
    PetscCall(PetscObjectTypeCompare((PetscObject)A,MATHERMITIANTRANSPOSEVIRTUAL,&flg));
    if (flg) PetscCall(BVMatMultHermitianTranspose(V,AT,Y));
    else PetscCall(BVMatMult(V,A,Y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_Randomized(SVD svd)
{
  PetscScalar    *w;
  PetscReal      res=1.0;
  PetscInt       i,k=0;
  Mat            A,U,V;

  PetscFunctionBegin;
  /* Form random matrix, G. Complete the initial basis with random vectors */
  PetscCall(BVSetActiveColumns(svd->V,svd->nini,svd->ncv));
  PetscCall(BVSetRandomNormal(svd->V));
  PetscCall(PetscCalloc1(svd->ncv,&w));

  /* Subspace Iteration */
  do {
    svd->its++;
    PetscCall(BVSetActiveColumns(svd->V,svd->nconv,svd->ncv));
    PetscCall(BVSetActiveColumns(svd->U,svd->nconv,svd->ncv));
    /* Form AG */
    PetscCall(BlockMatMult(svd->V,svd->A,svd->U,svd->AT));
    /* Orthogonalization Q=qr(AG)*/
    PetscCall(BVOrthogonalize(svd->U,NULL));
    /* Form B^*= AQ */
    PetscCall(BlockMatMult(svd->U,svd->AT,svd->V,svd->A));

    PetscCall(DSSetDimensions(svd->ds,svd->ncv,svd->nconv,svd->ncv));
    PetscCall(DSSVDSetDimensions(svd->ds,svd->ncv));
    PetscCall(DSGetMat(svd->ds,DS_MAT_A,&A));
    PetscCall(MatZeroEntries(A));
    PetscCall(BVOrthogonalize(svd->V,A));
    PetscCall(DSRestoreMat(svd->ds,DS_MAT_A,&A));
    PetscCall(DSSetState(svd->ds,DS_STATE_RAW));
    PetscCall(DSSolve(svd->ds,w,NULL));
    PetscCall(DSSort(svd->ds,w,NULL,NULL,NULL,NULL));
    PetscCall(DSSynchronize(svd->ds,w,NULL));
    PetscCall(DSGetMat(svd->ds,DS_MAT_U,&U));
    PetscCall(DSGetMat(svd->ds,DS_MAT_V,&V));
    PetscCall(BVMultInPlace(svd->U,V,svd->nconv,svd->ncv));
    PetscCall(BVMultInPlace(svd->V,U,svd->nconv,svd->ncv));
    PetscCall(DSRestoreMat(svd->ds,DS_MAT_U,&U));
    PetscCall(DSRestoreMat(svd->ds,DS_MAT_V,&V));
    /* Check convergence */
    k = 0;
    for (i=svd->nconv;i<svd->ncv;i++) {
      PetscCall(SVDRandomizedResidualNorm(svd,i,w[i],&res));
      svd->sigma[i] = PetscRealPart(w[i]);
      PetscCall((*svd->converged)(svd,svd->sigma[i],res,&svd->errest[i],svd->convergedctx));
      if (svd->errest[i] < svd->tol) k++;
      else break;
    }
    if (svd->conv == SVD_CONV_MAXIT && svd->its >= svd->max_it) {
      k = svd->nsv;
      for (i=0;i<svd->ncv;i++) svd->sigma[i] = PetscRealPart(w[i]);
    }
    PetscCall((*svd->stopping)(svd,svd->its,svd->max_it,svd->nconv+k,svd->nsv,&svd->reason,svd->stoppingctx));
    svd->nconv += k;
    PetscCall(SVDMonitor(svd,svd->its,svd->nconv,svd->sigma,svd->errest,svd->ncv));
  } while (svd->reason == SVD_CONVERGED_ITERATING);
  PetscCall(PetscFree(w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_Randomized(SVD svd)
{
  PetscFunctionBegin;
  svd->ops->setup          = SVDSetUp_Randomized;
  svd->ops->solve          = SVDSolve_Randomized;
  PetscFunctionReturn(PETSC_SUCCESS);
}
