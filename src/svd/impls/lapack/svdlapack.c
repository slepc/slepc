/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2021, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the LAPACK SVD subroutines
*/

#include <slepc/private/svdimpl.h>
#include <slepcblaslapack.h>

PetscErrorCode SVDSetUp_LAPACK(SVD svd)
{
  PetscInt       M,N,P=0;

  PetscFunctionBegin;
  PetscCall(MatGetSize(svd->A,&M,&N));
  if (!svd->isgeneralized) svd->ncv = N;
  else {
    PetscCall(MatGetSize(svd->OPb,&P,NULL));
    svd->ncv = PetscMin(M,PetscMin(N,P));
  }
  if (svd->mpd!=PETSC_DEFAULT) PetscCall(PetscInfo(svd,"Warning: parameter mpd ignored\n"));
  SVDCheckUnsupported(svd,SVD_FEATURE_STOPPING);
  if (svd->max_it==PETSC_DEFAULT) svd->max_it = 1;
  svd->leftbasis = PETSC_TRUE;
  PetscCall(SVDAllocateSolution(svd,0));
  PetscCall(DSSetType(svd->ds,svd->isgeneralized?DSGSVD:DSSVD));
  PetscCall(DSAllocate(svd->ds,PetscMax(N,PetscMax(M,P))));
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSolve_LAPACK(SVD svd)
{
  PetscInt          M,N,n,i,j,k,ld,lowu,lowv,highu,highv;
  Mat               Ar,mat;
  Vec               u,v;
  PetscScalar       *pU,*pV,*pu,*pv,*A,*w;
  const PetscScalar *pmat;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(svd->ds,&ld));
  PetscCall(MatCreateRedundantMatrix(svd->OP,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Ar));
  PetscCall(MatConvert(Ar,MATSEQDENSE,MAT_INITIAL_MATRIX,&mat));
  PetscCall(MatDestroy(&Ar));
  PetscCall(MatGetSize(mat,&M,&N));
  PetscCall(DSSetDimensions(svd->ds,M,0,0));
  PetscCall(DSSVDSetDimensions(svd->ds,N));
  PetscCall(MatDenseGetArrayRead(mat,&pmat));
  PetscCall(DSGetArray(svd->ds,DS_MAT_A,&A));
  for (i=0;i<M;i++)
    for (j=0;j<N;j++)
      A[i+j*ld] = pmat[i+j*M];
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_A,&A));
  PetscCall(MatDenseRestoreArrayRead(mat,&pmat));
  PetscCall(DSSetState(svd->ds,DS_STATE_RAW));

  n = PetscMin(M,N);
  PetscCall(PetscMalloc1(n,&w));
  PetscCall(DSSolve(svd->ds,w,NULL));
  PetscCall(DSSort(svd->ds,w,NULL,NULL,NULL,NULL));
  PetscCall(DSSynchronize(svd->ds,w,NULL));

  /* copy singular vectors */
  PetscCall(DSGetArray(svd->ds,DS_MAT_U,&pU));
  PetscCall(DSGetArray(svd->ds,DS_MAT_V,&pV));
  for (i=0;i<n;i++) {
    if (svd->which == SVD_SMALLEST) k = n - i - 1;
    else k = i;
    svd->sigma[k] = PetscRealPart(w[i]);
    PetscCall(BVGetColumn(svd->U,k,&u));
    PetscCall(BVGetColumn(svd->V,k,&v));
    PetscCall(VecGetOwnershipRange(u,&lowu,&highu));
    PetscCall(VecGetOwnershipRange(v,&lowv,&highv));
    PetscCall(VecGetArray(u,&pu));
    PetscCall(VecGetArray(v,&pv));
    if (M>=N) {
      for (j=lowu;j<highu;j++) pu[j-lowu] = pU[i*ld+j];
      for (j=lowv;j<highv;j++) pv[j-lowv] = pV[i*ld+j];
    } else {
      for (j=lowu;j<highu;j++) pu[j-lowu] = pV[i*ld+j];
      for (j=lowv;j<highv;j++) pv[j-lowv] = pU[i*ld+j];
    }
    PetscCall(VecRestoreArray(u,&pu));
    PetscCall(VecRestoreArray(v,&pv));
    PetscCall(BVRestoreColumn(svd->U,k,&u));
    PetscCall(BVRestoreColumn(svd->V,k,&v));
  }
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_U,&pU));
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_V,&pV));

  svd->nconv  = n;
  svd->its    = 1;
  svd->reason = SVD_CONVERGED_TOL;

  PetscCall(MatDestroy(&mat));
  PetscCall(PetscFree(w));
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSolve_LAPACK_GSVD(SVD svd)
{
  PetscInt          nsv,m,n,p,i,j,mlocal,plocal,ld,lowx,lowu,lowv,highx;
  Mat               Ar,A,Br,B;
  Vec               uv,x;
  PetscScalar       *Ads,*Bds,*U,*V,*X,*px,*puv,*w;
  const PetscScalar *pA,*pB;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(svd->ds,&ld));
  PetscCall(MatCreateRedundantMatrix(svd->OP,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Ar));
  PetscCall(MatConvert(Ar,MATSEQDENSE,MAT_INITIAL_MATRIX,&A));
  PetscCall(MatDestroy(&Ar));
  PetscCall(MatCreateRedundantMatrix(svd->OPb,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Br));
  PetscCall(MatConvert(Br,MATSEQDENSE,MAT_INITIAL_MATRIX,&B));
  PetscCall(MatDestroy(&Br));
  PetscCall(MatGetSize(A,&m,&n));
  PetscCall(MatGetLocalSize(svd->OP,&mlocal,NULL));
  PetscCall(MatGetLocalSize(svd->OPb,&plocal,NULL));
  PetscCall(MatGetSize(B,&p,NULL));
  PetscCall(DSSetDimensions(svd->ds,m,0,0));
  PetscCall(DSGSVDSetDimensions(svd->ds,n,p));
  PetscCall(MatDenseGetArrayRead(A,&pA));
  PetscCall(MatDenseGetArrayRead(B,&pB));
  PetscCall(DSGetArray(svd->ds,DS_MAT_A,&Ads));
  PetscCall(DSGetArray(svd->ds,DS_MAT_B,&Bds));
  for (j=0;j<n;j++) {
    for (i=0;i<m;i++) Ads[i+j*ld] = pA[i+j*m];
    for (i=0;i<p;i++) Bds[i+j*ld] = pB[i+j*p];
  }
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_B,&Bds));
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_A,&Ads));
  PetscCall(MatDenseRestoreArrayRead(B,&pB));
  PetscCall(MatDenseRestoreArrayRead(A,&pA));
  PetscCall(DSSetState(svd->ds,DS_STATE_RAW));

  nsv  = PetscMin(n,PetscMin(p,m));
  PetscCall(PetscMalloc1(nsv,&w));
  PetscCall(DSSolve(svd->ds,w,NULL));
  PetscCall(DSSort(svd->ds,w,NULL,NULL,NULL,NULL));
  PetscCall(DSSynchronize(svd->ds,w,NULL));
  PetscCall(DSGetDimensions(svd->ds,NULL,NULL,NULL,&nsv));

  /* copy singular vectors */
  PetscCall(MatGetOwnershipRange(svd->OP,&lowu,NULL));
  PetscCall(MatGetOwnershipRange(svd->OPb,&lowv,NULL));
  PetscCall(DSGetArray(svd->ds,DS_MAT_X,&X));
  PetscCall(DSGetArray(svd->ds,DS_MAT_U,&U));
  PetscCall(DSGetArray(svd->ds,DS_MAT_V,&V));
  for (j=0;j<nsv;j++) {
    svd->sigma[j] = PetscRealPart(w[j]);
    PetscCall(BVGetColumn(svd->V,j,&x));
    PetscCall(VecGetOwnershipRange(x,&lowx,&highx));
    PetscCall(VecGetArrayWrite(x,&px));
    for (i=lowx;i<highx;i++) px[i-lowx] = X[i+j*ld];
    PetscCall(VecRestoreArrayWrite(x,&px));
    PetscCall(BVRestoreColumn(svd->V,j,&x));
    PetscCall(BVGetColumn(svd->U,j,&uv));
    PetscCall(VecGetArrayWrite(uv,&puv));
    for (i=0;i<mlocal;i++) puv[i] = U[i+lowu+j*ld];
    for (i=0;i<plocal;i++) puv[i+mlocal] = V[i+lowv+j*ld];
    PetscCall(VecRestoreArrayWrite(uv,&puv));
    PetscCall(BVRestoreColumn(svd->U,j,&uv));
  }
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_X,&X));
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_U,&U));
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_V,&V));

  svd->nconv  = nsv;
  svd->its    = 1;
  svd->reason = SVD_CONVERGED_TOL;

  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(PetscFree(w));
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_LAPACK(SVD svd)
{
  PetscFunctionBegin;
  svd->ops->setup   = SVDSetUp_LAPACK;
  svd->ops->solve   = SVDSolve_LAPACK;
  svd->ops->solveg  = SVDSolve_LAPACK_GSVD;
  PetscFunctionReturn(0);
}
