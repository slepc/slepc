/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the LAPACK SVD subroutines
*/

#include <slepc/private/svdimpl.h>
#include <slepcblaslapack.h>

static PetscErrorCode SVDSetUp_LAPACK(SVD svd)
{
  PetscInt       M,N,P=0;

  PetscFunctionBegin;
  if (svd->nsv==0) svd->nsv = 1;
  PetscCall(MatGetSize(svd->A,&M,&N));
  if (!svd->isgeneralized) svd->ncv = N;
  else {
    PetscCall(MatGetSize(svd->OPb,&P,NULL));
    svd->ncv = PetscMin(M,PetscMin(N,P));
  }
  if (svd->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(svd,"Warning: parameter mpd ignored\n"));
  SVDCheckIgnored(svd,SVD_FEATURE_STOPPING);
  if (svd->max_it==PETSC_DETERMINE) svd->max_it = 1;
  svd->leftbasis = PETSC_TRUE;
  PetscCall(SVDAllocateSolution(svd,0));
  PetscCall(DSAllocate(svd->ds,PetscMax(N,PetscMax(M,P))));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_LAPACK(SVD svd)
{
  PetscInt          M,N,n,i,j,k,ld,lowu,lowv,highu,highv;
  Mat               A,Ar,mat;
  Vec               u,v;
  PetscScalar       *pU,*pV,*pu,*pv,*w;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(svd->ds,&ld));
  PetscCall(MatCreateRedundantMatrix(svd->OP,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Ar));
  PetscCall(MatConvert(Ar,MATSEQDENSE,MAT_INITIAL_MATRIX,&mat));
  PetscCall(MatDestroy(&Ar));
  PetscCall(MatGetSize(mat,&M,&N));
  PetscCall(DSSetDimensions(svd->ds,M,0,0));
  PetscCall(DSSVDSetDimensions(svd->ds,N));
  PetscCall(DSGetMat(svd->ds,DS_MAT_A,&A));
  PetscCall(MatCopy(mat,A,SAME_NONZERO_PATTERN));
  PetscCall(DSRestoreMat(svd->ds,DS_MAT_A,&A));
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
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_LAPACK_GSVD(SVD svd)
{
  PetscInt          nsv,m,n,p,i,j,mlocal,plocal,ld,lowx,lowu,lowv,highx;
  Mat               Ar,A,Ads,Br,B,Bds;
  Vec               uv,x;
  PetscScalar       *U,*V,*X,*px,*puv,*w;

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
  PetscCall(DSGetMat(svd->ds,DS_MAT_A,&Ads));
  PetscCall(MatCopy(A,Ads,SAME_NONZERO_PATTERN));
  PetscCall(DSRestoreMat(svd->ds,DS_MAT_A,&Ads));
  PetscCall(DSGetMat(svd->ds,DS_MAT_B,&Bds));
  PetscCall(MatCopy(B,Bds,SAME_NONZERO_PATTERN));
  PetscCall(DSRestoreMat(svd->ds,DS_MAT_B,&Bds));
  PetscCall(DSSetState(svd->ds,DS_STATE_RAW));

  nsv = PetscMin(n,PetscMin(p,m));
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
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_LAPACK_HSVD(SVD svd)
{
  PetscInt          M,N,n,i,j,k,ld,lowu,lowv,highu,highv;
  Mat               A,Ar,mat,D;
  Vec               u,v,vomega;
  PetscScalar       *pU,*pV,*pu,*pv,*w;
  PetscReal         *pD;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(svd->ds,&ld));
  PetscCall(MatCreateRedundantMatrix(svd->OP,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Ar));
  PetscCall(MatConvert(Ar,MATSEQDENSE,MAT_INITIAL_MATRIX,&mat));
  PetscCall(MatDestroy(&Ar));
  PetscCall(MatGetSize(mat,&M,&N));
  PetscCall(DSSetDimensions(svd->ds,M,0,0));
  PetscCall(DSHSVDSetDimensions(svd->ds,N));
  PetscCall(DSGetMat(svd->ds,DS_MAT_A,&A));
  PetscCall(MatCopy(mat,A,SAME_NONZERO_PATTERN));
  PetscCall(DSRestoreMat(svd->ds,DS_MAT_A,&A));
  PetscCall(DSGetMatAndColumn(svd->ds,DS_MAT_D,0,&D,&vomega));
  PetscCall(VecCopy(svd->omega,vomega));
  PetscCall(DSRestoreMatAndColumn(svd->ds,DS_MAT_D,0,&D,&vomega));
  PetscCall(DSSetState(svd->ds,DS_STATE_RAW));

  n = PetscMin(M,N);
  PetscCall(PetscMalloc1(n,&w));
  PetscCall(DSSolve(svd->ds,w,NULL));
  PetscCall(DSSort(svd->ds,w,NULL,NULL,NULL,NULL));
  PetscCall(DSSynchronize(svd->ds,w,NULL));

  /* copy singular vectors */
  PetscCall(DSGetArrayReal(svd->ds,DS_MAT_D,&pD));
  PetscCall(DSGetArray(svd->ds,DS_MAT_U,&pU));
  PetscCall(DSGetArray(svd->ds,DS_MAT_V,&pV));
  for (i=0;i<n;i++) {
    if (svd->which == SVD_SMALLEST) k = n - i - 1;
    else k = i;
    svd->sigma[k] = PetscRealPart(w[i]);
    svd->sign[k]  = pD[i];
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
  PetscCall(DSRestoreArrayReal(svd->ds,DS_MAT_D,&pD));
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_U,&pU));
  PetscCall(DSRestoreArray(svd->ds,DS_MAT_V,&pV));

  svd->nconv  = n;
  svd->its    = 1;
  svd->reason = SVD_CONVERGED_TOL;

  PetscCall(MatDestroy(&mat));
  PetscCall(PetscFree(w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSetDSType_LAPACK(SVD svd)
{
  PetscFunctionBegin;
  PetscCall(DSSetType(svd->ds,svd->OPb?DSGSVD:svd->omega?DSHSVD:DSSVD));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_LAPACK(SVD svd)
{
  PetscFunctionBegin;
  svd->ops->setup     = SVDSetUp_LAPACK;
  svd->ops->solve     = SVDSolve_LAPACK;
  svd->ops->solveg    = SVDSolve_LAPACK_GSVD;
  svd->ops->solveh    = SVDSolve_LAPACK_HSVD;
  svd->ops->setdstype = SVDSetDSType_LAPACK;
  PetscFunctionReturn(PETSC_SUCCESS);
}
