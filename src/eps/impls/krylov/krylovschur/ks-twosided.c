/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "krylovschur"

   Method: Two-sided Arnoldi with Krylov-Schur restart (for left eigenvectors)

   References:

       [1] I.N. Zwaan and M.E. Hochstenbach, "Krylov-Schur-type restarts
           for the two-sided Arnoldi method", SIAM J. Matrix Anal. Appl.
           38(2):297-321, 2017.

*/

#include <slepc/private/epsimpl.h>
#include "krylovschur.h"
#include <slepcblaslapack.h>

static PetscErrorCode EPSTwoSidedRQUpdate1(EPS eps,Mat M,PetscInt nv,PetscReal beta,PetscReal betat)
{
  PetscScalar       *T,*S,*A,*w;
  const PetscScalar *pM;
  Vec               u;
  PetscInt          ld,ncv=eps->ncv,i,l,nnv;
  PetscBLASInt      n_,ncv_,*p,one=1;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(PetscMalloc3(nv,&p,ncv*ncv,&A,ncv,&w));
  PetscCall(BVGetActiveColumns(eps->V,&l,&nnv));
  PetscCall(BVSetActiveColumns(eps->V,0,nv));
  PetscCall(BVSetActiveColumns(eps->W,0,nv));
  PetscCall(BVGetColumn(eps->V,nv,&u));
  PetscCall(BVDotVec(eps->W,u,w));
  PetscCall(BVRestoreColumn(eps->V,nv,&u));
  PetscCall(MatDenseGetArrayRead(M,&pM));
  PetscCall(PetscArraycpy(A,pM,ncv*ncv));
  PetscCall(MatDenseRestoreArrayRead(M,&pM));
  PetscCall(PetscBLASIntCast(nv,&n_));
  PetscCall(PetscBLASIntCast(ncv,&ncv_));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallLAPACKInfo("LAPACKgetrf",LAPACKgetrf_(&n_,&n_,A,&ncv_,p,&info));
  PetscCall(PetscLogFlops(2.0*n_*n_*n_/3.0));
  PetscCallLAPACKInfo("LAPACKgetrs",LAPACKgetrs_("N",&n_,&one,A,&ncv_,p,w,&ncv_,&info));
  PetscCall(PetscLogFlops(2.0*n_*n_-n_));
  PetscCall(BVMultColumn(eps->V,-1.0,1.0,nv,w));
  PetscCall(DSGetArray(eps->ds,DS_MAT_A,&S));
  for (i=0;i<nv;i++) S[(nv-1)*ld+i] += beta*w[i];
  PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&S));
  PetscCall(BVGetColumn(eps->W,nv,&u));
  PetscCall(BVDotVec(eps->V,u,w));
  PetscCall(BVRestoreColumn(eps->W,nv,&u));
  PetscCallLAPACKInfo("LAPACKgetrs",LAPACKgetrs_("C",&n_,&one,A,&ncv_,p,w,&ncv_,&info));
  PetscCall(PetscFPTrapPop());
  PetscCall(BVMultColumn(eps->W,-1.0,1.0,nv,w));
  PetscCall(DSGetArray(eps->ds,DS_MAT_B,&T));
  for (i=0;i<nv;i++) T[(nv-1)*ld+i] += betat*w[i];
  PetscCall(DSRestoreArray(eps->ds,DS_MAT_B,&T));
  PetscCall(PetscFree3(p,A,w));
  PetscCall(BVSetActiveColumns(eps->V,l,nnv));
  PetscCall(BVSetActiveColumns(eps->W,l,nnv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSTwoSidedRQUpdate2(EPS eps,Mat M,PetscInt k)
{
  PetscScalar    *Q,*pM,*w,zero=0.0,sone=1.0,*c,*A;
  PetscBLASInt   n_,ncv_,ld_;
  PetscReal      norm;
  PetscInt       l,nv,ncv=eps->ncv,ld,i,j;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(BVGetActiveColumns(eps->V,&l,&nv));
  PetscCall(BVSetActiveColumns(eps->V,0,nv));
  PetscCall(BVSetActiveColumns(eps->W,0,nv));
  PetscCall(PetscMalloc2(ncv*ncv,&w,ncv,&c));
  /* u = u - V*V'*u */
  PetscCall(BVOrthogonalizeColumn(eps->V,k,c,&norm,NULL));
  PetscCall(BVScaleColumn(eps->V,k,1.0/norm));
  PetscCall(DSGetArray(eps->ds,DS_MAT_A,&A));
  /* H = H + V'*u*b' */
  for (j=l;j<k;j++) {
    for (i=0;i<k;i++) A[i+j*ld] += c[i]*A[k+j*ld];
    A[k+j*ld] *= norm;
  }
  PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&A));
  PetscCall(BVOrthogonalizeColumn(eps->W,k,c,&norm,NULL));
  PetscCall(BVScaleColumn(eps->W,k,1.0/norm));
  PetscCall(DSGetArray(eps->ds,DS_MAT_B,&A));
  /* H = H + V'*u*b' */
  for (j=l;j<k;j++) {
    for (i=0;i<k;i++) A[i+j*ld] += c[i]*A[k+j*ld];
    A[k+j*ld] *= norm;
  }
  PetscCall(DSRestoreArray(eps->ds,DS_MAT_B,&A));

  /* M = Q'*M*Q */
  PetscCall(MatDenseGetArray(M,&pM));
  PetscCall(PetscBLASIntCast(ncv,&ncv_));
  PetscCall(PetscBLASIntCast(nv,&n_));
  PetscCall(PetscBLASIntCast(ld,&ld_));
  PetscCall(DSGetArray(eps->ds,DS_MAT_Q,&Q));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,pM,&ncv_,Q,&ld_,&zero,w,&ncv_));
  PetscCall(DSRestoreArray(eps->ds,DS_MAT_Q,&Q));
  PetscCall(DSGetArray(eps->ds,DS_MAT_Z,&Q));
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,Q,&ld_,w,&ncv_,&zero,pM,&ncv_));
  PetscCall(DSRestoreArray(eps->ds,DS_MAT_Z,&Q));
  PetscCall(MatDenseRestoreArray(M,&pM));
  PetscCall(PetscFree2(w,c));
  PetscCall(BVSetActiveColumns(eps->V,l,nv));
  PetscCall(BVSetActiveColumns(eps->W,l,nv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_TwoSided(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  Mat             M,U,Op,OpHT,S,T;
  PetscReal       norm,norm2,beta,betat;
  PetscInt        ld,l,nv,nvt,k,nconv,dsn,dsk;
  PetscBool       breakdownt,breakdown,breakdownl;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(EPSGetStartVector(eps,0,NULL));
  PetscCall(EPSGetLeftStartVector(eps,0,NULL));
  l = 0;
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,eps->ncv,eps->ncv,NULL,&M));

  PetscCall(STGetOperator(eps->st,&Op));
  PetscCall(MatCreateHermitianTranspose(Op,&OpHT));

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Arnoldi factorization for Op */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSGetMat(eps->ds,DS_MAT_A,&S));
    PetscCall(BVMatArnoldi(eps->V,Op,S,eps->nconv+l,&nv,&beta,&breakdown));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&S));

    /* Compute an nv-step Arnoldi factorization for Op' */
    nvt = nv;
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSGetMat(eps->ds,DS_MAT_B,&T));
    PetscCall(BVMatArnoldi(eps->W,OpHT,T,eps->nconv+l,&nvt,&betat,&breakdownt));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_B,&T));

    /* Make sure both factorizations have the same length */
    nv = PetscMin(nv,nvt);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    if (l==0) PetscCall(DSSetState(eps->ds,DS_STATE_INTERMEDIATE));
    else PetscCall(DSSetState(eps->ds,DS_STATE_RAW));
    breakdown = (breakdown || breakdownt)? PETSC_TRUE: PETSC_FALSE;

    /* Update M, modify Rayleigh quotients S and T */
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv+l,nv));
    PetscCall(BVSetActiveColumns(eps->W,eps->nconv+l,nv));
    PetscCall(BVMatProject(eps->V,NULL,eps->W,M));

    PetscCall(EPSTwoSidedRQUpdate1(eps,M,nv,beta,betat));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSUpdateExtraRow(eps->ds));

    /* Check convergence */
    PetscCall(BVNormColumn(eps->V,nv,NORM_2,&norm));
    PetscCall(BVNormColumn(eps->W,nv,NORM_2,&norm2));
    PetscCall(EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta*norm,betat*norm2,1.0,&k));
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
      PetscCall(DSGetTruncateSize(eps->ds,k,nv,&l));
    }
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(eps,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    /* Update the corresponding vectors V(:,idx) = V*Q(:,idx) */
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));
    PetscCall(BVSetActiveColumns(eps->W,eps->nconv,nv));
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&U));
    PetscCall(BVMultInPlace(eps->V,U,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&U));
    PetscCall(DSGetMat(eps->ds,DS_MAT_Z,&U));
    PetscCall(BVMultInPlace(eps->W,U,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Z,&U));
    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      PetscCall(BVCopyColumn(eps->V,nv,k+l));
      PetscCall(BVCopyColumn(eps->W,nv,k+l));
    }

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      if (breakdown || k==nv) {
        /* Start a new Arnoldi factorization */
        PetscCall(PetscInfo(eps,"Breakdown in Krylov-Schur method (it=%" PetscInt_FMT " norm=%g)\n",eps->its,(double)beta));
        if (k<eps->nev) {
          PetscCall(EPSGetStartVector(eps,k,&breakdown));
          PetscCall(EPSGetLeftStartVector(eps,k,&breakdownl));
          if (breakdown || breakdownl) {
            eps->reason = EPS_DIVERGED_BREAKDOWN;
            PetscCall(PetscInfo(eps,"Unable to generate more start vectors\n"));
          }
        }
      } else {
        PetscCall(DSGetDimensions(eps->ds,&dsn,NULL,&dsk,NULL));
        PetscCall(DSSetDimensions(eps->ds,dsn,k,dsk));
        PetscCall(DSTruncate(eps->ds,k+l,PETSC_FALSE));
      }
      PetscCall(EPSTwoSidedRQUpdate2(eps,M,k+l));
    }
    eps->nconv = k;
    PetscCall(EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv));
  }

  PetscCall(STRestoreOperator(eps->st,&Op));
  PetscCall(MatDestroy(&OpHT));

  PetscCall(DSTruncate(eps->ds,eps->nconv,PETSC_TRUE));
  PetscCall(MatDestroy(&M));
  PetscFunctionReturn(PETSC_SUCCESS);
}
