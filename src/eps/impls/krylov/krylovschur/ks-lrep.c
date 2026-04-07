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

   Method: thick-restarted Lanczos for Linear Response eigenvalue problems

   References:

       [1] Z. Teng, R.-C. Li, "Convergence analysis of Lanczos-type methods for the
           linear response eigenvalue problem", J. Comput. Appl. Math. 247, 2013.

       [2] F. Alvarruiz, B. Mellado-Pinto, J. E. Roman, "Restarted Lanczos methods
           for the linear response eigenvalue problem", in preparation, 2026.

*/
#include <slepc/private/epsimpl.h>
#include "krylovschur.h"

static PetscErrorCode Orthog_Teng(Vec x,BV U,BV V,PetscInt j,PetscScalar *h,PetscScalar *c)
{
  PetscInt i;

  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(U,0,j));
  PetscCall(BVSetActiveColumns(V,0,j));
  /* c = U^* x */
  PetscCall(BVDotVec(U,x,c));
  /* x = x-V*c */
  PetscCall(BVMultVec(V,-1.0,1.0,x,c));
  /* accumulate orthog coeffs into h */
  for (i=0;i<2*j;i++) h[i] += c[i];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Orthogonalize vector x against first j vectors in U and V
v is column j-1 of V */
static PetscErrorCode OrthogonalizeVector_Teng(Vec x,BV U,BV V,PetscInt j,Vec u,PetscReal *beta,PetscInt k,PetscScalar *h)
{
  PetscReal alpha;
  PetscInt  i,l;

  PetscFunctionBegin;
  PetscCall(PetscArrayzero(h,2*j));

  /* Local orthogonalization */
  l = j==k+1?0:j-2;  /* 1st column to orthogonalize against */
  PetscCall(VecDotRealPart(x,u,&alpha));
  for (i=l;i<j-1;i++) h[i] = beta[i];
  h[j-1] = alpha;
  /* x = x-V(:,l:j-1)*h(l:j-1) */
  PetscCall(BVSetActiveColumns(V,l,j));
  PetscCall(BVMultVec(V,-1.0,1.0,x,h+l));

  /* Full orthogonalization */
  PetscCall(Orthog_Teng(x,U,V,j,h,h+2*j));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLREPLanczos_Teng(EPS eps,Mat K,Mat M,BV U,BV V,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *min,PetscBool *breakdown)
{
  PetscInt       j,m = *min;
  Vec            u,v,uh,vh;
  PetscReal      beta0;
  PetscScalar    *hwork,lhwork[100],gamma;

  PetscFunctionBegin;
  if (4*m > 100) PetscCall(PetscMalloc1(4*m,&hwork));
  else hwork = lhwork;

  /* Normalize initial vector */
  if (k==0) {
    if (eps->nini==0) PetscCall(BVSetRandomColumn(V,0));
    PetscCall(BVGetColumn(U,0,&u));
    PetscCall(BVGetColumn(V,0,&v));
    PetscCall(MatMult(M,v,u));
    PetscCall(VecDot(u,v,&gamma));
    beta0 = PetscSqrtReal(PetscRealPart(gamma));
    if (beta0==0.0) {
      if (breakdown) *breakdown = PETSC_TRUE;
      *min = 1; m = 0;
    } else {
      PetscCall(VecScale(u,1.0/beta0));
      PetscCall(VecScale(v,1.0/beta0));
    }
    PetscCall(BVRestoreColumn(U,0,&u));
    PetscCall(BVRestoreColumn(V,0,&v));
  }

  for (j=k;j<m;j++) {
    /* j+1 columns (indices 0 to j) have been computed */
    PetscCall(BVGetColumn(U,j+1,&uh));
    PetscCall(BVGetColumn(V,j+1,&vh));
    PetscCall(BVGetColumn(U,j,&u));
    PetscCall(MatMult(K,u,vh));
    PetscCall(OrthogonalizeVector_Teng(vh,U,V,j+1,u,beta,k,hwork));
    alpha[j] = PetscRealPart(hwork[j]);
    PetscCall(MatMult(M,vh,uh));
    PetscCall(VecDot(uh,vh,&gamma));
    beta[j] = PetscSqrtReal(PetscRealPart(gamma));
    if (beta[j]==0.0) {
      if (breakdown) *breakdown = PETSC_TRUE;
      *min = j+1; m = j;
    } else {
      PetscCall(VecScale(uh,1.0/beta[j]));
      PetscCall(VecScale(vh,1.0/beta[j]));
    }
    PetscCall(BVRestoreColumn(U,j+1,&uh));
    PetscCall(BVRestoreColumn(V,j+1,&vh));
    PetscCall(BVRestoreColumn(U,j,&u));
  }
  if (4*m > 100) PetscCall(PetscFree(hwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_LREP_Teng(EPS eps)
{
  Mat         H;
  Vec         u,v,w;
  BV          U,V;
  IS          is[2];
  PetscInt    k;
  PetscScalar lambda;
  PetscBool   reduced;

  PetscFunctionBegin;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));
  PetscCall(SlepcCheckMatLREPReduced(H,&reduced));
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&V,&U));
  for (k=0;k<eps->nconv;k++) {
    PetscCall(BVGetColumn(V,k,&v));
    /* approx eigenvector is [eigr[k]*v; u] */
    lambda = eps->eigr[k];
    PetscCall(STBackTransform(eps->st,1,&lambda,&eps->eigi[k]));
    PetscCall(VecScale(v,lambda));
    PetscCall(BVRestoreColumn(V,k,&v));
  }
  if (!reduced) {
    /* the eigenvector [v;u] = J*[y;x] where [y;x] is the reduced eigenvector
       and J = 1/sqrt(2)[I I; I -I], i.e, v=1/sqrt(2)*(y+x) u=1/sqrt(2)*(y-x) */
    PetscCall(BVCreateVec(V,&w));
    for (k=0;k<eps->nconv;k++) {
      PetscCall(BVGetColumn(U,k,&u));
      PetscCall(BVGetColumn(V,k,&v));
      PetscCall(VecCopy(u,w));
      PetscCall(VecCopy(v,u));
      PetscCall(VecAXPY(u,-1.0,w));
      PetscCall(VecAXPY(v,1.0,w));
      PetscCall(BVRestoreColumn(U,k,&u));
      PetscCall(BVRestoreColumn(V,k,&v));
    }
    PetscCall(VecDestroy(&w));
  }
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&V,&U));
  /* Normalize eigenvectors */
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
  PetscCall(BVNormalize(eps->V,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSetUp_KrylovSchur_LREP(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscBool       flg;

  PetscFunctionBegin;
  PetscCheck((eps->problem_type==EPS_LREP),PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Problem type should be LREP");
  EPSCheckUnsupportedCondition(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_EXTRACTION | EPS_FEATURE_BALANCE,PETSC_TRUE," with LREP structure");
  PetscCall(EPSSetDimensions_Default(eps,&eps->nev,&eps->ncv,&eps->mpd));
  PetscCheck(eps->ncv<=eps->nev+eps->mpd,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv)*((eps->stop==EPS_STOP_THRESHOLD)?10:1);

  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&flg));
  PetscCheck(flg,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Krylov-Schur LREP only supports shift ST");
  if (!eps->which) eps->which = EPS_SMALLEST_MAGNITUDE;

  if (!ctx->keep) ctx->keep = 0.5;
  PetscCall(STSetStructured(eps->st,PETSC_FALSE));

  PetscCall(EPSAllocateSolution(eps,1));
  /* Teng */
  eps->ops->solve = EPSSolve_KrylovSchur_LREP_Teng;
  eps->ops->computevectors = EPSComputeVectors_LREP_Teng;
  PetscCall(DSSetType(eps->ds,DSHEP));
  PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
  PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
  PetscCall(DSAllocate(eps->ds,eps->ncv+1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_LREP_Teng(EPS eps)
{
  EPS_KRYLOVSCHUR   *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt          i,k,l,ld,nv,nconv=0,nevsave,ma,na,Ma,Na;
  Mat               H,Q,K,M,A,B;
  BV                U,V;
  IS                is[2];
  PetscReal         *a,*b,beta;
  PetscBool         reduced,breakdown=PETSC_FALSE;
  const PetscScalar scal[] = { 1.0, -1.0 };

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));

  /* Extract matrix blocks */
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));
  PetscCall(SlepcCheckMatLREPReduced(H,&reduced));
  if (reduced) {
    PetscCall(MatNestGetSubMat(H,0,1,&K));
    PetscCall(MatNestGetSubMat(H,1,0,&M));
  } else {
    PetscCall(MatNestGetSubMat(H,0,0,&A));
    PetscCall(MatNestGetSubMat(H,0,1,&B));
    PetscCall(MatGetSize(A,&Ma,&Na));
    PetscCall(MatGetLocalSize(A,&ma,&na));
    /* K = A-B */
    PetscCall(MatCreate(PetscObjectComm((PetscObject)A),&K));
    PetscCall(MatSetSizes(K,ma,na,Ma,Na));
    PetscCall(MatSetType(K,MATCOMPOSITE));
    PetscCall(MatCompositeAddMat(K,A));
    PetscCall(MatCompositeAddMat(K,B));
    PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatCompositeSetScalings(K,scal));
    /* M = A+B */
    PetscCall(MatCreate(PetscObjectComm((PetscObject)A),&M));
    PetscCall(MatSetSizes(M,ma,na,Ma,Na));
    PetscCall(MatSetType(M,MATCOMPOSITE));
    PetscCall(MatCompositeAddMat(M,A));
    PetscCall(MatCompositeAddMat(M,B));
    PetscCall(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY));
  }

  /* Get the split bases */
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&V,&U));

  nevsave  = eps->nev;
  eps->nev = (eps->nev+1)/2;
  l = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&a));
    b = a + ld;
    PetscCall(EPSLREPLanczos_Teng(eps,K,M,U,V,a,b,eps->nconv+l,&nv,&breakdown));
    beta = b[nv-1];
    PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&a));
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSSetState(eps->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(eps->ds));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

    /* Check convergence */
    for (i=0;i<nv;i++) eps->eigr[i] = PetscSqrtReal(PetscRealPart(eps->eigr[i]));
    PetscCall(EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta,0.0,1.0,&k));
    EPSSetCtxThreshold(eps,eps->eigr,eps->eigi,eps->errest,k,nv);
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(eps,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      PetscCheck(!breakdown,PetscObjectComm((PetscObject)eps),PETSC_ERR_CONV_FAILED,"Breakdown in LREP Krylov-Schur (beta=%g)",(double)beta);
      /* Prepare the Rayleigh quotient for restart */
      PetscCall(DSTruncate(eps->ds,k+l,PETSC_FALSE));
    }
    /* Update the corresponding vectors
       U(:,idx) = U*Q(:,idx),  V(:,idx) = V*Q(:,idx) */
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&Q));
    PetscCall(BVMultInPlace(U,Q,eps->nconv,k+l));
    PetscCall(BVMultInPlace(V,Q,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&Q));

    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      PetscCall(BVCopyColumn(eps->V,nv,k+l));
      if (eps->stop==EPS_STOP_THRESHOLD && nv-k<5) {  /* reallocate */
        eps->ncv = eps->mpd+k;
        PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&V,&U));
        PetscCall(EPSReallocateSolution(eps,eps->ncv+1));
        PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&V,&U));
        for (i=nv;i<eps->ncv;i++) eps->perm[i] = i;
        PetscCall(DSReallocate(eps->ds,eps->ncv+1));
        PetscCall(DSGetLeadingDimension(eps->ds,&ld));
      }
    }
    eps->nconv = k;
    PetscCall(EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv));
  }

  eps->nev = nevsave;

  PetscCall(DSTruncate(eps->ds,eps->nconv,PETSC_TRUE));
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&V,&U));
  if (!reduced) {
    PetscCall(MatDestroy(&K));
    PetscCall(MatDestroy(&M));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
