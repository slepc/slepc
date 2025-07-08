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

   Method: Hamiltonian Krylov-Schur

   References:

       [1] P. Benner et al., "A Hamiltonian Krylov-Schur-type method based on the
           symplectic Lanczos process", Linear Algebra Appl. 435(3), 2011.

       [2] D. Watkins, "The Matrix Eigenvalue Problem", SIAM, 2007.

*/
#include <slepc/private/epsimpl.h>
#include "krylovschur.h"

static PetscErrorCode Orthog_Hamilt(Vec x1,Vec x2,BV U1,BV U2,BV V1,BV V2,PetscInt j,PetscScalar *c,PetscScalar *d,PetscScalar *w,PetscReal *alpha,PetscBool *breakdown)
{
  PetscInt i;
  Vec      Jx1,Jx2;

  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(U1,0,j));
  PetscCall(BVSetActiveColumns(U2,0,j));
  PetscCall(BVSetActiveColumns(V1,0,j));
  PetscCall(BVSetActiveColumns(V2,0,j));
  PetscCall(VecDuplicate(x1,&Jx1));
  PetscCall(VecDuplicate(x2,&Jx2));
  PetscCall(VecCopy(x2,Jx1));
  PetscCall(VecCopy(x1,Jx2));
  PetscCall(VecScale(Jx2,-1.0));
  PetscCall(VecConjugate(Jx1));
  PetscCall(VecConjugate(Jx2));
  /* c = -V.'*J*x  computed as -conj(c) = V'*conj(J*x) */
  PetscCall(BVDotVecBegin(V1,Jx1,c));
  PetscCall(BVDotVecBegin(V2,Jx2,w));
  PetscCall(BVDotVecEnd(V1,Jx1,c));
  PetscCall(BVDotVecEnd(V2,Jx2,w));
  for (i=0;i<j;i++) c[i] = -PetscConj(c[i]+w[i]);
  /* e = U.'*J*x  computed as conj(c) = U'*conj(J*x) */
  PetscCall(BVDotVecBegin(U1,Jx1,d));
  PetscCall(BVDotVecBegin(U2,Jx2,w));
  PetscCall(BVDotVecEnd(U1,Jx1,d));
  PetscCall(BVDotVecEnd(U2,Jx2,w));
  for (i=0;i<j;i++) d[i] = PetscConj(d[i]+w[i]);
  /* x = x-U*c-V*d */
  PetscCall(BVMultVec(U1,-1.0,1.0,x1,c));
  PetscCall(BVMultVec(U2,-1.0,1.0,x2,c));
  PetscCall(BVMultVec(V1,-1.0,1.0,x1,d));
  PetscCall(BVMultVec(V2,-1.0,1.0,x2,d));
  *alpha += PetscRealPart(c[j-1]);
  PetscCall(VecDestroy(&Jx1));
  PetscCall(VecDestroy(&Jx2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* J-orthogonalize vector x against first j vectors in U and V */
static PetscErrorCode OrthogonalizeVector_Hamilt(Vec x1,Vec x2,BV U1,BV U2,BV V1,BV V2,PetscInt j,PetscReal *alpha,PetscReal *beta,PetscScalar *h,PetscBool *breakdown)
{
  PetscScalar p,p1,p2;
  Vec         u1,u2,v1,v2;

  PetscFunctionBegin;
  PetscCall(PetscArrayzero(h,3*j));
  /* compute p=-V(:,j).'*J*x */
  PetscCall(BVGetColumn(V1,j-1,&v1));
  PetscCall(BVGetColumn(V2,j-1,&v2));
  PetscCall(VecTDotBegin(x1,v2,&p1));
  PetscCall(VecTDotBegin(x2,v1,&p2));
  PetscCall(VecTDotEnd(x1,v2,&p1));
  PetscCall(VecTDotEnd(x2,v1,&p2));
  PetscCall(BVRestoreColumn(V1,j-1,&v1));
  PetscCall(BVRestoreColumn(V2,j-1,&v2));
  p = p1-p2;
  /* update x = x - U(:,j)*p */
  PetscCall(BVGetColumn(U1,j-1,&u1));
  PetscCall(BVGetColumn(U2,j-1,&u2));
  PetscCall(VecAXPY(x1,-p,u1));
  PetscCall(VecAXPY(x2,-p,u2));
  PetscCall(BVRestoreColumn(U1,j-1,&u1));
  PetscCall(BVRestoreColumn(U2,j-1,&u2));
  if (j>1) {   /* u = u - U(:,j-1)*bb(j-1) */
    PetscCall(BVGetColumn(U1,j-2,&u1));
    PetscCall(BVGetColumn(U2,j-2,&u2));
    PetscCall(VecAXPY(x1,-beta[j-2],u1));
    PetscCall(VecAXPY(x2,-beta[j-2],u2));
    PetscCall(BVRestoreColumn(U1,j-2,&u1));
    PetscCall(BVRestoreColumn(U2,j-2,&u2));
  }
  alpha[j-1] = PetscRealPart(p);
  /* J-orthogonalize x */
  PetscCall(Orthog_Hamilt(x1,x2,U1,U2,V1,V2,j,h,h+j,h+2*j,&alpha[j-1],breakdown));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Normalize u,v so that u.'*J*v=1 */
static PetscErrorCode NormalizeVector_Hamilt(Vec u1,Vec u2,Vec v1,Vec v2,PetscReal *rr,PetscReal *ss,PetscBool *breakdown)
{
  PetscScalar p,p1,p2;
  PetscReal   r,s;

  PetscFunctionBegin;
  PetscCall(VecTDotBegin(v1,u2,&p1));
  PetscCall(VecTDotBegin(v2,u1,&p2));
  PetscCall(VecTDotEnd(v1,u2,&p1));
  PetscCall(VecTDotEnd(v2,u1,&p2));
  p = p2-p1;
  if (breakdown) {
    *breakdown = p==0.0? PETSC_TRUE: PETSC_FALSE;
    if (*breakdown) PetscFunctionReturn(PETSC_SUCCESS);
  }
  r = PetscSqrtReal(PetscAbsScalar(p));
  s = PetscSign(PetscRealPart(p));
  PetscCall(VecScale(u1,1.0/r));
  PetscCall(VecScale(u2,1.0/r));
  PetscCall(VecScale(v1,s/r));
  PetscCall(VecScale(v2,s/r));
  if (rr) *rr = r;
  if (ss) *ss = s;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSHamiltonianKS(EPS eps,BV U,BV V,PetscReal *alpha,PetscReal *beta,PetscReal *omega,PetscInt k,PetscInt *M,PetscBool *breakdown)
{
  PetscInt       j,m = *M;
  Vec            u,v,u1,u2,v1,v2;
  Mat            H;
  IS             is[2];
  BV             U1,U2,V1,V2;
  PetscScalar    *hwork,lhwork[100];

  PetscFunctionBegin;
  if (3*m > 100) PetscCall(PetscMalloc1(3*m,&hwork));
  else hwork = lhwork;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));
  PetscCall(BVGetSplitRows(U,is[0],is[1],&U1,&U2));
  PetscCall(BVGetSplitRows(V,is[0],is[1],&V1,&V2));

  /* normalize initial vector */
  if (k==0) {
    if (eps->nini==0) PetscCall(BVSetRandomColumn(eps->V,0));
    PetscCall(BVGetColumn(U,0,&u));
    PetscCall(BVGetColumn(V,0,&v));
    PetscCall(STApply(eps->st,u,v));
    PetscCall(BVRestoreColumn(U,0,&u));
    PetscCall(BVRestoreColumn(V,0,&v));
    PetscCall(BVGetColumn(U1,0,&u1));
    PetscCall(BVGetColumn(U2,0,&u2));
    PetscCall(BVGetColumn(V1,0,&v1));
    PetscCall(BVGetColumn(V2,0,&v2));
    PetscCall(NormalizeVector_Hamilt(u1,u2,v1,v2,NULL,&omega[0],breakdown));
    PetscCheck(!*breakdown,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Breakdown in Hamiltonian Krylov-Schur");
    PetscCall(BVRestoreColumn(U1,0,&u1));
    PetscCall(BVRestoreColumn(U2,0,&u2));
    PetscCall(BVRestoreColumn(V1,0,&v1));
    PetscCall(BVRestoreColumn(V2,0,&v2));
  }

  for (j=k;j<m;j++) {
    PetscCall(BVGetColumn(U,j+1,&u));
    PetscCall(BVGetColumn(V,j,&v));
    PetscCall(STApply(eps->st,v,u));
    PetscCall(BVRestoreColumn(V,j,&v));
    PetscCall(BVGetColumn(U1,j+1,&u1));
    PetscCall(BVGetColumn(U2,j+1,&u2));
    PetscCall(OrthogonalizeVector_Hamilt(u1,u2,U1,U2,V1,V2,j+1,alpha,beta,hwork,breakdown));
    PetscCall(BVGetColumn(V,j+1,&v));
    PetscCall(STApply(eps->st,u,v));
    PetscCall(BVRestoreColumn(U,j+1,&u));
    PetscCall(BVRestoreColumn(V,j+1,&v));
    PetscCall(BVGetColumn(V1,j+1,&v1));
    PetscCall(BVGetColumn(V2,j+1,&v2));
    PetscCall(NormalizeVector_Hamilt(u1,u2,v1,v2,&beta[j],&omega[j+1],breakdown));
    PetscCheck(!*breakdown,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Breakdown in Hamiltonian Krylov-Schur");
    PetscCall(BVRestoreColumn(U1,j+1,&u1));
    PetscCall(BVRestoreColumn(U2,j+1,&u2));
    PetscCall(BVRestoreColumn(V1,j+1,&v1));
    PetscCall(BVRestoreColumn(V2,j+1,&v2));
  }
  PetscCall(BVRestoreSplitRows(U,is[0],is[1],&U1,&U2));
  PetscCall(BVRestoreSplitRows(V,is[0],is[1],&V1,&V2));
  if (3*m > 100) PetscCall(PetscFree(hwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_Hamilt(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,k,l,ld,nv,t,nconv=0,nevsave;
  Mat             Q,W,D;
  Vec             vomega,vomegaold;
  BV              U,V;
  PetscReal       *a,*b,*omega,beta,u_norm;
  PetscBool       breakdown=PETSC_FALSE;
  PetscComplex    eig;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));

  /* Get the split bases */
  PetscCall(BVSetActiveColumns(eps->V,eps->ncv/2+1,eps->ncv+2));
  PetscCall(BVGetSplit(eps->V,&U,&V));

  nevsave  = eps->nev;
  eps->nev = (eps->nev+1)/2;
  l = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Krylov factorization */
    nv = PetscMin(eps->nconv+eps->mpd/2,eps->ncv/2);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&a));
    b = a + ld;
    PetscCall(DSGetArrayReal(eps->ds,DS_MAT_D,&omega));
    PetscCall(EPSHamiltonianKS(eps,U,V,a,b,omega,eps->nconv+l,&nv,&breakdown));
    beta = b[nv-1];
    PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&a));
    PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_D,&omega));
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSSetState(eps->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(U,eps->nconv,nv));
    PetscCall(BVSetActiveColumns(V,eps->nconv,nv));

    /* Save a copy of the signature */
    PetscCall(DSGetMatAndColumn(eps->ds,DS_MAT_D,0,&D,&vomega));
    PetscCall(VecDuplicate(vomega,&vomegaold));
    PetscCall(VecCopy(vomega,vomegaold));
    PetscCall(DSRestoreMatAndColumn(eps->ds,DS_MAT_D,0,&D,&vomega));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(eps->ds));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

    /* Check convergence */
    PetscCall(DSGetDimensions(eps->ds,NULL,NULL,NULL,&t));
    PetscCall(BVNormColumn(U,nv,NORM_2,&u_norm));
    PetscCall(EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta,0.0,u_norm,&k));
    EPSSetCtxThreshold(eps,eps->eigr,eps->eigi,k);
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
      l = PetscMin(l,t);
      PetscCall(DSGetTruncateSize(eps->ds,k,t,&l));
    }
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(eps,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      PetscCheck(!breakdown,PetscObjectComm((PetscObject)eps),PETSC_ERR_CONV_FAILED,"Breakdown in Hamiltonian Krylov-Schur (beta=%g)",(double)beta);
      /* Prepare the Rayleigh quotient for restart */
      PetscCall(DSTruncate(eps->ds,k+l,PETSC_FALSE));
    }
    /* Update the corresponding vectors
       U(:,idx) = U*Omega*Q(:,idx)*Omega2
       V(:,idx) = V*Q(:,idx) */
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&Q));
    PetscCall(MatDuplicate(Q,MAT_COPY_VALUES,&W));
    PetscCall(DSGetMatAndColumn(eps->ds,DS_MAT_D,0,&D,&vomega));
    PetscCall(MatDiagonalScale(W,vomegaold,vomega));
    PetscCall(DSRestoreMatAndColumn(eps->ds,DS_MAT_D,0,&D,&vomega));
    PetscCall(VecDestroy(&vomegaold));
    PetscCall(BVMultInPlace(U,W,eps->nconv,k+l));
    PetscCall(MatDestroy(&W));
    PetscCall(BVMultInPlace(V,Q,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&Q));

    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      PetscCall(BVCopyColumn(U,nv,k+l));
      PetscCall(BVCopyColumn(V,nv,k+l));
      if (eps->stop==EPS_STOP_THRESHOLD && nv-k<5) {  /* reallocate */
        eps->ncv = eps->mpd+k;
        PetscCall(BVRestoreSplit(eps->V,&U,&V));
        PetscCall(EPSReallocateSolution(eps,eps->ncv+2));
        PetscCall(BVSetActiveColumns(eps->V,eps->ncv/2+1,eps->ncv+2));
        PetscCall(BVGetSplit(eps->V,&U,&V));
        for (i=nv;i<eps->ncv;i++) eps->perm[i] = i;
        PetscCall(DSReallocate(eps->ds,eps->ncv/2+1));
        PetscCall(DSGetLeadingDimension(eps->ds,&ld));
      }
    }
    eps->nconv = k;
    for (i=0;i<nv;i++) {
#if defined(PETSC_USE_COMPLEX)
      eig = PetscSqrtScalar(eps->eigr[i]);
#else
      eig = PetscSqrtComplex(PetscCMPLX(eps->eigr[i],eps->eigi[i]));
#endif
      eps->eigr[i] = PetscRealPartComplex(eig);
      eps->eigi[i] = PetscImaginaryPartComplex(eig);
    }
    PetscCall(EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv));
  }

  eps->nev = nevsave;

  PetscCall(DSTruncate(eps->ds,eps->nconv,PETSC_TRUE));
  PetscCall(BVRestoreSplit(eps->V,&U,&V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_Hamilt(EPS eps)
{
  Mat         X,W,D;
  Vec         vomega,vomegacopy;
  BV          U,V;
  PetscInt    n;

  PetscFunctionBegin;
  /* Get the split bases */
  PetscCall(BVSetActiveColumns(eps->V,eps->ncv/2+1,eps->ncv+2));
  PetscCall(BVGetSplit(eps->V,&U,&V));
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));

  /* Update bases again, needed due to 2x2 blocks */
  PetscCall(DSGetDimensions(eps->ds,&n,NULL,NULL,NULL));
  PetscCall(DSGetMatAndColumn(eps->ds,DS_MAT_D,0,&D,&vomega));
  PetscCall(VecDuplicate(vomega,&vomegacopy));
  PetscCall(VecCopy(vomega,vomegacopy));
  PetscCall(DSRestoreMatAndColumn(eps->ds,DS_MAT_D,0,&D,&vomega));
  PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));
  PetscCall(DSGetMat(eps->ds,DS_MAT_X,&X));
  PetscCall(MatDuplicate(X,MAT_COPY_VALUES,&W));
  PetscCall(MatDiagonalScale(W,vomegacopy,NULL));
  PetscCall(VecDestroy(&vomegacopy));
  PetscCall(BVSetActiveColumns(U,0,n));
  PetscCall(BVMultInPlace(U,W,0,n));
  PetscCall(MatDestroy(&W));
  PetscCall(BVSetActiveColumns(V,0,n));
  PetscCall(BVMultInPlace(V,X,0,n));
  PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&X));

  PetscCall(BVRestoreSplit(eps->V,&U,&V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSetUp_KrylovSchur_Hamilt(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscBool       flg,sinvert;

  PetscFunctionBegin;
  PetscCheck(eps->problem_type==EPS_HAMILT,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Problem type should be Hamiltonian");
  EPSCheckUnsupportedCondition(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_EXTRACTION | EPS_FEATURE_BALANCE,PETSC_TRUE," with Hamiltonian structure");
  PetscCall(EPSSetDimensions_Default(eps,&eps->nev,&eps->ncv,&eps->mpd));
  PetscCheck(eps->ncv<=eps->nev+eps->mpd,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv)*((eps->stop==EPS_STOP_THRESHOLD)?10:1);

  PetscCall(PetscObjectTypeCompareAny((PetscObject)eps->st,&flg,STSINVERT,STSHIFT,""));
  PetscCheck(flg,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Hamiltonian Krylov-Schur only supports shift and shift-and-invert ST");
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&sinvert));
  PetscCheck(!sinvert,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Hamiltonian Krylov-Schur does not currently support shift-and-invert");
  if (!eps->which) {
    if (sinvert) eps->which = EPS_TARGET_MAGNITUDE;
    else eps->which = EPS_LARGEST_MAGNITUDE;
  }

  if (!ctx->keep) ctx->keep = 0.5;
  PetscCall(STSetStructured(eps->st,PETSC_FALSE));

  PetscCall(EPSAllocateSolution(eps,2));
  eps->ops->solve = EPSSolve_KrylovSchur_Hamilt;
  eps->ops->computevectors = EPSComputeVectors_Hamilt;
  PetscCall(DSSetType(eps->ds,DSGHIEP));
  PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
  PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
  PetscCall(DSAllocate(eps->ds,eps->ncv/2+1));
  PetscFunctionReturn(PETSC_SUCCESS);
}
