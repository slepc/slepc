/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Two-sided variant of the NEPSLP solver.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include "slp.h"

typedef struct _n_nep_def_ctx *NEP_NEDEF_CTX;

struct _n_nep_def_ctx {
  PetscInt    n;
  PetscBool   ref;
  PetscScalar *eig;
  BV          V,W;
};

typedef struct {   /* context for two-sided solver */
  Mat         Ft;
  Mat         Jt;
  Vec         w;
} NEP_SLPTS_MATSHELL;

typedef struct {   /* context for non-equivalence deflation */
  NEP_NEDEF_CTX defctx;
  Mat           F;
  Mat           P;
  Mat           J;
  KSP           ksp;
  PetscBool     isJ;
  PetscScalar   lambda;
  Vec           w[2];
} NEP_NEDEF_MATSHELL;

static PetscErrorCode MatMult_SLPTS_Right(Mat M,Vec x,Vec y)
{
  NEP_SLPTS_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(MatMult(ctx->Jt,x,ctx->w));
  PetscCall(MatSolve(ctx->Ft,ctx->w,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_SLPTS_Left(Mat M,Vec x,Vec y)
{
  NEP_SLPTS_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(MatMultTranspose(ctx->Jt,x,ctx->w));
  PetscCall(MatSolveTranspose(ctx->Ft,ctx->w,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_SLPTS(Mat M)
{
  NEP_SLPTS_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(VecDestroy(&ctx->w));
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
static PetscErrorCode MatCreateVecs_SLPTS(Mat M,Vec *left,Vec *right)
{
  NEP_SLPTS_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  if (right) PetscCall(VecDuplicate(ctx->w,right));
  if (left) PetscCall(VecDuplicate(ctx->w,left));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode NEPSLPSetUpEPSMat(NEP nep,Mat F,Mat J,PetscBool left,Mat *M)
{
  Mat                Mshell;
  PetscInt           nloc,mloc;
  NEP_SLPTS_MATSHELL *shellctx;

  PetscFunctionBegin;
  /* Create mat shell */
  PetscCall(PetscNew(&shellctx));
  shellctx->Ft = F;
  shellctx->Jt = J;
  PetscCall(MatGetLocalSize(nep->function,&mloc,&nloc));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)nep),nloc,mloc,PETSC_DETERMINE,PETSC_DETERMINE,shellctx,&Mshell));
  if (left) PetscCall(MatShellSetOperation(Mshell,MATOP_MULT,(PetscErrorCodeFn*)MatMult_SLPTS_Left));
  else PetscCall(MatShellSetOperation(Mshell,MATOP_MULT,(PetscErrorCodeFn*)MatMult_SLPTS_Right));
  PetscCall(MatShellSetOperation(Mshell,MATOP_DESTROY,(PetscErrorCodeFn*)MatDestroy_SLPTS));
#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
  PetscCall(MatShellSetOperation(Mshell,MATOP_CREATE_VECS,(PetscErrorCodeFn*)MatCreateVecs_SLPTS));
#endif
  *M = Mshell;
  PetscCall(MatCreateVecs(nep->function,&shellctx->w,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Functions for deflation */
static PetscErrorCode NEPDeflationNEDestroy(NEP_NEDEF_CTX defctx)
{
  PetscFunctionBegin;
  if (!defctx) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscFree(defctx->eig));
  PetscCall(PetscFree(defctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationNECreate(NEP nep,BV V,BV W,PetscInt sz,NEP_NEDEF_CTX *defctx)
{
  NEP_NEDEF_CTX  op;

  PetscFunctionBegin;
  PetscCall(PetscNew(&op));
  *defctx = op;
  op->n   = 0;
  op->ref = PETSC_FALSE;
  PetscCall(PetscCalloc1(sz,&op->eig));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscCall(PetscObjectStateIncrease((PetscObject)W));
  op->V = V;
  op->W = W;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationNEComputeFunction(NEP nep,Mat M,PetscScalar lambda)
{
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  if (lambda==matctx->lambda) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(NEPComputeFunction(nep,lambda,matctx->F,matctx->P));
  if (matctx->isJ) PetscCall(NEPComputeJacobian(nep,lambda,matctx->J));
  if (matctx->ksp) PetscCall(NEP_KSPSetOperators(matctx->ksp,matctx->F,matctx->P));
  matctx->lambda = lambda;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_NEPDeflationNE(Mat M,Vec x,Vec r)
{
  Vec                t,tt;
  PetscScalar        *h,*alpha,lambda,*eig;
  PetscInt           i,k;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  if (matctx->defctx->n && !matctx->defctx->ref) {
    k = matctx->defctx->n;
    lambda = matctx->lambda;
    eig = matctx->defctx->eig;
    t = matctx->w[0];
    PetscCall(VecCopy(x,t));
    PetscCall(PetscMalloc2(k,&h,k,&alpha));
    for (i=0;i<k;i++) alpha[i] = (lambda-eig[i]-1.0)/(lambda-eig[i]);
    PetscCall(BVDotVec(matctx->defctx->V,t,h));
    for (i=0;i<k;i++) h[i] *= alpha[i];
    PetscCall(BVMultVec(matctx->defctx->W,-1.0,1.0,t,h));
    PetscCall(MatMult(matctx->isJ?matctx->J:matctx->F,t,r));
    if (matctx->isJ) {
      for (i=0;i<k;i++) h[i] *= (1.0/((lambda-eig[i])*(lambda-eig[i])))/alpha[i];
      tt = matctx->w[1];
      PetscCall(BVMultVec(matctx->defctx->W,-1.0,0.0,tt,h));
      PetscCall(MatMult(matctx->F,tt,t));
      PetscCall(VecAXPY(r,1.0,t));
    }
    PetscCall(PetscFree2(h,alpha));
  } else PetscCall(MatMult(matctx->isJ?matctx->J:matctx->F,x,r));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultTranspose_NEPDeflationNE(Mat M,Vec x,Vec r)
{
  Vec                t,tt;
  PetscScalar        *h,*alphaC,lambda,*eig;
  PetscInt           i,k;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  t    = matctx->w[0];
  PetscCall(VecCopy(x,t));
  if (matctx->defctx->n && !matctx->defctx->ref) {
    PetscCall(VecConjugate(t));
    k = matctx->defctx->n;
    lambda = matctx->lambda;
    eig = matctx->defctx->eig;
    PetscCall(PetscMalloc2(k,&h,k,&alphaC));
    for (i=0;i<k;i++) alphaC[i] = PetscConj((lambda-eig[i]-1.0)/(lambda-eig[i]));
    PetscCall(BVDotVec(matctx->defctx->W,t,h));
    for (i=0;i<k;i++) h[i] *= alphaC[i];
    PetscCall(BVMultVec(matctx->defctx->V,-1.0,1.0,t,h));
    PetscCall(VecConjugate(t));
    PetscCall(MatMultTranspose(matctx->isJ?matctx->J:matctx->F,t,r));
    if (matctx->isJ) {
      for (i=0;i<k;i++) h[i] *= PetscConj(1.0/((lambda-eig[i])*(lambda-eig[i])))/alphaC[i];
      tt = matctx->w[1];
      PetscCall(BVMultVec(matctx->defctx->V,-1.0,0.0,tt,h));
      PetscCall(VecConjugate(tt));
      PetscCall(MatMultTranspose(matctx->F,tt,t));
      PetscCall(VecAXPY(r,1.0,t));
    }
    PetscCall(PetscFree2(h,alphaC));
  } else PetscCall(MatMultTranspose(matctx->isJ?matctx->J:matctx->F,t,r));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatSolve_NEPDeflationNE(Mat M,Vec b,Vec x)
{
  PetscScalar        *h,*alpha,lambda,*eig;
  PetscInt           i,k;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  if (!matctx->ksp) {
    PetscCall(VecCopy(b,x));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCall(KSPSolve(matctx->ksp,b,x));
  if (matctx->defctx->n && !matctx->defctx->ref) {
    k = matctx->defctx->n;
    lambda = matctx->lambda;
    eig = matctx->defctx->eig;
    PetscCall(PetscMalloc2(k,&h,k,&alpha));
    PetscCall(BVDotVec(matctx->defctx->V,x,h));
    for (i=0;i<k;i++) alpha[i] = (lambda-eig[i]-1.0)/(lambda-eig[i]);
    for (i=0;i<k;i++) h[i] *= alpha[i]/(1.0-alpha[i]);
    PetscCall(BVMultVec(matctx->defctx->W,1.0,1.0,x,h));
    PetscCall(PetscFree2(h,alpha));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatSolveTranspose_NEPDeflationNE(Mat M,Vec b,Vec x)
{
  PetscScalar        *h,*alphaC,lambda,*eig;
  PetscInt           i,k;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  if (!matctx->ksp) {
    PetscCall(VecCopy(b,x));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCall(KSPSolveTranspose(matctx->ksp,b,x));
  if (matctx->defctx->n && !matctx->defctx->ref) {
    PetscCall(VecConjugate(x));
    k = matctx->defctx->n;
    lambda = matctx->lambda;
    eig = matctx->defctx->eig;
    PetscCall(PetscMalloc2(k,&h,k,&alphaC));
    PetscCall(BVDotVec(matctx->defctx->W,x,h));
    for (i=0;i<k;i++) alphaC[i] = PetscConj((lambda-eig[i]-1.0)/(lambda-eig[i]));
    for (i=0;i<k;i++) h[i] *= alphaC[i]/(1.0-alphaC[i]);
    PetscCall(BVMultVec(matctx->defctx->V,1.0,1.0,x,h));
    PetscCall(PetscFree2(h,alphaC));
    PetscCall(VecConjugate(x));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_NEPDeflationNE(Mat M)
{
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  PetscCall(VecDestroy(&matctx->w[0]));
  PetscCall(VecDestroy(&matctx->w[1]));
  PetscCall(PetscFree(matctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatCreateVecs_NEPDeflationNE(Mat M,Vec *right,Vec *left)
{
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  PetscCall(MatCreateVecs(matctx->F,right,left));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationNEFunctionCreate(NEP_NEDEF_CTX defctx,NEP nep,Mat F,Mat P,Mat J,KSP ksp,PetscBool isJ,Mat *Mshell)
{
  NEP_NEDEF_MATSHELL *matctx;
  PetscInt           nloc,mloc;

  PetscFunctionBegin;
  /* Create mat shell */
  PetscCall(PetscNew(&matctx));
  PetscCall(MatGetLocalSize(nep->function,&mloc,&nloc));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)nep),nloc,mloc,PETSC_DETERMINE,PETSC_DETERMINE,matctx,Mshell));
  matctx->F   = F;
  matctx->P   = P;
  matctx->J   = J;
  matctx->isJ = isJ;
  matctx->ksp = ksp;
  matctx->defctx = defctx;
  matctx->lambda = PETSC_MAX_REAL;
  PetscCall(MatCreateVecs(F,&matctx->w[0],NULL));
  PetscCall(VecDuplicate(matctx->w[0],&matctx->w[1]));
  PetscCall(MatShellSetOperation(*Mshell,MATOP_MULT,(PetscErrorCodeFn*)MatMult_NEPDeflationNE));
  PetscCall(MatShellSetOperation(*Mshell,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_NEPDeflationNE));
  PetscCall(MatShellSetOperation(*Mshell,MATOP_SOLVE,(PetscErrorCodeFn*)MatSolve_NEPDeflationNE));
  PetscCall(MatShellSetOperation(*Mshell,MATOP_SOLVE_TRANSPOSE,(PetscErrorCodeFn*)MatSolveTranspose_NEPDeflationNE));
  PetscCall(MatShellSetOperation(*Mshell,MATOP_DESTROY,(PetscErrorCodeFn*)MatDestroy_NEPDeflationNE));
  PetscCall(MatShellSetOperation(*Mshell,MATOP_CREATE_VECS,(PetscErrorCodeFn*)MatCreateVecs_NEPDeflationNE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationNERecoverEigenvectors(NEP_NEDEF_CTX defctx,Vec u,Vec w,PetscScalar lambda)
{
  PetscScalar    *h,*alpha,*eig;
  PetscInt       i,k;

  PetscFunctionBegin;
  if (w) PetscCall(VecConjugate(w));
  if (defctx->n && !defctx->ref) {
    eig = defctx->eig;
    k = defctx->n;
    PetscCall(PetscMalloc2(k,&h,k,&alpha));
    for (i=0;i<k;i++) alpha[i] = (lambda-eig[i]-1.0)/(lambda-eig[i]);
    PetscCall(BVDotVec(defctx->V,u,h));
    for (i=0;i<k;i++) h[i] *= alpha[i];
    PetscCall(BVMultVec(defctx->W,-1.0,1.0,u,h));
    PetscCall(VecNormalize(u,NULL));
    if (w) {
      PetscCall(BVDotVec(defctx->W,w,h));
      for (i=0;i<k;i++) alpha[i] = PetscConj((lambda-eig[i]-1.0)/(lambda-eig[i]));
      for (i=0;i<k;i++) h[i] *= alpha[i];
      PetscCall(BVMultVec(defctx->V,-1.0,1.0,w,h));
      PetscCall(VecNormalize(w,NULL));
    }
    PetscCall(PetscFree2(h,alpha));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationNELocking(NEP_NEDEF_CTX defctx,Vec u,Vec w,PetscScalar lambda)
{
  PetscInt       n;

  PetscFunctionBegin;
  n = defctx->n++;
  defctx->eig[n] = lambda;
  PetscCall(BVInsertVec(defctx->V,n,u));
  PetscCall(BVInsertVec(defctx->W,n,w));
  PetscCall(BVSetActiveColumns(defctx->V,0,defctx->n));
  PetscCall(BVSetActiveColumns(defctx->W,0,defctx->n));
  PetscCall(BVBiorthonormalizeColumn(defctx->V,defctx->W,n,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationNESetRefine(NEP_NEDEF_CTX defctx,PetscBool ref)
{
  PetscFunctionBegin;
  defctx->ref = ref;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPSolve_SLP_Twosided(NEP nep)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  Mat            mF,mJ,M,Mt;
  Vec            u,r,t,w;
  BV             X,Y;
  PetscScalar    sigma,lambda,mu,im=0.0,mu2,im2;
  PetscReal      resnorm,resl;
  PetscInt       nconv,nconv2,i;
  PetscBool      skip=PETSC_FALSE,lock=PETSC_FALSE;
  NEP_NEDEF_CTX  defctx=NULL;    /* Extended operator for deflation */

  PetscFunctionBegin;
  /* get initial approximation of eigenvalue and eigenvector */
  PetscCall(NEPGetDefaultShift(nep,&sigma));
  if (!nep->nini) PetscCall(BVSetRandomColumn(nep->V,0));
  PetscCall(BVSetRandomColumn(nep->W,0));
  lambda = sigma;
  if (!ctx->ksp) PetscCall(NEPSLPGetKSP(nep,&ctx->ksp));
  PetscCall(BVDuplicate(nep->V,&X));
  PetscCall(BVDuplicate(nep->W,&Y));
  PetscCall(NEPDeflationNECreate(nep,X,Y,nep->nev,&defctx));
  PetscCall(BVGetColumn(nep->V,0,&t));
  PetscCall(VecDuplicate(t,&u));
  PetscCall(VecDuplicate(t,&w));
  PetscCall(BVRestoreColumn(nep->V,0,&t));
  PetscCall(BVCopyVec(nep->V,0,u));
  PetscCall(BVCopyVec(nep->W,0,w));
  PetscCall(VecDuplicate(u,&r));
  PetscCall(NEPDeflationNEFunctionCreate(defctx,nep,nep->function,nep->function_pre?nep->function_pre:nep->function,NULL,ctx->ksp,PETSC_FALSE,&mF));
  PetscCall(NEPDeflationNEFunctionCreate(defctx,nep,nep->function,nep->function,nep->jacobian,NULL,PETSC_TRUE,&mJ));
  PetscCall(NEPSLPSetUpEPSMat(nep,mF,mJ,PETSC_FALSE,&M));
  PetscCall(NEPSLPSetUpEPSMat(nep,mF,mJ,PETSC_TRUE,&Mt));
  PetscCall(EPSSetOperators(ctx->eps,M,NULL));
  PetscCall(MatDestroy(&M));
  PetscCall(EPSSetOperators(ctx->epsts,Mt,NULL));
  PetscCall(MatDestroy(&Mt));

  /* Restart loop */
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;

    /* form residual,  r = T(lambda)*u (used in convergence test only) */
    PetscCall(NEPDeflationNEComputeFunction(nep,mF,lambda));
    PetscCall(MatMultTranspose(mF,w,r));
    PetscCall(VecNorm(r,NORM_2,&resl));
    PetscCall(MatMult(mF,u,r));

    /* convergence test */
    PetscCall(VecNorm(r,NORM_2,&resnorm));
    resnorm = PetscMax(resnorm,resl);
    PetscCall((*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx));
    nep->eigr[nep->nconv] = lambda;
    if (nep->errest[nep->nconv]<=nep->tol || nep->errest[nep->nconv]<=ctx->deftol) {
      if (nep->errest[nep->nconv]<=ctx->deftol && !defctx->ref && nep->nconv) {
        PetscCall(NEPDeflationNERecoverEigenvectors(defctx,u,w,lambda));
        PetscCall(VecConjugate(w));
        PetscCall(NEPDeflationNESetRefine(defctx,PETSC_TRUE));
        PetscCall(MatMultTranspose(mF,w,r));
        PetscCall(VecNorm(r,NORM_2,&resl));
        PetscCall(MatMult(mF,u,r));
        PetscCall(VecNorm(r,NORM_2,&resnorm));
        resnorm = PetscMax(resnorm,resl);
        PetscCall((*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx));
        if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
      } else if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
    }
    if (lock) {
      lock = PETSC_FALSE;
      skip = PETSC_TRUE;
      PetscCall(NEPDeflationNERecoverEigenvectors(defctx,u,w,lambda));
      PetscCall(NEPDeflationNELocking(defctx,u,w,lambda));
      PetscCall(NEPDeflationNESetRefine(defctx,PETSC_FALSE));
      PetscCall(BVInsertVec(nep->V,nep->nconv,u));
      PetscCall(BVInsertVec(nep->W,nep->nconv,w));
      PetscCall(VecConjugate(w));
      nep->nconv = nep->nconv + 1;
    }
    PetscCall((*nep->stopping)(nep,nep->its,nep->max_it,nep->nconv,nep->nev,&nep->reason,nep->stoppingctx));
    if (!skip || nep->reason>0) PetscCall(NEPMonitor(nep,nep->its,nep->nconv,nep->eigr,nep->eigi,nep->errest,(nep->reason>0)?nep->nconv:nep->nconv+1));

    if (nep->reason == NEP_CONVERGED_ITERATING) {
      if (!skip) {
        /* evaluate T(lambda) and T'(lambda) */
        PetscCall(NEPDeflationNEComputeFunction(nep,mF,lambda));
        PetscCall(NEPDeflationNEComputeFunction(nep,mJ,lambda));
        PetscCall(EPSSetInitialSpace(ctx->eps,1,&u));
        PetscCall(EPSSetInitialSpace(ctx->epsts,1,&w));

        /* compute new eigenvalue correction mu and eigenvector approximation u */
        PetscCall(EPSSolve(ctx->eps));
        PetscCall(EPSSolve(ctx->epsts));
        PetscCall(EPSGetConverged(ctx->eps,&nconv));
        PetscCall(EPSGetConverged(ctx->epsts,&nconv2));
        if (!nconv||!nconv2) {
          PetscCall(PetscInfo(nep,"iter=%" PetscInt_FMT ", inner iteration failed, stopping solve\n",nep->its));
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }
        PetscCall(EPSGetEigenpair(ctx->eps,0,&mu,&im,u,NULL));
        for (i=0;i<nconv2;i++) {
          PetscCall(EPSGetEigenpair(ctx->epsts,i,&mu2,&im2,w,NULL));
          if (SlepcAbsEigenvalue(mu-mu2,im-im2)/SlepcAbsEigenvalue(mu,im)<nep->tol*1000) break;
        }
        if (i==nconv2) {
          PetscCall(PetscInfo(nep,"iter=%" PetscInt_FMT ", inner iteration failed, stopping solve\n",nep->its));
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }

        mu = 1.0/mu;
        PetscCheck(PetscAbsScalar(im)<PETSC_MACHINE_EPSILON,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Complex eigenvalue approximation - not implemented in real scalars");
      } else {
        nep->its--;  /* do not count this as a full iteration */
        /* use second eigenpair computed in previous iteration */
        PetscCall(EPSGetConverged(ctx->eps,&nconv));
        if (nconv>=2 && nconv2>=2) {
          PetscCall(EPSGetEigenpair(ctx->eps,1,&mu,&im,u,NULL));
          PetscCall(EPSGetEigenpair(ctx->epsts,1,&mu2,&im2,w,NULL));
          mu = 1.0/mu;
        } else {
          PetscCall(BVSetRandomColumn(nep->V,nep->nconv));
          PetscCall(BVSetRandomColumn(nep->W,nep->nconv));
          PetscCall(BVCopyVec(nep->V,nep->nconv,u));
          PetscCall(BVCopyVec(nep->W,nep->nconv,w));
          mu = lambda-sigma;
        }
        skip = PETSC_FALSE;
      }
      /* correct eigenvalue */
      lambda = lambda - mu;
    }
  }
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&w));
  PetscCall(VecDestroy(&r));
  PetscCall(MatDestroy(&mF));
  PetscCall(MatDestroy(&mJ));
  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(NEPDeflationNEDestroy(defctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}
