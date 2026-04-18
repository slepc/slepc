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

       [2] H.-X. Zhong, H. Xu, "Weighted Golub-Kahan-Lanczos bidiagonalization
           algorithms", Elec. Trans. Numer. Anal. 47, 2017.

       [3] F. Alvarruiz, B. Mellado-Pinto, J. E. Roman, "Restarted Lanczos methods
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
  PetscBool      alloc=PETSC_FALSE;

  PetscFunctionBegin;
  if (4*m > 100) {
    PetscCall(PetscMalloc1(4*m,&hwork));
    alloc = PETSC_TRUE;
  } else hwork = lhwork;

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
  if (alloc) PetscCall(PetscFree(hwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* K-Orthogonalize vector vh against first j vectors in V using beta coeffs */
static PetscErrorCode OrthogonalizeVector_Zhong_v(Vec vh,BV V,PetscInt j,PetscReal *beta,PetscInt k,PetscScalar *h)
{
  PetscInt  i,l;

  PetscFunctionBegin;
  /* Local orthogonalization */
  l = j==k?0:j-1;  /* 1st column to orthogonalize against */
  for (i=l;i<j;i++) h[i] = beta[i];
  /* vh = vh-V[:,l:j-1]*h[l:j-1] */
  PetscCall(BVSetActiveColumns(V,l,j));
  PetscCall(BVMultVec(V,-1.0,1.0,vh,h+l));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* M-Orthogonalize vector uh against first j vectors in U. Full orthog */
static PetscErrorCode Orthogonalize_Zhong_u(Vec uh,BV U,BV MU,PetscInt j,PetscScalar *h)
{
  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(U,0,j));
  PetscCall(BVSetActiveColumns(MU,0,j));
  /* h=MU'*uh */
  PetscCall(BVDotVec(MU,uh,h));
  /* uh=uh-U*h */
  PetscCall(BVMultVec(U,-1.0,1.0,uh,h));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* M-Orthogonalize vector uh against first j vectors in U. Local+full orthog */
static PetscErrorCode OrthogonalizeVector_Zhong_u(Vec uh,BV U,BV MU,PetscInt j,PetscReal *alpha,PetscInt k,PetscScalar *h)
{
  Vec u;

  PetscFunctionBegin;
  /* Local orthogonalization: uh = uh-U[:,j-1]*alpha[j-1] */
  PetscCall(BVGetColumn(U,j-1,&u));
  PetscCall(VecAXPY(uh,-alpha[j-1],u));
  PetscCall(BVRestoreColumn(U,j-1,&u));
  /* Full orthogonalization */
  PetscCall(Orthogonalize_Zhong_u(uh,U,MU,j,h));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLREPLanczos_Zhong(EPS eps,Mat K,Mat M,BV U,BV V,BV MU,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *min,PetscBool *breakdown)
{
  PetscInt       j,m = *min;
  Vec            uh,vh,x;
  PetscReal      beta0;
  PetscScalar    *hwork,lhwork[100],gamma;
  PetscBool      alloc=PETSC_FALSE;

  PetscFunctionBegin;
  if (m > 100) {
    PetscCall(PetscMalloc1(m,&hwork));
    alloc = PETSC_TRUE;
  } else hwork = lhwork;

  /* Normalize initial vector */
  if (k==0) {
    if (eps->nini==0) PetscCall(BVSetRandomColumn(U,0));
    PetscCall(BVGetColumn(U,0,&uh));
    PetscCall(BVGetColumn(MU,0,&vh));
    PetscCall(MatMult(M,uh,vh));
    PetscCall(VecDot(uh,vh,&gamma));
    beta0 = PetscSqrtReal(PetscRealPart(gamma));
    if (beta0==0.0) {
      if (breakdown) *breakdown = PETSC_TRUE;
      *min = 1; m = 0;
    } else {
      PetscCall(VecScale(uh,1.0/beta0));
      PetscCall(VecScale(vh,1.0/beta0));
    }
    PetscCall(BVRestoreColumn(U,0,&uh));
    PetscCall(BVRestoreColumn(MU,0,&vh));
  }

  for (j=k;j<m;j++) {
    /* Compute column j of V, then column j+1 of U and MU */
    PetscCall(BVGetColumn(U,j+1,&uh));
    PetscCall(BVGetColumn(V,j,&vh));
    PetscCall(BVGetColumn(MU,j,&x));
    PetscCall(VecCopy(x,vh));
    PetscCall(BVRestoreColumn(MU,j,&x));
    PetscCall(OrthogonalizeVector_Zhong_v(vh,V,j,beta,k,hwork));
    PetscCall(MatMult(K,vh,uh));
    PetscCall(VecDot(uh,vh,&gamma));
    alpha[j] = PetscSqrtReal(PetscRealPart(gamma));
    if (alpha[j]==0.0) {
      if (breakdown) *breakdown = PETSC_TRUE;
      *min = j+1; m = j;
      PetscCall(BVRestoreColumn(U,j+1,&uh));
    } else {
      PetscCall(VecScale(uh,1.0/alpha[j]));
      PetscCall(VecScale(vh,1.0/alpha[j]));
    }
    PetscCall(BVRestoreColumn(V,j,&vh));
    if (breakdown && *breakdown) continue;

    PetscCall(OrthogonalizeVector_Zhong_u(uh,U,MU,j+1,alpha,k,hwork));
    PetscCall(BVGetColumn(MU,j+1,&vh));
    PetscCall(MatMult(M,uh,vh));
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
    PetscCall(BVRestoreColumn(MU,j+1,&vh));
  }
  if (alloc) PetscCall(PetscFree(hwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSConvergence_Zhong - convergence check based on SVDKrylovConvergence().
   FIXME: Code dulication. This is a copy of EPSConvergence_Gruning
*/
static PetscErrorCode EPSConvergence_Zhong(EPS eps,PetscBool getall,PetscInt kini,PetscInt nits,PetscInt *kout)
{
  PetscInt       k,marker,ld;
  PetscReal      *alpha,*beta,resnorm;
  PetscBool      extra;

  PetscFunctionBegin;
  *kout = 0;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(DSGetExtraRow(eps->ds,&extra));
  PetscCheck(extra,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Only implemented for DS with extra row");
  marker = -1;
  if (eps->trackall) getall = PETSC_TRUE;
  PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&alpha));
  beta = alpha + ld;
  for (k=kini;k<kini+nits;k++) {
    resnorm = PetscAbsReal(beta[k]);
    PetscCall((*eps->converged)(eps,eps->eigr[k],eps->eigi[k],resnorm,&eps->errest[k],eps->convergedctx));
    if (marker==-1 && eps->errest[k] >= eps->tol) marker = k;
    if (marker!=-1 && !getall) break;
  }
  PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&alpha));
  if (marker!=-1) k = marker;
  *kout = k;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSUnreduceVectors(EPS eps,BV U,BV V)
{
  PetscInt k;
  Vec      u,v,w;

  PetscFunctionBegin;
  /* The approximate eigenvector is [u+v; u-v], where [u; v] is the reduced eigenvector */
  PetscCall(BVCreateVec(V,&w));
  for (k=0;k<eps->nconv;k++) {
    PetscCall(BVGetColumn(U,k,&u));
    PetscCall(BVGetColumn(V,k,&v));
    PetscCall(VecCopy(v,w));
    PetscCall(VecCopy(u,v));
    PetscCall(VecAXPY(u,1.0,w));
    PetscCall(VecAXPY(v,-1.0,w));
    PetscCall(BVRestoreColumn(U,k,&u));
    PetscCall(BVRestoreColumn(V,k,&v));
  }
  PetscCall(VecDestroy(&w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_LREP_Teng(EPS eps)
{
  Mat         H;
  Vec         v;
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
  if (!reduced) PetscCall(EPSUnreduceVectors(eps,V,U));
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&V,&U));
  /* Normalize eigenvectors */
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
  PetscCall(BVNormalize(eps->V,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_LREP_Zhong(EPS eps)
{
  Mat         H;
  BV          U,V;
  IS          is[2];
  PetscBool   reduced;

  PetscFunctionBegin;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(SlepcCheckMatLREPReduced(H,&reduced));
  /* Approx eigenvector for the reduced form is [u; v] */
  if (!reduced) {
    PetscCall(MatNestGetISs(H,is,NULL));
    PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));
    PetscCall(EPSUnreduceVectors(eps, U, V));
    PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
  }
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
  switch (ctx->lrep) {
    case EPS_KRYLOVSCHUR_LREP_TENG:
      eps->ops->solve = EPSSolve_KrylovSchur_LREP_Teng;
      eps->ops->computevectors = EPSComputeVectors_LREP_Teng;
      PetscCall(DSSetType(eps->ds,DSHEP));
      PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      break;
    case EPS_KRYLOVSCHUR_LREP_ZHONG:
      eps->ops->solve = EPSSolve_KrylovSchur_LREP_Zhong;
      eps->ops->computevectors = EPSComputeVectors_LREP_Zhong;
      PetscCall(DSSetType(eps->ds,DSSVD));
      PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      break;
    default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Unexpected error");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCreateReducedMats(Mat H,Mat *K,Mat *M)
{
  PetscInt          ma,na,Ma,Na;
  Mat               A,B;
  const PetscScalar scal[] = { 1.0, -1.0 };

  PetscFunctionBegin;
  PetscCall(MatNestGetSubMat(H,0,0,&A));
  PetscCall(MatNestGetSubMat(H,0,1,&B));
  PetscCall(MatGetSize(A,&Ma,&Na));
  PetscCall(MatGetLocalSize(A,&ma,&na));
  /* K = A-B */
  PetscCall(MatCreate(PetscObjectComm((PetscObject)A),K));
  PetscCall(MatSetSizes(*K,ma,na,Ma,Na));
  PetscCall(MatSetType(*K,MATCOMPOSITE));
  PetscCall(MatCompositeAddMat(*K,A));
  PetscCall(MatCompositeAddMat(*K,B));
  PetscCall(MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatCompositeSetScalings(*K,scal));
  /* M = A+B */
  PetscCall(MatCreate(PetscObjectComm((PetscObject)A),M));
  PetscCall(MatSetSizes(*M,ma,na,Ma,Na));
  PetscCall(MatSetType(*M,MATCOMPOSITE));
  PetscCall(MatCompositeAddMat(*M,A));
  PetscCall(MatCompositeAddMat(*M,B));
  PetscCall(MatAssemblyBegin(*M,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(*M,MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_LREP_Teng(EPS eps)
{
  EPS_KRYLOVSCHUR   *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt          i,k,l,ld,nv,nconv=0,nevsave;
  Mat               H,Q,K,M;
  BV                U,V;
  IS                is[2];
  PetscReal         *a,*b,beta;
  PetscBool         reduced,breakdown=PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));

  /* Extract matrix blocks */
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));
  PetscCall(SlepcCheckMatLREPReduced(H,&reduced));
  if (reduced) {
    PetscCall(MatNestGetSubMat(H,0,1,&K));
    PetscCall(MatNestGetSubMat(H,1,0,&M));
  } else PetscCall(EPSCreateReducedMats(H,&K,&M));

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

PetscErrorCode EPSSolve_KrylovSchur_LREP_Zhong(EPS eps)
{
  EPS_KRYLOVSCHUR   *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt          i,k,l,ld,nv,nconv=0,nevsave;
  Mat               H,Q,Z,K,M;
  BV                U,V,MU;
  IS                is[2];
  PetscReal         *a,*b,beta;
  PetscBool         reduced,breakdown=PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));

  /* Extract matrix blocks */
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));
  PetscCall(SlepcCheckMatLREPReduced(H,&reduced));
  if (reduced) {
    PetscCall(MatNestGetSubMat(H,0,1,&K));
    PetscCall(MatNestGetSubMat(H,1,0,&M));
  } else PetscCall(EPSCreateReducedMats(H,&K,&M));

  /* Get the split bases */
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));

  /* Create MU */
  PetscCall(BVDuplicate(U,&MU));

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
    PetscCall(EPSLREPLanczos_Zhong(eps,K,M,U,V,MU,a,b,eps->nconv+l,&nv,&breakdown));
    beta = b[nv-1];
    PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&a));
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSSVDSetDimensions(eps->ds,nv));
    PetscCall(DSSetState(eps->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(U,eps->nconv,nv));
    PetscCall(BVSetActiveColumns(V,eps->nconv,nv));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(eps->ds));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

    /* Check convergence */
    PetscCall(EPSConvergence_Zhong(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,&k));
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
       U(:,idx) = U*Q(:,idx),  MU(:,idx) = MU*Q(:,idx),  V(:,idx) = V*Z(:,idx),  */
    PetscCall(DSGetMat(eps->ds,DS_MAT_U,&Z));
    PetscCall(DSGetMat(eps->ds,DS_MAT_V,&Q));
    PetscCall(BVMultInPlace(U,Q,eps->nconv,k+l));
    PetscCall(BVMultInPlace(MU,Q,eps->nconv,k+l));
    PetscCall(BVMultInPlace(V,Z,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_U,&Z));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_V,&Q));

    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      PetscCall(BVCopyColumn(U,nv,k+l));
      PetscCall(BVCopyColumn(MU,nv,k+l));
      if (eps->stop==EPS_STOP_THRESHOLD && nv-k<5) {  /* reallocate */
        eps->ncv = eps->mpd+k;
        PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
        PetscCall(EPSReallocateSolution(eps,eps->ncv+1));
        PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));
        PetscCall(BVResize(MU,eps->ncv+1,PETSC_TRUE));
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
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
  PetscCall(BVDestroy(&MU));
  if (!reduced) {
    PetscCall(MatDestroy(&K));
    PetscCall(MatDestroy(&M));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
