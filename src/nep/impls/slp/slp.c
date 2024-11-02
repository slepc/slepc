/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "slp"

   Method: Successive linear problems

   Algorithm:

       Newton-type iteration based on first order Taylor approximation.

   References:

       [1] A. Ruhe, "Algorithms for the nonlinear eigenvalue problem", SIAM J.
           Numer. Anal. 10(4):674-689, 1973.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <../src/nep/impls/nepdefl.h>
#include "slp.h"

typedef struct {
  NEP_EXT_OP extop;
  Vec        w;
} NEP_SLP_MATSHELL;

static PetscErrorCode NEPSetUp_SLP(NEP nep)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  PetscBool      flg;
  ST             st;

  PetscFunctionBegin;
  if (nep->ncv!=PETSC_DETERMINE) PetscCall(PetscInfo(nep,"Setting ncv = nev, ignoring user-provided value\n"));
  nep->ncv = nep->nev;
  if (nep->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(nep,"Setting mpd = nev, ignoring user-provided value\n"));
  nep->mpd = nep->nev;
  PetscCheck(nep->ncv<=nep->nev+nep->mpd,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (nep->max_it==PETSC_DETERMINE) nep->max_it = PetscMax(5000,2*nep->n/nep->ncv);
  if (!nep->which) nep->which = NEP_TARGET_MAGNITUDE;
  PetscCheck(nep->which==NEP_TARGET_MAGNITUDE,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver supports only target magnitude eigenvalues");
  NEPCheckUnsupported(nep,NEP_FEATURE_REGION);

  if (!ctx->eps) PetscCall(NEPSLPGetEPS(nep,&ctx->eps));
  PetscCall(EPSGetST(ctx->eps,&st));
  PetscCall(PetscObjectTypeCompareAny((PetscObject)st,&flg,STSINVERT,STCAYLEY,""));
  PetscCheck(!flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"SLP does not support spectral transformation");
  PetscCall(EPSSetDimensions(ctx->eps,1,PETSC_DECIDE,PETSC_DECIDE));
  PetscCall(EPSSetWhichEigenpairs(ctx->eps,EPS_LARGEST_MAGNITUDE));
  PetscCall(EPSSetTolerances(ctx->eps,SlepcDefaultTol(nep->tol)/10.0,nep->max_it));
  if (nep->tol==(PetscReal)PETSC_DETERMINE) nep->tol = SLEPC_DEFAULT_TOL;
  if (ctx->deftol==(PetscReal)PETSC_DETERMINE) ctx->deftol = nep->tol;

  if (nep->twosided) {
    nep->ops->solve = NEPSolve_SLP_Twosided;
    nep->ops->computevectors = NULL;
    if (!ctx->epsts) PetscCall(NEPSLPGetEPSLeft(nep,&ctx->epsts));
    PetscCall(EPSGetST(ctx->epsts,&st));
    PetscCall(PetscObjectTypeCompareAny((PetscObject)st,&flg,STSINVERT,STCAYLEY,""));
    PetscCheck(!flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"SLP does not support spectral transformation");
    PetscCall(EPSSetDimensions(ctx->epsts,1,PETSC_DECIDE,PETSC_DECIDE));
    PetscCall(EPSSetWhichEigenpairs(ctx->epsts,EPS_LARGEST_MAGNITUDE));
    PetscCall(EPSSetTolerances(ctx->epsts,SlepcDefaultTol(nep->tol)/10.0,nep->max_it));
  } else {
    nep->ops->solve = NEPSolve_SLP;
    nep->ops->computevectors = NEPComputeVectors_Schur;
  }
  PetscCall(NEPAllocateSolution(nep,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_SLP(Mat M,Vec x,Vec y)
{
  NEP_SLP_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(MatMult(ctx->extop->MJ,x,ctx->w));
  PetscCall(NEPDeflationFunctionSolve(ctx->extop,ctx->w,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_SLP(Mat M)
{
  NEP_SLP_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(VecDestroy(&ctx->w));
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
static PetscErrorCode MatCreateVecs_SLP(Mat M,Vec *left,Vec *right)
{
  NEP_SLP_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  if (right) PetscCall(VecDuplicate(ctx->w,right));
  if (left) PetscCall(VecDuplicate(ctx->w,left));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode NEPSLPSetUpLinearEP(NEP nep,NEP_EXT_OP extop,PetscScalar lambda,Vec u,PetscBool ini)
{
  NEP_SLP          *slpctx = (NEP_SLP*)nep->data;
  Mat              Mshell;
  PetscInt         nloc,mloc;
  NEP_SLP_MATSHELL *shellctx;

  PetscFunctionBegin;
  if (ini) {
    /* Create mat shell */
    PetscCall(PetscNew(&shellctx));
    shellctx->extop = extop;
    PetscCall(NEPDeflationCreateVec(extop,&shellctx->w));
    PetscCall(MatGetLocalSize(nep->function,&mloc,&nloc));
    nloc += extop->szd; mloc += extop->szd;
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)nep),nloc,mloc,PETSC_DETERMINE,PETSC_DETERMINE,shellctx,&Mshell));
    PetscCall(MatShellSetOperation(Mshell,MATOP_MULT,(void(*)(void))MatMult_SLP));
    PetscCall(MatShellSetOperation(Mshell,MATOP_DESTROY,(void(*)(void))MatDestroy_SLP));
#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
    PetscCall(MatShellSetOperation(Mshell,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_SLP));
#endif
    PetscCall(EPSSetOperators(slpctx->eps,Mshell,NULL));
    PetscCall(MatDestroy(&Mshell));
  }
  PetscCall(NEPDeflationSolveSetUp(extop,lambda));
  PetscCall(NEPDeflationComputeJacobian(extop,lambda,NULL));
  PetscCall(EPSSetInitialSpace(slpctx->eps,1,&u));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPSolve_SLP(NEP nep)
{
  NEP_SLP           *ctx = (NEP_SLP*)nep->data;
  Mat               F,H,A;
  Vec               uu,u,r;
  PetscScalar       sigma,lambda,mu,im;
  PetscReal         resnorm;
  PetscInt          nconv;
  PetscBool         skip=PETSC_FALSE,lock=PETSC_FALSE;
  NEP_EXT_OP        extop=NULL;    /* Extended operator for deflation */

  PetscFunctionBegin;
  /* get initial approximation of eigenvalue and eigenvector */
  PetscCall(NEPGetDefaultShift(nep,&sigma));
  if (!nep->nini) PetscCall(BVSetRandomColumn(nep->V,0));
  lambda = sigma;
  if (!ctx->ksp) PetscCall(NEPSLPGetKSP(nep,&ctx->ksp));
  PetscCall(NEPDeflationInitialize(nep,nep->V,ctx->ksp,PETSC_TRUE,nep->nev,&extop));
  PetscCall(NEPDeflationCreateVec(extop,&u));
  PetscCall(VecDuplicate(u,&r));
  PetscCall(BVGetColumn(nep->V,0,&uu));
  PetscCall(NEPDeflationCopyToExtendedVec(extop,uu,NULL,u,PETSC_FALSE));
  PetscCall(BVRestoreColumn(nep->V,0,&uu));

  /* Restart loop */
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;

    /* form residual,  r = T(lambda)*u (used in convergence test only) */
    PetscCall(NEPDeflationComputeFunction(extop,lambda,&F));
    PetscCall(MatMult(F,u,r));

    /* convergence test */
    PetscCall(VecNorm(r,NORM_2,&resnorm));
    PetscCall((*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx));
    nep->eigr[nep->nconv] = lambda;
    if (nep->errest[nep->nconv]<=nep->tol || nep->errest[nep->nconv]<=ctx->deftol) {
      if (nep->errest[nep->nconv]<=ctx->deftol && !extop->ref && nep->nconv) {
        PetscCall(NEPDeflationExtractEigenpair(extop,nep->nconv,u,lambda,nep->ds));
        PetscCall(NEPDeflationSetRefine(extop,PETSC_TRUE));
        PetscCall(MatMult(F,u,r));
        PetscCall(VecNorm(r,NORM_2,&resnorm));
        PetscCall((*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx));
        if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
      } else if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
    }

    if (lock) {
      PetscCall(NEPDeflationSetRefine(extop,PETSC_FALSE));
      nep->nconv = nep->nconv + 1;
      skip = PETSC_TRUE;
      lock = PETSC_FALSE;
      PetscCall(NEPDeflationLocking(extop,u,lambda));
    }
    PetscCall((*nep->stopping)(nep,nep->its,nep->max_it,nep->nconv,nep->nev,&nep->reason,nep->stoppingctx));
    if (!skip || nep->reason>0) PetscCall(NEPMonitor(nep,nep->its,nep->nconv,nep->eigr,nep->eigi,nep->errest,(nep->reason>0)?nep->nconv:nep->nconv+1));

    if (nep->reason == NEP_CONVERGED_ITERATING) {
      if (!skip) {
        /* evaluate T(lambda) and T'(lambda) */
        PetscCall(NEPSLPSetUpLinearEP(nep,extop,lambda,u,nep->its==1?PETSC_TRUE:PETSC_FALSE));
        /* compute new eigenvalue correction mu and eigenvector approximation u */
        PetscCall(EPSSolve(ctx->eps));
        PetscCall(EPSGetConverged(ctx->eps,&nconv));
        if (!nconv) {
          PetscCall(PetscInfo(nep,"iter=%" PetscInt_FMT ", inner iteration failed, stopping solve\n",nep->its));
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }
        PetscCall(EPSGetEigenpair(ctx->eps,0,&mu,&im,u,NULL));
        mu = 1.0/mu;
        PetscCheck(PetscAbsScalar(im)<PETSC_MACHINE_EPSILON,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Complex eigenvalue approximation - not implemented in real scalars");
      } else {
        nep->its--;  /* do not count this as a full iteration */
        /* use second eigenpair computed in previous iteration */
        PetscCall(EPSGetConverged(ctx->eps,&nconv));
        if (nconv>=2) {
          PetscCall(EPSGetEigenpair(ctx->eps,1,&mu,&im,u,NULL));
          mu = 1.0/mu;
        } else {
          PetscCall(NEPDeflationSetRandomVec(extop,u));
          mu = lambda-sigma;
        }
        skip = PETSC_FALSE;
      }
      /* correct eigenvalue */
      lambda = lambda - mu;
    }
  }
  PetscCall(NEPDeflationGetInvariantPair(extop,NULL,&H));
  PetscCall(DSSetType(nep->ds,DSNHEP));
  PetscCall(DSAllocate(nep->ds,PetscMax(nep->nconv,1)));
  PetscCall(DSSetDimensions(nep->ds,nep->nconv,0,nep->nconv));
  PetscCall(DSGetMat(nep->ds,DS_MAT_A,&A));
  PetscCall(MatCopy(H,A,SAME_NONZERO_PATTERN));
  PetscCall(DSRestoreMat(nep->ds,DS_MAT_A,&A));
  PetscCall(MatDestroy(&H));
  PetscCall(DSSolve(nep->ds,nep->eigr,nep->eigi));
  PetscCall(NEPDeflationReset(extop));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&r));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSetFromOptions_SLP(NEP nep,PetscOptionItems PetscOptionsObject)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  PetscBool      flg;
  PetscReal      r;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"NEP SLP Options");

    r = 0.0;
    PetscCall(PetscOptionsReal("-nep_slp_deflation_threshold","Tolerance used as a threshold for including deflated eigenpairs","NEPSLPSetDeflationThreshold",ctx->deftol,&r,&flg));
    if (flg) PetscCall(NEPSLPSetDeflationThreshold(nep,r));

  PetscOptionsHeadEnd();

  if (!ctx->eps) PetscCall(NEPSLPGetEPS(nep,&ctx->eps));
  PetscCall(EPSSetFromOptions(ctx->eps));
  if (nep->twosided) {
    if (!ctx->epsts) PetscCall(NEPSLPGetEPSLeft(nep,&ctx->epsts));
    PetscCall(EPSSetFromOptions(ctx->epsts));
  }
  if (!ctx->ksp) PetscCall(NEPSLPGetKSP(nep,&ctx->ksp));
  PetscCall(KSPSetFromOptions(ctx->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSLPSetDeflationThreshold_SLP(NEP nep,PetscReal deftol)
{
  NEP_SLP *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  if (deftol == (PetscReal)PETSC_DEFAULT || deftol == (PetscReal)PETSC_DECIDE) {
    ctx->deftol = PETSC_DETERMINE;
    nep->state  = NEP_STATE_INITIAL;
  } else {
    PetscCheck(deftol>0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of deftol. Must be > 0");
    ctx->deftol = deftol;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSLPSetDeflationThreshold - Sets the threshold value used to switch between
   deflated and non-deflated iteration.

   Logically Collective

   Input Parameters:
+  nep    - nonlinear eigenvalue solver
-  deftol - the threshold value

   Options Database Keys:
.  -nep_slp_deflation_threshold <deftol> - set the threshold

   Notes:
   Normally, the solver iterates on the extended problem in order to deflate
   previously converged eigenpairs. If this threshold is set to a nonzero value,
   then once the residual error is below this threshold the solver will
   continue the iteration without deflation. The intention is to be able to
   improve the current eigenpair further, despite having previous eigenpairs
   with somewhat bad precision.

   Level: advanced

.seealso: NEPSLPGetDeflationThreshold()
@*/
PetscErrorCode NEPSLPSetDeflationThreshold(NEP nep,PetscReal deftol)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,deftol,2);
  PetscTryMethod(nep,"NEPSLPSetDeflationThreshold_C",(NEP,PetscReal),(nep,deftol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSLPGetDeflationThreshold_SLP(NEP nep,PetscReal *deftol)
{
  NEP_SLP *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  *deftol = ctx->deftol;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSLPGetDeflationThreshold - Returns the threshold value that controls deflation.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  deftol - the threshold

   Level: advanced

.seealso: NEPSLPSetDeflationThreshold()
@*/
PetscErrorCode NEPSLPGetDeflationThreshold(NEP nep,PetscReal *deftol)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(deftol,2);
  PetscUseMethod(nep,"NEPSLPGetDeflationThreshold_C",(NEP,PetscReal*),(nep,deftol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSLPSetEPS_SLP(NEP nep,EPS eps)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)eps));
  PetscCall(EPSDestroy(&ctx->eps));
  ctx->eps = eps;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSLPSetEPS - Associate a linear eigensolver object (EPS) to the
   nonlinear eigenvalue solver.

   Collective

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  eps - the eigensolver object

   Level: advanced

.seealso: NEPSLPGetEPS()
@*/
PetscErrorCode NEPSLPSetEPS(NEP nep,EPS eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(nep,1,eps,2);
  PetscTryMethod(nep,"NEPSLPSetEPS_C",(NEP,EPS),(nep,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSLPGetEPS_SLP(NEP nep,EPS *eps)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  if (!ctx->eps) {
    PetscCall(EPSCreate(PetscObjectComm((PetscObject)nep),&ctx->eps));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->eps,(PetscObject)nep,1));
    PetscCall(EPSSetOptionsPrefix(ctx->eps,((PetscObject)nep)->prefix));
    PetscCall(EPSAppendOptionsPrefix(ctx->eps,"nep_slp_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->eps,((PetscObject)nep)->options));
  }
  *eps = ctx->eps;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSLPGetEPS - Retrieve the linear eigensolver object (EPS) associated
   to the nonlinear eigenvalue solver.

   Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  eps - the eigensolver object

   Level: advanced

.seealso: NEPSLPSetEPS()
@*/
PetscErrorCode NEPSLPGetEPS(NEP nep,EPS *eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(eps,2);
  PetscUseMethod(nep,"NEPSLPGetEPS_C",(NEP,EPS*),(nep,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSLPSetEPSLeft_SLP(NEP nep,EPS eps)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)eps));
  PetscCall(EPSDestroy(&ctx->epsts));
  ctx->epsts = eps;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSLPSetEPSLeft - Associate a linear eigensolver object (EPS) to the
   nonlinear eigenvalue solver, used to compute left eigenvectors in the
   two-sided variant of SLP.

   Collective

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  eps - the eigensolver object

   Level: advanced

.seealso: NEPSLPGetEPSLeft(), NEPSetTwoSided()
@*/
PetscErrorCode NEPSLPSetEPSLeft(NEP nep,EPS eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(nep,1,eps,2);
  PetscTryMethod(nep,"NEPSLPSetEPSLeft_C",(NEP,EPS),(nep,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSLPGetEPSLeft_SLP(NEP nep,EPS *eps)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  if (!ctx->epsts) {
    PetscCall(EPSCreate(PetscObjectComm((PetscObject)nep),&ctx->epsts));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->epsts,(PetscObject)nep,1));
    PetscCall(EPSSetOptionsPrefix(ctx->epsts,((PetscObject)nep)->prefix));
    PetscCall(EPSAppendOptionsPrefix(ctx->epsts,"nep_slp_left_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->epsts,((PetscObject)nep)->options));
  }
  *eps = ctx->epsts;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSLPGetEPSLeft - Retrieve the linear eigensolver object (EPS) associated
   to the nonlinear eigenvalue solver, used to compute left eigenvectors in the
   two-sided variant of SLP.

   Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  eps - the eigensolver object

   Level: advanced

.seealso: NEPSLPSetEPSLeft(), NEPSetTwoSided()
@*/
PetscErrorCode NEPSLPGetEPSLeft(NEP nep,EPS *eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(eps,2);
  PetscUseMethod(nep,"NEPSLPGetEPSLeft_C",(NEP,EPS*),(nep,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSLPSetKSP_SLP(NEP nep,KSP ksp)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)ksp));
  PetscCall(KSPDestroy(&ctx->ksp));
  ctx->ksp   = ksp;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSLPSetKSP - Associate a linear solver object (KSP) to the nonlinear
   eigenvalue solver.

   Collective

   Input Parameters:
+  nep - eigenvalue solver
-  ksp - the linear solver object

   Level: advanced

.seealso: NEPSLPGetKSP()
@*/
PetscErrorCode NEPSLPSetKSP(NEP nep,KSP ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  PetscCheckSameComm(nep,1,ksp,2);
  PetscTryMethod(nep,"NEPSLPSetKSP_C",(NEP,KSP),(nep,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSLPGetKSP_SLP(NEP nep,KSP *ksp)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    PetscCall(KSPCreate(PetscObjectComm((PetscObject)nep),&ctx->ksp));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->ksp,(PetscObject)nep,1));
    PetscCall(KSPSetOptionsPrefix(ctx->ksp,((PetscObject)nep)->prefix));
    PetscCall(KSPAppendOptionsPrefix(ctx->ksp,"nep_slp_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->ksp,((PetscObject)nep)->options));
    PetscCall(KSPSetErrorIfNotConverged(ctx->ksp,PETSC_TRUE));
    PetscCall(KSPSetTolerances(ctx->ksp,SlepcDefaultTol(nep->tol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
  }
  *ksp = ctx->ksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSLPGetKSP - Retrieve the linear solver object (KSP) associated with
   the nonlinear eigenvalue solver.

   Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  ksp - the linear solver object

   Level: advanced

.seealso: NEPSLPSetKSP()
@*/
PetscErrorCode NEPSLPGetKSP(NEP nep,KSP *ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(ksp,2);
  PetscUseMethod(nep,"NEPSLPGetKSP_C",(NEP,KSP*),(nep,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPView_SLP(NEP nep,PetscViewer viewer)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (ctx->deftol) PetscCall(PetscViewerASCIIPrintf(viewer,"  deflation threshold: %g\n",(double)ctx->deftol));
    if (!ctx->eps) PetscCall(NEPSLPGetEPS(nep,&ctx->eps));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(EPSView(ctx->eps,viewer));
    if (nep->twosided) {
      if (!ctx->epsts) PetscCall(NEPSLPGetEPSLeft(nep,&ctx->epsts));
      PetscCall(EPSView(ctx->epsts,viewer));
    }
    if (!ctx->ksp) PetscCall(NEPSLPGetKSP(nep,&ctx->ksp));
    PetscCall(KSPView(ctx->ksp,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPReset_SLP(NEP nep)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  PetscCall(EPSReset(ctx->eps));
  if (nep->twosided) PetscCall(EPSReset(ctx->epsts));
  PetscCall(KSPReset(ctx->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDestroy_SLP(NEP nep)
{
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  PetscCall(KSPDestroy(&ctx->ksp));
  PetscCall(EPSDestroy(&ctx->eps));
  PetscCall(EPSDestroy(&ctx->epsts));
  PetscCall(PetscFree(nep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetDeflationThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetDeflationThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetEPSLeft_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetEPSLeft_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetKSP_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetKSP_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_SLP(NEP nep)
{
  NEP_SLP        *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  nep->data = (void*)ctx;

  nep->useds  = PETSC_TRUE;
  ctx->deftol = PETSC_DETERMINE;

  nep->ops->solve          = NEPSolve_SLP;
  nep->ops->setup          = NEPSetUp_SLP;
  nep->ops->setfromoptions = NEPSetFromOptions_SLP;
  nep->ops->reset          = NEPReset_SLP;
  nep->ops->destroy        = NEPDestroy_SLP;
  nep->ops->view           = NEPView_SLP;
  nep->ops->computevectors = NEPComputeVectors_Schur;

  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetDeflationThreshold_C",NEPSLPSetDeflationThreshold_SLP));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetDeflationThreshold_C",NEPSLPGetDeflationThreshold_SLP));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetEPS_C",NEPSLPSetEPS_SLP));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetEPS_C",NEPSLPGetEPS_SLP));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetEPSLeft_C",NEPSLPSetEPSLeft_SLP));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetEPSLeft_C",NEPSLPGetEPSLeft_SLP));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetKSP_C",NEPSLPSetKSP_SLP));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetKSP_C",NEPSLPGetKSP_SLP));
  PetscFunctionReturn(PETSC_SUCCESS);
}
