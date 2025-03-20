/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "rii"

   Method: Residual inverse iteration

   Algorithm:

       Simple residual inverse iteration with varying shift.

   References:

       [1] A. Neumaier, "Residual inverse iteration for the nonlinear
           eigenvalue problem", SIAM J. Numer. Anal. 22(5):914-923, 1985.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <../src/nep/impls/nepdefl.h>

typedef struct {
  PetscInt  max_inner_it;     /* maximum number of Newton iterations */
  PetscInt  lag;              /* interval to rebuild preconditioner */
  PetscBool cctol;            /* constant correction tolerance */
  PetscBool herm;             /* whether the Hermitian version of the scalar equation must be used */
  PetscReal deftol;           /* tolerance for the deflation (threshold) */
  KSP       ksp;              /* linear solver object */
} NEP_RII;

static PetscErrorCode NEPSetUp_RII(NEP nep)
{
  PetscFunctionBegin;
  if (nep->ncv!=PETSC_DETERMINE) PetscCall(PetscInfo(nep,"Setting ncv = nev, ignoring user-provided value\n"));
  nep->ncv = nep->nev;
  if (nep->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(nep,"Setting mpd = nev, ignoring user-provided value\n"));
  nep->mpd = nep->nev;
  if (nep->max_it==PETSC_DETERMINE) nep->max_it = PetscMax(5000,2*nep->n/nep->ncv);
  if (!nep->which) nep->which = NEP_TARGET_MAGNITUDE;
  PetscCheck(nep->which==NEP_TARGET_MAGNITUDE,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver supports only target magnitude eigenvalues");
  NEPCheckUnsupported(nep,NEP_FEATURE_REGION | NEP_FEATURE_TWOSIDED);
  PetscCall(NEPAllocateSolution(nep,0));
  PetscCall(NEPSetWorkVecs(nep,2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSolve_RII(NEP nep)
{
  NEP_RII            *ctx = (NEP_RII*)nep->data;
  Mat                T,Tp,H,A;
  Vec                uu,u,r,delta,t;
  PetscScalar        lambda,lambda2,sigma,a1,a2,corr;
  PetscReal          nrm,resnorm=1.0,ktol=0.1,perr,rtol;
  PetscBool          skip=PETSC_FALSE,lock=PETSC_FALSE;
  PetscInt           inner_its,its=0;
  NEP_EXT_OP         extop=NULL;
  KSPConvergedReason kspreason;

  PetscFunctionBegin;
  /* get initial approximation of eigenvalue and eigenvector */
  PetscCall(NEPGetDefaultShift(nep,&sigma));
  lambda = sigma;
  if (!nep->nini) {
    PetscCall(BVSetRandomColumn(nep->V,0));
    PetscCall(BVNormColumn(nep->V,0,NORM_2,&nrm));
    PetscCall(BVScaleColumn(nep->V,0,1.0/nrm));
  }
  if (!ctx->ksp) PetscCall(NEPRIIGetKSP(nep,&ctx->ksp));
  PetscCall(NEPDeflationInitialize(nep,nep->V,ctx->ksp,PETSC_FALSE,nep->nev,&extop));
  PetscCall(NEPDeflationCreateVec(extop,&u));
  PetscCall(VecDuplicate(u,&r));
  PetscCall(VecDuplicate(u,&delta));
  PetscCall(BVGetColumn(nep->V,0,&uu));
  PetscCall(NEPDeflationCopyToExtendedVec(extop,uu,NULL,u,PETSC_FALSE));
  PetscCall(BVRestoreColumn(nep->V,0,&uu));

  /* prepare linear solver */
  PetscCall(NEPDeflationSolveSetUp(extop,sigma));
  PetscCall(KSPGetTolerances(ctx->ksp,&rtol,NULL,NULL,NULL));

  PetscCall(VecCopy(u,r));
  PetscCall(NEPDeflationFunctionSolve(extop,r,u));
  PetscCall(VecNorm(u,NORM_2,&nrm));
  PetscCall(VecScale(u,1.0/nrm));

  /* Restart loop */
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    its++;

    /* Use Newton's method to compute nonlinear Rayleigh functional. Current eigenvalue
       estimate as starting value. */
    inner_its=0;
    lambda2 = lambda;
    do {
      PetscCall(NEPDeflationComputeFunction(extop,lambda,&T));
      PetscCall(MatMult(T,u,r));
      if (!ctx->herm) {
        PetscCall(NEPDeflationFunctionSolve(extop,r,delta));
        PetscCall(KSPGetConvergedReason(ctx->ksp,&kspreason));
        if (kspreason<0) PetscCall(PetscInfo(nep,"iter=%" PetscInt_FMT ", linear solve failed\n",nep->its));
        t = delta;
      } else t = r;
      PetscCall(VecDot(t,u,&a1));
      PetscCall(NEPDeflationComputeJacobian(extop,lambda,&Tp));
      PetscCall(MatMult(Tp,u,r));
      if (!ctx->herm) {
        PetscCall(NEPDeflationFunctionSolve(extop,r,delta));
        PetscCall(KSPGetConvergedReason(ctx->ksp,&kspreason));
        if (kspreason<0) PetscCall(PetscInfo(nep,"iter=%" PetscInt_FMT ", linear solve failed\n",nep->its));
        t = delta;
      } else t = r;
      PetscCall(VecDot(t,u,&a2));
      corr = a1/a2;
      lambda = lambda - corr;
      inner_its++;
    } while (PetscAbsScalar(corr)/PetscAbsScalar(lambda)>PETSC_SQRT_MACHINE_EPSILON && inner_its<ctx->max_inner_it);

    /* form residual,  r = T(lambda)*u */
    PetscCall(NEPDeflationComputeFunction(extop,lambda,&T));
    PetscCall(MatMult(T,u,r));

    /* convergence test */
    perr = nep->errest[nep->nconv];
    PetscCall(VecNorm(r,NORM_2,&resnorm));
    PetscCall((*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx));
    nep->eigr[nep->nconv] = lambda;
    if (its>1 && (nep->errest[nep->nconv]<=nep->tol || nep->errest[nep->nconv]<=ctx->deftol)) {
      if (nep->errest[nep->nconv]<=ctx->deftol && !extop->ref && nep->nconv) {
        PetscCall(NEPDeflationExtractEigenpair(extop,nep->nconv,u,lambda,nep->ds));
        PetscCall(NEPDeflationSetRefine(extop,PETSC_TRUE));
        PetscCall(MatMult(T,u,r));
        PetscCall(VecNorm(r,NORM_2,&resnorm));
        PetscCall((*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx));
        if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
      } else if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
    }
    if (lock) {
      PetscCall(NEPDeflationSetRefine(extop,PETSC_FALSE));
      nep->nconv = nep->nconv + 1;
      PetscCall(NEPDeflationLocking(extop,u,lambda));
      nep->its += its;
      skip = PETSC_TRUE;
      lock = PETSC_FALSE;
    }
    PetscCall((*nep->stopping)(nep,nep->its+its,nep->max_it,nep->nconv,nep->nev,&nep->reason,nep->stoppingctx));
    if (!skip || nep->reason>0) PetscCall(NEPMonitor(nep,nep->its+its,nep->nconv,nep->eigr,nep->eigi,nep->errest,(nep->reason>0)?nep->nconv:nep->nconv+1));

    if (nep->reason == NEP_CONVERGED_ITERATING) {
      if (!skip) {
        /* update preconditioner and set adaptive tolerance */
        if (ctx->lag && !(its%ctx->lag) && its>=2*ctx->lag && perr && nep->errest[nep->nconv]>.5*perr) PetscCall(NEPDeflationSolveSetUp(extop,lambda2));
        if (!ctx->cctol) {
          ktol = PetscMax(ktol/2.0,rtol);
          PetscCall(KSPSetTolerances(ctx->ksp,ktol,PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
        }

        /* eigenvector correction: delta = T(sigma)\r */
        PetscCall(NEPDeflationFunctionSolve(extop,r,delta));
        PetscCall(KSPGetConvergedReason(ctx->ksp,&kspreason));
        if (kspreason<0) {
          PetscCall(PetscInfo(nep,"iter=%" PetscInt_FMT ", linear solve failed, stopping solve\n",nep->its));
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }

        /* update eigenvector: u = u - delta */
        PetscCall(VecAXPY(u,-1.0,delta));

        /* normalize eigenvector */
        PetscCall(VecNormalize(u,NULL));
      } else {
        its = -1;
        PetscCall(NEPDeflationSetRandomVec(extop,u));
        PetscCall(NEPDeflationSolveSetUp(extop,sigma));
        PetscCall(VecCopy(u,r));
        PetscCall(NEPDeflationFunctionSolve(extop,r,u));
        PetscCall(VecNorm(u,NORM_2,&nrm));
        PetscCall(VecScale(u,1.0/nrm));
        lambda = sigma;
        skip = PETSC_FALSE;
      }
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
  PetscCall(VecDestroy(&delta));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSetFromOptions_RII(NEP nep,PetscOptionItems PetscOptionsObject)
{
  NEP_RII        *ctx = (NEP_RII*)nep->data;
  PetscBool      flg;
  PetscInt       i;
  PetscReal      r;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"NEP RII Options");

    i = 0;
    PetscCall(PetscOptionsInt("-nep_rii_max_it","Maximum number of Newton iterations for updating Rayleigh functional","NEPRIISetMaximumIterations",ctx->max_inner_it,&i,&flg));
    if (flg) PetscCall(NEPRIISetMaximumIterations(nep,i));

    PetscCall(PetscOptionsBool("-nep_rii_const_correction_tol","Constant correction tolerance for the linear solver","NEPRIISetConstCorrectionTol",ctx->cctol,&ctx->cctol,NULL));

    PetscCall(PetscOptionsBool("-nep_rii_hermitian","Use Hermitian version of the scalar nonlinear equation","NEPRIISetHermitian",ctx->herm,&ctx->herm,NULL));

    i = 0;
    PetscCall(PetscOptionsInt("-nep_rii_lag_preconditioner","Interval to rebuild preconditioner","NEPRIISetLagPreconditioner",ctx->lag,&i,&flg));
    if (flg) PetscCall(NEPRIISetLagPreconditioner(nep,i));

    r = 0.0;
    PetscCall(PetscOptionsReal("-nep_rii_deflation_threshold","Tolerance used as a threshold for including deflated eigenpairs","NEPRIISetDeflationThreshold",ctx->deftol,&r,&flg));
    if (flg) PetscCall(NEPRIISetDeflationThreshold(nep,r));

  PetscOptionsHeadEnd();

  if (!ctx->ksp) PetscCall(NEPRIIGetKSP(nep,&ctx->ksp));
  PetscCall(KSPSetFromOptions(ctx->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIISetMaximumIterations_RII(NEP nep,PetscInt its)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  if (its==PETSC_DEFAULT || its==PETSC_DECIDE) ctx->max_inner_it = 10;
  else {
    PetscCheck(its>0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Number of iterations must be >0");
    ctx->max_inner_it = its;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIISetMaximumIterations - Sets the maximum number of inner iterations to be
   used in the RII solver. These are the Newton iterations related to the computation
   of the nonlinear Rayleigh functional.

   Logically Collective

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  its - maximum inner iterations

   Level: advanced

.seealso: NEPRIIGetMaximumIterations()
@*/
PetscErrorCode NEPRIISetMaximumIterations(NEP nep,PetscInt its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,its,2);
  PetscTryMethod(nep,"NEPRIISetMaximumIterations_C",(NEP,PetscInt),(nep,its));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIIGetMaximumIterations_RII(NEP nep,PetscInt *its)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *its = ctx->max_inner_it;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIIGetMaximumIterations - Gets the maximum number of inner iterations of RII.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  its - maximum inner iterations

   Level: advanced

.seealso: NEPRIISetMaximumIterations()
@*/
PetscErrorCode NEPRIIGetMaximumIterations(NEP nep,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(its,2);
  PetscUseMethod(nep,"NEPRIIGetMaximumIterations_C",(NEP,PetscInt*),(nep,its));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIISetLagPreconditioner_RII(NEP nep,PetscInt lag)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  PetscCheck(lag>=0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Lag must be non-negative");
  ctx->lag = lag;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIISetLagPreconditioner - Determines when the preconditioner is rebuilt in the
   nonlinear solve.

   Logically Collective

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  lag - 0 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is
          computed within the nonlinear iteration, 2 means every second time
          the Jacobian is built, etc.

   Options Database Keys:
.  -nep_rii_lag_preconditioner <lag> - the lag value

   Notes:
   The default is 1.
   The preconditioner is ALWAYS built in the first iteration of a nonlinear solve.

   Level: intermediate

.seealso: NEPRIIGetLagPreconditioner()
@*/
PetscErrorCode NEPRIISetLagPreconditioner(NEP nep,PetscInt lag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,lag,2);
  PetscTryMethod(nep,"NEPRIISetLagPreconditioner_C",(NEP,PetscInt),(nep,lag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIIGetLagPreconditioner_RII(NEP nep,PetscInt *lag)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *lag = ctx->lag;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIIGetLagPreconditioner - Indicates how often the preconditioner is rebuilt.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  lag - the lag parameter

   Level: intermediate

.seealso: NEPRIISetLagPreconditioner()
@*/
PetscErrorCode NEPRIIGetLagPreconditioner(NEP nep,PetscInt *lag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(lag,2);
  PetscUseMethod(nep,"NEPRIIGetLagPreconditioner_C",(NEP,PetscInt*),(nep,lag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIISetConstCorrectionTol_RII(NEP nep,PetscBool cct)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ctx->cctol = cct;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIISetConstCorrectionTol - Sets a flag to keep the tolerance used
   in the linear solver constant.

   Logically Collective

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  cct - a boolean value

   Options Database Keys:
.  -nep_rii_const_correction_tol <bool> - set the boolean flag

   Notes:
   By default, an exponentially decreasing tolerance is set in the KSP used
   within the nonlinear iteration, so that each Newton iteration requests
   better accuracy than the previous one. The constant correction tolerance
   flag stops this behaviour.

   Level: intermediate

.seealso: NEPRIIGetConstCorrectionTol()
@*/
PetscErrorCode NEPRIISetConstCorrectionTol(NEP nep,PetscBool cct)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,cct,2);
  PetscTryMethod(nep,"NEPRIISetConstCorrectionTol_C",(NEP,PetscBool),(nep,cct));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIIGetConstCorrectionTol_RII(NEP nep,PetscBool *cct)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *cct = ctx->cctol;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIIGetConstCorrectionTol - Returns the constant tolerance flag.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  cct - the value of the constant tolerance flag

   Level: intermediate

.seealso: NEPRIISetConstCorrectionTol()
@*/
PetscErrorCode NEPRIIGetConstCorrectionTol(NEP nep,PetscBool *cct)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(cct,2);
  PetscUseMethod(nep,"NEPRIIGetConstCorrectionTol_C",(NEP,PetscBool*),(nep,cct));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIISetHermitian_RII(NEP nep,PetscBool herm)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ctx->herm = herm;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIISetHermitian - Sets a flag to indicate if the Hermitian version of the
   scalar nonlinear equation must be used by the solver.

   Logically Collective

   Input Parameters:
+  nep  - nonlinear eigenvalue solver
-  herm - a boolean value

   Options Database Keys:
.  -nep_rii_hermitian <bool> - set the boolean flag

   Notes:
   By default, the scalar nonlinear equation x'*inv(T(sigma))*T(z)*x=0 is solved
   at each step of the nonlinear iteration. When this flag is set the simpler
   form x'*T(z)*x=0 is used, which is supposed to be valid only for Hermitian
   problems.

   Level: intermediate

.seealso: NEPRIIGetHermitian()
@*/
PetscErrorCode NEPRIISetHermitian(NEP nep,PetscBool herm)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,herm,2);
  PetscTryMethod(nep,"NEPRIISetHermitian_C",(NEP,PetscBool),(nep,herm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIIGetHermitian_RII(NEP nep,PetscBool *herm)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *herm = ctx->herm;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIIGetHermitian - Returns the flag about using the Hermitian version of
   the scalar nonlinear equation.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  herm - the value of the hermitian flag

   Level: intermediate

.seealso: NEPRIISetHermitian()
@*/
PetscErrorCode NEPRIIGetHermitian(NEP nep,PetscBool *herm)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(herm,2);
  PetscUseMethod(nep,"NEPRIIGetHermitian_C",(NEP,PetscBool*),(nep,herm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIISetDeflationThreshold_RII(NEP nep,PetscReal deftol)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ctx->deftol = deftol;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIISetDeflationThreshold - Sets the threshold value used to switch between
   deflated and non-deflated iteration.

   Logically Collective

   Input Parameters:
+  nep    - nonlinear eigenvalue solver
-  deftol - the threshold value

   Options Database Keys:
.  -nep_rii_deflation_threshold <deftol> - set the threshold

   Notes:
   Normally, the solver iterates on the extended problem in order to deflate
   previously converged eigenpairs. If this threshold is set to a nonzero value,
   then once the residual error is below this threshold the solver will
   continue the iteration without deflation. The intention is to be able to
   improve the current eigenpair further, despite having previous eigenpairs
   with somewhat bad precision.

   Level: advanced

.seealso: NEPRIIGetDeflationThreshold()
@*/
PetscErrorCode NEPRIISetDeflationThreshold(NEP nep,PetscReal deftol)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,deftol,2);
  PetscTryMethod(nep,"NEPRIISetDeflationThreshold_C",(NEP,PetscReal),(nep,deftol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIIGetDeflationThreshold_RII(NEP nep,PetscReal *deftol)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *deftol = ctx->deftol;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIIGetDeflationThreshold - Returns the threshold value that controls deflation.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  deftol - the threshold

   Level: advanced

.seealso: NEPRIISetDeflationThreshold()
@*/
PetscErrorCode NEPRIIGetDeflationThreshold(NEP nep,PetscReal *deftol)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(deftol,2);
  PetscUseMethod(nep,"NEPRIIGetDeflationThreshold_C",(NEP,PetscReal*),(nep,deftol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIISetKSP_RII(NEP nep,KSP ksp)
{
  NEP_RII        *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)ksp));
  PetscCall(KSPDestroy(&ctx->ksp));
  ctx->ksp = ksp;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIISetKSP - Associate a linear solver object (KSP) to the nonlinear
   eigenvalue solver.

   Collective

   Input Parameters:
+  nep - eigenvalue solver
-  ksp - the linear solver object

   Level: advanced

.seealso: NEPRIIGetKSP()
@*/
PetscErrorCode NEPRIISetKSP(NEP nep,KSP ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  PetscCheckSameComm(nep,1,ksp,2);
  PetscTryMethod(nep,"NEPRIISetKSP_C",(NEP,KSP),(nep,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPRIIGetKSP_RII(NEP nep,KSP *ksp)
{
  NEP_RII        *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    PetscCall(KSPCreate(PetscObjectComm((PetscObject)nep),&ctx->ksp));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->ksp,(PetscObject)nep,1));
    PetscCall(KSPSetOptionsPrefix(ctx->ksp,((PetscObject)nep)->prefix));
    PetscCall(KSPAppendOptionsPrefix(ctx->ksp,"nep_rii_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->ksp,((PetscObject)nep)->options));
    PetscCall(KSPSetErrorIfNotConverged(ctx->ksp,PETSC_TRUE));
    PetscCall(KSPSetTolerances(ctx->ksp,SlepcDefaultTol(nep->tol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
  }
  *ksp = ctx->ksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRIIGetKSP - Retrieve the linear solver object (KSP) associated with
   the nonlinear eigenvalue solver.

   Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  ksp - the linear solver object

   Level: advanced

.seealso: NEPRIISetKSP()
@*/
PetscErrorCode NEPRIIGetKSP(NEP nep,KSP *ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(ksp,2);
  PetscUseMethod(nep,"NEPRIIGetKSP_C",(NEP,KSP*),(nep,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPView_RII(NEP nep,PetscViewer viewer)
{
  NEP_RII        *ctx = (NEP_RII*)nep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum number of inner iterations: %" PetscInt_FMT "\n",ctx->max_inner_it));
    if (ctx->cctol) PetscCall(PetscViewerASCIIPrintf(viewer,"  using a constant tolerance for the linear solver\n"));
    if (ctx->herm) PetscCall(PetscViewerASCIIPrintf(viewer,"  using the Hermitian version of the scalar nonlinear equation\n"));
    if (ctx->lag) PetscCall(PetscViewerASCIIPrintf(viewer,"  updating the preconditioner every %" PetscInt_FMT " iterations\n",ctx->lag));
    if (ctx->deftol) PetscCall(PetscViewerASCIIPrintf(viewer,"  deflation threshold: %g\n",(double)ctx->deftol));
    if (!ctx->ksp) PetscCall(NEPRIIGetKSP(nep,&ctx->ksp));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(KSPView(ctx->ksp,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPReset_RII(NEP nep)
{
  NEP_RII        *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  PetscCall(KSPReset(ctx->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDestroy_RII(NEP nep)
{
  NEP_RII        *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  PetscCall(KSPDestroy(&ctx->ksp));
  PetscCall(PetscFree(nep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetMaximumIterations_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetMaximumIterations_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetLagPreconditioner_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetLagPreconditioner_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetConstCorrectionTol_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetConstCorrectionTol_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetHermitian_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetHermitian_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetDeflationThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetDeflationThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetKSP_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetKSP_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_RII(NEP nep)
{
  NEP_RII        *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  nep->data = (void*)ctx;
  ctx->max_inner_it = 10;
  ctx->lag          = 1;
  ctx->cctol        = PETSC_FALSE;
  ctx->herm         = PETSC_FALSE;
  ctx->deftol       = 0.0;

  nep->useds = PETSC_TRUE;

  nep->ops->solve          = NEPSolve_RII;
  nep->ops->setup          = NEPSetUp_RII;
  nep->ops->setfromoptions = NEPSetFromOptions_RII;
  nep->ops->reset          = NEPReset_RII;
  nep->ops->destroy        = NEPDestroy_RII;
  nep->ops->view           = NEPView_RII;
  nep->ops->computevectors = NEPComputeVectors_Schur;

  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetMaximumIterations_C",NEPRIISetMaximumIterations_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetMaximumIterations_C",NEPRIIGetMaximumIterations_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetLagPreconditioner_C",NEPRIISetLagPreconditioner_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetLagPreconditioner_C",NEPRIIGetLagPreconditioner_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetConstCorrectionTol_C",NEPRIISetConstCorrectionTol_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetConstCorrectionTol_C",NEPRIIGetConstCorrectionTol_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetHermitian_C",NEPRIISetHermitian_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetHermitian_C",NEPRIIGetHermitian_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetDeflationThreshold_C",NEPRIISetDeflationThreshold_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetDeflationThreshold_C",NEPRIIGetDeflationThreshold_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetKSP_C",NEPRIISetKSP_RII));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetKSP_C",NEPRIIGetKSP_RII));
  PetscFunctionReturn(PETSC_SUCCESS);
}
