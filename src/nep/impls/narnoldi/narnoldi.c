/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "narnoldi"

   Method: Nonlinear Arnoldi

   Algorithm:

       Arnoldi for nonlinear eigenproblems.

   References:

       [1] H. Voss, "An Arnoldi method for nonlinear eigenvalue problems",
           BIT 44:387-401, 2004.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <../src/nep/impls/nepdefl.h>

typedef struct {
  PetscInt lag;             /* interval to rebuild preconditioner */
  KSP      ksp;             /* linear solver object */
} NEP_NARNOLDI;

static PetscErrorCode NEPSetUp_NArnoldi(NEP nep)
{
  PetscFunctionBegin;
  PetscCall(NEPSetDimensions_Default(nep,nep->nev,&nep->ncv,&nep->mpd));
  PetscCheck(nep->ncv<=nep->nev+nep->mpd,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (nep->max_it==PETSC_DETERMINE) nep->max_it = nep->nev*nep->ncv;
  if (!nep->which) nep->which = NEP_TARGET_MAGNITUDE;
  PetscCheck(nep->which==NEP_TARGET_MAGNITUDE,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver supports only target magnitude eigenvalues");
  NEPCheckUnsupported(nep,NEP_FEATURE_CALLBACK | NEP_FEATURE_REGION | NEP_FEATURE_TWOSIDED);
  PetscCall(NEPAllocateSolution(nep,0));
  PetscCall(NEPSetWorkVecs(nep,3));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSolve_NArnoldi(NEP nep)
{
  NEP_NARNOLDI       *ctx = (NEP_NARNOLDI*)nep->data;
  Mat                T,H,A;
  Vec                f,r,u,uu;
  PetscScalar        *X,lambda=0.0,lambda2=0.0,*eigr,sigma;
  PetscReal          beta,resnorm=0.0,nrm,perr=0.0;
  PetscInt           n;
  PetscBool          breakdown,skip=PETSC_FALSE;
  BV                 Vext;
  DS                 ds;
  NEP_EXT_OP         extop=NULL;
  SlepcSC            sc;
  KSPConvergedReason kspreason;

  PetscFunctionBegin;
  /* get initial space and shift */
  PetscCall(NEPGetDefaultShift(nep,&sigma));
  if (!nep->nini) {
    PetscCall(BVSetRandomColumn(nep->V,0));
    PetscCall(BVNormColumn(nep->V,0,NORM_2,&nrm));
    PetscCall(BVScaleColumn(nep->V,0,1.0/nrm));
    n = 1;
  } else n = nep->nini;

  if (!ctx->ksp) PetscCall(NEPNArnoldiGetKSP(nep,&ctx->ksp));
  PetscCall(NEPDeflationInitialize(nep,nep->V,ctx->ksp,PETSC_FALSE,nep->nev,&extop));
  PetscCall(NEPDeflationCreateBV(extop,nep->ncv,&Vext));

  /* prepare linear solver */
  PetscCall(NEPDeflationSolveSetUp(extop,sigma));

  PetscCall(BVGetColumn(Vext,0,&f));
  PetscCall(VecDuplicate(f,&r));
  PetscCall(VecDuplicate(f,&u));
  PetscCall(BVGetColumn(nep->V,0,&uu));
  PetscCall(NEPDeflationCopyToExtendedVec(extop,uu,NULL,f,PETSC_FALSE));
  PetscCall(BVRestoreColumn(nep->V,0,&uu));
  PetscCall(VecCopy(f,r));
  PetscCall(NEPDeflationFunctionSolve(extop,r,f));
  PetscCall(VecNorm(f,NORM_2,&nrm));
  PetscCall(VecScale(f,1.0/nrm));
  PetscCall(BVRestoreColumn(Vext,0,&f));

  PetscCall(DSCreate(PetscObjectComm((PetscObject)nep),&ds));
  PetscCall(DSSetType(ds,DSNEP));
  PetscCall(DSNEPSetFN(ds,nep->nt,nep->f));
  PetscCall(DSAllocate(ds,nep->ncv));
  PetscCall(DSGetSlepcSC(ds,&sc));
  sc->comparison    = nep->sc->comparison;
  sc->comparisonctx = nep->sc->comparisonctx;
  PetscCall(DSSetFromOptions(ds));

  /* build projected matrices for initial space */
  PetscCall(DSSetDimensions(ds,n,0,0));
  PetscCall(NEPDeflationProjectOperator(extop,Vext,ds,0,n));

  PetscCall(PetscMalloc1(nep->ncv,&eigr));

  /* Restart loop */
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;

    /* solve projected problem */
    PetscCall(DSSetDimensions(ds,n,0,0));
    PetscCall(DSSetState(ds,DS_STATE_RAW));
    PetscCall(DSSolve(ds,eigr,NULL));
    PetscCall(DSSynchronize(ds,eigr,NULL));
    if (nep->its>1) lambda2 = lambda;
    lambda = eigr[0];
    nep->eigr[nep->nconv] = lambda;

    /* compute Ritz vector, x = V*s */
    PetscCall(DSGetArray(ds,DS_MAT_X,&X));
    PetscCall(BVSetActiveColumns(Vext,0,n));
    PetscCall(BVMultVec(Vext,1.0,0.0,u,X));
    PetscCall(DSRestoreArray(ds,DS_MAT_X,&X));

    /* compute the residual, r = T(lambda)*x */
    PetscCall(NEPDeflationComputeFunction(extop,lambda,&T));
    PetscCall(MatMult(T,u,r));

    /* convergence test */
    PetscCall(VecNorm(r,NORM_2,&resnorm));
    if (nep->its>1) perr=nep->errest[nep->nconv];
    PetscCall((*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx));
    if (nep->errest[nep->nconv]<=nep->tol) {
      nep->nconv = nep->nconv + 1;
      PetscCall(NEPDeflationLocking(extop,u,lambda));
      skip = PETSC_TRUE;
    }
    PetscCall((*nep->stopping)(nep,nep->its,nep->max_it,nep->nconv,nep->nev,&nep->reason,nep->stoppingctx));
    if (!skip || nep->reason>0) PetscCall(NEPMonitor(nep,nep->its,nep->nconv,nep->eigr,nep->eigi,nep->errest,(nep->reason>0)?nep->nconv:nep->nconv+1));

    if (nep->reason == NEP_CONVERGED_ITERATING) {
      if (!skip) {
        if (n>=nep->ncv) {
          nep->reason = NEP_DIVERGED_SUBSPACE_EXHAUSTED;
          break;
        }
        if (ctx->lag && !(nep->its%ctx->lag) && nep->its>=2*ctx->lag && perr && nep->errest[nep->nconv]>.5*perr) PetscCall(NEPDeflationSolveSetUp(extop,lambda2));

        /* continuation vector: f = T(sigma)\r */
        PetscCall(BVGetColumn(Vext,n,&f));
        PetscCall(NEPDeflationFunctionSolve(extop,r,f));
        PetscCall(BVRestoreColumn(Vext,n,&f));
        PetscCall(KSPGetConvergedReason(ctx->ksp,&kspreason));
        if (kspreason<0) {
          PetscCall(PetscInfo(nep,"iter=%" PetscInt_FMT ", linear solve failed, stopping solve\n",nep->its));
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }

        /* orthonormalize */
        PetscCall(BVOrthonormalizeColumn(Vext,n,PETSC_FALSE,&beta,&breakdown));
        if (breakdown || beta==0.0) {
          PetscCall(PetscInfo(nep,"iter=%" PetscInt_FMT ", orthogonalization failed, stopping solve\n",nep->its));
          nep->reason = NEP_DIVERGED_BREAKDOWN;
          break;
        }

        /* update projected matrices */
        PetscCall(DSSetDimensions(ds,n+1,0,0));
        PetscCall(NEPDeflationProjectOperator(extop,Vext,ds,n,n+1));
        n++;
      } else {
        nep->its--;  /* do not count this as a full iteration */
        PetscCall(BVGetColumn(Vext,0,&f));
        PetscCall(NEPDeflationSetRandomVec(extop,f));
        PetscCall(NEPDeflationSolveSetUp(extop,sigma));
        PetscCall(VecCopy(f,r));
        PetscCall(NEPDeflationFunctionSolve(extop,r,f));
        PetscCall(VecNorm(f,NORM_2,&nrm));
        PetscCall(VecScale(f,1.0/nrm));
        PetscCall(BVRestoreColumn(Vext,0,&f));
        n = 1;
        PetscCall(DSSetDimensions(ds,n,0,0));
        PetscCall(NEPDeflationProjectOperator(extop,Vext,ds,n-1,n));
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
  PetscCall(BVDestroy(&Vext));
  PetscCall(DSDestroy(&ds));
  PetscCall(PetscFree(eigr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNArnoldiSetLagPreconditioner_NArnoldi(NEP nep,PetscInt lag)
{
  NEP_NARNOLDI *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  PetscCheck(lag>=0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Lag must be non-negative");
  ctx->lag = lag;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNArnoldiSetLagPreconditioner - Determines when the preconditioner is rebuilt in the
   nonlinear solve.

   Logically Collective

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  lag - 0 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is
          computed within the nonlinear iteration, 2 means every second time
          the Jacobian is built, etc.

   Options Database Keys:
.  -nep_narnoldi_lag_preconditioner <lag> - the lag value

   Notes:
   The default is 1.
   The preconditioner is ALWAYS built in the first iteration of a nonlinear solve.

   Level: intermediate

.seealso: NEPNArnoldiGetLagPreconditioner()
@*/
PetscErrorCode NEPNArnoldiSetLagPreconditioner(NEP nep,PetscInt lag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,lag,2);
  PetscTryMethod(nep,"NEPNArnoldiSetLagPreconditioner_C",(NEP,PetscInt),(nep,lag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNArnoldiGetLagPreconditioner_NArnoldi(NEP nep,PetscInt *lag)
{
  NEP_NARNOLDI *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  *lag = ctx->lag;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNArnoldiGetLagPreconditioner - Indicates how often the preconditioner is rebuilt.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  lag - the lag parameter

   Level: intermediate

.seealso: NEPNArnoldiSetLagPreconditioner()
@*/
PetscErrorCode NEPNArnoldiGetLagPreconditioner(NEP nep,PetscInt *lag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(lag,2);
  PetscUseMethod(nep,"NEPNArnoldiGetLagPreconditioner_C",(NEP,PetscInt*),(nep,lag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSetFromOptions_NArnoldi(NEP nep,PetscOptionItems PetscOptionsObject)
{
  PetscInt       i;
  PetscBool      flg;
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"NEP N-Arnoldi Options");
    i = 0;
    PetscCall(PetscOptionsInt("-nep_narnoldi_lag_preconditioner","Interval to rebuild preconditioner","NEPNArnoldiSetLagPreconditioner",ctx->lag,&i,&flg));
    if (flg) PetscCall(NEPNArnoldiSetLagPreconditioner(nep,i));

  PetscOptionsHeadEnd();

  if (!ctx->ksp) PetscCall(NEPNArnoldiGetKSP(nep,&ctx->ksp));
  PetscCall(KSPSetFromOptions(ctx->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNArnoldiSetKSP_NArnoldi(NEP nep,KSP ksp)
{
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)ksp));
  PetscCall(KSPDestroy(&ctx->ksp));
  ctx->ksp = ksp;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNArnoldiSetKSP - Associate a linear solver object (KSP) to the nonlinear
   eigenvalue solver.

   Collective

   Input Parameters:
+  nep - eigenvalue solver
-  ksp - the linear solver object

   Level: advanced

.seealso: NEPNArnoldiGetKSP()
@*/
PetscErrorCode NEPNArnoldiSetKSP(NEP nep,KSP ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  PetscCheckSameComm(nep,1,ksp,2);
  PetscTryMethod(nep,"NEPNArnoldiSetKSP_C",(NEP,KSP),(nep,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNArnoldiGetKSP_NArnoldi(NEP nep,KSP *ksp)
{
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    PetscCall(KSPCreate(PetscObjectComm((PetscObject)nep),&ctx->ksp));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->ksp,(PetscObject)nep,1));
    PetscCall(KSPSetOptionsPrefix(ctx->ksp,((PetscObject)nep)->prefix));
    PetscCall(KSPAppendOptionsPrefix(ctx->ksp,"nep_narnoldi_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->ksp,((PetscObject)nep)->options));
    PetscCall(KSPSetErrorIfNotConverged(ctx->ksp,PETSC_TRUE));
    PetscCall(KSPSetTolerances(ctx->ksp,SlepcDefaultTol(nep->tol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
  }
  *ksp = ctx->ksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNArnoldiGetKSP - Retrieve the linear solver object (KSP) associated with
   the nonlinear eigenvalue solver.

   Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  ksp - the linear solver object

   Level: advanced

.seealso: NEPNArnoldiSetKSP()
@*/
PetscErrorCode NEPNArnoldiGetKSP(NEP nep,KSP *ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(ksp,2);
  PetscUseMethod(nep,"NEPNArnoldiGetKSP_C",(NEP,KSP*),(nep,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPView_NArnoldi(NEP nep,PetscViewer viewer)
{
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (ctx->lag) PetscCall(PetscViewerASCIIPrintf(viewer,"  updating the preconditioner every %" PetscInt_FMT " iterations\n",ctx->lag));
    if (!ctx->ksp) PetscCall(NEPNArnoldiGetKSP(nep,&ctx->ksp));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(KSPView(ctx->ksp,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPReset_NArnoldi(NEP nep)
{
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  PetscCall(KSPReset(ctx->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDestroy_NArnoldi(NEP nep)
{
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  PetscCall(KSPDestroy(&ctx->ksp));
  PetscCall(PetscFree(nep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiSetLagPreconditioner_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiGetLagPreconditioner_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiSetKSP_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiGetKSP_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_NArnoldi(NEP nep)
{
  NEP_NARNOLDI   *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  nep->data = (void*)ctx;
  ctx->lag  = 1;

  nep->useds = PETSC_TRUE;

  nep->ops->solve          = NEPSolve_NArnoldi;
  nep->ops->setup          = NEPSetUp_NArnoldi;
  nep->ops->setfromoptions = NEPSetFromOptions_NArnoldi;
  nep->ops->reset          = NEPReset_NArnoldi;
  nep->ops->destroy        = NEPDestroy_NArnoldi;
  nep->ops->view           = NEPView_NArnoldi;
  nep->ops->computevectors = NEPComputeVectors_Schur;

  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiSetLagPreconditioner_C",NEPNArnoldiSetLagPreconditioner_NArnoldi));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiGetLagPreconditioner_C",NEPNArnoldiGetLagPreconditioner_NArnoldi));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiSetKSP_C",NEPNArnoldiSetKSP_NArnoldi));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiGetKSP_C",NEPNArnoldiGetKSP_NArnoldi));
  PetscFunctionReturn(PETSC_SUCCESS);
}
