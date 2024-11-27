/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "arnoldi"

   Method: Explicitly Restarted Arnoldi

   Algorithm:

       Arnoldi method with explicit restart and deflation.

   References:

       [1] "Arnoldi Methods in SLEPc", SLEPc Technical Report STR-4,
           available at https://slepc.upv.es.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/

typedef struct {
  PetscBool delayed;
} EPS_ARNOLDI;

static PetscErrorCode EPSSetUp_Arnoldi(EPS eps)
{
  PetscFunctionBegin;
  EPSCheckDefinite(eps);
  EPSCheckNotStructured(eps);
  PetscCall(EPSSetDimensions_Default(eps,&eps->nev,&eps->ncv,&eps->mpd));
  PetscCheck(eps->ncv<=eps->nev+eps->mpd,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  PetscCheck(eps->which!=EPS_ALL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_THRESHOLD | EPS_FEATURE_TWOSIDED);

  PetscCall(EPSAllocateSolution(eps,1));
  PetscCall(EPS_SetInnerProduct(eps));
  PetscCall(DSSetType(eps->ds,DSNHEP));
  if (eps->extraction==EPS_REFINED || eps->extraction==EPS_REFINED_HARMONIC) PetscCall(DSSetRefined(eps->ds,PETSC_TRUE));
  PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
  PetscCall(DSAllocate(eps->ds,eps->ncv+1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_Arnoldi(EPS eps)
{
  PetscInt           k,nv,ld;
  Mat                U,Op,H;
  PetscScalar        *Harray;
  PetscReal          beta,gamma=1.0;
  PetscBool          breakdown,harmonic,refined;
  BVOrthogRefineType orthog_ref;
  EPS_ARNOLDI        *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(DSGetRefined(eps->ds,&refined));
  harmonic = (eps->extraction==EPS_HARMONIC || eps->extraction==EPS_REFINED_HARMONIC)?PETSC_TRUE:PETSC_FALSE;
  PetscCall(BVGetOrthogonalization(eps->V,NULL,&orthog_ref,NULL,NULL));

  /* Get the starting Arnoldi vector */
  PetscCall(EPSGetStartVector(eps,0,NULL));

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Arnoldi factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,0));
    if (!arnoldi->delayed) {
      PetscCall(STGetOperator(eps->st,&Op));
      PetscCall(DSGetMat(eps->ds,DS_MAT_A,&H));
      PetscCall(BVMatArnoldi(eps->V,Op,H,eps->nconv,&nv,&beta,&breakdown));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&H));
      PetscCall(STRestoreOperator(eps->st,&Op));
    } else if (orthog_ref == BV_ORTHOG_REFINE_NEVER) {
      PetscCall(DSGetArray(eps->ds,DS_MAT_A,&Harray));
      PetscCall(EPSDelayedArnoldi1(eps,Harray,ld,eps->nconv,&nv,&beta,&breakdown));
      PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&Harray));
    } else {
      PetscCall(DSGetArray(eps->ds,DS_MAT_A,&Harray));
      PetscCall(EPSDelayedArnoldi(eps,Harray,ld,eps->nconv,&nv,&beta,&breakdown));
      PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&Harray));
    }
    PetscCall(DSSetState(eps->ds,DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));

    /* Compute translation of Krylov decomposition if harmonic extraction used */
    if (harmonic) PetscCall(DSTranslateHarmonic(eps->ds,eps->target,beta,PETSC_FALSE,NULL,&gamma));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(eps->ds));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

    /* Check convergence */
    PetscCall(EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta,0.0,gamma,&k));
    if (refined) {
      PetscCall(DSGetMat(eps->ds,DS_MAT_X,&U));
      PetscCall(BVMultInPlace(eps->V,U,eps->nconv,k+1));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&U));
      PetscCall(BVOrthonormalizeColumn(eps->V,k,PETSC_FALSE,NULL,NULL));
    } else {
      PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&U));
      PetscCall(BVMultInPlace(eps->V,U,eps->nconv,PetscMin(k+1,nv)));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&U));
    }
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));
    if (eps->reason == EPS_CONVERGED_ITERATING && breakdown) {
      PetscCall(PetscInfo(eps,"Breakdown in Arnoldi method (it=%" PetscInt_FMT " norm=%g)\n",eps->its,(double)beta));
      PetscCall(EPSGetStartVector(eps,k,&breakdown));
      if (breakdown) {
        eps->reason = EPS_DIVERGED_BREAKDOWN;
        PetscCall(PetscInfo(eps,"Unable to generate more start vectors\n"));
      }
    }
    eps->nconv = k;
    PetscCall(EPSMonitor(eps,eps->its,eps->nconv,eps->eigr,eps->eigi,eps->errest,nv));
  }
  PetscCall(DSTruncate(eps->ds,eps->nconv,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_Arnoldi(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  PetscBool      set,val;
  EPS_ARNOLDI    *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS Arnoldi Options");

    PetscCall(PetscOptionsBool("-eps_arnoldi_delayed","Use delayed reorthogonalization","EPSArnoldiSetDelayed",arnoldi->delayed,&val,&set));
    if (set) PetscCall(EPSArnoldiSetDelayed(eps,val));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSArnoldiSetDelayed_Arnoldi(EPS eps,PetscBool delayed)
{
  EPS_ARNOLDI *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  arnoldi->delayed = delayed;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSArnoldiSetDelayed - Activates or deactivates delayed reorthogonalization
   in the Arnoldi iteration.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  delayed - boolean flag

   Options Database Key:
.  -eps_arnoldi_delayed - Activates delayed reorthogonalization in Arnoldi

   Note:
   Delayed reorthogonalization is an aggressive optimization for the Arnoldi
   eigensolver than may provide better scalability, but sometimes makes the
   solver converge less than the default algorithm.

   Level: advanced

.seealso: EPSArnoldiGetDelayed()
@*/
PetscErrorCode EPSArnoldiSetDelayed(EPS eps,PetscBool delayed)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,delayed,2);
  PetscTryMethod(eps,"EPSArnoldiSetDelayed_C",(EPS,PetscBool),(eps,delayed));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSArnoldiGetDelayed_Arnoldi(EPS eps,PetscBool *delayed)
{
  EPS_ARNOLDI *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  *delayed = arnoldi->delayed;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSArnoldiGetDelayed - Gets the type of reorthogonalization used during the Arnoldi
   iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  delayed - boolean flag indicating if delayed reorthogonalization has been enabled

   Level: advanced

.seealso: EPSArnoldiSetDelayed()
@*/
PetscErrorCode EPSArnoldiGetDelayed(EPS eps,PetscBool *delayed)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(delayed,2);
  PetscUseMethod(eps,"EPSArnoldiGetDelayed_C",(EPS,PetscBool*),(eps,delayed));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_Arnoldi(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSArnoldiSetDelayed_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSArnoldiGetDelayed_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_Arnoldi(EPS eps,PetscViewer viewer)
{
  PetscBool      isascii;
  EPS_ARNOLDI    *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii && arnoldi->delayed) PetscCall(PetscViewerASCIIPrintf(viewer,"  using delayed reorthogonalization\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_Arnoldi(EPS eps)
{
  EPS_ARNOLDI    *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->useds = PETSC_TRUE;

  eps->ops->solve          = EPSSolve_Arnoldi;
  eps->ops->setup          = EPSSetUp_Arnoldi;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->setfromoptions = EPSSetFromOptions_Arnoldi;
  eps->ops->destroy        = EPSDestroy_Arnoldi;
  eps->ops->view           = EPSView_Arnoldi;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_Schur;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSArnoldiSetDelayed_C",EPSArnoldiSetDelayed_Arnoldi));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSArnoldiGetDelayed_C",EPSArnoldiGetDelayed_Arnoldi));
  PetscFunctionReturn(PETSC_SUCCESS);
}
