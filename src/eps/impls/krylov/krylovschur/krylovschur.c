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

   Method: Krylov-Schur

   Algorithm:

       Single-vector Krylov-Schur method for non-symmetric problems,
       including harmonic extraction.

   References:

       [1] "Krylov-Schur Methods in SLEPc", SLEPc Technical Report STR-7,
           available at https://slepc.upv.es.

       [2] G.W. Stewart, "A Krylov-Schur Algorithm for Large Eigenproblems",
           SIAM J. Matrix Anal. App. 23(3):601-614, 2001.

       [3] "Practical Implementation of Harmonic Krylov-Schur", SLEPc Technical
            Report STR-9, available at https://slepc.upv.es.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include "krylovschur.h"

PetscErrorCode EPSGetArbitraryValues(EPS eps,PetscScalar *rr,PetscScalar *ri)
{
  PetscInt       i,newi,ld,n,l;
  Vec            xr=eps->work[0],xi=eps->work[1];
  PetscScalar    re,im,*Zr,*Zi,*X;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(DSGetDimensions(eps->ds,&n,&l,NULL,NULL));
  for (i=l;i<n;i++) {
    re = eps->eigr[i];
    im = eps->eigi[i];
    PetscCall(STBackTransform(eps->st,1,&re,&im));
    newi = i;
    PetscCall(DSVectors(eps->ds,DS_MAT_X,&newi,NULL));
    PetscCall(DSGetArray(eps->ds,DS_MAT_X,&X));
    Zr = X+i*ld;
    if (newi==i+1) Zi = X+newi*ld;
    else Zi = NULL;
    PetscCall(EPSComputeRitzVector(eps,Zr,Zi,eps->V,xr,xi));
    PetscCall(DSRestoreArray(eps->ds,DS_MAT_X,&X));
    PetscCall((*eps->arbitrary)(re,im,xr,xi,rr+i,ri+i,eps->arbitraryctx));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUp_KrylovSchur_Filter(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscBool       estimaterange=PETSC_TRUE;
  PetscReal       rleft,rright;
  Mat             A;

  PetscFunctionBegin;
  EPSCheckHermitianCondition(eps,PETSC_TRUE," with polynomial filter");
  EPSCheckStandardCondition(eps,PETSC_TRUE," with polynomial filter");
  PetscCheck(eps->intb<PETSC_MAX_REAL || eps->inta>PETSC_MIN_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"The defined computational interval should have at least one of their sides bounded");
  EPSCheckUnsupportedCondition(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_EXTRACTION | EPS_FEATURE_THRESHOLD,PETSC_TRUE," with polynomial filter");
  if (eps->tol==(PetscReal)PETSC_DETERMINE) eps->tol = SLEPC_DEFAULT_TOL*1e-2;  /* use tighter tolerance */
  PetscCall(STFilterSetInterval(eps->st,eps->inta,eps->intb));
  if (!ctx->estimatedrange) {
    PetscCall(STFilterGetRange(eps->st,&rleft,&rright));
    estimaterange = (!rleft && !rright)? PETSC_TRUE: PETSC_FALSE;
  }
  if (estimaterange) { /* user did not set a range */
    PetscCall(STGetMatrix(eps->st,0,&A));
    PetscCall(MatEstimateSpectralRange_EPS(A,&rleft,&rright));
    PetscCall(PetscInfo(eps,"Setting eigenvalue range to [%g,%g]\n",(double)rleft,(double)rright));
    PetscCall(STFilterSetRange(eps->st,rleft,rright));
    ctx->estimatedrange = PETSC_TRUE;
  }
  if (eps->ncv==PETSC_DETERMINE && eps->nev==0) eps->nev = 40;  /* user did not provide nev estimation */
  PetscCall(EPSSetDimensions_Default(eps,&eps->nev,&eps->ncv,&eps->mpd));
  PetscCheck(eps->ncv<=eps->nev+eps->mpd,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUp_KrylovSchur(EPS eps)
{
  PetscReal         eta;
  PetscBool         isfilt=PETSC_FALSE;
  BVOrthogType      otype;
  BVOrthogBlockType obtype;
  EPS_KRYLOVSCHUR   *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  enum { EPS_KS_DEFAULT,EPS_KS_SYMM,EPS_KS_SLICE,EPS_KS_FILTER,EPS_KS_INDEF,EPS_KS_TWOSIDED } variant;

  PetscFunctionBegin;
  if (eps->which==EPS_ALL) {  /* default values in case of spectrum slicing or polynomial filter  */
    PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt));
    if (isfilt) PetscCall(EPSSetUp_KrylovSchur_Filter(eps));
    else PetscCall(EPSSetUp_KrylovSchur_Slice(eps));
  } else if (eps->isstructured) {
    PetscCall(EPSSetUp_KrylovSchur_BSE(eps));
    PetscFunctionReturn(PETSC_SUCCESS);
  } else {
    PetscCall(EPSSetDimensions_Default(eps,&eps->nev,&eps->ncv,&eps->mpd));
    PetscCheck(eps->ncv<=eps->nev+eps->mpd,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
    if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv)*((eps->stop==EPS_STOP_THRESHOLD)?10:1);
    if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  }
  PetscCheck(ctx->lock || eps->mpd>=eps->ncv,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Should not use mpd parameter in non-locking variant");

  EPSCheckDefiniteCondition(eps,eps->arbitrary," with arbitrary selection of eigenpairs");

  PetscCheck(eps->extraction==EPS_RITZ || eps->extraction==EPS_HARMONIC,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");

  if (!ctx->keep) ctx->keep = 0.5;

  PetscCall(EPSAllocateSolution(eps,1));
  PetscCall(EPS_SetInnerProduct(eps));
  if (eps->arbitrary) PetscCall(EPSSetWorkVecs(eps,2));
  else if (eps->ishermitian && !eps->ispositive) PetscCall(EPSSetWorkVecs(eps,1));

  /* dispatch solve method */
  if (eps->ishermitian) {
    if (eps->which==EPS_ALL) {
      EPSCheckDefiniteCondition(eps,eps->which==EPS_ALL," with spectrum slicing");
      variant = isfilt? EPS_KS_FILTER: EPS_KS_SLICE;
    } else if (eps->isgeneralized && !eps->ispositive) {
      variant = EPS_KS_INDEF;
    } else {
      switch (eps->extraction) {
        case EPS_RITZ:     variant = EPS_KS_SYMM; break;
        case EPS_HARMONIC: variant = EPS_KS_DEFAULT; break;
        default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");
      }
    }
  } else if (eps->twosided) {
    variant = EPS_KS_TWOSIDED;
  } else {
    switch (eps->extraction) {
      case EPS_RITZ:     variant = EPS_KS_DEFAULT; break;
      case EPS_HARMONIC: variant = EPS_KS_DEFAULT; break;
      default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");
    }
  }
  switch (variant) {
    case EPS_KS_DEFAULT:
      eps->ops->solve = EPSSolve_KrylovSchur_Default;
      eps->ops->computevectors = EPSComputeVectors_Schur;
      PetscCall(DSSetType(eps->ds,DSNHEP));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      break;
    case EPS_KS_SYMM:
    case EPS_KS_FILTER:
      eps->ops->solve = EPSSolve_KrylovSchur_Default;
      eps->ops->computevectors = EPSComputeVectors_Hermitian;
      PetscCall(DSSetType(eps->ds,DSHEP));
      PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      break;
    case EPS_KS_SLICE:
      eps->ops->solve = EPSSolve_KrylovSchur_Slice;
      eps->ops->computevectors = EPSComputeVectors_Slice;
      break;
    case EPS_KS_INDEF:
      eps->ops->solve = EPSSolve_KrylovSchur_Indefinite;
      eps->ops->computevectors = EPSComputeVectors_Indefinite;
      PetscCall(DSSetType(eps->ds,DSGHIEP));
      PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      /* force reorthogonalization for pseudo-Lanczos */
      PetscCall(BVGetOrthogonalization(eps->V,&otype,NULL,&eta,&obtype));
      PetscCall(BVSetOrthogonalization(eps->V,otype,BV_ORTHOG_REFINE_ALWAYS,eta,obtype));
      break;
    case EPS_KS_TWOSIDED:
      eps->ops->solve = EPSSolve_KrylovSchur_TwoSided;
      eps->ops->computevectors = EPSComputeVectors_Schur;
      PetscCall(DSSetType(eps->ds,DSNHEPTS));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      break;
    default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Unexpected error");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUpSort_KrylovSchur(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  SlepcSC         sc;
  PetscBool       isfilt;

  PetscFunctionBegin;
  PetscCall(EPSSetUpSort_Default(eps));
  if (eps->which==EPS_ALL) {
    PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt));
    if (isfilt) {
      PetscCall(DSGetSlepcSC(eps->ds,&sc));
      sc->rg            = NULL;
      sc->comparison    = SlepcCompareLargestReal;
      sc->comparisonctx = NULL;
      sc->map           = NULL;
      sc->mapobj        = NULL;
    } else {
      if (!ctx->global && ctx->sr->numEigs>0) {
        PetscCall(DSGetSlepcSC(eps->ds,&sc));
        sc->rg            = NULL;
        sc->comparison    = SlepcCompareLargestMagnitude;
        sc->comparisonctx = NULL;
        sc->map           = NULL;
        sc->mapobj        = NULL;
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_Default(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,j,*pj,k,l,nv,ld,nconv;
  Mat             U,Op,H,T;
  PetscScalar     *g;
  PetscReal       beta,gamma=1.0;
  PetscBool       breakdown,harmonic,hermitian;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  harmonic = (eps->extraction==EPS_HARMONIC || eps->extraction==EPS_REFINED_HARMONIC)?PETSC_TRUE:PETSC_FALSE;
  hermitian = (eps->ishermitian && !harmonic)?PETSC_TRUE:PETSC_FALSE;
  if (harmonic) PetscCall(PetscMalloc1(ld,&g));
  if (eps->arbitrary) pj = &j;
  else pj = NULL;

  /* Get the starting Arnoldi vector */
  PetscCall(EPSGetStartVector(eps,0,NULL));
  l = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Arnoldi factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(STGetOperator(eps->st,&Op));
    if (hermitian) {
      PetscCall(DSGetMat(eps->ds,DS_MAT_T,&T));
      PetscCall(BVMatLanczos(eps->V,Op,T,eps->nconv+l,&nv,&beta,&breakdown));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_T,&T));
    } else {
      PetscCall(DSGetMat(eps->ds,DS_MAT_A,&H));
      PetscCall(BVMatArnoldi(eps->V,Op,H,eps->nconv+l,&nv,&beta,&breakdown));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&H));
    }
    PetscCall(STRestoreOperator(eps->st,&Op));
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSSetState(eps->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));

    /* Compute translation of Krylov decomposition if harmonic extraction used */
    if (PetscUnlikely(harmonic)) PetscCall(DSTranslateHarmonic(eps->ds,eps->target,beta,PETSC_FALSE,g,&gamma));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    if (PetscUnlikely(eps->arbitrary)) {
      PetscCall(EPSGetArbitraryValues(eps,eps->rr,eps->ri));
      j=1;
    }
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,eps->rr,eps->ri,pj));
    PetscCall(DSUpdateExtraRow(eps->ds));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

    /* Check convergence */
    PetscCall(EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta,0.0,gamma,&k));
    EPSSetCtxThreshold(eps,eps->eigr,eps->eigi,k);
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
      if (!hermitian) PetscCall(DSGetTruncateSize(eps->ds,k,nv,&l));
    }
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(eps,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      if (PetscUnlikely(breakdown || k==nv)) {
        /* Start a new Arnoldi factorization */
        PetscCall(PetscInfo(eps,"Breakdown in Krylov-Schur method (it=%" PetscInt_FMT " norm=%g)\n",eps->its,(double)beta));
        if (k<eps->nev) {
          PetscCall(EPSGetStartVector(eps,k,&breakdown));
          if (breakdown) {
            eps->reason = EPS_DIVERGED_BREAKDOWN;
            PetscCall(PetscInfo(eps,"Unable to generate more start vectors\n"));
          }
        }
      } else {
        /* Undo translation of Krylov decomposition */
        if (PetscUnlikely(harmonic)) {
          PetscCall(DSSetDimensions(eps->ds,nv,k,l));
          PetscCall(DSTranslateHarmonic(eps->ds,0.0,beta,PETSC_TRUE,g,&gamma));
          /* gamma u^ = u - U*g~ */
          PetscCall(BVSetActiveColumns(eps->V,0,nv));
          PetscCall(BVMultColumn(eps->V,-1.0,1.0,nv,g));
          PetscCall(BVScaleColumn(eps->V,nv,1.0/gamma));
          PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));
          PetscCall(DSSetDimensions(eps->ds,nv,k,nv));
        }
        /* Prepare the Rayleigh quotient for restart */
        PetscCall(DSTruncate(eps->ds,k+l,PETSC_FALSE));
      }
    }
    /* Update the corresponding vectors V(:,idx) = V*Q(:,idx) */
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&U));
    PetscCall(BVMultInPlace(eps->V,U,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&U));

    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      PetscCall(BVCopyColumn(eps->V,nv,k+l));  /* copy restart vector from the last column */
      if (eps->stop==EPS_STOP_THRESHOLD && nv-k<5) {  /* reallocate */
        eps->ncv = eps->mpd+k;
        PetscCall(EPSReallocateSolution(eps,eps->ncv+1));
        for (i=nv;i<eps->ncv;i++) eps->perm[i] = i;
        PetscCall(DSReallocate(eps->ds,eps->ncv+1));
        PetscCall(DSGetLeadingDimension(eps->ds,&ld));
      }
    }

    eps->nconv = k;
    PetscCall(EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv));
  }

  if (harmonic) PetscCall(PetscFree(g));
  PetscCall(DSTruncate(eps->ds,eps->nconv,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurSetRestart_KrylovSchur(EPS eps,PetscReal keep)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (keep==(PetscReal)PETSC_DEFAULT || keep==(PetscReal)PETSC_DECIDE) ctx->keep = 0.5;
  else {
    PetscCheck(keep>=0.1 && keep<=0.9,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument %g must be in the range [0.1,0.9]",(double)keep);
    ctx->keep = keep;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurSetRestart - Sets the restart parameter for the Krylov-Schur
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -eps_krylovschur_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: EPSKrylovSchurGetRestart()
@*/
PetscErrorCode EPSKrylovSchurSetRestart(EPS eps,PetscReal keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,keep,2);
  PetscTryMethod(eps,"EPSKrylovSchurSetRestart_C",(EPS,PetscReal),(eps,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetRestart_KrylovSchur(EPS eps,PetscReal *keep)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *keep = ctx->keep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetRestart - Gets the restart parameter used in the
   Krylov-Schur method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: EPSKrylovSchurSetRestart()
@*/
PetscErrorCode EPSKrylovSchurGetRestart(EPS eps,PetscReal *keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(keep,2);
  PetscUseMethod(eps,"EPSKrylovSchurGetRestart_C",(EPS,PetscReal*),(eps,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurSetLocking_KrylovSchur(EPS eps,PetscBool lock)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurSetLocking - Choose between locking and non-locking variants of
   the Krylov-Schur method.

   Logically Collective

   Input Parameters:
+  eps  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -eps_krylovschur_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: EPSKrylovSchurGetLocking()
@*/
PetscErrorCode EPSKrylovSchurSetLocking(EPS eps,PetscBool lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,lock,2);
  PetscTryMethod(eps,"EPSKrylovSchurSetLocking_C",(EPS,PetscBool),(eps,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetLocking_KrylovSchur(EPS eps,PetscBool *lock)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetLocking - Gets the locking flag used in the Krylov-Schur
   method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: EPSKrylovSchurSetLocking()
@*/
PetscErrorCode EPSKrylovSchurGetLocking(EPS eps,PetscBool *lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(lock,2);
  PetscUseMethod(eps,"EPSKrylovSchurGetLocking_C",(EPS,PetscBool*),(eps,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurSetPartitions_KrylovSchur(EPS eps,PetscInt npart)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscMPIInt     size;
  PetscInt        newnpart;

  PetscFunctionBegin;
  if (npart == PETSC_DEFAULT || npart == PETSC_DECIDE) {
    newnpart = 1;
  } else {
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)eps),&size));
    PetscCheck(npart>0 && npart<=size,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of npart");
    newnpart = npart;
  }
  if (ctx->npart!=newnpart) {
    if (ctx->npart>1) {
      PetscCall(PetscSubcommDestroy(&ctx->subc));
      if (ctx->commset) {
        PetscCallMPI(MPI_Comm_free(&ctx->commrank));
        ctx->commset = PETSC_FALSE;
      }
    }
    PetscCall(EPSDestroy(&ctx->eps));
    ctx->npart = newnpart;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurSetPartitions - Sets the number of partitions for the
   case of doing spectrum slicing for a computational interval with the
   communicator split in several sub-communicators.

   Logically Collective

   Input Parameters:
+  eps   - the eigenproblem solver context
-  npart - number of partitions

   Options Database Key:
.  -eps_krylovschur_partitions <npart> - Sets the number of partitions

   Notes:
   By default, npart=1 so all processes in the communicator participate in
   the processing of the whole interval. If npart>1 then the interval is
   divided into npart subintervals, each of them being processed by a
   subset of processes.

   The interval is split proportionally unless the separation points are
   specified with EPSKrylovSchurSetSubintervals().

   Level: advanced

.seealso: EPSKrylovSchurSetSubintervals(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurSetPartitions(EPS eps,PetscInt npart)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,npart,2);
  PetscTryMethod(eps,"EPSKrylovSchurSetPartitions_C",(EPS,PetscInt),(eps,npart));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetPartitions_KrylovSchur(EPS eps,PetscInt *npart)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *npart  = ctx->npart;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetPartitions - Gets the number of partitions of the
   communicator in case of spectrum slicing.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  npart - number of partitions

   Level: advanced

.seealso: EPSKrylovSchurSetPartitions()
@*/
PetscErrorCode EPSKrylovSchurGetPartitions(EPS eps,PetscInt *npart)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(npart,2);
  PetscUseMethod(eps,"EPSKrylovSchurGetPartitions_C",(EPS,PetscInt*),(eps,npart));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurSetDetectZeros_KrylovSchur(EPS eps,PetscBool detect)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  ctx->detect = detect;
  eps->state  = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurSetDetectZeros - Sets a flag to enforce detection of
   zeros during the factorizations throughout the spectrum slicing computation.

   Logically Collective

   Input Parameters:
+  eps    - the eigenproblem solver context
-  detect - check for zeros

   Options Database Key:
.  -eps_krylovschur_detect_zeros - Check for zeros; this takes an optional
   bool value (0/1/no/yes/true/false)

   Notes:
   A zero in the factorization indicates that a shift coincides with an eigenvalue.

   This flag is turned off by default, and may be necessary in some cases,
   especially when several partitions are being used. This feature currently
   requires an external package for factorizations with support for zero
   detection, e.g. MUMPS.

   Level: advanced

.seealso: EPSKrylovSchurSetPartitions(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurSetDetectZeros(EPS eps,PetscBool detect)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,detect,2);
  PetscTryMethod(eps,"EPSKrylovSchurSetDetectZeros_C",(EPS,PetscBool),(eps,detect));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetDetectZeros_KrylovSchur(EPS eps,PetscBool *detect)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *detect = ctx->detect;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetDetectZeros - Gets the flag that enforces zero detection
   in spectrum slicing.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  detect - whether zeros detection is enforced during factorizations

   Level: advanced

.seealso: EPSKrylovSchurSetDetectZeros()
@*/
PetscErrorCode EPSKrylovSchurGetDetectZeros(EPS eps,PetscBool *detect)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(detect,2);
  PetscUseMethod(eps,"EPSKrylovSchurGetDetectZeros_C",(EPS,PetscBool*),(eps,detect));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurSetDimensions_KrylovSchur(EPS eps,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (nev != PETSC_CURRENT) {
    PetscCheck(nev>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
    ctx->nev = nev;
  }
  if (ncv == PETSC_DETERMINE) {
    ctx->ncv = PETSC_DETERMINE;
  } else if (ncv != PETSC_CURRENT) {
    PetscCheck(ncv>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    ctx->ncv = ncv;
  }
  if (mpd == PETSC_DETERMINE) {
    ctx->mpd = PETSC_DETERMINE;
  } else if (mpd != PETSC_CURRENT) {
    PetscCheck(mpd>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    ctx->mpd = mpd;
  }
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurSetDimensions - Sets the dimensions used for each subsolve
   step in case of doing spectrum slicing for a computational interval.
   The meaning of the parameters is the same as in EPSSetDimensions().

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the subsolve
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Key:
+  -eps_krylovschur_nev <nev> - Sets the number of eigenvalues
.  -eps_krylovschur_ncv <ncv> - Sets the dimension of the subspace
-  -eps_krylovschur_mpd <mpd> - Sets the maximum projected dimension

   Note:
   Use PETSC_DETERMINE for ncv and mpd to assign a default value. For any
   of the arguments, use PETSC_CURRENT to preserve the current value.

   Level: advanced

.seealso: EPSKrylovSchurGetDimensions(), EPSSetDimensions(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurSetDimensions(EPS eps,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nev,2);
  PetscValidLogicalCollectiveInt(eps,ncv,3);
  PetscValidLogicalCollectiveInt(eps,mpd,4);
  PetscTryMethod(eps,"EPSKrylovSchurSetDimensions_C",(EPS,PetscInt,PetscInt,PetscInt),(eps,nev,ncv,mpd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetDimensions_KrylovSchur(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (nev) *nev = ctx->nev;
  if (ncv) *ncv = ctx->ncv;
  if (mpd) *mpd = ctx->mpd;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetDimensions - Gets the dimensions used for each subsolve
   step in case of doing spectrum slicing for a computational interval.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the subsolve
-  mpd - the maximum dimension allowed for the projected problem

   Level: advanced

.seealso: EPSKrylovSchurSetDimensions()
@*/
PetscErrorCode EPSKrylovSchurGetDimensions(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSKrylovSchurGetDimensions_C",(EPS,PetscInt*,PetscInt*,PetscInt*),(eps,nev,ncv,mpd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurSetSubintervals_KrylovSchur(EPS eps,PetscReal* subint)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i;

  PetscFunctionBegin;
  PetscCheck(subint[0]==eps->inta && subint[ctx->npart]==eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"First and last values must match the endpoints of EPSSetInterval()");
  for (i=0;i<ctx->npart;i++) PetscCheck(subint[i]<=subint[i+1],PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Array must contain values in ascending order");
  if (ctx->subintervals) PetscCall(PetscFree(ctx->subintervals));
  PetscCall(PetscMalloc1(ctx->npart+1,&ctx->subintervals));
  for (i=0;i<ctx->npart+1;i++) ctx->subintervals[i] = subint[i];
  ctx->subintset = PETSC_TRUE;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurSetSubintervals - Sets the points that delimit the
   subintervals to be used in spectrum slicing with several partitions.

   Logically Collective

   Input Parameters:
+  eps    - the eigenproblem solver context
-  subint - array of real values specifying subintervals

   Notes:
   This function must be called after EPSKrylovSchurSetPartitions(). For npart
   partitions, the argument subint must contain npart+1 real values sorted in
   ascending order, subint_0, subint_1, ..., subint_npart, where the first
   and last values must coincide with the interval endpoints set with
   EPSSetInterval().

   The subintervals are then defined by two consecutive points [subint_0,subint_1],
   [subint_1,subint_2], and so on.

   Level: advanced

.seealso: EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubintervals(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurSetSubintervals(EPS eps,PetscReal subint[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(subint,2);
  PetscTryMethod(eps,"EPSKrylovSchurSetSubintervals_C",(EPS,PetscReal*),(eps,subint));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetSubintervals_KrylovSchur(EPS eps,PetscReal **subint)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i;

  PetscFunctionBegin;
  if (!ctx->subintset) {
    PetscCheck(eps->state,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
    PetscCheck(ctx->sr,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  }
  PetscCall(PetscMalloc1(ctx->npart+1,subint));
  for (i=0;i<=ctx->npart;i++) (*subint)[i] = ctx->subintervals[i];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSKrylovSchurGetSubintervals - Returns the points that delimit the
   subintervals used in spectrum slicing with several partitions.

   Not Collective

   Input Parameter:
.  eps    - the eigenproblem solver context

   Output Parameter:
.  subint - array of real values specifying subintervals

   Notes:
   If the user passed values with EPSKrylovSchurSetSubintervals(), then the
   same values are returned. Otherwise, the values computed internally are
   obtained.

   This function is only available for spectrum slicing runs.

   The returned array has length npart+1 (see EPSKrylovSchurGetPartitions())
   and should be freed by the user.

   Fortran Notes:
   The calling sequence from Fortran is
.vb
   EPSKrylovSchurGetSubintervals(eps,subint,ierr)
   double precision subint(npart+1) output
.ve

   Level: advanced

.seealso: EPSKrylovSchurSetSubintervals(), EPSKrylovSchurGetPartitions(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurGetSubintervals(EPS eps,PetscReal **subint)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(subint,2);
  PetscUseMethod(eps,"EPSKrylovSchurGetSubintervals_C",(EPS,PetscReal**),(eps,subint));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetInertias_KrylovSchur(EPS eps,PetscInt *n,PetscReal **shifts,PetscInt **inertias)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,numsh;
  EPS_SR          sr = ctx->sr;

  PetscFunctionBegin;
  PetscCheck(eps->state,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  PetscCheck(ctx->sr,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  switch (eps->state) {
  case EPS_STATE_INITIAL:
    break;
  case EPS_STATE_SETUP:
    numsh = ctx->npart+1;
    if (n) *n = numsh;
    if (shifts) {
      PetscCall(PetscMalloc1(numsh,shifts));
      (*shifts)[0] = eps->inta;
      if (ctx->npart==1) (*shifts)[1] = eps->intb;
      else for (i=1;i<numsh;i++) (*shifts)[i] = ctx->subintervals[i];
    }
    if (inertias) {
      PetscCall(PetscMalloc1(numsh,inertias));
      (*inertias)[0] = (sr->dir==1)?sr->inertia0:sr->inertia1;
      if (ctx->npart==1) (*inertias)[1] = (sr->dir==1)?sr->inertia1:sr->inertia0;
      else for (i=1;i<numsh;i++) (*inertias)[i] = (*inertias)[i-1]+ctx->nconv_loc[i-1];
    }
    break;
  case EPS_STATE_SOLVED:
  case EPS_STATE_EIGENVECTORS:
    numsh = ctx->nshifts;
    if (n) *n = numsh;
    if (shifts) {
      PetscCall(PetscMalloc1(numsh,shifts));
      for (i=0;i<numsh;i++) (*shifts)[i] = ctx->shifts[i];
    }
    if (inertias) {
      PetscCall(PetscMalloc1(numsh,inertias));
      for (i=0;i<numsh;i++) (*inertias)[i] = ctx->inertias[i];
    }
    break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSKrylovSchurGetInertias - Gets the values of the shifts and their
   corresponding inertias in case of doing spectrum slicing for a
   computational interval.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  n        - number of shifts, including the endpoints of the interval
.  shifts   - the values of the shifts used internally in the solver
-  inertias - the values of the inertia in each shift

   Notes:
   If called after EPSSolve(), all shifts used internally by the solver are
   returned (including both endpoints and any intermediate ones). If called
   before EPSSolve() and after EPSSetUp() then only the information of the
   endpoints of subintervals is available.

   This function is only available for spectrum slicing runs.

   The returned arrays should be freed by the user. Can pass NULL in any of
   the two arrays if not required.

   Fortran Notes:
   The calling sequence from Fortran is
.vb
   EPSKrylovSchurGetInertias(eps,n,shifts,inertias,ierr)
   integer n
   double precision shifts(*)
   integer inertias(*)
.ve
   The arrays should be at least of length n. The value of n can be determined
   by an initial call
.vb
   EPSKrylovSchurGetInertias(eps,n,PETSC_NULL_REAL_ARRAY,PETSC_NULL_INTEGER_ARRAY,ierr)
.ve

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetSubintervals()
@*/
PetscErrorCode EPSKrylovSchurGetInertias(EPS eps,PetscInt *n,PetscReal *shifts[],PetscInt *inertias[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(n,2);
  PetscUseMethod(eps,"EPSKrylovSchurGetInertias_C",(EPS,PetscInt*,PetscReal**,PetscInt**),(eps,n,shifts,inertias));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetSubcommInfo_KrylovSchur(EPS eps,PetscInt *k,PetscInt *n,Vec *v)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  EPS_SR          sr = ((EPS_KRYLOVSCHUR*)ctx->eps->data)->sr;

  PetscFunctionBegin;
  PetscCheck(eps->state,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  PetscCheck(ctx->sr,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  if (k) *k = (ctx->npart==1)? 0: ctx->subc->color;
  if (n) *n = sr->numEigs;
  if (v) PetscCall(BVCreateVec(sr->V,v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetSubcommInfo - Gets information related to the case of
   doing spectrum slicing for a computational interval with multiple
   communicators.

   Collective on the subcommunicator (if v is given)

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  k - index of the subinterval for the calling process
.  n - number of eigenvalues found in the k-th subinterval
-  v - a vector owned by processes in the subcommunicator with dimensions
       compatible for locally computed eigenvectors (or NULL)

   Notes:
   This function is only available for spectrum slicing runs.

   The returned Vec should be destroyed by the user.

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubcommPairs()
@*/
PetscErrorCode EPSKrylovSchurGetSubcommInfo(EPS eps,PetscInt *k,PetscInt *n,Vec *v)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSKrylovSchurGetSubcommInfo_C",(EPS,PetscInt*,PetscInt*,Vec*),(eps,k,n,v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetSubcommPairs_KrylovSchur(EPS eps,PetscInt i,PetscScalar *eig,Vec v)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  EPS_SR          sr = ((EPS_KRYLOVSCHUR*)ctx->eps->data)->sr;

  PetscFunctionBegin;
  EPSCheckSolved(eps,1);
  PetscCheck(ctx->sr,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  PetscCheck(i>=0 && i<sr->numEigs,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  if (eig) *eig = sr->eigr[sr->perm[i]];
  if (v) PetscCall(BVCopyVec(sr->V,sr->perm[i],v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetSubcommPairs - Gets the i-th eigenpair stored
   internally in the subcommunicator to which the calling process belongs.

   Collective on the subcommunicator (if v is given)

   Input Parameters:
+  eps - the eigenproblem solver context
-  i   - index of the solution

   Output Parameters:
+  eig - the eigenvalue
-  v   - the eigenvector

   Notes:
   It is allowed to pass NULL for v if the eigenvector is not required.
   Otherwise, the caller must provide a valid Vec objects, i.e.,
   it must be created by the calling program with EPSKrylovSchurGetSubcommInfo().

   The index i should be a value between 0 and n-1, where n is the number of
   vectors in the local subinterval, see EPSKrylovSchurGetSubcommInfo().

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubcommInfo(), EPSKrylovSchurGetSubcommMats()
@*/
PetscErrorCode EPSKrylovSchurGetSubcommPairs(EPS eps,PetscInt i,PetscScalar *eig,Vec v)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (v) PetscValidLogicalCollectiveInt(v,i,2);
  PetscUseMethod(eps,"EPSKrylovSchurGetSubcommPairs_C",(EPS,PetscInt,PetscScalar*,Vec),(eps,i,eig,v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetSubcommMats_KrylovSchur(EPS eps,Mat *A,Mat *B)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  PetscCheck(ctx->sr,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  PetscCheck(eps->state,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  PetscCall(EPSGetOperators(ctx->eps,A,B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetSubcommMats - Gets the eigenproblem matrices stored
   internally in the subcommunicator to which the calling process belongs.

   Collective on the subcommunicator

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  A  - the matrix associated with the eigensystem
-  B  - the second matrix in the case of generalized eigenproblems

   Notes:
   This is the analog of EPSGetOperators(), but returns the matrices distributed
   differently (in the subcommunicator rather than in the parent communicator).

   These matrices should not be modified by the user.

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubcommInfo()
@*/
PetscErrorCode EPSKrylovSchurGetSubcommMats(EPS eps,Mat *A,Mat *B)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscTryMethod(eps,"EPSKrylovSchurGetSubcommMats_C",(EPS,Mat*,Mat*),(eps,A,B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurUpdateSubcommMats_KrylovSchur(EPS eps,PetscScalar a,PetscScalar ap,Mat Au,PetscScalar b,PetscScalar bp, Mat Bu,MatStructure str,PetscBool globalup)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data,*subctx;
  Mat             A,B=NULL,Ag,Bg=NULL;
  PetscBool       reuse=PETSC_TRUE;

  PetscFunctionBegin;
  PetscCheck(ctx->sr,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  PetscCheck(eps->state,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  PetscCall(EPSGetOperators(eps,&Ag,&Bg));
  PetscCall(EPSGetOperators(ctx->eps,&A,&B));

  PetscCall(MatScale(A,a));
  if (Au) PetscCall(MatAXPY(A,ap,Au,str));
  if (B) PetscCall(MatScale(B,b));
  if (Bu) PetscCall(MatAXPY(B,bp,Bu,str));
  PetscCall(EPSSetOperators(ctx->eps,A,B));

  /* Update stored matrix state */
  subctx = (EPS_KRYLOVSCHUR*)ctx->eps->data;
  PetscCall(MatGetState(A,&subctx->Astate));
  if (B) PetscCall(MatGetState(B,&subctx->Bstate));

  /* Update matrices in the parent communicator if requested by user */
  if (globalup) {
    if (ctx->npart>1) {
      if (!ctx->isrow) {
        PetscCall(MatGetOwnershipIS(Ag,&ctx->isrow,&ctx->iscol));
        reuse = PETSC_FALSE;
      }
      if (str==DIFFERENT_NONZERO_PATTERN || str==UNKNOWN_NONZERO_PATTERN) reuse = PETSC_FALSE;
      if (ctx->submata && !reuse) PetscCall(MatDestroyMatrices(1,&ctx->submata));
      PetscCall(MatCreateSubMatrices(A,1,&ctx->isrow,&ctx->iscol,(reuse)?MAT_REUSE_MATRIX:MAT_INITIAL_MATRIX,&ctx->submata));
      PetscCall(MatCreateMPIMatConcatenateSeqMat(((PetscObject)Ag)->comm,ctx->submata[0],PETSC_DECIDE,MAT_REUSE_MATRIX,&Ag));
      if (B) {
        if (ctx->submatb && !reuse) PetscCall(MatDestroyMatrices(1,&ctx->submatb));
        PetscCall(MatCreateSubMatrices(B,1,&ctx->isrow,&ctx->iscol,(reuse)?MAT_REUSE_MATRIX:MAT_INITIAL_MATRIX,&ctx->submatb));
        PetscCall(MatCreateMPIMatConcatenateSeqMat(((PetscObject)Bg)->comm,ctx->submatb[0],PETSC_DECIDE,MAT_REUSE_MATRIX,&Bg));
      }
    }
    PetscCall(MatGetState(Ag,&ctx->Astate));
    if (Bg) PetscCall(MatGetState(Bg,&ctx->Bstate));
  }
  PetscCall(EPSSetOperators(eps,Ag,Bg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurUpdateSubcommMats - Update the eigenproblem matrices stored
   internally in the subcommunicator to which the calling process belongs.

   Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  s   - scalar that multiplies the existing A matrix
.  a   - scalar used in the axpy operation on A
.  Au  - matrix used in the axpy operation on A
.  t   - scalar that multiplies the existing B matrix
.  b   - scalar used in the axpy operation on B
.  Bu  - matrix used in the axpy operation on B
.  str - structure flag
-  globalup - flag indicating if global matrices must be updated

   Notes:
   This function modifies the eigenproblem matrices at the subcommunicator level,
   and optionally updates the global matrices in the parent communicator. The updates
   are expressed as A <-- s*A + a*Au,  B <-- t*B + b*Bu.

   It is possible to update one of the matrices, or both.

   The matrices Au and Bu must be equal in all subcommunicators.

   The str flag is passed to the MatAXPY() operations to perform the updates.

   If globalup is true, communication is carried out to reconstruct the updated
   matrices in the parent communicator. The user must be warned that if global
   matrices are not in sync with subcommunicator matrices, the errors computed
   by EPSComputeError() will be wrong even if the computed solution is correct
   (the synchronization may be done only once at the end).

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubcommMats()
@*/
PetscErrorCode EPSKrylovSchurUpdateSubcommMats(EPS eps,PetscScalar s,PetscScalar a,Mat Au,PetscScalar t,PetscScalar b,Mat Bu,MatStructure str,PetscBool globalup)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveScalar(eps,s,2);
  PetscValidLogicalCollectiveScalar(eps,a,3);
  if (Au) PetscValidHeaderSpecific(Au,MAT_CLASSID,4);
  PetscValidLogicalCollectiveScalar(eps,t,5);
  PetscValidLogicalCollectiveScalar(eps,b,6);
  if (Bu) PetscValidHeaderSpecific(Bu,MAT_CLASSID,7);
  PetscValidLogicalCollectiveEnum(eps,str,8);
  PetscValidLogicalCollectiveBool(eps,globalup,9);
  PetscTryMethod(eps,"EPSKrylovSchurUpdateSubcommMats_C",(EPS,PetscScalar,PetscScalar,Mat,PetscScalar,PetscScalar,Mat,MatStructure,PetscBool),(eps,s,a,Au,t,b,Bu,str,globalup));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSKrylovSchurGetChildEPS(EPS eps,EPS *childeps)
{
  EPS_KRYLOVSCHUR  *ctx=(EPS_KRYLOVSCHUR*)eps->data,*ctx_local;
  Mat              A,B=NULL,Ar=NULL,Br=NULL;
  PetscMPIInt      rank;
  PetscObjectState Astate,Bstate=0;
  PetscObjectId    Aid,Bid=0;
  STType           sttype;
  PetscInt         nmat;
  const char       *prefix;
  MPI_Comm         child;

  PetscFunctionBegin;
  PetscCall(EPSGetOperators(eps,&A,&B));
  if (ctx->npart==1) {
    if (!ctx->eps) PetscCall(EPSCreate(((PetscObject)eps)->comm,&ctx->eps));
    PetscCall(EPSGetOptionsPrefix(eps,&prefix));
    PetscCall(EPSSetOptionsPrefix(ctx->eps,prefix));
    PetscCall(EPSSetOperators(ctx->eps,A,B));
  } else {
    PetscCall(MatGetState(A,&Astate));
    PetscCall(PetscObjectGetId((PetscObject)A,&Aid));
    if (B) {
      PetscCall(MatGetState(B,&Bstate));
      PetscCall(PetscObjectGetId((PetscObject)B,&Bid));
    }
    if (!ctx->subc) {
      /* Create context for subcommunicators */
      PetscCall(PetscSubcommCreate(PetscObjectComm((PetscObject)eps),&ctx->subc));
      PetscCall(PetscSubcommSetNumber(ctx->subc,ctx->npart));
      PetscCall(PetscSubcommSetType(ctx->subc,PETSC_SUBCOMM_CONTIGUOUS));
      PetscCall(PetscSubcommGetChild(ctx->subc,&child));

      /* Duplicate matrices */
      PetscCall(MatCreateRedundantMatrix(A,0,child,MAT_INITIAL_MATRIX,&Ar));
      ctx->Astate = Astate;
      ctx->Aid = Aid;
      PetscCall(MatPropagateSymmetryOptions(A,Ar));
      if (B) {
        PetscCall(MatCreateRedundantMatrix(B,0,child,MAT_INITIAL_MATRIX,&Br));
        ctx->Bstate = Bstate;
        ctx->Bid = Bid;
        PetscCall(MatPropagateSymmetryOptions(B,Br));
      }
    } else {
      PetscCall(PetscSubcommGetChild(ctx->subc,&child));
      if (ctx->Astate != Astate || (B && ctx->Bstate != Bstate) || ctx->Aid != Aid || (B && ctx->Bid != Bid)) {
        PetscCall(STGetNumMatrices(ctx->eps->st,&nmat));
        if (nmat) PetscCall(EPSGetOperators(ctx->eps,&Ar,&Br));
        PetscCall(MatCreateRedundantMatrix(A,0,child,MAT_INITIAL_MATRIX,&Ar));
        ctx->Astate = Astate;
        ctx->Aid = Aid;
        PetscCall(MatPropagateSymmetryOptions(A,Ar));
        if (B) {
          PetscCall(MatCreateRedundantMatrix(B,0,child,MAT_INITIAL_MATRIX,&Br));
          ctx->Bstate = Bstate;
          ctx->Bid = Bid;
          PetscCall(MatPropagateSymmetryOptions(B,Br));
        }
        PetscCall(EPSSetOperators(ctx->eps,Ar,Br));
        PetscCall(MatDestroy(&Ar));
        PetscCall(MatDestroy(&Br));
      }
    }

    /* Create auxiliary EPS */
    if (!ctx->eps) {
      PetscCall(EPSCreate(child,&ctx->eps));
      PetscCall(EPSGetOptionsPrefix(eps,&prefix));
      PetscCall(EPSSetOptionsPrefix(ctx->eps,prefix));
      PetscCall(EPSSetOperators(ctx->eps,Ar,Br));
      PetscCall(MatDestroy(&Ar));
      PetscCall(MatDestroy(&Br));
    }
    /* Create subcommunicator grouping processes with same rank */
    if (!ctx->commset) {
      PetscCallMPI(MPI_Comm_rank(child,&rank));
      PetscCallMPI(MPI_Comm_split(((PetscObject)eps)->comm,rank,ctx->subc->color,&ctx->commrank));
      ctx->commset = PETSC_TRUE;
    }
  }
  PetscCall(EPSSetType(ctx->eps,((PetscObject)eps)->type_name));
  PetscCall(STGetType(eps->st,&sttype));
  PetscCall(STSetType(ctx->eps->st,sttype));

  ctx_local = (EPS_KRYLOVSCHUR*)ctx->eps->data;
  ctx_local->npart = ctx->npart;
  ctx_local->global = PETSC_FALSE;
  ctx_local->eps = eps;
  ctx_local->subc = ctx->subc;
  ctx_local->commrank = ctx->commrank;
  *childeps = ctx->eps;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetKSP_KrylovSchur(EPS eps,KSP *ksp)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  ST              st;
  PetscBool       isfilt;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt));
  PetscCheck(eps->which==EPS_ALL && !isfilt,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations with spectrum slicing");
  PetscCall(EPSKrylovSchurGetChildEPS(eps,&ctx->eps));
  PetscCall(EPSGetST(ctx->eps,&st));
  PetscCall(STGetOperator(st,NULL));
  PetscCall(STGetKSP(st,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetKSP - Retrieve the linear solver object associated with the
   internal EPS object in case of doing spectrum slicing for a computational interval.

   Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  ksp - the internal KSP object

   Notes:
   When invoked to compute all eigenvalues in an interval with spectrum
   slicing, EPSKRYLOVSCHUR creates another EPS object internally that is
   used to compute eigenvalues by chunks near selected shifts. This function
   allows access to the KSP object associated to this internal EPS object.

   This function is only available for spectrum slicing runs. In case of
   having more than one partition, the returned KSP will be different
   in MPI processes belonging to different partitions. Hence, if required,
   EPSKrylovSchurSetPartitions() must be called BEFORE this function.

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions()
@*/
PetscErrorCode EPSKrylovSchurGetKSP(EPS eps,KSP *ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSKrylovSchurGetKSP_C",(EPS,KSP*),(eps,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurSetBSEType_KrylovSchur(EPS eps,EPSKrylovSchurBSEType bse)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  switch (bse) {
    case EPS_KRYLOVSCHUR_BSE_SHAO:
    case EPS_KRYLOVSCHUR_BSE_GRUNING:
    case EPS_KRYLOVSCHUR_BSE_PROJECTEDBSE:
      if (ctx->bse != bse) {
        ctx->bse = bse;
        eps->state = EPS_STATE_INITIAL;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid BSE type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurSetBSEType - Sets the method to be used for BSE structured eigenproblems
   in the Krylov-Schur solver.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  bse - the BSE method

   Options Database Key:
.  -eps_krylovschur_bse_type - Sets the BSE type (either 'shao', 'gruning', or 'projectedbse')

   Level: advanced

.seealso: EPSKrylovSchurGetBSEType(), EPSKrylovSchurBSEType, MatCreateBSE()
@*/
PetscErrorCode EPSKrylovSchurSetBSEType(EPS eps,EPSKrylovSchurBSEType bse)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,bse,2);
  PetscTryMethod(eps,"EPSKrylovSchurSetBSEType_C",(EPS,EPSKrylovSchurBSEType),(eps,bse));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSKrylovSchurGetBSEType_KrylovSchur(EPS eps,EPSKrylovSchurBSEType *bse)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *bse = ctx->bse;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSKrylovSchurGetBSEType - Gets the method used for BSE structured eigenproblems
   in the Krylov-Schur solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  bse - the BSE method

   Level: advanced

.seealso: EPSKrylovSchurSetBSEType(), EPSKrylovSchurBSEType, MatCreateBSE()
@*/
PetscErrorCode EPSKrylovSchurGetBSEType(EPS eps,EPSKrylovSchurBSEType *bse)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(bse,2);
  PetscUseMethod(eps,"EPSKrylovSchurGetBSEType_C",(EPS,EPSKrylovSchurBSEType*),(eps,bse));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_KrylovSchur(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  EPS_KRYLOVSCHUR       *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscBool             flg,lock,b,f1,f2,f3,isfilt;
  PetscReal             keep;
  PetscInt              i,j,k;
  KSP                   ksp;
  EPSKrylovSchurBSEType bse;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS Krylov-Schur Options");

    PetscCall(PetscOptionsReal("-eps_krylovschur_restart","Proportion of vectors kept after restart","EPSKrylovSchurSetRestart",0.5,&keep,&flg));
    if (flg) PetscCall(EPSKrylovSchurSetRestart(eps,keep));

    PetscCall(PetscOptionsBool("-eps_krylovschur_locking","Choose between locking and non-locking variants","EPSKrylovSchurSetLocking",PETSC_TRUE,&lock,&flg));
    if (flg) PetscCall(EPSKrylovSchurSetLocking(eps,lock));

    i = ctx->npart;
    PetscCall(PetscOptionsInt("-eps_krylovschur_partitions","Number of partitions of the communicator for spectrum slicing","EPSKrylovSchurSetPartitions",ctx->npart,&i,&flg));
    if (flg) PetscCall(EPSKrylovSchurSetPartitions(eps,i));

    b = ctx->detect;
    PetscCall(PetscOptionsBool("-eps_krylovschur_detect_zeros","Check zeros during factorizations at subinterval boundaries","EPSKrylovSchurSetDetectZeros",ctx->detect,&b,&flg));
    if (flg) PetscCall(EPSKrylovSchurSetDetectZeros(eps,b));

    i = 1;
    j = k = PETSC_DECIDE;
    PetscCall(PetscOptionsInt("-eps_krylovschur_nev","Number of eigenvalues to compute in each subsolve (only for spectrum slicing)","EPSKrylovSchurSetDimensions",40,&i,&f1));
    PetscCall(PetscOptionsInt("-eps_krylovschur_ncv","Number of basis vectors in each subsolve (only for spectrum slicing)","EPSKrylovSchurSetDimensions",80,&j,&f2));
    PetscCall(PetscOptionsInt("-eps_krylovschur_mpd","Maximum dimension of projected problem in each subsolve (only for spectrum slicing)","EPSKrylovSchurSetDimensions",80,&k,&f3));
    if (f1 || f2 || f3) PetscCall(EPSKrylovSchurSetDimensions(eps,i,j,k));

    PetscCall(PetscOptionsEnum("-eps_krylovschur_bse_type","Method for BSE structured eigenproblems","EPSKrylovSchurSetBSEType",EPSKrylovSchurBSETypes,(PetscEnum)ctx->bse,(PetscEnum*)&bse,&flg));
    if (flg) PetscCall(EPSKrylovSchurSetBSEType(eps,bse));

  PetscOptionsHeadEnd();

  /* set options of child KSP in spectrum slicing */
  if (eps->which==EPS_ALL) {
    if (!eps->st) PetscCall(EPSGetST(eps,&eps->st));
    PetscCall(EPSSetDefaultST(eps));
    PetscCall(STSetFromOptions(eps->st));  /* need to advance this to check ST type */
    PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt));
    if (!isfilt) {
      PetscCall(EPSKrylovSchurGetKSP_KrylovSchur(eps,&ksp));
      PetscCall(KSPSetFromOptions(ksp));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_KrylovSchur(EPS eps,PetscViewer viewer)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscBool       isascii,isfilt;
  KSP             ksp;
  PetscViewer     sviewer;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*ctx->keep)));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-"));
    if (eps->problem_type==EPS_BSE) PetscCall(PetscViewerASCIIPrintf(viewer,"  BSE method: %s\n",EPSKrylovSchurBSETypes[ctx->bse]));
    if (eps->which==EPS_ALL) {
      PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt));
      if (isfilt) PetscCall(PetscViewerASCIIPrintf(viewer,"  using filtering to extract all eigenvalues in an interval\n"));
      else {
        PetscCall(PetscViewerASCIIPrintf(viewer,"  doing spectrum slicing with nev=%" PetscInt_FMT ", ncv=%" PetscInt_FMT ", mpd=%" PetscInt_FMT "\n",ctx->nev,ctx->ncv,ctx->mpd));
        if (ctx->npart>1) {
          PetscCall(PetscViewerASCIIPrintf(viewer,"  multi-communicator spectrum slicing with %" PetscInt_FMT " partitions\n",ctx->npart));
          if (ctx->detect) PetscCall(PetscViewerASCIIPrintf(viewer,"  detecting zeros when factorizing at subinterval boundaries\n"));
        }
        /* view child KSP */
        PetscCall(EPSKrylovSchurGetKSP_KrylovSchur(eps,&ksp));
        PetscCall(PetscViewerASCIIPushTab(viewer));
        if (ctx->npart>1 && ctx->subc) {
          PetscCall(PetscViewerGetSubViewer(viewer,ctx->subc->child,&sviewer));
          if (!ctx->subc->color) PetscCall(KSPView(ksp,sviewer));
          PetscCall(PetscViewerFlush(sviewer));
          PetscCall(PetscViewerRestoreSubViewer(viewer,ctx->subc->child,&sviewer));
          /* extra call needed because of the two calls to PetscViewerASCIIPushSynchronized() in PetscViewerGetSubViewer() */
          PetscCall(PetscViewerASCIIPopSynchronized(viewer));
        } else PetscCall(KSPView(ksp,viewer));
        PetscCall(PetscViewerASCIIPopTab(viewer));
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_KrylovSchur(EPS eps)
{
  PetscBool      isfilt;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt));
  if (eps->which==EPS_ALL && !isfilt) PetscCall(EPSDestroy_KrylovSchur_Slice(eps));
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetPartitions_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetPartitions_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetDetectZeros_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetDetectZeros_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetDimensions_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetDimensions_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetSubintervals_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubintervals_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetInertias_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommInfo_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommPairs_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommMats_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurUpdateSubcommMats_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetKSP_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetBSEType_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetBSEType_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_KrylovSchur(EPS eps)
{
  PetscBool      isfilt;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt));
  if (eps->which==EPS_ALL && !isfilt) PetscCall(EPSReset_KrylovSchur_Slice(eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetDefaultST_KrylovSchur(EPS eps)
{
  PetscFunctionBegin;
  if (eps->which==EPS_ALL) {
    if (!((PetscObject)eps->st)->type_name) PetscCall(STSetType(eps->st,STSINVERT));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_KrylovSchur(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data   = (void*)ctx;
  ctx->lock   = PETSC_TRUE;
  ctx->nev    = 1;
  ctx->ncv    = PETSC_DETERMINE;
  ctx->mpd    = PETSC_DETERMINE;
  ctx->npart  = 1;
  ctx->detect = PETSC_FALSE;
  ctx->global = PETSC_TRUE;

  eps->useds = PETSC_TRUE;

  /* solve and computevectors determined at setup */
  eps->ops->setup          = EPSSetUp_KrylovSchur;
  eps->ops->setupsort      = EPSSetUpSort_KrylovSchur;
  eps->ops->setfromoptions = EPSSetFromOptions_KrylovSchur;
  eps->ops->destroy        = EPSDestroy_KrylovSchur;
  eps->ops->reset          = EPSReset_KrylovSchur;
  eps->ops->view           = EPSView_KrylovSchur;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_KrylovSchur;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetRestart_C",EPSKrylovSchurSetRestart_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetRestart_C",EPSKrylovSchurGetRestart_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetLocking_C",EPSKrylovSchurSetLocking_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetLocking_C",EPSKrylovSchurGetLocking_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetPartitions_C",EPSKrylovSchurSetPartitions_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetPartitions_C",EPSKrylovSchurGetPartitions_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetDetectZeros_C",EPSKrylovSchurSetDetectZeros_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetDetectZeros_C",EPSKrylovSchurGetDetectZeros_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetDimensions_C",EPSKrylovSchurSetDimensions_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetDimensions_C",EPSKrylovSchurGetDimensions_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetSubintervals_C",EPSKrylovSchurSetSubintervals_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubintervals_C",EPSKrylovSchurGetSubintervals_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetInertias_C",EPSKrylovSchurGetInertias_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommInfo_C",EPSKrylovSchurGetSubcommInfo_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommPairs_C",EPSKrylovSchurGetSubcommPairs_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommMats_C",EPSKrylovSchurGetSubcommMats_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurUpdateSubcommMats_C",EPSKrylovSchurUpdateSubcommMats_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetKSP_C",EPSKrylovSchurGetKSP_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetBSEType_C",EPSKrylovSchurSetBSEType_KrylovSchur));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetBSEType_C",EPSKrylovSchurGetBSEType_KrylovSchur));
  PetscFunctionReturn(PETSC_SUCCESS);
}
