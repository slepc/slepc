/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to eigensolvers in Elemental.
*/

#include <slepc/private/epsimpl.h>    /*I "slepceps.h" I*/
#include <petsc/private/petscelemental.h>

typedef struct {
  Mat Ae,Be;        /* converted matrices */
} EPS_Elemental;

static PetscErrorCode EPSSetUp_Elemental(EPS eps)
{
  EPS_Elemental  *ctx = (EPS_Elemental*)eps->data;
  Mat            A,B;
  PetscInt       nmat;
  PetscBool      isshift;
  PetscScalar    shift;

  PetscFunctionBegin;
  EPSCheckHermitianDefinite(eps);
  EPSCheckNotStructured(eps);
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isshift));
  PetscCheck(isshift,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support spectral transformations");
  if (eps->nev==0) eps->nev = 1;
  eps->ncv = eps->n;
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = 1;
  if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  PetscCheck(eps->which!=EPS_ALL || eps->inta==eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support interval computation");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_STOPPING);
  PetscCall(EPSAllocateSolution(eps,0));

  /* convert matrices */
  PetscCall(MatDestroy(&ctx->Ae));
  PetscCall(MatDestroy(&ctx->Be));
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  PetscCall(MatConvert(A,MATELEMENTAL,MAT_INITIAL_MATRIX,&ctx->Ae));
  if (nmat>1) {
    PetscCall(STGetMatrix(eps->st,1,&B));
    PetscCall(MatConvert(B,MATELEMENTAL,MAT_INITIAL_MATRIX,&ctx->Be));
  }
  PetscCall(STGetShift(eps->st,&shift));
  if (shift != 0.0) {
    if (nmat>1) PetscCall(MatAXPY(ctx->Ae,-shift,ctx->Be,SAME_NONZERO_PATTERN));
    else PetscCall(MatShift(ctx->Ae,-shift));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_Elemental(EPS eps)
{
  EPS_Elemental  *ctx = (EPS_Elemental*)eps->data;
  Mat            A = ctx->Ae,B = ctx->Be,Q,V;
  Mat_Elemental  *a = (Mat_Elemental*)A->data,*b,*q;
  PetscInt       i,rrank,ridx,erow;

  PetscFunctionBegin;
  El::DistMatrix<PetscReal,El::VR,El::STAR> w(*a->grid);
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&Q));
  q = (Mat_Elemental*)Q->data;

  if (B) {
    b = (Mat_Elemental*)B->data;
    El::HermitianGenDefEig(El::AXBX,El::LOWER,*a->emat,*b->emat,w,*q->emat);
  } else El::HermitianEig(El::LOWER,*a->emat,w,*q->emat);

  for (i=0;i<eps->ncv;i++) {
    P2RO(A,0,i,&rrank,&ridx);
    RO2E(A,0,rrank,ridx,&erow);
    eps->eigr[i] = w.Get(erow,0);
  }
  PetscCall(BVGetMat(eps->V,&V));
  PetscCall(MatConvert(Q,MATDENSE,MAT_REUSE_MATRIX,&V));
  PetscCall(BVRestoreMat(eps->V,&V));
  PetscCall(MatDestroy(&Q));

  eps->nconv  = eps->ncv;
  eps->its    = 1;
  eps->reason = EPS_CONVERGED_TOL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_Elemental(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_Elemental(EPS eps)
{
  EPS_Elemental  *ctx = (EPS_Elemental*)eps->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->Ae));
  PetscCall(MatDestroy(&ctx->Be));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_Elemental(EPS eps)
{
  EPS_Elemental  *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->categ = EPS_CATEGORY_OTHER;

  eps->ops->solve          = EPSSolve_Elemental;
  eps->ops->setup          = EPSSetUp_Elemental;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->destroy        = EPSDestroy_Elemental;
  eps->ops->reset          = EPSReset_Elemental;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_NoFactor;
  PetscFunctionReturn(PETSC_SUCCESS);
}
