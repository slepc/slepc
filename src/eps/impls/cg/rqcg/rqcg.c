/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "rqcg"

   Method: Rayleigh Quotient Conjugate Gradient

   Algorithm:

       Conjugate Gradient minimization of the Rayleigh quotient with
       periodic Rayleigh-Ritz acceleration.

   References:

       [1] L. Bergamaschi et al., "Parallel preconditioned conjugate gradient
           optimization of the Rayleigh quotient for the solution of sparse
           eigenproblems", Appl. Math. Comput. 175(2):1694-1715, 2006.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/

static PetscErrorCode EPSSolve_RQCG(EPS);

typedef struct {
  PetscInt nrest;         /* user-provided reset parameter */
  PetscInt allocsize;     /* number of columns of work BV's allocated at setup */
  BV       AV,W,P,G;
} EPS_RQCG;

static PetscErrorCode EPSSetUp_RQCG(EPS eps)
{
  PetscInt       nmat;
  EPS_RQCG       *ctx = (EPS_RQCG*)eps->data;

  PetscFunctionBegin;
  EPSCheckHermitianDefinite(eps);
  EPSCheckNotStructured(eps);
  PetscCall(EPSSetDimensions_Default(eps,&eps->nev,&eps->ncv,&eps->mpd));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) eps->which = EPS_SMALLEST_REAL;
  PetscCheck(eps->which==EPS_SMALLEST_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver supports only smallest real eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_EXTRACTION | EPS_FEATURE_THRESHOLD);
  EPSCheckIgnored(eps,EPS_FEATURE_BALANCE);

  if (!ctx->nrest) ctx->nrest = 20;

  PetscCall(EPSAllocateSolution(eps,0));
  PetscCall(EPS_SetInnerProduct(eps));

  PetscCall(STGetNumMatrices(eps->st,&nmat));
  if (!ctx->allocsize) {
    ctx->allocsize = eps->mpd;
    PetscCall(BVDuplicateResize(eps->V,eps->mpd,&ctx->AV));
    if (nmat>1) PetscCall(BVDuplicate(ctx->AV,&ctx->W));
    PetscCall(BVDuplicate(ctx->AV,&ctx->P));
    PetscCall(BVDuplicate(ctx->AV,&ctx->G));
  } else if (ctx->allocsize!=eps->mpd) {
    ctx->allocsize = eps->mpd;
    PetscCall(BVResize(ctx->AV,eps->mpd,PETSC_FALSE));
    if (nmat>1) PetscCall(BVResize(ctx->W,eps->mpd,PETSC_FALSE));
    PetscCall(BVResize(ctx->P,eps->mpd,PETSC_FALSE));
    PetscCall(BVResize(ctx->G,eps->mpd,PETSC_FALSE));
  }
  PetscCall(DSSetType(eps->ds,DSHEP));
  PetscCall(DSAllocate(eps->ds,eps->ncv));
  PetscCall(EPSSetWorkVecs(eps,1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_RQCG(EPS eps)
{
  EPS_RQCG       *ctx = (EPS_RQCG*)eps->data;
  PetscInt       i,j,k,ld,nv,ncv = eps->ncv,kini,nmat;
  PetscScalar    *C,*gamma,g,pap,pbp,pbx,pax,nu,mu,alpha,beta;
  PetscReal      resnorm,a,b,c,d,disc,t;
  PetscBool      reset;
  Mat            A,B,Q,Q1;
  Vec            v,av,bv,p,w=eps->work[0];

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  if (nmat>1) PetscCall(STGetMatrix(eps->st,1,&B));
  else B = NULL;
  PetscCall(PetscMalloc1(eps->mpd,&gamma));

  kini = eps->nini;
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;
    nv = PetscMin(eps->nconv+eps->mpd,ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,0));
    for (;kini<nv;kini++) { /* Generate more initial vectors if necessary */
      PetscCall(BVSetRandomColumn(eps->V,kini));
      PetscCall(BVOrthonormalizeColumn(eps->V,kini,PETSC_TRUE,NULL,NULL));
    }
    reset = (eps->its>1 && (eps->its-1)%ctx->nrest==0)? PETSC_TRUE: PETSC_FALSE;

    if (reset) {
      /* Prevent BVDotVec below to use B-product, restored at the end */
      PetscCall(BVSetMatrix(eps->V,NULL,PETSC_FALSE));

      /* Compute Rayleigh quotient */
      PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));
      PetscCall(BVSetActiveColumns(ctx->AV,0,nv-eps->nconv));
      PetscCall(BVMatMult(eps->V,A,ctx->AV));
      PetscCall(DSGetArray(eps->ds,DS_MAT_A,&C));
      for (i=eps->nconv;i<nv;i++) {
        PetscCall(BVSetActiveColumns(eps->V,eps->nconv,i+1));
        PetscCall(BVGetColumn(ctx->AV,i-eps->nconv,&av));
        PetscCall(BVDotVec(eps->V,av,C+eps->nconv+i*ld));
        PetscCall(BVRestoreColumn(ctx->AV,i-eps->nconv,&av));
        for (j=eps->nconv;j<i-1;j++) C[i+j*ld] = PetscConj(C[j+i*ld]);
      }
      PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&C));
      PetscCall(DSSetState(eps->ds,DS_STATE_RAW));

      /* Solve projected problem */
      PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
      PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
      PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

      /* Update vectors V(:,idx) = V * Y(:,idx) */
      PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&Q));
      PetscCall(BVMultInPlace(eps->V,Q,eps->nconv,nv));
      PetscCall(MatDenseGetSubMatrix(Q,eps->nconv,PETSC_DECIDE,eps->nconv,PETSC_DECIDE,&Q1));
      PetscCall(BVMultInPlace(ctx->AV,Q1,0,nv-eps->nconv));
      PetscCall(MatDenseRestoreSubMatrix(Q,&Q1));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&Q));
      if (B) PetscCall(BVSetMatrix(eps->V,B,PETSC_FALSE));
    } else {
      /* No need to do Rayleigh-Ritz, just take diag(V'*A*V) */
      for (i=eps->nconv;i<nv;i++) {
        PetscCall(BVGetColumn(eps->V,i,&v));
        PetscCall(BVGetColumn(ctx->AV,i-eps->nconv,&av));
        PetscCall(MatMult(A,v,av));
        PetscCall(VecDot(av,v,eps->eigr+i));
        PetscCall(BVRestoreColumn(eps->V,i,&v));
        PetscCall(BVRestoreColumn(ctx->AV,i-eps->nconv,&av));
      }
    }

    /* Compute gradient and check convergence */
    k = -1;
    for (i=eps->nconv;i<nv;i++) {
      PetscCall(BVGetColumn(eps->V,i,&v));
      PetscCall(BVGetColumn(ctx->AV,i-eps->nconv,&av));
      PetscCall(BVGetColumn(ctx->G,i-eps->nconv,&p));
      if (B) {
        PetscCall(BVGetColumn(ctx->W,i-eps->nconv,&bv));
        PetscCall(MatMult(B,v,bv));
        PetscCall(VecWAXPY(p,-eps->eigr[i],bv,av));
        PetscCall(BVRestoreColumn(ctx->W,i-eps->nconv,&bv));
      } else PetscCall(VecWAXPY(p,-eps->eigr[i],v,av));
      PetscCall(BVRestoreColumn(eps->V,i,&v));
      PetscCall(BVRestoreColumn(ctx->AV,i-eps->nconv,&av));
      PetscCall(VecNorm(p,NORM_2,&resnorm));
      PetscCall(BVRestoreColumn(ctx->G,i-eps->nconv,&p));
      PetscCall((*eps->converged)(eps,eps->eigr[i],0.0,resnorm,&eps->errest[i],eps->convergedctx));
      if (k==-1 && eps->errest[i] >= eps->tol) k = i;
    }
    if (k==-1) k = nv;
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));

    /* The next lines are necessary to avoid DS zeroing eigr */
    PetscCall(DSGetArray(eps->ds,DS_MAT_A,&C));
    for (i=eps->nconv;i<k;i++) C[i+i*ld] = eps->eigr[i];
    PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&C));

    if (eps->reason == EPS_CONVERGED_ITERATING) {

      /* Search direction */
      for (i=0;i<nv-eps->nconv;i++) {
        PetscCall(BVGetColumn(ctx->G,i,&v));
        PetscCall(STApply(eps->st,v,w));
        PetscCall(VecDot(w,v,&g));
        PetscCall(BVRestoreColumn(ctx->G,i,&v));
        beta = (!reset && eps->its>1)? g/gamma[i]: 0.0;
        gamma[i] = g;
        PetscCall(BVGetColumn(ctx->P,i,&v));
        PetscCall(VecAXPBY(v,1.0,beta,w));
        if (i+eps->nconv>0) {
          PetscCall(BVSetActiveColumns(eps->V,0,i+eps->nconv));
          PetscCall(BVOrthogonalizeVec(eps->V,v,NULL,NULL,NULL));
        }
        PetscCall(BVRestoreColumn(ctx->P,i,&v));
      }

      /* Minimization problem */
      for (i=eps->nconv;i<nv;i++) {
        PetscCall(BVGetColumn(eps->V,i,&v));
        PetscCall(BVGetColumn(ctx->AV,i-eps->nconv,&av));
        PetscCall(BVGetColumn(ctx->P,i-eps->nconv,&p));
        PetscCall(VecDot(av,v,&nu));
        PetscCall(VecDot(av,p,&pax));
        PetscCall(MatMult(A,p,w));
        PetscCall(VecDot(w,p,&pap));
        if (B) {
          PetscCall(BVGetColumn(ctx->W,i-eps->nconv,&bv));
          PetscCall(VecDot(bv,v,&mu));
          PetscCall(VecDot(bv,p,&pbx));
          PetscCall(BVRestoreColumn(ctx->W,i-eps->nconv,&bv));
          PetscCall(MatMult(B,p,w));
          PetscCall(VecDot(w,p,&pbp));
        } else {
          PetscCall(VecDot(v,v,&mu));
          PetscCall(VecDot(v,p,&pbx));
          PetscCall(VecDot(p,p,&pbp));
        }
        PetscCall(BVRestoreColumn(ctx->AV,i-eps->nconv,&av));
        a = PetscRealPart(pap*pbx-pax*pbp);
        b = PetscRealPart(nu*pbp-mu*pap);
        c = PetscRealPart(mu*pax-nu*pbx);
        t = PetscMax(PetscMax(PetscAbsReal(a),PetscAbsReal(b)),PetscAbsReal(c));
        if (t!=0.0) { a /= t; b /= t; c /= t; }
        disc = b*b-4.0*a*c;
        d = PetscSqrtReal(PetscAbsReal(disc));
        if (b>=0.0 && a!=0.0) alpha = (b+d)/(2.0*a);
        else if (b!=d) alpha = 2.0*c/(b-d);
        else alpha = 0;
        /* Next iterate */
        if (alpha!=0.0) PetscCall(VecAXPY(v,alpha,p));
        PetscCall(BVRestoreColumn(eps->V,i,&v));
        PetscCall(BVRestoreColumn(ctx->P,i-eps->nconv,&p));
        PetscCall(BVOrthonormalizeColumn(eps->V,i,PETSC_TRUE,NULL,NULL));
      }
    }

    PetscCall(EPSMonitor(eps,eps->its,k,eps->eigr,eps->eigi,eps->errest,nv));
    eps->nconv = k;
  }

  PetscCall(PetscFree(gamma));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSRQCGSetReset_RQCG(EPS eps,PetscInt nrest)
{
  EPS_RQCG *ctx = (EPS_RQCG*)eps->data;

  PetscFunctionBegin;
  if (nrest==PETSC_DEFAULT || nrest==PETSC_DECIDE) {
    ctx->nrest = 0;
    eps->state = EPS_STATE_INITIAL;
  } else {
    PetscCheck(nrest>0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Reset parameter must be >0");
    ctx->nrest = nrest;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSRQCGSetReset - Sets the reset parameter of the RQCG iteration. Every
   nrest iterations, the solver performs a Rayleigh-Ritz projection step.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  nrest - the number of iterations between resets

   Options Database Key:
.  -eps_rqcg_reset - Sets the reset parameter

   Level: advanced

.seealso: EPSRQCGGetReset()
@*/
PetscErrorCode EPSRQCGSetReset(EPS eps,PetscInt nrest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nrest,2);
  PetscTryMethod(eps,"EPSRQCGSetReset_C",(EPS,PetscInt),(eps,nrest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSRQCGGetReset_RQCG(EPS eps,PetscInt *nrest)
{
  EPS_RQCG *ctx = (EPS_RQCG*)eps->data;

  PetscFunctionBegin;
  *nrest = ctx->nrest;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSRQCGGetReset - Gets the reset parameter used in the RQCG method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  nrest - the reset parameter

   Level: advanced

.seealso: EPSRQCGSetReset()
@*/
PetscErrorCode EPSRQCGGetReset(EPS eps,PetscInt *nrest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(nrest,2);
  PetscUseMethod(eps,"EPSRQCGGetReset_C",(EPS,PetscInt*),(eps,nrest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_RQCG(EPS eps)
{
  EPS_RQCG       *ctx = (EPS_RQCG*)eps->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&ctx->AV));
  PetscCall(BVDestroy(&ctx->W));
  PetscCall(BVDestroy(&ctx->P));
  PetscCall(BVDestroy(&ctx->G));
  ctx->allocsize = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_RQCG(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  PetscBool      flg;
  PetscInt       nrest;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS RQCG Options");

    PetscCall(PetscOptionsInt("-eps_rqcg_reset","Reset parameter","EPSRQCGSetReset",20,&nrest,&flg));
    if (flg) PetscCall(EPSRQCGSetReset(eps,nrest));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_RQCG(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSRQCGSetReset_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSRQCGGetReset_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_RQCG(EPS eps,PetscViewer viewer)
{
  EPS_RQCG       *ctx = (EPS_RQCG*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) PetscCall(PetscViewerASCIIPrintf(viewer,"  reset every %" PetscInt_FMT " iterations\n",ctx->nrest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_RQCG(EPS eps)
{
  EPS_RQCG       *rqcg;

  PetscFunctionBegin;
  PetscCall(PetscNew(&rqcg));
  eps->data = (void*)rqcg;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_RQCG;
  eps->ops->setup          = EPSSetUp_RQCG;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->setfromoptions = EPSSetFromOptions_RQCG;
  eps->ops->destroy        = EPSDestroy_RQCG;
  eps->ops->reset          = EPSReset_RQCG;
  eps->ops->view           = EPSView_RQCG;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_GMRES;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSRQCGSetReset_C",EPSRQCGSetReset_RQCG));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSRQCGGetReset_C",EPSRQCGGetReset_RQCG));
  PetscFunctionReturn(PETSC_SUCCESS);
}
