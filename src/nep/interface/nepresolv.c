/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to resolvent T^{-1}(z) = sum_i (z-lambda_i)^{-1} x_i y_i'
*/

#include <slepc/private/nepimpl.h>       /*I "slepcnep.h" I*/

typedef struct {
  NEP              nep;
  RG               rg;
  PetscScalar      omega;
  PetscScalar      *nfactor;         /* normalization factors y_i'*T'(lambda_i)*x_i */
  PetscBool        *nfactor_avail;
  PetscScalar      *dots;            /* inner products y_i'*v */
  PetscBool        *dots_avail;
  PetscObjectId    vid;
  PetscObjectState vstate;
} NEP_RESOLVENT_MATSHELL;

static PetscErrorCode MatMult_Resolvent(Mat M,Vec v,Vec r)
{
  NEP_RESOLVENT_MATSHELL *ctx;
  NEP                    nep;
  PetscInt               i,inside=1;
  PetscScalar            alpha;
  Vec                    x,y,z,w;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  nep = ctx->nep;
  w = nep->work[0];
  z = nep->work[1];
  if (((PetscObject)v)->id != ctx->vid || ((PetscObject)v)->state != ctx->vstate) {
    PetscCall(PetscArrayzero(ctx->dots_avail,ctx->nep->nconv));
    PetscCall(PetscObjectGetId((PetscObject)v,&ctx->vid));
    PetscCall(VecGetState(v,&ctx->vstate));
  }
  PetscCall(VecSet(r,0.0));
  for (i=0;i<nep->nconv;i++) {
    if (ctx->rg) PetscCall(RGCheckInside(ctx->rg,1,&nep->eigr[i],&nep->eigi[i],&inside));
    if (inside>=0) {
      PetscCall(BVGetColumn(nep->V,i,&x));
      PetscCall(BVGetColumn(nep->W,i,&y));
      PetscCall(NEPApplyJacobian(nep,nep->eigr[i],x,z,w,NULL));
      if (!ctx->dots_avail[i]) {
        PetscCall(VecDot(v,y,&ctx->dots[i]));
        ctx->dots_avail[i] = PETSC_TRUE;
      }
      if (!ctx->nfactor_avail[i]) {
        PetscCall(VecDot(w,y,&ctx->nfactor[i]));
        ctx->nfactor_avail[i] = PETSC_TRUE;
      }
      alpha = ctx->dots[i]/(ctx->nfactor[i]*(ctx->omega-nep->eigr[i]));
      PetscCall(VecAXPY(r,alpha,x));
      PetscCall(BVRestoreColumn(nep->V,i,&x));
      PetscCall(BVRestoreColumn(nep->W,i,&y));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_Resolvent(Mat M)
{
  NEP_RESOLVENT_MATSHELL *ctx;

  PetscFunctionBegin;
  if (M) {
    PetscCall(MatShellGetContext(M,&ctx));
    PetscCall(PetscFree4(ctx->nfactor,ctx->nfactor_avail,ctx->dots,ctx->dots_avail));
    PetscCall(PetscFree(ctx));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPApplyResolvent - Applies the resolvent $T^{-1}(z)$ to a given vector.

   Collective

   Input Parameters:
+  nep   - the nonlinear eigensolver context
.  rg    - optional region
.  omega - value where the resolvent must be evaluated
-  v     - input vector

   Output Parameter:
.  r     - result vector

   Note:
   The resolvent $T^{-1}(z) = \sum_i (z-\lambda_i)^{-1} x_i y_i^*$ is evaluated at
   $z=\omega$ and the matrix-vector multiplication $r = T^{-1}(\omega) v$ is computed.
   Vectors $x_i$ and $y_i$ are right and left eigenvectors, respectively, normalized
   so that $y_i^*T'(\lambda_i)x_i=1$. The sum contains only eigenvectors that have
   been previously computed with `NEPSolve()`, and if a region `rg` is given then only
   those corresponding to eigenvalues inside the region are considered.

   Level: intermediate

.seealso: [](ch:nep), `NEPGetLeftEigenvector()`, `NEPSolve()`
@*/
PetscErrorCode NEPApplyResolvent(NEP nep,RG rg,PetscScalar omega,Vec v,Vec r)
{
  NEP_RESOLVENT_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,omega,3);
  PetscValidHeaderSpecific(v,VEC_CLASSID,4);
  PetscValidHeaderSpecific(r,VEC_CLASSID,5);
  NEPCheckSolved(nep,1);

  PetscCall(PetscLogEventBegin(NEP_Resolvent,nep,0,0,0));
  if (!nep->resolvent) {
    PetscCall(PetscNew(&ctx));
    ctx->nep = nep;
    PetscCall(PetscCalloc4(nep->nconv,&ctx->nfactor,nep->nconv,&ctx->nfactor_avail,nep->nconv,&ctx->dots,nep->nconv,&ctx->dots_avail));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)nep),nep->nloc,nep->nloc,nep->n,nep->n,ctx,&nep->resolvent));
    PetscCall(MatShellSetOperation(nep->resolvent,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Resolvent));
    PetscCall(MatShellSetOperation(nep->resolvent,MATOP_DESTROY,(PetscErrorCodeFn*)MatDestroy_Resolvent));
  } else PetscCall(MatShellGetContext(nep->resolvent,&ctx));
  PetscCall(NEPComputeVectors(nep));
  PetscCall(NEPSetWorkVecs(nep,2));
  ctx->rg    = rg;
  ctx->omega = omega;
  PetscCall(MatMult(nep->resolvent,v,r));
  PetscCall(PetscLogEventEnd(NEP_Resolvent,nep,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}
