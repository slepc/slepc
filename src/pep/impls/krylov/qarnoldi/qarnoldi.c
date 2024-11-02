/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc quadratic eigensolver: "qarnoldi"

   Method: Q-Arnoldi

   Algorithm:

       Quadratic Arnoldi with Krylov-Schur type restart.

   References:

       [1] K. Meerbergen, "The Quadratic Arnoldi method for the solution
           of the quadratic eigenvalue problem", SIAM J. Matrix Anal.
           Appl. 30(4):1462-1482, 2008.
*/

#include <slepc/private/pepimpl.h>    /*I "slepcpep.h" I*/
#include <petscblaslapack.h>

typedef struct {
  PetscReal keep;         /* restart parameter */
  PetscBool lock;         /* locking/non-locking variant */
} PEP_QARNOLDI;

static PetscErrorCode PEPSetUp_QArnoldi(PEP pep)
{
  PEP_QARNOLDI   *ctx = (PEP_QARNOLDI*)pep->data;
  PetscBool      flg;

  PetscFunctionBegin;
  PEPCheckQuadratic(pep);
  PEPCheckShiftSinvert(pep);
  PetscCall(PEPSetDimensions_Default(pep,pep->nev,&pep->ncv,&pep->mpd));
  PetscCheck(ctx->lock || pep->mpd>=pep->ncv,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Should not use mpd parameter in non-locking variant");
  if (pep->max_it==PETSC_DETERMINE) pep->max_it = PetscMax(100,4*pep->n/pep->ncv);
  if (!pep->which) PetscCall(PEPSetWhichEigenpairs_Default(pep));
  PetscCheck(pep->which!=PEP_ALL,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues");

  PetscCall(STGetTransform(pep->st,&flg));
  PetscCheck(flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Solver requires the ST transformation flag set, see STSetTransform()");

  /* set default extraction */
  if (!pep->extract) pep->extract = PEP_EXTRACT_NONE;
  PEPCheckUnsupported(pep,PEP_FEATURE_NONMONOMIAL | PEP_FEATURE_EXTRACT);

  if (!ctx->keep) ctx->keep = 0.5;

  PetscCall(PEPAllocateSolution(pep,0));
  PetscCall(PEPSetWorkVecs(pep,4));

  PetscCall(DSSetType(pep->ds,DSNHEP));
  PetscCall(DSSetExtraRow(pep->ds,PETSC_TRUE));
  PetscCall(DSAllocate(pep->ds,pep->ncv+1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPExtractVectors_QArnoldi(PEP pep)
{
  PetscInt       k=pep->nconv;
  Mat            X,X0;

  PetscFunctionBegin;
  if (pep->nconv==0) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DSVectors(pep->ds,DS_MAT_X,NULL,NULL));

  /* update vectors V = V*X */
  PetscCall(DSGetMat(pep->ds,DS_MAT_X,&X));
  PetscCall(MatDenseGetSubMatrix(X,0,k,0,k,&X0));
  PetscCall(BVMultInPlace(pep->V,X0,0,k));
  PetscCall(MatDenseRestoreSubMatrix(X,&X0));
  PetscCall(DSRestoreMat(pep->ds,DS_MAT_X,&X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute a step of Classical Gram-Schmidt orthogonalization
*/
static PetscErrorCode PEPQArnoldiCGS(PEP pep,PetscScalar *H,PetscBLASInt ldh,PetscScalar *h,PetscBLASInt j,BV V,Vec t,Vec v,Vec w,PetscReal *onorm,PetscReal *norm,PetscScalar *work)
{
  PetscBLASInt   ione = 1,j_1 = j+1;
  PetscReal      x,y;
  PetscScalar    dot,one = 1.0,zero = 0.0;

  PetscFunctionBegin;
  /* compute norm of v and w */
  if (onorm) {
    PetscCall(VecNorm(v,NORM_2,&x));
    PetscCall(VecNorm(w,NORM_2,&y));
    *onorm = SlepcAbs(x,y);
  }

  /* orthogonalize: compute h */
  PetscCall(BVDotVec(V,v,h));
  PetscCall(BVDotVec(V,w,work));
  if (j>0) PetscCallBLAS("BLASgemv",BLASgemv_("C",&j_1,&j,&one,H,&ldh,work,&ione,&one,h,&ione));
  PetscCall(VecDot(w,t,&dot));
  h[j] += dot;

  /* orthogonalize: update v and w */
  PetscCall(BVMultVec(V,-1.0,1.0,v,h));
  if (j>0) {
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&j_1,&j,&one,H,&ldh,h,&ione,&zero,work,&ione));
    PetscCall(BVMultVec(V,-1.0,1.0,w,work));
  }
  PetscCall(VecAXPY(w,-h[j],t));

  /* compute norm of v and w */
  if (norm) {
    PetscCall(VecNorm(v,NORM_2,&x));
    PetscCall(VecNorm(w,NORM_2,&y));
    *norm = SlepcAbs(x,y);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute a run of Q-Arnoldi iterations
*/
static PetscErrorCode PEPQArnoldi(PEP pep,Mat A,PetscInt k,PetscInt *M,Vec v,Vec w,PetscReal *beta,PetscBool *breakdown,PetscScalar *work)
{
  PetscInt           i,j,l,m = *M,ldh;
  PetscBLASInt       jj,ldhh;
  Vec                t = pep->work[2],u = pep->work[3];
  BVOrthogRefineType refinement;
  PetscReal          norm=0.0,onorm,eta;
  PetscScalar        *H,*c = work + m;

  PetscFunctionBegin;
  *beta = 0.0;
  PetscCall(MatDenseGetArray(A,&H));
  PetscCall(MatDenseGetLDA(A,&ldh));
  PetscCall(BVGetOrthogonalization(pep->V,NULL,&refinement,&eta,NULL));
  PetscCall(BVInsertVec(pep->V,k,v));
  for (j=k;j<m;j++) {
    /* apply operator */
    PetscCall(VecCopy(w,t));
    if (pep->Dr) PetscCall(VecPointwiseMult(v,v,pep->Dr));
    PetscCall(STMatMult(pep->st,0,v,u));
    PetscCall(VecCopy(t,v));
    if (pep->Dr) PetscCall(VecPointwiseMult(t,t,pep->Dr));
    PetscCall(STMatMult(pep->st,1,t,w));
    PetscCall(VecAXPY(u,pep->sfactor,w));
    PetscCall(STMatSolve(pep->st,u,w));
    PetscCall(VecScale(w,-1.0/(pep->sfactor*pep->sfactor)));
    if (pep->Dr) PetscCall(VecPointwiseDivide(w,w,pep->Dr));
    PetscCall(VecCopy(v,t));
    PetscCall(BVSetActiveColumns(pep->V,0,j+1));

    /* orthogonalize */
    PetscCall(PetscBLASIntCast(j,&jj));
    PetscCall(PetscBLASIntCast(ldh,&ldhh));
    switch (refinement) {
      case BV_ORTHOG_REFINE_NEVER:
        PetscCall(PEPQArnoldiCGS(pep,H,ldhh,H+ldh*j,jj,pep->V,t,v,w,NULL,&norm,work));
        *breakdown = PETSC_FALSE;
        break;
      case BV_ORTHOG_REFINE_ALWAYS:
        PetscCall(PEPQArnoldiCGS(pep,H,ldhh,H+ldh*j,jj,pep->V,t,v,w,NULL,NULL,work));
        PetscCall(PEPQArnoldiCGS(pep,H,ldhh,c,jj,pep->V,t,v,w,&onorm,&norm,work));
        for (i=0;i<=j;i++) H[ldh*j+i] += c[i];
        if (norm < eta * onorm) *breakdown = PETSC_TRUE;
        else *breakdown = PETSC_FALSE;
        break;
      case BV_ORTHOG_REFINE_IFNEEDED:
        PetscCall(PEPQArnoldiCGS(pep,H,ldhh,H+ldh*j,jj,pep->V,t,v,w,&onorm,&norm,work));
        /* ||q|| < eta ||h|| */
        l = 1;
        while (l<3 && norm < eta * onorm) {
          l++;
          onorm = norm;
          PetscCall(PEPQArnoldiCGS(pep,H,ldhh,c,jj,pep->V,t,v,w,NULL,&norm,work));
          for (i=0;i<=j;i++) H[ldh*j+i] += c[i];
        }
        if (norm < eta * onorm) *breakdown = PETSC_TRUE;
        else *breakdown = PETSC_FALSE;
        break;
    }
    PetscCall(VecScale(v,1.0/norm));
    PetscCall(VecScale(w,1.0/norm));

    H[j+1+ldh*j] = norm;
    if (j<m-1) PetscCall(BVInsertVec(pep->V,j+1,v));
  }
  *beta = norm;
  PetscCall(MatDenseRestoreArray(A,&H));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSolve_QArnoldi(PEP pep)
{
  PEP_QARNOLDI   *ctx = (PEP_QARNOLDI*)pep->data;
  PetscInt       j,k,l,lwork,nv,nconv;
  Vec            v=pep->work[0],w=pep->work[1];
  Mat            Q,S;
  PetscScalar    *work;
  PetscReal      beta,norm,x,y;
  PetscBool      breakdown=PETSC_FALSE,sinv;

  PetscFunctionBegin;
  lwork = 7*pep->ncv;
  PetscCall(PetscMalloc1(lwork,&work));
  PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv));
  PetscCall(RGPushScale(pep->rg,sinv?pep->sfactor:1.0/pep->sfactor));
  PetscCall(STScaleShift(pep->st,sinv?pep->sfactor:1.0/pep->sfactor));

  /* Get the starting Arnoldi vector */
  for (j=0;j<2;j++) {
    if (j>=pep->nini) PetscCall(BVSetRandomColumn(pep->V,j));
  }
  PetscCall(BVCopyVec(pep->V,0,v));
  PetscCall(BVCopyVec(pep->V,1,w));
  PetscCall(VecNorm(v,NORM_2,&x));
  PetscCall(VecNorm(w,NORM_2,&y));
  norm = SlepcAbs(x,y);
  PetscCall(VecScale(v,1.0/norm));
  PetscCall(VecScale(w,1.0/norm));

  /* clean projected matrix (including the extra-arrow) */
  PetscCall(DSSetDimensions(pep->ds,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(DSGetMat(pep->ds,DS_MAT_A,&S));
  PetscCall(MatZeroEntries(S));
  PetscCall(DSRestoreMat(pep->ds,DS_MAT_A,&S));

   /* Restart loop */
  l = 0;
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;

    /* Compute an nv-step Arnoldi factorization */
    nv = PetscMin(pep->nconv+pep->mpd,pep->ncv);
    PetscCall(DSGetMat(pep->ds,DS_MAT_A,&S));
    PetscCall(PEPQArnoldi(pep,S,pep->nconv+l,&nv,v,w,&beta,&breakdown,work));
    PetscCall(DSRestoreMat(pep->ds,DS_MAT_A,&S));
    PetscCall(DSSetDimensions(pep->ds,nv,pep->nconv,pep->nconv+l));
    PetscCall(DSSetState(pep->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(pep->V,pep->nconv,nv));

    /* Solve projected problem */
    PetscCall(DSSolve(pep->ds,pep->eigr,pep->eigi));
    PetscCall(DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(pep->ds));
    PetscCall(DSSynchronize(pep->ds,pep->eigr,pep->eigi));

    /* Check convergence */
    PetscCall(PEPKrylovConvergence(pep,PETSC_FALSE,pep->nconv,nv-pep->nconv,beta,&k));
    PetscCall((*pep->stopping)(pep,pep->its,pep->max_it,k,pep->nev,&pep->reason,pep->stoppingctx));
    nconv = k;

    /* Update l */
    if (pep->reason != PEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
      PetscCall(DSGetTruncateSize(pep->ds,k,nv,&l));
    }
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(pep,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    if (pep->reason == PEP_CONVERGED_ITERATING) {
      if (PetscUnlikely(breakdown)) {
        /* Stop if breakdown */
        PetscCall(PetscInfo(pep,"Breakdown Quadratic Arnoldi method (it=%" PetscInt_FMT " norm=%g)\n",pep->its,(double)beta));
        pep->reason = PEP_DIVERGED_BREAKDOWN;
      } else {
        /* Prepare the Rayleigh quotient for restart */
        PetscCall(DSTruncate(pep->ds,k+l,PETSC_FALSE));
      }
    }
    /* Update the corresponding vectors V(:,idx) = V*Q(:,idx) */
    PetscCall(DSGetMat(pep->ds,DS_MAT_Q,&Q));
    PetscCall(BVMultInPlace(pep->V,Q,pep->nconv,k+l));
    PetscCall(DSRestoreMat(pep->ds,DS_MAT_Q,&Q));

    pep->nconv = k;
    PetscCall(PEPMonitor(pep,pep->its,nconv,pep->eigr,pep->eigi,pep->errest,nv));
  }
  PetscCall(BVSetActiveColumns(pep->V,0,pep->nconv));
  for (j=0;j<pep->nconv;j++) {
    pep->eigr[j] *= pep->sfactor;
    pep->eigi[j] *= pep->sfactor;
  }

  PetscCall(STScaleShift(pep->st,sinv?1.0/pep->sfactor:pep->sfactor));
  PetscCall(RGPopScale(pep->rg));

  PetscCall(DSTruncate(pep->ds,pep->nconv,PETSC_TRUE));
  PetscCall(PetscFree(work));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPQArnoldiSetRestart_QArnoldi(PEP pep,PetscReal keep)
{
  PEP_QARNOLDI *ctx = (PEP_QARNOLDI*)pep->data;

  PetscFunctionBegin;
  if (keep==(PetscReal)PETSC_DEFAULT || keep==(PetscReal)PETSC_DECIDE) ctx->keep = 0.5;
  else {
    PetscCheck(keep>=0.1 && keep<=0.9,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument must be in the range [0.1,0.9]");
    ctx->keep = keep;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPQArnoldiSetRestart - Sets the restart parameter for the Q-Arnoldi
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective

   Input Parameters:
+  pep  - the eigenproblem solver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -pep_qarnoldi_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: PEPQArnoldiGetRestart()
@*/
PetscErrorCode PEPQArnoldiSetRestart(PEP pep,PetscReal keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,keep,2);
  PetscTryMethod(pep,"PEPQArnoldiSetRestart_C",(PEP,PetscReal),(pep,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPQArnoldiGetRestart_QArnoldi(PEP pep,PetscReal *keep)
{
  PEP_QARNOLDI *ctx = (PEP_QARNOLDI*)pep->data;

  PetscFunctionBegin;
  *keep = ctx->keep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPQArnoldiGetRestart - Gets the restart parameter used in the Q-Arnoldi method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: PEPQArnoldiSetRestart()
@*/
PetscErrorCode PEPQArnoldiGetRestart(PEP pep,PetscReal *keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(keep,2);
  PetscUseMethod(pep,"PEPQArnoldiGetRestart_C",(PEP,PetscReal*),(pep,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPQArnoldiSetLocking_QArnoldi(PEP pep,PetscBool lock)
{
  PEP_QARNOLDI *ctx = (PEP_QARNOLDI*)pep->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPQArnoldiSetLocking - Choose between locking and non-locking variants of
   the Q-Arnoldi method.

   Logically Collective

   Input Parameters:
+  pep  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -pep_qarnoldi_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: PEPQArnoldiGetLocking()
@*/
PetscErrorCode PEPQArnoldiSetLocking(PEP pep,PetscBool lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,lock,2);
  PetscTryMethod(pep,"PEPQArnoldiSetLocking_C",(PEP,PetscBool),(pep,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPQArnoldiGetLocking_QArnoldi(PEP pep,PetscBool *lock)
{
  PEP_QARNOLDI *ctx = (PEP_QARNOLDI*)pep->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPQArnoldiGetLocking - Gets the locking flag used in the Q-Arnoldi method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: PEPQArnoldiSetLocking()
@*/
PetscErrorCode PEPQArnoldiGetLocking(PEP pep,PetscBool *lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(lock,2);
  PetscUseMethod(pep,"PEPQArnoldiGetLocking_C",(PEP,PetscBool*),(pep,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetFromOptions_QArnoldi(PEP pep,PetscOptionItems PetscOptionsObject)
{
  PetscBool      flg,lock;
  PetscReal      keep;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"PEP Q-Arnoldi Options");

    PetscCall(PetscOptionsReal("-pep_qarnoldi_restart","Proportion of vectors kept after restart","PEPQArnoldiSetRestart",0.5,&keep,&flg));
    if (flg) PetscCall(PEPQArnoldiSetRestart(pep,keep));

    PetscCall(PetscOptionsBool("-pep_qarnoldi_locking","Choose between locking and non-locking variants","PEPQArnoldiSetLocking",PETSC_FALSE,&lock,&flg));
    if (flg) PetscCall(PEPQArnoldiSetLocking(pep,lock));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPView_QArnoldi(PEP pep,PetscViewer viewer)
{
  PEP_QARNOLDI   *ctx = (PEP_QARNOLDI*)pep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*ctx->keep)));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPDestroy_QArnoldi(PEP pep)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(pep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPQArnoldiSetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPQArnoldiGetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPQArnoldiSetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPQArnoldiGetLocking_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_QArnoldi(PEP pep)
{
  PEP_QARNOLDI   *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  pep->data = (void*)ctx;

  pep->lineariz = PETSC_TRUE;
  ctx->lock     = PETSC_TRUE;

  pep->ops->solve          = PEPSolve_QArnoldi;
  pep->ops->setup          = PEPSetUp_QArnoldi;
  pep->ops->setfromoptions = PEPSetFromOptions_QArnoldi;
  pep->ops->destroy        = PEPDestroy_QArnoldi;
  pep->ops->view           = PEPView_QArnoldi;
  pep->ops->backtransform  = PEPBackTransform_Default;
  pep->ops->computevectors = PEPComputeVectors_Default;
  pep->ops->extractvectors = PEPExtractVectors_QArnoldi;
  pep->ops->setdefaultst   = PEPSetDefaultST_Transform;

  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPQArnoldiSetRestart_C",PEPQArnoldiSetRestart_QArnoldi));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPQArnoldiGetRestart_C",PEPQArnoldiGetRestart_QArnoldi));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPQArnoldiSetLocking_C",PEPQArnoldiSetLocking_QArnoldi));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPQArnoldiGetLocking_C",PEPQArnoldiGetLocking_QArnoldi));
  PetscFunctionReturn(PETSC_SUCCESS);
}
