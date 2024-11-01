/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "lanczos"

   Method: Explicitly restarted Lanczos

   Algorithm:

       Golub-Kahan-Lanczos bidiagonalization with explicit restart.

   References:

       [1] G.H. Golub and W. Kahan, "Calculating the singular values
           and pseudo-inverse of a matrix", SIAM J. Numer. Anal. Ser.
           B 2:205-224, 1965.

       [2] V. Hernandez, J.E. Roman, and A. Tomas, "A robust and
           efficient parallel SVD solver based on restarted Lanczos
           bidiagonalization", Elec. Trans. Numer. Anal. 31:68-85,
           2008.
*/

#include <slepc/private/svdimpl.h>                /*I "slepcsvd.h" I*/

typedef struct {
  PetscBool oneside;
} SVD_LANCZOS;

static PetscErrorCode SVDSetUp_Lanczos(SVD svd)
{
  SVD_LANCZOS    *lanczos = (SVD_LANCZOS*)svd->data;
  PetscInt       N;

  PetscFunctionBegin;
  SVDCheckStandard(svd);
  SVDCheckDefinite(svd);
  PetscCall(MatGetSize(svd->A,NULL,&N));
  PetscCall(SVDSetDimensions_Default(svd));
  PetscCheck(svd->ncv<=svd->nsv+svd->mpd,PetscObjectComm((PetscObject)svd),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (svd->max_it==PETSC_DETERMINE) svd->max_it = PetscMax(N/svd->ncv,100)*((svd->stop==SVD_STOP_THRESHOLD)?10:1);
  svd->leftbasis = PetscNot(lanczos->oneside);
  PetscCall(SVDAllocateSolution(svd,1));
  PetscCall(DSSetType(svd->ds,DSSVD));
  PetscCall(DSSetCompact(svd->ds,PETSC_TRUE));
  PetscCall(DSSetExtraRow(svd->ds,PETSC_TRUE));
  PetscCall(DSAllocate(svd->ds,svd->ncv+1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode SVDTwoSideLanczos(SVD svd,PetscReal *alpha,PetscReal *beta,BV V,BV U,PetscInt k,PetscInt *n,PetscBool *breakdown)
{
  PetscInt       i;
  Vec            u,v;
  PetscBool      lindep=PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(BVGetColumn(svd->V,k,&v));
  PetscCall(BVGetColumn(svd->U,k,&u));
  PetscCall(MatMult(svd->A,v,u));
  PetscCall(BVRestoreColumn(svd->V,k,&v));
  PetscCall(BVRestoreColumn(svd->U,k,&u));
  PetscCall(BVOrthonormalizeColumn(svd->U,k,PETSC_FALSE,alpha+k,&lindep));
  if (PetscUnlikely(lindep)) {
    *n = k;
    if (breakdown) *breakdown = lindep;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  for (i=k+1;i<*n;i++) {
    PetscCall(BVGetColumn(svd->V,i,&v));
    PetscCall(BVGetColumn(svd->U,i-1,&u));
    PetscCall(MatMult(svd->AT,u,v));
    PetscCall(BVRestoreColumn(svd->V,i,&v));
    PetscCall(BVRestoreColumn(svd->U,i-1,&u));
    PetscCall(BVOrthonormalizeColumn(svd->V,i,PETSC_FALSE,beta+i-1,&lindep));
    if (PetscUnlikely(lindep)) {
      *n = i;
      break;
    }
    PetscCall(BVGetColumn(svd->V,i,&v));
    PetscCall(BVGetColumn(svd->U,i,&u));
    PetscCall(MatMult(svd->A,v,u));
    PetscCall(BVRestoreColumn(svd->V,i,&v));
    PetscCall(BVRestoreColumn(svd->U,i,&u));
    PetscCall(BVOrthonormalizeColumn(svd->U,i,PETSC_FALSE,alpha+i,&lindep));
    if (PetscUnlikely(lindep)) {
      *n = i;
      break;
    }
  }

  if (!lindep) {
    PetscCall(BVGetColumn(svd->V,*n,&v));
    PetscCall(BVGetColumn(svd->U,*n-1,&u));
    PetscCall(MatMult(svd->AT,u,v));
    PetscCall(BVRestoreColumn(svd->V,*n,&v));
    PetscCall(BVRestoreColumn(svd->U,*n-1,&u));
    PetscCall(BVOrthogonalizeColumn(svd->V,*n,NULL,beta+*n-1,&lindep));
  }
  if (breakdown) *breakdown = lindep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDOneSideLanczos(SVD svd,PetscReal *alpha,PetscReal *beta,BV V,Vec u,Vec u_1,PetscInt k,PetscInt n,PetscScalar* work)
{
  PetscInt       i,bvl,bvk;
  PetscReal      a,b;
  Vec            z,temp;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(V,&bvl,&bvk));
  PetscCall(BVGetColumn(V,k,&z));
  PetscCall(MatMult(svd->A,z,u));
  PetscCall(BVRestoreColumn(V,k,&z));

  for (i=k+1;i<n;i++) {
    PetscCall(BVGetColumn(V,i,&z));
    PetscCall(MatMult(svd->AT,u,z));
    PetscCall(BVRestoreColumn(V,i,&z));
    PetscCall(VecNormBegin(u,NORM_2,&a));
    PetscCall(BVSetActiveColumns(V,0,i));
    PetscCall(BVDotColumnBegin(V,i,work));
    PetscCall(VecNormEnd(u,NORM_2,&a));
    PetscCall(BVDotColumnEnd(V,i,work));
    PetscCall(VecScale(u,1.0/a));
    PetscCall(BVMultColumn(V,-1.0/a,1.0/a,i,work));

    /* h = V^* z, z = z - V h  */
    PetscCall(BVDotColumn(V,i,work));
    PetscCall(BVMultColumn(V,-1.0,1.0,i,work));
    PetscCall(BVNormColumn(V,i,NORM_2,&b));
    PetscCheck(PetscAbsReal(b)>10*PETSC_MACHINE_EPSILON,PetscObjectComm((PetscObject)svd),PETSC_ERR_PLIB,"Recurrence generated a zero vector; use a two-sided variant");
    PetscCall(BVScaleColumn(V,i,1.0/b));

    PetscCall(BVGetColumn(V,i,&z));
    PetscCall(MatMult(svd->A,z,u_1));
    PetscCall(BVRestoreColumn(V,i,&z));
    PetscCall(VecAXPY(u_1,-b,u));
    alpha[i-1] = a;
    beta[i-1] = b;
    temp = u;
    u = u_1;
    u_1 = temp;
  }

  PetscCall(BVGetColumn(V,n,&z));
  PetscCall(MatMult(svd->AT,u,z));
  PetscCall(BVRestoreColumn(V,n,&z));
  PetscCall(VecNormBegin(u,NORM_2,&a));
  PetscCall(BVDotColumnBegin(V,n,work));
  PetscCall(VecNormEnd(u,NORM_2,&a));
  PetscCall(BVDotColumnEnd(V,n,work));
  PetscCall(VecScale(u,1.0/a));
  PetscCall(BVMultColumn(V,-1.0/a,1.0/a,n,work));

  /* h = V^* z, z = z - V h  */
  PetscCall(BVDotColumn(V,n,work));
  PetscCall(BVMultColumn(V,-1.0,1.0,n,work));
  PetscCall(BVNormColumn(V,i,NORM_2,&b));

  alpha[n-1] = a;
  beta[n-1] = b;
  PetscCall(BVSetActiveColumns(V,bvl,bvk));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SVDKrylovConvergence - Implements the loop that checks for convergence
   in Krylov methods.

   Input Parameters:
     svd     - the solver; some error estimates are updated in svd->errest
     getall  - whether all residuals must be computed
     kini    - initial value of k (the loop variable)
     nits    - number of iterations of the loop
     normr   - norm of triangular factor of qr([A;B]), used only in GSVD

   Output Parameter:
     kout  - the first index where the convergence test failed
*/
PetscErrorCode SVDKrylovConvergence(SVD svd,PetscBool getall,PetscInt kini,PetscInt nits,PetscReal normr,PetscInt *kout)
{
  PetscInt       k,marker,ld;
  PetscReal      *alpha,*beta,*betah,resnorm;
  PetscBool      extra;

  PetscFunctionBegin;
  if (PetscUnlikely(svd->conv == SVD_CONV_MAXIT && svd->its >= svd->max_it)) *kout = svd->nsv;
  else {
    PetscCall(DSGetLeadingDimension(svd->ds,&ld));
    PetscCall(DSGetExtraRow(svd->ds,&extra));
    PetscCheck(extra,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Only implemented for DS with extra row");
    marker = -1;
    if (svd->trackall) getall = PETSC_TRUE;
    PetscCall(DSGetArrayReal(svd->ds,DS_MAT_T,&alpha));
    beta = alpha + ld;
    betah = alpha + 2*ld;
    for (k=kini;k<kini+nits;k++) {
      if (svd->isgeneralized) resnorm = SlepcAbs(beta[k],betah[k])*normr;
      else resnorm = PetscAbsReal(beta[k]);
      PetscCall((*svd->converged)(svd,svd->sigma[k],resnorm,&svd->errest[k],svd->convergedctx));
      if (marker==-1 && svd->errest[k] >= svd->tol) marker = k;
      if (marker!=-1 && !getall) break;
    }
    PetscCall(DSRestoreArrayReal(svd->ds,DS_MAT_T,&alpha));
    if (marker!=-1) k = marker;
    *kout = k;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_Lanczos(SVD svd)
{
  SVD_LANCZOS    *lanczos = (SVD_LANCZOS*)svd->data;
  PetscReal      *alpha,*beta;
  PetscScalar    *swork,*w,*P,*aux1,*aux2;
  PetscInt       i,k,j,nv,ld;
  Vec            u=NULL,u_1=NULL;
  Mat            U,V;

  PetscFunctionBegin;
  /* allocate working space */
  PetscCall(DSGetLeadingDimension(svd->ds,&ld));
  PetscCall(PetscMalloc2(ld,&w,svd->ncv,&swork));

  if (lanczos->oneside) {
    PetscCall(MatCreateVecs(svd->A,NULL,&u));
    PetscCall(MatCreateVecs(svd->A,NULL,&u_1));
  }

  /* normalize start vector */
  if (!svd->nini) {
    PetscCall(BVSetRandomColumn(svd->V,0));
    PetscCall(BVOrthonormalizeColumn(svd->V,0,PETSC_TRUE,NULL,NULL));
  }

  while (svd->reason == SVD_CONVERGED_ITERATING) {
    svd->its++;

    /* inner loop */
    nv = PetscMin(svd->nconv+svd->mpd,svd->ncv);
    PetscCall(DSGetArrayReal(svd->ds,DS_MAT_T,&alpha));
    beta = alpha + ld;
    if (lanczos->oneside) PetscCall(SVDOneSideLanczos(svd,alpha,beta,svd->V,u,u_1,svd->nconv,nv,swork));
    else {
      PetscCall(SVDTwoSideLanczos(svd,alpha,beta,svd->V,svd->U,svd->nconv,&nv,NULL));
      PetscCall(BVSetActiveColumns(svd->U,svd->nconv,nv));
    }
    PetscCall(DSRestoreArrayReal(svd->ds,DS_MAT_T,&alpha));
    PetscCall(BVSetActiveColumns(svd->V,svd->nconv,nv));

    /* compute SVD of bidiagonal matrix */
    PetscCall(DSSetDimensions(svd->ds,nv,svd->nconv,0));
    PetscCall(DSSVDSetDimensions(svd->ds,nv));
    PetscCall(DSSetState(svd->ds,DS_STATE_INTERMEDIATE));
    PetscCall(DSSolve(svd->ds,w,NULL));
    PetscCall(DSSort(svd->ds,w,NULL,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(svd->ds));
    PetscCall(DSSynchronize(svd->ds,w,NULL));
    for (i=svd->nconv;i<nv;i++) svd->sigma[i] = PetscRealPart(w[i]);

    /* check convergence */
    PetscCall(SVDKrylovConvergence(svd,PETSC_FALSE,svd->nconv,nv-svd->nconv,1.0,&k));
    SVDSetCtxThreshold(svd,svd->sigma,k);
    PetscCall((*svd->stopping)(svd,svd->its,svd->max_it,k,svd->nsv,&svd->reason,svd->stoppingctx));

    /* compute restart vector */
    if (svd->reason == SVD_CONVERGED_ITERATING) {
      if (k<nv) {
        PetscCall(DSGetArray(svd->ds,DS_MAT_V,&P));
        for (j=svd->nconv;j<nv;j++) swork[j-svd->nconv] = PetscConj(P[j+k*ld]);
        PetscCall(DSRestoreArray(svd->ds,DS_MAT_V,&P));
        PetscCall(BVMultColumn(svd->V,1.0,0.0,nv,swork));
      } else {
        /* all approximations have converged, generate a new initial vector */
        PetscCall(BVSetRandomColumn(svd->V,nv));
        PetscCall(BVOrthonormalizeColumn(svd->V,nv,PETSC_FALSE,NULL,NULL));
      }
    }

    /* compute converged singular vectors */
    PetscCall(DSGetMat(svd->ds,DS_MAT_V,&V));
    PetscCall(BVMultInPlace(svd->V,V,svd->nconv,k));
    PetscCall(DSRestoreMat(svd->ds,DS_MAT_V,&V));
    if (!lanczos->oneside) {
      PetscCall(DSGetMat(svd->ds,DS_MAT_U,&U));
      PetscCall(BVMultInPlace(svd->U,U,svd->nconv,k));
      PetscCall(DSRestoreMat(svd->ds,DS_MAT_U,&U));
    }

    if (svd->reason == SVD_CONVERGED_ITERATING) {
      PetscCall(BVCopyColumn(svd->V,nv,k));  /* copy restart vector from the last column */
      if (svd->stop==SVD_STOP_THRESHOLD && nv-k<5) {  /* reallocate */
        svd->ncv = svd->mpd+k;
        PetscCall(SVDReallocateSolution(svd,svd->ncv+1));
        for (i=nv;i<svd->ncv;i++) svd->perm[i] = i;
        PetscCall(DSReallocate(svd->ds,svd->ncv+1));
        aux1 = w;
        aux2 = swork;
        PetscCall(PetscMalloc2(svd->ncv+1,&w,svd->ncv,&swork));
        PetscCall(PetscArraycpy(w,aux1,ld));
        PetscCall(PetscFree2(aux1,aux2));
        PetscCall(DSGetLeadingDimension(svd->ds,&ld));
      }
    }

    svd->nconv = k;
    PetscCall(SVDMonitor(svd,svd->its,svd->nconv,svd->sigma,svd->errest,nv));
  }

  /* free working space */
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&u_1));
  PetscCall(PetscFree2(w,swork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSetFromOptions_Lanczos(SVD svd,PetscOptionItems *PetscOptionsObject)
{
  PetscBool      set,val;
  SVD_LANCZOS    *lanczos = (SVD_LANCZOS*)svd->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"SVD Lanczos Options");

    PetscCall(PetscOptionsBool("-svd_lanczos_oneside","Use one-side reorthogonalization","SVDLanczosSetOneSide",lanczos->oneside,&val,&set));
    if (set) PetscCall(SVDLanczosSetOneSide(svd,val));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDLanczosSetOneSide_Lanczos(SVD svd,PetscBool oneside)
{
  SVD_LANCZOS *lanczos = (SVD_LANCZOS*)svd->data;

  PetscFunctionBegin;
  if (lanczos->oneside != oneside) {
    lanczos->oneside = oneside;
    svd->state = SVD_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDLanczosSetOneSide - Indicate if the variant of the Lanczos method
   to be used is one-sided or two-sided.

   Logically Collective

   Input Parameters:
+  svd     - singular value solver
-  oneside - boolean flag indicating if the method is one-sided or not

   Options Database Key:
.  -svd_lanczos_oneside <boolean> - Indicates the boolean flag

   Note:
   By default, a two-sided variant is selected, which is sometimes slightly
   more robust. However, the one-sided variant is faster because it avoids
   the orthogonalization associated to left singular vectors. It also saves
   the memory required for storing such vectors.

   Level: advanced

.seealso: SVDTRLanczosSetOneSide()
@*/
PetscErrorCode SVDLanczosSetOneSide(SVD svd,PetscBool oneside)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,oneside,2);
  PetscTryMethod(svd,"SVDLanczosSetOneSide_C",(SVD,PetscBool),(svd,oneside));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDLanczosGetOneSide_Lanczos(SVD svd,PetscBool *oneside)
{
  SVD_LANCZOS *lanczos = (SVD_LANCZOS*)svd->data;

  PetscFunctionBegin;
  *oneside = lanczos->oneside;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDLanczosGetOneSide - Gets if the variant of the Lanczos method
   to be used is one-sided or two-sided.

   Not Collective

   Input Parameters:
.  svd     - singular value solver

   Output Parameters:
.  oneside - boolean flag indicating if the method is one-sided or not

   Level: advanced

.seealso: SVDLanczosSetOneSide()
@*/
PetscErrorCode SVDLanczosGetOneSide(SVD svd,PetscBool *oneside)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(oneside,2);
  PetscUseMethod(svd,"SVDLanczosGetOneSide_C",(SVD,PetscBool*),(svd,oneside));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDDestroy_Lanczos(SVD svd)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(svd->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDLanczosSetOneSide_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDLanczosGetOneSide_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDView_Lanczos(SVD svd,PetscViewer viewer)
{
  SVD_LANCZOS    *lanczos = (SVD_LANCZOS*)svd->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) PetscCall(PetscViewerASCIIPrintf(viewer,"  %s-sided reorthogonalization\n",lanczos->oneside? "one": "two"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_Lanczos(SVD svd)
{
  SVD_LANCZOS    *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  svd->data = (void*)ctx;

  svd->ops->setup          = SVDSetUp_Lanczos;
  svd->ops->solve          = SVDSolve_Lanczos;
  svd->ops->destroy        = SVDDestroy_Lanczos;
  svd->ops->setfromoptions = SVDSetFromOptions_Lanczos;
  svd->ops->view           = SVDView_Lanczos;
  svd->ops->computevectors = SVDComputeVectors_Left;
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDLanczosSetOneSide_C",SVDLanczosSetOneSide_Lanczos));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDLanczosGetOneSide_C",SVDLanczosGetOneSide_Lanczos));
  PetscFunctionReturn(PETSC_SUCCESS);
}
