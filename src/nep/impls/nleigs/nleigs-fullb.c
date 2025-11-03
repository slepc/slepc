/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Full basis for the linearization of the rational approximation of non-linear eigenproblems
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include "nleigs.h"

static PetscErrorCode MatMult_FullBasis_Sinvert(Mat M,Vec x,Vec y)
{
  NEP_NLEIGS        *ctx;
  NEP               nep;
  const PetscScalar *px;
  PetscScalar       *beta,*s,*xi,*t,*py,sigma;
  PetscInt          nmat,d,i,k,m;
  Vec               xx,xxx,yy,yyy,w,ww,www;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&nep));
  ctx = (NEP_NLEIGS*)nep->data;
  beta = ctx->beta; s = ctx->s; xi = ctx->xi;
  sigma = ctx->shifts[0];
  nmat = ctx->nmat;
  d = nmat-1;
  m = nep->nloc;
  PetscCall(PetscMalloc1(ctx->nmat,&t));
  xx = ctx->w[0]; xxx = ctx->w[1]; yy = ctx->w[2]; yyy=ctx->w[3];
  w = nep->work[0]; ww = nep->work[1]; www = nep->work[2];
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  PetscCall(VecPlaceArray(xx,px+(d-1)*m));
  PetscCall(VecPlaceArray(xxx,px+(d-2)*m));
  PetscCall(VecPlaceArray(yy,py+(d-2)*m));
  PetscCall(VecCopy(xxx,yy));
  PetscCall(VecAXPY(yy,beta[d-1]/xi[d-2],xx));
  PetscCall(VecScale(yy,1.0/(s[d-2]-sigma)));
  PetscCall(VecResetArray(xx));
  PetscCall(VecResetArray(xxx));
  PetscCall(VecResetArray(yy));
  for (i=d-3;i>=0;i--) {
    PetscCall(VecPlaceArray(xx,px+(i+1)*m));
    PetscCall(VecPlaceArray(xxx,px+i*m));
    PetscCall(VecPlaceArray(yy,py+i*m));
    PetscCall(VecPlaceArray(yyy,py+(i+1)*m));
    PetscCall(VecCopy(xxx,yy));
    PetscCall(VecAXPY(yy,beta[i+1]/xi[i],xx));
    PetscCall(VecAXPY(yy,-beta[i+1]*(1.0-sigma/xi[i]),yyy));
    PetscCall(VecScale(yy,1.0/(s[i]-sigma)));
    PetscCall(VecResetArray(xx));
    PetscCall(VecResetArray(xxx));
    PetscCall(VecResetArray(yy));
    PetscCall(VecResetArray(yyy));
  }
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    PetscCall(VecZeroEntries(w));
    for (k=0;k<nep->nt;k++) {
      PetscCall(VecZeroEntries(ww));
      PetscCall(VecPlaceArray(xx,px+(d-1)*m));
      PetscCall(VecAXPY(ww,-ctx->coeffD[k+nep->nt*d]/beta[d],xx));
      PetscCall(VecResetArray(xx));
      for (i=0;i<d-1;i++) {
        PetscCall(VecPlaceArray(yy,py+i*m));
        PetscCall(VecAXPY(ww,-ctx->coeffD[nep->nt*i+k],yy));
        PetscCall(VecResetArray(yy));
      }
      PetscCall(MatMult(nep->A[k],ww,www));
      PetscCall(VecAXPY(w,1.0,www));
    }
  } else {
    PetscCall(VecPlaceArray(xx,px+(d-1)*m));
    PetscCall(MatMult(ctx->D[d],xx,w));
    PetscCall(VecScale(w,-1.0/beta[d]));
    PetscCall(VecResetArray(xx));
    for (i=0;i<d-1;i++) {
      PetscCall(VecPlaceArray(yy,py+i*m));
      PetscCall(MatMult(ctx->D[i],yy,ww));
      PetscCall(VecResetArray(yy));
      PetscCall(VecAXPY(w,-1.0,ww));
    }
  }
  PetscCall(VecPlaceArray(yy,py+(d-1)*m));
  PetscCall(KSPSolve(ctx->ksp[0],w,yy));
  PetscCall(NEPNLEIGSEvalNRTFunct(nep,d-1,sigma,t));
  for (i=0;i<d-1;i++) {
    PetscCall(VecPlaceArray(yyy,py+i*m));
    PetscCall(VecAXPY(yyy,t[i],yy));
    PetscCall(VecResetArray(yyy));
  }
  PetscCall(VecScale(yy,t[d-1]));
  PetscCall(VecResetArray(yy));
  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscCall(PetscFree(t));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultTranspose_FullBasis_Sinvert(Mat M,Vec x,Vec y)
{
  NEP_NLEIGS        *ctx;
  NEP               nep;
  const PetscScalar *px;
  PetscScalar       *beta,*s,*xi,*t,*py,sigma;
  PetscInt          nmat,d,i,k,m;
  Vec               xx,yy,yyy,w,z0;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&nep));
  ctx = (NEP_NLEIGS*)nep->data;
  beta = ctx->beta; s = ctx->s; xi = ctx->xi;
  sigma = ctx->shifts[0];
  nmat = ctx->nmat;
  d = nmat-1;
  m = nep->nloc;
  PetscCall(PetscMalloc1(ctx->nmat,&t));
  xx = ctx->w[0]; yy = ctx->w[1]; yyy=ctx->w[2];
  w = nep->work[0]; z0 = nep->work[1];
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  PetscCall(NEPNLEIGSEvalNRTFunct(nep,d,sigma,t));
  PetscCall(VecPlaceArray(xx,px+(d-1)*m));
  PetscCall(VecCopy(xx,w));
  PetscCall(VecScale(w,t[d-1]));
  PetscCall(VecResetArray(xx));
  for (i=0;i<d-1;i++) {
    PetscCall(VecPlaceArray(xx,px+i*m));
    PetscCall(VecAXPY(w,t[i],xx));
    PetscCall(VecResetArray(xx));
  }
  PetscCall(KSPSolveTranspose(ctx->ksp[0],w,z0));

  PetscCall(VecPlaceArray(yy,py));
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    PetscCall(VecZeroEntries(yy));
    for (k=0;k<nep->nt;k++) {
      PetscCall(MatMult(nep->A[k],z0,w));
      PetscCall(VecAXPY(yy,ctx->coeffD[k],w));
    }
  } else PetscCall(MatMultTranspose(ctx->D[0],z0,yy));
  PetscCall(VecPlaceArray(xx,px));
  PetscCall(VecAXPY(yy,-1.0,xx));
  PetscCall(VecResetArray(xx));
  PetscCall(VecScale(yy,-1.0/(s[0]-sigma)));
  PetscCall(VecResetArray(yy));
  for (i=2;i<d;i++) {
    PetscCall(VecPlaceArray(yy,py+(i-1)*m));
    if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
      PetscCall(VecZeroEntries(yy));
      for (k=0;k<nep->nt;k++) {
        PetscCall(MatMult(nep->A[k],z0,w));
        PetscCall(VecAXPY(yy,ctx->coeffD[k+(i-1)*nep->nt],w));
      }
    } else PetscCall(MatMultTranspose(ctx->D[i-1],z0,yy));
    PetscCall(VecPlaceArray(yyy,py+(i-2)*m));
    PetscCall(VecAXPY(yy,beta[i-1]*(1.0-sigma/xi[i-2]),yyy));
    PetscCall(VecResetArray(yyy));
    PetscCall(VecPlaceArray(xx,px+(i-1)*m));
    PetscCall(VecAXPY(yy,-1.0,xx));
    PetscCall(VecResetArray(xx));
    PetscCall(VecScale(yy,-1.0/(s[i-1]-sigma)));
    PetscCall(VecResetArray(yy));
  }
  PetscCall(VecPlaceArray(yy,py+(d-1)*m));
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    PetscCall(VecZeroEntries(yy));
    for (k=0;k<nep->nt;k++) {
      PetscCall(MatMult(nep->A[k],z0,w));
      PetscCall(VecAXPY(yy,ctx->coeffD[k+d*nep->nt],w));
    }
  } else PetscCall(MatMultTranspose(ctx->D[d],z0,yy));
  PetscCall(VecScale(yy,-1.0/beta[d]));
  PetscCall(VecPlaceArray(yyy,py+(d-2)*m));
  PetscCall(VecAXPY(yy,beta[d-1]/xi[d-2],yyy));
  PetscCall(VecResetArray(yyy));
  PetscCall(VecResetArray(yy));

  for (i=d-2;i>0;i--) {
    PetscCall(VecPlaceArray(yyy,py+(i-1)*m));
    PetscCall(VecPlaceArray(yy,py+i*m));
    PetscCall(VecAXPY(yy,beta[i]/xi[i-1],yyy));
    PetscCall(VecResetArray(yyy));
    PetscCall(VecResetArray(yy));
  }

  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscCall(PetscFree(t));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BackTransform_FullBasis(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  NEP            nep;

  PetscFunctionBegin;
  PetscCall(STShellGetContext(st,&nep));
  PetscCall(NEPNLEIGSBackTransform((PetscObject)nep,n,eigr,eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode Apply_FullBasis(ST st,Vec x,Vec y)
{
  NEP            nep;
  NEP_NLEIGS     *ctx;

  PetscFunctionBegin;
  PetscCall(STShellGetContext(st,&nep));
  ctx = (NEP_NLEIGS*)nep->data;
  PetscCall(MatMult(ctx->A,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode ApplyTranspose_FullBasis(ST st,Vec x,Vec y)
{
  NEP            nep;
  NEP_NLEIGS     *ctx;

  PetscFunctionBegin;
  PetscCall(STShellGetContext(st,&nep));
  ctx = (NEP_NLEIGS*)nep->data;
  PetscCall(MatMultTranspose(ctx->A,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPSetUp_NLEIGS_FullBasis(NEP nep)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  ST             st;
  Mat            Q;
  PetscInt       i=0,deg=ctx->nmat-1;
  PetscBool      trackall,istrivial,ks;
  PetscScalar    *epsarray,*neparray;
  Vec            veps,w=NULL;
  EPSWhich       which;

  PetscFunctionBegin;
  PetscCheck(ctx->nshifts==0,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The full-basis option is not supported with rational Krylov");
  if (!ctx->eps) PetscCall(NEPNLEIGSGetEPS(nep,&ctx->eps));
  PetscCall(EPSGetST(ctx->eps,&st));
  PetscCall(EPSSetTarget(ctx->eps,nep->target));
  PetscCall(STSetDefaultShift(st,nep->target));
  if (!((PetscObject)ctx->eps)->type_name) PetscCall(EPSSetType(ctx->eps,EPSKRYLOVSCHUR));
  else {
    PetscCall(PetscObjectTypeCompare((PetscObject)ctx->eps,EPSKRYLOVSCHUR,&ks));
    PetscCheck(ks,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Full-basis option only implemented for Krylov-Schur");
  }
  PetscCall(STSetType(st,STSHELL));
  PetscCall(STShellSetContext(st,nep));
  PetscCall(STShellSetBackTransform(st,BackTransform_FullBasis));
  PetscCall(KSPGetOperators(ctx->ksp[0],&Q,NULL));
  PetscCall(MatCreateVecsEmpty(Q,&ctx->w[0],&ctx->w[1]));
  PetscCall(MatCreateVecsEmpty(Q,&ctx->w[2],&ctx->w[3]));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)nep),deg*nep->nloc,deg*nep->nloc,deg*nep->n,deg*nep->n,nep,&ctx->A));
  PetscCall(MatShellSetOperation(ctx->A,MATOP_MULT,(PetscErrorCodeFn*)MatMult_FullBasis_Sinvert));
  PetscCall(MatShellSetOperation(ctx->A,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_FullBasis_Sinvert));
  PetscCall(STShellSetApply(st,Apply_FullBasis));
  PetscCall(STShellSetApplyTranspose(st,ApplyTranspose_FullBasis));
  PetscCall(EPSSetOperators(ctx->eps,ctx->A,NULL));
  PetscCall(EPSSetProblemType(ctx->eps,EPS_NHEP));
  switch (nep->which) {
    case NEP_TARGET_MAGNITUDE:   which = EPS_TARGET_MAGNITUDE; break;
    case NEP_TARGET_REAL:        which = EPS_TARGET_REAL; break;
    case NEP_TARGET_IMAGINARY:   which = EPS_TARGET_IMAGINARY; break;
    case NEP_WHICH_USER:         which = EPS_WHICH_USER;
      PetscCall(EPSSetEigenvalueComparison(ctx->eps,nep->sc->comparison,nep->sc->comparisonctx));
      break;
    default: SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Should set a target selection in NEPSetWhichEigenpairs()");
  }
  PetscCall(EPSSetWhichEigenpairs(ctx->eps,which));
  PetscCall(RGIsTrivial(nep->rg,&istrivial));
  if (!istrivial) PetscCall(EPSSetRG(ctx->eps,nep->rg));
  PetscCall(EPSSetDimensions(ctx->eps,nep->nev,nep->ncv,nep->mpd));
  PetscCall(EPSSetTolerances(ctx->eps,SlepcDefaultTol(nep->tol),nep->max_it));
  PetscCall(EPSSetTwoSided(ctx->eps,nep->twosided));
  /* Transfer the trackall option from pep to eps */
  PetscCall(NEPGetTrackAll(nep,&trackall));
  PetscCall(EPSSetTrackAll(ctx->eps,trackall));

  /* process initial vector */
  if (nep->nini<0) {
    PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)ctx->eps),deg*nep->nloc,deg*nep->n,&veps));
    PetscCall(VecGetArray(veps,&epsarray));
    for (i=0;i<deg;i++) {
      if (i<-nep->nini) {
        PetscCall(VecGetArray(nep->IS[i],&neparray));
        PetscCall(PetscArraycpy(epsarray+i*nep->nloc,neparray,nep->nloc));
        PetscCall(VecRestoreArray(nep->IS[i],&neparray));
      } else {
        if (!w) PetscCall(VecDuplicate(nep->IS[0],&w));
        PetscCall(VecSetRandom(w,NULL));
        PetscCall(VecGetArray(w,&neparray));
        PetscCall(PetscArraycpy(epsarray+i*nep->nloc,neparray,nep->nloc));
        PetscCall(VecRestoreArray(w,&neparray));
      }
    }
    PetscCall(VecRestoreArray(veps,&epsarray));
    PetscCall(EPSSetInitialSpace(ctx->eps,1,&veps));
    PetscCall(VecDestroy(&veps));
    PetscCall(VecDestroy(&w));
    PetscCall(SlepcBasisDestroy_Private(&nep->nini,&nep->IS));
  }

  PetscCall(EPSSetUp(ctx->eps));
  PetscCall(EPSGetDimensions(ctx->eps,NULL,&nep->ncv,&nep->mpd));
  PetscCall(EPSGetTolerances(ctx->eps,NULL,&nep->max_it));
  PetscCall(NEPAllocateSolution(nep,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   NEPNLEIGSExtract_None - Extracts the first block of the basis
   and normalizes the columns.
*/
static PetscErrorCode NEPNLEIGSExtract_None(NEP nep,EPS eps)
{
  PetscInt          i,k,m,d;
  const PetscScalar *px;
  PetscScalar       sigma=nep->target,*b;
  Mat               A;
  Vec               xxr,xxi=NULL,w,t,xx;
  PetscReal         norm;
  NEP_NLEIGS        *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  d = ctx->nmat-1;
  PetscCall(EPSGetOperators(eps,&A,NULL));
  PetscCall(MatCreateVecs(A,&xxr,NULL));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDuplicate(xxr,&xxi));
#endif
  w = nep->work[0];
  for (i=0;i<nep->nconv;i++) {
    PetscCall(EPSGetEigenvector(eps,i,xxr,xxi));
    PetscCall(VecGetArrayRead(xxr,&px));
    PetscCall(VecPlaceArray(w,px));
    PetscCall(BVInsertVec(nep->V,i,w));
    PetscCall(BVNormColumn(nep->V,i,NORM_2,&norm));
    PetscCall(BVScaleColumn(nep->V,i,1.0/norm));
    PetscCall(VecResetArray(w));
    PetscCall(VecRestoreArrayRead(xxr,&px));
  }
  if (nep->twosided) {
    PetscCall(PetscMalloc1(ctx->nmat,&b));
    PetscCall(NEPNLEIGSEvalNRTFunct(nep,d,sigma,b));
    m = nep->nloc;
    xx = ctx->w[0];
    w = nep->work[0]; t = nep->work[1];
    for (k=0;k<nep->nconv;k++) {
      PetscCall(EPSGetLeftEigenvector(eps,k,xxr,xxi));
      PetscCall(VecGetArrayRead(xxr,&px));
      PetscCall(VecPlaceArray(xx,px+(d-1)*m));
      PetscCall(VecCopy(xx,w));
      PetscCall(VecScale(w,PetscConj(b[d-1])));
      PetscCall(VecResetArray(xx));
      for (i=0;i<d-1;i++) {
        PetscCall(VecPlaceArray(xx,px+i*m));
        PetscCall(VecAXPY(w,PetscConj(b[i]),xx));
        PetscCall(VecResetArray(xx));
      }
      PetscCall(VecConjugate(w));
      PetscCall(KSPSolveTranspose(ctx->ksp[0],w,t));
      PetscCall(VecConjugate(t));
      PetscCall(BVInsertVec(nep->W,k,t));
      PetscCall(BVNormColumn(nep->W,k,NORM_2,&norm));
      PetscCall(BVScaleColumn(nep->W,k,1.0/norm));
      PetscCall(VecRestoreArrayRead(xxr,&px));
    }
    PetscCall(PetscFree(b));
  }
  PetscCall(VecDestroy(&xxr));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDestroy(&xxi));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPSolve_NLEIGS_FullBasis(NEP nep)
{
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt       i;
  PetscScalar    eigi=0.0;

  PetscFunctionBegin;
  PetscCall(EPSSolve(ctx->eps));
  PetscCall(EPSGetConverged(ctx->eps,&nep->nconv));
  PetscCall(EPSGetIterationNumber(ctx->eps,&nep->its));
  PetscCall(EPSGetConvergedReason(ctx->eps,(EPSConvergedReason*)&nep->reason));

  /* recover eigenvalues */
  for (i=0;i<nep->nconv;i++) {
    PetscCall(EPSGetEigenpair(ctx->eps,i,&nep->eigr[i],&eigi,NULL,NULL));
#if !defined(PETSC_USE_COMPLEX)
    PetscCheck(eigi==0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Complex value requires complex arithmetic");
#endif
  }
  PetscCall(NEPNLEIGSExtract_None(nep,ctx->eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPNLEIGSSetEPS_NLEIGS(NEP nep,EPS eps)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)eps));
  PetscCall(EPSDestroy(&ctx->eps));
  ctx->eps = eps;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSSetEPS - Associate an eigensolver object (`EPS`) to the NLEIGS solver.

   Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  eps - the linear eigensolver context

   Note:
   By default, the linear eigensolver is integrated within the NLEIGS method.
   This `EPS` object is used only in the case that the explicit basis
   has been selected with `NEPNLEIGSSetFullBasis()`.

   Level: advanced

.seealso: [](ch:nep), `NEPNLEIGS`, `NEPNLEIGSGetEPS()`, `NEPNLEIGSSetFullBasis()`
@*/
PetscErrorCode NEPNLEIGSSetEPS(NEP nep,EPS eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(nep,1,eps,2);
  PetscTryMethod(nep,"NEPNLEIGSSetEPS_C",(NEP,EPS),(nep,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSMonitor_NLEIGS(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  NEP            nep = (NEP)ctx;
  PetscInt       i,nv = PetscMin(nest,nep->ncv);

  PetscFunctionBegin;
  for (i=0;i<nv;i++) {
    nep->eigr[i]   = eigr[i];
    nep->eigi[i]   = eigi[i];
    nep->errest[i] = errest[i];
  }
  PetscCall(NEPNLEIGSBackTransform((PetscObject)nep,nv,nep->eigr,nep->eigi));
  PetscCall(NEPMonitor(nep,its,nconv,nep->eigr,nep->eigi,nep->errest,nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPNLEIGSGetEPS_NLEIGS(NEP nep,EPS *eps)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (!ctx->eps) {
    PetscCall(EPSCreate(PetscObjectComm((PetscObject)nep),&ctx->eps));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->eps,(PetscObject)nep,1));
    PetscCall(EPSSetOptionsPrefix(ctx->eps,((PetscObject)nep)->prefix));
    PetscCall(EPSAppendOptionsPrefix(ctx->eps,"nep_nleigs_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->eps,((PetscObject)nep)->options));
    PetscCall(EPSMonitorSet(ctx->eps,EPSMonitor_NLEIGS,nep,NULL));
  }
  *eps = ctx->eps;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSGetEPS - Retrieve the linear eigensolver object (`EPS`) associated
   to the nonlinear eigenvalue solver.

   Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  eps - the linear eigensolver context

   Note:
   By default, the linear eigensolver is integrated within the NLEIGS method.
   This `EPS` object is used only in the case that the explicit basis
   has been selected with `NEPNLEIGSSetFullBasis()`.

   Level: advanced

.seealso: [](ch:nep), `NEPNLEIGS`, `NEPNLEIGSSetEPS()`, `NEPNLEIGSSetFullBasis()`
@*/
PetscErrorCode NEPNLEIGSGetEPS(NEP nep,EPS *eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(eps,2);
  PetscUseMethod(nep,"NEPNLEIGSGetEPS_C",(NEP,EPS*),(nep,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}
