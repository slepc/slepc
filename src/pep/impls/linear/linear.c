/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Explicit linearization for polynomial eigenproblems
*/

#include <slepc/private/pepimpl.h>         /*I "slepcpep.h" I*/
#include "linear.h"

static PetscErrorCode MatMult_Linear_Shift(Mat M,Vec x,Vec y)
{
  PEP_LINEAR        *ctx;
  PEP               pep;
  const PetscScalar *px;
  PetscScalar       *py,a,sigma=0.0;
  PetscInt          nmat,deg,i,m;
  Vec               x1,x2,x3,y1,aux;
  PetscReal         *ca,*cb,*cg;
  PetscBool         flg;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  pep = ctx->pep;
  PetscCall(STGetTransform(pep->st,&flg));
  if (!flg) PetscCall(STGetShift(pep->st,&sigma));
  nmat = pep->nmat;
  deg = nmat-1;
  m = pep->nloc;
  ca = pep->pbc;
  cb = pep->pbc+nmat;
  cg = pep->pbc+2*nmat;
  x1=ctx->w[0];x2=ctx->w[1];x3=ctx->w[2];y1=ctx->w[3];aux=ctx->w[4];

  PetscCall(VecSet(y,0.0));
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  a = 1.0;

  /* first block */
  PetscCall(VecPlaceArray(x2,px));
  PetscCall(VecPlaceArray(x3,px+m));
  PetscCall(VecPlaceArray(y1,py));
  PetscCall(VecAXPY(y1,cb[0]-sigma,x2));
  PetscCall(VecAXPY(y1,ca[0],x3));
  PetscCall(VecResetArray(x2));
  PetscCall(VecResetArray(x3));
  PetscCall(VecResetArray(y1));

  /* inner blocks */
  for (i=1;i<deg-1;i++) {
    PetscCall(VecPlaceArray(x1,px+(i-1)*m));
    PetscCall(VecPlaceArray(x2,px+i*m));
    PetscCall(VecPlaceArray(x3,px+(i+1)*m));
    PetscCall(VecPlaceArray(y1,py+i*m));
    PetscCall(VecAXPY(y1,cg[i],x1));
    PetscCall(VecAXPY(y1,cb[i]-sigma,x2));
    PetscCall(VecAXPY(y1,ca[i],x3));
    PetscCall(VecResetArray(x1));
    PetscCall(VecResetArray(x2));
    PetscCall(VecResetArray(x3));
    PetscCall(VecResetArray(y1));
  }

  /* last block */
  PetscCall(VecPlaceArray(y1,py+(deg-1)*m));
  for (i=0;i<deg;i++) {
    PetscCall(VecPlaceArray(x1,px+i*m));
    PetscCall(STMatMult(pep->st,i,x1,aux));
    PetscCall(VecAXPY(y1,a,aux));
    PetscCall(VecResetArray(x1));
    a *= pep->sfactor;
  }
  PetscCall(VecCopy(y1,aux));
  PetscCall(STMatSolve(pep->st,aux,y1));
  PetscCall(VecScale(y1,-ca[deg-1]/a));
  PetscCall(VecPlaceArray(x1,px+(deg-2)*m));
  PetscCall(VecPlaceArray(x2,px+(deg-1)*m));
  PetscCall(VecAXPY(y1,cg[deg-1],x1));
  PetscCall(VecAXPY(y1,cb[deg-1]-sigma,x2));
  PetscCall(VecResetArray(x1));
  PetscCall(VecResetArray(x2));
  PetscCall(VecResetArray(y1));

  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_Linear_Sinvert(Mat M,Vec x,Vec y)
{
  PEP_LINEAR        *ctx;
  PEP               pep;
  const PetscScalar *px;
  PetscScalar       *py,a,sigma,t=1.0,tp=0.0,tt;
  PetscInt          nmat,deg,i,m;
  Vec               x1,y1,y2,y3,aux,aux2;
  PetscReal         *ca,*cb,*cg;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  pep = ctx->pep;
  nmat = pep->nmat;
  deg = nmat-1;
  m = pep->nloc;
  ca = pep->pbc;
  cb = pep->pbc+nmat;
  cg = pep->pbc+2*nmat;
  x1=ctx->w[0];y1=ctx->w[1];y2=ctx->w[2];y3=ctx->w[3];aux=ctx->w[4];aux2=ctx->w[5];
  PetscCall(EPSGetTarget(ctx->eps,&sigma));
  PetscCall(VecSet(y,0.0));
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  a = pep->sfactor;

  /* first block */
  PetscCall(VecPlaceArray(x1,px));
  PetscCall(VecPlaceArray(y1,py+m));
  PetscCall(VecCopy(x1,y1));
  PetscCall(VecScale(y1,1.0/ca[0]));
  PetscCall(VecResetArray(x1));
  PetscCall(VecResetArray(y1));

  /* second block */
  if (deg>2) {
    PetscCall(VecPlaceArray(x1,px+m));
    PetscCall(VecPlaceArray(y1,py+m));
    PetscCall(VecPlaceArray(y2,py+2*m));
    PetscCall(VecCopy(x1,y2));
    PetscCall(VecAXPY(y2,sigma-cb[1],y1));
    PetscCall(VecScale(y2,1.0/ca[1]));
    PetscCall(VecResetArray(x1));
    PetscCall(VecResetArray(y1));
    PetscCall(VecResetArray(y2));
  }

  /* inner blocks */
  for (i=2;i<deg-1;i++) {
    PetscCall(VecPlaceArray(x1,px+i*m));
    PetscCall(VecPlaceArray(y1,py+(i-1)*m));
    PetscCall(VecPlaceArray(y2,py+i*m));
    PetscCall(VecPlaceArray(y3,py+(i+1)*m));
    PetscCall(VecCopy(x1,y3));
    PetscCall(VecAXPY(y3,sigma-cb[i],y2));
    PetscCall(VecAXPY(y3,-cg[i],y1));
    PetscCall(VecScale(y3,1.0/ca[i]));
    PetscCall(VecResetArray(x1));
    PetscCall(VecResetArray(y1));
    PetscCall(VecResetArray(y2));
    PetscCall(VecResetArray(y3));
  }

  /* last block */
  PetscCall(VecPlaceArray(y1,py));
  for (i=0;i<deg-2;i++) {
    PetscCall(VecPlaceArray(y2,py+(i+1)*m));
    PetscCall(STMatMult(pep->st,i+1,y2,aux));
    PetscCall(VecAXPY(y1,a,aux));
    PetscCall(VecResetArray(y2));
    a *= pep->sfactor;
  }
  i = deg-2;
  PetscCall(VecPlaceArray(y2,py+(i+1)*m));
  PetscCall(VecPlaceArray(y3,py+i*m));
  PetscCall(VecCopy(y2,aux2));
  PetscCall(VecAXPY(aux2,cg[i+1]/ca[i+1],y3));
  PetscCall(STMatMult(pep->st,i+1,aux2,aux));
  PetscCall(VecAXPY(y1,a,aux));
  PetscCall(VecResetArray(y2));
  PetscCall(VecResetArray(y3));
  a *= pep->sfactor;
  i = deg-1;
  PetscCall(VecPlaceArray(x1,px+i*m));
  PetscCall(VecPlaceArray(y3,py+i*m));
  PetscCall(VecCopy(x1,aux2));
  PetscCall(VecAXPY(aux2,sigma-cb[i],y3));
  PetscCall(VecScale(aux2,1.0/ca[i]));
  PetscCall(STMatMult(pep->st,i+1,aux2,aux));
  PetscCall(VecAXPY(y1,a,aux));
  PetscCall(VecResetArray(x1));
  PetscCall(VecResetArray(y3));

  PetscCall(VecCopy(y1,aux));
  PetscCall(STMatSolve(pep->st,aux,y1));
  PetscCall(VecScale(y1,-1.0));

  /* final update */
  for (i=1;i<deg;i++) {
    PetscCall(VecPlaceArray(y2,py+i*m));
    tt = t;
    t = ((sigma-cb[i-1])*t-cg[i-1]*tp)/ca[i-1]; /* i-th basis polynomial */
    tp = tt;
    PetscCall(VecAXPY(y2,t,y1));
    PetscCall(VecResetArray(y2));
  }
  PetscCall(VecResetArray(y1));

  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BackTransform_Linear(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  PEP_LINEAR     *ctx;
  ST             stctx;

  PetscFunctionBegin;
  PetscCall(STShellGetContext(st,&ctx));
  PetscCall(PEPGetST(ctx->pep,&stctx));
  PetscCall(STBackTransform(stctx,n,eigr,eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Dummy backtransform operation
 */
static PetscErrorCode BackTransform_Skip(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscFunctionBegin;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode Apply_Linear(ST st,Vec x,Vec y)
{
  PEP_LINEAR     *ctx;

  PetscFunctionBegin;
  PetscCall(STShellGetContext(st,&ctx));
  PetscCall(MatMult(ctx->A,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetUp_Linear(PEP pep)
{
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;
  ST             st;
  PetscInt       i=0,deg=pep->nmat-1;
  EPSWhich       which = EPS_LARGEST_MAGNITUDE;
  EPSProblemType ptype;
  PetscBool      trackall,istrivial,transf,sinv,ks;
  PetscScalar    sigma,*epsarray,*peparray;
  Vec            veps,w=NULL;
  /* function tables */
  PetscErrorCode (*fcreate[][2])(MPI_Comm,PEP_LINEAR*,Mat*) = {
    { MatCreateExplicit_Linear_NA, MatCreateExplicit_Linear_NB },
    { MatCreateExplicit_Linear_SA, MatCreateExplicit_Linear_SB },
    { MatCreateExplicit_Linear_HA, MatCreateExplicit_Linear_HB },
  };

  PetscFunctionBegin;
  PEPCheckShiftSinvert(pep);
  PEPCheckUnsupported(pep,PEP_FEATURE_STOPPING);
  PEPCheckIgnored(pep,PEP_FEATURE_CONVERGENCE);
  PetscCall(STGetTransform(pep->st,&transf));
  PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv));
  if (!pep->which) PetscCall(PEPSetWhichEigenpairs_Default(pep));
  PetscCheck(pep->which!=PEP_ALL,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues");
  PetscCall(STSetUp(pep->st));
  if (!ctx->eps) PetscCall(PEPLinearGetEPS(pep,&ctx->eps));
  PetscCall(EPSGetST(ctx->eps,&st));
  if (!transf && !ctx->usereps) PetscCall(EPSSetTarget(ctx->eps,pep->target));
  if (sinv && !transf && !ctx->usereps) PetscCall(STSetDefaultShift(st,pep->target));
  /* compute scale factor if not set by user */
  PetscCall(PEPComputeScaleFactor(pep));

  if (ctx->explicitmatrix) {
    PEPCheckQuadraticCondition(pep,PETSC_TRUE," (with explicit matrix)");
    PEPCheckUnsupportedCondition(pep,PEP_FEATURE_NONMONOMIAL,PETSC_TRUE," (with explicit matrix)");
    PetscCheck(!transf,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Explicit matrix option is not implemented with st-transform flag active");
    PetscCheck(pep->scale!=PEP_SCALE_DIAGONAL && pep->scale!=PEP_SCALE_BOTH,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Diagonal scaling not allowed in PEPLINEAR with explicit matrices");
    if (sinv && !transf) PetscCall(STSetType(st,STSINVERT));
    PetscCall(RGPushScale(pep->rg,1.0/pep->sfactor));
    PetscCall(STGetMatrixTransformed(pep->st,0,&ctx->K));
    PetscCall(STGetMatrixTransformed(pep->st,1,&ctx->C));
    PetscCall(STGetMatrixTransformed(pep->st,2,&ctx->M));
    ctx->sfactor = pep->sfactor;
    ctx->dsfactor = pep->dsfactor;

    PetscCall(MatDestroy(&ctx->A));
    PetscCall(MatDestroy(&ctx->B));
    PetscCall(VecDestroy(&ctx->w[0]));
    PetscCall(VecDestroy(&ctx->w[1]));
    PetscCall(VecDestroy(&ctx->w[2]));
    PetscCall(VecDestroy(&ctx->w[3]));

    switch (pep->problem_type) {
      case PEP_GENERAL:    i = 0; break;
      case PEP_HERMITIAN:
      case PEP_HYPERBOLIC: i = 1; break;
      case PEP_GYROSCOPIC: i = 2; break;
    }

    PetscCall((*fcreate[i][0])(PetscObjectComm((PetscObject)pep),ctx,&ctx->A));
    PetscCall((*fcreate[i][1])(PetscObjectComm((PetscObject)pep),ctx,&ctx->B));

  } else {   /* implicit matrix */
    PetscCheck(pep->problem_type==PEP_GENERAL,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Must use the explicit matrix option if problem type is not general");
    if (!((PetscObject)ctx->eps)->type_name) PetscCall(EPSSetType(ctx->eps,EPSKRYLOVSCHUR));
    else {
      PetscCall(PetscObjectTypeCompare((PetscObject)ctx->eps,EPSKRYLOVSCHUR,&ks));
      PetscCheck(ks,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Implicit matrix option only implemented for Krylov-Schur");
    }
    PetscCheck(ctx->alpha==1.0 && ctx->beta==0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Implicit matrix option does not support setting alpha,beta parameters of the linearization");
    PetscCall(STSetType(st,STSHELL));
    PetscCall(STShellSetContext(st,ctx));
    if (!transf) PetscCall(STShellSetBackTransform(st,BackTransform_Linear));
    else PetscCall(STShellSetBackTransform(st,BackTransform_Skip));
    PetscCall(MatCreateVecsEmpty(pep->A[0],&ctx->w[0],&ctx->w[1]));
    PetscCall(MatCreateVecsEmpty(pep->A[0],&ctx->w[2],&ctx->w[3]));
    PetscCall(MatCreateVecs(pep->A[0],&ctx->w[4],&ctx->w[5]));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)pep),deg*pep->nloc,deg*pep->nloc,deg*pep->n,deg*pep->n,ctx,&ctx->A));
    if (sinv && !transf) PetscCall(MatShellSetOperation(ctx->A,MATOP_MULT,(void(*)(void))MatMult_Linear_Sinvert));
    else PetscCall(MatShellSetOperation(ctx->A,MATOP_MULT,(void(*)(void))MatMult_Linear_Shift));
    PetscCall(STShellSetApply(st,Apply_Linear));
    ctx->pep = pep;

    PetscCall(PEPBasisCoefficients(pep,pep->pbc));
    if (!transf) {
      PetscCall(PetscMalloc1(pep->nmat,&pep->solvematcoeffs));
      if (sinv) PetscCall(PEPEvaluateBasis(pep,pep->target,0,pep->solvematcoeffs,NULL));
      else {
        for (i=0;i<deg;i++) pep->solvematcoeffs[i] = 0.0;
        pep->solvematcoeffs[deg] = 1.0;
      }
      PetscCall(STScaleShift(pep->st,1.0/pep->sfactor));
      PetscCall(RGPushScale(pep->rg,1.0/pep->sfactor));
    }
    if (pep->sfactor!=1.0) {
      for (i=0;i<pep->nmat;i++) {
        pep->pbc[pep->nmat+i] /= pep->sfactor;
        pep->pbc[2*pep->nmat+i] /= pep->sfactor*pep->sfactor;
      }
    }
  }

  PetscCall(EPSSetOperators(ctx->eps,ctx->A,ctx->B));
  PetscCall(EPSGetProblemType(ctx->eps,&ptype));
  if (!ptype) {
    if (ctx->explicitmatrix) PetscCall(EPSSetProblemType(ctx->eps,EPS_GNHEP));
    else PetscCall(EPSSetProblemType(ctx->eps,EPS_NHEP));
  }
  if (!ctx->usereps) {
    if (transf) which = EPS_LARGEST_MAGNITUDE;
    else {
      switch (pep->which) {
        case PEP_LARGEST_MAGNITUDE:  which = EPS_LARGEST_MAGNITUDE; break;
        case PEP_SMALLEST_MAGNITUDE: which = EPS_SMALLEST_MAGNITUDE; break;
        case PEP_LARGEST_REAL:       which = EPS_LARGEST_REAL; break;
        case PEP_SMALLEST_REAL:      which = EPS_SMALLEST_REAL; break;
        case PEP_LARGEST_IMAGINARY:  which = EPS_LARGEST_IMAGINARY; break;
        case PEP_SMALLEST_IMAGINARY: which = EPS_SMALLEST_IMAGINARY; break;
        case PEP_TARGET_MAGNITUDE:   which = EPS_TARGET_MAGNITUDE; break;
        case PEP_TARGET_REAL:        which = EPS_TARGET_REAL; break;
        case PEP_TARGET_IMAGINARY:   which = EPS_TARGET_IMAGINARY; break;
        case PEP_ALL:                which = EPS_ALL; break;
        case PEP_WHICH_USER:         which = EPS_WHICH_USER;
          PetscCall(EPSSetEigenvalueComparison(ctx->eps,pep->sc->comparison,pep->sc->comparisonctx));
          break;
      }
    }
    PetscCall(EPSSetWhichEigenpairs(ctx->eps,which));

    PetscCall(EPSSetDimensions(ctx->eps,pep->nev,pep->ncv,pep->mpd));
    PetscCall(EPSSetTolerances(ctx->eps,SlepcDefaultTol(pep->tol),pep->max_it));
  }
  PetscCall(RGIsTrivial(pep->rg,&istrivial));
  if (!istrivial) {
    PetscCheck(!transf,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"PEPLINEAR does not support a nontrivial region with st-transform");
    PetscCall(EPSSetRG(ctx->eps,pep->rg));
  }
  /* Transfer the trackall option from pep to eps */
  PetscCall(PEPGetTrackAll(pep,&trackall));
  PetscCall(EPSSetTrackAll(ctx->eps,trackall));

  /* temporary change of target */
  if (pep->sfactor!=1.0) {
    PetscCall(EPSGetTarget(ctx->eps,&sigma));
    PetscCall(EPSSetTarget(ctx->eps,sigma/pep->sfactor));
  }

  /* process initial vector */
  if (pep->nini<0) {
    PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)ctx->eps),deg*pep->nloc,deg*pep->n,&veps));
    PetscCall(VecGetArray(veps,&epsarray));
    for (i=0;i<deg;i++) {
      if (i<-pep->nini) {
        PetscCall(VecGetArray(pep->IS[i],&peparray));
        PetscCall(PetscArraycpy(epsarray+i*pep->nloc,peparray,pep->nloc));
        PetscCall(VecRestoreArray(pep->IS[i],&peparray));
      } else {
        if (!w) PetscCall(VecDuplicate(pep->IS[0],&w));
        PetscCall(VecSetRandom(w,NULL));
        PetscCall(VecGetArray(w,&peparray));
        PetscCall(PetscArraycpy(epsarray+i*pep->nloc,peparray,pep->nloc));
        PetscCall(VecRestoreArray(w,&peparray));
      }
    }
    PetscCall(VecRestoreArray(veps,&epsarray));
    PetscCall(EPSSetInitialSpace(ctx->eps,1,&veps));
    PetscCall(VecDestroy(&veps));
    PetscCall(VecDestroy(&w));
    PetscCall(SlepcBasisDestroy_Private(&pep->nini,&pep->IS));
  }

  PetscCall(EPSSetUp(ctx->eps));
  PetscCall(EPSGetDimensions(ctx->eps,NULL,&pep->ncv,&pep->mpd));
  PetscCall(EPSGetTolerances(ctx->eps,NULL,&pep->max_it));
  PetscCall(PEPAllocateSolution(pep,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPLinearExtract_Residual - Auxiliary routine that copies the solution of the
   linear eigenproblem to the PEP object. The eigenvector of the generalized
   problem is supposed to be
                               z = [  x  ]
                                   [ l*x ]
   The eigenvector is taken from z(1:n) or z(n+1:2*n) depending on the explicitly
   computed residual norm.
   Finally, x is normalized so that ||x||_2 = 1.
*/
static PetscErrorCode PEPLinearExtract_Residual(PEP pep,EPS eps)
{
  PetscInt          i,k;
  const PetscScalar *px;
  PetscScalar       *er=pep->eigr,*ei=pep->eigi;
  PetscReal         rn1,rn2;
  Vec               xr,xi=NULL,wr;
  Mat               A;
#if !defined(PETSC_USE_COMPLEX)
  Vec               wi;
  const PetscScalar *py;
#endif

  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PEPSetWorkVecs(pep,2));
#else
  PetscCall(PEPSetWorkVecs(pep,4));
#endif
  PetscCall(EPSGetOperators(eps,&A,NULL));
  PetscCall(MatCreateVecs(A,&xr,NULL));
  PetscCall(MatCreateVecsEmpty(pep->A[0],&wr,NULL));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDuplicate(xr,&xi));
  PetscCall(VecDuplicateEmpty(wr,&wi));
#endif
  for (i=0;i<pep->nconv;i++) {
    PetscCall(EPSGetEigenpair(eps,i,NULL,NULL,xr,xi));
#if !defined(PETSC_USE_COMPLEX)
    if (ei[i]!=0.0) {   /* complex conjugate pair */
      PetscCall(VecGetArrayRead(xr,&px));
      PetscCall(VecGetArrayRead(xi,&py));
      PetscCall(VecPlaceArray(wr,px));
      PetscCall(VecPlaceArray(wi,py));
      PetscCall(VecNormalizeComplex(wr,wi,PETSC_TRUE,NULL));
      PetscCall(PEPComputeResidualNorm_Private(pep,er[i],ei[i],wr,wi,pep->work,&rn1));
      PetscCall(BVInsertVec(pep->V,i,wr));
      PetscCall(BVInsertVec(pep->V,i+1,wi));
      for (k=1;k<pep->nmat-1;k++) {
        PetscCall(VecResetArray(wr));
        PetscCall(VecResetArray(wi));
        PetscCall(VecPlaceArray(wr,px+k*pep->nloc));
        PetscCall(VecPlaceArray(wi,py+k*pep->nloc));
        PetscCall(VecNormalizeComplex(wr,wi,PETSC_TRUE,NULL));
        PetscCall(PEPComputeResidualNorm_Private(pep,er[i],ei[i],wr,wi,pep->work,&rn2));
        if (rn1>rn2) {
          PetscCall(BVInsertVec(pep->V,i,wr));
          PetscCall(BVInsertVec(pep->V,i+1,wi));
          rn1 = rn2;
        }
      }
      PetscCall(VecResetArray(wr));
      PetscCall(VecResetArray(wi));
      PetscCall(VecRestoreArrayRead(xr,&px));
      PetscCall(VecRestoreArrayRead(xi,&py));
      i++;
    } else   /* real eigenvalue */
#endif
    {
      PetscCall(VecGetArrayRead(xr,&px));
      PetscCall(VecPlaceArray(wr,px));
      PetscCall(VecNormalizeComplex(wr,NULL,PETSC_FALSE,NULL));
      PetscCall(PEPComputeResidualNorm_Private(pep,er[i],ei[i],wr,NULL,pep->work,&rn1));
      PetscCall(BVInsertVec(pep->V,i,wr));
      for (k=1;k<pep->nmat-1;k++) {
        PetscCall(VecResetArray(wr));
        PetscCall(VecPlaceArray(wr,px+k*pep->nloc));
        PetscCall(VecNormalizeComplex(wr,NULL,PETSC_FALSE,NULL));
        PetscCall(PEPComputeResidualNorm_Private(pep,er[i],ei[i],wr,NULL,pep->work,&rn2));
        if (rn1>rn2) {
          PetscCall(BVInsertVec(pep->V,i,wr));
          rn1 = rn2;
        }
      }
      PetscCall(VecResetArray(wr));
      PetscCall(VecRestoreArrayRead(xr,&px));
    }
  }
  PetscCall(VecDestroy(&wr));
  PetscCall(VecDestroy(&xr));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDestroy(&wi));
  PetscCall(VecDestroy(&xi));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPLinearExtract_None - Same as PEPLinearExtract_Norm but always takes
   the first block.
*/
static PetscErrorCode PEPLinearExtract_None(PEP pep,EPS eps)
{
  PetscInt          i;
  const PetscScalar *px;
  Mat               A;
  Vec               xr,xi=NULL,w;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar       *ei=pep->eigi;
#endif

  PetscFunctionBegin;
  PetscCall(EPSGetOperators(eps,&A,NULL));
  PetscCall(MatCreateVecs(A,&xr,NULL));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDuplicate(xr,&xi));
#endif
  PetscCall(MatCreateVecsEmpty(pep->A[0],&w,NULL));
  for (i=0;i<pep->nconv;i++) {
    PetscCall(EPSGetEigenvector(eps,i,xr,xi));
#if !defined(PETSC_USE_COMPLEX)
    if (ei[i]!=0.0) {   /* complex conjugate pair */
      PetscCall(VecGetArrayRead(xr,&px));
      PetscCall(VecPlaceArray(w,px));
      PetscCall(BVInsertVec(pep->V,i,w));
      PetscCall(VecResetArray(w));
      PetscCall(VecRestoreArrayRead(xr,&px));
      PetscCall(VecGetArrayRead(xi,&px));
      PetscCall(VecPlaceArray(w,px));
      PetscCall(BVInsertVec(pep->V,i+1,w));
      PetscCall(VecResetArray(w));
      PetscCall(VecRestoreArrayRead(xi,&px));
      i++;
    } else   /* real eigenvalue */
#endif
    {
      PetscCall(VecGetArrayRead(xr,&px));
      PetscCall(VecPlaceArray(w,px));
      PetscCall(BVInsertVec(pep->V,i,w));
      PetscCall(VecResetArray(w));
      PetscCall(VecRestoreArrayRead(xr,&px));
    }
  }
  PetscCall(VecDestroy(&w));
  PetscCall(VecDestroy(&xr));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDestroy(&xi));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPLinearExtract_Norm - Auxiliary routine that copies the solution of the
   linear eigenproblem to the PEP object. The eigenvector of the generalized
   problem is supposed to be
                               z = [  x  ]
                                   [ l*x ]
   If |l|<1.0, the eigenvector is taken from z(1:n), otherwise from z(n+1:2*n).
   Finally, x is normalized so that ||x||_2 = 1.
*/
static PetscErrorCode PEPLinearExtract_Norm(PEP pep,EPS eps)
{
  PetscInt          i,offset;
  const PetscScalar *px;
  PetscScalar       *er=pep->eigr;
  Mat               A;
  Vec               xr,xi=NULL,w;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar       *ei=pep->eigi;
#endif

  PetscFunctionBegin;
  PetscCall(EPSGetOperators(eps,&A,NULL));
  PetscCall(MatCreateVecs(A,&xr,NULL));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDuplicate(xr,&xi));
#endif
  PetscCall(MatCreateVecsEmpty(pep->A[0],&w,NULL));
  for (i=0;i<pep->nconv;i++) {
    PetscCall(EPSGetEigenpair(eps,i,NULL,NULL,xr,xi));
    if (SlepcAbsEigenvalue(er[i],ei[i])>1.0) offset = (pep->nmat-2)*pep->nloc;
    else offset = 0;
#if !defined(PETSC_USE_COMPLEX)
    if (ei[i]!=0.0) {   /* complex conjugate pair */
      PetscCall(VecGetArrayRead(xr,&px));
      PetscCall(VecPlaceArray(w,px+offset));
      PetscCall(BVInsertVec(pep->V,i,w));
      PetscCall(VecResetArray(w));
      PetscCall(VecRestoreArrayRead(xr,&px));
      PetscCall(VecGetArrayRead(xi,&px));
      PetscCall(VecPlaceArray(w,px+offset));
      PetscCall(BVInsertVec(pep->V,i+1,w));
      PetscCall(VecResetArray(w));
      PetscCall(VecRestoreArrayRead(xi,&px));
      i++;
    } else /* real eigenvalue */
#endif
    {
      PetscCall(VecGetArrayRead(xr,&px));
      PetscCall(VecPlaceArray(w,px+offset));
      PetscCall(BVInsertVec(pep->V,i,w));
      PetscCall(VecResetArray(w));
      PetscCall(VecRestoreArrayRead(xr,&px));
    }
  }
  PetscCall(VecDestroy(&w));
  PetscCall(VecDestroy(&xr));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDestroy(&xi));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPExtractVectors_Linear(PEP pep)
{
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  switch (pep->extract) {
  case PEP_EXTRACT_NONE:
    PetscCall(PEPLinearExtract_None(pep,ctx->eps));
    break;
  case PEP_EXTRACT_NORM:
    PetscCall(PEPLinearExtract_Norm(pep,ctx->eps));
    break;
  case PEP_EXTRACT_RESIDUAL:
    PetscCall(PEPLinearExtract_Residual(pep,ctx->eps));
    break;
  case PEP_EXTRACT_STRUCTURED:
    SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Extraction not implemented in this solver");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSolve_Linear(PEP pep)
{
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;
  PetscScalar    sigma;
  PetscBool      flg;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCall(EPSSolve(ctx->eps));
  PetscCall(EPSGetConverged(ctx->eps,&pep->nconv));
  PetscCall(EPSGetIterationNumber(ctx->eps,&pep->its));
  PetscCall(EPSGetConvergedReason(ctx->eps,(EPSConvergedReason*)&pep->reason));

  /* recover eigenvalues */
  for (i=0;i<pep->nconv;i++) {
    PetscCall(EPSGetEigenpair(ctx->eps,i,&pep->eigr[i],&pep->eigi[i],NULL,NULL));
    pep->eigr[i] *= pep->sfactor;
    pep->eigi[i] *= pep->sfactor;
  }

  /* restore target */
  PetscCall(EPSGetTarget(ctx->eps,&sigma));
  PetscCall(EPSSetTarget(ctx->eps,sigma*pep->sfactor));

  PetscCall(STGetTransform(pep->st,&flg));
  if (flg) PetscTryTypeMethod(pep,backtransform);
  if (pep->sfactor!=1.0) {
    /* Restore original values */
    for (i=0;i<pep->nmat;i++) {
      pep->pbc[pep->nmat+i] *= pep->sfactor;
      pep->pbc[2*pep->nmat+i] *= pep->sfactor*pep->sfactor;
    }
    if (!flg && !ctx->explicitmatrix) PetscCall(STScaleShift(pep->st,pep->sfactor));
  }
  if (ctx->explicitmatrix || !flg) PetscCall(RGPopScale(pep->rg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSMonitor_Linear(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PEP            pep = (PEP)ctx;

  PetscFunctionBegin;
  PetscCall(PEPMonitor(pep,its,nconv,eigr,eigi,errest,nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetFromOptions_Linear(PEP pep,PetscOptionItems PetscOptionsObject)
{
  PetscBool      set,val;
  PetscInt       k;
  PetscReal      array[2]={0,0};
  PetscBool      flg;
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"PEP Linear Options");

    k = 2;
    PetscCall(PetscOptionsRealArray("-pep_linear_linearization","Parameters of the linearization","PEPLinearSetLinearization",array,&k,&flg));
    if (flg) PetscCall(PEPLinearSetLinearization(pep,array[0],array[1]));

    PetscCall(PetscOptionsBool("-pep_linear_explicitmatrix","Use explicit matrix in linearization","PEPLinearSetExplicitMatrix",ctx->explicitmatrix,&val,&set));
    if (set) PetscCall(PEPLinearSetExplicitMatrix(pep,val));

  PetscOptionsHeadEnd();

  if (!ctx->eps) PetscCall(PEPLinearGetEPS(pep,&ctx->eps));
  PetscCall(EPSSetFromOptions(ctx->eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPLinearSetLinearization_Linear(PEP pep,PetscReal alpha,PetscReal beta)
{
  PEP_LINEAR *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  PetscCheck(beta!=0.0 || alpha!=0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"Parameters alpha and beta cannot be zero simultaneously");
  ctx->alpha = alpha;
  ctx->beta  = beta;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPLinearSetLinearization - Set the coefficients that define
   the linearization of a quadratic eigenproblem.

   Logically Collective

   Input Parameters:
+  pep   - polynomial eigenvalue solver
.  alpha - first parameter of the linearization
-  beta  - second parameter of the linearization

   Options Database Key:
.  -pep_linear_linearization <alpha,beta> - Sets the coefficients

   Notes:
   Cannot pass zero for both alpha and beta. The default values are
   alpha=1 and beta=0.

   Level: advanced

.seealso: PEPLinearGetLinearization()
@*/
PetscErrorCode PEPLinearSetLinearization(PEP pep,PetscReal alpha,PetscReal beta)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,alpha,2);
  PetscValidLogicalCollectiveReal(pep,beta,3);
  PetscTryMethod(pep,"PEPLinearSetLinearization_C",(PEP,PetscReal,PetscReal),(pep,alpha,beta));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPLinearGetLinearization_Linear(PEP pep,PetscReal *alpha,PetscReal *beta)
{
  PEP_LINEAR *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  if (alpha) *alpha = ctx->alpha;
  if (beta)  *beta  = ctx->beta;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPLinearGetLinearization - Returns the coefficients that define
   the linearization of a quadratic eigenproblem.

   Not Collective

   Input Parameter:
.  pep  - polynomial eigenvalue solver

   Output Parameters:
+  alpha - the first parameter of the linearization
-  beta  - the second parameter of the linearization

   Level: advanced

.seealso: PEPLinearSetLinearization()
@*/
PetscErrorCode PEPLinearGetLinearization(PEP pep,PetscReal *alpha,PetscReal *beta)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscUseMethod(pep,"PEPLinearGetLinearization_C",(PEP,PetscReal*,PetscReal*),(pep,alpha,beta));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPLinearSetExplicitMatrix_Linear(PEP pep,PetscBool explicitmatrix)
{
  PEP_LINEAR *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  if (ctx->explicitmatrix != explicitmatrix) {
    ctx->explicitmatrix = explicitmatrix;
    pep->state = PEP_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPLinearSetExplicitMatrix - Indicate if the matrices A and B for the
   linearization of the problem must be built explicitly.

   Logically Collective

   Input Parameters:
+  pep         - polynomial eigenvalue solver
-  explicitmat - boolean flag indicating if the matrices are built explicitly

   Options Database Key:
.  -pep_linear_explicitmatrix <boolean> - Indicates the boolean flag

   Level: advanced

.seealso: PEPLinearGetExplicitMatrix()
@*/
PetscErrorCode PEPLinearSetExplicitMatrix(PEP pep,PetscBool explicitmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,explicitmat,2);
  PetscTryMethod(pep,"PEPLinearSetExplicitMatrix_C",(PEP,PetscBool),(pep,explicitmat));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPLinearGetExplicitMatrix_Linear(PEP pep,PetscBool *explicitmat)
{
  PEP_LINEAR *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  *explicitmat = ctx->explicitmatrix;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPLinearGetExplicitMatrix - Returns the flag indicating if the matrices
   A and B for the linearization are built explicitly.

   Not Collective

   Input Parameter:
.  pep  - polynomial eigenvalue solver

   Output Parameter:
.  explicitmat - the mode flag

   Level: advanced

.seealso: PEPLinearSetExplicitMatrix()
@*/
PetscErrorCode PEPLinearGetExplicitMatrix(PEP pep,PetscBool *explicitmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(explicitmat,2);
  PetscUseMethod(pep,"PEPLinearGetExplicitMatrix_C",(PEP,PetscBool*),(pep,explicitmat));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPLinearSetEPS_Linear(PEP pep,EPS eps)
{
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)eps));
  PetscCall(EPSDestroy(&ctx->eps));
  ctx->eps     = eps;
  ctx->usereps = PETSC_TRUE;
  pep->state   = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPLinearSetEPS - Associate an eigensolver object (EPS) to the
   polynomial eigenvalue solver.

   Collective

   Input Parameters:
+  pep - polynomial eigenvalue solver
-  eps - the eigensolver object

   Level: advanced

.seealso: PEPLinearGetEPS()
@*/
PetscErrorCode PEPLinearSetEPS(PEP pep,EPS eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(pep,1,eps,2);
  PetscTryMethod(pep,"PEPLinearSetEPS_C",(PEP,EPS),(pep,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPLinearGetEPS_Linear(PEP pep,EPS *eps)
{
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  if (!ctx->eps) {
    PetscCall(EPSCreate(PetscObjectComm((PetscObject)pep),&ctx->eps));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->eps,(PetscObject)pep,1));
    PetscCall(EPSSetOptionsPrefix(ctx->eps,((PetscObject)pep)->prefix));
    PetscCall(EPSAppendOptionsPrefix(ctx->eps,"pep_linear_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->eps,((PetscObject)pep)->options));
    PetscCall(EPSMonitorSet(ctx->eps,EPSMonitor_Linear,pep,NULL));
  }
  *eps = ctx->eps;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPLinearGetEPS - Retrieve the eigensolver object (EPS) associated
   to the polynomial eigenvalue solver.

   Collective

   Input Parameter:
.  pep - polynomial eigenvalue solver

   Output Parameter:
.  eps - the eigensolver object

   Level: advanced

.seealso: PEPLinearSetEPS()
@*/
PetscErrorCode PEPLinearGetEPS(PEP pep,EPS *eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(eps,2);
  PetscUseMethod(pep,"PEPLinearGetEPS_C",(PEP,EPS*),(pep,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPView_Linear(PEP pep,PetscViewer viewer)
{
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (!ctx->eps) PetscCall(PEPLinearGetEPS(pep,&ctx->eps));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %s matrices\n",ctx->explicitmatrix? "explicit": "implicit"));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  linearization parameters: alpha=%g beta=%g\n",(double)ctx->alpha,(double)ctx->beta));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(EPSView(ctx->eps,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPReset_Linear(PEP pep)
{
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  if (!ctx->eps) PetscCall(EPSReset(ctx->eps));
  PetscCall(MatDestroy(&ctx->A));
  PetscCall(MatDestroy(&ctx->B));
  PetscCall(VecDestroy(&ctx->w[0]));
  PetscCall(VecDestroy(&ctx->w[1]));
  PetscCall(VecDestroy(&ctx->w[2]));
  PetscCall(VecDestroy(&ctx->w[3]));
  PetscCall(VecDestroy(&ctx->w[4]));
  PetscCall(VecDestroy(&ctx->w[5]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPDestroy_Linear(PEP pep)
{
  PEP_LINEAR     *ctx = (PEP_LINEAR*)pep->data;

  PetscFunctionBegin;
  PetscCall(EPSDestroy(&ctx->eps));
  PetscCall(PetscFree(pep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearSetLinearization_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearGetLinearization_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearSetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearGetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearSetExplicitMatrix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearGetExplicitMatrix_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_Linear(PEP pep)
{
  PEP_LINEAR     *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  pep->data = (void*)ctx;

  pep->lineariz       = PETSC_TRUE;
  ctx->explicitmatrix = PETSC_FALSE;
  ctx->alpha          = 1.0;
  ctx->beta           = 0.0;

  pep->ops->solve          = PEPSolve_Linear;
  pep->ops->setup          = PEPSetUp_Linear;
  pep->ops->setfromoptions = PEPSetFromOptions_Linear;
  pep->ops->destroy        = PEPDestroy_Linear;
  pep->ops->reset          = PEPReset_Linear;
  pep->ops->view           = PEPView_Linear;
  pep->ops->backtransform  = PEPBackTransform_Default;
  pep->ops->computevectors = PEPComputeVectors_Default;
  pep->ops->extractvectors = PEPExtractVectors_Linear;

  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearSetLinearization_C",PEPLinearSetLinearization_Linear));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearGetLinearization_C",PEPLinearGetLinearization_Linear));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearSetEPS_C",PEPLinearSetEPS_Linear));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearGetEPS_C",PEPLinearGetEPS_Linear));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearSetExplicitMatrix_C",PEPLinearSetExplicitMatrix_Linear));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPLinearGetExplicitMatrix_C",PEPLinearGetExplicitMatrix_Linear));
  PetscFunctionReturn(PETSC_SUCCESS);
}
