/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   A polynomial filter based on a truncated Chebyshev series.
*/

#include <slepc/private/stimpl.h>
#include "filter.h"

#define DAMPING_LANCZOS_POW 2

/*
   Gateway to Chebyshev for evaluating y=p(A)*x
   (Clenshaw with a vector)
*/
static PetscErrorCode MatMult_Chebyshev(Mat A,Vec x,Vec y)
{
  ST             st;
  ST_FILTER      *ctx;
  CHEBYSHEV_CTX  cheby;
  PetscInt       degree,i,w0=0,w1=1;
  PetscReal      s1,s2;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&st));
  ctx = (ST_FILTER*)st->data;
  cheby = (CHEBYSHEV_CTX)ctx->data;
  degree = ctx->polyDegree;
  s1 = 4/(ctx->right - ctx->left);
  s2 = -2*(ctx->right+ctx->left)/(ctx->right-ctx->left);

  PetscCall(VecSet(st->work[w0],0.0));
  PetscCall(VecCopy(x,y));
  PetscCall(VecScale(y,cheby->damping_coeffs[degree]*cheby->coeffs[degree]));

  for (i=0;i<degree;i++) {
    /* swap work vectors, save a copy of current y */
    w0 = 1-w0;
    w1 = 1-w1;
    PetscCall(VecCopy(y,st->work[w0]));
    /* compute next y */
    if (i == degree-1) {
      s1 = s1/2;
      s2 = s2/2;
    }
    PetscCall(VecAXPBYPCZ(y,cheby->damping_coeffs[degree-i-1]*cheby->coeffs[degree-i-1],-1,s2,x,st->work[w1]));
    PetscCall(MatMult(ctx->T,st->work[w0],st->work[w1]));
    PetscCall(VecAXPY(y,s1,st->work[w1]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Gateway to Chebyshev for evaluating C=p(A)*B
   (Clenshaw with a matrix)
*/
static PetscErrorCode MatMatMult_Chebyshev(Mat A,Mat B,Mat C,void *pctx)
{
  ST             st;
  ST_FILTER      *ctx;
  CHEBYSHEV_CTX  cheby;
  PetscInt       degree,i,m1,m2,w0=0,w1=1;
  PetscReal      s1,s2;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&st));
  ctx = (ST_FILTER*)st->data;
  cheby = (CHEBYSHEV_CTX)ctx->data;
  degree = ctx->polyDegree;
  s1 = 4/(ctx->right - ctx->left);
  s2 = -2*(ctx->right+ctx->left)/(ctx->right-ctx->left);

  if (ctx->nW) {  /* check if work matrices must be resized */
    PetscCall(MatGetSize(B,NULL,&m1));
    PetscCall(MatGetSize(ctx->W[0],NULL,&m2));
    if (m1!=m2) {
      PetscCall(MatDestroyMatrices(ctx->nW,&ctx->W));
      ctx->nW = 0;
    }
  }
  if (!ctx->nW) {  /* allocate work matrices */
    ctx->nW = 2;
    PetscCall(PetscMalloc1(ctx->nW,&ctx->W));
    for (i=0;i<ctx->nW;i++) PetscCall(MatDuplicate(B,MAT_DO_NOT_COPY_VALUES,&ctx->W[i]));
  }

  PetscCall(MatScale(ctx->W[w0],0));
  PetscCall(MatCopy(B,C,SAME_NONZERO_PATTERN));
  PetscCall(MatScale(C,cheby->damping_coeffs[degree]*cheby->coeffs[degree]));

  for (i=0;i<degree;i++) {
    /* swap work matrices, save a copy of current C */
    w0 = 1-w0;
    w1 = 1-w1;
    PetscCall(MatCopy(C,ctx->W[w0],SAME_NONZERO_PATTERN));
    if (i == degree-1) {
      s1 = s1/2;
      s2 = s2/2;
    }
    PetscCall(MatScale(C,s2));
    PetscCall(MatAXPY(C,-1,ctx->W[w1],SAME_NONZERO_PATTERN));
    PetscCall(MatAXPY(C,cheby->damping_coeffs[degree-i-1]*cheby->coeffs[degree-i-1],B,SAME_NONZERO_PATTERN));
    PetscCall(MatMatMult(ctx->T,ctx->W[w0],MAT_REUSE_MATRIX,PETSC_CURRENT,&ctx->W[w1]));
    PetscCall(MatAXPY(C,s1,ctx->W[w1],SAME_NONZERO_PATTERN));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode Chebyshev_Compute_Coefficients(ST st)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  CHEBYSHEV_CTX cheby = (CHEBYSHEV_CTX)ctx->data;
  PetscInt      i,degree;
  PetscReal     s_inta,s_intb,acosa,acosb,t1,t2,s;

  PetscFunctionBegin;
  degree = ctx->polyDegree;
  PetscCall(PetscFree(cheby->coeffs));
  PetscCall(PetscMalloc1(degree+1,&cheby->coeffs));
  s_inta = PetscMax((2.0*ctx->inta - (ctx->right + ctx->left))/(ctx->right - ctx->left),-1.0);
  s_intb = PetscMin((2.0*ctx->intb - (ctx->right + ctx->left))/(ctx->right - ctx->left),1.0);
  acosa = PetscAcosReal(s_inta);
  acosb = PetscAcosReal(s_intb);
  t1 = (acosa - acosb)/PETSC_PI;
  t2 = t1;
  cheby->coeffs[0] = t1;
  for (i=0;i<degree;i++) {
    cheby->coeffs[i+1] = 2.0*(PetscSinReal((i+1)*acosa) - PetscSinReal((i+1)*acosb)) / ((i+1)*PETSC_PI);
    t1 = t1 + cheby->coeffs[i+1]*PetscCosReal((i+1)*acosa);
    t2 = t2 + cheby->coeffs[i+1]*PetscCosReal((i+1)*acosb);
  }
  if (s_inta == -1.0) s = 0.5/PetscAbs(t1);
  else if (s_intb == 1.0) s = 0.5/PetscMin(PetscAbs(t1),PetscAbs(t2));
  else s = 0.5/PetscMin(PetscAbs(t1),PetscAbs(t2));
  for (i=0;i<degree+1;i++) cheby->coeffs[i] = s*cheby->coeffs[i];
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode Chebyshev_Compute_Damping_Coefficients(ST st)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  CHEBYSHEV_CTX cheby = (CHEBYSHEV_CTX)ctx->data;
  PetscInt      i,degree;
  PetscReal     invdeg;

  PetscFunctionBegin;
  degree = ctx->polyDegree;
  PetscCall(PetscFree(cheby->damping_coeffs));
  PetscCall(PetscMalloc1(degree+1,&cheby->damping_coeffs));
  cheby->damping_coeffs[0] = 1;
  switch (ctx->damping) {
    case ST_FILTER_DAMPING_NONE:
      for (i=1;i<=degree;i++) cheby->damping_coeffs[i] = 1;
      break;
    case ST_FILTER_DAMPING_LANCZOS:
      for (i=1;i<=degree;i++) cheby->damping_coeffs[i] = PetscPowRealInt(PetscSinReal(i*PETSC_PI/(degree+1))/(i*PETSC_PI/(degree+1)),DAMPING_LANCZOS_POW);
      break;
    case ST_FILTER_DAMPING_JACKSON:
      invdeg = 1/((PetscReal)degree+2);
      for (i=1;i<=degree;i++) cheby->damping_coeffs[i] = invdeg*((degree+2-i)*PetscCosReal(PETSC_PI*i*invdeg)+PetscSinReal(PETSC_PI*i*invdeg)/PetscTanReal(PETSC_PI*invdeg));
      break;
    case ST_FILTER_DAMPING_FEJER:
      for (i=1;i<=degree;i++) cheby->damping_coeffs[i] = 1-i/(PetscReal)(degree+1);
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONG,"Invalid Chebyshev filter damping type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Computes the coefficients and damping coefficients for the Chebyshev filter

   Creates the shifted (and scaled) matrix and the filter P(z).
   G is a shell matrix whose MatMult() and MatMatMult() applies the filter.
*/
static PetscErrorCode STComputeOperator_Filter_Chebyshev(ST st,Mat *G)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  PetscInt      n,m,N,M;

  PetscFunctionBegin;
  PetscCall(STSetWorkVecs(st,2));
  PetscCall(Chebyshev_Compute_Coefficients(st));
  PetscCall(Chebyshev_Compute_Damping_Coefficients(st));
  /* create shell matrix*/
  PetscCall(MatDestroy(&ctx->T));
  PetscCall(MatDuplicate(st->A[0],MAT_COPY_VALUES,&ctx->T));
  if (!*G) {
    PetscCall(MatGetSize(ctx->T,&N,&M));
    PetscCall(MatGetLocalSize(ctx->T,&n,&m));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)st),n,m,N,M,st,G));
    PetscCall(MatShellSetOperation(*G,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Chebyshev));
    PetscCall(MatShellSetMatProductOperation(*G,MATPRODUCT_AB,NULL,MatMatMult_Chebyshev,NULL,MATDENSE,MATDENSE));
    PetscCall(MatShellSetMatProductOperation(*G,MATPRODUCT_AB,NULL,MatMatMult_Chebyshev,NULL,MATDENSECUDA,MATDENSECUDA));
    PetscCall(MatShellSetMatProductOperation(*G,MATPRODUCT_AB,NULL,MatMatMult_Chebyshev,NULL,MATDENSEHIP,MATDENSEHIP));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode Chebyshev_Clenshaw_Damping_Real(ST st,PetscReal x,PetscReal *px)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  CHEBYSHEV_CTX cheby = (CHEBYSHEV_CTX)ctx->data;
  PetscReal     t1=0,t2=0;
  PetscInt      i,degree;

  PetscFunctionBegin;
  degree = ctx->polyDegree;
  PetscCheck(cheby->coeffs,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_NULL,"Chebyshev coefficients are not computed");
  PetscCheck(cheby->damping_coeffs,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_NULL,"Chebyshev damping coefficients are not computed");
  x = 2*x;
  for (i=degree;i>1;i=i-2) {
    t2 = cheby->damping_coeffs[i]*cheby->coeffs[i] + x*t1 - t2;
    t1 = cheby->damping_coeffs[i-1]*cheby->coeffs[i-1] + x*t2 - t1;
  }
  if (degree % 2) {
    t2 = cheby->damping_coeffs[1]*cheby->coeffs[1] + x*t1 - t2;
    t1 = t2;
  }
  *px = cheby->damping_coeffs[0]*cheby->coeffs[0] + .5*x*t1 - t2;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Computes the threshold for the Chebyshev filter

   The threshold is P(scale(ctx->inta)), where ctx->inta is the lower bound of the interval.
   It should be the same value if the upper bound was used instead.
*/
static PetscErrorCode STFilterGetThreshold_Filter_Chebyshev(ST st,PetscReal *gamma)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  PetscReal     c,e;

  PetscFunctionBegin;
  /* scale ctx->inta */
  c = (ctx->right + ctx->left) / 2.0;
  e = (ctx->right - ctx->left) / 2.0;
  PetscCall(Chebyshev_Clenshaw_Damping_Real(st,(ctx->inta-c)/e,gamma));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STDestroy_Filter_Chebyshev(ST st)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  CHEBYSHEV_CTX cheby = (CHEBYSHEV_CTX)ctx->data;

  PetscFunctionBegin;
  PetscCall(PetscFree(cheby->coeffs));
  PetscCall(PetscFree(cheby->damping_coeffs));
  PetscCall(PetscFree(cheby));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STCreate_Filter_Chebyshev(ST st)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  CHEBYSHEV_CTX cheby;

  PetscFunctionBegin;
  PetscCall(PetscNew(&cheby));
  ctx->data = (void*)cheby;

  ctx->computeoperator = STComputeOperator_Filter_Chebyshev;
  ctx->getthreshold    = STFilterGetThreshold_Filter_Chebyshev;
  ctx->destroy         = STDestroy_Filter_Chebyshev;
  PetscFunctionReturn(PETSC_SUCCESS);
}
