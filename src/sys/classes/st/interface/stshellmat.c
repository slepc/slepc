/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Subroutines that implement various operations of the matrix associated with
   the shift-and-invert technique for eigenvalue problems
*/

#include <slepc/private/stimpl.h>

typedef struct {
  PetscScalar alpha;
  PetscScalar *coeffs;
  ST          st;
  Vec         z;
  PetscInt    nmat;
  PetscInt    *matIdx;
} ST_MATSHELL;

PetscErrorCode STMatShellShift(Mat A,PetscScalar alpha)
{
  ST_MATSHELL    *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  ctx->alpha = alpha;
  PetscCall(PetscObjectStateIncrease((PetscObject)A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  For i=0:nmat-1 computes y = (sum_i (coeffs[i]*alpha^i*st->A[idx[i]]))x
  If null coeffs computes with coeffs[i]=1.0
*/
static PetscErrorCode MatMult_Shell(Mat A,Vec x,Vec y)
{
  ST_MATSHELL    *ctx;
  ST             st;
  PetscInt       i;
  PetscScalar    t=1.0,c;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  st = ctx->st;
  PetscCall(MatMult(st->A[ctx->matIdx[0]],x,y));
  if (ctx->coeffs && ctx->coeffs[0]!=1.0) PetscCall(VecScale(y,ctx->coeffs[0]));
  if (ctx->alpha!=0.0) {
    for (i=1;i<ctx->nmat;i++) {
      PetscCall(MatMult(st->A[ctx->matIdx[i]],x,ctx->z));
      t *= ctx->alpha;
      c = (ctx->coeffs)?t*ctx->coeffs[i]:t;
      PetscCall(VecAXPY(y,c,ctx->z));
    }
    if (ctx->nmat==1) PetscCall(VecAXPY(y,ctx->alpha,x)); /* y = (A + alpha*I) x */
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultTranspose_Shell(Mat A,Vec x,Vec y)
{
  ST_MATSHELL    *ctx;
  ST             st;
  PetscInt       i;
  PetscScalar    t=1.0,c;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  st = ctx->st;
  PetscCall(MatMultTranspose(st->A[ctx->matIdx[0]],x,y));
  if (ctx->coeffs && ctx->coeffs[0]!=1.0) PetscCall(VecScale(y,ctx->coeffs[0]));
  if (ctx->alpha!=0.0) {
    for (i=1;i<ctx->nmat;i++) {
      PetscCall(MatMultTranspose(st->A[ctx->matIdx[i]],x,ctx->z));
      t *= ctx->alpha;
      c = (ctx->coeffs)?t*ctx->coeffs[i]:t;
      PetscCall(VecAXPY(y,c,ctx->z));
    }
    if (ctx->nmat==1) PetscCall(VecAXPY(y,ctx->alpha,x)); /* y = (A + alpha*I) x */
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode MatMultHermitianTranspose_Shell(Mat A,Vec x,Vec y)
{
  ST_MATSHELL    *ctx;
  ST             st;
  PetscInt       i;
  PetscScalar    t=1.0,c;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  st = ctx->st;
  PetscCall(MatMultHermitianTranspose(st->A[ctx->matIdx[0]],x,y));
  if (ctx->coeffs && ctx->coeffs[0]!=1.0) PetscCall(VecScale(y,PetscConj(ctx->coeffs[0])));
  if (ctx->alpha!=0.0) {
    for (i=1;i<ctx->nmat;i++) {
      PetscCall(MatMultHermitianTranspose(st->A[ctx->matIdx[i]],x,ctx->z));
      t *= ctx->alpha;
      c = (ctx->coeffs)?t*ctx->coeffs[i]:t;
      PetscCall(VecAXPY(y,PetscConj(c),ctx->z));
    }
    if (ctx->nmat==1) PetscCall(VecAXPY(y,PetscConj(ctx->alpha),x)); /* y = (A + alpha*I) x */
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode MatGetDiagonal_Shell(Mat A,Vec diag)
{
  ST_MATSHELL    *ctx;
  ST             st;
  Vec            diagb;
  PetscInt       i;
  PetscScalar    t=1.0,c;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  st = ctx->st;
  PetscCall(MatGetDiagonal(st->A[ctx->matIdx[0]],diag));
  if (ctx->coeffs && ctx->coeffs[0]!=1.0) PetscCall(VecScale(diag,ctx->coeffs[0]));
  if (ctx->alpha!=0.0) {
    if (ctx->nmat==1) PetscCall(VecShift(diag,ctx->alpha)); /* y = (A + alpha*I) x */
    else {
      PetscCall(VecDuplicate(diag,&diagb));
      for (i=1;i<ctx->nmat;i++) {
        PetscCall(MatGetDiagonal(st->A[ctx->matIdx[i]],diagb));
        t *= ctx->alpha;
        c = (ctx->coeffs)?t*ctx->coeffs[i]:t;
        PetscCall(VecAYPX(diag,c,diagb));
      }
      PetscCall(VecDestroy(&diagb));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_Shell(Mat A)
{
  ST_MATSHELL    *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(VecDestroy(&ctx->z));
  PetscCall(PetscFree(ctx->matIdx));
  PetscCall(PetscFree(ctx->coeffs));
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STMatShellCreate(ST st,PetscScalar alpha,PetscInt nmat,PetscInt *matIdx,PetscScalar *coeffs,Mat *mat)
{
  PetscInt       n,m,N,M,i;
  PetscBool      has=PETSC_FALSE,hasA,hasB;
  ST_MATSHELL    *ctx;

  PetscFunctionBegin;
  PetscCall(MatGetSize(st->A[0],&M,&N));
  PetscCall(MatGetLocalSize(st->A[0],&m,&n));
  PetscCall(PetscNew(&ctx));
  ctx->st = st;
  ctx->alpha = alpha;
  ctx->nmat = matIdx?nmat:st->nmat;
  PetscCall(PetscMalloc1(ctx->nmat,&ctx->matIdx));
  if (matIdx) {
    for (i=0;i<ctx->nmat;i++) ctx->matIdx[i] = matIdx[i];
  } else {
    ctx->matIdx[0] = 0;
    if (ctx->nmat>1) ctx->matIdx[1] = 1;
  }
  if (coeffs) {
    PetscCall(PetscMalloc1(ctx->nmat,&ctx->coeffs));
    for (i=0;i<ctx->nmat;i++) ctx->coeffs[i] = coeffs[i];
  }
  PetscCall(MatCreateVecs(st->A[0],&ctx->z,NULL));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)st),m,n,M,N,(void*)ctx,mat));
  PetscCall(MatShellSetOperation(*mat,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Shell));
  PetscCall(MatShellSetOperation(*mat,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_Shell));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(MatShellSetOperation(*mat,MATOP_MULT_HERMITIAN_TRANSPOSE,(PetscErrorCodeFn*)MatMultHermitianTranspose_Shell));
#else
  PetscCall(MatShellSetOperation(*mat,MATOP_MULT_HERMITIAN_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_Shell));
#endif
  PetscCall(MatShellSetOperation(*mat,MATOP_DESTROY,(PetscErrorCodeFn*)MatDestroy_Shell));

  PetscCall(MatHasOperation(st->A[ctx->matIdx[0]],MATOP_GET_DIAGONAL,&hasA));
  if (st->nmat>1) {
    has = hasA;
    for (i=1;i<ctx->nmat;i++) {
      PetscCall(MatHasOperation(st->A[ctx->matIdx[i]],MATOP_GET_DIAGONAL,&hasB));
      has = (has && hasB)? PETSC_TRUE: PETSC_FALSE;
    }
  }
  if ((hasA && st->nmat==1) || has) PetscCall(MatShellSetOperation(*mat,MATOP_GET_DIAGONAL,(PetscErrorCodeFn*)MatGetDiagonal_Shell));
  PetscFunctionReturn(PETSC_SUCCESS);
}
