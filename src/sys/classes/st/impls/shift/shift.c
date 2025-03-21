/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Shift spectral transformation, applies (A + sigma I) as operator, or
   inv(B)(A + sigma B) for generalized problems
*/

#include <slepc/private/stimpl.h>
#include <petsc/private/matimpl.h>  /* Mat_MPIDense */

/*
   Special STApply() for the BSE structured matrix

       H = [ R  C; -C^H -R^T ].

   Assumes that H is a MATNEST and x,y are VECNEST.

   Only the upper part of the product y=H*x is computed.

       y1 = R*x1+C*x2

   The bottom part of the input vector x2 is computed as
   either conj(x1) or -conj(x1), where the sign is given by
   s in the context SlepcMatStruct.
   The bottom part of the output vector y2 is not referenced.
*/
static PetscErrorCode STApply_Shift_BSE(ST st,Vec x,Vec y)
{
  Mat            H,R,C;
  Vec            x1,x2,y1;
  PetscContainer container;
  SlepcMatStruct mctx;

  PetscFunctionBegin;
  H = st->T[0];
  PetscCall(PetscObjectQuery((PetscObject)H,"SlepcMatStruct",(PetscObject*)&container));
  PetscCall(PetscContainerGetPointer(container,(void**)&mctx));
  PetscCall(MatNestGetSubMat(H,0,0,&R));
  PetscCall(MatNestGetSubMat(H,0,1,&C));
  PetscCall(VecNestGetSubVec(x,0,&x1));
  PetscCall(VecNestGetSubVec(x,1,&x2));
  PetscCall(VecNestGetSubVec(y,0,&y1));

  /* x2 = +/-conj(x1) */
  PetscCall(VecCopy(x1,x2));
  PetscCall(VecConjugate(x2));
  if (mctx->s==-1.0) PetscCall(VecScale(x2,-1.0));

  PetscCall(MatMult(C,x2,y1));
  PetscCall(MatMultAdd(R,x1,y1,y1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Specialized version that avoids communication for multiplication by R.
   It needs to access internal data structures of MATMPIDENSE.
*/
static PetscErrorCode STApply_Shift_BSE_Dense(ST st,Vec x,Vec y)
{
  Mat            H,R,C;
  Vec            x1,x2,y1;
  Mat_MPIDense   *mdnR,*mdnC;
  PetscContainer container;
  SlepcMatStruct mctx;

  PetscFunctionBegin;
  H = st->T[0];
  PetscCall(PetscObjectQuery((PetscObject)H,"SlepcMatStruct",(PetscObject*)&container));
  PetscCall(PetscContainerGetPointer(container,(void**)&mctx));
  PetscCall(MatNestGetSubMat(H,0,0,&R));
  PetscCall(MatNestGetSubMat(H,0,1,&C));
  mdnR = (Mat_MPIDense*)R->data;
  mdnC = (Mat_MPIDense*)C->data;
  PetscCall(VecNestGetSubVec(x,0,&x1));
  PetscCall(VecNestGetSubVec(x,1,&x2));
  PetscCall(VecNestGetSubVec(y,0,&y1));

  /* x2 = +/-conj(x1) */
  PetscCall(VecCopy(x1,x2));
  PetscCall(VecConjugate(x2));
  if (mctx->s==-1.0) PetscCall(VecScale(x2,-1.0));

  PetscCall(MatMult(C,x2,y1));
  /* PetscCall(MatMultAdd(R,x1,y1,y1)); */
  PetscCall(VecConjugate(mdnC->lvec));
  if (mctx->s==-1.0) PetscCall(VecScale(mdnC->lvec,-1.0));
  PetscCall(PetscLogEventBegin(MAT_MultAdd,R,mdnC->lvec,y1,y1));
  PetscCall(VecLockReadPush(mdnC->lvec));
  PetscCall((*mdnR->A->ops->multadd)(mdnR->A,mdnC->lvec,y1,y1));
  PetscCall(VecLockReadPop(mdnC->lvec));
  PetscCall(PetscLogEventEnd(MAT_MultAdd,R,mdnC->lvec,y1,y1));
  PetscCall(PetscObjectStateIncrease((PetscObject)y1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STBackTransform_Shift(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscInt j;

  PetscFunctionBegin;
  for (j=0;j<n;j++) {
    eigr[j] += st->sigma;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STPostSolve_Shift(ST st)
{
  PetscFunctionBegin;
  if (st->matmode == ST_MATMODE_INPLACE) {
    if (st->nmat>1) PetscCall(MatAXPY(st->A[0],st->sigma,st->A[1],st->str));
    else PetscCall(MatShift(st->A[0],st->sigma));
    st->Astate[0] = ((PetscObject)st->A[0])->state;
    st->state   = ST_STATE_INITIAL;
    st->opready = PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Operator (shift):
               Op               P         M
   if nmat=1:  A-sI             NULL      A-sI
   if nmat=2:  B^-1 (A-sB)      B         A-sB
*/
static PetscErrorCode STComputeOperator_Shift(ST st)
{
  PetscFunctionBegin;
  st->usesksp = (st->nmat>1)? PETSC_TRUE: PETSC_FALSE;
  PetscCall(PetscObjectReference((PetscObject)st->A[1]));
  PetscCall(MatDestroy(&st->T[1]));
  st->T[1] = st->A[1];
  PetscCall(STMatMAXPY_Private(st,-st->sigma,0.0,0,NULL,PetscNot(st->state==ST_STATE_UPDATED),PETSC_FALSE,&st->T[0]));
  if (st->nmat>1) PetscCall(PetscObjectReference((PetscObject)st->T[1]));
  PetscCall(MatDestroy(&st->P));
  st->P = (st->nmat>1)? st->T[1]: NULL;
  st->M = st->T[0];
  if (st->nmat>1 && st->Psplit) {  /* build custom preconditioner from the split matrices */
    PetscCall(MatDestroy(&st->Pmat));
    PetscCall(PetscObjectReference((PetscObject)st->Psplit[1]));
    st->Pmat = st->Psplit[1];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STSetUp_Shift(ST st)
{
  PetscInt       k,nc,nmat=st->nmat;
  PetscScalar    *coeffs=NULL;
  PetscBool      denseR,denseC;
  Mat            H,R,C;

  PetscFunctionBegin;
  if (nmat>1) PetscCall(STSetWorkVecs(st,1));
  if (nmat>2) {  /* set-up matrices for polynomial eigenproblems */
    if (st->transform) {
      nc = (nmat*(nmat+1))/2;
      PetscCall(PetscMalloc1(nc,&coeffs));
      /* Compute coeffs */
      PetscCall(STCoeffs_Monomial(st,coeffs));
      /* T[n] = A_n */
      k = nmat-1;
      PetscCall(PetscObjectReference((PetscObject)st->A[k]));
      PetscCall(MatDestroy(&st->T[k]));
      st->T[k] = st->A[k];
      for (k=0;k<nmat-1;k++) PetscCall(STMatMAXPY_Private(st,nmat>2?st->sigma:-st->sigma,0.0,k,coeffs?coeffs+((nmat-k)*(nmat-k-1))/2:NULL,PetscNot(st->state==ST_STATE_UPDATED),PETSC_FALSE,&st->T[k]));
      PetscCall(PetscFree(coeffs));
      PetscCall(PetscObjectReference((PetscObject)st->T[nmat-1]));
      PetscCall(MatDestroy(&st->P));
      st->P = st->T[nmat-1];
      if (st->Psplit) {  /* build custom preconditioner from the split matrices */
        PetscCall(STMatMAXPY_Private(st,st->sigma,0.0,nmat-1,coeffs?coeffs:NULL,PETSC_TRUE,PETSC_TRUE,&st->Pmat));
      }
      PetscCall(ST_KSPSetOperators(st,st->P,st->Pmat?st->Pmat:st->P));
    } else {
      for (k=0;k<nmat;k++) {
        PetscCall(PetscObjectReference((PetscObject)st->A[k]));
        PetscCall(MatDestroy(&st->T[k]));
        st->T[k] = st->A[k];
      }
    }
  }
  if (st->P) PetscCall(KSPSetUp(st->ksp));
  if (st->structured) {
    H = st->T[0];
    PetscCall(MatNestGetSubMat(H,0,0,&R));
    PetscCall(MatNestGetSubMat(H,0,1,&C));
    PetscCall(PetscObjectTypeCompareAny((PetscObject)R,&denseR,MATMPIDENSE,MATMPIDENSECUDA,MATMPIDENSEHIP,""));
    PetscCall(PetscObjectTypeCompareAny((PetscObject)C,&denseC,MATMPIDENSE,MATMPIDENSECUDA,MATMPIDENSEHIP,""));
    st->ops->apply = (denseR && denseC)? STApply_Shift_BSE_Dense: STApply_Shift_BSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STSetShift_Shift(ST st,PetscScalar newshift)
{
  PetscInt       k,nc,nmat=PetscMax(st->nmat,2);
  PetscScalar    *coeffs=NULL;

  PetscFunctionBegin;
  if (st->transform) {
    if (st->matmode == ST_MATMODE_COPY && nmat>2) {
      nc = (nmat*(nmat+1))/2;
      PetscCall(PetscMalloc1(nc,&coeffs));
      /* Compute coeffs */
      PetscCall(STCoeffs_Monomial(st,coeffs));
    }
    for (k=0;k<nmat-1;k++) PetscCall(STMatMAXPY_Private(st,nmat>2?newshift:-newshift,nmat>2?st->sigma:-st->sigma,k,coeffs?coeffs+((nmat-k)*(nmat-k-1))/2:NULL,PETSC_FALSE,PETSC_FALSE,&st->T[k]));
    if (st->matmode == ST_MATMODE_COPY && nmat>2) PetscCall(PetscFree(coeffs));
    if (st->nmat<=2) st->M = st->T[0];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode STCreate_Shift(ST st)
{
  PetscFunctionBegin;
  st->usesksp = PETSC_TRUE;

  st->ops->apply           = STApply_Generic;
  st->ops->applytrans      = STApplyTranspose_Generic;
  st->ops->applyhermtrans  = STApplyHermitianTranspose_Generic;
  st->ops->backtransform   = STBackTransform_Shift;
  st->ops->setshift        = STSetShift_Shift;
  st->ops->getbilinearform = STGetBilinearForm_Default;
  st->ops->setup           = STSetUp_Shift;
  st->ops->computeoperator = STComputeOperator_Shift;
  st->ops->postsolve       = STPostSolve_Shift;
  st->ops->setdefaultksp   = STSetDefaultKSP_Default;
  PetscFunctionReturn(PETSC_SUCCESS);
}
