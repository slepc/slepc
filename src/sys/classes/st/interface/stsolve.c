/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   ST interface routines, callable by users
*/

#include <slepc/private/stimpl.h>            /*I "slepcst.h" I*/

PetscErrorCode STApply_Generic(ST st,Vec x,Vec y)
{
  PetscFunctionBegin;
  if (st->M && st->P) {
    PetscCall(MatMult(st->M,x,st->work[0]));
    PetscCall(STMatSolve(st,st->work[0],y));
  } else if (st->M) PetscCall(MatMult(st->M,x,y));
  else PetscCall(STMatSolve(st,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STApply - Applies the spectral transformation operator to a vector, for
   instance $y=(A - \sigma B)^{-1} Bx$ in the case of the shift-and-invert transformation
   and generalized eigenproblem.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector

   Level: developer

.seealso: [](ch:st), `STApplyTranspose()`, `STApplyHermitianTranspose()`, `STApplyMat()`
@*/
PetscErrorCode STApply(ST st,Vec x,Vec y)
{
  Mat            Op;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  PetscCheck(x!=y,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  PetscCall(VecSetErrorIfLocked(y,3));
  PetscCall(STGetOperator_Private(st,&Op));
  PetscCall(MatMult(Op,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STApplyMat_Generic(ST st,Mat B,Mat C)
{
  Mat            work;

  PetscFunctionBegin;
  if (st->M && st->P) {
    PetscCall(MatMatMult(st->M,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&work));
    PetscCall(STMatMatSolve(st,work,C));
    PetscCall(MatDestroy(&work));
  } else if (st->M) PetscCall(MatMatMult(st->M,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&C));
  else PetscCall(STMatMatSolve(st,B,C));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STApplyMat - Applies the spectral transformation operator to a matrix, for
   instance $Y = (A - \sigma B)^{-1} B X$ in the case of the shift-and-invert
   transformation and generalized eigenproblem.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  X  - input matrix

   Output Parameter:
.  Y - output matrix

   Level: developer

.seealso: [](ch:st), `STApply()`
@*/
PetscErrorCode STApplyMat(ST st,Mat X,Mat Y)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(X,MAT_CLASSID,2);
  PetscValidHeaderSpecific(Y,MAT_CLASSID,3);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  PetscCheck(X!=Y,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"X and Y must be different matrices");
  PetscUseTypeMethod(st,applymat,X,Y);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STApplyTranspose_Generic(ST st,Vec x,Vec y)
{
  PetscFunctionBegin;
  if (st->M && st->P) {
    PetscCall(STMatSolveTranspose(st,x,st->work[0]));
    PetscCall(MatMultTranspose(st->M,st->work[0],y));
  } else if (st->M) PetscCall(MatMultTranspose(st->M,x,y));
  else PetscCall(STMatSolveTranspose(st,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STApplyTranspose - Applies the transpose of the operator to a vector, for
   instance $y = B^T(A - \sigma B)^{-T} x$ in the case of the shift-and-invert
   transformation and generalized eigenproblem.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector

   Level: developer

.seealso: [](ch:st), `STApply()`, `STApplyHermitianTranspose()`
@*/
PetscErrorCode STApplyTranspose(ST st,Vec x,Vec y)
{
  Mat            Op;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  PetscCheck(x!=y,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  PetscCall(VecSetErrorIfLocked(y,3));
  PetscCall(STGetOperator_Private(st,&Op));
  PetscCall(MatMultTranspose(Op,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STApplyHermitianTranspose_Generic(ST st,Vec x,Vec y)
{
  PetscFunctionBegin;
  if (st->M && st->P) {
    PetscCall(STMatSolveHermitianTranspose(st,x,st->work[0]));
    PetscCall(MatMultHermitianTranspose(st->M,st->work[0],y));
  } else if (st->M) PetscCall(MatMultHermitianTranspose(st->M,x,y));
  else PetscCall(STMatSolveHermitianTranspose(st,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STApplyHermitianTranspose - Applies the Hermitian-transpose of the operator
   to a vector, for instance $y=B^*(A - \sigma B)^{-*}x$ in the case of the shift-and-invert
   transformation and generalized eigenproblem.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector

   Note:
   Currently implemented via `STApplyTranspose()` with appropriate conjugation.

   Level: developer

.seealso: [](ch:st), `STApply()`, `STApplyTranspose()`
@*/
PetscErrorCode STApplyHermitianTranspose(ST st,Vec x,Vec y)
{
  Mat            Op;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  PetscCheck(x!=y,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  PetscCall(VecSetErrorIfLocked(y,3));
  PetscCall(STGetOperator_Private(st,&Op));
  PetscCall(MatMultHermitianTranspose(Op,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetBilinearForm - Returns the matrix used in the bilinear form with a
   generalized eigenproblem with semi-definite $B$.

   Logically Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  B - output matrix

   Notes:
   The output matrix `B` must be destroyed after use. It will be `NULL` in
   case of standard eigenproblems.

   Level: developer

.seealso: [](ch:st), `BVSetMatrix()`
@*/
PetscErrorCode STGetBilinearForm(ST st,Mat *B)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  PetscAssertPointer(B,2);
  STCheckMatrices(st,1);
  PetscUseTypeMethod(st,getbilinearform,B);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STGetBilinearForm_Default(ST st,Mat *B)
{
  PetscFunctionBegin;
  if (st->nmat==1) *B = NULL;
  else {
    *B = st->A[1];
    PetscCall(PetscObjectReference((PetscObject)*B));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_STOperator(Mat Op,Vec x,Vec y)
{
  ST             st;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(Op,&st));
  PetscCall(STSetUp(st));
  PetscCall(PetscLogEventBegin(ST_Apply,st,x,y,0));
  if (st->D) { /* with balancing */
    PetscCall(VecPointwiseDivide(st->wb,x,st->D));
    PetscUseTypeMethod(st,apply,st->wb,y);
    PetscCall(VecPointwiseMult(y,y,st->D));
  } else PetscUseTypeMethod(st,apply,x,y);
  PetscCall(PetscLogEventEnd(ST_Apply,st,x,y,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultTranspose_STOperator(Mat Op,Vec x,Vec y)
{
  ST             st;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(Op,&st));
  PetscCall(STSetUp(st));
  PetscCall(PetscLogEventBegin(ST_ApplyTranspose,st,x,y,0));
  if (st->D) { /* with balancing */
    PetscCall(VecPointwiseMult(st->wb,x,st->D));
    PetscUseTypeMethod(st,applytrans,st->wb,y);
    PetscCall(VecPointwiseDivide(y,y,st->D));
  } else PetscUseTypeMethod(st,applytrans,x,y);
  PetscCall(PetscLogEventEnd(ST_ApplyTranspose,st,x,y,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode MatMultHermitianTranspose_STOperator(Mat Op,Vec x,Vec y)
{
  ST             st;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(Op,&st));
  PetscCall(STSetUp(st));
  PetscCall(PetscLogEventBegin(ST_ApplyHermitianTranspose,st,x,y,0));
  if (st->ops->applyhermtrans) {
    if (st->D) { /* with balancing */
      PetscCall(VecConjugate(st->D));
      PetscCall(VecPointwiseMult(st->wb,x,st->D));
      PetscUseTypeMethod(st,applyhermtrans,st->wb,y);
      PetscCall(VecPointwiseDivide(y,y,st->D));
      PetscCall(VecConjugate(st->D));
    } else PetscUseTypeMethod(st,applyhermtrans,x,y);
  } else {
    if (!st->wht) PetscCall(MatCreateVecs(st->A[0],&st->wht,NULL));
    PetscCall(VecCopy(x,st->wht));
    PetscCall(VecConjugate(st->wht));
    if (st->D) { /* with balancing */
      PetscCall(VecPointwiseMult(st->wb,st->wht,st->D));
      PetscUseTypeMethod(st,applytrans,st->wb,y);
      PetscCall(VecPointwiseDivide(y,y,st->D));
    } else PetscUseTypeMethod(st,applytrans,st->wht,y);
    PetscCall(VecConjugate(y));
  }
  PetscCall(PetscLogEventEnd(ST_ApplyHermitianTranspose,st,x,y,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode MatMatMult_STOperator(Mat Op,Mat B,Mat C,void *ctx)
{
  ST             st;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(Op,&st));
  PetscCall(STSetUp(st));
  PetscCall(PetscLogEventBegin(ST_Apply,st,B,C,0));
  PetscCall(STApplyMat_Generic(st,B,C));
  PetscCall(PetscLogEventEnd(ST_Apply,st,B,C,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STGetOperator_Private(ST st,Mat *Op)
{
  PetscInt       m,n,M,N;
  Vec            v;
  VecType        vtype;

  PetscFunctionBegin;
  if (!st->Op) {
    if (Op) *Op = NULL;
    /* create the shell matrix */
    PetscCall(MatGetLocalSize(st->A[0],&m,&n));
    PetscCall(MatGetSize(st->A[0],&M,&N));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)st),m,n,M,N,st,&st->Op));
    PetscCall(MatShellSetOperation(st->Op,MATOP_MULT,(PetscErrorCodeFn*)MatMult_STOperator));
    PetscCall(MatShellSetOperation(st->Op,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_STOperator));
#if defined(PETSC_USE_COMPLEX)
    PetscCall(MatShellSetOperation(st->Op,MATOP_MULT_HERMITIAN_TRANSPOSE,(PetscErrorCodeFn*)MatMultHermitianTranspose_STOperator));
#else
    PetscCall(MatShellSetOperation(st->Op,MATOP_MULT_HERMITIAN_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_STOperator));
#endif
    if (!st->D && st->ops->apply==STApply_Generic) {
      PetscCall(MatShellSetMatProductOperation(st->Op,MATPRODUCT_AB,NULL,MatMatMult_STOperator,NULL,MATDENSE,MATDENSE));
      PetscCall(MatShellSetMatProductOperation(st->Op,MATPRODUCT_AB,NULL,MatMatMult_STOperator,NULL,MATDENSECUDA,MATDENSECUDA));
      PetscCall(MatShellSetMatProductOperation(st->Op,MATPRODUCT_AB,NULL,MatMatMult_STOperator,NULL,MATDENSEHIP,MATDENSEHIP));
    }
    /* make sure the shell matrix generates a vector of the same type as the problem matrices */
    PetscCall(MatCreateVecs(st->A[0],&v,NULL));
    PetscCall(VecGetType(v,&vtype));
    PetscCall(MatShellSetVecType(st->Op,vtype));
    PetscCall(VecDestroy(&v));
    /* build the operator matrices */
    PetscCall(STComputeOperator(st));
  }
  if (Op) *Op = st->Op;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetOperator - Returns a shell matrix that represents the operator of the
   spectral transformation.

   Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  Op - operator matrix

   Notes:
   The operator is defined in linear eigenproblems only, not in polynomial ones,
   so the call will fail if more than 2 matrices were passed in `STSetMatrices()`.

   The returned shell matrix is essentially a wrapper to the `STApply()` and
   `STApplyTranspose()` operations. The operator can often be expressed as

      $$Op = D K^{-1} M D^{-1}$$

   where $D$ is the balancing matrix, and $M$ and $K$ are two matrices corresponding
   to the numerator and denominator for spectral transformations that represent
   a rational matrix function. In the case of `STSHELL`, the inner part $K^{-1}M$
   is replaced by the user-provided operation from `STShellSetApply()`.

   The preconditioner matrix $K$ typically depends on the value of the shift, and
   its inverse is handled via an internal `KSP` object. Normal usage does not
   require explicitly calling `STGetOperator()`, but it can be used to force the
   creation of $K$ and $M$, and then $K$ is passed to the `KSP`. This is useful for
   setting options associated with the `PCFactor` (to set MUMPS options, for instance).

   The returned matrix must NOT be destroyed by the user. Instead, when no
   longer needed it must be returned with `STRestoreOperator()`. In particular,
   this is required before modifying the `ST` matrices or the shift.

   A `NULL` pointer can be passed in `Op` in case the matrix is not required but we
   want to force its creation. In this case, `STRestoreOperator()` should not be
   called.

   Level: advanced

.seealso: [](ch:st), `STApply()`, `STApplyTranspose()`, `STSetBalanceMatrix()`, `STShellSetApply()`, `STGetKSP()`, `STSetShift()`, `STRestoreOperator()`, `STSetMatrices()`
@*/
PetscErrorCode STGetOperator(ST st,Mat *Op)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  STCheckNotSeized(st,1);
  PetscCheck(st->nmat<=2,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONGSTATE,"The operator is not defined in polynomial eigenproblems");
  PetscCall(STGetOperator_Private(st,Op));
  if (Op) st->opseized = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STRestoreOperator - Restore the previously seized operator matrix.

   Logically Collective

   Input Parameters:
+  st - the spectral transformation context
-  Op - operator matrix

   Notes:
   The arguments must match the corresponding call to `STGetOperator()`.

   Level: advanced

.seealso: [](ch:st), `STGetOperator()`
@*/
PetscErrorCode STRestoreOperator(ST st,Mat *Op)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(Op,2);
  PetscValidHeaderSpecific(*Op,MAT_CLASSID,2);
  PetscCheck(st->opseized,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONGSTATE,"Must be called after STGetOperator()");
  *Op = NULL;
  st->opseized = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   STComputeOperator - Computes the matrices that constitute the operator

      Op = D*inv(K)*M*inv(D).

   K and M are computed here (D is user-provided) from the system matrices
   and the shift sigma (whenever these are changed, this function recomputes
   K and M). This is used only in linear eigenproblems (nmat<3).

   K is the "preconditioner matrix": it is the denominator in rational operators,
   e.g. (A-sigma*B) in shift-and-invert. In non-rational transformations such
   as STFILTER, K=NULL which means identity. After computing K, it is passed to
   the internal KSP object via KSPSetOperators.

   M is the numerator in rational operators. If unused it is set to NULL (e.g.
   in STPRECOND).

   STSHELL does not compute anything here, but sets the flag as if it was ready.
*/
PetscErrorCode STComputeOperator(ST st)
{
  PC             pc;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  if (!st->opready && st->ops->computeoperator) {
    PetscCall(PetscInfo(st,"Building the operator matrices\n"));
    STCheckMatrices(st,1);
    if (!st->T) PetscCall(PetscCalloc1(PetscMax(2,st->nmat),&st->T));
    PetscCall(PetscLogEventBegin(ST_ComputeOperator,st,0,0,0));
    PetscUseTypeMethod(st,computeoperator);
    PetscCall(PetscLogEventEnd(ST_ComputeOperator,st,0,0,0));
    if (st->usesksp) {
      if (!st->ksp) PetscCall(STGetKSP(st,&st->ksp));
      if (st->P) {
        PetscCall(STSetDefaultKSP(st));
        PetscCall(ST_KSPSetOperators(st,st->P,st->Pmat?st->Pmat:st->P));
      } else {
        /* STPRECOND defaults to PCNONE if st->P is empty */
        PetscCall(KSPGetPC(st->ksp,&pc));
        PetscCall(PCSetType(pc,PCNONE));
      }
    }
  }
  st->opready = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetUp - Prepares for the use of a spectral transformation.

   Collective

   Input Parameter:
.  st - the spectral transformation context

   Level: advanced

.seealso: [](ch:st), `STCreate()`, `STApply()`, `STDestroy()`
@*/
PetscErrorCode STSetUp(ST st)
{
  PetscInt       i,n,k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  switch (st->state) {
    case ST_STATE_INITIAL:
      PetscCall(PetscInfo(st,"Setting up new ST\n"));
      if (!((PetscObject)st)->type_name) PetscCall(STSetType(st,STSHIFT));
      break;
    case ST_STATE_SETUP:
      PetscFunctionReturn(PETSC_SUCCESS);
    case ST_STATE_UPDATED:
      PetscCall(PetscInfo(st,"Setting up updated ST\n"));
      break;
  }
  PetscCall(PetscLogEventBegin(ST_SetUp,st,0,0,0));
  if (st->state!=ST_STATE_UPDATED) {
    if (!(st->nmat<3 && st->opready)) {
      if (st->T) {
        for (i=0;i<PetscMax(2,st->nmat);i++) PetscCall(MatDestroy(&st->T[i]));
      }
      PetscCall(MatDestroy(&st->P));
    }
  }
  if (st->D) {
    PetscCall(MatGetLocalSize(st->A[0],NULL,&n));
    PetscCall(VecGetLocalSize(st->D,&k));
    PetscCheck(n==k,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_SIZ,"Balance matrix has wrong dimension %" PetscInt_FMT " (should be %" PetscInt_FMT ")",k,n);
    if (!st->wb) PetscCall(VecDuplicate(st->D,&st->wb));
  }
  if (st->nmat<3 && st->transform) PetscCall(STComputeOperator(st));
  else {
    if (!st->T) PetscCall(PetscCalloc1(PetscMax(2,st->nmat),&st->T));
  }
  PetscTryTypeMethod(st,setup);
  st->state = ST_STATE_SETUP;
  PetscCall(PetscLogEventEnd(ST_SetUp,st,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Computes coefficients for the transformed polynomial,
   and stores the result in argument S.

   alpha - value of the parameter of the transformed polynomial
   beta - value of the previous shift (only used in inplace mode)
   k - index of first matrix included in the computation
   coeffs - coefficients of the expansion
   initial - true if this is the first time
   precond - whether the preconditioner matrix must be computed
*/
PetscErrorCode STMatMAXPY_Private(ST st,PetscScalar alpha,PetscScalar beta,PetscInt k,PetscScalar *coeffs,PetscBool initial,PetscBool precond,Mat *S)
{
  PetscInt       *matIdx=NULL,nmat,i,ini=-1;
  PetscScalar    t=1.0,ta,gamma;
  PetscBool      nz=PETSC_FALSE;
  Mat            *A=precond?st->Psplit:st->A;
  MatStructure   str=precond?st->strp:st->str;

  PetscFunctionBegin;
  nmat = st->nmat-k;
  switch (st->matmode) {
  case ST_MATMODE_INPLACE:
    PetscCheck(st->nmat<=2,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"ST_MATMODE_INPLACE not supported for polynomial eigenproblems");
    PetscCheck(!precond,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"ST_MATMODE_INPLACE not supported for split preconditioner");
    if (initial) {
      PetscCall(PetscObjectReference((PetscObject)A[0]));
      *S = A[0];
      gamma = alpha;
    } else gamma = alpha-beta;
    if (gamma != 0.0) {
      if (st->nmat>1) PetscCall(MatAXPY(*S,gamma,A[1],str));
      else PetscCall(MatShift(*S,gamma));
    }
    break;
  case ST_MATMODE_SHELL:
    PetscCheck(!precond,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"ST_MATMODE_SHELL not supported for split preconditioner");
    if (initial) {
      if (st->nmat>2) {
        PetscCall(PetscMalloc1(nmat,&matIdx));
        for (i=0;i<nmat;i++) matIdx[i] = k+i;
      }
      PetscCall(STMatShellCreate(st,alpha,nmat,matIdx,coeffs,S));
      if (st->nmat>2) PetscCall(PetscFree(matIdx));
    } else PetscCall(STMatShellShift(*S,alpha));
    break;
  case ST_MATMODE_COPY:
    if (coeffs) {
      for (i=0;i<nmat && ini==-1;i++) {
        if (coeffs[i]!=0.0) ini = i;
        else t *= alpha;
      }
      if (coeffs[ini] != 1.0) nz = PETSC_TRUE;
      for (i=ini+1;i<nmat&&!nz;i++) if (coeffs[i]!=0.0) nz = PETSC_TRUE;
    } else { nz = PETSC_TRUE; ini = 0; }
    if ((alpha == 0.0 || !nz) && t==1.0) {
      PetscCall(PetscObjectReference((PetscObject)A[k+ini]));
      PetscCall(MatDestroy(S));
      *S = A[k+ini];
    } else {
      if (*S && *S!=A[k+ini]) {
        PetscCall(MatSetOption(*S,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE));
        PetscCall(MatCopy(A[k+ini],*S,DIFFERENT_NONZERO_PATTERN));
      } else {
        PetscCall(MatDestroy(S));
        PetscCall(MatDuplicate(A[k+ini],MAT_COPY_VALUES,S));
        PetscCall(MatSetOption(*S,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE));
      }
      if (coeffs && coeffs[ini]!=1.0) PetscCall(MatScale(*S,coeffs[ini]));
      for (i=ini+k+1;i<PetscMax(2,st->nmat);i++) {
        t *= alpha;
        ta = t;
        if (coeffs) ta *= coeffs[i-k];
        if (ta!=0.0) {
          if (st->nmat>1) PetscCall(MatAXPY(*S,ta,A[i],str));
          else PetscCall(MatShift(*S,ta));
        }
      }
    }
  }
  PetscCall(MatSetOption(*S,MAT_SYMMETRIC,st->asymm));
  PetscCall(MatSetOption(*S,MAT_HERMITIAN,(PetscImaginaryPart(st->sigma)==0.0)?st->aherm:PETSC_FALSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Computes the values of the coefficients required by STMatMAXPY_Private
   for the case of monomial basis.
*/
PetscErrorCode STCoeffs_Monomial(ST st, PetscScalar *coeffs)
{
  PetscInt  k,i,ini,inip;

  PetscFunctionBegin;
  /* Compute binomial coefficients */
  ini = (st->nmat*(st->nmat-1))/2;
  for (i=0;i<st->nmat;i++) coeffs[ini+i]=1.0;
  for (k=st->nmat-1;k>=1;k--) {
    inip = ini+1;
    ini = (k*(k-1))/2;
    coeffs[ini] = 1.0;
    for (i=1;i<k;i++) coeffs[ini+i] = coeffs[ini+i-1]+coeffs[inip+i-1];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STPostSolve - Optional post-solve phase, intended for any actions that must
   be performed on the `ST` object after the eigensolver has finished.

   Collective

   Input Parameter:
.  st  - the spectral transformation context

   Level: developer

.seealso: [](ch:st), `EPSSolve()`
@*/
PetscErrorCode STPostSolve(ST st)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  PetscTryTypeMethod(st,postsolve);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STBackTransform - Back-transformation phase, intended for
   spectral transformations which require to transform the computed
   eigenvalues back to the original eigenvalue problem.

   Not Collective

   Input Parameters:
+  st   - the spectral transformation context
.  n    - number of eigenvalues
.  eigr - real part of a computed eigenvalues
-  eigi - imaginary part of a computed eigenvalues

   Level: developer

.seealso: [](ch:st), `STIsInjective()`
@*/
PetscErrorCode STBackTransform(ST st,PetscInt n,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  PetscTryTypeMethod(st,backtransform,n,eigr,eigi);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STIsInjective - Returns whether this spectral transformation is injective or not
   (that is, if it corresponds to a one-to-one mapping). If not, then it
   does not make sense to call `STBackTransform()`.

   Not Collective

   Input Parameter:
.  st   - the spectral transformation context

   Output Parameter:
.  is - the answer

   Note:
   In case of non-injective transformations such as `STFILTER`, the eigenvalues
   of the original eigenproblem are computed via the Rayleigh quotient
   $\rho(A,x) = \frac{x^*Ax}{x^*x}$ for each computed eigenvector $x$.

   Level: developer

.seealso: [](ch:st), `STBackTransform()`
@*/
PetscErrorCode STIsInjective(ST st,PetscBool* is)
{
  PetscBool      shell;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  PetscAssertPointer(is,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)st,STSHELL,&shell));
  if (shell) PetscCall(STIsInjective_Shell(st,is));
  else *is = st->ops->backtransform? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatSetUp - Build the preconditioner matrix used in `STMatSolve()`.

   Collective

   Input Parameters:
+  st     - the spectral transformation context
.  sigma  - the shift
-  coeffs - the coefficients (may be `NULL`)

   Note:
   This function is not intended to be called by end users, but by SLEPc
   solvers that use `ST`. It builds the internal matrix for the preconditioner as
   $$
   P=\begin{cases}
     \sum_{i=0}^{\mathtt{nmat}-1}\mathtt{coeffs[i]}\,\sigma^i A_i, & \text{if coefficients given}\\
     \sum_{i=0}^{\mathtt{nmat}-1}\sigma^i A_i, & \text{otherwise}
     \end{cases}
   $$
   then calls `KSPSetUp()`.

   Level: developer

.seealso: [](ch:st), `STMatSolve()`
@*/
PetscErrorCode STMatSetUp(ST st,PetscScalar sigma,PetscScalar coeffs[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveScalar(st,sigma,2);
  STCheckMatrices(st,1);

  PetscCall(PetscLogEventBegin(ST_MatSetUp,st,0,0,0));
  PetscCall(STMatMAXPY_Private(st,sigma,0.0,0,coeffs,PETSC_TRUE,PETSC_FALSE,&st->P));
  if (st->Psplit) PetscCall(STMatMAXPY_Private(st,sigma,0.0,0,coeffs,PETSC_TRUE,PETSC_TRUE,&st->Pmat));
  PetscCall(ST_KSPSetOperators(st,st->P,st->Pmat?st->Pmat:st->P));
  PetscCall(KSPSetUp(st->ksp));
  PetscCall(PetscLogEventEnd(ST_MatSetUp,st,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetWorkVecs - Sets a number of work vectors into the `ST` object.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  nw - number of work vectors to allocate

   Developer Notes:
   This is `SLEPC_EXTERN` because it may be required by shell `ST`s.

   Level: developer

.seealso: [](ch:st), `STMatCreateVecs()`
@*/
PetscErrorCode STSetWorkVecs(ST st,PetscInt nw)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,nw,2);
  PetscCheck(nw>0,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"nw must be > 0: nw = %" PetscInt_FMT,nw);
  if (st->nwork < nw) {
    PetscCall(VecDestroyVecs(st->nwork,&st->work));
    st->nwork = nw;
    PetscCall(PetscMalloc1(nw,&st->work));
    for (i=0;i<nw;i++) PetscCall(STMatCreateVecs(st,&st->work[i],NULL));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
