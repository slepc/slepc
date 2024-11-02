/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Implements the Cayley spectral transform
*/

#include <slepc/private/stimpl.h>          /*I "slepcst.h" I*/

typedef struct {
  PetscScalar nu;
  PetscBool   nu_set;
} ST_CAYLEY;

static PetscErrorCode MatMult_Cayley(Mat B,Vec x,Vec y)
{
  ST             st;
  ST_CAYLEY      *ctx;
  PetscScalar    nu;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&st));
  ctx = (ST_CAYLEY*)st->data;
  nu = ctx->nu;

  if (st->matmode == ST_MATMODE_INPLACE) { nu = nu + st->sigma; }

  if (st->nmat>1) {
    /* generalized eigenproblem: y = (A + tB)x */
    PetscCall(MatMult(st->A[0],x,y));
    PetscCall(MatMult(st->A[1],x,st->work[1]));
    PetscCall(VecAXPY(y,nu,st->work[1]));
  } else {
    /* standard eigenproblem: y = (A + tI)x */
    PetscCall(MatMult(st->A[0],x,y));
    PetscCall(VecAXPY(y,nu,x));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultTranspose_Cayley(Mat B,Vec x,Vec y)
{
  ST             st;
  ST_CAYLEY      *ctx;
  PetscScalar    nu;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&st));
  ctx = (ST_CAYLEY*)st->data;
  nu = ctx->nu;

  if (st->matmode == ST_MATMODE_INPLACE) { nu = nu + st->sigma; }
  nu = PetscConj(nu);

  if (st->nmat>1) {
    /* generalized eigenproblem: y = (A + tB)x */
    PetscCall(MatMultTranspose(st->A[0],x,y));
    PetscCall(MatMultTranspose(st->A[1],x,st->work[1]));
    PetscCall(VecAXPY(y,nu,st->work[1]));
  } else {
    /* standard eigenproblem: y = (A + tI)x */
    PetscCall(MatMultTranspose(st->A[0],x,y));
    PetscCall(VecAXPY(y,nu,x));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STGetBilinearForm_Cayley(ST st,Mat *B)
{
  PetscFunctionBegin;
  PetscCall(STSetUp(st));
  *B = st->T[0];
  PetscCall(PetscObjectReference((PetscObject)*B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STBackTransform_Cayley(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  ST_CAYLEY   *ctx = (ST_CAYLEY*)st->data;
  PetscInt    j;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar t,i,r;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  for (j=0;j<n;j++) {
    if (eigi[j] == 0.0) eigr[j] = (ctx->nu + eigr[j] * st->sigma) / (eigr[j] - 1.0);
    else {
      r = eigr[j];
      i = eigi[j];
      r = st->sigma * (r * r + i * i - r) + ctx->nu * (r - 1);
      i = - st->sigma * i - ctx->nu * i;
      t = i * i + r * (r - 2.0) + 1.0;
      eigr[j] = r / t;
      eigi[j] = i / t;
    }
  }
#else
  for (j=0;j<n;j++) {
    eigr[j] = (ctx->nu + eigr[j] * st->sigma) / (eigr[j] - 1.0);
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STPostSolve_Cayley(ST st)
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
   Operator (cayley):
               Op                  P         M
   if nmat=1:  (A-sI)^-1 (A+tI)    A-sI      A+tI
   if nmat=2:  (A-sB)^-1 (A+tB)    A-sB      A+tI
*/
static PetscErrorCode STComputeOperator_Cayley(ST st)
{
  PetscInt       n,m;
  ST_CAYLEY      *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  /* if the user did not set the shift, use the target value */
  if (!st->sigma_set) st->sigma = st->defsigma;

  if (!ctx->nu_set) ctx->nu = st->sigma;
  PetscCheck(ctx->nu!=0.0 || st->sigma!=0.0,PetscObjectComm((PetscObject)st),PETSC_ERR_USER_INPUT,"Values of shift and antishift cannot be zero simultaneously");
  PetscCheck(ctx->nu!=-st->sigma,PetscObjectComm((PetscObject)st),PETSC_ERR_USER_INPUT,"It is not allowed to set the antishift equal to minus the shift (the target)");

  /* T[0] = A+nu*B */
  if (st->matmode==ST_MATMODE_INPLACE) {
    PetscCall(MatGetLocalSize(st->A[0],&n,&m));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)st),n,m,PETSC_DETERMINE,PETSC_DETERMINE,st,&st->T[0]));
    PetscCall(MatShellSetOperation(st->T[0],MATOP_MULT,(void(*)(void))MatMult_Cayley));
    PetscCall(MatShellSetOperation(st->T[0],MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTranspose_Cayley));
  } else PetscCall(STMatMAXPY_Private(st,ctx->nu,0.0,0,NULL,PetscNot(st->state==ST_STATE_UPDATED),PETSC_FALSE,&st->T[0]));
  st->M = st->T[0];

  /* T[1] = A-sigma*B */
  PetscCall(STMatMAXPY_Private(st,-st->sigma,0.0,0,NULL,PetscNot(st->state==ST_STATE_UPDATED),PETSC_FALSE,&st->T[1]));
  PetscCall(PetscObjectReference((PetscObject)st->T[1]));
  PetscCall(MatDestroy(&st->P));
  st->P = st->T[1];
  if (st->Psplit) {  /* build custom preconditioner from the split matrices */
    PetscCall(STMatMAXPY_Private(st,-st->sigma,0.0,0,NULL,PETSC_TRUE,PETSC_TRUE,&st->Pmat));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STSetUp_Cayley(ST st)
{
  PetscFunctionBegin;
  PetscCheck(st->nmat<=2,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"Cayley transform cannot be used in polynomial eigenproblems");
  PetscCall(STSetWorkVecs(st,2));
  PetscCall(KSPSetUp(st->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STSetShift_Cayley(ST st,PetscScalar newshift)
{
  ST_CAYLEY      *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  PetscCheck(newshift!=0.0 || (ctx->nu_set && ctx->nu!=0.0),PetscObjectComm((PetscObject)st),PETSC_ERR_USER_INPUT,"Values of shift and antishift cannot be zero simultaneously");
  PetscCheck(ctx->nu!=-newshift,PetscObjectComm((PetscObject)st),PETSC_ERR_USER_INPUT,"It is not allowed to set the shift equal to minus the antishift");

  if (!ctx->nu_set) {
    if (st->matmode!=ST_MATMODE_INPLACE) PetscCall(STMatMAXPY_Private(st,newshift,ctx->nu,0,NULL,PETSC_FALSE,PETSC_FALSE,&st->T[0]));
    ctx->nu = newshift;
  }
  PetscCall(STMatMAXPY_Private(st,-newshift,-st->sigma,0,NULL,PETSC_FALSE,PETSC_FALSE,&st->T[1]));
  if (st->P!=st->T[1]) {
    PetscCall(PetscObjectReference((PetscObject)st->T[1]));
    PetscCall(MatDestroy(&st->P));
    st->P = st->T[1];
  }
  if (st->Psplit) {  /* build custom preconditioner from the split matrices */
    PetscCall(STMatMAXPY_Private(st,-newshift,-st->sigma,0,NULL,PETSC_FALSE,PETSC_TRUE,&st->Pmat));
  }
  PetscCall(ST_KSPSetOperators(st,st->P,st->Pmat?st->Pmat:st->P));
  PetscCall(KSPSetUp(st->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STSetFromOptions_Cayley(ST st,PetscOptionItems PetscOptionsObject)
{
  PetscScalar    nu;
  PetscBool      flg;
  ST_CAYLEY      *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"ST Cayley Options");

    PetscCall(PetscOptionsScalar("-st_cayley_antishift","Value of the antishift","STCayleySetAntishift",ctx->nu,&nu,&flg));
    if (flg) PetscCall(STCayleySetAntishift(st,nu));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STCayleySetAntishift_Cayley(ST st,PetscScalar newshift)
{
  ST_CAYLEY *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  if (ctx->nu != newshift) {
    STCheckNotSeized(st,1);
    if (st->state && st->matmode!=ST_MATMODE_INPLACE) PetscCall(STMatMAXPY_Private(st,newshift,ctx->nu,0,NULL,PETSC_FALSE,PETSC_FALSE,&st->T[0]));
    ctx->nu = newshift;
  }
  ctx->nu_set = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STCayleySetAntishift - Sets the value of the anti-shift for the Cayley
   spectral transformation.

   Logically Collective

   Input Parameters:
+  st  - the spectral transformation context
-  nu  - the anti-shift

   Options Database Key:
.  -st_cayley_antishift - Sets the value of the anti-shift

   Level: intermediate

   Note:
   In the generalized Cayley transform, the operator can be expressed as
   OP = inv(A - sigma B)*(A + nu B). This function sets the value of nu.
   Use STSetShift() for setting sigma. The value nu=-sigma is not allowed.

.seealso: STSetShift(), STCayleyGetAntishift()
@*/
PetscErrorCode STCayleySetAntishift(ST st,PetscScalar nu)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveScalar(st,nu,2);
  PetscTryMethod(st,"STCayleySetAntishift_C",(ST,PetscScalar),(st,nu));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STCayleyGetAntishift_Cayley(ST st,PetscScalar *nu)
{
  ST_CAYLEY *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  *nu = ctx->nu;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STCayleyGetAntishift - Gets the value of the anti-shift used in the Cayley
   spectral transformation.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  nu  - the anti-shift

   Level: intermediate

.seealso: STGetShift(), STCayleySetAntishift()
@*/
PetscErrorCode STCayleyGetAntishift(ST st,PetscScalar *nu)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(nu,2);
  PetscUseMethod(st,"STCayleyGetAntishift_C",(ST,PetscScalar*),(st,nu));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STView_Cayley(ST st,PetscViewer viewer)
{
  char           str[50];
  ST_CAYLEY      *ctx = (ST_CAYLEY*)st->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->nu,PETSC_FALSE));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  antishift: %s\n",str));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STDestroy_Cayley(ST st)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(st->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STCayleySetAntishift_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STCayleyGetAntishift_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode STCreate_Cayley(ST st)
{
  ST_CAYLEY      *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  st->data = (void*)ctx;

  st->usesksp = PETSC_TRUE;

  st->ops->apply           = STApply_Generic;
  st->ops->applytrans      = STApplyTranspose_Generic;
  st->ops->applyhermtrans  = STApplyHermitianTranspose_Generic;
  st->ops->backtransform   = STBackTransform_Cayley;
  st->ops->setshift        = STSetShift_Cayley;
  st->ops->getbilinearform = STGetBilinearForm_Cayley;
  st->ops->setup           = STSetUp_Cayley;
  st->ops->computeoperator = STComputeOperator_Cayley;
  st->ops->setfromoptions  = STSetFromOptions_Cayley;
  st->ops->postsolve       = STPostSolve_Cayley;
  st->ops->destroy         = STDestroy_Cayley;
  st->ops->view            = STView_Cayley;
  st->ops->checknullspace  = STCheckNullSpace_Default;
  st->ops->setdefaultksp   = STSetDefaultKSP_Default;

  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STCayleySetAntishift_C",STCayleySetAntishift_Cayley));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STCayleyGetAntishift_C",STCayleyGetAntishift_Cayley));
  PetscFunctionReturn(PETSC_SUCCESS);
}
