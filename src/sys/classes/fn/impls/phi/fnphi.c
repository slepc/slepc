/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Phi functions
      phi_0(x) = exp(x)
      phi_1(x) = (exp(x)-1)/x
      phi_k(x) = (phi_{k-1}(x)-1/(k-1)!)/x
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/

typedef struct {
  PetscInt k;    /* index of the phi-function, defaults to k=1 */
  Mat      H;    /* auxiliary matrix of order m+k */
  Mat      F;    /* auxiliary matrix to store exp(H) */
} FN_PHI;

#define MAX_INDEX 10

static const PetscReal rfactorial[MAX_INDEX+2] = { 1, 1, 0.5, 1.0/6, 1.0/24, 1.0/120, 1.0/720, 1.0/5040, 1.0/40320, 1.0/362880, 1.0/3628800, 1.0/39916800 };

static PetscErrorCode FNEvaluateFunction_Phi(FN fn,PetscScalar x,PetscScalar *y)
{
  FN_PHI      *ctx = (FN_PHI*)fn->data;
  PetscInt    i;
  PetscScalar phi[MAX_INDEX+1];

  PetscFunctionBegin;
  if (x==0.0) *y = rfactorial[ctx->k];
  else {
    phi[0] = PetscExpScalar(x);
    for (i=1;i<=ctx->k;i++) phi[i] = (phi[i-1]-rfactorial[i-1])/x;
    *y = phi[ctx->k];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateDerivative_Phi(FN fn,PetscScalar x,PetscScalar *y)
{
  FN_PHI      *ctx = (FN_PHI*)fn->data;
  PetscInt    i;
  PetscScalar phi[MAX_INDEX+2];

  PetscFunctionBegin;
  if (x==0.0) *y = rfactorial[ctx->k+1];
  else {
    phi[0] = PetscExpScalar(x);
    for (i=1;i<=ctx->k+1;i++) phi[i] = (phi[i-1]-rfactorial[i-1])/x;
    *y = phi[ctx->k] - phi[ctx->k+1]*(PetscReal)ctx->k;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMatVec_Phi(FN fn,Mat A,Vec v)
{
  FN_PHI            *ctx = (FN_PHI*)fn->data;
  PetscInt          i,j,m,n,nh;
  PetscScalar       *Ha,*va,sfactor=1.0;
  const PetscScalar *Aa,*Fa;

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,&m,NULL));
  n = m+ctx->k;
  if (ctx->H) {
    PetscCall(MatGetSize(ctx->H,&nh,NULL));
    if (n!=nh) {
      PetscCall(MatDestroy(&ctx->H));
      PetscCall(MatDestroy(&ctx->F));
    }
  }
  if (!ctx->H) {
    PetscCall(MatCreateDense(PETSC_COMM_SELF,n,n,n,n,NULL,&ctx->H));
    PetscCall(MatCreateDense(PETSC_COMM_SELF,n,n,n,n,NULL,&ctx->F));
  }
  PetscCall(MatDenseGetArray(ctx->H,&Ha));
  PetscCall(MatDenseGetArrayRead(A,&Aa));
  for (j=0;j<m;j++) PetscCall(PetscArraycpy(Ha+j*n,Aa+j*m,m));
  PetscCall(MatDenseRestoreArrayRead(A,&Aa));
  if (ctx->k) {
    for (j=0;j<m;j++) for (i=m;i<n;i++) Ha[i+j*n] = 0.0;
    for (j=m;j<n;j++) for (i=0;i<n;i++) Ha[i+j*n] = 0.0;
    Ha[0+m*n] = fn->alpha;
    for (j=m+1;j<n;j++) Ha[j-1+j*n] = fn->alpha;
  }
  PetscCall(MatDenseRestoreArray(ctx->H,&Ha));

  PetscCall(FNEvaluateFunctionMat_Exp_Higham(fn,ctx->H,ctx->F));

  PetscCall(MatDenseGetArrayRead(ctx->F,&Fa));
  PetscCall(VecGetArray(v,&va));
  if (ctx->k) {
    sfactor = PetscPowScalarInt(fn->alpha,-ctx->k);
    for (i=0;i<m;i++) va[i] = sfactor*Fa[i+(n-1)*n];
  } else {
    for (i=0;i<m;i++) va[i] = sfactor*Fa[i+0*n];
  }
  PetscCall(VecRestoreArray(v,&va));
  PetscCall(MatDenseRestoreArrayRead(ctx->F,&Fa));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNPhiSetIndex_Phi(FN fn,PetscInt k)
{
  FN_PHI         *ctx = (FN_PHI*)fn->data;

  PetscFunctionBegin;
  PetscCheck(k>=0,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"Index cannot be negative");
  PetscCheck(k<=MAX_INDEX,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"Phi functions only implemented for k<=%d",MAX_INDEX);
  if (k!=ctx->k) {
    ctx->k = k;
    PetscCall(MatDestroy(&ctx->H));
    PetscCall(MatDestroy(&ctx->F));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNPhiSetIndex - Sets the index of the phi-function.

   Logically Collective

   Input Parameters:
+  fn - the math function context
-  k  - the index

   Notes:
   The phi-functions are defined as follows. The default is k=1.
.vb
      phi_0(x) = exp(x)
      phi_1(x) = (exp(x)-1)/x
      phi_k(x) = (phi_{k-1}(x)-1/(k-1)!)/x
.ve

   Level: intermediate

.seealso: FNPhiGetIndex()
@*/
PetscErrorCode FNPhiSetIndex(FN fn,PetscInt k)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveInt(fn,k,2);
  PetscTryMethod(fn,"FNPhiSetIndex_C",(FN,PetscInt),(fn,k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNPhiGetIndex_Phi(FN fn,PetscInt *k)
{
  FN_PHI *ctx = (FN_PHI*)fn->data;

  PetscFunctionBegin;
  *k = ctx->k;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNPhiGetIndex - Gets the index of the phi-function.

   Not Collective

   Input Parameter:
.  fn - the math function context

   Output Parameter:
.  k  - the index

   Level: intermediate

.seealso: FNPhiSetIndex()
@*/
PetscErrorCode FNPhiGetIndex(FN fn,PetscInt *k)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscAssertPointer(k,2);
  PetscUseMethod(fn,"FNPhiGetIndex_C",(FN,PetscInt*),(fn,k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNView_Phi(FN fn,PetscViewer viewer)
{
  FN_PHI         *ctx = (FN_PHI*)fn->data;
  PetscBool      isascii;
  char           str[50],strx[50];

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  phi_%" PetscInt_FMT ": ",ctx->k));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    if (fn->beta!=(PetscScalar)1.0) {
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->beta,PETSC_TRUE));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%s*",str));
    }
    if (fn->alpha==(PetscScalar)1.0) PetscCall(PetscSNPrintf(strx,sizeof(strx),"x"));
    else {
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->alpha,PETSC_TRUE));
      PetscCall(PetscSNPrintf(strx,sizeof(strx),"(%s*x)",str));
    }
    if (!ctx->k) PetscCall(PetscViewerASCIIPrintf(viewer,"exp(%s)\n",strx));
    else if (ctx->k==1) PetscCall(PetscViewerASCIIPrintf(viewer,"(exp(%s)-1)/%s\n",strx,strx));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"(phi_%" PetscInt_FMT "(%s)-1/%" PetscInt_FMT "!)/%s\n",ctx->k-1,strx,ctx->k-1,strx));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNSetFromOptions_Phi(FN fn,PetscOptionItems PetscOptionsObject)
{
  FN_PHI         *ctx = (FN_PHI*)fn->data;
  PetscInt       k;
  PetscBool      flag;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"FN Phi Options");

    PetscCall(PetscOptionsInt("-fn_phi_index","Index of the phi-function","FNPhiSetIndex",ctx->k,&k,&flag));
    if (flag) PetscCall(FNPhiSetIndex(fn,k));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNDuplicate_Phi(FN fn,MPI_Comm comm,FN *newfn)
{
  FN_PHI *ctx = (FN_PHI*)fn->data,*ctx2 = (FN_PHI*)(*newfn)->data;

  PetscFunctionBegin;
  ctx2->k = ctx->k;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNDestroy_Phi(FN fn)
{
  FN_PHI         *ctx = (FN_PHI*)fn->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->H));
  PetscCall(MatDestroy(&ctx->F));
  PetscCall(PetscFree(fn->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNPhiSetIndex_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNPhiGetIndex_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode FNCreate_Phi(FN fn)
{
  FN_PHI         *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  fn->data = (void*)ctx;
  ctx->k   = 1;

  fn->ops->evaluatefunction          = FNEvaluateFunction_Phi;
  fn->ops->evaluatederivative        = FNEvaluateDerivative_Phi;
  fn->ops->evaluatefunctionmatvec[0] = FNEvaluateFunctionMatVec_Phi;
  fn->ops->setfromoptions            = FNSetFromOptions_Phi;
  fn->ops->view                      = FNView_Phi;
  fn->ops->duplicate                 = FNDuplicate_Phi;
  fn->ops->destroy                   = FNDestroy_Phi;
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNPhiSetIndex_C",FNPhiSetIndex_Phi));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNPhiGetIndex_C",FNPhiGetIndex_Phi));
  PetscFunctionReturn(PETSC_SUCCESS);
}
