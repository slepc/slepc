/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   A function that is obtained by combining two other functions (either by
   addition, multiplication, division or composition)

      addition:          f(x) = f1(x)+f2(x)
      multiplication:    f(x) = f1(x)*f2(x)
      division:          f(x) = f1(x)/f2(x)      f(A) = f2(A)\f1(A)
      composition:       f(x) = f2(f1(x))
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/

typedef struct {
  FN            f1,f2;    /* functions */
  FNCombineType comb;     /* how the functions are combined */
} FN_COMBINE;

static PetscErrorCode FNEvaluateFunction_Combine(FN fn,PetscScalar x,PetscScalar *y)
{
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscScalar    a,b;

  PetscFunctionBegin;
  PetscCall(FNEvaluateFunction(ctx->f1,x,&a));
  switch (ctx->comb) {
    case FN_COMBINE_ADD:
      PetscCall(FNEvaluateFunction(ctx->f2,x,&b));
      *y = a+b;
      break;
    case FN_COMBINE_MULTIPLY:
      PetscCall(FNEvaluateFunction(ctx->f2,x,&b));
      *y = a*b;
      break;
    case FN_COMBINE_DIVIDE:
      PetscCall(FNEvaluateFunction(ctx->f2,x,&b));
      PetscCheck(b!=0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Function not defined in the requested value");
      *y = a/b;
      break;
    case FN_COMBINE_COMPOSE:
      PetscCall(FNEvaluateFunction(ctx->f2,a,y));
      break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateDerivative_Combine(FN fn,PetscScalar x,PetscScalar *yp)
{
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscScalar    a,b,ap,bp;

  PetscFunctionBegin;
  switch (ctx->comb) {
    case FN_COMBINE_ADD:
      PetscCall(FNEvaluateDerivative(ctx->f1,x,&ap));
      PetscCall(FNEvaluateDerivative(ctx->f2,x,&bp));
      *yp = ap+bp;
      break;
    case FN_COMBINE_MULTIPLY:
      PetscCall(FNEvaluateDerivative(ctx->f1,x,&ap));
      PetscCall(FNEvaluateDerivative(ctx->f2,x,&bp));
      PetscCall(FNEvaluateFunction(ctx->f1,x,&a));
      PetscCall(FNEvaluateFunction(ctx->f2,x,&b));
      *yp = ap*b+a*bp;
      break;
    case FN_COMBINE_DIVIDE:
      PetscCall(FNEvaluateDerivative(ctx->f1,x,&ap));
      PetscCall(FNEvaluateDerivative(ctx->f2,x,&bp));
      PetscCall(FNEvaluateFunction(ctx->f1,x,&a));
      PetscCall(FNEvaluateFunction(ctx->f2,x,&b));
      PetscCheck(b!=0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Derivative not defined in the requested value");
      *yp = (ap*b-a*bp)/(b*b);
      break;
    case FN_COMBINE_COMPOSE:
      PetscCall(FNEvaluateFunction(ctx->f1,x,&a));
      PetscCall(FNEvaluateDerivative(ctx->f1,x,&ap));
      PetscCall(FNEvaluateDerivative(ctx->f2,a,yp));
      *yp *= ap;
      break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Combine(FN fn,Mat A,Mat B)
{
  FN_COMBINE   *ctx = (FN_COMBINE*)fn->data;
  Mat          W,Z,F;
  PetscBool    iscuda;

  PetscFunctionBegin;
  PetscCall(FN_AllocateWorkMat(fn,A,&W));
  switch (ctx->comb) {
    case FN_COMBINE_ADD:
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f1,A,W,PETSC_FALSE));
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f2,A,B,PETSC_FALSE));
      PetscCall(MatAXPY(B,1.0,W,SAME_NONZERO_PATTERN));
      break;
    case FN_COMBINE_MULTIPLY:
      PetscCall(FN_AllocateWorkMat(fn,A,&Z));
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f1,A,W,PETSC_FALSE));
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f2,A,Z,PETSC_FALSE));
      PetscCall(MatMatMult(W,Z,MAT_REUSE_MATRIX,PETSC_DEFAULT,&B));
      PetscCall(FN_FreeWorkMat(fn,&Z));
      break;
    case FN_COMBINE_DIVIDE:
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f2,A,W,PETSC_FALSE));
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f1,A,B,PETSC_FALSE));
      PetscCall(PetscObjectTypeCompare((PetscObject)A,MATSEQDENSECUDA,&iscuda));
      PetscCall(MatGetFactor(W,iscuda?MATSOLVERCUDA:MATSOLVERPETSC,MAT_FACTOR_LU,&F));
      PetscCall(MatLUFactorSymbolic(F,W,NULL,NULL,NULL));
      PetscCall(MatLUFactorNumeric(F,W,NULL));
      PetscCall(MatMatSolve(F,B,B));
      PetscCall(MatDestroy(&F));
      break;
    case FN_COMBINE_COMPOSE:
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f1,A,W,PETSC_FALSE));
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f2,W,B,PETSC_FALSE));
      break;
  }
  PetscCall(FN_FreeWorkMat(fn,&W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMatVec_Combine(FN fn,Mat A,Vec v)
{
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscBool      iscuda;
  Mat            Z,F;
  Vec            w;

  PetscFunctionBegin;
  switch (ctx->comb) {
    case FN_COMBINE_ADD:
      PetscCall(VecDuplicate(v,&w));
      PetscCall(FNEvaluateFunctionMatVec_Private(ctx->f1,A,w,PETSC_FALSE));
      PetscCall(FNEvaluateFunctionMatVec_Private(ctx->f2,A,v,PETSC_FALSE));
      PetscCall(VecAXPY(v,1.0,w));
      PetscCall(VecDestroy(&w));
      break;
    case FN_COMBINE_MULTIPLY:
      PetscCall(VecDuplicate(v,&w));
      PetscCall(FN_AllocateWorkMat(fn,A,&Z));
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f1,A,Z,PETSC_FALSE));
      PetscCall(FNEvaluateFunctionMatVec_Private(ctx->f2,A,w,PETSC_FALSE));
      PetscCall(MatMult(Z,w,v));
      PetscCall(FN_FreeWorkMat(fn,&Z));
      PetscCall(VecDestroy(&w));
      break;
    case FN_COMBINE_DIVIDE:
      PetscCall(VecDuplicate(v,&w));
      PetscCall(FN_AllocateWorkMat(fn,A,&Z));
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f2,A,Z,PETSC_FALSE));
      PetscCall(FNEvaluateFunctionMatVec_Private(ctx->f1,A,w,PETSC_FALSE));
      PetscCall(PetscObjectTypeCompare((PetscObject)A,MATSEQDENSECUDA,&iscuda));
      PetscCall(MatGetFactor(Z,iscuda?MATSOLVERCUDA:MATSOLVERPETSC,MAT_FACTOR_LU,&F));
      PetscCall(MatLUFactorSymbolic(F,Z,NULL,NULL,NULL));
      PetscCall(MatLUFactorNumeric(F,Z,NULL));
      PetscCall(MatSolve(F,w,v));
      PetscCall(MatDestroy(&F));
      PetscCall(FN_FreeWorkMat(fn,&Z));
      PetscCall(VecDestroy(&w));
      break;
    case FN_COMBINE_COMPOSE:
      PetscCall(FN_AllocateWorkMat(fn,A,&Z));
      PetscCall(FNEvaluateFunctionMat_Private(ctx->f1,A,Z,PETSC_FALSE));
      PetscCall(FNEvaluateFunctionMatVec_Private(ctx->f2,Z,v,PETSC_FALSE));
      PetscCall(FN_FreeWorkMat(fn,&Z));
      break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNView_Combine(FN fn,PetscViewer viewer)
{
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscBool      isascii;
  char           str[50];

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (fn->alpha!=(PetscScalar)1.0 || fn->beta!=(PetscScalar)1.0) {
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->alpha,PETSC_TRUE));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  scaling factors: inner %s",str));
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->beta,PETSC_TRUE));
      PetscCall(PetscViewerASCIIPrintf(viewer,"outer %s\n",str));
    }
    switch (ctx->comb) {
      case FN_COMBINE_ADD:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  two added functions f1+f2\n"));
        break;
      case FN_COMBINE_MULTIPLY:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  two multiplied functions f1*f2\n"));
        break;
      case FN_COMBINE_DIVIDE:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  a quotient of two functions f1/f2\n"));
        break;
      case FN_COMBINE_COMPOSE:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  two composed functions f2(f1(.))\n"));
        break;
    }
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(FNView(ctx->f1,viewer));
    PetscCall(FNView(ctx->f2,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNCombineSetChildren_Combine(FN fn,FNCombineType comb,FN f1,FN f2)
{
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;

  PetscFunctionBegin;
  ctx->comb = comb;
  PetscCall(PetscObjectReference((PetscObject)f1));
  PetscCall(FNDestroy(&ctx->f1));
  ctx->f1 = f1;
  PetscCall(PetscObjectReference((PetscObject)f2));
  PetscCall(FNDestroy(&ctx->f2));
  ctx->f2 = f2;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNCombineSetChildren - Sets the two child functions that constitute this
   combined function, and the way they must be combined.

   Logically Collective

   Input Parameters:
+  fn   - the math function context
.  comb - how to combine the functions (addition, multiplication, division or composition)
.  f1   - first function
-  f2   - second function

   Level: intermediate

.seealso: [](sec:fn), `FNCOMBINE`, `FNCombineGetChildren()`
@*/
PetscErrorCode FNCombineSetChildren(FN fn,FNCombineType comb,FN f1,FN f2)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveEnum(fn,comb,2);
  PetscValidHeaderSpecific(f1,FN_CLASSID,3);
  PetscValidHeaderSpecific(f2,FN_CLASSID,4);
  PetscTryMethod(fn,"FNCombineSetChildren_C",(FN,FNCombineType,FN,FN),(fn,comb,f1,f2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNCombineGetChildren_Combine(FN fn,FNCombineType *comb,FN *f1,FN *f2)
{
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;

  PetscFunctionBegin;
  if (comb) *comb = ctx->comb;
  if (f1) {
    if (!ctx->f1) PetscCall(FNCreate(PetscObjectComm((PetscObject)fn),&ctx->f1));
    *f1 = ctx->f1;
  }
  if (f2) {
    if (!ctx->f2) PetscCall(FNCreate(PetscObjectComm((PetscObject)fn),&ctx->f2));
    *f2 = ctx->f2;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNCombineGetChildren - Gets the two child functions that constitute this
   combined function, and the way they are combined.

   Not Collective

   Input Parameter:
.  fn   - the math function context

   Output Parameters:
+  comb - how to combine the functions (addition, multiplication, division or composition)
.  f1   - first function
-  f2   - second function

   Level: intermediate

.seealso: [](sec:fn), `FNCOMBINE`, `FNCombineSetChildren()`
@*/
PetscErrorCode FNCombineGetChildren(FN fn,FNCombineType *comb,FN *f1,FN *f2)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscUseMethod(fn,"FNCombineGetChildren_C",(FN,FNCombineType*,FN*,FN*),(fn,comb,f1,f2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNDuplicate_Combine(FN fn,MPI_Comm comm,FN *newfn)
{
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data,*ctx2 = (FN_COMBINE*)(*newfn)->data;

  PetscFunctionBegin;
  ctx2->comb = ctx->comb;
  PetscCall(FNDuplicate(ctx->f1,comm,&ctx2->f1));
  PetscCall(FNDuplicate(ctx->f2,comm,&ctx2->f2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNDestroy_Combine(FN fn)
{
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;

  PetscFunctionBegin;
  PetscCall(FNDestroy(&ctx->f1));
  PetscCall(FNDestroy(&ctx->f2));
  PetscCall(PetscFree(fn->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNCombineSetChildren_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNCombineGetChildren_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   FNCOMBINE - FNCOMBINE = "combine" - A mathematical function that is defined
   by combining two previously defined functions, e.g., by addition. Use
   `FNCombineSetChildren()` to specify which functions must be combined and how.

   Level: beginner

.seealso: [](sec:fn), `FN`, `FNType`, `FNSetType()`, `FNCombineSetChildren()`
M*/

SLEPC_EXTERN PetscErrorCode FNCreate_Combine(FN fn)
{
  FN_COMBINE     *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  fn->data = (void*)ctx;

  fn->ops->evaluatefunction          = FNEvaluateFunction_Combine;
  fn->ops->evaluatederivative        = FNEvaluateDerivative_Combine;
  fn->ops->evaluatefunctionmat[0]    = FNEvaluateFunctionMat_Combine;
  fn->ops->evaluatefunctionmatvec[0] = FNEvaluateFunctionMatVec_Combine;
#if defined(PETSC_HAVE_CUDA)
  fn->ops->evaluatefunctionmatcuda[0]    = FNEvaluateFunctionMat_Combine;
  fn->ops->evaluatefunctionmatveccuda[0] = FNEvaluateFunctionMatVec_Combine;
#endif
  fn->ops->view                      = FNView_Combine;
  fn->ops->duplicate                 = FNDuplicate_Combine;
  fn->ops->destroy                   = FNDestroy_Combine;
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNCombineSetChildren_C",FNCombineSetChildren_Combine));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNCombineGetChildren_C",FNCombineGetChildren_Combine));
  PetscFunctionReturn(PETSC_SUCCESS);
}
