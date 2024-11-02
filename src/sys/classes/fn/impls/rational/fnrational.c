/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Rational function  r(x) = p(x)/q(x), where p(x) and q(x) are polynomials
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/

typedef struct {
  PetscScalar *pcoeff;    /* numerator coefficients */
  PetscInt    np;         /* length of array pcoeff, p(x) has degree np-1 */
  PetscScalar *qcoeff;    /* denominator coefficients */
  PetscInt    nq;         /* length of array qcoeff, q(x) has degree nq-1 */
} FN_RATIONAL;

static PetscErrorCode FNEvaluateFunction_Rational(FN fn,PetscScalar x,PetscScalar *y)
{
  FN_RATIONAL *ctx = (FN_RATIONAL*)fn->data;
  PetscInt    i;
  PetscScalar p,q;

  PetscFunctionBegin;
  if (!ctx->np) p = 1.0;
  else {
    p = ctx->pcoeff[0];
    for (i=1;i<ctx->np;i++)
      p = ctx->pcoeff[i]+x*p;
  }
  if (!ctx->nq) *y = p;
  else {
    q = ctx->qcoeff[0];
    for (i=1;i<ctx->nq;i++)
      q = ctx->qcoeff[i]+x*q;
    PetscCheck(q!=0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Function not defined in the requested value");
    *y = p/q;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Horner evaluation of P=p(A)
   d = degree of polynomial;   coeff = coefficients of polynomial;    W = workspace
*/
static PetscErrorCode EvaluatePoly(Mat A,Mat P,Mat W,PetscInt d,PetscScalar *coeff)
{
  PetscInt j;

  PetscFunctionBegin;
  PetscCall(MatZeroEntries(P));
  if (!d) PetscCall(MatShift(P,1.0));
  else {
    PetscCall(MatShift(P,coeff[0]));
    for (j=1;j<d;j++) {
      PetscCall(MatMatMult(P,A,MAT_REUSE_MATRIX,PETSC_DEFAULT,&W));
      PetscCall(MatCopy(W,P,SAME_NONZERO_PATTERN));
      PetscCall(MatShift(P,coeff[j]));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Rational(FN fn,Mat A,Mat B)
{
  FN_RATIONAL *ctx = (FN_RATIONAL*)fn->data;
  Mat         P,Q,W,F;
  PetscBool   iscuda;

  PetscFunctionBegin;
  if (A==B) PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&P));
  else P = B;
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&W));

  PetscCall(EvaluatePoly(A,P,W,ctx->np,ctx->pcoeff));
  if (ctx->nq) {
    PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&Q));
    PetscCall(EvaluatePoly(A,Q,W,ctx->nq,ctx->qcoeff));
    PetscCall(PetscObjectTypeCompare((PetscObject)A,MATSEQDENSECUDA,&iscuda));
    PetscCall(MatGetFactor(Q,iscuda?MATSOLVERCUDA:MATSOLVERPETSC,MAT_FACTOR_LU,&F));
    PetscCall(MatLUFactorSymbolic(F,Q,NULL,NULL,NULL));
    PetscCall(MatLUFactorNumeric(F,Q,NULL));
    PetscCall(MatMatSolve(F,P,P));
    PetscCall(MatDestroy(&F));
    PetscCall(MatDestroy(&Q));
  }

  if (A==B) {
    PetscCall(MatCopy(P,B,SAME_NONZERO_PATTERN));
    PetscCall(MatDestroy(&P));
  }
  PetscCall(MatDestroy(&W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMatVec_Rational(FN fn,Mat A,Vec v)
{
  FN_RATIONAL *ctx = (FN_RATIONAL*)fn->data;
  Mat         P,Q,W,F;
  Vec         b;
  PetscBool   iscuda;

  PetscFunctionBegin;
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&P));
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&W));

  PetscCall(EvaluatePoly(A,P,W,ctx->np,ctx->pcoeff));
  if (ctx->nq) {
    PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&Q));
    PetscCall(EvaluatePoly(A,Q,W,ctx->nq,ctx->qcoeff));
    PetscCall(PetscObjectTypeCompare((PetscObject)A,MATSEQDENSECUDA,&iscuda));
    PetscCall(MatGetFactor(Q,iscuda?MATSOLVERCUDA:MATSOLVERPETSC,MAT_FACTOR_LU,&F));
    PetscCall(MatLUFactorSymbolic(F,Q,NULL,NULL,NULL));
    PetscCall(MatLUFactorNumeric(F,Q,NULL));
    PetscCall(MatCreateVecs(P,&b,NULL));
    PetscCall(MatGetColumnVector(P,b,0));
    PetscCall(MatSolve(F,b,v));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&F));
    PetscCall(MatDestroy(&Q));
  } else PetscCall(MatGetColumnVector(P,v,0));

  PetscCall(MatDestroy(&P));
  PetscCall(MatDestroy(&W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateDerivative_Rational(FN fn,PetscScalar x,PetscScalar *yp)
{
  FN_RATIONAL *ctx = (FN_RATIONAL*)fn->data;
  PetscInt    i;
  PetscScalar p,q,pp,qp;

  PetscFunctionBegin;
  if (!ctx->np) {
    p = 1.0;
    pp = 0.0;
  } else {
    p = ctx->pcoeff[0];
    pp = 0.0;
    for (i=1;i<ctx->np;i++) {
      pp = p+x*pp;
      p = ctx->pcoeff[i]+x*p;
    }
  }
  if (!ctx->nq) *yp = pp;
  else {
    q = ctx->qcoeff[0];
    qp = 0.0;
    for (i=1;i<ctx->nq;i++) {
      qp = q+x*qp;
      q = ctx->qcoeff[i]+x*q;
    }
    PetscCheck(q!=0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Derivative not defined in the requested value");
    *yp = (pp*q-p*qp)/(q*q);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNView_Rational(FN fn,PetscViewer viewer)
{
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscBool      isascii;
  PetscInt       i;
  char           str[50];

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (fn->alpha!=(PetscScalar)1.0 || fn->beta!=(PetscScalar)1.0) {
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->alpha,PETSC_FALSE));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  scale factors: alpha=%s,",str));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->beta,PETSC_FALSE));
      PetscCall(PetscViewerASCIIPrintf(viewer," beta=%s\n",str));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    }
    if (!ctx->nq) {
      if (!ctx->np) PetscCall(PetscViewerASCIIPrintf(viewer,"  constant: 1.0\n"));
      else if (ctx->np==1) {
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->pcoeff[0],PETSC_FALSE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"  constant: %s\n",str));
      } else {
        PetscCall(PetscViewerASCIIPrintf(viewer,"  polynomial: "));
        PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
        for (i=0;i<ctx->np-1;i++) {
          PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->pcoeff[i],PETSC_TRUE));
          PetscCall(PetscViewerASCIIPrintf(viewer,"%s*x^%1" PetscInt_FMT,str,ctx->np-i-1));
        }
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->pcoeff[ctx->np-1],PETSC_TRUE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"%s\n",str));
        PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
      }
    } else if (!ctx->np) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"  inverse polynomial: 1 / ("));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
      for (i=0;i<ctx->nq-1;i++) {
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->qcoeff[i],PETSC_TRUE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"%s*x^%1" PetscInt_FMT,str,ctx->nq-i-1));
      }
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->qcoeff[ctx->nq-1],PETSC_TRUE));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%s)\n",str));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    } else {
      PetscCall(PetscViewerASCIIPrintf(viewer,"  rational function: ("));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
      for (i=0;i<ctx->np-1;i++) {
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->pcoeff[i],PETSC_TRUE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"%s*x^%1" PetscInt_FMT,str,ctx->np-i-1));
      }
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->pcoeff[ctx->np-1],PETSC_TRUE));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%s) / (",str));
      for (i=0;i<ctx->nq-1;i++) {
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->qcoeff[i],PETSC_TRUE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"%s*x^%1" PetscInt_FMT,str,ctx->nq-i-1));
      }
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->qcoeff[ctx->nq-1],PETSC_TRUE));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%s)\n",str));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNRationalSetNumerator_Rational(FN fn,PetscInt np,PetscScalar *pcoeff)
{
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCheck(np>=0,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"Argument np cannot be negative");
  ctx->np = np;
  PetscCall(PetscFree(ctx->pcoeff));
  if (np) {
    PetscCall(PetscMalloc1(np,&ctx->pcoeff));
    for (i=0;i<np;i++) ctx->pcoeff[i] = pcoeff[i];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNRationalSetNumerator - Sets the parameters defining the numerator of the
   rational function.

   Logically Collective

   Input Parameters:
+  fn     - the math function context
.  np     - number of coefficients
-  pcoeff - coefficients (array of scalar values)

   Notes:
   Let the rational function r(x) = p(x)/q(x), where p(x) and q(x) are polynomials.
   This function provides the coefficients of the numerator p(x).
   Hence, p(x) is of degree np-1.
   If np is zero, then the numerator is assumed to be p(x)=1.

   In polynomials, high order coefficients are stored in the first positions
   of the array, e.g. to represent x^2-3 use {1,0,-3}.

   Level: intermediate

.seealso: FNRationalSetDenominator(), FNRationalGetNumerator()
@*/
PetscErrorCode FNRationalSetNumerator(FN fn,PetscInt np,PetscScalar pcoeff[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveInt(fn,np,2);
  if (np) PetscAssertPointer(pcoeff,3);
  PetscTryMethod(fn,"FNRationalSetNumerator_C",(FN,PetscInt,PetscScalar*),(fn,np,pcoeff));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNRationalGetNumerator_Rational(FN fn,PetscInt *np,PetscScalar *pcoeff[])
{
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (np) *np = ctx->np;
  if (pcoeff) {
    if (!ctx->np) *pcoeff = NULL;
    else {
      PetscCall(PetscMalloc1(ctx->np,pcoeff));
      for (i=0;i<ctx->np;i++) (*pcoeff)[i] = ctx->pcoeff[i];
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   FNRationalGetNumerator - Gets the parameters that define the numerator of the
   rational function.

   Not Collective

   Input Parameter:
.  fn     - the math function context

   Output Parameters:
+  np     - number of coefficients
-  pcoeff - coefficients (array of scalar values, length nq)

   Notes:
   The values passed by user with FNRationalSetNumerator() are returned (or null
   pointers otherwise).
   The pcoeff array should be freed by the user when no longer needed.

   Level: intermediate

.seealso: FNRationalSetNumerator()
@*/
PetscErrorCode FNRationalGetNumerator(FN fn,PetscInt *np,PetscScalar *pcoeff[]) PeNS
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscUseMethod(fn,"FNRationalGetNumerator_C",(FN,PetscInt*,PetscScalar**),(fn,np,pcoeff));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNRationalSetDenominator_Rational(FN fn,PetscInt nq,PetscScalar *qcoeff)
{
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCheck(nq>=0,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"Argument nq cannot be negative");
  ctx->nq = nq;
  PetscCall(PetscFree(ctx->qcoeff));
  if (nq) {
    PetscCall(PetscMalloc1(nq,&ctx->qcoeff));
    for (i=0;i<nq;i++) ctx->qcoeff[i] = qcoeff[i];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNRationalSetDenominator - Sets the parameters defining the denominator of the
   rational function.

   Logically Collective

   Input Parameters:
+  fn     - the math function context
.  nq     - number of coefficients
-  qcoeff - coefficients (array of scalar values)

   Notes:
   Let the rational function r(x) = p(x)/q(x), where p(x) and q(x) are polynomials.
   This function provides the coefficients of the denominator q(x).
   Hence, q(x) is of degree nq-1.
   If nq is zero, then the function is assumed to be polynomial, r(x) = p(x).

   In polynomials, high order coefficients are stored in the first positions
   of the array, e.g. to represent x^2-3 use {1,0,-3}.

   Level: intermediate

.seealso: FNRationalSetNumerator(), FNRationalGetDenominator()
@*/
PetscErrorCode FNRationalSetDenominator(FN fn,PetscInt nq,PetscScalar qcoeff[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveInt(fn,nq,2);
  if (nq) PetscAssertPointer(qcoeff,3);
  PetscTryMethod(fn,"FNRationalSetDenominator_C",(FN,PetscInt,PetscScalar*),(fn,nq,qcoeff));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNRationalGetDenominator_Rational(FN fn,PetscInt *nq,PetscScalar *qcoeff[])
{
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (nq) *nq = ctx->nq;
  if (qcoeff) {
    if (!ctx->nq) *qcoeff = NULL;
    else {
      PetscCall(PetscMalloc1(ctx->nq,qcoeff));
      for (i=0;i<ctx->nq;i++) (*qcoeff)[i] = ctx->qcoeff[i];
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   FNRationalGetDenominator - Gets the parameters that define the denominator of the
   rational function.

   Not Collective

   Input Parameter:
.  fn     - the math function context

   Output Parameters:
+  nq     - number of coefficients
-  qcoeff - coefficients (array of scalar values, length nq)

   Notes:
   The values passed by user with FNRationalSetDenominator() are returned (or a null
   pointer otherwise).
   The qcoeff array should be freed by the user when no longer needed.

   Level: intermediate

.seealso: FNRationalSetDenominator()
@*/
PetscErrorCode FNRationalGetDenominator(FN fn,PetscInt *nq,PetscScalar *qcoeff[]) PeNS
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscUseMethod(fn,"FNRationalGetDenominator_C",(FN,PetscInt*,PetscScalar**),(fn,nq,qcoeff));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNSetFromOptions_Rational(FN fn,PetscOptionItems PetscOptionsObject)
{
#define PARMAX 10
  PetscScalar    array[PARMAX];
  PetscInt       i,k;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"FN Rational Options");

    k = PARMAX;
    for (i=0;i<k;i++) array[i] = 0;
    PetscCall(PetscOptionsScalarArray("-fn_rational_numerator","Numerator coefficients (one or more scalar values separated with a comma without spaces)","FNRationalSetNumerator",array,&k,&flg));
    if (flg) PetscCall(FNRationalSetNumerator(fn,k,array));

    k = PARMAX;
    for (i=0;i<k;i++) array[i] = 0;
    PetscCall(PetscOptionsScalarArray("-fn_rational_denominator","Denominator coefficients (one or more scalar values separated with a comma without spaces)","FNRationalSetDenominator",array,&k,&flg));
    if (flg) PetscCall(FNRationalSetDenominator(fn,k,array));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNDuplicate_Rational(FN fn,MPI_Comm comm,FN *newfn)
{
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data,*ctx2 = (FN_RATIONAL*)(*newfn)->data;
  PetscInt       i;

  PetscFunctionBegin;
  ctx2->np = ctx->np;
  if (ctx->np) {
    PetscCall(PetscMalloc1(ctx->np,&ctx2->pcoeff));
    for (i=0;i<ctx->np;i++) ctx2->pcoeff[i] = ctx->pcoeff[i];
  }
  ctx2->nq = ctx->nq;
  if (ctx->nq) {
    PetscCall(PetscMalloc1(ctx->nq,&ctx2->qcoeff));
    for (i=0;i<ctx->nq;i++) ctx2->qcoeff[i] = ctx->qcoeff[i];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNDestroy_Rational(FN fn)
{
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;

  PetscFunctionBegin;
  PetscCall(PetscFree(ctx->pcoeff));
  PetscCall(PetscFree(ctx->qcoeff));
  PetscCall(PetscFree(fn->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNRationalSetNumerator_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNRationalGetNumerator_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNRationalSetDenominator_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNRationalGetDenominator_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode FNCreate_Rational(FN fn)
{
  FN_RATIONAL    *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  fn->data = (void*)ctx;

  fn->ops->evaluatefunction          = FNEvaluateFunction_Rational;
  fn->ops->evaluatederivative        = FNEvaluateDerivative_Rational;
  fn->ops->evaluatefunctionmat[0]    = FNEvaluateFunctionMat_Rational;
  fn->ops->evaluatefunctionmatvec[0] = FNEvaluateFunctionMatVec_Rational;
#if defined(PETSC_HAVE_CUDA)
  fn->ops->evaluatefunctionmatcuda[0]    = FNEvaluateFunctionMat_Rational;
  fn->ops->evaluatefunctionmatveccuda[0] = FNEvaluateFunctionMatVec_Rational;
#endif
  fn->ops->setfromoptions            = FNSetFromOptions_Rational;
  fn->ops->view                      = FNView_Rational;
  fn->ops->duplicate                 = FNDuplicate_Rational;
  fn->ops->destroy                   = FNDestroy_Rational;
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNRationalSetNumerator_C",FNRationalSetNumerator_Rational));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNRationalGetNumerator_C",FNRationalGetNumerator_Rational));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNRationalSetDenominator_C",FNRationalSetDenominator_Rational));
  PetscCall(PetscObjectComposeFunction((PetscObject)fn,"FNRationalGetDenominator_C",FNRationalGetDenominator_Rational));
  PetscFunctionReturn(PETSC_SUCCESS);
}
