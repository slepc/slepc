/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Inverse square root function  x^(-1/2)
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/
#include <slepcblaslapack.h>

static PetscErrorCode FNEvaluateFunction_Invsqrt(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscFunctionBegin;
  PetscCheck(x!=0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Function not defined in the requested value");
#if !defined(PETSC_USE_COMPLEX)
  PetscCheck(x>0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Function not defined in the requested value");
#endif
  *y = 1.0/PetscSqrtScalar(x);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateDerivative_Invsqrt(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscFunctionBegin;
  PetscCheck(x!=0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Derivative not defined in the requested value");
#if !defined(PETSC_USE_COMPLEX)
  PetscCheck(x>0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Derivative not defined in the requested value");
#endif
  *y = -1.0/(2.0*PetscPowScalarReal(x,1.5));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Invsqrt_Schur(FN fn,Mat A,Mat B)
{
  PetscBLASInt   n=0,ld,*ipiv,info;
  PetscScalar    *Ba,*Wa;
  PetscInt       m;
  Mat            W;

  PetscFunctionBegin;
  PetscCall(FN_AllocateWorkMat(fn,A,&W));
  if (A!=B) PetscCall(MatCopy(A,B,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseGetArray(B,&Ba));
  PetscCall(MatDenseGetArray(W,&Wa));
  /* compute B = sqrtm(A) */
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  ld = n;
  PetscCall(FNSqrtmSchur(fn,n,Ba,n,PETSC_FALSE));
  /* compute B = A\B */
  PetscCall(PetscMalloc1(ld,&ipiv));
  PetscCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&n,Wa,&ld,ipiv,Ba,&ld,&info));
  SlepcCheckLapackInfo("gesv",info);
  PetscCall(PetscLogFlops(2.0*n*n*n/3.0+2.0*n*n*n));
  PetscCall(PetscFree(ipiv));
  PetscCall(MatDenseRestoreArray(W,&Wa));
  PetscCall(MatDenseRestoreArray(B,&Ba));
  PetscCall(FN_FreeWorkMat(fn,&W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMatVec_Invsqrt_Schur(FN fn,Mat A,Vec v)
{
  PetscBLASInt   n=0,ld,*ipiv,info,one=1;
  PetscScalar    *Ba,*Wa;
  PetscInt       m;
  Mat            B,W;

  PetscFunctionBegin;
  PetscCall(FN_AllocateWorkMat(fn,A,&B));
  PetscCall(FN_AllocateWorkMat(fn,A,&W));
  PetscCall(MatDenseGetArray(B,&Ba));
  PetscCall(MatDenseGetArray(W,&Wa));
  /* compute B_1 = sqrtm(A)*e_1 */
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  ld = n;
  PetscCall(FNSqrtmSchur(fn,n,Ba,n,PETSC_TRUE));
  /* compute B_1 = A\B_1 */
  PetscCall(PetscMalloc1(ld,&ipiv));
  PetscCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&one,Wa,&ld,ipiv,Ba,&ld,&info));
  SlepcCheckLapackInfo("gesv",info);
  PetscCall(PetscFree(ipiv));
  PetscCall(MatDenseRestoreArray(W,&Wa));
  PetscCall(MatDenseRestoreArray(B,&Ba));
  PetscCall(MatGetColumnVector(B,v,0));
  PetscCall(FN_FreeWorkMat(fn,&W));
  PetscCall(FN_FreeWorkMat(fn,&B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Invsqrt_DBP(FN fn,Mat A,Mat B)
{
  PetscBLASInt   n=0;
  PetscScalar    *T;
  PetscInt       m;

  PetscFunctionBegin;
  if (A!=B) PetscCall(MatCopy(A,B,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseGetArray(B,&T));
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(FNSqrtmDenmanBeavers(fn,n,T,n,PETSC_TRUE));
  PetscCall(MatDenseRestoreArray(B,&T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Invsqrt_NS(FN fn,Mat A,Mat B)
{
  PetscBLASInt   n=0;
  PetscScalar    *T;
  PetscInt       m;

  PetscFunctionBegin;
  if (A!=B) PetscCall(MatCopy(A,B,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseGetArray(B,&T));
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(FNSqrtmNewtonSchulz(fn,n,T,n,PETSC_TRUE));
  PetscCall(MatDenseRestoreArray(B,&T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Invsqrt_Sadeghi(FN fn,Mat A,Mat B)
{
  PetscBLASInt   n=0,ld,*ipiv,info;
  PetscScalar    *Ba,*Wa;
  PetscInt       m;
  Mat            W;

  PetscFunctionBegin;
  PetscCall(FN_AllocateWorkMat(fn,A,&W));
  if (A!=B) PetscCall(MatCopy(A,B,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseGetArray(B,&Ba));
  PetscCall(MatDenseGetArray(W,&Wa));
  /* compute B = sqrtm(A) */
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  ld = n;
  PetscCall(FNSqrtmSadeghi(fn,n,Ba,n));
  /* compute B = A\B */
  PetscCall(PetscMalloc1(ld,&ipiv));
  PetscCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&n,Wa,&ld,ipiv,Ba,&ld,&info));
  SlepcCheckLapackInfo("gesv",info);
  PetscCall(PetscLogFlops(2.0*n*n*n/3.0+2.0*n*n*n));
  PetscCall(PetscFree(ipiv));
  PetscCall(MatDenseRestoreArray(W,&Wa));
  PetscCall(MatDenseRestoreArray(B,&Ba));
  PetscCall(FN_FreeWorkMat(fn,&W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_CUDA)
PetscErrorCode FNEvaluateFunctionMat_Invsqrt_NS_CUDA(FN fn,Mat A,Mat B)
{
  PetscBLASInt   n=0;
  PetscScalar    *Ba;
  PetscInt       m;

  PetscFunctionBegin;
  if (A!=B) PetscCall(MatCopy(A,B,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseCUDAGetArray(B,&Ba));
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(FNSqrtmNewtonSchulz_CUDA(fn,n,Ba,n,PETSC_TRUE));
  PetscCall(MatDenseCUDARestoreArray(B,&Ba));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_MAGMA)
#include <slepcmagma.h>

PetscErrorCode FNEvaluateFunctionMat_Invsqrt_DBP_CUDAm(FN fn,Mat A,Mat B)
{
  PetscBLASInt   n=0;
  PetscScalar    *T;
  PetscInt       m;

  PetscFunctionBegin;
  if (A!=B) PetscCall(MatCopy(A,B,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseCUDAGetArray(B,&T));
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(FNSqrtmDenmanBeavers_CUDAm(fn,n,T,n,PETSC_TRUE));
  PetscCall(MatDenseCUDARestoreArray(B,&T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FNEvaluateFunctionMat_Invsqrt_Sadeghi_CUDAm(FN fn,Mat A,Mat B)
{
  PetscBLASInt   n=0,ld,*ipiv;
  PetscScalar    *Ba,*Wa;
  PetscInt       m;
  Mat            W;

  PetscFunctionBegin;
  PetscCall(FN_AllocateWorkMat(fn,A,&W));
  if (A!=B) PetscCall(MatCopy(A,B,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseCUDAGetArray(B,&Ba));
  PetscCall(MatDenseCUDAGetArray(W,&Wa));
  /* compute B = sqrtm(A) */
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  ld = n;
  PetscCall(FNSqrtmSadeghi_CUDAm(fn,n,Ba,n));
  /* compute B = A\B */
  PetscCall(SlepcMagmaInit());
  PetscCall(PetscMalloc1(ld,&ipiv));
  PetscCallMAGMA(magma_xgesv_gpu,n,n,Wa,ld,ipiv,Ba,ld);
  PetscCall(PetscLogFlops(2.0*n*n*n/3.0+2.0*n*n*n));
  PetscCall(PetscFree(ipiv));
  PetscCall(MatDenseCUDARestoreArray(W,&Wa));
  PetscCall(MatDenseCUDARestoreArray(B,&Ba));
  PetscCall(FN_FreeWorkMat(fn,&W));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* PETSC_HAVE_MAGMA */
#endif /* PETSC_HAVE_CUDA */

static PetscErrorCode FNView_Invsqrt(FN fn,PetscViewer viewer)
{
  PetscBool      isascii;
  char           str[50];
  const char     *methodname[] = {
                  "Schur method for inv(A)*sqrtm(A)",
                  "Denman-Beavers (product form)",
                  "Newton-Schulz iteration",
                  "Sadeghi iteration"
  };
  const int      nmeth=PETSC_STATIC_ARRAY_LENGTH(methodname);

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (fn->beta==(PetscScalar)1.0) {
      if (fn->alpha==(PetscScalar)1.0) PetscCall(PetscViewerASCIIPrintf(viewer,"  inverse square root: x^(-1/2)\n"));
      else {
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->alpha,PETSC_TRUE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"  inverse square root: (%s*x)^(-1/2)\n",str));
      }
    } else {
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->beta,PETSC_TRUE));
      if (fn->alpha==(PetscScalar)1.0) PetscCall(PetscViewerASCIIPrintf(viewer,"  inverse square root: %s*x^(-1/2)\n",str));
      else {
        PetscCall(PetscViewerASCIIPrintf(viewer,"  inverse square root: %s",str));
        PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->alpha,PETSC_TRUE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"*(%s*x)^(-1/2)\n",str));
        PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
      }
    }
    if (fn->method<nmeth) PetscCall(PetscViewerASCIIPrintf(viewer,"  computing matrix functions with: %s\n",methodname[fn->method]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   FNINVSQRT - FNINVSQRT = "invsqrt" - The inverse square root function $f(x)=x^{-\frac{1}{2}}$.

   Level: beginner

.seealso: [](sec:fn), `FN`, `FNType`, `FNSetType()`
M*/

SLEPC_EXTERN PetscErrorCode FNCreate_Invsqrt(FN fn)
{
  PetscFunctionBegin;
  fn->ops->evaluatefunction          = FNEvaluateFunction_Invsqrt;
  fn->ops->evaluatederivative        = FNEvaluateDerivative_Invsqrt;
  fn->ops->evaluatefunctionmat[0]    = FNEvaluateFunctionMat_Invsqrt_Schur;
  fn->ops->evaluatefunctionmat[1]    = FNEvaluateFunctionMat_Invsqrt_DBP;
  fn->ops->evaluatefunctionmat[2]    = FNEvaluateFunctionMat_Invsqrt_NS;
  fn->ops->evaluatefunctionmat[3]    = FNEvaluateFunctionMat_Invsqrt_Sadeghi;
#if defined(PETSC_HAVE_CUDA)
  fn->ops->evaluatefunctionmatcuda[2] = FNEvaluateFunctionMat_Invsqrt_NS_CUDA;
#if defined(PETSC_HAVE_MAGMA)
  fn->ops->evaluatefunctionmatcuda[1] = FNEvaluateFunctionMat_Invsqrt_DBP_CUDAm;
  fn->ops->evaluatefunctionmatcuda[3] = FNEvaluateFunctionMat_Invsqrt_Sadeghi_CUDAm;
#endif /* PETSC_HAVE_MAGMA */
#endif /* PETSC_HAVE_CUDA */
  fn->ops->evaluatefunctionmatvec[0] = FNEvaluateFunctionMatVec_Invsqrt_Schur;
  fn->ops->view                      = FNView_Invsqrt;
  PetscFunctionReturn(PETSC_SUCCESS);
}
