/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV private kernels that use the BLAS
*/

#include <slepc/private/bvimpl.h>
#include <slepcblaslapack.h>

#define BLOCKSIZE 64

/*
    C := alpha*A*B + beta*C

    A is mxk (ld=lda), B is kxn (ld=ldb), C is mxn (ld=ldc)
*/
PetscErrorCode BVMult_BLAS_Private(BV bv,PetscInt m_,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *A,PetscInt lda_,const PetscScalar *B,PetscInt ldb_,PetscScalar beta,PetscScalar *C,PetscInt ldc_)
{
  PetscBLASInt   m,n,k,lda,ldb,ldc;
#if defined(PETSC_HAVE_FBLASLAPACK) || defined(PETSC_HAVE_F2CBLASLAPACK)
  PetscBLASInt   l,bs=BLOCKSIZE;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(m_,&m));
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(k_,&k));
  PetscCall(PetscBLASIntCast(lda_,&lda));
  PetscCall(PetscBLASIntCast(ldb_,&ldb));
  PetscCall(PetscBLASIntCast(ldc_,&ldc));
#if defined(PETSC_HAVE_FBLASLAPACK) || defined(PETSC_HAVE_F2CBLASLAPACK)
  l = m % bs;
  if (l) PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&k,&alpha,(PetscScalar*)A,&lda,(PetscScalar*)B,&ldb,&beta,C,&ldc));
  for (;l<m;l+=bs) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&bs,&n,&k,&alpha,(PetscScalar*)A+l,&lda,(PetscScalar*)B,&ldb,&beta,C+l,&ldc));
  }
#else
  if (m) PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&m,&n,&k,&alpha,(PetscScalar*)A,&lda,(PetscScalar*)B,&ldb,&beta,C,&ldc));
#endif
  PetscCall(PetscLogFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    y := alpha*A*x + beta*y

    A is nxk (ld=lda)
*/
PetscErrorCode BVMultVec_BLAS_Private(BV bv,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *A,PetscInt lda_,const PetscScalar *x,PetscScalar beta,PetscScalar *y)
{
  PetscBLASInt   n,k,lda,one=1;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(k_,&k));
  PetscCall(PetscBLASIntCast(lda_,&lda));
  if (n) PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&k,&alpha,A,&lda,x,&one,&beta,y,&one));
  PetscCall(PetscLogFlops(2.0*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    A(:,s:e-1) := A*B(:,s:e-1)

    A is mxk (ld=lda), B is kxn (ld=ldb), n=e-s
*/
PetscErrorCode BVMultInPlace_BLAS_Private(BV bv,PetscInt m_,PetscInt k_,PetscInt s,PetscInt e,PetscScalar *A,PetscInt lda_,const PetscScalar *B,PetscInt ldb_,PetscBool btrans)
{
  PetscScalar    *pb,zero=0.0,one=1.0;
  PetscBLASInt   m,n,k,l,lda,ldb,bs=BLOCKSIZE;
  PetscInt       j,n_=e-s;
  const char     *bt;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(m_,&m));
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(k_,&k));
  PetscCall(PetscBLASIntCast(lda_,&lda));
  PetscCall(PetscBLASIntCast(ldb_,&ldb));
  PetscCall(BVAllocateWork_Private(bv,bs*n_));
  if (PetscUnlikely(btrans)) {
    pb = (PetscScalar*)B+s;
    bt = "C";
  } else {
    pb = (PetscScalar*)B+s*ldb;
    bt = "N";
  }
  l = m % bs;
  if (l) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N",bt,&l,&n,&k,&one,A,&lda,pb,&ldb,&zero,bv->work,&l));
    for (j=0;j<n;j++) PetscCall(PetscArraycpy(A+(s+j)*lda,bv->work+j*l,l));
  }
  for (;l<m;l+=bs) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N",bt,&bs,&n,&k,&one,A+l,&lda,pb,&ldb,&zero,bv->work,&bs));
    for (j=0;j<n;j++) PetscCall(PetscArraycpy(A+(s+j)*lda+l,bv->work+j*bs,bs));
  }
  PetscCall(PetscLogFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    V := V*B

    V is mxn (ld=m), B is nxn (ld=k)
*/
PetscErrorCode BVMultInPlace_Vecs_Private(BV bv,PetscInt m_,PetscInt n_,PetscInt k_,Vec *V,const PetscScalar *B,PetscBool btrans)
{
  PetscScalar       zero=0.0,one=1.0,*out,*pout;
  const PetscScalar *pin;
  PetscBLASInt      m = 0,n,k,l,bs=BLOCKSIZE;
  PetscInt          j;
  const char        *bt;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(m_,&m));
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(k_,&k));
  PetscCall(BVAllocateWork_Private(bv,2*bs*n_));
  out = bv->work+bs*n_;
  if (btrans) bt = "C";
  else bt = "N";
  l = m % bs;
  if (l) {
    for (j=0;j<n;j++) {
      PetscCall(VecGetArrayRead(V[j],&pin));
      PetscCall(PetscArraycpy(bv->work+j*l,pin,l));
      PetscCall(VecRestoreArrayRead(V[j],&pin));
    }
    PetscCallBLAS("BLASgemm",BLASgemm_("N",bt,&l,&n,&n,&one,bv->work,&l,(PetscScalar*)B,&k,&zero,out,&l));
    for (j=0;j<n;j++) {
      PetscCall(VecGetArray(V[j],&pout));
      PetscCall(PetscArraycpy(pout,out+j*l,l));
      PetscCall(VecRestoreArray(V[j],&pout));
    }
  }
  for (;l<m;l+=bs) {
    for (j=0;j<n;j++) {
      PetscCall(VecGetArrayRead(V[j],&pin));
      PetscCall(PetscArraycpy(bv->work+j*bs,pin+l,bs));
      PetscCall(VecRestoreArrayRead(V[j],&pin));
    }
    PetscCallBLAS("BLASgemm",BLASgemm_("N",bt,&bs,&n,&n,&one,bv->work,&bs,(PetscScalar*)B,&k,&zero,out,&bs));
    for (j=0;j<n;j++) {
      PetscCall(VecGetArray(V[j],&pout));
      PetscCall(PetscArraycpy(pout+l,out+j*bs,bs));
      PetscCall(VecRestoreArray(V[j],&pout));
    }
  }
  PetscCall(PetscLogFlops(2.0*n*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    B := alpha*A + beta*B

    A,B are nxk
*/
PetscErrorCode BVAXPY_BLAS_Private(BV bv,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *A,PetscInt lda_,PetscScalar beta,PetscScalar *B,PetscInt ldb_)
{
  PetscBLASInt   m,one=1;
  PetscInt       i,j;

  PetscFunctionBegin;
  if (lda_==n_ && ldb_==n_) {
    PetscCall(PetscBLASIntCast(n_*k_,&m));
    if (beta!=(PetscScalar)1.0) PetscCallBLAS("BLASscal",BLASscal_(&m,&beta,B,&one));
    PetscCallBLAS("BLASaxpy",BLASaxpy_(&m,&alpha,A,&one,B,&one));
  } else {
    if (beta!=(PetscScalar)1.0) {
      for (j=0;j<k_;j++) {
        for (i=0;i<n_;i++) {
          B[i+j*ldb_] = alpha*A[i+j*lda_] + beta*B[i+j*ldb_];
        }
      }
    } else {
      for (j=0;j<k_;j++) {
        for (i=0;i<n_;i++) {
          B[i+j*ldb_] += alpha*A[i+j*lda_];
        }
      }
    }
  }
  PetscCall(PetscLogFlops((beta==(PetscScalar)1.0)?2.0*n_*k_:3.0*n_*k_));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    C := A'*B

    A' is mxk (ld=lda), B is kxn (ld=ldb), C is mxn (ld=ldc)
*/
PetscErrorCode BVDot_BLAS_Private(BV bv,PetscInt m_,PetscInt n_,PetscInt k_,const PetscScalar *A,PetscInt lda_,const PetscScalar *B,PetscInt ldb_,PetscScalar *C,PetscInt ldc_,PetscBool mpi)
{
  PetscScalar    zero=0.0,one=1.0,*CC;
  PetscBLASInt   m,n,k,lda,ldb,ldc,j;
  PetscMPIInt    len;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(m_,&m));
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(k_,&k));
  PetscCall(PetscBLASIntCast(lda_,&lda));
  PetscCall(PetscBLASIntCast(ldb_,&ldb));
  PetscCall(PetscBLASIntCast(ldc_,&ldc));
  if (mpi) {
    if (ldc==m) {
      PetscCall(BVAllocateWork_Private(bv,m*n));
      if (k) PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&m,&n,&k,&one,(PetscScalar*)A,&lda,(PetscScalar*)B,&ldb,&zero,bv->work,&ldc));
      else PetscCall(PetscArrayzero(bv->work,m*n));
      PetscCall(PetscMPIIntCast(m*n,&len));
      PetscCallMPI(MPIU_Allreduce(bv->work,C,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
    } else {
      PetscCall(BVAllocateWork_Private(bv,2*m*n));
      CC = bv->work+m*n;
      if (k) PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&m,&n,&k,&one,(PetscScalar*)A,&lda,(PetscScalar*)B,&ldb,&zero,bv->work,&m));
      else PetscCall(PetscArrayzero(bv->work,m*n));
      PetscCall(PetscMPIIntCast(m*n,&len));
      PetscCallMPI(MPIU_Allreduce(bv->work,CC,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
      for (j=0;j<n;j++) PetscCall(PetscArraycpy(C+j*ldc,CC+j*m,m));
    }
  } else {
    if (k) PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&m,&n,&k,&one,(PetscScalar*)A,&lda,(PetscScalar*)B,&ldb,&zero,C,&ldc));
  }
  PetscCall(PetscLogFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    y := A'*x

    A is nxk (ld=lda)
*/
PetscErrorCode BVDotVec_BLAS_Private(BV bv,PetscInt n_,PetscInt k_,const PetscScalar *A,PetscInt lda_,const PetscScalar *x,PetscScalar *y,PetscBool mpi)
{
  PetscScalar    zero=0.0,done=1.0;
  PetscBLASInt   n,k,lda,one=1;
  PetscMPIInt    len;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(k_,&k));
  PetscCall(PetscBLASIntCast(lda_,&lda));
  if (mpi) {
    PetscCall(BVAllocateWork_Private(bv,k));
    if (n) PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&k,&done,A,&lda,x,&one,&zero,bv->work,&one));
    else PetscCall(PetscArrayzero(bv->work,k));
    PetscCall(PetscMPIIntCast(k,&len));
    PetscCallMPI(MPIU_Allreduce(bv->work,y,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
  } else {
    if (n) PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&k,&done,A,&lda,x,&one,&zero,y,&one));
  }
  PetscCall(PetscLogFlops(2.0*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Scale n scalars
*/
PetscErrorCode BVScale_BLAS_Private(BV bv,PetscInt n_,PetscScalar *A,PetscScalar alpha)
{
  PetscBLASInt   n,one=1;

  PetscFunctionBegin;
  if (PetscUnlikely(alpha == (PetscScalar)0.0)) PetscCall(PetscArrayzero(A,n_));
  else if (alpha!=(PetscScalar)1.0) {
    PetscCall(PetscBLASIntCast(n_,&n));
    PetscCallBLAS("BLASscal",BLASscal_(&n,&alpha,A,&one));
    PetscCall(PetscLogFlops(1.0*n));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
