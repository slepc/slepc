/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   CUDA-related code common to several BV impls
*/

#include <slepc/private/bvimpl.h>
#include <slepccupmblas.h>

#define BLOCKSIZE 64

/*
    C := alpha*A*B + beta*C
*/
PetscErrorCode BVMult_BLAS_CUDA(BV,PetscInt m_,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_B,PetscInt ldb_,PetscScalar beta,PetscScalar *d_C,PetscInt ldc_)
{
  PetscCuBLASInt    m=0,n=0,k=0,lda=0,ldb=0,ldc=0;
  cublasHandle_t    cublasv2handle;

  PetscFunctionBegin;
  PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
  PetscCall(PetscCuBLASIntCast(m_,&m));
  PetscCall(PetscCuBLASIntCast(n_,&n));
  PetscCall(PetscCuBLASIntCast(k_,&k));
  PetscCall(PetscCuBLASIntCast(lda_,&lda));
  PetscCall(PetscCuBLASIntCast(ldb_,&ldb));
  PetscCall(PetscCuBLASIntCast(ldc_,&ldc));
  PetscCall(PetscLogGpuTimeBegin());
  PetscCallCUBLAS(cublasXgemm(cublasv2handle,CUBLAS_OP_N,CUBLAS_OP_N,m,n,k,&alpha,d_A,lda,d_B,ldb,&beta,d_C,ldc));
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    y := alpha*A*x + beta*y
*/
PetscErrorCode BVMultVec_BLAS_CUDA(BV,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_x,PetscScalar beta,PetscScalar *d_y)
{
  PetscCuBLASInt    n=0,k=0,lda=0,one=1;
  cublasHandle_t    cublasv2handle;

  PetscFunctionBegin;
  PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
  PetscCall(PetscCuBLASIntCast(n_,&n));
  PetscCall(PetscCuBLASIntCast(k_,&k));
  PetscCall(PetscCuBLASIntCast(lda_,&lda));
  PetscCall(PetscLogGpuTimeBegin());
  PetscCallCUBLAS(cublasXgemv(cublasv2handle,CUBLAS_OP_N,n,k,&alpha,d_A,lda,d_x,one,&beta,d_y,one));
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops(2.0*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    A(:,s:e-1) := A*B(:,s:e-1)
*/
PetscErrorCode BVMultInPlace_BLAS_CUDA(BV,PetscInt m_,PetscInt k_,PetscInt s,PetscInt e,PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_B,PetscInt ldb_,PetscBool btrans)
{
  const PetscScalar *d_B1;
  PetscScalar       *d_work,sone=1.0,szero=0.0;
  PetscCuBLASInt    m=0,n=0,k=0,l=0,lda=0,ldb=0,bs=BLOCKSIZE;
  size_t            freemem,totmem;
  cublasHandle_t    cublasv2handle;
  cublasOperation_t bt;

  PetscFunctionBegin;
  PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
  PetscCall(PetscCuBLASIntCast(m_,&m));
  PetscCall(PetscCuBLASIntCast(e-s,&n));
  PetscCall(PetscCuBLASIntCast(k_,&k));
  PetscCall(PetscCuBLASIntCast(lda_,&lda));
  PetscCall(PetscCuBLASIntCast(ldb_,&ldb));
  PetscCall(PetscLogGpuTimeBegin());
  if (PetscUnlikely(btrans)) {
    d_B1 = d_B+s;
    bt   = CUBLAS_OP_C;
  } else {
    d_B1 = d_B+s*ldb;
    bt   = CUBLAS_OP_N;
  }
  /* try to allocate the whole matrix */
  PetscCallCUDA(cudaMemGetInfo(&freemem,&totmem));
  if (freemem>=lda*n*sizeof(PetscScalar)) {
    PetscCallCUDA(cudaMalloc((void**)&d_work,lda*n*sizeof(PetscScalar)));
    PetscCallCUBLAS(cublasXgemm(cublasv2handle,CUBLAS_OP_N,bt,m,n,k,&sone,d_A,lda,d_B1,ldb,&szero,d_work,lda));
    PetscCallCUDA(cudaMemcpy2D(d_A+s*lda,lda*sizeof(PetscScalar),d_work,lda*sizeof(PetscScalar),m*sizeof(PetscScalar),n,cudaMemcpyDeviceToDevice));
  } else {
    PetscCall(PetscCuBLASIntCast(freemem/(m*sizeof(PetscScalar)),&bs));
    PetscCallCUDA(cudaMalloc((void**)&d_work,bs*n*sizeof(PetscScalar)));
    PetscCall(PetscCuBLASIntCast(m % bs,&l));
    if (l) {
      PetscCallCUBLAS(cublasXgemm(cublasv2handle,CUBLAS_OP_N,bt,l,n,k,&sone,d_A,lda,d_B1,ldb,&szero,d_work,l));
      PetscCallCUDA(cudaMemcpy2D(d_A+s*lda,lda*sizeof(PetscScalar),d_work,l*sizeof(PetscScalar),l*sizeof(PetscScalar),n,cudaMemcpyDeviceToDevice));
    }
    for (;l<m;l+=bs) {
      PetscCallCUBLAS(cublasXgemm(cublasv2handle,CUBLAS_OP_N,bt,bs,n,k,&sone,d_A+l,lda,d_B1,ldb,&szero,d_work,bs));
      PetscCallCUDA(cudaMemcpy2D(d_A+l+s*lda,lda*sizeof(PetscScalar),d_work,bs*sizeof(PetscScalar),bs*sizeof(PetscScalar),n,cudaMemcpyDeviceToDevice));
    }
  }
  PetscCall(PetscLogGpuTimeEnd());
  PetscCallCUDA(cudaFree(d_work));
  PetscCall(PetscLogGpuFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    B := alpha*A + beta*B
*/
PetscErrorCode BVAXPY_BLAS_CUDA(BV,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *d_A,PetscInt lda_,PetscScalar beta,PetscScalar *d_B,PetscInt ldb_)
{
  PetscCuBLASInt n=0,k=0,lda=0,ldb=0;
  cublasHandle_t cublasv2handle;

  PetscFunctionBegin;
  PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
  PetscCall(PetscCuBLASIntCast(n_,&n));
  PetscCall(PetscCuBLASIntCast(k_,&k));
  PetscCall(PetscCuBLASIntCast(lda_,&lda));
  PetscCall(PetscCuBLASIntCast(ldb_,&ldb));
  PetscCall(PetscLogGpuTimeBegin());
  PetscCallCUBLAS(cublasXgeam(cublasv2handle,CUBLAS_OP_N,CUBLAS_OP_N,n,k,&alpha,d_A,lda,&beta,d_B,ldb,d_B,ldb));
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops((beta==(PetscScalar)1.0)?2.0*n*k:3.0*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    C := A'*B

    C is a CPU array
*/
PetscErrorCode BVDot_BLAS_CUDA(BV bv,PetscInt m_,PetscInt n_,PetscInt k_,const PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_B,PetscInt ldb_,PetscScalar *C,PetscInt ldc_,PetscBool mpi)
{
  PetscScalar       *d_work,sone=1.0,szero=0.0,*CC;
  PetscInt          j;
  PetscCuBLASInt    m=0,n=0,k=0,lda=0,ldb=0,ldc=0;
  PetscMPIInt       len;
  cublasHandle_t    cublasv2handle;

  PetscFunctionBegin;
  PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
  PetscCall(PetscCuBLASIntCast(m_,&m));
  PetscCall(PetscCuBLASIntCast(n_,&n));
  PetscCall(PetscCuBLASIntCast(k_,&k));
  PetscCall(PetscCuBLASIntCast(lda_,&lda));
  PetscCall(PetscCuBLASIntCast(ldb_,&ldb));
  PetscCall(PetscCuBLASIntCast(ldc_,&ldc));
  PetscCallCUDA(cudaMalloc((void**)&d_work,m*n*sizeof(PetscScalar)));
  if (mpi) {
    if (ldc==m) {
      PetscCall(BVAllocateWork_Private(bv,m*n));
      if (k) {
        PetscCall(PetscLogGpuTimeBegin());
        PetscCallCUBLAS(cublasXgemm(cublasv2handle,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&sone,d_A,lda,d_B,ldb,&szero,d_work,ldc));
        PetscCall(PetscLogGpuTimeEnd());
        PetscCallCUDA(cudaMemcpy(bv->work,d_work,m*n*sizeof(PetscScalar),cudaMemcpyDeviceToHost));
        PetscCall(PetscLogGpuToCpu(m*n*sizeof(PetscScalar)));
      } else PetscCall(PetscArrayzero(bv->work,m*n));
      PetscCall(PetscMPIIntCast(m*n,&len));
      PetscCallMPI(MPIU_Allreduce(bv->work,C,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
    } else {
      PetscCall(BVAllocateWork_Private(bv,2*m*n));
      CC = bv->work+m*n;
      if (k) {
        PetscCall(PetscLogGpuTimeBegin());
        PetscCallCUBLAS(cublasXgemm(cublasv2handle,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&sone,d_A,lda,d_B,ldb,&szero,d_work,m));
        PetscCall(PetscLogGpuTimeEnd());
        PetscCallCUDA(cudaMemcpy(bv->work,d_work,m*n*sizeof(PetscScalar),cudaMemcpyDeviceToHost));
        PetscCall(PetscLogGpuToCpu(m*n*sizeof(PetscScalar)));
      } else PetscCall(PetscArrayzero(bv->work,m*n));
      PetscCall(PetscMPIIntCast(m*n,&len));
      PetscCallMPI(MPIU_Allreduce(bv->work,CC,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
      for (j=0;j<n;j++) PetscCall(PetscArraycpy(C+j*ldc,CC+j*m,m));
    }
  } else {
    if (k) {
      PetscCall(BVAllocateWork_Private(bv,m*n));
      PetscCall(PetscLogGpuTimeBegin());
      PetscCallCUBLAS(cublasXgemm(cublasv2handle,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&sone,d_A,lda,d_B,ldb,&szero,d_work,m));
      PetscCall(PetscLogGpuTimeEnd());
      PetscCallCUDA(cudaMemcpy(bv->work,d_work,m*n*sizeof(PetscScalar),cudaMemcpyDeviceToHost));
      PetscCall(PetscLogGpuToCpu(m*n*sizeof(PetscScalar)));
      for (j=0;j<n;j++) PetscCall(PetscArraycpy(C+j*ldc,bv->work+j*m,m));
    }
  }
  PetscCallCUDA(cudaFree(d_work));
  PetscCall(PetscLogGpuFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    y := A'*x

    y is a CPU array, if NULL bv->buffer is used as a workspace
*/
PetscErrorCode BVDotVec_BLAS_CUDA(BV bv,PetscInt n_,PetscInt k_,const PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_x,PetscScalar *y,PetscBool mpi)
{
  PetscScalar       *d_work,szero=0.0,sone=1.0,*yy;
  PetscCuBLASInt    n=0,k=0,lda=0,one=1;
  PetscMPIInt       len;
  cublasHandle_t    cublasv2handle;

  PetscFunctionBegin;
  PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
  PetscCall(PetscCuBLASIntCast(n_,&n));
  PetscCall(PetscCuBLASIntCast(k_,&k));
  PetscCall(PetscCuBLASIntCast(lda_,&lda));
  if (!y) PetscCall(VecCUDAGetArrayWrite(bv->buffer,&d_work));
  else PetscCallCUDA(cudaMalloc((void**)&d_work,k*sizeof(PetscScalar)));
  if (mpi) {
    PetscCall(BVAllocateWork_Private(bv,k));
    if (n) {
      PetscCall(PetscLogGpuTimeBegin());
      PetscCallCUBLAS(cublasXgemv(cublasv2handle,CUBLAS_OP_C,n,k,&sone,d_A,lda,d_x,one,&szero,d_work,one));
      PetscCall(PetscLogGpuTimeEnd());
      PetscCallCUDA(cudaMemcpy(bv->work,d_work,k*sizeof(PetscScalar),cudaMemcpyDeviceToHost));
      PetscCall(PetscLogGpuToCpu(k*sizeof(PetscScalar)));
    } else PetscCall(PetscArrayzero(bv->work,k));
    /* reduction */
    PetscCall(PetscMPIIntCast(k,&len));
    if (!y) {
      if (use_gpu_aware_mpi) {  /* case 1: reduce on GPU using a temporary buffer */
        PetscCallCUDA(cudaMalloc((void**)&yy,k*sizeof(PetscScalar)));
        PetscCallMPI(MPIU_Allreduce(d_work,yy,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
        PetscCallCUDA(cudaMemcpy(d_work,yy,k*sizeof(PetscScalar),cudaMemcpyDeviceToDevice));
        PetscCallCUDA(cudaFree(yy));
      } else {  /* case 2: reduce on CPU, copy result back to GPU */
        PetscCall(BVAllocateWork_Private(bv,2*k));
        yy = bv->work+k;
        PetscCallCUDA(cudaMemcpy(bv->work,d_work,k*sizeof(PetscScalar),cudaMemcpyDeviceToHost));
        PetscCall(PetscLogGpuToCpu(k*sizeof(PetscScalar)));
        PetscCallMPI(MPIU_Allreduce(bv->work,yy,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
        PetscCallCUDA(cudaMemcpy(d_work,yy,k*sizeof(PetscScalar),cudaMemcpyHostToDevice));
        PetscCall(PetscLogCpuToGpu(k*sizeof(PetscScalar)));
      }
      PetscCall(VecCUDARestoreArrayWrite(bv->buffer,&d_work));
    } else {  /* case 3: user-provided array y, reduce on CPU */
      PetscCallCUDA(cudaFree(d_work));
      PetscCallMPI(MPIU_Allreduce(bv->work,y,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
    }
  } else {
    if (n) {
      PetscCall(PetscLogGpuTimeBegin());
      PetscCallCUBLAS(cublasXgemv(cublasv2handle,CUBLAS_OP_C,n,k,&sone,d_A,lda,d_x,one,&szero,d_work,one));
      PetscCall(PetscLogGpuTimeEnd());
    }
    if (!y) PetscCall(VecCUDARestoreArrayWrite(bv->buffer,&d_work));
    else {
      PetscCallCUDA(cudaMemcpy(y,d_work,k*sizeof(PetscScalar),cudaMemcpyDeviceToHost));
      PetscCall(PetscLogGpuToCpu(k*sizeof(PetscScalar)));
      PetscCallCUDA(cudaFree(d_work));
    }
  }
  PetscCall(PetscLogGpuFlops(2.0*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Scale n scalars
*/
PetscErrorCode BVScale_BLAS_CUDA(BV,PetscInt n_,PetscScalar *d_A,PetscScalar alpha)
{
  PetscCuBLASInt n=0,one=1;
  cublasHandle_t cublasv2handle;

  PetscFunctionBegin;
  PetscCall(PetscCuBLASIntCast(n_,&n));
  if (PetscUnlikely(alpha == (PetscScalar)0.0)) PetscCallCUDA(cudaMemset(d_A,0,n*sizeof(PetscScalar)));
  else if (alpha != (PetscScalar)1.0) {
    PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallCUBLAS(cublasXscal(cublasv2handle,n,&alpha,d_A,one));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(PetscLogGpuFlops(1.0*n));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Compute 2-norm of vector consisting of n scalars
*/
PetscErrorCode BVNorm_BLAS_CUDA(BV,PetscInt n_,const PetscScalar *d_A,PetscReal *nrm)
{
  PetscCuBLASInt n=0,one=1;
  cublasHandle_t cublasv2handle;

  PetscFunctionBegin;
  PetscCall(PetscCuBLASIntCast(n_,&n));
  PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
  PetscCall(PetscLogGpuTimeBegin());
  PetscCallCUBLAS(cublasXnrm2(cublasv2handle,n,d_A,one,nrm));
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops(2.0*n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Normalize the columns of A
*/
PetscErrorCode BVNormalize_BLAS_CUDA(BV,PetscInt m_,PetscInt n_,PetscScalar *d_A,PetscInt lda_,PetscScalar *eigi)
{
  PetscInt       j,k;
  PetscReal      nrm,nrm1;
  PetscScalar    alpha;
  PetscCuBLASInt m=0,one=1;
  cublasHandle_t cublasv2handle;

  PetscFunctionBegin;
  PetscCall(PetscCuBLASIntCast(m_,&m));
  PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
  PetscCall(PetscLogGpuTimeBegin());
  for (j=0;j<n_;j++) {
    k = 1;
#if !defined(PETSC_USE_COMPLEX)
    if (eigi && eigi[j] != 0.0) k = 2;
#endif
    PetscCallCUBLAS(cublasXnrm2(cublasv2handle,m,d_A+j*lda_,one,&nrm));
    if (k==2) {
      PetscCallCUBLAS(cublasXnrm2(cublasv2handle,m,d_A+(j+1)*lda_,one,&nrm1));
      nrm = SlepcAbs(nrm,nrm1);
    }
    alpha = 1.0/nrm;
    PetscCallCUBLAS(cublasXscal(cublasv2handle,m,&alpha,d_A+j*lda_,one));
    if (k==2) {
      PetscCallCUBLAS(cublasXscal(cublasv2handle,m,&alpha,d_A+(j+1)*lda_,one));
      j++;
    }
  }
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops(3.0*m*n_));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_CleanCoefficients_CUDA - Sets to zero all entries of column j of the bv buffer
*/
PetscErrorCode BV_CleanCoefficients_CUDA(BV bv,PetscInt j,PetscScalar *h)
{
  PetscScalar    *d_hh,*d_a;
  PetscInt       i;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecCUDAGetArray(bv->buffer,&d_a));
    PetscCall(PetscLogGpuTimeBegin());
    d_hh = d_a + j*(bv->nc+bv->m);
    PetscCallCUDA(cudaMemset(d_hh,0,(bv->nc+j)*sizeof(PetscScalar)));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(VecCUDARestoreArray(bv->buffer,&d_a));
  } else { /* cpu memory */
    for (i=0;i<bv->nc+j;i++) h[i] = 0.0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_AddCoefficients_CUDA - Add the contents of the scratch (0-th column) of the bv buffer
   into column j of the bv buffer
 */
PetscErrorCode BV_AddCoefficients_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscScalar *c)
{
  PetscScalar    *d_h,*d_c,sone=1.0;
  PetscInt       i;
  PetscCuBLASInt idx=0,one=1;
  cublasHandle_t cublasv2handle;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
    PetscCall(VecCUDAGetArray(bv->buffer,&d_c));
    d_h = d_c + j*(bv->nc+bv->m);
    PetscCall(PetscCuBLASIntCast(bv->nc+j,&idx));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallCUBLAS(cublasXaxpy(cublasv2handle,idx,&sone,d_c,one,d_h,one));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(PetscLogGpuFlops(1.0*(bv->nc+j)));
    PetscCall(VecCUDARestoreArray(bv->buffer,&d_c));
  } else { /* cpu memory */
    for (i=0;i<bv->nc+j;i++) h[i] += c[i];
    PetscCall(PetscLogFlops(1.0*(bv->nc+j)));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SetValue_CUDA - Sets value in row j (counted after the constraints) of column k
   of the coefficients array
*/
PetscErrorCode BV_SetValue_CUDA(BV bv,PetscInt j,PetscInt k,PetscScalar *h,PetscScalar value)
{
  PetscScalar    *d_h,*a;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecCUDAGetArray(bv->buffer,&a));
    PetscCall(PetscLogGpuTimeBegin());
    d_h = a + k*(bv->nc+bv->m) + bv->nc+j;
    PetscCallCUDA(cudaMemcpy(d_h,&value,sizeof(PetscScalar),cudaMemcpyHostToDevice));
    PetscCall(PetscLogCpuToGpu(sizeof(PetscScalar)));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(VecCUDARestoreArray(bv->buffer,&a));
  } else { /* cpu memory */
    h[bv->nc+j] = value;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SquareSum_CUDA - Returns the value h'*h, where h represents the contents of the
   coefficients array (up to position j)
*/
PetscErrorCode BV_SquareSum_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscReal *sum)
{
  const PetscScalar *d_h;
  PetscScalar       dot;
  PetscInt          i;
  PetscCuBLASInt    idx=0,one=1;
  cublasHandle_t    cublasv2handle;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(PetscCUBLASGetHandle(&cublasv2handle));
    PetscCall(VecCUDAGetArrayRead(bv->buffer,&d_h));
    PetscCall(PetscCuBLASIntCast(bv->nc+j,&idx));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallCUBLAS(cublasXdot(cublasv2handle,idx,d_h,one,d_h,one,&dot));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(PetscLogGpuFlops(2.0*(bv->nc+j)));
    *sum = PetscRealPart(dot);
    PetscCall(VecCUDARestoreArrayRead(bv->buffer,&d_h));
  } else { /* cpu memory */
    *sum = 0.0;
    for (i=0;i<bv->nc+j;i++) *sum += PetscRealPart(h[i]*PetscConj(h[i]));
    PetscCall(PetscLogFlops(2.0*(bv->nc+j)));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* pointwise multiplication */
static __global__ void PointwiseMult_kernel(PetscInt xcount,PetscScalar *a,const PetscScalar *b,PetscInt n)
{
  PetscInt x;

  x = xcount*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x;
  if (x<n) a[x] *= PetscRealPart(b[x]);
}

/* pointwise division */
static __global__ void PointwiseDiv_kernel(PetscInt xcount,PetscScalar *a,const PetscScalar *b,PetscInt n)
{
  PetscInt x;

  x = xcount*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x;
  if (x<n) a[x] /= PetscRealPart(b[x]);
}

/*
   BV_ApplySignature_CUDA - Computes the pointwise product h*omega, where h represents
   the contents of the coefficients array (up to position j) and omega is the signature;
   if inverse=TRUE then the operation is h/omega
*/
PetscErrorCode BV_ApplySignature_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscBool inverse)
{
  PetscScalar       *d_h;
  const PetscScalar *d_omega,*omega;
  PetscInt          i,xcount;
  dim3              blocks3d, threads3d;

  PetscFunctionBegin;
  if (!(bv->nc+j)) PetscFunctionReturn(PETSC_SUCCESS);
  if (!h) {
    PetscCall(VecCUDAGetArray(bv->buffer,&d_h));
    PetscCall(VecCUDAGetArrayRead(bv->omega,&d_omega));
    PetscCall(SlepcKernelSetGrid1D(bv->nc+j,&blocks3d,&threads3d,&xcount));
    PetscCall(PetscLogGpuTimeBegin());
    if (inverse) {
      for (i=0;i<xcount;i++) PointwiseDiv_kernel<<<blocks3d,threads3d>>>(i,d_h,d_omega,bv->nc+j);
    } else {
      for (i=0;i<xcount;i++) PointwiseMult_kernel<<<blocks3d,threads3d>>>(i,d_h,d_omega,bv->nc+j);
    }
    PetscCallCUDA(cudaGetLastError());
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(PetscLogGpuFlops(1.0*(bv->nc+j)));
    PetscCall(VecCUDARestoreArrayRead(bv->omega,&d_omega));
    PetscCall(VecCUDARestoreArray(bv->buffer,&d_h));
  } else {
    PetscCall(VecGetArrayRead(bv->omega,&omega));
    if (inverse) for (i=0;i<bv->nc+j;i++) h[i] /= PetscRealPart(omega[i]);
    else for (i=0;i<bv->nc+j;i++) h[i] *= PetscRealPart(omega[i]);
    PetscCall(VecRestoreArrayRead(bv->omega,&omega));
    PetscCall(PetscLogFlops(1.0*(bv->nc+j)));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SquareRoot_CUDA - Returns the square root of position j (counted after the constraints)
   of the coefficients array
*/
PetscErrorCode BV_SquareRoot_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscReal *beta)
{
  const PetscScalar *d_h;
  PetscScalar       hh;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecCUDAGetArrayRead(bv->buffer,&d_h));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallCUDA(cudaMemcpy(&hh,d_h+bv->nc+j,sizeof(PetscScalar),cudaMemcpyDeviceToHost));
    PetscCall(PetscLogGpuToCpu(sizeof(PetscScalar)));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(BV_SafeSqrt(bv,hh,beta));
    PetscCall(VecCUDARestoreArrayRead(bv->buffer,&d_h));
  } else PetscCall(BV_SafeSqrt(bv,h[bv->nc+j],beta));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_StoreCoefficients_CUDA - Copy the contents of the coefficients array to an array dest
   provided by the caller (only values from l to j are copied)
*/
PetscErrorCode BV_StoreCoefficients_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscScalar *dest)
{
  const PetscScalar *d_h,*d_a;
  PetscInt          i;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecCUDAGetArrayRead(bv->buffer,&d_a));
    PetscCall(PetscLogGpuTimeBegin());
    d_h = d_a + j*(bv->nc+bv->m)+bv->nc;
    PetscCallCUDA(cudaMemcpy(dest-bv->l,d_h,(j-bv->l)*sizeof(PetscScalar),cudaMemcpyDeviceToHost));
    PetscCall(PetscLogGpuToCpu((j-bv->l)*sizeof(PetscScalar)));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(VecCUDARestoreArrayRead(bv->buffer,&d_a));
  } else {
    for (i=bv->l;i<j;i++) dest[i-bv->l] = h[bv->nc+i];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
