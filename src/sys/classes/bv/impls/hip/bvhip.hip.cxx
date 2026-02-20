/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   HIP-related code common to several BV impls
*/

#include <slepc/private/bvimpl.h>
#include <slepccupmblas.h>

#define BLOCKSIZE 64

/*
    C := alpha*A*B + beta*C
*/
PetscErrorCode BVMult_BLAS_HIP(BV,PetscInt m_,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_B,PetscInt ldb_,PetscScalar beta,PetscScalar *d_C,PetscInt ldc_)
{
  PetscHipBLASInt    m=0,n=0,k=0,lda=0,ldb=0,ldc=0;
  hipblasHandle_t    hipblashandle;

  PetscFunctionBegin;
  PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
  PetscCall(PetscHipBLASIntCast(m_,&m));
  PetscCall(PetscHipBLASIntCast(n_,&n));
  PetscCall(PetscHipBLASIntCast(k_,&k));
  PetscCall(PetscHipBLASIntCast(lda_,&lda));
  PetscCall(PetscHipBLASIntCast(ldb_,&ldb));
  PetscCall(PetscHipBLASIntCast(ldc_,&ldc));
  PetscCall(PetscLogGpuTimeBegin());
  PetscCallHIPBLAS(hipblasXgemm(hipblashandle,HIPBLAS_OP_N,HIPBLAS_OP_N,m,n,k,&alpha,d_A,lda,d_B,ldb,&beta,d_C,ldc));
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    y := alpha*A*x + beta*y
*/
PetscErrorCode BVMultVec_BLAS_HIP(BV,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_x,PetscScalar beta,PetscScalar *d_y)
{
  PetscHipBLASInt    n=0,k=0,lda=0,one=1;
  hipblasHandle_t    hipblashandle;

  PetscFunctionBegin;
  PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
  PetscCall(PetscHipBLASIntCast(n_,&n));
  PetscCall(PetscHipBLASIntCast(k_,&k));
  PetscCall(PetscHipBLASIntCast(lda_,&lda));
  PetscCall(PetscLogGpuTimeBegin());
  PetscCallHIPBLAS(hipblasXgemv(hipblashandle,HIPBLAS_OP_N,n,k,&alpha,d_A,lda,d_x,one,&beta,d_y,one));
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops(2.0*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    A(:,s:e-1) := A*B(:,s:e-1)
*/
PetscErrorCode BVMultInPlace_BLAS_HIP(BV,PetscInt m_,PetscInt k_,PetscInt s,PetscInt e,PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_B,PetscInt ldb_,PetscBool btrans)
{
  const PetscScalar  *d_B1;
  PetscScalar        *d_work,sone=1.0,szero=0.0;
  PetscHipBLASInt    m=0,n=0,k=0,l=0,lda=0,ldb=0,bs=BLOCKSIZE;
  size_t             freemem,totmem;
  hipblasHandle_t    hipblashandle;
  hipblasOperation_t bt;

  PetscFunctionBegin;
  PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
  PetscCall(PetscHipBLASIntCast(m_,&m));
  PetscCall(PetscHipBLASIntCast(e-s,&n));
  PetscCall(PetscHipBLASIntCast(k_,&k));
  PetscCall(PetscHipBLASIntCast(lda_,&lda));
  PetscCall(PetscHipBLASIntCast(ldb_,&ldb));
  PetscCall(PetscLogGpuTimeBegin());
  if (PetscUnlikely(btrans)) {
    d_B1 = d_B+s;
    bt   = HIPBLAS_OP_C;
  } else {
    d_B1 = d_B+s*ldb;
    bt   = HIPBLAS_OP_N;
  }
  /* try to allocate the whole matrix */
  PetscCallHIP(hipMemGetInfo(&freemem,&totmem));
  if (freemem>=lda*n*sizeof(PetscScalar)) {
    PetscCallHIP(hipMalloc((void**)&d_work,lda*n*sizeof(PetscScalar)));
    PetscCallHIPBLAS(hipblasXgemm(hipblashandle,HIPBLAS_OP_N,bt,m,n,k,&sone,d_A,lda,d_B1,ldb,&szero,d_work,lda));
    PetscCallHIP(hipMemcpy2D(d_A+s*lda,lda*sizeof(PetscScalar),d_work,lda*sizeof(PetscScalar),m*sizeof(PetscScalar),n,hipMemcpyDeviceToDevice));
  } else {
    PetscCall(PetscHipBLASIntCast(freemem/(m*sizeof(PetscScalar)),&bs));
    PetscCallHIP(hipMalloc((void**)&d_work,bs*n*sizeof(PetscScalar)));
    PetscCall(PetscHipBLASIntCast(m % bs,&l));
    if (l) {
      PetscCallHIPBLAS(hipblasXgemm(hipblashandle,HIPBLAS_OP_N,bt,l,n,k,&sone,d_A,lda,d_B1,ldb,&szero,d_work,l));
      PetscCallHIP(hipMemcpy2D(d_A+s*lda,lda*sizeof(PetscScalar),d_work,l*sizeof(PetscScalar),l*sizeof(PetscScalar),n,hipMemcpyDeviceToDevice));
    }
    for (;l<m;l+=bs) {
      PetscCallHIPBLAS(hipblasXgemm(hipblashandle,HIPBLAS_OP_N,bt,bs,n,k,&sone,d_A+l,lda,d_B1,ldb,&szero,d_work,bs));
      PetscCallHIP(hipMemcpy2D(d_A+l+s*lda,lda*sizeof(PetscScalar),d_work,bs*sizeof(PetscScalar),bs*sizeof(PetscScalar),n,hipMemcpyDeviceToDevice));
    }
  }
  PetscCall(PetscLogGpuTimeEnd());
  PetscCallHIP(hipFree(d_work));
  PetscCall(PetscLogGpuFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    B := alpha*A + beta*B
*/
PetscErrorCode BVAXPY_BLAS_HIP(BV,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *d_A,PetscInt lda_,PetscScalar beta,PetscScalar *d_B,PetscInt ldb_)
{
  PetscHipBLASInt n=0,k=0,lda=0,ldb=0;
  hipblasHandle_t hipblashandle;

  PetscFunctionBegin;
  PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
  PetscCall(PetscHipBLASIntCast(n_,&n));
  PetscCall(PetscHipBLASIntCast(k_,&k));
  PetscCall(PetscHipBLASIntCast(lda_,&lda));
  PetscCall(PetscHipBLASIntCast(ldb_,&ldb));
  PetscCall(PetscLogGpuTimeBegin());
  PetscCallHIPBLAS(hipblasXgeam(hipblashandle,HIPBLAS_OP_N,HIPBLAS_OP_N,n,k,&alpha,d_A,lda,&beta,d_B,ldb,d_B,ldb));
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops((beta==(PetscScalar)1.0)?2.0*n*k:3.0*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    C := A'*B

    C is a CPU array
*/
PetscErrorCode BVDot_BLAS_HIP(BV bv,PetscInt m_,PetscInt n_,PetscInt k_,const PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_B,PetscInt ldb_,PetscScalar *C,PetscInt ldc_,PetscBool mpi)
{
  PetscScalar       *d_work,sone=1.0,szero=0.0,*CC;
  PetscInt          j;
  PetscHipBLASInt   m=0,n=0,k=0,lda=0,ldb=0,ldc=0;
  PetscMPIInt       len;
  hipblasHandle_t   hipblashandle;

  PetscFunctionBegin;
  PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
  PetscCall(PetscHipBLASIntCast(m_,&m));
  PetscCall(PetscHipBLASIntCast(n_,&n));
  PetscCall(PetscHipBLASIntCast(k_,&k));
  PetscCall(PetscHipBLASIntCast(lda_,&lda));
  PetscCall(PetscHipBLASIntCast(ldb_,&ldb));
  PetscCall(PetscHipBLASIntCast(ldc_,&ldc));
  PetscCallHIP(hipMalloc((void**)&d_work,m*n*sizeof(PetscScalar)));
  if (mpi) {
    if (ldc==m) {
      PetscCall(BVAllocateWork_Private(bv,m*n));
      if (k) {
        PetscCall(PetscLogGpuTimeBegin());
        PetscCallHIPBLAS(hipblasXgemm(hipblashandle,HIPBLAS_OP_C,HIPBLAS_OP_N,m,n,k,&sone,d_A,lda,d_B,ldb,&szero,d_work,ldc));
        PetscCall(PetscLogGpuTimeEnd());
        PetscCallHIP(hipMemcpy(bv->work,d_work,m*n*sizeof(PetscScalar),hipMemcpyDeviceToHost));
        PetscCall(PetscLogGpuToCpu(m*n*sizeof(PetscScalar)));
      } else PetscCall(PetscArrayzero(bv->work,m*n));
      PetscCall(PetscMPIIntCast(m*n,&len));
      PetscCallMPI(MPIU_Allreduce(bv->work,C,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
    } else {
      PetscCall(BVAllocateWork_Private(bv,2*m*n));
      CC = bv->work+m*n;
      if (k) {
        PetscCall(PetscLogGpuTimeBegin());
        PetscCallHIPBLAS(hipblasXgemm(hipblashandle,HIPBLAS_OP_C,HIPBLAS_OP_N,m,n,k,&sone,d_A,lda,d_B,ldb,&szero,d_work,m));
        PetscCall(PetscLogGpuTimeEnd());
        PetscCallHIP(hipMemcpy(bv->work,d_work,m*n*sizeof(PetscScalar),hipMemcpyDeviceToHost));
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
      PetscCallHIPBLAS(hipblasXgemm(hipblashandle,HIPBLAS_OP_C,HIPBLAS_OP_N,m,n,k,&sone,d_A,lda,d_B,ldb,&szero,d_work,m));
      PetscCall(PetscLogGpuTimeEnd());
      PetscCallHIP(hipMemcpy(bv->work,d_work,m*n*sizeof(PetscScalar),hipMemcpyDeviceToHost));
      PetscCall(PetscLogGpuToCpu(m*n*sizeof(PetscScalar)));
      for (j=0;j<n;j++) PetscCall(PetscArraycpy(C+j*ldc,bv->work+j*m,m));
    }
  }
  PetscCallHIP(hipFree(d_work));
  PetscCall(PetscLogGpuFlops(2.0*m*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    y := A'*x

    y is a CPU array, if NULL bv->buffer is used as a workspace
*/
PetscErrorCode BVDotVec_BLAS_HIP(BV bv,PetscInt n_,PetscInt k_,const PetscScalar *d_A,PetscInt lda_,const PetscScalar *d_x,PetscScalar *y,PetscBool mpi)
{
  PetscScalar       *d_work,szero=0.0,sone=1.0,*yy;
  PetscHipBLASInt   n=0,k=0,lda=0,one=1;
  PetscMPIInt       len;
  hipblasHandle_t   hipblashandle;

  PetscFunctionBegin;
  PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
  PetscCall(PetscHipBLASIntCast(n_,&n));
  PetscCall(PetscHipBLASIntCast(k_,&k));
  PetscCall(PetscHipBLASIntCast(lda_,&lda));
  if (!y) PetscCall(VecHIPGetArrayWrite(bv->buffer,&d_work));
  else PetscCallHIP(hipMalloc((void**)&d_work,k*sizeof(PetscScalar)));
  if (mpi) {
    PetscCall(BVAllocateWork_Private(bv,k));
    if (n) {
      PetscCall(PetscLogGpuTimeBegin());
      PetscCallHIPBLAS(hipblasXgemv(hipblashandle,HIPBLAS_OP_C,n,k,&sone,d_A,lda,d_x,one,&szero,d_work,one));
      PetscCall(PetscLogGpuTimeEnd());
      PetscCallHIP(hipMemcpy(bv->work,d_work,k*sizeof(PetscScalar),hipMemcpyDeviceToHost));
      PetscCall(PetscLogGpuToCpu(k*sizeof(PetscScalar)));
    } else PetscCall(PetscArrayzero(bv->work,k));
    /* reduction */
    PetscCall(PetscMPIIntCast(k,&len));
    if (!y) {
      if (use_gpu_aware_mpi) {  /* case 1: reduce on GPU using a temporary buffer */
        PetscCallHIP(hipMalloc((void**)&yy,k*sizeof(PetscScalar)));
        PetscCallMPI(MPIU_Allreduce(d_work,yy,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
        PetscCallHIP(hipMemcpy(d_work,yy,k*sizeof(PetscScalar),hipMemcpyDeviceToDevice));
        PetscCallHIP(hipFree(yy));
      } else {  /* case 2: reduce on CPU, copy result back to GPU */
        PetscCall(BVAllocateWork_Private(bv,2*k));
        yy = bv->work+k;
        PetscCallHIP(hipMemcpy(bv->work,d_work,k*sizeof(PetscScalar),hipMemcpyDeviceToHost));
        PetscCall(PetscLogGpuToCpu(k*sizeof(PetscScalar)));
        PetscCallMPI(MPIU_Allreduce(bv->work,yy,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
        PetscCallHIP(hipMemcpy(d_work,yy,k*sizeof(PetscScalar),hipMemcpyHostToDevice));
        PetscCall(PetscLogCpuToGpu(k*sizeof(PetscScalar)));
      }
      PetscCall(VecHIPRestoreArrayWrite(bv->buffer,&d_work));
    } else {  /* case 3: user-provided array y, reduce on CPU */
      PetscCallHIP(hipFree(d_work));
      PetscCallMPI(MPIU_Allreduce(bv->work,y,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
    }
  } else {
    if (n) {
      PetscCall(PetscLogGpuTimeBegin());
      PetscCallHIPBLAS(hipblasXgemv(hipblashandle,HIPBLAS_OP_C,n,k,&sone,d_A,lda,d_x,one,&szero,d_work,one));
      PetscCall(PetscLogGpuTimeEnd());
    }
    if (!y) PetscCall(VecHIPRestoreArrayWrite(bv->buffer,&d_work));
    else {
      PetscCallHIP(hipMemcpy(y,d_work,k*sizeof(PetscScalar),hipMemcpyDeviceToHost));
      PetscCall(PetscLogGpuToCpu(k*sizeof(PetscScalar)));
      PetscCallHIP(hipFree(d_work));
    }
  }
  PetscCall(PetscLogGpuFlops(2.0*n*k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Scale n scalars
*/
PetscErrorCode BVScale_BLAS_HIP(BV,PetscInt n_,PetscScalar *d_A,PetscScalar alpha)
{
  PetscHipBLASInt n=0,one=1;
  hipblasHandle_t hipblashandle;

  PetscFunctionBegin;
  PetscCall(PetscHipBLASIntCast(n_,&n));
  if (PetscUnlikely(alpha == (PetscScalar)0.0)) PetscCallHIP(hipMemset(d_A,0,n*sizeof(PetscScalar)));
  else if (alpha != (PetscScalar)1.0) {
    PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallHIPBLAS(hipblasXscal(hipblashandle,n,&alpha,d_A,one));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(PetscLogGpuFlops(1.0*n));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Compute 2-norm of vector consisting of n scalars
*/
PetscErrorCode BVNorm_BLAS_HIP(BV,PetscInt n_,const PetscScalar *d_A,PetscReal *nrm)
{
  PetscHipBLASInt n=0,one=1;
  hipblasHandle_t hipblashandle;

  PetscFunctionBegin;
  PetscCall(PetscHipBLASIntCast(n_,&n));
  PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
  PetscCall(PetscLogGpuTimeBegin());
  PetscCallHIPBLAS(hipblasXnrm2(hipblashandle,n,d_A,one,nrm));
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops(2.0*n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Normalize the columns of A
*/
PetscErrorCode BVNormalize_BLAS_HIP(BV,PetscInt m_,PetscInt n_,PetscScalar *d_A,PetscInt lda_,PetscScalar *eigi)
{
  PetscInt        j,k;
  PetscReal       nrm,nrm1;
  PetscScalar     alpha;
  PetscHipBLASInt m=0,one=1;
  hipblasHandle_t hipblashandle;

  PetscFunctionBegin;
  PetscCall(PetscHipBLASIntCast(m_,&m));
  PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
  PetscCall(PetscLogGpuTimeBegin());
  for (j=0;j<n_;j++) {
    k = 1;
#if !defined(PETSC_USE_COMPLEX)
    if (eigi && eigi[j] != 0.0) k = 2;
#endif
    PetscCallHIPBLAS(hipblasXnrm2(hipblashandle,m,d_A+j*lda_,one,&nrm));
    if (k==2) {
      PetscCallHIPBLAS(hipblasXnrm2(hipblashandle,m,d_A+(j+1)*lda_,one,&nrm1));
      nrm = SlepcAbs(nrm,nrm1);
    }
    alpha = 1.0/nrm;
    PetscCallHIPBLAS(hipblasXscal(hipblashandle,m,&alpha,d_A+j*lda_,one));
    if (k==2) {
      PetscCallHIPBLAS(hipblasXscal(hipblashandle,m,&alpha,d_A+(j+1)*lda_,one));
      j++;
    }
  }
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogGpuFlops(3.0*m*n_));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_CleanCoefficients_HIP - Sets to zero all entries of column j of the bv buffer
*/
PetscErrorCode BV_CleanCoefficients_HIP(BV bv,PetscInt j,PetscScalar *h)
{
  PetscScalar    *d_hh,*d_a;
  PetscInt       i;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecHIPGetArray(bv->buffer,&d_a));
    PetscCall(PetscLogGpuTimeBegin());
    d_hh = d_a + j*(bv->nc+bv->m);
    PetscCallHIP(hipMemset(d_hh,0,(bv->nc+j)*sizeof(PetscScalar)));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(VecHIPRestoreArray(bv->buffer,&d_a));
  } else { /* cpu memory */
    for (i=0;i<bv->nc+j;i++) h[i] = 0.0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_AddCoefficients_HIP - Add the contents of the scratch (0-th column) of the bv buffer
   into column j of the bv buffer
 */
PetscErrorCode BV_AddCoefficients_HIP(BV bv,PetscInt j,PetscScalar *h,PetscScalar *c)
{
  PetscScalar     *d_h,*d_c,sone=1.0;
  PetscInt        i;
  PetscHipBLASInt idx=0,one=1;
  hipblasHandle_t hipblashandle;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
    PetscCall(VecHIPGetArray(bv->buffer,&d_c));
    d_h = d_c + j*(bv->nc+bv->m);
    PetscCall(PetscHipBLASIntCast(bv->nc+j,&idx));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallHIPBLAS(hipblasXaxpy(hipblashandle,idx,&sone,d_c,one,d_h,one));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(PetscLogGpuFlops(1.0*(bv->nc+j)));
    PetscCall(VecHIPRestoreArray(bv->buffer,&d_c));
  } else { /* cpu memory */
    for (i=0;i<bv->nc+j;i++) h[i] += c[i];
    PetscCall(PetscLogFlops(1.0*(bv->nc+j)));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SetValue_HIP - Sets value in row j (counted after the constraints) of column k
   of the coefficients array
*/
PetscErrorCode BV_SetValue_HIP(BV bv,PetscInt j,PetscInt k,PetscScalar *h,PetscScalar value)
{
  PetscScalar    *d_h,*a;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecHIPGetArray(bv->buffer,&a));
    PetscCall(PetscLogGpuTimeBegin());
    d_h = a + k*(bv->nc+bv->m) + bv->nc+j;
    PetscCallHIP(hipMemcpy(d_h,&value,sizeof(PetscScalar),hipMemcpyHostToDevice));
    PetscCall(PetscLogCpuToGpu(sizeof(PetscScalar)));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(VecHIPRestoreArray(bv->buffer,&a));
  } else { /* cpu memory */
    h[bv->nc+j] = value;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SquareSum_HIP - Returns the value h'*h, where h represents the contents of the
   coefficients array (up to position j)
*/
PetscErrorCode BV_SquareSum_HIP(BV bv,PetscInt j,PetscScalar *h,PetscReal *sum)
{
  const PetscScalar *d_h;
  PetscScalar       dot;
  PetscInt          i;
  PetscHipBLASInt   idx=0,one=1;
  hipblasHandle_t   hipblashandle;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(PetscHIPBLASGetHandle(&hipblashandle));
    PetscCall(VecHIPGetArrayRead(bv->buffer,&d_h));
    PetscCall(PetscHipBLASIntCast(bv->nc+j,&idx));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallHIPBLAS(hipblasXdot(hipblashandle,idx,d_h,one,d_h,one,&dot));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(PetscLogGpuFlops(2.0*(bv->nc+j)));
    *sum = PetscRealPart(dot);
    PetscCall(VecHIPRestoreArrayRead(bv->buffer,&d_h));
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
   BV_ApplySignature_HIP - Computes the pointwise product h*omega, where h represents
   the contents of the coefficients array (up to position j) and omega is the signature;
   if inverse=TRUE then the operation is h/omega
*/
PetscErrorCode BV_ApplySignature_HIP(BV bv,PetscInt j,PetscScalar *h,PetscBool inverse)
{
  PetscScalar       *d_h;
  const PetscScalar *d_omega,*omega;
  PetscInt          i,xcount;
  dim3              blocks3d, threads3d;

  PetscFunctionBegin;
  if (!(bv->nc+j)) PetscFunctionReturn(PETSC_SUCCESS);
  if (!h) {
    PetscCall(VecHIPGetArray(bv->buffer,&d_h));
    PetscCall(VecHIPGetArrayRead(bv->omega,&d_omega));
    PetscCall(SlepcKernelSetGrid1D(bv->nc+j,&blocks3d,&threads3d,&xcount));
    PetscCall(PetscLogGpuTimeBegin());
    if (inverse) {
      for (i=0;i<xcount;i++) PointwiseDiv_kernel<<<blocks3d,threads3d,0,0>>>(i,d_h,d_omega,bv->nc+j);
    } else {
      for (i=0;i<xcount;i++) PointwiseMult_kernel<<<blocks3d,threads3d,0,0>>>(i,d_h,d_omega,bv->nc+j);
    }
    PetscCallHIP(hipGetLastError());
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(PetscLogGpuFlops(1.0*(bv->nc+j)));
    PetscCall(VecHIPRestoreArrayRead(bv->omega,&d_omega));
    PetscCall(VecHIPRestoreArray(bv->buffer,&d_h));
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
   BV_SquareRoot_HIP - Returns the square root of position j (counted after the constraints)
   of the coefficients array
*/
PetscErrorCode BV_SquareRoot_HIP(BV bv,PetscInt j,PetscScalar *h,PetscReal *beta)
{
  const PetscScalar *d_h;
  PetscScalar       hh;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecHIPGetArrayRead(bv->buffer,&d_h));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallHIP(hipMemcpy(&hh,d_h+bv->nc+j,sizeof(PetscScalar),hipMemcpyDeviceToHost));
    PetscCall(PetscLogGpuToCpu(sizeof(PetscScalar)));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(BV_SafeSqrt(bv,hh,beta));
    PetscCall(VecHIPRestoreArrayRead(bv->buffer,&d_h));
  } else PetscCall(BV_SafeSqrt(bv,h[bv->nc+j],beta));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_StoreCoefficients_HIP - Copy the contents of the coefficients array to an array dest
   provided by the caller (only values from l to j are copied)
*/
PetscErrorCode BV_StoreCoefficients_HIP(BV bv,PetscInt j,PetscScalar *h,PetscScalar *dest)
{
  const PetscScalar *d_h,*d_a;
  PetscInt          i;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecHIPGetArrayRead(bv->buffer,&d_a));
    PetscCall(PetscLogGpuTimeBegin());
    d_h = d_a + j*(bv->nc+bv->m)+bv->nc;
    PetscCallHIP(hipMemcpy(dest-bv->l,d_h,(j-bv->l)*sizeof(PetscScalar),hipMemcpyDeviceToHost));
    PetscCall(PetscLogGpuToCpu((j-bv->l)*sizeof(PetscScalar)));
    PetscCall(PetscLogGpuTimeEnd());
    PetscCall(VecHIPRestoreArrayRead(bv->buffer,&d_a));
  } else {
    for (i=bv->l;i<j;i++) dest[i-bv->l] = h[bv->nc+i];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
