/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "cyclic" (CUDA implementation)
*/
#include <slepc/private/svdimpl.h>
#include <petscdevice_cuda.h>
#include "../src/svd/impls/cyclic/cyclic.h"

/* block alignment can differ between ranks, so copies must not log collective VecCopy() events */
static PetscErrorCode Copy_Cyclic_CUDA(PetscScalar *dest,const PetscScalar *src,PetscInt n)
{
  cudaStream_t stream;

  PetscFunctionBegin;
  if (n) {
    PetscCall(PetscGetCurrentCUDAStream(&stream));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallCUDA(cudaMemcpyAsync(dest,src,n*sizeof(*src),cudaMemcpyDeviceToDevice,stream));
    PetscCall(PetscLogGpuTimeEnd());
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMult_Cyclic_CUDA(Mat B,Vec x,Vec y)
{
  SVD_CYCLIC_SHELL  *ctx;
  const PetscScalar *d_px,*read;
  PetscScalar       *d_py,*write;
  PetscInt          m,n;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatGetLocalSize(ctx->A,&m,&n));
  PetscCall(VecCUDAGetArrayRead(x,&d_px));
  PetscCall(VecCUDAGetArrayWrite(y,&d_py));
  /* dense matrix columns may misalign either block; retain the same Vec objects on every rank */
  if (((PETSC_UINTPTR_T)(d_px))%16) {
    PetscCall(VecCUDAGetArrayWrite(ctx->x1,&write));
    PetscCall(Copy_Cyclic_CUDA(write,d_px,m));
    PetscCall(VecCUDARestoreArrayWrite(ctx->x1,&write));
  } else PetscCall(VecCUDAPlaceArray(ctx->x1,d_px));
  if (((PETSC_UINTPTR_T)(d_px+m))%16) {
    PetscCall(VecCUDAGetArrayWrite(ctx->x2,&write));
    PetscCall(Copy_Cyclic_CUDA(write,d_px+m,n));
    PetscCall(VecCUDARestoreArrayWrite(ctx->x2,&write));
  } else PetscCall(VecCUDAPlaceArray(ctx->x2,d_px+m));
  if (!(((PETSC_UINTPTR_T)(d_py))%16)) PetscCall(VecCUDAPlaceArray(ctx->y1,d_py));
  if (!(((PETSC_UINTPTR_T)(d_py+m))%16)) PetscCall(VecCUDAPlaceArray(ctx->y2,d_py+m));
  PetscCall(MatMult(ctx->A,ctx->x2,ctx->y1));
  PetscCall(MatMult(ctx->AT,ctx->x1,ctx->y2));
  if (((PETSC_UINTPTR_T)(d_py))%16) {
    PetscCall(VecCUDAGetArrayRead(ctx->y1,&read));
    PetscCall(Copy_Cyclic_CUDA(d_py,read,m));
    PetscCall(VecCUDARestoreArrayRead(ctx->y1,&read));
  } else PetscCall(VecCUDAResetArray(ctx->y1));
  if (((PETSC_UINTPTR_T)(d_py+m))%16) {
    PetscCall(VecCUDAGetArrayRead(ctx->y2,&read));
    PetscCall(Copy_Cyclic_CUDA(d_py+m,read,n));
    PetscCall(VecCUDARestoreArrayRead(ctx->y2,&read));
  } else PetscCall(VecCUDAResetArray(ctx->y2));
  if (!(((PETSC_UINTPTR_T)(d_px))%16)) PetscCall(VecCUDAResetArray(ctx->x1));
  if (!(((PETSC_UINTPTR_T)(d_px+m))%16)) PetscCall(VecCUDAResetArray(ctx->x2));
  PetscCall(VecCUDARestoreArrayRead(x,&d_px));
  PetscCall(VecCUDARestoreArrayWrite(y,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMult_ECross_CUDA(Mat B,Vec x,Vec y)
{
  SVD_CYCLIC_SHELL  *ctx;
  const PetscScalar *d_px,*read;
  PetscScalar       *d_py,*write;
  PetscInt          mn,m,n;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatGetLocalSize(ctx->A,NULL,&n));
  PetscCall(VecGetLocalSize(y,&mn));
  m = mn-n;
  PetscCall(VecCUDAGetArrayRead(x,&d_px));
  PetscCall(VecCUDAGetArrayWrite(y,&d_py));
  PetscCall(VecCUDAPlaceArray(ctx->x1,d_px));
  PetscCall(VecCUDAPlaceArray(ctx->y1,d_py));
  PetscCall(VecCopy(ctx->x1,ctx->y1));
  if (((PETSC_UINTPTR_T)(d_px+m))%16) {
    PetscCall(VecCUDAGetArrayWrite(ctx->x2,&write));
    PetscCall(Copy_Cyclic_CUDA(write,d_px+m,n));
    PetscCall(VecCUDARestoreArrayWrite(ctx->x2,&write));
  } else PetscCall(VecCUDAPlaceArray(ctx->x2,d_px+m));
  if (!(((PETSC_UINTPTR_T)(d_py+m))%16)) PetscCall(VecCUDAPlaceArray(ctx->y2,d_py+m));
  PetscCall(MatMult(ctx->A,ctx->x2,ctx->w));
  PetscCall(MatMult(ctx->AT,ctx->w,ctx->y2));
  if (((PETSC_UINTPTR_T)(d_py+m))%16) {
    PetscCall(VecCUDAGetArrayRead(ctx->y2,&read));
    PetscCall(Copy_Cyclic_CUDA(d_py+m,read,n));
    PetscCall(VecCUDARestoreArrayRead(ctx->y2,&read));
  } else PetscCall(VecCUDAResetArray(ctx->y2));
  if (!(((PETSC_UINTPTR_T)(d_px+m))%16)) PetscCall(VecCUDAResetArray(ctx->x2));
  PetscCall(VecCUDAResetArray(ctx->x1));
  PetscCall(VecCUDAResetArray(ctx->y1));
  PetscCall(VecCUDARestoreArrayRead(x,&d_px));
  PetscCall(VecCUDARestoreArrayWrite(y,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}
