/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "cyclic" (HIP implementation)
*/
#include <slepc/private/svdimpl.h>
#include <petscdevice_hip.h>
#include "../src/svd/impls/cyclic/cyclic.h"

/* block alignment can differ between ranks, so copies must not log collective VecCopy() events */
static PetscErrorCode Copy_Cyclic_HIP(PetscScalar *dest,const PetscScalar *src,PetscInt n)
{
  hipStream_t stream;

  PetscFunctionBegin;
  if (n) {
    PetscCall(PetscGetCurrentHIPStream(&stream));
    PetscCall(PetscLogGpuTimeBegin());
    PetscCallHIP(hipMemcpyAsync(dest,src,n*sizeof(*src),hipMemcpyDeviceToDevice,stream));
    PetscCall(PetscLogGpuTimeEnd());
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMult_Cyclic_HIP(Mat B,Vec x,Vec y)
{
  SVD_CYCLIC_SHELL  *ctx;
  const PetscScalar *d_px,*read;
  PetscScalar       *d_py,*write;
  PetscInt          m,n;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatGetLocalSize(ctx->A,&m,&n));
  PetscCall(VecHIPGetArrayRead(x,&d_px));
  PetscCall(VecHIPGetArrayWrite(y,&d_py));
  /* dense matrix columns may misalign either block; retain the same Vec objects on every rank */
  if (((PETSC_UINTPTR_T)(d_px))%16) {
    PetscCall(VecHIPGetArrayWrite(ctx->x1,&write));
    PetscCall(Copy_Cyclic_HIP(write,d_px,m));
    PetscCall(VecHIPRestoreArrayWrite(ctx->x1,&write));
  } else PetscCall(VecHIPPlaceArray(ctx->x1,d_px));
  if (((PETSC_UINTPTR_T)(d_px+m))%16) {
    PetscCall(VecHIPGetArrayWrite(ctx->x2,&write));
    PetscCall(Copy_Cyclic_HIP(write,d_px+m,n));
    PetscCall(VecHIPRestoreArrayWrite(ctx->x2,&write));
  } else PetscCall(VecHIPPlaceArray(ctx->x2,d_px+m));
  if (!(((PETSC_UINTPTR_T)(d_py))%16)) PetscCall(VecHIPPlaceArray(ctx->y1,d_py));
  if (!(((PETSC_UINTPTR_T)(d_py+m))%16)) PetscCall(VecHIPPlaceArray(ctx->y2,d_py+m));
  PetscCall(MatMult(ctx->A,ctx->x2,ctx->y1));
  PetscCall(MatMult(ctx->AT,ctx->x1,ctx->y2));
  if (((PETSC_UINTPTR_T)(d_py))%16) {
    PetscCall(VecHIPGetArrayRead(ctx->y1,&read));
    PetscCall(Copy_Cyclic_HIP(d_py,read,m));
    PetscCall(VecHIPRestoreArrayRead(ctx->y1,&read));
  } else PetscCall(VecHIPResetArray(ctx->y1));
  if (((PETSC_UINTPTR_T)(d_py+m))%16) {
    PetscCall(VecHIPGetArrayRead(ctx->y2,&read));
    PetscCall(Copy_Cyclic_HIP(d_py+m,read,n));
    PetscCall(VecHIPRestoreArrayRead(ctx->y2,&read));
  } else PetscCall(VecHIPResetArray(ctx->y2));
  if (!(((PETSC_UINTPTR_T)(d_px))%16)) PetscCall(VecHIPResetArray(ctx->x1));
  if (!(((PETSC_UINTPTR_T)(d_px+m))%16)) PetscCall(VecHIPResetArray(ctx->x2));
  PetscCall(VecHIPRestoreArrayRead(x,&d_px));
  PetscCall(VecHIPRestoreArrayWrite(y,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMult_ECross_HIP(Mat B,Vec x,Vec y)
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
  PetscCall(VecHIPGetArrayRead(x,&d_px));
  PetscCall(VecHIPGetArrayWrite(y,&d_py));
  PetscCall(VecHIPPlaceArray(ctx->x1,d_px));
  PetscCall(VecHIPPlaceArray(ctx->y1,d_py));
  PetscCall(VecCopy(ctx->x1,ctx->y1));
  if (((PETSC_UINTPTR_T)(d_px+m))%16) {
    PetscCall(VecHIPGetArrayWrite(ctx->x2,&write));
    PetscCall(Copy_Cyclic_HIP(write,d_px+m,n));
    PetscCall(VecHIPRestoreArrayWrite(ctx->x2,&write));
  } else PetscCall(VecHIPPlaceArray(ctx->x2,d_px+m));
  if (!(((PETSC_UINTPTR_T)(d_py+m))%16)) PetscCall(VecHIPPlaceArray(ctx->y2,d_py+m));
  PetscCall(MatMult(ctx->A,ctx->x2,ctx->w));
  PetscCall(MatMult(ctx->AT,ctx->w,ctx->y2));
  if (((PETSC_UINTPTR_T)(d_py+m))%16) {
    PetscCall(VecHIPGetArrayRead(ctx->y2,&read));
    PetscCall(Copy_Cyclic_HIP(d_py+m,read,n));
    PetscCall(VecHIPRestoreArrayRead(ctx->y2,&read));
  } else PetscCall(VecHIPResetArray(ctx->y2));
  if (!(((PETSC_UINTPTR_T)(d_px+m))%16)) PetscCall(VecHIPResetArray(ctx->x2));
  PetscCall(VecHIPResetArray(ctx->x1));
  PetscCall(VecHIPResetArray(ctx->y1));
  PetscCall(VecHIPRestoreArrayRead(x,&d_px));
  PetscCall(VecHIPRestoreArrayWrite(y,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}
