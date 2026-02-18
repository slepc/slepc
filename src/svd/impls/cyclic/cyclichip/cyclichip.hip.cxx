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
#include "../src/svd/impls/cyclic/cyclic.h"

PetscErrorCode MatMult_Cyclic_HIP(Mat B,Vec x,Vec y)
{
  SVD_CYCLIC_SHELL  *ctx;
  const PetscScalar *d_px;
  PetscScalar       *d_py;
  PetscInt          m;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatGetLocalSize(ctx->A,&m,NULL));
  PetscCall(VecHIPGetArrayRead(x,&d_px));
  PetscCall(VecHIPGetArrayWrite(y,&d_py));
  PetscCall(VecHIPPlaceArray(ctx->x1,d_px));
  PetscCall(VecHIPPlaceArray(ctx->x2,d_px+m));
  PetscCall(VecHIPPlaceArray(ctx->y1,d_py));
  PetscCall(VecHIPPlaceArray(ctx->y2,d_py+m));
  if (!ctx->misaligned) {
    PetscCall(MatMult(ctx->A,ctx->x2,ctx->y1));
    PetscCall(MatMult(ctx->AT,ctx->x1,ctx->y2));
  } else { /* prevent HIP errors when bottom part is misaligned */
    PetscCall(VecCopy(ctx->x2,ctx->wx2));
    PetscCall(MatMult(ctx->A,ctx->wx2,ctx->y1));
    PetscCall(MatMult(ctx->AT,ctx->x1,ctx->wy2));
    PetscCall(VecCopy(ctx->wy2,ctx->y2));
  }
  PetscCall(VecHIPResetArray(ctx->x1));
  PetscCall(VecHIPResetArray(ctx->x2));
  PetscCall(VecHIPResetArray(ctx->y1));
  PetscCall(VecHIPResetArray(ctx->y2));
  PetscCall(VecHIPRestoreArrayRead(x,&d_px));
  PetscCall(VecHIPRestoreArrayWrite(y,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMult_ECross_HIP(Mat B,Vec x,Vec y)
{
  SVD_CYCLIC_SHELL  *ctx;
  const PetscScalar *d_px;
  PetscScalar       *d_py;
  PetscInt          mn,m,n;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatGetLocalSize(ctx->A,NULL,&n));
  PetscCall(VecGetLocalSize(y,&mn));
  m = mn-n;
  PetscCall(VecHIPGetArrayRead(x,&d_px));
  PetscCall(VecHIPGetArrayWrite(y,&d_py));
  PetscCall(VecHIPPlaceArray(ctx->x1,d_px));
  PetscCall(VecHIPPlaceArray(ctx->x2,d_px+m));
  PetscCall(VecHIPPlaceArray(ctx->y1,d_py));
  PetscCall(VecHIPPlaceArray(ctx->y2,d_py+m));
  PetscCall(VecCopy(ctx->x1,ctx->y1));
  if (!ctx->misaligned) {
    PetscCall(MatMult(ctx->A,ctx->x2,ctx->w));
    PetscCall(MatMult(ctx->AT,ctx->w,ctx->y2));
  } else { /* prevent HIP errors when bottom part is misaligned */
    PetscCall(VecCopy(ctx->x2,ctx->wx2));
    PetscCall(MatMult(ctx->A,ctx->wx2,ctx->w));
    PetscCall(MatMult(ctx->AT,ctx->w,ctx->wy2));
    PetscCall(VecCopy(ctx->wy2,ctx->y2));
  }
  PetscCall(VecHIPResetArray(ctx->x1));
  PetscCall(VecHIPResetArray(ctx->x2));
  PetscCall(VecHIPResetArray(ctx->y1));
  PetscCall(VecHIPResetArray(ctx->y2));
  PetscCall(VecHIPRestoreArrayRead(x,&d_px));
  PetscCall(VecHIPRestoreArrayWrite(y,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}
