/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "cyclic"

   Method: Uses a Hermitian eigensolver for H(A) = [ 0  A ; A^T 0 ]
*/

#include <slepc/private/svdimpl.h>                /*I "slepcsvd.h" I*/
#include <slepc/private/bvimpl.h>
#include "cyclic.h"

static PetscErrorCode MatMult_Cyclic(Mat B,Vec x,Vec y)
{
  SVD_CYCLIC_SHELL  *ctx;
  const PetscScalar *px;
  PetscScalar       *py;
  PetscInt          m;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatGetLocalSize(ctx->A,&m,NULL));
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArrayWrite(y,&py));
  PetscCall(VecPlaceArray(ctx->x1,px));
  PetscCall(VecPlaceArray(ctx->x2,px+m));
  PetscCall(VecPlaceArray(ctx->y1,py));
  PetscCall(VecPlaceArray(ctx->y2,py+m));
  PetscCall(MatMult(ctx->A,ctx->x2,ctx->y1));
  PetscCall(MatMult(ctx->AT,ctx->x1,ctx->y2));
  PetscCall(VecResetArray(ctx->x1));
  PetscCall(VecResetArray(ctx->x2));
  PetscCall(VecResetArray(ctx->y1));
  PetscCall(VecResetArray(ctx->y2));
  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArrayWrite(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatGetDiagonal_Cyclic(Mat B,Vec diag)
{
  PetscFunctionBegin;
  PetscCall(VecSet(diag,0.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_Cyclic(Mat B)
{
  SVD_CYCLIC_SHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(VecDestroy(&ctx->x1));
  PetscCall(VecDestroy(&ctx->x2));
  PetscCall(VecDestroy(&ctx->y1));
  PetscCall(VecDestroy(&ctx->y2));
  if (ctx->misaligned) {
    PetscCall(VecDestroy(&ctx->wx2));
    PetscCall(VecDestroy(&ctx->wy2));
  }
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Builds cyclic matrix   C = | 0   A |
                              | AT  0 |
*/
static PetscErrorCode SVDCyclicGetCyclicMat(SVD svd,Mat A,Mat AT,Mat *C)
{
  SVD_CYCLIC       *cyclic = (SVD_CYCLIC*)svd->data;
  SVD_CYCLIC_SHELL *ctx;
  PetscInt         i,M,N,m,n,Istart,Iend;
  VecType          vtype;
  Mat              Zm,Zn;
#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
  PetscBool        gpu;
  const PetscInt   *ranges;
  PetscMPIInt      size;
#endif

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));

  if (cyclic->explicitmatrix) {
    PetscCheck(svd->expltrans,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Cannot use explicit cyclic matrix with implicit transpose");
    PetscCall(MatCreate(PetscObjectComm((PetscObject)svd),&Zm));
    PetscCall(MatSetSizes(Zm,m,m,M,M));
    PetscCall(MatSetFromOptions(Zm));
    PetscCall(MatGetOwnershipRange(Zm,&Istart,&Iend));
    for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(Zm,i,i,0.0,INSERT_VALUES));
    PetscCall(MatAssemblyBegin(Zm,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zm,MAT_FINAL_ASSEMBLY));
    PetscCall(MatCreate(PetscObjectComm((PetscObject)svd),&Zn));
    PetscCall(MatSetSizes(Zn,n,n,N,N));
    PetscCall(MatSetFromOptions(Zn));
    PetscCall(MatGetOwnershipRange(Zn,&Istart,&Iend));
    for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(Zn,i,i,0.0,INSERT_VALUES));
    PetscCall(MatAssemblyBegin(Zn,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zn,MAT_FINAL_ASSEMBLY));
    PetscCall(MatCreateTile(1.0,Zm,1.0,A,1.0,AT,1.0,Zn,C));
    PetscCall(MatDestroy(&Zm));
    PetscCall(MatDestroy(&Zn));
  } else {
    PetscCall(PetscNew(&ctx));
    ctx->A       = A;
    ctx->AT      = AT;
    ctx->swapped = svd->swapped;
    PetscCall(MatCreateVecsEmpty(A,&ctx->x2,&ctx->x1));
    PetscCall(MatCreateVecsEmpty(A,&ctx->y2,&ctx->y1));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)svd),m+n,m+n,M+N,M+N,ctx,C));
    PetscCall(MatShellSetOperation(*C,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Cyclic));
    PetscCall(MatShellSetOperation(*C,MATOP_DESTROY,(void(*)(void))MatDestroy_Cyclic));
#if defined(PETSC_HAVE_CUDA)
    PetscCall(PetscObjectTypeCompareAny((PetscObject)(svd->swapped?AT:A),&gpu,MATSEQAIJCUSPARSE,MATMPIAIJCUSPARSE,""));
    if (gpu) PetscCall(MatShellSetOperation(*C,MATOP_MULT,(void(*)(void))MatMult_Cyclic_CUDA));
    else
#elif defined(PETSC_HAVE_HIP)
    PetscCall(PetscObjectTypeCompareAny((PetscObject)(svd->swapped?AT:A),&gpu,MATSEQAIJHIPSPARSE,MATMPIAIJHIPSPARSE,""));
    if (gpu) PetscCall(MatShellSetOperation(*C,MATOP_MULT,(void(*)(void))MatMult_Cyclic_HIP));
    else
#endif
      PetscCall(MatShellSetOperation(*C,MATOP_MULT,(void(*)(void))MatMult_Cyclic));
    PetscCall(MatGetVecType(A,&vtype));
    PetscCall(MatSetVecType(*C,vtype));
#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
    if (gpu) {
      /* check alignment of bottom block */
      PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ctx->x1),&size));
      PetscCall(VecGetOwnershipRanges(ctx->x1,&ranges));
      for (i=0;i<size;i++) {
        ctx->misaligned = (((ranges[i+1]-ranges[i])*sizeof(PetscScalar))%16)? PETSC_TRUE: PETSC_FALSE;
        if (ctx->misaligned) break;
      }
      if (ctx->misaligned) {  /* create work vectors for MatMult */
        PetscCall(VecDuplicate(ctx->x2,&ctx->wx2));
        PetscCall(VecDuplicate(ctx->y2,&ctx->wy2));
      }
    }
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_ECross(Mat B,Vec x,Vec y)
{
  SVD_CYCLIC_SHELL  *ctx;
  const PetscScalar *px;
  PetscScalar       *py;
  PetscInt          mn,m,n;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatGetLocalSize(ctx->A,NULL,&n));
  PetscCall(VecGetLocalSize(y,&mn));
  m = mn-n;
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArrayWrite(y,&py));
  PetscCall(VecPlaceArray(ctx->x1,px));
  PetscCall(VecPlaceArray(ctx->x2,px+m));
  PetscCall(VecPlaceArray(ctx->y1,py));
  PetscCall(VecPlaceArray(ctx->y2,py+m));
  PetscCall(VecCopy(ctx->x1,ctx->y1));
  PetscCall(MatMult(ctx->A,ctx->x2,ctx->w));
  PetscCall(MatMult(ctx->AT,ctx->w,ctx->y2));
  PetscCall(VecResetArray(ctx->x1));
  PetscCall(VecResetArray(ctx->x2));
  PetscCall(VecResetArray(ctx->y1));
  PetscCall(VecResetArray(ctx->y2));
  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArrayWrite(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatGetDiagonal_ECross(Mat B,Vec d)
{
  SVD_CYCLIC_SHELL  *ctx;
  PetscScalar       *pd;
  PetscMPIInt       len;
  PetscInt          mn,m,n,N,i,j,start,end,ncols;
  PetscScalar       *work1,*work2,*diag;
  const PetscInt    *cols;
  const PetscScalar *vals;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatGetLocalSize(ctx->A,NULL,&n));
  PetscCall(VecGetLocalSize(d,&mn));
  m = mn-n;
  PetscCall(VecGetArrayWrite(d,&pd));
  PetscCall(VecPlaceArray(ctx->y1,pd));
  PetscCall(VecSet(ctx->y1,1.0));
  PetscCall(VecResetArray(ctx->y1));
  PetscCall(VecPlaceArray(ctx->y2,pd+m));
  if (!ctx->diag) {
    /* compute diagonal from rows and store in ctx->diag */
    PetscCall(VecDuplicate(ctx->y2,&ctx->diag));
    PetscCall(MatGetSize(ctx->A,NULL,&N));
    PetscCall(PetscCalloc2(N,&work1,N,&work2));
    if (ctx->swapped) {
      PetscCall(MatGetOwnershipRange(ctx->AT,&start,&end));
      for (i=start;i<end;i++) {
        PetscCall(MatGetRow(ctx->AT,i,&ncols,NULL,&vals));
        for (j=0;j<ncols;j++) work1[i] += vals[j]*vals[j];
        PetscCall(MatRestoreRow(ctx->AT,i,&ncols,NULL,&vals));
      }
    } else {
      PetscCall(MatGetOwnershipRange(ctx->A,&start,&end));
      for (i=start;i<end;i++) {
        PetscCall(MatGetRow(ctx->A,i,&ncols,&cols,&vals));
        for (j=0;j<ncols;j++) work1[cols[j]] += vals[j]*vals[j];
        PetscCall(MatRestoreRow(ctx->A,i,&ncols,&cols,&vals));
      }
    }
    PetscCall(PetscMPIIntCast(N,&len));
    PetscCallMPI(MPIU_Allreduce(work1,work2,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)B)));
    PetscCall(VecGetOwnershipRange(ctx->diag,&start,&end));
    PetscCall(VecGetArrayWrite(ctx->diag,&diag));
    for (i=start;i<end;i++) diag[i-start] = work2[i];
    PetscCall(VecRestoreArrayWrite(ctx->diag,&diag));
    PetscCall(PetscFree2(work1,work2));
  }
  PetscCall(VecCopy(ctx->diag,ctx->y2));
  PetscCall(VecResetArray(ctx->y2));
  PetscCall(VecRestoreArrayWrite(d,&pd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_ECross(Mat B)
{
  SVD_CYCLIC_SHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(VecDestroy(&ctx->x1));
  PetscCall(VecDestroy(&ctx->x2));
  PetscCall(VecDestroy(&ctx->y1));
  PetscCall(VecDestroy(&ctx->y2));
  PetscCall(VecDestroy(&ctx->diag));
  PetscCall(VecDestroy(&ctx->w));
  if (ctx->misaligned) {
    PetscCall(VecDestroy(&ctx->wx2));
    PetscCall(VecDestroy(&ctx->wy2));
  }
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Builds extended cross product matrix   C = | I_m   0  |
                                              |  0  AT*A |
   t is an auxiliary Vec used to take the dimensions of the upper block
*/
static PetscErrorCode SVDCyclicGetECrossMat(SVD svd,Mat A,Mat AT,Mat *C,Vec t)
{
  SVD_CYCLIC       *cyclic = (SVD_CYCLIC*)svd->data;
  SVD_CYCLIC_SHELL *ctx;
  PetscInt         i,M,N,m,n,Istart,Iend;
  VecType          vtype;
  Mat              Id,Zm,Zn,ATA;
#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
  PetscBool        gpu;
  const PetscInt   *ranges;
  PetscMPIInt      size;
#endif

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,NULL,&N));
  PetscCall(MatGetLocalSize(A,NULL,&n));
  PetscCall(VecGetSize(t,&M));
  PetscCall(VecGetLocalSize(t,&m));

  if (cyclic->explicitmatrix) {
    PetscCheck(svd->expltrans,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Cannot use explicit cyclic matrix with implicit transpose");
    PetscCall(MatCreateConstantDiagonal(PetscObjectComm((PetscObject)svd),m,m,M,M,1.0,&Id));
    PetscCall(MatCreate(PetscObjectComm((PetscObject)svd),&Zm));
    PetscCall(MatSetSizes(Zm,m,n,M,N));
    PetscCall(MatSetFromOptions(Zm));
    PetscCall(MatGetOwnershipRange(Zm,&Istart,&Iend));
    for (i=Istart;i<Iend;i++) {
      if (i<N) PetscCall(MatSetValue(Zm,i,i,0.0,INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(Zm,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zm,MAT_FINAL_ASSEMBLY));
    PetscCall(MatCreate(PetscObjectComm((PetscObject)svd),&Zn));
    PetscCall(MatSetSizes(Zn,n,m,N,M));
    PetscCall(MatSetFromOptions(Zn));
    PetscCall(MatGetOwnershipRange(Zn,&Istart,&Iend));
    for (i=Istart;i<Iend;i++) {
      if (i<m) PetscCall(MatSetValue(Zn,i,i,0.0,INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(Zn,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Zn,MAT_FINAL_ASSEMBLY));
    PetscCall(MatProductCreate(AT,A,NULL,&ATA));
    PetscCall(MatProductSetType(ATA,MATPRODUCT_AB));
    PetscCall(MatProductSetFromOptions(ATA));
    PetscCall(MatProductSymbolic(ATA));
    PetscCall(MatProductNumeric(ATA));
    PetscCall(MatCreateTile(1.0,Id,1.0,Zm,1.0,Zn,1.0,ATA,C));
    PetscCall(MatDestroy(&Id));
    PetscCall(MatDestroy(&Zm));
    PetscCall(MatDestroy(&Zn));
    PetscCall(MatDestroy(&ATA));
  } else {
    PetscCall(PetscNew(&ctx));
    ctx->A       = A;
    ctx->AT      = AT;
    ctx->swapped = svd->swapped;
    PetscCall(VecDuplicateEmpty(t,&ctx->x1));
    PetscCall(VecDuplicateEmpty(t,&ctx->y1));
    PetscCall(MatCreateVecsEmpty(A,&ctx->x2,NULL));
    PetscCall(MatCreateVecsEmpty(A,&ctx->y2,NULL));
    PetscCall(MatCreateVecs(A,NULL,&ctx->w));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)svd),m+n,m+n,M+N,M+N,ctx,C));
    PetscCall(MatShellSetOperation(*C,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_ECross));
    PetscCall(MatShellSetOperation(*C,MATOP_DESTROY,(void(*)(void))MatDestroy_ECross));
#if defined(PETSC_HAVE_CUDA)
    PetscCall(PetscObjectTypeCompareAny((PetscObject)(svd->swapped?AT:A),&gpu,MATSEQAIJCUSPARSE,MATMPIAIJCUSPARSE,""));
    if (gpu) PetscCall(MatShellSetOperation(*C,MATOP_MULT,(void(*)(void))MatMult_ECross_CUDA));
    else
#elif defined(PETSC_HAVE_HIP)
    PetscCall(PetscObjectTypeCompareAny((PetscObject)(svd->swapped?AT:A),&gpu,MATSEQAIJHIPSPARSE,MATMPIAIJHIPSPARSE,""));
    if (gpu) PetscCall(MatShellSetOperation(*C,MATOP_MULT,(void(*)(void))MatMult_ECross_HIP));
    else
#endif
      PetscCall(MatShellSetOperation(*C,MATOP_MULT,(void(*)(void))MatMult_ECross));
    PetscCall(MatGetVecType(A,&vtype));
    PetscCall(MatSetVecType(*C,vtype));
#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
    if (gpu) {
      /* check alignment of bottom block */
      PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ctx->x1),&size));
      PetscCall(VecGetOwnershipRanges(ctx->x1,&ranges));
      for (i=0;i<size;i++) {
        ctx->misaligned = (((ranges[i+1]-ranges[i])*sizeof(PetscScalar))%16)? PETSC_TRUE: PETSC_FALSE;
        if (ctx->misaligned) break;
      }
      if (ctx->misaligned) {  /* create work vectors for MatMult */
        PetscCall(VecDuplicate(ctx->x2,&ctx->wx2));
        PetscCall(VecDuplicate(ctx->y2,&ctx->wy2));
      }
    }
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Convergence test relative to the norm of R (used in GSVD only) */
static PetscErrorCode EPSConv_Cyclic(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  SVD svd = (SVD)ctx;

  PetscFunctionBegin;
  *errest = res/PetscMax(svd->nrma,svd->nrmb);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSetUp_Cyclic(SVD svd)
{
  SVD_CYCLIC        *cyclic = (SVD_CYCLIC*)svd->data;
  PetscInt          M,N,m,n,p,k,i,isl,offset,nev,ncv,mpd,maxit;
  PetscReal         tol;
  const PetscScalar *isa,*oa;
  PetscScalar       *va;
  EPSProblemType    ptype;
  PetscBool         trackall,issinv;
  Vec               v,t;
  ST                st;
  Mat               Omega;
  MatType           Atype;

  PetscFunctionBegin;
  if (svd->nsv==0) svd->nsv = 1;
  PetscCall(MatGetSize(svd->A,&M,&N));
  PetscCall(MatGetLocalSize(svd->A,&m,&n));
  if (!cyclic->eps) PetscCall(SVDCyclicGetEPS(svd,&cyclic->eps));
  PetscCall(MatDestroy(&cyclic->C));
  PetscCall(MatDestroy(&cyclic->D));
  if (svd->isgeneralized) {
    if (svd->which==SVD_SMALLEST) {  /* alternative pencil */
      PetscCall(MatCreateVecs(svd->B,NULL,&t));
      PetscCall(SVDCyclicGetCyclicMat(svd,svd->B,svd->BT,&cyclic->C));
      PetscCall(SVDCyclicGetECrossMat(svd,svd->A,svd->AT,&cyclic->D,t));
    } else {
      PetscCall(MatCreateVecs(svd->A,NULL,&t));
      PetscCall(SVDCyclicGetCyclicMat(svd,svd->A,svd->AT,&cyclic->C));
      PetscCall(SVDCyclicGetECrossMat(svd,svd->B,svd->BT,&cyclic->D,t));
    }
    PetscCall(VecDestroy(&t));
    PetscCall(EPSSetOperators(cyclic->eps,cyclic->C,cyclic->D));
    PetscCall(EPSGetProblemType(cyclic->eps,&ptype));
    if (!ptype) PetscCall(EPSSetProblemType(cyclic->eps,EPS_GHEP));
  } else if (svd->ishyperbolic) {
    PetscCall(SVDCyclicGetCyclicMat(svd,svd->A,svd->AT,&cyclic->C));
    PetscCall(MatCreateVecs(cyclic->C,&v,NULL));
    PetscCall(VecSet(v,1.0));
    PetscCall(VecGetArrayRead(svd->omega,&oa));
    PetscCall(VecGetArray(v,&va));
    if (svd->swapped) PetscCall(PetscArraycpy(va+m,oa,n));
    else PetscCall(PetscArraycpy(va,oa,m));
    PetscCall(VecRestoreArrayRead(svd->omega,&oa));
    PetscCall(VecRestoreArray(v,&va));
    PetscCall(MatGetType(svd->OP,&Atype));
    PetscCall(MatCreate(PetscObjectComm((PetscObject)svd),&Omega));
    PetscCall(MatSetSizes(Omega,m+n,m+n,M+N,M+N));
    PetscCall(MatSetType(Omega,Atype));
    PetscCall(MatDiagonalSet(Omega,v,INSERT_VALUES));
    PetscCall(EPSSetOperators(cyclic->eps,cyclic->C,Omega));
    PetscCall(EPSSetProblemType(cyclic->eps,EPS_GHIEP));
    PetscCall(MatDestroy(&Omega));
    PetscCall(VecDestroy(&v));
  } else {
    PetscCall(SVDCyclicGetCyclicMat(svd,svd->A,svd->AT,&cyclic->C));
    PetscCall(EPSSetOperators(cyclic->eps,cyclic->C,NULL));
    PetscCall(EPSSetProblemType(cyclic->eps,EPS_HEP));
  }
  if (!cyclic->usereps) {
    if (svd->which == SVD_LARGEST) {
      PetscCall(EPSGetST(cyclic->eps,&st));
      PetscCall(PetscObjectTypeCompare((PetscObject)st,STSINVERT,&issinv));
      if (issinv) PetscCall(EPSSetWhichEigenpairs(cyclic->eps,EPS_TARGET_MAGNITUDE));
      else if (svd->ishyperbolic) PetscCall(EPSSetWhichEigenpairs(cyclic->eps,EPS_LARGEST_MAGNITUDE));
      else PetscCall(EPSSetWhichEigenpairs(cyclic->eps,EPS_LARGEST_REAL));
    } else {
      if (svd->isgeneralized) {  /* computes sigma^{-1} via alternative pencil */
        PetscCall(EPSSetWhichEigenpairs(cyclic->eps,EPS_LARGEST_REAL));
      } else {
        if (svd->ishyperbolic) PetscCall(EPSSetWhichEigenpairs(cyclic->eps,EPS_TARGET_MAGNITUDE));
        else PetscCall(EPSSetEigenvalueComparison(cyclic->eps,SlepcCompareSmallestPosReal,NULL));
        PetscCall(EPSSetTarget(cyclic->eps,0.0));
      }
    }
    PetscCall(EPSGetDimensions(cyclic->eps,&nev,&ncv,&mpd));
    PetscCheck(nev==1 || nev>=2*svd->nsv,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_WRONG,"The number of requested eigenvalues %" PetscInt_FMT " must be at least 2*%" PetscInt_FMT,nev,svd->nsv);
    nev = PetscMax(nev,2*svd->nsv);
    if (ncv==PETSC_DETERMINE && svd->ncv!=PETSC_DETERMINE) ncv = PetscMax(3*svd->nsv,svd->ncv);
    if (mpd==PETSC_DETERMINE && svd->mpd!=PETSC_DETERMINE) mpd = svd->mpd;
    PetscCall(EPSSetDimensions(cyclic->eps,nev,ncv,mpd));
    PetscCall(EPSGetTolerances(cyclic->eps,&tol,&maxit));
    if (tol==(PetscReal)PETSC_DETERMINE) tol = svd->tol==(PetscReal)PETSC_DETERMINE? SLEPC_DEFAULT_TOL/10.0: svd->tol;
    if (maxit==PETSC_DETERMINE) maxit = svd->max_it;
    PetscCall(EPSSetTolerances(cyclic->eps,tol,maxit));
    switch (svd->conv) {
    case SVD_CONV_ABS:
      PetscCall(EPSSetConvergenceTest(cyclic->eps,EPS_CONV_ABS));break;
    case SVD_CONV_REL:
      PetscCall(EPSSetConvergenceTest(cyclic->eps,EPS_CONV_REL));break;
    case SVD_CONV_NORM:
      if (svd->isgeneralized) {
        if (!svd->nrma) PetscCall(MatNorm(svd->OP,NORM_INFINITY,&svd->nrma));
        if (!svd->nrmb) PetscCall(MatNorm(svd->OPb,NORM_INFINITY,&svd->nrmb));
        PetscCall(EPSSetConvergenceTestFunction(cyclic->eps,EPSConv_Cyclic,svd,NULL));
      } else {
        PetscCall(EPSSetConvergenceTest(cyclic->eps,EPS_CONV_NORM));break;
      }
      break;
    case SVD_CONV_MAXIT:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Maxit convergence test not supported in this solver");
    case SVD_CONV_USER:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"User-defined convergence test not supported in this solver");
    }
  }
  SVDCheckUnsupported(svd,SVD_FEATURE_STOPPING);
  /* Transfer the trackall option from svd to eps */
  PetscCall(SVDGetTrackAll(svd,&trackall));
  PetscCall(EPSSetTrackAll(cyclic->eps,trackall));
  /* Transfer the initial subspace from svd to eps */
  if (svd->nini<0 || svd->ninil<0) {
    for (i=0;i<-PetscMin(svd->nini,svd->ninil);i++) {
      PetscCall(MatCreateVecs(cyclic->C,&v,NULL));
      PetscCall(VecGetArrayWrite(v,&va));
      if (svd->isgeneralized) PetscCall(MatGetLocalSize(svd->B,&p,NULL));
      k = (svd->isgeneralized && svd->which==SVD_SMALLEST)? p: m;  /* size of upper block row */
      if (i<-svd->ninil) {
        PetscCall(VecGetArrayRead(svd->ISL[i],&isa));
        if (svd->isgeneralized) {
          PetscCall(VecGetLocalSize(svd->ISL[i],&isl));
          PetscCheck(isl==m+p,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Size mismatch for left initial vector");
          offset = (svd->which==SVD_SMALLEST)? m: 0;
          PetscCall(PetscArraycpy(va,isa+offset,k));
        } else {
          PetscCall(VecGetLocalSize(svd->ISL[i],&isl));
          PetscCheck(isl==k,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Size mismatch for left initial vector");
          PetscCall(PetscArraycpy(va,isa,k));
        }
        PetscCall(VecRestoreArrayRead(svd->IS[i],&isa));
      } else PetscCall(PetscArrayzero(&va,k));
      if (i<-svd->nini) {
        PetscCall(VecGetLocalSize(svd->IS[i],&isl));
        PetscCheck(isl==n,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Size mismatch for right initial vector");
        PetscCall(VecGetArrayRead(svd->IS[i],&isa));
        PetscCall(PetscArraycpy(va+k,isa,n));
        PetscCall(VecRestoreArrayRead(svd->IS[i],&isa));
      } else PetscCall(PetscArrayzero(va+k,n));
      PetscCall(VecRestoreArrayWrite(v,&va));
      PetscCall(VecDestroy(&svd->IS[i]));
      svd->IS[i] = v;
    }
    svd->nini = PetscMin(svd->nini,svd->ninil);
    PetscCall(EPSSetInitialSpace(cyclic->eps,-svd->nini,svd->IS));
    PetscCall(SlepcBasisDestroy_Private(&svd->nini,&svd->IS));
    PetscCall(SlepcBasisDestroy_Private(&svd->ninil,&svd->ISL));
  }
  PetscCall(EPSSetUp(cyclic->eps));
  PetscCall(EPSGetDimensions(cyclic->eps,NULL,&svd->ncv,&svd->mpd));
  svd->ncv = PetscMin(svd->ncv,PetscMin(M,N));
  PetscCall(EPSGetTolerances(cyclic->eps,NULL,&svd->max_it));
  if (svd->tol==(PetscReal)PETSC_DETERMINE) svd->tol = SLEPC_DEFAULT_TOL;

  svd->leftbasis = PETSC_TRUE;
  PetscCall(SVDAllocateSolution(svd,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCyclicCheckEigenvalue(SVD svd,PetscScalar er,PetscScalar ei,PetscReal *sigma,PetscBool *isreal)
{
  PetscFunctionBegin;
  if (svd->ishyperbolic && PetscDefined(USE_COMPLEX) && PetscAbsReal(PetscImaginaryPart(er))>10*PetscAbsReal(PetscRealPart(er))) {
    *sigma = PetscImaginaryPart(er);
    if (isreal) *isreal = PETSC_FALSE;
  } else if (svd->ishyperbolic && !PetscDefined(USE_COMPLEX) && PetscAbsScalar(ei)>10*PetscAbsScalar(er)) {
    *sigma = PetscRealPart(ei);
    if (isreal) *isreal = PETSC_FALSE;
  } else {
    *sigma = PetscRealPart(er);
    if (isreal) *isreal = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_Cyclic(SVD svd)
{
  SVD_CYCLIC     *cyclic = (SVD_CYCLIC*)svd->data;
  PetscInt       i,j,nconv;
  PetscScalar    er,ei;
  PetscReal      sigma;

  PetscFunctionBegin;
  PetscCall(EPSSolve(cyclic->eps));
  PetscCall(EPSGetConverged(cyclic->eps,&nconv));
  PetscCall(EPSGetIterationNumber(cyclic->eps,&svd->its));
  PetscCall(EPSGetConvergedReason(cyclic->eps,(EPSConvergedReason*)&svd->reason));
  for (i=0,j=0;i<nconv;i++) {
    PetscCall(EPSGetEigenvalue(cyclic->eps,i,&er,&ei));
    PetscCall(SVDCyclicCheckEigenvalue(svd,er,ei,&sigma,NULL));
    if (sigma>0.0) {
      if (svd->isgeneralized && svd->which==SVD_SMALLEST) svd->sigma[j] = 1.0/sigma;
      else svd->sigma[j] = sigma;
      j++;
    }
  }
  svd->nconv = j;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDComputeVectors_Cyclic_Standard(SVD svd)
{
  SVD_CYCLIC        *cyclic = (SVD_CYCLIC*)svd->data;
  PetscInt          i,j,m,nconv;
  PetscScalar       er,ei;
  PetscReal         sigma;
  const PetscScalar *px;
  Vec               x,x1,x2;

  PetscFunctionBegin;
  PetscCall(MatCreateVecs(cyclic->C,&x,NULL));
  PetscCall(MatGetLocalSize(svd->A,&m,NULL));
  PetscCall(MatCreateVecsEmpty(svd->A,&x2,&x1));
  PetscCall(EPSGetConverged(cyclic->eps,&nconv));
  for (i=0,j=0;i<nconv;i++) {
    PetscCall(EPSGetEigenpair(cyclic->eps,i,&er,&ei,x,NULL));
    PetscCall(SVDCyclicCheckEigenvalue(svd,er,ei,&sigma,NULL));
    if (sigma<0.0) continue;
    PetscCall(VecGetArrayRead(x,&px));
    PetscCall(VecPlaceArray(x1,px));
    PetscCall(VecPlaceArray(x2,px+m));
    PetscCall(BVInsertVec(svd->U,j,x1));
    PetscCall(BVScaleColumn(svd->U,j,PETSC_SQRT2));
    PetscCall(BVInsertVec(svd->V,j,x2));
    PetscCall(BVScaleColumn(svd->V,j,PETSC_SQRT2));
    PetscCall(VecResetArray(x1));
    PetscCall(VecResetArray(x2));
    PetscCall(VecRestoreArrayRead(x,&px));
    j++;
  }
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&x1));
  PetscCall(VecDestroy(&x2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDComputeVectors_Cyclic_Generalized(SVD svd)
{
  SVD_CYCLIC        *cyclic = (SVD_CYCLIC*)svd->data;
  PetscInt          i,j,m,p,nconv;
  PetscScalar       *dst,er,ei;
  PetscReal         sigma;
  const PetscScalar *src,*px;
  Vec               u,v,x,x1,x2,uv;

  PetscFunctionBegin;
  PetscCall(MatGetLocalSize(svd->A,&m,NULL));
  PetscCall(MatGetLocalSize(svd->B,&p,NULL));
  PetscCall(MatCreateVecs(cyclic->C,&x,NULL));
  if (svd->which==SVD_SMALLEST) PetscCall(MatCreateVecsEmpty(svd->B,&x1,&x2));
  else PetscCall(MatCreateVecsEmpty(svd->A,&x2,&x1));
  PetscCall(MatCreateVecs(svd->A,NULL,&u));
  PetscCall(MatCreateVecs(svd->B,NULL,&v));
  PetscCall(EPSGetConverged(cyclic->eps,&nconv));
  for (i=0,j=0;i<nconv;i++) {
    PetscCall(EPSGetEigenpair(cyclic->eps,i,&er,&ei,x,NULL));
    PetscCall(SVDCyclicCheckEigenvalue(svd,er,ei,&sigma,NULL));
    if (sigma<0.0) continue;
    if (svd->which==SVD_SMALLEST) {
      /* evec_i = 1/sqrt(2)*[ v_i; w_i ],  w_i = x_i/c_i */
      PetscCall(VecGetArrayRead(x,&px));
      PetscCall(VecPlaceArray(x2,px));
      PetscCall(VecPlaceArray(x1,px+p));
      PetscCall(VecCopy(x2,v));
      PetscCall(VecScale(v,PETSC_SQRT2));  /* v_i = sqrt(2)*evec_i_1 */
      PetscCall(VecScale(x1,PETSC_SQRT2)); /* w_i = sqrt(2)*evec_i_2 */
      PetscCall(MatMult(svd->A,x1,u));     /* A*w_i = u_i */
      PetscCall(VecScale(x1,1.0/PetscSqrtScalar(1.0+sigma*sigma)));  /* x_i = w_i*c_i */
      PetscCall(BVInsertVec(svd->V,j,x1));
      PetscCall(VecResetArray(x2));
      PetscCall(VecResetArray(x1));
      PetscCall(VecRestoreArrayRead(x,&px));
    } else {
      /* evec_i = 1/sqrt(2)*[ u_i; w_i ],  w_i = x_i/s_i */
      PetscCall(VecGetArrayRead(x,&px));
      PetscCall(VecPlaceArray(x1,px));
      PetscCall(VecPlaceArray(x2,px+m));
      PetscCall(VecCopy(x1,u));
      PetscCall(VecScale(u,PETSC_SQRT2));  /* u_i = sqrt(2)*evec_i_1 */
      PetscCall(VecScale(x2,PETSC_SQRT2)); /* w_i = sqrt(2)*evec_i_2 */
      PetscCall(MatMult(svd->B,x2,v));     /* B*w_i = v_i */
      PetscCall(VecScale(x2,1.0/PetscSqrtScalar(1.0+sigma*sigma)));  /* x_i = w_i*s_i */
      PetscCall(BVInsertVec(svd->V,j,x2));
      PetscCall(VecResetArray(x1));
      PetscCall(VecResetArray(x2));
      PetscCall(VecRestoreArrayRead(x,&px));
    }
    /* copy [u;v] to U[j] */
    PetscCall(BVGetColumn(svd->U,j,&uv));
    PetscCall(VecGetArrayWrite(uv,&dst));
    PetscCall(VecGetArrayRead(u,&src));
    PetscCall(PetscArraycpy(dst,src,m));
    PetscCall(VecRestoreArrayRead(u,&src));
    PetscCall(VecGetArrayRead(v,&src));
    PetscCall(PetscArraycpy(dst+m,src,p));
    PetscCall(VecRestoreArrayRead(v,&src));
    PetscCall(VecRestoreArrayWrite(uv,&dst));
    PetscCall(BVRestoreColumn(svd->U,j,&uv));
    j++;
  }
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&x1));
  PetscCall(VecDestroy(&x2));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_COMPLEX)
/* VecMaxAbs: returns the entry of x that has max(abs(x(i))), using w as a workspace vector */
static PetscErrorCode VecMaxAbs(Vec x,Vec w,PetscScalar *v)
{
  PetscMPIInt       size,rank,root;
  const PetscScalar *xx;
  const PetscInt    *ranges;
  PetscReal         val;
  PetscInt          p;

  PetscFunctionBegin;
  PetscCall(VecCopy(x,w));
  PetscCall(VecAbs(w));
  PetscCall(VecMax(w,&p,&val));
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)x),&size));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)x),&rank));
  PetscCall(VecGetOwnershipRanges(x,&ranges));
  for (root=0;root<size;root++) if (p>=ranges[root] && p<ranges[root+1]) break;
  if (rank==root) {
    PetscCall(VecGetArrayRead(x,&xx));
    *v = xx[p-ranges[root]];
    PetscCall(VecRestoreArrayRead(x,&xx));
  }
  PetscCallMPI(MPI_Bcast(v,1,MPIU_SCALAR,root,PetscObjectComm((PetscObject)x)));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode SVDComputeVectors_Cyclic_Hyperbolic(SVD svd)
{
  SVD_CYCLIC        *cyclic = (SVD_CYCLIC*)svd->data;
  PetscInt          i,j,m,n,nconv;
  PetscScalar       er,ei;
  PetscReal         sigma,nrm;
  PetscBool         isreal;
  const PetscScalar *px;
  Vec               u,x,xi=NULL,x1,x2,x1i=NULL,x2i;
  BV                U=NULL,V=NULL;
#if !defined(PETSC_USE_COMPLEX)
  const PetscScalar *pxi;
  PetscReal         nrmr,nrmi;
#else
  PetscScalar       alpha;
#endif

  PetscFunctionBegin;
  PetscCall(MatCreateVecs(cyclic->C,&x,svd->ishyperbolic?&xi:NULL));
  PetscCall(MatGetLocalSize(svd->A,&m,NULL));
  PetscCall(MatCreateVecsEmpty(svd->OP,&x2,&x1));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(MatCreateVecs(svd->OP,&x2i,&x1i));
#else
  PetscCall(MatCreateVecsEmpty(svd->OP,&x2i,&x1i));
#endif
  /* set-up Omega-normalization of U */
  U = svd->swapped? svd->V: svd->U;
  V = svd->swapped? svd->U: svd->V;
  PetscCall(BVGetSizes(U,&n,NULL,NULL));
  PetscCall(BV_SetMatrixDiagonal(U,svd->omega,svd->A));
  PetscCall(EPSGetConverged(cyclic->eps,&nconv));
  for (i=0,j=0;i<nconv;i++) {
    PetscCall(EPSGetEigenpair(cyclic->eps,i,&er,&ei,x,xi));
    PetscCall(SVDCyclicCheckEigenvalue(svd,er,ei,&sigma,&isreal));
    if (sigma<0.0) continue;
    PetscCall(VecGetArrayRead(x,&px));
    if (svd->swapped) {
      PetscCall(VecPlaceArray(x2,px));
      PetscCall(VecPlaceArray(x1,px+m));
    } else {
      PetscCall(VecPlaceArray(x1,px));
      PetscCall(VecPlaceArray(x2,px+n));
    }
#if defined(PETSC_USE_COMPLEX)
    PetscCall(BVInsertVec(U,j,x1));
    PetscCall(BVInsertVec(V,j,x2));
    if (!isreal) {
      PetscCall(VecMaxAbs(x1,x1i,&alpha));
      PetscCall(BVScaleColumn(U,j,PetscAbsScalar(alpha)/alpha));
      PetscCall(BVScaleColumn(V,j,PetscAbsScalar(alpha)/(alpha*PETSC_i)));
    }
#else
    PetscCall(VecGetArrayRead(xi,&pxi));
    if (svd->swapped) {
      PetscCall(VecPlaceArray(x2i,pxi));
      PetscCall(VecPlaceArray(x1i,pxi+m));
    } else {
      PetscCall(VecPlaceArray(x1i,pxi));
      PetscCall(VecPlaceArray(x2i,pxi+n));
    }
    PetscCall(VecNorm(x2,NORM_2,&nrmr));
    PetscCall(VecNorm(x2i,NORM_2,&nrmi));
    if (nrmi>nrmr) {
      if (isreal) {
        PetscCall(BVInsertVec(U,j,x1i));
        PetscCall(BVInsertVec(V,j,x2i));
      } else {
        PetscCall(BVInsertVec(U,j,x1));
        PetscCall(BVInsertVec(V,j,x2i));
      }
    } else {
      if (isreal) {
        PetscCall(BVInsertVec(U,j,x1));
        PetscCall(BVInsertVec(V,j,x2));
      } else {
        PetscCall(BVInsertVec(U,j,x1i));
        PetscCall(BVScaleColumn(U,j,-1.0));
        PetscCall(BVInsertVec(V,j,x2));
      }
    }
    PetscCall(VecResetArray(x1i));
    PetscCall(VecResetArray(x2i));
    PetscCall(VecRestoreArrayRead(xi,&pxi));
#endif
    PetscCall(VecResetArray(x1));
    PetscCall(VecResetArray(x2));
    PetscCall(VecRestoreArrayRead(x,&px));
    PetscCall(BVGetColumn(U,j,&u));
    PetscCall(VecPointwiseMult(u,u,svd->omega));
    PetscCall(BVRestoreColumn(U,j,&u));
    PetscCall(BVNormColumn(U,j,NORM_2,&nrm));
    PetscCall(BVScaleColumn(U,j,1.0/PetscAbs(nrm)));
    svd->sign[j] = PetscSign(nrm);
    PetscCall(BVNormColumn(V,j,NORM_2,&nrm));
    PetscCall(BVScaleColumn(V,j,1.0/nrm));
    j++;
  }
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&x1));
  PetscCall(VecDestroy(&x2));
  PetscCall(VecDestroy(&xi));
  PetscCall(VecDestroy(&x1i));
  PetscCall(VecDestroy(&x2i));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDComputeVectors_Cyclic(SVD svd)
{
  PetscFunctionBegin;
  switch (svd->problem_type) {
    case SVD_STANDARD:
      PetscCall(SVDComputeVectors_Cyclic_Standard(svd));
      break;
    case SVD_GENERALIZED:
      PetscCall(SVDComputeVectors_Cyclic_Generalized(svd));
      break;
    case SVD_HYPERBOLIC:
      PetscCall(SVDComputeVectors_Cyclic_Hyperbolic(svd));
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_WRONG,"Unknown singular value problem type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSMonitor_Cyclic(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscInt       i,j;
  SVD            svd = (SVD)ctx;
  PetscScalar    er,ei;
  PetscReal      sigma;
  ST             st;

  PetscFunctionBegin;
  nconv = 0;
  PetscCall(EPSGetST(eps,&st));
  for (i=0,j=0;i<PetscMin(nest,svd->ncv);i++) {
    er = eigr[i]; ei = eigi[i];
    PetscCall(STBackTransform(st,1,&er,&ei));
    PetscCall(SVDCyclicCheckEigenvalue(svd,er,ei,&sigma,NULL));
    if (sigma>0.0) {
      svd->sigma[j]  = sigma;
      svd->errest[j] = errest[i];
      if (errest[i] && errest[i] < svd->tol) nconv++;
      j++;
    }
  }
  nest = j;
  PetscCall(SVDMonitor(svd,its,nconv,svd->sigma,svd->errest,nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSetFromOptions_Cyclic(SVD svd,PetscOptionItems *PetscOptionsObject)
{
  PetscBool      set,val;
  SVD_CYCLIC     *cyclic = (SVD_CYCLIC*)svd->data;
  ST             st;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"SVD Cyclic Options");

    PetscCall(PetscOptionsBool("-svd_cyclic_explicitmatrix","Use cyclic explicit matrix","SVDCyclicSetExplicitMatrix",cyclic->explicitmatrix,&val,&set));
    if (set) PetscCall(SVDCyclicSetExplicitMatrix(svd,val));

  PetscOptionsHeadEnd();

  if (!cyclic->eps) PetscCall(SVDCyclicGetEPS(svd,&cyclic->eps));
  if (!cyclic->explicitmatrix && !cyclic->usereps) {
    /* use as default an ST with shell matrix and Jacobi */
    PetscCall(EPSGetST(cyclic->eps,&st));
    PetscCall(STSetMatMode(st,ST_MATMODE_SHELL));
  }
  PetscCall(EPSSetFromOptions(cyclic->eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCyclicSetExplicitMatrix_Cyclic(SVD svd,PetscBool explicitmat)
{
  SVD_CYCLIC *cyclic = (SVD_CYCLIC*)svd->data;

  PetscFunctionBegin;
  if (cyclic->explicitmatrix != explicitmat) {
    cyclic->explicitmatrix = explicitmat;
    svd->state = SVD_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDCyclicSetExplicitMatrix - Indicate if the eigensolver operator
   H(A) = [ 0  A ; A^T 0 ] must be computed explicitly.

   Logically Collective

   Input Parameters:
+  svd         - singular value solver
-  explicitmat - boolean flag indicating if H(A) is built explicitly

   Options Database Key:
.  -svd_cyclic_explicitmatrix <boolean> - Indicates the boolean flag

   Level: advanced

.seealso: SVDCyclicGetExplicitMatrix()
@*/
PetscErrorCode SVDCyclicSetExplicitMatrix(SVD svd,PetscBool explicitmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,explicitmat,2);
  PetscTryMethod(svd,"SVDCyclicSetExplicitMatrix_C",(SVD,PetscBool),(svd,explicitmat));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCyclicGetExplicitMatrix_Cyclic(SVD svd,PetscBool *explicitmat)
{
  SVD_CYCLIC *cyclic = (SVD_CYCLIC*)svd->data;

  PetscFunctionBegin;
  *explicitmat = cyclic->explicitmatrix;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDCyclicGetExplicitMatrix - Returns the flag indicating if H(A) is built explicitly.

   Not Collective

   Input Parameter:
.  svd  - singular value solver

   Output Parameter:
.  explicitmat - the mode flag

   Level: advanced

.seealso: SVDCyclicSetExplicitMatrix()
@*/
PetscErrorCode SVDCyclicGetExplicitMatrix(SVD svd,PetscBool *explicitmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(explicitmat,2);
  PetscUseMethod(svd,"SVDCyclicGetExplicitMatrix_C",(SVD,PetscBool*),(svd,explicitmat));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCyclicSetEPS_Cyclic(SVD svd,EPS eps)
{
  SVD_CYCLIC      *cyclic = (SVD_CYCLIC*)svd->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)eps));
  PetscCall(EPSDestroy(&cyclic->eps));
  cyclic->eps     = eps;
  cyclic->usereps = PETSC_TRUE;
  svd->state      = SVD_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDCyclicSetEPS - Associate an eigensolver object (EPS) to the
   singular value solver.

   Collective

   Input Parameters:
+  svd - singular value solver
-  eps - the eigensolver object

   Level: advanced

.seealso: SVDCyclicGetEPS()
@*/
PetscErrorCode SVDCyclicSetEPS(SVD svd,EPS eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(svd,1,eps,2);
  PetscTryMethod(svd,"SVDCyclicSetEPS_C",(SVD,EPS),(svd,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCyclicGetEPS_Cyclic(SVD svd,EPS *eps)
{
  SVD_CYCLIC     *cyclic = (SVD_CYCLIC*)svd->data;

  PetscFunctionBegin;
  if (!cyclic->eps) {
    PetscCall(EPSCreate(PetscObjectComm((PetscObject)svd),&cyclic->eps));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)cyclic->eps,(PetscObject)svd,1));
    PetscCall(EPSSetOptionsPrefix(cyclic->eps,((PetscObject)svd)->prefix));
    PetscCall(EPSAppendOptionsPrefix(cyclic->eps,"svd_cyclic_"));
    PetscCall(PetscObjectSetOptions((PetscObject)cyclic->eps,((PetscObject)svd)->options));
    PetscCall(EPSSetWhichEigenpairs(cyclic->eps,EPS_LARGEST_REAL));
    PetscCall(EPSMonitorSet(cyclic->eps,EPSMonitor_Cyclic,svd,NULL));
  }
  *eps = cyclic->eps;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDCyclicGetEPS - Retrieve the eigensolver object (EPS) associated
   to the singular value solver.

   Collective

   Input Parameter:
.  svd - singular value solver

   Output Parameter:
.  eps - the eigensolver object

   Level: advanced

.seealso: SVDCyclicSetEPS()
@*/
PetscErrorCode SVDCyclicGetEPS(SVD svd,EPS *eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(eps,2);
  PetscUseMethod(svd,"SVDCyclicGetEPS_C",(SVD,EPS*),(svd,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDView_Cyclic(SVD svd,PetscViewer viewer)
{
  SVD_CYCLIC     *cyclic = (SVD_CYCLIC*)svd->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (!cyclic->eps) PetscCall(SVDCyclicGetEPS(svd,&cyclic->eps));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %s matrix\n",cyclic->explicitmatrix?"explicit":"implicit"));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(EPSView(cyclic->eps,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDReset_Cyclic(SVD svd)
{
  SVD_CYCLIC     *cyclic = (SVD_CYCLIC*)svd->data;

  PetscFunctionBegin;
  PetscCall(EPSReset(cyclic->eps));
  PetscCall(MatDestroy(&cyclic->C));
  PetscCall(MatDestroy(&cyclic->D));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDDestroy_Cyclic(SVD svd)
{
  SVD_CYCLIC     *cyclic = (SVD_CYCLIC*)svd->data;

  PetscFunctionBegin;
  PetscCall(EPSDestroy(&cyclic->eps));
  PetscCall(PetscFree(svd->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCyclicSetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCyclicGetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCyclicSetExplicitMatrix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCyclicGetExplicitMatrix_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_Cyclic(SVD svd)
{
  SVD_CYCLIC     *cyclic;

  PetscFunctionBegin;
  PetscCall(PetscNew(&cyclic));
  svd->data                = (void*)cyclic;
  svd->ops->solve          = SVDSolve_Cyclic;
  svd->ops->solveg         = SVDSolve_Cyclic;
  svd->ops->solveh         = SVDSolve_Cyclic;
  svd->ops->setup          = SVDSetUp_Cyclic;
  svd->ops->setfromoptions = SVDSetFromOptions_Cyclic;
  svd->ops->destroy        = SVDDestroy_Cyclic;
  svd->ops->reset          = SVDReset_Cyclic;
  svd->ops->view           = SVDView_Cyclic;
  svd->ops->computevectors = SVDComputeVectors_Cyclic;
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCyclicSetEPS_C",SVDCyclicSetEPS_Cyclic));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCyclicGetEPS_C",SVDCyclicGetEPS_Cyclic));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCyclicSetExplicitMatrix_C",SVDCyclicSetExplicitMatrix_Cyclic));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCyclicGetExplicitMatrix_C",SVDCyclicGetExplicitMatrix_Cyclic));
  PetscFunctionReturn(PETSC_SUCCESS);
}
