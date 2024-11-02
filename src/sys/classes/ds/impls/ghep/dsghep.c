/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/dsimpl.h>
#include <slepcblaslapack.h>

static PetscErrorCode DSAllocate_GHEP(DS ds,PetscInt ld)
{
  PetscFunctionBegin;
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_B));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Q));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSView_GHEP(DS ds,PetscViewer viewer)
{
  PetscViewerFormat format;

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DSViewMat(ds,viewer,DS_MAT_A));
  PetscCall(DSViewMat(ds,viewer,DS_MAT_B));
  if (ds->state>DS_STATE_INTERMEDIATE) PetscCall(DSViewMat(ds,viewer,DS_MAT_Q));
  if (ds->omat[DS_MAT_X]) PetscCall(DSViewMat(ds,viewer,DS_MAT_X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_GHEP(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscScalar       *Z;
  const PetscScalar *Q;
  PetscInt          ld = ds->ld;

  PetscFunctionBegin;
  PetscCheck(!rnorm,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
  switch (mat) {
    case DS_MAT_X:
    case DS_MAT_Y:
      if (j) {
        PetscCall(MatDenseGetArray(ds->omat[mat],&Z));
        if (ds->state>=DS_STATE_CONDENSED) {
          PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Q],&Q));
          PetscCall(PetscArraycpy(Z+(*j)*ld,Q+(*j)*ld,ld));
          PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_Q],&Q));
        } else {
          PetscCall(PetscArrayzero(Z+(*j)*ld,ld));
          Z[(*j)+(*j)*ld] = 1.0;
        }
        PetscCall(MatDenseRestoreArray(ds->omat[mat],&Z));
      } else {
        if (ds->state>=DS_STATE_CONDENSED) PetscCall(MatCopy(ds->omat[DS_MAT_Q],ds->omat[mat],SAME_NONZERO_PATTERN));
        else PetscCall(DSSetIdentity(ds,mat));
      }
      break;
    case DS_MAT_U:
    case DS_MAT_V:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_GHEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscInt       n,l,i,*perm,ld=ds->ld;
  PetscScalar    *A;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(PETSC_SUCCESS);
  n = ds->n;
  l = ds->l;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  perm = ds->perm;
  for (i=l;i<n;i++) wr[i] = A[i+i*ld];
  if (rr) PetscCall(DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_FALSE));
  else PetscCall(DSSortEigenvalues_Private(ds,wr,NULL,perm,PETSC_FALSE));
  for (i=l;i<n;i++) A[i+i*ld] = wr[perm[i]];
  for (i=l;i<n;i++) wr[i] = A[i+i*ld];
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(DSPermuteColumns_Private(ds,l,n,n,DS_MAT_Q,perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_GHEP(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscScalar    *work,*A,*B,*Q;
  PetscBLASInt   itype = 1,*iwork,info,n1,liwork,ld,lrwork=0,lwork;
  PetscInt       off,i;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork,*rr;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n-ds->l,&n1));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(5*ds->n+3,&liwork));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscBLASIntCast(ds->n*ds->n+2*ds->n,&lwork));
  PetscCall(PetscBLASIntCast(2*ds->n*ds->n+5*ds->n+1+n1,&lrwork));
#else
  PetscCall(PetscBLASIntCast(2*ds->n*ds->n+6*ds->n+1,&lwork));
#endif
  PetscCall(DSAllocateWork_Private(ds,lwork,lrwork,liwork));
  work = ds->work;
  iwork = ds->iwork;
  off = ds->l+ds->l*ld;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
#if defined(PETSC_USE_COMPLEX)
  rr = ds->rwork;
  rwork = ds->rwork+n1;
  PetscCall(PetscBLASIntCast(ds->lrwork-n1,&lrwork));
  PetscCallBLAS("LAPACKsygvd",LAPACKsygvd_(&itype,"V","U",&n1,A+off,&ld,B+off,&ld,rr,work,&lwork,rwork,&lrwork,iwork,&liwork,&info));
  for (i=0;i<n1;i++) wr[ds->l+i] = rr[i];
#else
  PetscCallBLAS("LAPACKsygvd",LAPACKsygvd_(&itype,"V","U",&n1,A+off,&ld,B+off,&ld,wr+ds->l,work,&lwork,iwork,&liwork,&info));
#endif
  SlepcCheckLapackInfo("sygvd",info);
  PetscCall(PetscArrayzero(Q+ds->l*ld,n1*ld));
  for (i=ds->l;i<ds->n;i++) PetscCall(PetscArraycpy(Q+ds->l+i*ld,A+ds->l+i*ld,n1));
  PetscCall(PetscArrayzero(B+ds->l*ld,n1*ld));
  PetscCall(PetscArrayzero(A+ds->l*ld,n1*ld));
  for (i=ds->l;i<ds->n;i++) {
    if (wi) wi[i] = 0.0;
    B[i+i*ld] = 1.0;
    A[i+i*ld] = wr[i];
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_GHEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscScalar    *A,*B,*Q;
  PetscInt       ld=ds->ld,l=ds->l,k;
  PetscMPIInt    n,rank,off=0,size,ldn;

  PetscFunctionBegin;
  k = 2*(ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += (ds->n-l)*ld;
  if (eigr) k += (ds->n-l);
  PetscCall(DSAllocateWork_Private(ds,k,0,0));
  PetscCall(PetscMPIIntCast(k*sizeof(PetscScalar),&size));
  PetscCall(PetscMPIIntCast(ds->n-l,&n));
  PetscCall(PetscMPIIntCast(ld*(ds->n-l),&ldn));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  if (ds->state>DS_STATE_RAW) PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank));
  if (!rank) {
    PetscCallMPI(MPI_Pack(A+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    PetscCallMPI(MPI_Pack(B+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) PetscCallMPI(MPI_Pack(Q+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (eigr) PetscCallMPI(MPI_Pack(eigr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
  }
  PetscCallMPI(MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds)));
  if (rank) {
    PetscCallMPI(MPI_Unpack(ds->work,size,&off,A+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    PetscCallMPI(MPI_Unpack(ds->work,size,&off,B+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) PetscCallMPI(MPI_Unpack(ds->work,size,&off,Q+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (eigr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  if (ds->state>DS_STATE_RAW) PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSHermitian_GHEP(DS ds,DSMatType m,PetscBool *flg)
{
  PetscFunctionBegin;
  if (m==DS_MAT_A || m==DS_MAT_B) *flg = PETSC_TRUE;
  else *flg = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSGHEP - Dense Generalized Hermitian Eigenvalue Problem.

   Level: beginner

   Notes:
   The problem is expressed as A*X = B*X*Lambda, where both A and B are
   real symmetric (or complex Hermitian) and B is positive-definite. Lambda
   is a diagonal matrix whose diagonal elements are the arguments of DSSolve().
   After solve, A is overwritten with Lambda, and B is overwritten with I.

   No intermediate state is implemented, nor compact storage.

   Used DS matrices:
+  DS_MAT_A - first problem matrix
.  DS_MAT_B - second problem matrix
-  DS_MAT_Q - matrix of B-orthogonal eigenvectors, which is equal to X

   Implemented methods:
.  0 - Divide and Conquer (_sygvd)

.seealso: DSCreate(), DSSetType(), DSType
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_GHEP(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate      = DSAllocate_GHEP;
  ds->ops->view          = DSView_GHEP;
  ds->ops->vectors       = DSVectors_GHEP;
  ds->ops->solve[0]      = DSSolve_GHEP;
  ds->ops->sort          = DSSort_GHEP;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize   = DSSynchronize_GHEP;
#endif
  ds->ops->hermitian     = DSHermitian_GHEP;
  PetscFunctionReturn(PETSC_SUCCESS);
}
