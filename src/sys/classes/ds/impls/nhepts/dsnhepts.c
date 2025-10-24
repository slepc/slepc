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

typedef struct {
  PetscScalar *wr,*wi;     /* eigenvalues of B */
} DS_NHEPTS;

static PetscErrorCode DSAllocate_NHEPTS(DS ds,PetscInt ld)
{
  DS_NHEPTS      *ctx = (DS_NHEPTS*)ds->data;

  PetscFunctionBegin;
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_B));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Q));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Z));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscCall(PetscMalloc1(ld,&ctx->wr));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscMalloc1(ld,&ctx->wi));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSView_NHEPTS(DS ds,PetscViewer viewer)
{
  PetscViewerFormat format;

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DSViewMat(ds,viewer,DS_MAT_A));
  PetscCall(DSViewMat(ds,viewer,DS_MAT_B));
  if (ds->state>DS_STATE_INTERMEDIATE) {
    PetscCall(DSViewMat(ds,viewer,DS_MAT_Q));
    PetscCall(DSViewMat(ds,viewer,DS_MAT_Z));
  }
  if (ds->omat[DS_MAT_X]) PetscCall(DSViewMat(ds,viewer,DS_MAT_X));
  if (ds->omat[DS_MAT_Y]) PetscCall(DSViewMat(ds,viewer,DS_MAT_Y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_NHEPTS_Eigen_Some(DS ds,PetscInt *k,PetscReal *rnorm,PetscBool left)
{
  PetscInt          i;
  PetscBLASInt      mm=1,mout,info,ld,n,*select,inc=1,cols=1,zero=0;
  PetscScalar       sone=1.0,szero=0.0;
  PetscReal         norm,done=1.0;
  PetscBool         iscomplex = PETSC_FALSE;
  PetscScalar       *X,*Y;
  const PetscScalar *A,*Q;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(DSAllocateWork_Private(ds,0,0,ld));
  select = ds->iwork;
  for (i=0;i<n;i++) select[i] = (PetscBLASInt)PETSC_FALSE;

  /* compute k-th eigenvector Y of A */
  PetscCall(MatDenseGetArrayRead(ds->omat[left?DS_MAT_B:DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[left?DS_MAT_Y:DS_MAT_X],&X));
  Y = X+(*k)*ld;
  select[*k] = (PetscBLASInt)PETSC_TRUE;
#if !defined(PETSC_USE_COMPLEX)
  if ((*k)<n-1 && A[(*k)+1+(*k)*ld]!=0.0) iscomplex = PETSC_TRUE;
  mm = iscomplex? 2: 1;
  if (iscomplex) select[(*k)+1] = (PetscBLASInt)PETSC_TRUE;
  PetscCall(DSAllocateWork_Private(ds,3*ld,0,0));
  PetscCallBLAS("LAPACKtrevc",LAPACKtrevc_("R","S",select,&n,(PetscScalar*)A,&ld,Y,&ld,Y,&ld,&mm,&mout,ds->work,&info));
#else
  PetscCall(DSAllocateWork_Private(ds,2*ld,ld,0));
  PetscCallBLAS("LAPACKtrevc",LAPACKtrevc_("R","S",select,&n,(PetscScalar*)A,&ld,Y,&ld,Y,&ld,&mm,&mout,ds->work,ds->rwork,&info));
#endif
  SlepcCheckLapackInfo("trevc",info);
  PetscCheck(mout==mm,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Inconsistent arguments");
  PetscCall(MatDenseRestoreArrayRead(ds->omat[left?DS_MAT_B:DS_MAT_A],&A));

  /* accumulate and normalize eigenvectors */
  if (ds->state>=DS_STATE_CONDENSED) {
    PetscCall(MatDenseGetArrayRead(ds->omat[left?DS_MAT_Z:DS_MAT_Q],&Q));
    PetscCall(PetscArraycpy(ds->work,Y,mout*ld));
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,Q,&ld,ds->work,&inc,&szero,Y,&inc));
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,Q,&ld,ds->work+ld,&inc,&szero,Y+ld,&inc));
#endif
    PetscCall(MatDenseRestoreArrayRead(ds->omat[left?DS_MAT_Z:DS_MAT_Q],&Q));
    cols = 1;
    norm = BLASnrm2_(&n,Y,&inc);
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) {
      norm = SlepcAbsEigenvalue(norm,BLASnrm2_(&n,Y+ld,&inc));
      cols = 2;
    }
#endif
    PetscCallBLAS("LAPACKlascl",LAPACKlascl_("G",&zero,&zero,&norm,&done,&n,&cols,Y,&ld,&info));
    SlepcCheckLapackInfo("lascl",info);
  }

  /* set output arguments */
  if (iscomplex) (*k)++;
  if (rnorm) {
    if (iscomplex) *rnorm = SlepcAbsEigenvalue(Y[n-1],Y[n-1+ld]);
    else *rnorm = PetscAbsScalar(Y[n-1]);
  }
  PetscCall(MatDenseRestoreArray(ds->omat[left?DS_MAT_Y:DS_MAT_X],&X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_NHEPTS_Eigen_All(DS ds,PetscBool left)
{
  PetscInt          i;
  PetscBLASInt      n,ld,mout,info,inc=1,cols,zero=0;
  PetscBool         iscomplex;
  PetscScalar       *X;
  const PetscScalar *A;
  PetscReal         norm,done=1.0;
  const char        *back;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(MatDenseGetArrayRead(ds->omat[left?DS_MAT_B:DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayWrite(ds->omat[left?DS_MAT_Y:DS_MAT_X],&X));
  if (ds->state>=DS_STATE_CONDENSED) {
    /* DSSolve() has been called, backtransform with matrix Q */
    back = "B";
    PetscCall(MatCopy(ds->omat[left?DS_MAT_Z:DS_MAT_Q],ds->omat[left?DS_MAT_Y:DS_MAT_X],SAME_NONZERO_PATTERN));
  } else back = "A";
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(DSAllocateWork_Private(ds,3*ld,0,0));
  PetscCallBLAS("LAPACKtrevc",LAPACKtrevc_("R",back,NULL,&n,(PetscScalar*)A,&ld,X,&ld,X,&ld,&n,&mout,ds->work,&info));
#else
  PetscCall(DSAllocateWork_Private(ds,2*ld,ld,0));
  PetscCallBLAS("LAPACKtrevc",LAPACKtrevc_("R",back,NULL,&n,(PetscScalar*)A,&ld,X,&ld,X,&ld,&n,&mout,ds->work,ds->rwork,&info));
#endif
  SlepcCheckLapackInfo("trevc",info);

  /* normalize eigenvectors */
  for (i=0;i<n;i++) {
    iscomplex = (i<n-1 && A[i+1+i*ld]!=0.0)? PETSC_TRUE: PETSC_FALSE;
    cols = 1;
    norm = BLASnrm2_(&n,X+i*ld,&inc);
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) {
      norm = SlepcAbsEigenvalue(norm,BLASnrm2_(&n,X+(i+1)*ld,&inc));
      cols = 2;
    }
#endif
    PetscCallBLAS("LAPACKlascl",LAPACKlascl_("G",&zero,&zero,&norm,&done,&n,&cols,X+i*ld,&ld,&info));
    SlepcCheckLapackInfo("lascl",info);
    if (iscomplex) i++;
  }
  PetscCall(MatDenseRestoreArrayRead(ds->omat[left?DS_MAT_B:DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[left?DS_MAT_Y:DS_MAT_X],&X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_NHEPTS(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_X:
      PetscCheck(!ds->refined,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
      if (j) PetscCall(DSVectors_NHEPTS_Eigen_Some(ds,j,rnorm,PETSC_FALSE));
      else PetscCall(DSVectors_NHEPTS_Eigen_All(ds,PETSC_FALSE));
      break;
    case DS_MAT_Y:
      PetscCheck(!ds->refined,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
      if (j) PetscCall(DSVectors_NHEPTS_Eigen_Some(ds,j,rnorm,PETSC_TRUE));
      else PetscCall(DSVectors_NHEPTS_Eigen_All(ds,PETSC_TRUE));
      break;
    case DS_MAT_U:
    case DS_MAT_V:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_NHEPTS(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  DS_NHEPTS      *ctx = (DS_NHEPTS*)ds->data;
  PetscInt       i,j,cont,id=0,*p,*idx,*idx2;
  PetscReal      s,t;
#if defined(PETSC_USE_COMPLEX)
  Mat            A,U;
#endif

  PetscFunctionBegin;
  PetscCheck(!rr || wr==rr,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
  PetscCall(PetscMalloc3(ds->ld,&idx,ds->ld,&idx2,ds->ld,&p));
  PetscCall(DSSort_NHEP_Total(ds,DS_MAT_A,DS_MAT_Q,wr,wi));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(DSGetMat(ds,DS_MAT_B,&A));
  PetscCall(MatConjugate(A));
  PetscCall(DSRestoreMat(ds,DS_MAT_B,&A));
  PetscCall(DSGetMat(ds,DS_MAT_Z,&U));
  PetscCall(MatConjugate(U));
  PetscCall(DSRestoreMat(ds,DS_MAT_Z,&U));
  for (i=0;i<ds->n;i++) ctx->wr[i] = PetscConj(ctx->wr[i]);
#endif
  PetscCall(DSSort_NHEP_Total(ds,DS_MAT_B,DS_MAT_Z,ctx->wr,ctx->wi));
  /* check correct eigenvalue correspondence */
  cont = 0;
  for (i=0;i<ds->n;i++) {
    if (SlepcAbsEigenvalue(ctx->wr[i]-wr[i],ctx->wi[i]-wi[i])>PETSC_SQRT_MACHINE_EPSILON) {idx2[cont] = i; idx[cont++] = i;}
    p[i] = -1;
  }
  if (cont) {
    for (i=0;i<cont;i++) {
      t = PETSC_MAX_REAL;
      for (j=0;j<cont;j++) if (idx2[j]!=-1 && (s=SlepcAbsEigenvalue(ctx->wr[idx[j]]-wr[idx[i]],ctx->wi[idx[j]]-wi[idx[i]]))<t) { id = j; t = s; }
      p[idx[i]] = idx[id];
      idx2[id] = -1;
    }
    for (i=0;i<ds->n;i++) if (p[i]==-1) p[i] = i;
    PetscCall(DSSortWithPermutation_NHEP_Private(ds,p,DS_MAT_B,DS_MAT_Z,ctx->wr,ctx->wi));
  }
#if defined(PETSC_USE_COMPLEX)
  PetscCall(DSGetMat(ds,DS_MAT_B,&A));
  PetscCall(MatConjugate(A));
  PetscCall(DSRestoreMat(ds,DS_MAT_B,&A));
  PetscCall(DSGetMat(ds,DS_MAT_Z,&U));
  PetscCall(MatConjugate(U));
  PetscCall(DSRestoreMat(ds,DS_MAT_Z,&U));
#endif
  PetscCall(PetscFree3(idx,idx2,p));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSUpdateExtraRow_NHEPTS(DS ds)
{
  PetscInt          i;
  PetscBLASInt      n,ld,incx=1;
  PetscScalar       *A,*x,*y,one=1.0,zero=0.0;
  const PetscScalar *Q;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(DSAllocateWork_Private(ds,2*ld,0,0));
  x = ds->work;
  y = ds->work+ld;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Q],&Q));
  for (i=0;i<n;i++) x[i] = PetscConj(A[n+i*ld]);
  PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&one,Q,&ld,x,&incx,&zero,y,&incx));
  for (i=0;i<n;i++) A[n+i*ld] = PetscConj(y[i]);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&A));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Z],&Q));
  for (i=0;i<n;i++) x[i] = PetscConj(A[n+i*ld]);
  PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&one,Q,&ld,x,&incx,&zero,y,&incx));
  for (i=0;i<n;i++) A[n+i*ld] = PetscConj(y[i]);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&A));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_Z],&Q));
  ds->k = n;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_NHEPTS(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  DS_NHEPTS      *ctx = (DS_NHEPTS*)ds->data;

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscAssertPointer(wi,3);
#endif
  PetscCall(DSSolve_NHEP_Private(ds,DS_MAT_A,DS_MAT_Q,wr,wi));
  PetscCall(DSSolve_NHEP_Private(ds,DS_MAT_B,DS_MAT_Z,ctx->wr,ctx->wi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_NHEPTS(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscInt       ld=ds->ld,l=ds->l,k;
  PetscMPIInt    n,rank,off=0,size,ldn;
  DS_NHEPTS      *ctx = (DS_NHEPTS*)ds->data;
  PetscScalar    *A,*B,*Q,*Z;

  PetscFunctionBegin;
  k = 2*(ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += 2*(ds->n-l)*ld;
  if (eigr) k += ds->n-l;
  if (eigi) k += ds->n-l;
  if (ctx->wr) k += ds->n-l;
  if (ctx->wi) k += ds->n-l;
  PetscCall(DSAllocateWork_Private(ds,k,0,0));
  PetscCall(PetscMPIIntCast(k*sizeof(PetscScalar),&size));
  PetscCall(PetscMPIIntCast(ds->n-l,&n));
  PetscCall(PetscMPIIntCast(ld*(ds->n-l),&ldn));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  if (ds->state>DS_STATE_RAW) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Z],&Z));
  }
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank));
  if (!rank) {
    PetscCallMPI(MPI_Pack(A+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    PetscCallMPI(MPI_Pack(B+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) {
      PetscCallMPI(MPI_Pack(Q+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Pack(Z+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    }
    if (eigr) PetscCallMPI(MPI_Pack(eigr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
#if !defined(PETSC_USE_COMPLEX)
    if (eigi) PetscCallMPI(MPI_Pack(eigi+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
#endif
    if (ctx->wr) PetscCallMPI(MPI_Pack(ctx->wr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (ctx->wi) PetscCallMPI(MPI_Pack(ctx->wi+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
  }
  PetscCallMPI(MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds)));
  if (rank) {
    PetscCallMPI(MPI_Unpack(ds->work,size,&off,A+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    PetscCallMPI(MPI_Unpack(ds->work,size,&off,B+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) {
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,Q+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,Z+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    }
    if (eigr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
#if !defined(PETSC_USE_COMPLEX)
    if (eigi) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigi+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
#endif
    if (ctx->wr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,ctx->wr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (ctx->wi) PetscCallMPI(MPI_Unpack(ds->work,size,&off,ctx->wi+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  if (ds->state>DS_STATE_RAW) {
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Z],&Z));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSGetTruncateSize_NHEPTS(DS ds,PetscInt l,PetscInt n,PetscInt *k)
{
#if !defined(PETSC_USE_COMPLEX)
  const PetscScalar *A,*B;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&B));
  if (A[l+(*k)+(l+(*k)-1)*ds->ld] != 0.0 || B[l+(*k)+(l+(*k)-1)*ds->ld] != 0.0) {
    if (l+(*k)<n-1) (*k)++;
    else (*k)--;
  }
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSTruncate_NHEPTS(DS ds,PetscInt n,PetscBool trim)
{
  PetscInt    i,ld=ds->ld,l=ds->l;
  PetscScalar *A,*B;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
#if defined(PETSC_USE_DEBUG)
  /* make sure diagonal 2x2 block is not broken */
  PetscCheck(ds->state<DS_STATE_CONDENSED || n==0 || n==ds->n || A[n+(n-1)*ld]==0.0 || B[n+(n-1)*ld]==0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"The given size would break a 2x2 block, call DSGetTruncateSize() first");
#endif
  if (trim) {
    if (ds->extrarow) {   /* clean extra row */
      for (i=l;i<ds->n;i++) { A[ds->n+i*ld] = 0.0; B[ds->n+i*ld] = 0.0; }
    }
    ds->l = 0;
    ds->k = 0;
    ds->n = n;
    ds->t = ds->n;   /* truncated length equal to the new dimension */
  } else {
    if (ds->extrarow && ds->k==ds->n) {
      /* copy entries of extra row to the new position, then clean last row */
      for (i=l;i<n;i++) { A[n+i*ld] = A[ds->n+i*ld]; B[n+i*ld] = B[ds->n+i*ld]; }
      for (i=l;i<ds->n;i++) { A[ds->n+i*ld] = 0.0; B[ds->n+i*ld] = 0.0; }
    }
    ds->k = (ds->extrarow)? n: 0;
    ds->t = ds->n;   /* truncated length equal to previous dimension */
    ds->n = n;
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSDestroy_NHEPTS(DS ds)
{
  DS_NHEPTS      *ctx = (DS_NHEPTS*)ds->data;

  PetscFunctionBegin;
  if (ctx->wr) PetscCall(PetscFree(ctx->wr));
  if (ctx->wi) PetscCall(PetscFree(ctx->wi));
  PetscCall(PetscFree(ds->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSMatGetSize_NHEPTS(DS ds,DSMatType t,PetscInt *rows,PetscInt *cols)
{
  PetscFunctionBegin;
  *rows = ((t==DS_MAT_A || t==DS_MAT_B) && ds->extrarow)? ds->n+1: ds->n;
  *cols = ds->n;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSNHEPTS - Dense Non-Hermitian Eigenvalue Problem (special variant intended
   for two-sided Krylov solvers).

   Notes:
   Two related problems are solved, $AX = X\Lambda$ and $BY = Y\Lambda^*$, where $A$ and
   $B$ are supposed to come from the Arnoldi factorizations of a certain matrix and its
   (conjugate) transpose, respectively. Hence, in exact arithmetic the columns of $Y$
   are equal to the left eigenvectors of $A$. $\Lambda$ is a diagonal matrix whose diagonal
   elements are the arguments of `DSSolve()`. After solve, $A$ is overwritten with the
   upper quasi-triangular matrix $T$ of the (real) Schur form, $AQ = QT$, and similarly
   another (real) Schur relation is computed, $BZ = ZS$, overwriting $B$.

   In the intermediate state $A$ and $B$ are reduced to upper Hessenberg form.

   When left eigenvectors `DS_MAT_Y` are requested, right eigenvectors of $B$ are returned,
   while `DS_MAT_X` contains right eigenvectors of $A$.

   Used DS matrices:
+  `DS_MAT_A` - first problem matrix obtained from Arnoldi
.  `DS_MAT_B` - second problem matrix obtained from Arnoldi on the transpose
.  `DS_MAT_Q` - orthogonal/unitary transformation that reduces $A$ to Hessenberg form
   (intermediate step) or matrix of orthogonal Schur vectors of $A$
-  `DS_MAT_Z` - orthogonal/unitary transformation that reduces $B$ to Hessenberg form
   (intermediate step) or matrix of orthogonal Schur vectors of $B$

   Implemented methods:
.  0 - Implicit QR (`_hseqr`)

   Level: beginner

.seealso: [](sec:ds), `DSCreate()`, `DSSetType()`, `DSType`
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_NHEPTS(DS ds)
{
  DS_NHEPTS      *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  ds->data = (void*)ctx;

  ds->ops->allocate        = DSAllocate_NHEPTS;
  ds->ops->view            = DSView_NHEPTS;
  ds->ops->vectors         = DSVectors_NHEPTS;
  ds->ops->solve[0]        = DSSolve_NHEPTS;
  ds->ops->sort            = DSSort_NHEPTS;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize     = DSSynchronize_NHEPTS;
#endif
  ds->ops->gettruncatesize = DSGetTruncateSize_NHEPTS;
  ds->ops->truncate        = DSTruncate_NHEPTS;
  ds->ops->update          = DSUpdateExtraRow_NHEPTS;
  ds->ops->destroy         = DSDestroy_NHEPTS;
  ds->ops->matgetsize      = DSMatGetSize_NHEPTS;
  PetscFunctionReturn(PETSC_SUCCESS);
}
