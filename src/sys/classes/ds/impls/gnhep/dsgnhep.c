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

/*
  1) Patterns of A and B
      DS_STATE_RAW:       DS_STATE_INTERM/CONDENSED
       0       n-1              0       n-1
      -------------            -------------
    0 |* * * * * *|          0 |* * * * * *|
      |* * * * * *|            |  * * * * *|
      |* * * * * *|            |    * * * *|
      |* * * * * *|            |    * * * *|
      |* * * * * *|            |        * *|
  n-1 |* * * * * *|        n-1 |          *|
      -------------            -------------

  2) Moreover, P and Q are assumed to be the identity in DS_STATE_INTERMEDIATE.
*/

static PetscErrorCode CleanDenseSchur(PetscInt n,PetscInt k,PetscScalar *S,PetscInt ldS,PetscScalar *T,PetscInt ldT,PetscScalar *X,PetscInt ldX,PetscScalar *Y,PetscInt ldY);

static PetscErrorCode DSAllocate_GNHEP(DS ds,PetscInt ld)
{
  PetscFunctionBegin;
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_B));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Z));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Q));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSView_GNHEP(DS ds,PetscViewer viewer)
{
  PetscViewerFormat format;

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DSViewMat(ds,viewer,DS_MAT_A));
  PetscCall(DSViewMat(ds,viewer,DS_MAT_B));
  if (ds->state>DS_STATE_INTERMEDIATE) {
    PetscCall(DSViewMat(ds,viewer,DS_MAT_Z));
    PetscCall(DSViewMat(ds,viewer,DS_MAT_Q));
  }
  if (ds->omat[DS_MAT_X]) PetscCall(DSViewMat(ds,viewer,DS_MAT_X));
  if (ds->omat[DS_MAT_Y]) PetscCall(DSViewMat(ds,viewer,DS_MAT_Y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_GNHEP_Eigen_Some(DS ds,PetscInt *k,PetscReal *rnorm,PetscBool left)
{
  PetscInt       i;
  PetscBLASInt   n,ld,mout,info,*select,mm,inc=1,cols=1,zero=0;
  PetscScalar    *X,*Y,*XY,*Z,*Q,*A,*B,fone=1.0,fzero=0.0;
  PetscReal      norm,done=1.0;
  PetscBool      iscomplex = PETSC_FALSE;
  const char     *side;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  if (left) {
    X = NULL;
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Y],&Y));
    side = "L";
  } else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
    Y = NULL;
    side = "R";
  }
  XY = left? Y: X;
  PetscCall(DSAllocateWork_Private(ds,0,0,ld));
  select = ds->iwork;
  for (i=0;i<n;i++) select[i] = (PetscBLASInt)PETSC_FALSE;
  if (ds->state <= DS_STATE_INTERMEDIATE) {
    PetscCall(DSSetIdentity(ds,DS_MAT_Q));
    PetscCall(DSSetIdentity(ds,DS_MAT_Z));
  }
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Z],&Z));
  PetscCall(CleanDenseSchur(n,0,A,ld,B,ld,Q,ld,Z,ld));
  if (ds->state < DS_STATE_CONDENSED) PetscCall(DSSetState(ds,DS_STATE_CONDENSED));

  /* compute k-th eigenvector */
  select[*k] = (PetscBLASInt)PETSC_TRUE;
#if defined(PETSC_USE_COMPLEX)
  mm = 1;
  PetscCall(DSAllocateWork_Private(ds,2*ld,2*ld,0));
  PetscCallBLAS("LAPACKtgevc",LAPACKtgevc_(side,"S",select,&n,A,&ld,B,&ld,PetscSafePointerPlusOffset(Y,(*k)*ld),&ld,PetscSafePointerPlusOffset(X,(*k)*ld),&ld,&mm,&mout,ds->work,ds->rwork,&info));
#else
  if ((*k)<n-1 && (A[ld*(*k)+(*k)+1] != 0.0 || B[ld*(*k)+(*k)+1] != 0.0)) iscomplex = PETSC_TRUE;
  mm = iscomplex? 2: 1;
  if (iscomplex) select[(*k)+1] = (PetscBLASInt)PETSC_TRUE;
  PetscCall(DSAllocateWork_Private(ds,6*ld,0,0));
  PetscCallBLAS("LAPACKtgevc",LAPACKtgevc_(side,"S",select,&n,A,&ld,B,&ld,PetscSafePointerPlusOffset(Y,(*k)*ld),&ld,PetscSafePointerPlusOffset(X,(*k)*ld),&ld,&mm,&mout,ds->work,&info));
#endif
  SlepcCheckLapackInfo("tgevc",info);
  PetscCheck(select[*k] && mout==mm,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong arguments in call to Lapack xTGEVC");
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));

  /* accumulate and normalize eigenvectors */
  PetscCall(PetscArraycpy(ds->work,XY+(*k)*ld,mm*ld));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&mm,&n,&fone,left?Z:Q,&ld,ds->work,&ld,&fzero,XY+(*k)*ld,&ld));
  norm = BLASnrm2_(&n,XY+(*k)*ld,&inc);
#if !defined(PETSC_USE_COMPLEX)
  if (iscomplex) {
    norm = SlepcAbsEigenvalue(norm,BLASnrm2_(&n,XY+(*k+1)*ld,&inc));
    cols = 2;
  }
#endif
  PetscCallBLAS("LAPACKlascl",LAPACKlascl_("G",&zero,&zero,&norm,&done,&n,&cols,XY+(*k)*ld,&ld,&info));
  SlepcCheckLapackInfo("lascl",info);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Z],&Z));

  /* set output arguments */
  if (rnorm) {
    if (iscomplex) *rnorm = SlepcAbsEigenvalue(XY[n-1+(*k)*ld],XY[n-1+(*k+1)*ld]);
    else *rnorm = PetscAbsScalar(XY[n-1+(*k)*ld]);
  }
  if (iscomplex) (*k)++;
  PetscCall(MatDenseRestoreArray(ds->omat[left?DS_MAT_Y:DS_MAT_X],&XY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_GNHEP_Eigen_All(DS ds,PetscBool left)
{
  PetscInt       i;
  PetscBLASInt   n,ld,mout,info,inc = 1;
  PetscBool      iscomplex;
  PetscScalar    *X,*Y,*XY,*Q,*Z,*A,*B,tmp;
  PetscReal      norm;
  const char     *side,*back;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  if (left) {
    X = NULL;
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Y],&Y));
    side = "L";
  } else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
    Y = NULL;
    side = "R";
  }
  XY = left? Y: X;
  if (ds->state <= DS_STATE_INTERMEDIATE) {
    PetscCall(DSSetIdentity(ds,DS_MAT_Q));
    PetscCall(DSSetIdentity(ds,DS_MAT_Z));
  }
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Z],&Z));
  PetscCall(CleanDenseSchur(n,0,A,ld,B,ld,Q,ld,Z,ld));
  if (ds->state>=DS_STATE_CONDENSED) {
    /* DSSolve() has been called, backtransform with matrix Q */
    back = "B";
    PetscCall(PetscArraycpy(left?Y:X,left?Z:Q,ld*ld));
  } else {
    back = "A";
    PetscCall(DSSetState(ds,DS_STATE_CONDENSED));
  }
#if defined(PETSC_USE_COMPLEX)
  PetscCall(DSAllocateWork_Private(ds,2*ld,2*ld,0));
  PetscCallBLAS("LAPACKtgevc",LAPACKtgevc_(side,back,NULL,&n,A,&ld,B,&ld,Y,&ld,X,&ld,&n,&mout,ds->work,ds->rwork,&info));
#else
  PetscCall(DSAllocateWork_Private(ds,6*ld,0,0));
  PetscCallBLAS("LAPACKtgevc",LAPACKtgevc_(side,back,NULL,&n,A,&ld,B,&ld,Y,&ld,X,&ld,&n,&mout,ds->work,&info));
#endif
  SlepcCheckLapackInfo("tgevc",info);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Z],&Z));

  /* normalize eigenvectors */
  for (i=0;i<n;i++) {
    iscomplex = (i<n-1 && (A[i+1+i*ld]!=0.0 || B[i+1+i*ld]!=0.0))? PETSC_TRUE: PETSC_FALSE;
    norm = BLASnrm2_(&n,XY+i*ld,&inc);
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) {
      tmp = BLASnrm2_(&n,XY+(i+1)*ld,&inc);
      norm = SlepcAbsEigenvalue(norm,tmp);
    }
#endif
    tmp = 1.0 / norm;
    PetscCallBLAS("BLASscal",BLASscal_(&n,&tmp,XY+i*ld,&inc));
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) PetscCallBLAS("BLASscal",BLASscal_(&n,&tmp,XY+(i+1)*ld,&inc));
#endif
    if (iscomplex) i++;
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArray(ds->omat[left?DS_MAT_Y:DS_MAT_X],&XY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_GNHEP(DS ds,DSMatType mat,PetscInt *k,PetscReal *rnorm)
{
  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_X:
    case DS_MAT_Y:
      if (k) PetscCall(DSVectors_GNHEP_Eigen_Some(ds,k,rnorm,mat == DS_MAT_Y?PETSC_TRUE:PETSC_FALSE));
      else PetscCall(DSVectors_GNHEP_Eigen_All(ds,mat == DS_MAT_Y?PETSC_TRUE:PETSC_FALSE));
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_GNHEP_Arbitrary(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscInt       i;
  PetscBLASInt   info,n,ld,mout,lwork,liwork,*iwork,*selection,zero_=0,true_=1;
  PetscScalar    *S,*T,*Q,*Z,*work,*beta;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
#if !defined(PETSC_USE_COMPLEX)
  lwork = 4*n+16;
#else
  lwork = 1;
#endif
  liwork = 1;
  PetscCall(DSAllocateWork_Private(ds,lwork+2*n,0,liwork+n));
  beta      = ds->work;
  work      = ds->work + n;
  PetscCall(PetscBLASIntCast(ds->lwork-n,&lwork));
  selection = ds->iwork;
  iwork     = ds->iwork + n;
  PetscCall(PetscBLASIntCast(ds->liwork-n,&liwork));
  /* Compute the selected eigenvalue to be in the leading position */
  PetscCall(DSSortEigenvalues_Private(ds,rr,ri,ds->perm,PETSC_FALSE));
  PetscCall(PetscArrayzero(selection,n));
  for (i=0; i<*k; i++) selection[ds->perm[i]] = 1;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&S));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&T));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Z],&Z));
#if !defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKtgsen",LAPACKtgsen_(&zero_,&true_,&true_,selection,&n,S,&ld,T,&ld,wr,wi,beta,Z,&ld,Q,&ld,&mout,NULL,NULL,NULL,work,&lwork,iwork,&liwork,&info));
#else
  PetscCallBLAS("LAPACKtgsen",LAPACKtgsen_(&zero_,&true_,&true_,selection,&n,S,&ld,T,&ld,wr,beta,Z,&ld,Q,&ld,&mout,NULL,NULL,NULL,work,&lwork,iwork,&liwork,&info));
#endif
  SlepcCheckLapackInfo("tgsen",info);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&S));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&T));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Z],&Z));
  *k = mout;
  for (i=0;i<n;i++) {
    if (beta[i]==0.0) wr[i] = (PetscRealPart(wr[i])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wr[i] /= beta[i];
#if !defined(PETSC_USE_COMPLEX)
    if (beta[i]==0.0) wi[i] = (wi[i]>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wi[i] /= beta[i];
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_GNHEP_Total(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscScalar    re;
  PetscInt       i,j,pos,result;
  PetscBLASInt   ifst,ilst,info,n,ld,one=1;
  PetscScalar    *S,*T,*Z,*Q;
#if !defined(PETSC_USE_COMPLEX)
  PetscBLASInt   lwork;
  PetscScalar    *work,a,safmin,scale1,scale2,im;
#endif

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&S));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&T));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Z],&Z));
#if !defined(PETSC_USE_COMPLEX)
  lwork = -1;
  PetscCallBLAS("LAPACKtgexc",LAPACKtgexc_(&one,&one,&ld,NULL,&ld,NULL,&ld,NULL,&ld,NULL,&ld,&one,&one,&a,&lwork,&info));
  SlepcCheckLapackInfo("tgexc",info);
  safmin = LAPACKlamch_("S");
  PetscCall(PetscBLASIntCast((PetscInt)a,&lwork));
  PetscCall(DSAllocateWork_Private(ds,lwork,0,0));
  work = ds->work;
#endif
  /* selection sort */
  for (i=ds->l;i<n-1;i++) {
    re = wr[i];
#if !defined(PETSC_USE_COMPLEX)
    im = wi[i];
#endif
    pos = 0;
    j = i+1; /* j points to the next eigenvalue */
#if !defined(PETSC_USE_COMPLEX)
    if (im != 0) j=i+2;
#endif
    /* find minimum eigenvalue */
    for (;j<n;j++) {
#if !defined(PETSC_USE_COMPLEX)
      PetscCall(SlepcSCCompare(ds->sc,re,im,wr[j],wi[j],&result));
#else
      PetscCall(SlepcSCCompare(ds->sc,re,0.0,wr[j],0.0,&result));
#endif
      if (result > 0) {
        re = wr[j];
#if !defined(PETSC_USE_COMPLEX)
        im = wi[j];
#endif
        pos = j;
      }
#if !defined(PETSC_USE_COMPLEX)
      if (wi[j] != 0) j++;
#endif
    }
    if (pos) {
      /* interchange blocks */
      PetscCall(PetscBLASIntCast(pos+1,&ifst));
      PetscCall(PetscBLASIntCast(i+1,&ilst));
#if !defined(PETSC_USE_COMPLEX)
      PetscCallBLAS("LAPACKtgexc",LAPACKtgexc_(&one,&one,&n,S,&ld,T,&ld,Z,&ld,Q,&ld,&ifst,&ilst,work,&lwork,&info));
#else
      PetscCallBLAS("LAPACKtgexc",LAPACKtgexc_(&one,&one,&n,S,&ld,T,&ld,Z,&ld,Q,&ld,&ifst,&ilst,&info));
#endif
      SlepcCheckLapackInfo("tgexc",info);
      /* recover original eigenvalues from T and S matrices */
      for (j=i;j<n;j++) {
#if !defined(PETSC_USE_COMPLEX)
        if (j<n-1 && S[j*ld+j+1] != 0.0) {
          /* complex conjugate eigenvalue */
          PetscCallBLAS("LAPACKlag2",LAPACKlag2_(S+j*ld+j,&ld,T+j*ld+j,&ld,&safmin,&scale1,&scale2,&re,&a,&im));
          wr[j] = re / scale1;
          wi[j] = im / scale1;
          wr[j+1] = a / scale2;
          wi[j+1] = -wi[j];
          j++;
        } else
#endif
        {
          if (T[j*ld+j] == 0.0) wr[j] = (PetscRealPart(S[j*ld+j])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
          else wr[j] = S[j*ld+j] / T[j*ld+j];
#if !defined(PETSC_USE_COMPLEX)
          wi[j] = 0.0;
#endif
        }
      }
    }
#if !defined(PETSC_USE_COMPLEX)
    if (wi[i] != 0.0) i++;
#endif
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&S));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&T));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Z],&Z));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_GNHEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscFunctionBegin;
  if (!rr || wr == rr) PetscCall(DSSort_GNHEP_Total(ds,wr,wi));
  else PetscCall(DSSort_GNHEP_Arbitrary(ds,wr,wi,rr,ri,k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSUpdateExtraRow_GNHEP(DS ds)
{
  PetscInt          i;
  PetscBLASInt      n,ld,incx=1;
  PetscScalar       *A,*B,*x,*y,one=1.0,zero=0.0;
  const PetscScalar *Q;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Q],&Q));
  PetscCall(DSAllocateWork_Private(ds,2*ld,0,0));
  x = ds->work;
  y = ds->work+ld;
  for (i=0;i<n;i++) x[i] = PetscConj(A[n+i*ld]);
  PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&one,Q,&ld,x,&incx,&zero,y,&incx));
  for (i=0;i<n;i++) A[n+i*ld] = PetscConj(y[i]);
  for (i=0;i<n;i++) x[i] = PetscConj(B[n+i*ld]);
  PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&one,Q,&ld,x,&incx,&zero,y,&incx));
  for (i=0;i<n;i++) B[n+i*ld] = PetscConj(y[i]);
  ds->k = n;
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_Q],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Write zeros from the column k to n in the lower triangular part of the
   matrices S and T, and inside 2-by-2 diagonal blocks of T in order to
   make (S,T) a valid Schur decompositon.
*/
static PetscErrorCode CleanDenseSchur(PetscInt n,PetscInt k,PetscScalar *S,PetscInt ldS,PetscScalar *T,PetscInt ldT,PetscScalar *X,PetscInt ldX,PetscScalar *Y,PetscInt ldY)
{
  PetscInt       i;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
  PetscScalar    s;
#else
  PetscBLASInt   ldS_,ldT_,n_i,n_i_2,one=1,n_,i_2,i_;
  PetscScalar    b11,b22,sr,cr,sl,cl;
#endif

  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  for (i=k; i<n; i++) {
    /* Some functions need the diagonal elements in T be real */
    if (T && PetscImaginaryPart(T[ldT*i+i]) != 0.0) {
      s = PetscConj(T[ldT*i+i])/PetscAbsScalar(T[ldT*i+i]);
      for (j=0;j<=i;j++) {
        T[ldT*i+j] *= s;
        S[ldS*i+j] *= s;
      }
      T[ldT*i+i] = PetscRealPart(T[ldT*i+i]);
      if (X) for (j=0;j<n;j++) X[ldX*i+j] *= s;
    }
    j = i+1;
    if (j<n) {
      S[ldS*i+j] = 0.0;
      if (T) T[ldT*i+j] = 0.0;
    }
  }
#else
  PetscCall(PetscBLASIntCast(ldS,&ldS_));
  PetscCall(PetscBLASIntCast(ldT,&ldT_));
  PetscCall(PetscBLASIntCast(n,&n_));
  for (i=k;i<n-1;i++) {
    if (S[ldS*i+i+1] != 0.0) {
      /* Check if T(i+1,i) and T(i,i+1) are zero */
      if (T[ldT*(i+1)+i] != 0.0 || T[ldT*i+i+1] != 0.0) {
        /* Check if T(i+1,i) and T(i,i+1) are negligible */
        if (PetscAbs(T[ldT*(i+1)+i])+PetscAbs(T[ldT*i+i+1]) < (PetscAbs(T[ldT*i+i])+PetscAbs(T[ldT*(i+1)+i+1]))*PETSC_MACHINE_EPSILON) {
          T[ldT*i+i+1] = 0.0;
          T[ldT*(i+1)+i] = 0.0;
        } else {
          /* If one of T(i+1,i) or T(i,i+1) is negligible, we make zero the other element */
          if (PetscAbs(T[ldT*i+i+1]) < (PetscAbs(T[ldT*i+i])+PetscAbs(T[ldT*(i+1)+i+1])+PetscAbs(T[ldT*(i+1)+i]))*PETSC_MACHINE_EPSILON) {
            PetscCallBLAS("LAPACKlasv2",LAPACKlasv2_(&T[ldT*i+i],&T[ldT*(i+1)+i],&T[ldT*(i+1)+i+1],&b22,&b11,&sl,&cl,&sr,&cr));
          } else if (PetscAbs(T[ldT*(i+1)+i]) < (PetscAbs(T[ldT*i+i])+PetscAbs(T[ldT*(i+1)+i+1])+PetscAbs(T[ldT*i+i+1]))*PETSC_MACHINE_EPSILON) {
            PetscCallBLAS("LAPACKlasv2",LAPACKlasv2_(&T[ldT*i+i],&T[ldT*i+i+1],&T[ldT*(i+1)+i+1],&b22,&b11,&sr,&cr,&sl,&cl));
          } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unsupported format. Call DSSolve before this function");
          PetscCall(PetscBLASIntCast(n-i,&n_i));
          n_i_2 = n_i - 2;
          PetscCall(PetscBLASIntCast(i+2,&i_2));
          PetscCall(PetscBLASIntCast(i,&i_));
          if (b11 < 0.0) {
            cr = -cr; sr = -sr;
            b11 = -b11; b22 = -b22;
          }
          PetscCallBLAS("BLASrot",BLASrot_(&n_i,&S[ldS*i+i],&ldS_,&S[ldS*i+i+1],&ldS_,&cl,&sl));
          PetscCallBLAS("BLASrot",BLASrot_(&i_2,&S[ldS*i],&one,&S[ldS*(i+1)],&one,&cr,&sr));
          PetscCallBLAS("BLASrot",BLASrot_(&n_i_2,&T[ldT*(i+2)+i],&ldT_,&T[ldT*(i+2)+i+1],&ldT_,&cl,&sl));
          PetscCallBLAS("BLASrot",BLASrot_(&i_,&T[ldT*i],&one,&T[ldT*(i+1)],&one,&cr,&sr));
          if (X) PetscCallBLAS("BLASrot",BLASrot_(&n_,&X[ldX*i],&one,&X[ldX*(i+1)],&one,&cr,&sr));
          if (Y) PetscCallBLAS("BLASrot",BLASrot_(&n_,&Y[ldY*i],&one,&Y[ldY*(i+1)],&one,&cl,&sl));
          T[ldT*i+i] = b11; T[ldT*i+i+1] = 0.0;
          T[ldT*(i+1)+i] = 0.0; T[ldT*(i+1)+i+1] = b22;
        }
      }
      i++;
    }
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_GNHEP(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscScalar    *work,*beta,a;
  PetscInt       i;
  PetscBLASInt   lwork,info,n,ld,iaux;
  PetscScalar    *A,*B,*Z,*Q;
  PetscBool      usegges3=(ds->method==1)?PETSC_TRUE:PETSC_FALSE;

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscAssertPointer(wi,3);
#endif
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  lwork = -1;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Z],&Z));
#if !defined(PETSC_USE_COMPLEX)
  if (usegges3) PetscCallBLAS("LAPACKgges3",LAPACKgges3_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,wi,NULL,Z,&ld,Q,&ld,&a,&lwork,NULL,&info));
  else PetscCallBLAS("LAPACKgges",LAPACKgges_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,wi,NULL,Z,&ld,Q,&ld,&a,&lwork,NULL,&info));
  PetscCall(PetscBLASIntCast((PetscInt)a,&lwork));
  PetscCall(DSAllocateWork_Private(ds,lwork+ld,0,0));
  beta = ds->work;
  work = beta+ds->n;
  PetscCall(PetscBLASIntCast(ds->lwork-ds->n,&lwork));
  if (usegges3) PetscCallBLAS("LAPACKgges3",LAPACKgges3_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,wi,beta,Z,&ld,Q,&ld,work,&lwork,NULL,&info));
  else PetscCallBLAS("LAPACKgges",LAPACKgges_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,wi,beta,Z,&ld,Q,&ld,work,&lwork,NULL,&info));
#else
  if (usegges3) PetscCallBLAS("LAPACKgges3",LAPACKgges3_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,NULL,Z,&ld,Q,&ld,&a,&lwork,NULL,NULL,&info));
  else PetscCallBLAS("LAPACKgges",LAPACKgges_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,NULL,Z,&ld,Q,&ld,&a,&lwork,NULL,NULL,&info));
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork));
  PetscCall(DSAllocateWork_Private(ds,lwork+ld,8*ld,0));
  beta = ds->work;
  work = beta+ds->n;
  PetscCall(PetscBLASIntCast(ds->lwork-ds->n,&lwork));
  if (usegges3) PetscCallBLAS("LAPACKgges3",LAPACKgges3_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,beta,Z,&ld,Q,&ld,work,&lwork,ds->rwork,NULL,&info));
  else PetscCallBLAS("LAPACKgges",LAPACKgges_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,beta,Z,&ld,Q,&ld,work,&lwork,ds->rwork,NULL,&info));
#endif
  SlepcCheckLapackInfo(usegges3?"gges3":"gges",info);
  for (i=0;i<n;i++) {
    if (beta[i]==0.0) wr[i] = (PetscRealPart(wr[i])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wr[i] /= beta[i];
#if !defined(PETSC_USE_COMPLEX)
    if (beta[i]==0.0) wi[i] = (wi[i]>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wi[i] /= beta[i];
#else
    if (wi) wi[i] = 0.0;
#endif
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Z],&Z));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_GNHEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscInt       ld=ds->ld,l=ds->l,k;
  PetscMPIInt    n,rank,off=0,size,ldn;
  PetscScalar    *A,*B,*Q,*Z;

  PetscFunctionBegin;
  k = 2*(ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += 2*(ds->n-l)*ld;
  if (eigr) k += (ds->n-l);
  if (eigi) k += (ds->n-l);
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

static PetscErrorCode DSTruncate_GNHEP(DS ds,PetscInt n,PetscBool trim)
{
  PetscInt    i,ld=ds->ld,l=ds->l;
  PetscScalar *A,*B;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
#if defined(PETSC_USE_DEBUG)
  /* make sure diagonal 2x2 block is not broken */
  PetscCheck(ds->state<DS_STATE_CONDENSED || n==0 || n==ds->n || (A[n+(n-1)*ld]==0.0 && B[n+(n-1)*ld]==0.0),PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"The given size would break a 2x2 block, call DSGetTruncateSize() first");
#endif
  if (trim) {
    if (ds->extrarow) {   /* clean extra row */
      for (i=l;i<ds->n;i++) A[ds->n+i*ld] = 0.0;
      for (i=l;i<ds->n;i++) B[ds->n+i*ld] = 0.0;
    }
    ds->l = 0;
    ds->k = 0;
    ds->n = n;
    ds->t = ds->n;   /* truncated length equal to the new dimension */
  } else {
    if (ds->extrarow && ds->k==ds->n) {
      /* copy entries of extra row to the new position, then clean last row */
      for (i=l;i<n;i++) A[n+i*ld] = A[ds->n+i*ld];
      for (i=l;i<ds->n;i++) A[ds->n+i*ld] = 0.0;
      for (i=l;i<n;i++) B[n+i*ld] = B[ds->n+i*ld];
      for (i=l;i<ds->n;i++) B[ds->n+i*ld] = 0.0;
    }
    ds->k = (ds->extrarow)? n: 0;
    ds->t = ds->n;   /* truncated length equal to previous dimension */
    ds->n = n;
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSGNHEP - Dense Generalized Non-Hermitian Eigenvalue Problem.

   Level: beginner

   Notes:
   The problem is expressed as A*X = B*X*Lambda, where (A,B) is the input
   matrix pencil. Lambda is a diagonal matrix whose diagonal elements are the
   arguments of DSSolve(). After solve, (A,B) is overwritten with the
   generalized (real) Schur form (S,T) = (Z'*A*Q,Z'*B*Q), with the first
   matrix being upper quasi-triangular and the second one triangular.

   Used DS matrices:
+  DS_MAT_A - first problem matrix
.  DS_MAT_B - second problem matrix
.  DS_MAT_Q - first orthogonal/unitary transformation that reduces to
   generalized (real) Schur form
-  DS_MAT_Z - second orthogonal/unitary transformation that reduces to
   generalized (real) Schur form

   Implemented methods:
+  0 - QZ iteration (_gges)
-  1 - blocked QZ iteration (_gges3, if available)

.seealso: DSCreate(), DSSetType(), DSType
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_GNHEP(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate        = DSAllocate_GNHEP;
  ds->ops->view            = DSView_GNHEP;
  ds->ops->vectors         = DSVectors_GNHEP;
  ds->ops->solve[0]        = DSSolve_GNHEP;
#if !defined(SLEPC_MISSING_LAPACK_GGES3)
  ds->ops->solve[1]        = DSSolve_GNHEP;
#endif
  ds->ops->sort            = DSSort_GNHEP;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize     = DSSynchronize_GNHEP;
#endif
  ds->ops->gettruncatesize = DSGetTruncateSize_Default;
  ds->ops->truncate        = DSTruncate_GNHEP;
  ds->ops->update          = DSUpdateExtraRow_GNHEP;
  PetscFunctionReturn(PETSC_SUCCESS);
}
