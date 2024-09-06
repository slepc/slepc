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

static PetscErrorCode DSAllocate_GHIEP(DS ds,PetscInt ld)
{
  PetscFunctionBegin;
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_B));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Q));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_T));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_D));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DSSwitchFormat_GHIEP(DS ds,PetscBool tocompact)
{
  PetscReal      *T,*S;
  PetscScalar    *A,*B;
  PetscInt       i,n,ld;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&S));
  n = ds->n;
  ld = ds->ld;
  if (tocompact) { /* switch from dense (arrow) to compact storage */
    PetscCall(PetscArrayzero(T,n));
    PetscCall(PetscArrayzero(T+ld,n));
    PetscCall(PetscArrayzero(T+2*ld,n));
    PetscCall(PetscArrayzero(S,n));
    for (i=0;i<n-1;i++) {
      T[i]    = PetscRealPart(A[i+i*ld]);
      T[ld+i] = PetscRealPart(A[i+1+i*ld]);
      S[i]    = PetscRealPart(B[i+i*ld]);
    }
    T[n-1] = PetscRealPart(A[n-1+(n-1)*ld]);
    S[n-1] = PetscRealPart(B[n-1+(n-1)*ld]);
    for (i=ds->l;i<ds->k;i++) T[2*ld+i] = PetscRealPart(A[ds->k+i*ld]);
  } else { /* switch from compact (arrow) to dense storage */
    for (i=0;i<n;i++) {
      PetscCall(PetscArrayzero(A+i*ld,n));
      PetscCall(PetscArrayzero(B+i*ld,n));
    }
    for (i=0;i<n-1;i++) {
      A[i+i*ld]     = T[i];
      A[i+1+i*ld]   = T[ld+i];
      A[i+(i+1)*ld] = T[ld+i];
      B[i+i*ld]     = S[i];
    }
    A[n-1+(n-1)*ld] = T[n-1];
    B[n-1+(n-1)*ld] = S[n-1];
    for (i=ds->l;i<ds->k;i++) {
      A[ds->k+i*ld] = T[2*ld+i];
      A[i+ds->k*ld] = T[2*ld+i];
    }
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&S));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSView_GHIEP(DS ds,PetscViewer viewer)
{
  PetscViewerFormat format;
  PetscInt          i,j;
  PetscReal         *T,*S,value;
  const char        *methodname[] = {
                     "QR + Inverse Iteration",
                     "HZ method",
                     "QR"
  };
  const int         nmeth=PETSC_STATIC_ARRAY_LENGTH(methodname);

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
    if (ds->method<nmeth) PetscCall(PetscViewerASCIIPrintf(viewer,"solving the problem with: %s\n",methodname[ds->method]));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSGetArrayReal(ds,DS_MAT_D,&S));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%% Size = %" PetscInt_FMT " %" PetscInt_FMT "\n",ds->n,ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = zeros(%" PetscInt_FMT ",3);\n",3*ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = [\n"));
      for (i=0;i<ds->n;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,i+1,(double)T[i]));
      for (i=0;i<ds->n-1;i++) {
        if (T[i+ds->ld] !=0 && i!=ds->k-1) {
          PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+2,i+1,(double)T[i+ds->ld]));
          PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,i+2,(double)T[i+ds->ld]));
        }
      }
      for (i = ds->l;i<ds->k;i++) {
        if (T[i+2*ds->ld]) {
          PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",ds->k+1,i+1,(double)T[i+2*ds->ld]));
          PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,ds->k+1,(double)T[i+2*ds->ld]));
        }
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(zzz);\n",DSMatName[DS_MAT_A]));

      PetscCall(PetscViewerASCIIPrintf(viewer,"%% Size = %" PetscInt_FMT " %" PetscInt_FMT "\n",ds->n,ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"omega = zeros(%" PetscInt_FMT ",3);\n",3*ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"omega = [\n"));
      for (i=0;i<ds->n;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,i+1,(double)S[i]));
      PetscCall(PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(omega);\n",DSMatName[DS_MAT_B]));

    } else {
      PetscCall(PetscViewerASCIIPrintf(viewer,"T\n"));
      for (i=0;i<ds->n;i++) {
        for (j=0;j<ds->n;j++) {
          if (i==j) value = T[i];
          else if (i==j+1 || j==i+1) value = T[PetscMin(i,j)+ds->ld];
          else if ((i<ds->k && j==ds->k) || (i==ds->k && j<ds->k)) value = T[PetscMin(i,j)+2*ds->ld];
          else value = 0.0;
          PetscCall(PetscViewerASCIIPrintf(viewer," %18.16e ",(double)value));
        }
        PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"omega\n"));
      for (i=0;i<ds->n;i++) {
        for (j=0;j<ds->n;j++) {
          if (i==j) value = S[i];
          else value = 0.0;
          PetscCall(PetscViewerASCIIPrintf(viewer," %18.16e ",(double)value));
        }
        PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      }
    }
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerFlush(viewer));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&S));
  } else {
    PetscCall(DSViewMat(ds,viewer,DS_MAT_A));
    PetscCall(DSViewMat(ds,viewer,DS_MAT_B));
  }
  if (ds->state>DS_STATE_INTERMEDIATE) PetscCall(DSViewMat(ds,viewer,DS_MAT_Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_GHIEP_Eigen_Some(DS ds,PetscInt *idx,PetscReal *rnorm)
{
  PetscReal         b[4],M[4],*T,*S,d1,d2,s1,s2,e,scal1,scal2,wr1,wr2,wi,ep,norm;
  PetscScalar       *X,Y[4],alpha,szero=0.0;
  const PetscScalar *A,*B,*Q;
  PetscInt          k;
  PetscBLASInt      two=2,n_,ld,one=1;
#if !defined(PETSC_USE_COMPLEX)
  PetscBLASInt      four=4;
#endif

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&S));
  k = *idx;
  PetscCall(PetscBLASIntCast(ds->n,&n_));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  if (k < ds->n-1) e = (ds->compact)?T[k+ld]:PetscRealPart(A[(k+1)+ld*k]);
  else e = 0.0;
  if (e == 0.0) { /* Real */
    if (ds->state>=DS_STATE_CONDENSED) PetscCall(PetscArraycpy(X+k*ld,Q+k*ld,ld));
    else {
      PetscCall(PetscArrayzero(X+k*ds->ld,ds->ld));
      X[k+k*ds->ld] = 1.0;
    }
    if (rnorm) *rnorm = PetscAbsScalar(X[ds->n-1+k*ld]);
  } else { /* 2x2 block */
    if (ds->compact) {
      s1 = S[k];
      d1 = T[k];
      s2 = S[k+1];
      d2 = T[k+1];
    } else {
      s1 = PetscRealPart(B[k*ld+k]);
      d1 = PetscRealPart(A[k+k*ld]);
      s2 = PetscRealPart(B[(k+1)*ld+k+1]);
      d2 = PetscRealPart(A[k+1+(k+1)*ld]);
    }
    M[0] = d1; M[1] = e; M[2] = e; M[3]= d2;
    b[0] = s1; b[1] = 0.0; b[2] = 0.0; b[3] = s2;
    ep = LAPACKlamch_("S");
    /* Compute eigenvalues of the block */
    PetscCallBLAS("LAPACKlag2",LAPACKlag2_(M,&two,b,&two,&ep,&scal1,&scal2,&wr1,&wr2,&wi));
    PetscCheck(wi!=0.0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Real block in DSVectors_GHIEP");
    /* Complex eigenvalues */
    PetscCheck(scal1>=ep,PETSC_COMM_SELF,PETSC_ERR_FP,"Nearly infinite eigenvalue");
    wr1 /= scal1;
    wi  /= scal1;
#if !defined(PETSC_USE_COMPLEX)
    if (SlepcAbs(s1*d1-wr1,wi)<SlepcAbs(s2*d2-wr1,wi)) {
      Y[0] = wr1-s2*d2; Y[1] = s2*e; Y[2] = wi; Y[3] = 0.0;
    } else {
      Y[0] = s1*e; Y[1] = wr1-s1*d1; Y[2] = 0.0; Y[3] = wi;
    }
    norm = BLASnrm2_(&four,Y,&one);
    norm = 1.0/norm;
    if (ds->state >= DS_STATE_CONDENSED) {
      alpha = norm;
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&two,&two,&alpha,Q+k*ld,&ld,Y,&two,&szero,X+k*ld,&ld));
      if (rnorm) *rnorm = SlepcAbsEigenvalue(X[ds->n-1+k*ld],X[ds->n-1+(k+1)*ld]);
    } else {
      PetscCall(PetscArrayzero(X+k*ld,2*ld));
      X[k*ld+k]       = Y[0]*norm;
      X[k*ld+k+1]     = Y[1]*norm;
      X[(k+1)*ld+k]   = Y[2]*norm;
      X[(k+1)*ld+k+1] = Y[3]*norm;
    }
#else
    if (SlepcAbs(s1*d1-wr1,wi)<SlepcAbs(s2*d2-wr1,wi)) {
      Y[0] = PetscCMPLX(wr1-s2*d2,wi);
      Y[1] = s2*e;
    } else {
      Y[0] = s1*e;
      Y[1] = PetscCMPLX(wr1-s1*d1,wi);
    }
    norm = BLASnrm2_(&two,Y,&one);
    norm = 1.0/norm;
    if (ds->state >= DS_STATE_CONDENSED) {
      alpha = norm;
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&two,&alpha,Q+k*ld,&ld,Y,&one,&szero,X+k*ld,&one));
      if (rnorm) *rnorm = PetscAbsScalar(X[ds->n-1+k*ld]);
    } else {
      PetscCall(PetscArrayzero(X+k*ld,2*ld));
      X[k*ld+k]   = Y[0]*norm;
      X[k*ld+k+1] = Y[1]*norm;
    }
    X[(k+1)*ld+k]   = PetscConj(X[k*ld+k]);
    X[(k+1)*ld+k+1] = PetscConj(X[k*ld+k+1]);
#endif
    (*idx)++;
  }
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&S));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_GHIEP(DS ds,DSMatType mat,PetscInt *k,PetscReal *rnorm)
{
  PetscScalar       *Z;
  const PetscScalar *A,*Q;
  PetscInt          i;
  PetscReal         e,*T,*S;

  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_X:
    case DS_MAT_Y:
      if (k) PetscCall(DSVectors_GHIEP_Eigen_Some(ds,k,rnorm));
      else {
        PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
        PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Q],&Q));
        PetscCall(MatDenseGetArray(ds->omat[mat],&Z));
        PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
        for (i=0; i<ds->n; i++) {
          e = (ds->compact)?T[i+ds->ld]:PetscRealPart(A[(i+1)+ds->ld*i]);
          if (e == 0.0) { /* real */
            if (ds->state >= DS_STATE_CONDENSED) PetscCall(PetscArraycpy(Z+i*ds->ld,Q+i*ds->ld,ds->ld));
            else {
              PetscCall(PetscArrayzero(Z+i*ds->ld,ds->ld));
              Z[i+i*ds->ld] = 1.0;
            }
          } else {
            PetscCall(DSVectors_GHIEP_Eigen_Some(ds,&i,rnorm));
            PetscCall(DSGetArrayReal(ds,DS_MAT_D,&S));
            S[i]=S[i-1]=1.0;
            PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&S));
          }
        }
        PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
        PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_Q],&Q));
        PetscCall(MatDenseRestoreArray(ds->omat[mat],&Z));
        PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
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

/*
  Extract the eigenvalues contained in the block-diagonal of the indefinite problem.
  Only the index range n0..n1 is processed.
*/
PetscErrorCode DSGHIEPComplexEigs(DS ds,PetscInt n0,PetscInt n1,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt          k,ld;
  PetscBLASInt      two=2;
  const PetscScalar *A,*B;
  PetscReal         *D,*T,b[4],M[4],d1,d2,s1,s2,e,scal1,scal2,ep,wr1,wr2,wi1;

  PetscFunctionBegin;
  ld = ds->ld;
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&D));
  for (k=n0;k<n1;k++) {
    if (k < n1-1) e = (ds->compact)?T[ld+k]:PetscRealPart(A[(k+1)+ld*k]);
    else e = 0.0;
    if (e==0.0) { /* real eigenvalue */
      wr[k] = (ds->compact)?T[k]/D[k]:A[k+k*ld]/B[k+k*ld];
#if !defined(PETSC_USE_COMPLEX)
      wi[k] = 0.0 ;
#endif
    } else { /* diagonal block */
      if (ds->compact) {
        s1 = D[k];
        d1 = T[k];
        s2 = D[k+1];
        d2 = T[k+1];
      } else {
        s1 = PetscRealPart(B[k*ld+k]);
        d1 = PetscRealPart(A[k+k*ld]);
        s2 = PetscRealPart(B[(k+1)*ld+k+1]);
        d2 = PetscRealPart(A[k+1+(k+1)*ld]);
      }
      M[0] = d1; M[1] = e; M[2] = e; M[3]= d2;
      b[0] = s1; b[1] = 0.0; b[2] = 0.0; b[3] = s2;
      ep = LAPACKlamch_("S");
      /* Compute eigenvalues of the block */
      PetscCallBLAS("LAPACKlag2",LAPACKlag2_(M,&two,b,&two,&ep,&scal1,&scal2,&wr1,&wr2,&wi1));
      PetscCheck(scal1>=ep,PETSC_COMM_SELF,PETSC_ERR_FP,"Nearly infinite eigenvalue");
      if (wi1==0.0) { /* Real eigenvalues */
        PetscCheck(scal2>=ep,PETSC_COMM_SELF,PETSC_ERR_FP,"Nearly infinite eigenvalue");
        wr[k] = wr1/scal1; wr[k+1] = wr2/scal2;
#if !defined(PETSC_USE_COMPLEX)
        wi[k] = wi[k+1] = 0.0;
#endif
      } else { /* Complex eigenvalues */
#if !defined(PETSC_USE_COMPLEX)
        wr[k]   = wr1/scal1;
        wr[k+1] = wr[k];
        wi[k]   = wi1/scal1;
        wi[k+1] = -wi[k];
#else
        wr[k]   = PetscCMPLX(wr1,wi1)/scal1;
        wr[k+1] = PetscConj(wr[k]);
#endif
      }
      k++;
    }
  }
#if defined(PETSC_USE_COMPLEX)
  if (wi) {
    for (k=n0;k<n1;k++) wi[k] = 0.0;
  }
#endif
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&D));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_GHIEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscInt       n,i,*perm;
  PetscReal      *d,*e,*s;

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscAssertPointer(wi,3);
#endif
  n = ds->n;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d + ds->ld;
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&s));
  PetscCall(DSAllocateWork_Private(ds,ds->ld,ds->ld,0));
  perm = ds->perm;
  if (!rr) {
    rr = wr;
    ri = wi;
  }
  PetscCall(DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_TRUE));
  if (!ds->compact) PetscCall(DSSwitchFormat_GHIEP(ds,PETSC_TRUE));
  PetscCall(PetscArraycpy(ds->work,wr,n));
  for (i=ds->l;i<n;i++) wr[i] = *(ds->work+perm[i]);
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscArraycpy(ds->work,wi,n));
  for (i=ds->l;i<n;i++) wi[i] = *(ds->work+perm[i]);
#endif
  PetscCall(PetscArraycpy(ds->rwork,s,n));
  for (i=ds->l;i<n;i++) s[i] = *(ds->rwork+perm[i]);
  PetscCall(PetscArraycpy(ds->rwork,d,n));
  for (i=ds->l;i<n;i++) d[i] = *(ds->rwork+perm[i]);
  PetscCall(PetscArraycpy(ds->rwork,e,n-1));
  PetscCall(PetscArrayzero(e+ds->l,n-1-ds->l));
  for (i=ds->l;i<n-1;i++) {
    if (perm[i]<n-1) e[i] = *(ds->rwork+perm[i]);
  }
  if (!ds->compact) PetscCall(DSSwitchFormat_GHIEP(ds,PETSC_FALSE));
  PetscCall(DSPermuteColumns_Private(ds,ds->l,n,n,DS_MAT_Q,perm));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSUpdateExtraRow_GHIEP(DS ds)
{
  PetscInt          i;
  PetscBLASInt      n,ld,incx=1;
  PetscScalar       *A,*x,*y,one=1.0,zero=0.0;
  const PetscScalar *Q;
  PetscReal         *T,*b,*r,beta;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Q],&Q));
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    b = T+ld;
    r = T+2*ld;
    beta = b[n-1];   /* in compact, we assume all entries are zero except the last one */
    for (i=0;i<n;i++) r[i] = PetscRealPart(beta*Q[n-1+i*ld]);
    ds->k = n;
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  } else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    PetscCall(DSAllocateWork_Private(ds,2*ld,0,0));
    x = ds->work;
    y = ds->work+ld;
    for (i=0;i<n;i++) x[i] = PetscConj(A[n+i*ld]);
    PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&one,Q,&ld,x,&incx,&zero,y,&incx));
    for (i=0;i<n;i++) A[n+i*ld] = PetscConj(y[i]);
    ds->k = n;
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_Q],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Get eigenvectors with inverse iteration.
  The system matrix is in Hessenberg form.
*/
PetscErrorCode DSGHIEPInverseIteration(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt          i,off;
  PetscBLASInt      *select,*infoC,ld,n1,mout,info;
  const PetscScalar *A,*B;
  PetscScalar       *H,*X;
  PetscReal         *s,*d,*e;
#if defined(PETSC_USE_COMPLEX)
  PetscInt          j;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(ds->n-ds->l,&n1));
  PetscCall(DSAllocateWork_Private(ds,ld*ld+2*ld,ld,2*ld));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_W));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_W],&H));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&s));
  e = d + ld;
  select = ds->iwork;
  infoC = ds->iwork + ld;
  off = ds->l+ds->l*ld;
  if (ds->compact) {
    H[off] = d[ds->l]*s[ds->l];
    H[off+ld] = e[ds->l]*s[ds->l];
    for (i=ds->l+1;i<ds->n-1;i++) {
      H[i+(i-1)*ld] = e[i-1]*s[i];
      H[i+i*ld] = d[i]*s[i];
      H[i+(i+1)*ld] = e[i]*s[i];
    }
    H[ds->n-1+(ds->n-2)*ld] = e[ds->n-2]*s[ds->n-1];
    H[ds->n-1+(ds->n-1)*ld] = d[ds->n-1]*s[ds->n-1];
  } else {
    s[ds->l]  = PetscRealPart(B[off]);
    H[off]    = A[off]*s[ds->l];
    H[off+ld] = A[off+ld]*s[ds->l];
    for (i=ds->l+1;i<ds->n-1;i++) {
      s[i] = PetscRealPart(B[i+i*ld]);
      H[i+(i-1)*ld] = A[i+(i-1)*ld]*s[i];
      H[i+i*ld]     = A[i+i*ld]*s[i];
      H[i+(i+1)*ld] = A[i+(i+1)*ld]*s[i];
    }
    s[ds->n-1] = PetscRealPart(B[ds->n-1+(ds->n-1)*ld]);
    H[ds->n-1+(ds->n-2)*ld] = A[ds->n-1+(ds->n-2)*ld]*s[ds->n-1];
    H[ds->n-1+(ds->n-1)*ld] = A[ds->n-1+(ds->n-1)*ld]*s[ds->n-1];
  }
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_X));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
  for (i=0;i<n1;i++) select[i] = 1;
#if !defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKhsein",LAPACKhsein_("R","N","N",select,&n1,H+off,&ld,wr+ds->l,wi+ds->l,NULL,&ld,X+off,&ld,&n1,&mout,ds->work,NULL,infoC,&info));
#else
  PetscCallBLAS("LAPACKhsein",LAPACKhsein_("R","N","N",select,&n1,H+off,&ld,wr+ds->l,NULL,&ld,X+off,&ld,&n1,&mout,ds->work,ds->rwork,NULL,infoC,&info));

  /* Separate real and imaginary part of complex eigenvectors */
  for (j=ds->l;j<ds->n;j++) {
    if (PetscAbsReal(PetscImaginaryPart(wr[j])) > PetscAbsScalar(wr[j])*PETSC_SQRT_MACHINE_EPSILON) {
      for (i=ds->l;i<ds->n;i++) {
        X[i+(j+1)*ds->ld] = PetscImaginaryPart(X[i+j*ds->ld]);
        X[i+j*ds->ld] = PetscRealPart(X[i+j*ds->ld]);
      }
      j++;
    }
  }
#endif
  SlepcCheckLapackInfo("hsein",info);
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_W],&H));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));
  PetscCall(DSGHIEPOrthogEigenv(ds,DS_MAT_X,wr,wi,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Undo 2x2 blocks that have real eigenvalues.
*/
PetscErrorCode DSGHIEPRealBlocks(DS ds)
{
  PetscInt       i;
  PetscReal      e,d1,d2,s1,s2,ss1,ss2,t,dd,ss;
  PetscReal      maxy,ep,scal1,scal2,snorm;
  PetscReal      *T,*D,b[4],M[4],wr1,wr2,wi;
  PetscScalar    *A,*B,*Q,Y[4],sone=1.0,szero=0.0;
  PetscBLASInt   m,two=2,ld;
  PetscBool      isreal;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(ds->n-ds->l,&m));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&D));
  PetscCall(DSAllocateWork_Private(ds,2*m,0,0));
  for (i=ds->l;i<ds->n-1;i++) {
    e = (ds->compact)?T[ld+i]:PetscRealPart(A[(i+1)+ld*i]);
    if (e != 0.0) { /* 2x2 block */
      if (ds->compact) {
        s1 = D[i];
        d1 = T[i];
        s2 = D[i+1];
        d2 = T[i+1];
      } else {
        s1 = PetscRealPart(B[i*ld+i]);
        d1 = PetscRealPart(A[i*ld+i]);
        s2 = PetscRealPart(B[(i+1)*ld+i+1]);
        d2 = PetscRealPart(A[(i+1)*ld+i+1]);
      }
      isreal = PETSC_FALSE;
      if (s1==s2) { /* apply a Jacobi rotation to compute the eigendecomposition */
        dd = d1-d2;
        if (2*PetscAbsReal(e) <= dd) {
          t = 2*e/dd;
          t = t/(1 + PetscSqrtReal(1+t*t));
        } else {
          t = dd/(2*e);
          ss = (t>=0)?1.0:-1.0;
          t = ss/(PetscAbsReal(t)+PetscSqrtReal(1+t*t));
        }
        Y[0] = 1/PetscSqrtReal(1 + t*t); Y[3] = Y[0]; /* c */
        Y[1] = Y[0]*t; Y[2] = -Y[1]; /* s */
        wr1 = d1+t*e; wr2 = d2-t*e;
        ss1 = s1; ss2 = s2;
        isreal = PETSC_TRUE;
      } else {
        ss1 = 1.0; ss2 = 1.0,
        M[0] = d1; M[1] = e; M[2] = e; M[3]= d2;
        b[0] = s1; b[1] = 0.0; b[2] = 0.0; b[3] = s2;
        ep = LAPACKlamch_("S");

        /* Compute eigenvalues of the block */
        PetscCallBLAS("LAPACKlag2",LAPACKlag2_(M,&two,b,&two,&ep,&scal1,&scal2,&wr1,&wr2,&wi));
        if (wi==0.0) { /* Real eigenvalues */
          isreal = PETSC_TRUE;
          PetscCheck(scal1>=ep && scal2>=ep,PETSC_COMM_SELF,PETSC_ERR_FP,"Nearly infinite eigenvalue");
          wr1 /= scal1;
          wr2 /= scal2;
          if (PetscAbsReal(s1*d1-wr1)<PetscAbsReal(s2*d2-wr1)) {
            Y[0] = wr1-s2*d2;
            Y[1] = s2*e;
          } else {
            Y[0] = s1*e;
            Y[1] = wr1-s1*d1;
          }
          /* normalize with a signature*/
          maxy = PetscMax(PetscAbsScalar(Y[0]),PetscAbsScalar(Y[1]));
          scal1 = PetscRealPart(Y[0])/maxy;
          scal2 = PetscRealPart(Y[1])/maxy;
          snorm = scal1*scal1*s1 + scal2*scal2*s2;
          if (snorm<0) { ss1 = -1.0; snorm = -snorm; }
          snorm = maxy*PetscSqrtReal(snorm);
          Y[0] = Y[0]/snorm;
          Y[1] = Y[1]/snorm;
          if (PetscAbsReal(s1*d1-wr2)<PetscAbsReal(s2*d2-wr2)) {
            Y[2] = wr2-s2*d2;
            Y[3] = s2*e;
          } else {
            Y[2] = s1*e;
            Y[3] = wr2-s1*d1;
          }
          maxy = PetscMax(PetscAbsScalar(Y[2]),PetscAbsScalar(Y[3]));
          scal1 = PetscRealPart(Y[2])/maxy;
          scal2 = PetscRealPart(Y[3])/maxy;
          snorm = scal1*scal1*s1 + scal2*scal2*s2;
          if (snorm<0) { ss2 = -1.0; snorm = -snorm; }
          snorm = maxy*PetscSqrtReal(snorm); Y[2] = Y[2]/snorm; Y[3] = Y[3]/snorm;
        }
        wr1 *= ss1; wr2 *= ss2;
      }
      if (isreal) {
        if (ds->compact) {
          D[i]    = ss1;
          T[i]    = wr1;
          D[i+1]  = ss2;
          T[i+1]  = wr2;
          T[ld+i] = 0.0;
        } else {
          B[i*ld+i]       = ss1;
          A[i*ld+i]       = wr1;
          B[(i+1)*ld+i+1] = ss2;
          A[(i+1)*ld+i+1] = wr2;
          A[(i+1)+ld*i]   = 0.0;
          A[i+ld*(i+1)]   = 0.0;
        }
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&m,&two,&two,&sone,Q+ds->l+i*ld,&ld,Y,&two,&szero,ds->work,&m));
        PetscCall(PetscArraycpy(Q+ds->l+i*ld,ds->work,m));
        PetscCall(PetscArraycpy(Q+ds->l+(i+1)*ld,ds->work+m,m));
      }
      i++;
    }
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&D));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_GHIEP_QR_II(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt          i,off;
  PetscBLASInt      n1,ld,one=1,info,lwork;
  const PetscScalar *A,*B;
  PetscScalar       *H,*Q;
  PetscReal         *d,*e,*s;
#if defined(PETSC_USE_COMPLEX)
  PetscInt          j;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscAssertPointer(wi,3);
#endif
  PetscCall(PetscBLASIntCast(ds->n-ds->l,&n1));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  off = ds->l + ds->l*ld;
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&s));
  e = d + ld;
#if defined(PETSC_USE_DEBUG)
  /* Check signature */
  for (i=0;i<ds->n;i++) {
    PetscReal de = (ds->compact)?s[i]:PetscRealPart(B[i*ld+i]);
    PetscCheck(de==1.0 || de==-1.0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Diagonal elements of the signature matrix must be 1 or -1");
  }
#endif

  /* Quick return if possible */
  if (n1 == 1) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
    for (i=0;i<=ds->l;i++) Q[i+i*ld] = 1.0;
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
    PetscCall(DSGHIEPComplexEigs(ds,0,ds->l,wr,wi));
    if (!ds->compact) {
      d[ds->l] = PetscRealPart(A[off]);
      s[ds->l] = PetscRealPart(B[off]);
    }
    PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
    PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
    wr[ds->l] = d[ds->l]/s[ds->l];
    if (wi) wi[ds->l] = 0.0;
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscCall(DSAllocateWork_Private(ds,ld*ld,2*ld,ld*2));
  lwork = ld*ld;

  /* Reduce to pseudotriadiagonal form */
  PetscCall(DSIntermediate_GHIEP(ds));

  /* Compute Eigenvalues (QR) */
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_W));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_W],&H));
  if (ds->compact) {
    H[off]    = d[ds->l]*s[ds->l];
    H[off+ld] = e[ds->l]*s[ds->l];
    for (i=ds->l+1;i<ds->n-1;i++) {
      H[i+(i-1)*ld] = e[i-1]*s[i];
      H[i+i*ld]     = d[i]*s[i];
      H[i+(i+1)*ld] = e[i]*s[i];
    }
    H[ds->n-1+(ds->n-2)*ld] = e[ds->n-2]*s[ds->n-1];
    H[ds->n-1+(ds->n-1)*ld] = d[ds->n-1]*s[ds->n-1];
  } else {
    s[ds->l]  = PetscRealPart(B[off]);
    H[off]    = A[off]*s[ds->l];
    H[off+ld] = A[off+ld]*s[ds->l];
    for (i=ds->l+1;i<ds->n-1;i++) {
      s[i] = PetscRealPart(B[i+i*ld]);
      H[i+(i-1)*ld] = A[i+(i-1)*ld]*s[i];
      H[i+i*ld]     = A[i+i*ld]*s[i];
      H[i+(i+1)*ld] = A[i+(i+1)*ld]*s[i];
    }
    s[ds->n-1] = PetscRealPart(B[ds->n-1+(ds->n-1)*ld]);
    H[ds->n-1+(ds->n-2)*ld] = A[ds->n-1+(ds->n-2)*ld]*s[ds->n-1];
    H[ds->n-1+(ds->n-1)*ld] = A[ds->n-1+(ds->n-1)*ld]*s[ds->n-1];
  }

#if !defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKhseqr",LAPACKhseqr_("E","N",&n1,&one,&n1,H+off,&ld,wr+ds->l,wi+ds->l,NULL,&ld,ds->work,&lwork,&info));
#else
  PetscCallBLAS("LAPACKhseqr",LAPACKhseqr_("E","N",&n1,&one,&n1,H+off,&ld,wr+ds->l,NULL,&ld,ds->work,&lwork,&info));
  for (i=ds->l;i<ds->n;i++) if (PetscAbsReal(PetscImaginaryPart(wr[i]))<10*PETSC_MACHINE_EPSILON) wr[i] = PetscRealPart(wr[i]);
  /* Sort to have consecutive conjugate pairs */
  for (i=ds->l;i<ds->n;i++) {
    j=i+1;
    while (j<ds->n && (PetscAbsScalar(wr[i]-PetscConj(wr[j]))>PetscAbsScalar(wr[i])*PETSC_SQRT_MACHINE_EPSILON)) j++;
    if (j==ds->n) {
      PetscCheck(PetscAbsReal(PetscImaginaryPart(wr[i]))<PetscAbsScalar(wr[i])*PETSC_SQRT_MACHINE_EPSILON,PETSC_COMM_SELF,PETSC_ERR_LIB,"Found complex without conjugate pair");
      wr[i] = PetscRealPart(wr[i]);
    } else { /* complex eigenvalue */
      wr[j] = wr[i+1];
      if (PetscImaginaryPart(wr[i])<0) wr[i] = PetscConj(wr[i]);
      wr[i+1] = PetscConj(wr[i]);
      i++;
    }
  }
#endif
  SlepcCheckLapackInfo("hseqr",info);
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_W],&H));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));

  /* Compute Eigenvectors with Inverse Iteration */
  PetscCall(DSGHIEPInverseIteration(ds,wr,wi));

  /* Recover eigenvalues from diagonal */
  PetscCall(DSGHIEPComplexEigs(ds,0,ds->l,wr,wi));
#if defined(PETSC_USE_COMPLEX)
  if (wi) {
    for (i=ds->l;i<ds->n;i++) wi[i] = 0.0;
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_GHIEP_QR(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt          i,j,off,nwu=0,n,lw,lwr,nwru=0;
  PetscBLASInt      n_,ld,info,lwork,ilo,ihi;
  const PetscScalar *A,*B;
  PetscScalar       *H,*Q,*X;
  PetscReal         *d,*s,*scale,nrm,*rcde,*rcdv;
#if defined(PETSC_USE_COMPLEX)
  PetscInt          k;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscAssertPointer(wi,3);
#endif
  n = ds->n-ds->l;
  PetscCall(PetscBLASIntCast(n,&n_));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  off = ds->l + ds->l*ld;
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&s));
#if defined(PETSC_USE_DEBUG)
  /* Check signature */
  for (i=0;i<ds->n;i++) {
    PetscReal de = (ds->compact)?s[i]:PetscRealPart(B[i*ld+i]);
    PetscCheck(de==1.0 || de==-1.0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Diagonal elements of the signature matrix must be 1 or -1");
  }
#endif

  /* Quick return if possible */
  if (n_ == 1) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
    for (i=0;i<=ds->l;i++) Q[i+i*ld] = 1.0;
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
    PetscCall(DSGHIEPComplexEigs(ds,0,ds->l,wr,wi));
    if (!ds->compact) {
      d[ds->l] = PetscRealPart(A[off]);
      s[ds->l] = PetscRealPart(B[off]);
    }
    PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
    PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
    wr[ds->l] = d[ds->l]/s[ds->l];
    if (wi) wi[ds->l] = 0.0;
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  lw = 14*ld+ld*ld;
  lwr = 7*ld;
  PetscCall(DSAllocateWork_Private(ds,lw,lwr,0));
  scale = ds->rwork+nwru;
  nwru += ld;
  rcde = ds->rwork+nwru;
  nwru += ld;
  rcdv = ds->rwork+nwru;

  /* Form pseudo-symmetric matrix */
  H =  ds->work+nwu;
  nwu += n*n;
  PetscCall(PetscArrayzero(H,n*n));
  if (ds->compact) {
    for (i=0;i<n-1;i++) {
      H[i+i*n]     = s[ds->l+i]*d[ds->l+i];
      H[i+1+i*n]   = s[ds->l+i+1]*d[ld+ds->l+i];
      H[i+(i+1)*n] = s[ds->l+i]*d[ld+ds->l+i];
    }
    H[n-1+(n-1)*n] = s[ds->l+n-1]*d[ds->l+n-1];
    for (i=0;i<ds->k-ds->l;i++) {
      H[ds->k-ds->l+i*n] = s[ds->k]*d[2*ld+ds->l+i];
      H[i+(ds->k-ds->l)*n] = s[i+ds->l]*d[2*ld+ds->l+i];
    }
  } else {
    for (j=0;j<n;j++) {
      for (i=0;i<n;i++) H[i+j*n] = B[off+i+i*ld]*A[off+i+j*ld];
    }
  }

  /* Compute eigenpairs */
  PetscCall(PetscBLASIntCast(lw-nwu,&lwork));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_X));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_X],&X));
#if !defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgeevx",LAPACKgeevx_("B","N","V","N",&n_,H,&n_,wr+ds->l,wi+ds->l,NULL,&ld,X+off,&ld,&ilo,&ihi,scale,&nrm,rcde,rcdv,ds->work+nwu,&lwork,NULL,&info));
#else
  PetscCallBLAS("LAPACKgeevx",LAPACKgeevx_("B","N","V","N",&n_,H,&n_,wr+ds->l,NULL,&ld,X+off,&ld,&ilo,&ihi,scale,&nrm,rcde,rcdv,ds->work+nwu,&lwork,ds->rwork+nwru,&info));

  /* Sort to have consecutive conjugate pairs
     Separate real and imaginary part of complex eigenvectors*/
  for (i=ds->l;i<ds->n;i++) {
    j=i+1;
    while (j<ds->n && (PetscAbsScalar(wr[i]-PetscConj(wr[j]))>PetscAbsScalar(wr[i])*PETSC_SQRT_MACHINE_EPSILON)) j++;
    if (j==ds->n) {
      PetscCheck(PetscAbsReal(PetscImaginaryPart(wr[i]))<PetscAbsScalar(wr[i])*PETSC_SQRT_MACHINE_EPSILON,PETSC_COMM_SELF,PETSC_ERR_LIB,"Found complex without conjugate pair");
      wr[i]=PetscRealPart(wr[i]); /* real eigenvalue */
      for (k=ds->l;k<ds->n;k++) {
        X[k+i*ds->ld] = PetscRealPart(X[k+i*ds->ld]);
      }
    } else { /* complex eigenvalue */
      if (j!=i+1) {
        wr[j] = wr[i+1];
        PetscCall(PetscArraycpy(X+j*ds->ld,X+(i+1)*ds->ld,ds->ld));
      }
      if (PetscImaginaryPart(wr[i])<0) {
        wr[i] = PetscConj(wr[i]);
        for (k=ds->l;k<ds->n;k++) {
          X[k+(i+1)*ds->ld] = -PetscImaginaryPart(X[k+i*ds->ld]);
          X[k+i*ds->ld] = PetscRealPart(X[k+i*ds->ld]);
        }
      } else {
        for (k=ds->l;k<ds->n;k++) {
          X[k+(i+1)*ds->ld] = PetscImaginaryPart(X[k+i*ds->ld]);
          X[k+i*ds->ld] = PetscRealPart(X[k+i*ds->ld]);
        }
      }
      wr[i+1] = PetscConj(wr[i]);
      i++;
    }
  }
#endif
  SlepcCheckLapackInfo("geevx",info);
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_X],&X));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));

  /* Compute real s-orthonormal basis */
  PetscCall(DSGHIEPOrthogEigenv(ds,DS_MAT_X,wr,wi,PETSC_FALSE));

  /* Recover eigenvalues from diagonal */
  PetscCall(DSGHIEPComplexEigs(ds,0,ds->l,wr,wi));
#if defined(PETSC_USE_COMPLEX)
  if (wi) {
    for (i=ds->l;i<ds->n;i++) wi[i] = 0.0;
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSGetTruncateSize_GHIEP(DS ds,PetscInt l,PetscInt n,PetscInt *k)
{
  PetscReal *T;

  PetscFunctionBegin;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  if (T[l+(*k)-1+ds->ld] !=0.0) {
    if (l+(*k)<n-1) (*k)++;
    else (*k)--;
  }
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSTruncate_GHIEP(DS ds,PetscInt n,PetscBool trim)
{
  PetscInt    i,ld=ds->ld,l=ds->l;
  PetscScalar *A;
  PetscReal   *T,*b,*r,*omega;

  PetscFunctionBegin;
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSGetArrayReal(ds,DS_MAT_D,&omega));
#if defined(PETSC_USE_DEBUG)
    /* make sure diagonal 2x2 block is not broken */
    PetscCheck(ds->state<DS_STATE_CONDENSED || n==0 || n==ds->n || T[n-1+ld]==0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"The given size would break a 2x2 block, call DSGetTruncateSize() first");
#endif
  }
  if (trim) {
    if (!ds->compact && ds->extrarow) {   /* clean extra row */
      PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
      for (i=l;i<ds->n;i++) A[ds->n+i*ld] = 0.0;
      PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
    }
    ds->l = 0;
    ds->k = 0;
    ds->n = n;
    ds->t = ds->n;   /* truncated length equal to the new dimension */
  } else {
    if (!ds->compact && ds->extrarow && ds->k==ds->n) {
      /* copy entries of extra row to the new position, then clean last row */
      PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
      for (i=l;i<n;i++) A[n+i*ld] = A[ds->n+i*ld];
      for (i=l;i<ds->n;i++) A[ds->n+i*ld] = 0.0;
      PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
    }
    if (ds->compact) {
      b = T+ld;
      r = T+2*ld;
      b[n-1] = r[n-1];
      b[n] = b[ds->n];
      omega[n] = omega[ds->n];
    }
    ds->k = (ds->extrarow)? n: 0;
    ds->t = ds->n;   /* truncated length equal to previous dimension */
    ds->n = n;
  }
  if (ds->compact) {
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&omega));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_GHIEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscScalar    *A,*B,*Q;
  PetscReal      *T,*D;
  PetscInt       ld=ds->ld,l=ds->l,k=0,kr=0;
  PetscMPIInt    n,rank,off=0,size,ldn,ld3,ld_;

  PetscFunctionBegin;
  if (ds->compact) kr = 4*ld;
  else k = 2*(ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += (ds->n-l)*ld;
  if (eigr) k += (ds->n-l);
  if (eigi) k += (ds->n-l);
  PetscCall(DSAllocateWork_Private(ds,k+kr,0,0));
  PetscCall(PetscMPIIntCast(k*sizeof(PetscScalar)+kr*sizeof(PetscReal),&size));
  PetscCall(PetscMPIIntCast(ds->n-l,&n));
  PetscCall(PetscMPIIntCast(ld*(ds->n-l),&ldn));
  PetscCall(PetscMPIIntCast(ld*3,&ld3));
  PetscCall(PetscMPIIntCast(ld,&ld_));
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSGetArrayReal(ds,DS_MAT_D,&D));
  } else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  }
  if (ds->state>DS_STATE_RAW) PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank));
  if (!rank) {
    if (ds->compact) {
      PetscCallMPI(MPI_Pack(T,ld3,MPIU_REAL,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Pack(D,ld_,MPIU_REAL,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    } else {
      PetscCallMPI(MPI_Pack(A+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Pack(B+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    }
    if (ds->state>DS_STATE_RAW) PetscCallMPI(MPI_Pack(Q+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (eigr) PetscCallMPI(MPI_Pack(eigr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
#if !defined(PETSC_USE_COMPLEX)
    if (eigi) PetscCallMPI(MPI_Pack(eigi+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
#endif
  }
  PetscCallMPI(MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds)));
  if (rank) {
    if (ds->compact) {
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,T,ld3,MPIU_REAL,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,D,ld_,MPIU_REAL,PetscObjectComm((PetscObject)ds)));
    } else {
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,A+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,B+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    }
    if (ds->state>DS_STATE_RAW) PetscCallMPI(MPI_Unpack(ds->work,size,&off,Q+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (eigr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
#if !defined(PETSC_USE_COMPLEX)
    if (eigi) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigi+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
#endif
  }
  if (ds->compact) {
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&D));
  } else {
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  }
  if (ds->state>DS_STATE_RAW) PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSHermitian_GHIEP(DS ds,DSMatType m,PetscBool *flg)
{
  PetscFunctionBegin;
  if ((m==DS_MAT_A && !ds->extrarow) || m==DS_MAT_B) *flg = PETSC_TRUE;
  else *flg = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSGHIEP - Dense Generalized Hermitian Indefinite Eigenvalue Problem.

   Level: beginner

   Notes:
   The problem is expressed as A*X = B*X*Lambda, where both A and B are
   real symmetric (or complex Hermitian) and possibly indefinite. Lambda
   is a diagonal matrix whose diagonal elements are the arguments of DSSolve().
   After solve, A is overwritten with Lambda. Note that in the case of real
   scalars, A is overwritten with a real representation of Lambda, i.e.,
   complex conjugate eigenvalue pairs are stored as a 2x2 block in the
   quasi-diagonal matrix.

   In the intermediate state A is reduced to tridiagonal form and B is
   transformed into a signature matrix. In compact storage format, these
   matrices are stored in T and D, respectively.

   Used DS matrices:
+  DS_MAT_A - first problem matrix
.  DS_MAT_B - second problem matrix
.  DS_MAT_T - symmetric tridiagonal matrix of the reduced pencil
.  DS_MAT_D - diagonal matrix (signature) of the reduced pencil
-  DS_MAT_Q - pseudo-orthogonal transformation that reduces (A,B) to
   tridiagonal-diagonal form (intermediate step) or a real basis of eigenvectors

   Implemented methods:
+  0 - QR iteration plus inverse iteration for the eigenvectors
.  1 - HZ iteration
-  2 - QR iteration plus pseudo-orthogonalization for the eigenvectors

   References:
.  1. - C. Campos and J. E. Roman, "Restarted Q-Arnoldi-type methods exploiting
   symmetry in quadratic eigenvalue problems", BIT Numer. Math. 56(4):1213-1236, 2016.

.seealso: DSCreate(), DSSetType(), DSType
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_GHIEP(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate        = DSAllocate_GHIEP;
  ds->ops->view            = DSView_GHIEP;
  ds->ops->vectors         = DSVectors_GHIEP;
  ds->ops->solve[0]        = DSSolve_GHIEP_QR_II;
  ds->ops->solve[1]        = DSSolve_GHIEP_HZ;
  ds->ops->solve[2]        = DSSolve_GHIEP_QR;
  ds->ops->sort            = DSSort_GHIEP;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize     = DSSynchronize_GHIEP;
#endif
  ds->ops->gettruncatesize = DSGetTruncateSize_GHIEP;
  ds->ops->truncate        = DSTruncate_GHIEP;
  ds->ops->update          = DSUpdateExtraRow_GHIEP;
  ds->ops->hermitian       = DSHermitian_GHIEP;
  PetscFunctionReturn(PETSC_SUCCESS);
}
