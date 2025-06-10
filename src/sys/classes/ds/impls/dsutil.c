/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Utility subroutines common to several impls
*/

#include <slepc/private/dsimpl.h>      /*I "slepcds.h" I*/
#include <slepcblaslapack.h>

/*
   Compute the (real) Schur form of A. At the end, A is (quasi-)triangular and Q
   contains the unitary matrix of Schur vectors. Eigenvalues are returned in wr,wi
*/
PetscErrorCode DSSolve_NHEP_Private(DS ds,DSMatType mA,DSMatType mQ,PetscScalar *wr,PetscScalar *wi)
{
  PetscScalar    *work,*tau,*A,*Q;
  PetscInt       i,j;
  PetscBLASInt   ilo,lwork,info,n,k,ld;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ds->omat[mA],&A));
  PetscCall(MatDenseGetArray(ds->omat[mQ],&Q));
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(ds->l+1,&ilo));
  PetscCall(PetscBLASIntCast(ds->k,&k));
  PetscCall(DSAllocateWork_Private(ds,ld+6*ld,0,0));
  tau  = ds->work;
  work = ds->work+ld;
  lwork = 6*ld;

  /* initialize orthogonal matrix */
  PetscCall(PetscArrayzero(Q,ld*ld));
  for (i=0;i<n;i++) Q[i+i*ld] = 1.0;
  if (n==1) { /* quick return */
    wr[0] = A[0];
    if (wi) wi[0] = 0.0;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* reduce to upper Hessenberg form */
  if (ds->state<DS_STATE_INTERMEDIATE) {
    PetscCallBLAS("LAPACKgehrd",LAPACKgehrd_(&n,&ilo,&n,A,&ld,tau,work,&lwork,&info));
    SlepcCheckLapackInfo("gehrd",info);
    for (j=0;j<n-1;j++) {
      for (i=j+2;i<n;i++) {
        Q[i+j*ld] = A[i+j*ld];
        A[i+j*ld] = 0.0;
      }
    }
    PetscCallBLAS("LAPACKorghr",LAPACKorghr_(&n,&ilo,&n,Q,&ld,tau,work,&lwork,&info));
    SlepcCheckLapackInfo("orghr",info);
  }

  /* compute the (real) Schur form */
#if !defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKhseqr",LAPACKhseqr_("S","V",&n,&ilo,&n,A,&ld,wr,wi,Q,&ld,work,&lwork,&info));
  for (j=0;j<ds->l;j++) {
    if (j==n-1 || A[j+1+j*ld] == 0.0) {
      /* real eigenvalue */
      wr[j] = A[j+j*ld];
      wi[j] = 0.0;
    } else {
      /* complex eigenvalue */
      wr[j] = A[j+j*ld];
      wr[j+1] = A[j+j*ld];
      wi[j] = PetscSqrtReal(PetscAbsReal(A[j+1+j*ld]))*PetscSqrtReal(PetscAbsReal(A[j+(j+1)*ld]));
      wi[j+1] = -wi[j];
      j++;
    }
  }
#else
  PetscCallBLAS("LAPACKhseqr",LAPACKhseqr_("S","V",&n,&ilo,&n,A,&ld,wr,Q,&ld,work,&lwork,&info));
  if (wi) for (i=ds->l;i<n;i++) wi[i] = 0.0;
#endif
  SlepcCheckLapackInfo("hseqr",info);
  PetscCall(MatDenseRestoreArray(ds->omat[mA],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[mQ],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Sort a Schur form represented by the (quasi-)triangular matrix T and
   the unitary matrix Q, and return the sorted eigenvalues in wr,wi
*/
PetscErrorCode DSSort_NHEP_Total(DS ds,DSMatType mT,DSMatType mQ,PetscScalar *wr,PetscScalar *wi)
{
  PetscScalar    re,*T,*Q;
  PetscInt       i,j,pos,result;
  PetscBLASInt   ifst,ilst,info,n,ld;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    *work,im;
#endif

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ds->omat[mT],&T));
  PetscCall(MatDenseGetArray(ds->omat[mQ],&Q));
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(DSAllocateWork_Private(ds,ld,0,0));
  work = ds->work;
#endif
  /* selection sort */
  for (i=ds->l;i<n-1;i++) {
    re = wr[i];
#if !defined(PETSC_USE_COMPLEX)
    im = wi[i];
#endif
    pos = 0;
    j=i+1; /* j points to the next eigenvalue */
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
      PetscCallBLAS("LAPACKtrexc",LAPACKtrexc_("V",&n,T,&ld,Q,&ld,&ifst,&ilst,work,&info));
#else
      PetscCallBLAS("LAPACKtrexc",LAPACKtrexc_("V",&n,T,&ld,Q,&ld,&ifst,&ilst,&info));
#endif
      SlepcCheckLapackInfo("trexc",info);
      /* recover original eigenvalues from T matrix */
      for (j=i;j<n;j++) {
        wr[j] = T[j+j*ld];
#if !defined(PETSC_USE_COMPLEX)
        if (j<n-1 && T[j+1+j*ld] != 0.0) {
          /* complex conjugate eigenvalue */
          wi[j] = PetscSqrtReal(PetscAbsReal(T[j+1+j*ld]))*PetscSqrtReal(PetscAbsReal(T[j+(j+1)*ld]));
          wr[j+1] = wr[j];
          wi[j+1] = -wi[j];
          j++;
        } else wi[j] = 0.0;
#endif
      }
    }
#if !defined(PETSC_USE_COMPLEX)
    if (wi[i] != 0) i++;
#endif
  }
  PetscCall(MatDenseRestoreArray(ds->omat[mT],&T));
  PetscCall(MatDenseRestoreArray(ds->omat[mQ],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Reorder a Schur form represented by T,Q according to a permutation perm,
   and return the sorted eigenvalues in wr,wi
*/
PetscErrorCode DSSortWithPermutation_NHEP_Private(DS ds,PetscInt *perm,DSMatType mT,DSMatType mQ,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt       i,j,pos,inc=1;
  PetscBLASInt   ifst,ilst,info,n,ld;
  PetscScalar    *T,*Q;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    *work;
#endif

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ds->omat[mT],&T));
  PetscCall(MatDenseGetArray(ds->omat[mQ],&Q));
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(DSAllocateWork_Private(ds,ld,0,0));
  work = ds->work;
#endif
  for (i=ds->l;i<n-1;i++) {
    pos = perm[i];
#if !defined(PETSC_USE_COMPLEX)
    inc = (pos<n-1 && T[pos+1+pos*ld] != 0.0)? 2: 1;
#endif
    if (pos!=i) {
#if !defined(PETSC_USE_COMPLEX)
      PetscCheck((T[pos+(pos-1)*ld]==0.0 || perm[i+1]==pos-1) && (pos==n-1 || (T[pos+1+pos*ld]==0.0 || perm[i+1]==pos+1)),PETSC_COMM_SELF,PETSC_ERR_FP,"Invalid permutation due to a 2x2 block at position %" PetscInt_FMT,pos);
#endif
      /* interchange blocks */
      PetscCall(PetscBLASIntCast(pos+1,&ifst));
      PetscCall(PetscBLASIntCast(i+1,&ilst));
#if !defined(PETSC_USE_COMPLEX)
      PetscCallBLAS("LAPACKtrexc",LAPACKtrexc_("V",&n,T,&ld,Q,&ld,&ifst,&ilst,work,&info));
#else
      PetscCallBLAS("LAPACKtrexc",LAPACKtrexc_("V",&n,T,&ld,Q,&ld,&ifst,&ilst,&info));
#endif
      SlepcCheckLapackInfo("trexc",info);
      for (j=i+1;j<n;j++) {
        if (perm[j]>=i && perm[j]<pos) perm[j]+=inc;
      }
      perm[i] = i;
      if (inc==2) perm[i+1] = i+1;
    }
    if (inc==2) i++;
  }
  /* recover original eigenvalues from T matrix */
  for (j=ds->l;j<n;j++) {
    wr[j] = T[j+j*ld];
#if !defined(PETSC_USE_COMPLEX)
    if (j<n-1 && T[j+1+j*ld] != 0.0) {
      /* complex conjugate eigenvalue */
      wi[j] = PetscSqrtReal(PetscAbsReal(T[j+1+j*ld]))*PetscSqrtReal(PetscAbsReal(T[j+(j+1)*ld]));
      wr[j+1] = wr[j];
      wi[j+1] = -wi[j];
      j++;
    } else wi[j] = 0.0;
#endif
  }
  PetscCall(MatDenseRestoreArray(ds->omat[mT],&T));
  PetscCall(MatDenseRestoreArray(ds->omat[mQ],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}
