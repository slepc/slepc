/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/dsimpl.h>       /*I "slepcds.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  PetscInt  d;              /* polynomial degree */
  PetscReal *pbc;           /* polynomial basis coefficients */
} DS_PEP;

static PetscErrorCode DSAllocate_PEP(DS ds,PetscInt ld)
{
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCheck(ctx->d,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"DSPEP requires specifying the polynomial degree via DSPEPSetDegree()");
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_X));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Y));
  for (i=0;i<=ctx->d;i++) PetscCall(DSAllocateMat_Private(ds,DSMatExtra[i]));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld*ctx->d,&ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSView_PEP(DS ds,PetscViewer viewer)
{
  DS_PEP            *ctx = (DS_PEP*)ds->data;
  PetscViewerFormat format;
  PetscInt          i;

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO) PetscFunctionReturn(PETSC_SUCCESS);
  if (format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"polynomial degree: %" PetscInt_FMT "\n",ctx->d));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  for (i=0;i<=ctx->d;i++) PetscCall(DSViewMat(ds,viewer,DSMatExtra[i]));
  if (ds->state>DS_STATE_INTERMEDIATE) PetscCall(DSViewMat(ds,viewer,DS_MAT_X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_PEP(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscFunctionBegin;
  PetscCheck(!rnorm,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
  switch (mat) {
    case DS_MAT_X:
      break;
    case DS_MAT_Y:
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_PEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *kout)
{
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       n,i,*perm,told;
  PetscScalar    *A;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(PETSC_SUCCESS);
  n = ds->n*ctx->d;
  perm = ds->perm;
  for (i=0;i<n;i++) perm[i] = i;
  told = ds->t;
  ds->t = n;  /* force the sorting routines to consider d*n eigenvalues */
  if (rr) PetscCall(DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_FALSE));
  else PetscCall(DSSortEigenvalues_Private(ds,wr,wi,perm,PETSC_FALSE));
  ds->t = told;  /* restore value of t */
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  for (i=0;i<n;i++) A[i]  = wr[perm[i]];
  for (i=0;i<n;i++) wr[i] = A[i];
  for (i=0;i<n;i++) A[i]  = wi[perm[i]];
  for (i=0;i<n;i++) wi[i] = A[i];
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(DSPermuteColumnsTwo_Private(ds,0,n,ds->n,DS_MAT_X,DS_MAT_Y,perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_PEP_QZ(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  DS_PEP            *ctx = (DS_PEP*)ds->data;
  PetscInt          i,j,k,off;
  PetscScalar       *A,*B,*W,*X,*U,*Y,*work,*beta,a;
  const PetscScalar *Ed,*Ei;
  PetscReal         *ca,*cb,*cg,norm,done=1.0;
  PetscBLASInt      info,n,ld,ldd,nd,lwork,one=1,zero=0,cols;
  PetscBool         useggev3=(ds->method==1)?PETSC_TRUE:PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n*ctx->d,&nd));
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(ds->ld*ctx->d,&ldd));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_B));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_W));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_U));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));

  /* build matrices A and B of the linearization */
  PetscCall(MatDenseGetArrayRead(ds->omat[DSMatExtra[ctx->d]],&Ed));
  PetscCall(PetscArrayzero(A,ldd*ldd));
  if (!ctx->pbc) { /* monomial basis */
    for (i=0;i<nd-ds->n;i++) A[i+(i+ds->n)*ldd] = 1.0;
    for (i=0;i<ctx->d;i++) {
      PetscCall(MatDenseGetArrayRead(ds->omat[DSMatExtra[i]],&Ei));
      off = i*ds->n*ldd+(ctx->d-1)*ds->n;
      for (j=0;j<ds->n;j++) PetscCall(PetscArraycpy(A+off+j*ldd,Ei+j*ds->ld,ds->n));
      PetscCall(MatDenseRestoreArrayRead(ds->omat[DSMatExtra[i]],&Ei));
    }
  } else {
    ca = ctx->pbc;
    cb = ca+ctx->d+1;
    cg = cb+ctx->d+1;
    for (i=0;i<ds->n;i++) {
      A[i+(i+ds->n)*ldd] = ca[0];
      A[i+i*ldd] = cb[0];
    }
    for (;i<nd-ds->n;i++) {
      j = i/ds->n;
      A[i+(i+ds->n)*ldd] = ca[j];
      A[i+i*ldd] = cb[j];
      A[i+(i-ds->n)*ldd] = cg[j];
    }
    for (i=0;i<ctx->d-2;i++) {
      PetscCall(MatDenseGetArrayRead(ds->omat[DSMatExtra[i]],&Ei));
      off = i*ds->n*ldd+(ctx->d-1)*ds->n;
      for (j=0;j<ds->n;j++)
        for (k=0;k<ds->n;k++)
          A[off+j*ldd+k] = Ei[j*ds->ld+k]*ca[ctx->d-1];
      PetscCall(MatDenseRestoreArrayRead(ds->omat[DSMatExtra[i]],&Ei));
    }
    PetscCall(MatDenseGetArrayRead(ds->omat[DSMatExtra[i]],&Ei));
    off = i*ds->n*ldd+(ctx->d-1)*ds->n;
    for (j=0;j<ds->n;j++)
      for (k=0;k<ds->n;k++)
        A[off+j*ldd+k] = Ei[j*ds->ld+k]*ca[ctx->d-1]-Ed[j*ds->ld+k]*cg[ctx->d-1];
    PetscCall(MatDenseRestoreArrayRead(ds->omat[DSMatExtra[i]],&Ei));
    i++;
    PetscCall(MatDenseGetArrayRead(ds->omat[DSMatExtra[i]],&Ei));
    off = i*ds->n*ldd+(ctx->d-1)*ds->n;
    for (j=0;j<ds->n;j++)
      for (k=0;k<ds->n;k++)
        A[off+j*ldd+k] = Ei[j*ds->ld+k]*ca[ctx->d-1]-Ed[j*ds->ld+k]*cb[ctx->d-1];
    PetscCall(MatDenseRestoreArrayRead(ds->omat[DSMatExtra[i]],&Ei));
  }
  PetscCall(PetscArrayzero(B,ldd*ldd));
  for (i=0;i<nd-ds->n;i++) B[i+i*ldd] = 1.0;
  off = (ctx->d-1)*ds->n*(ldd+1);
  for (j=0;j<ds->n;j++) {
    for (i=0;i<ds->n;i++) B[off+i+j*ldd] = -Ed[i+j*ds->ld];
  }
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DSMatExtra[ctx->d]],&Ed));

  /* solve generalized eigenproblem */
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_W],&W));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_U],&U));
  lwork = -1;
#if defined(PETSC_USE_COMPLEX)
  if (useggev3) PetscCallBLAS("LAPACKggev3",LAPACKggev3_("V","V",&nd,A,&ldd,B,&ldd,wr,NULL,U,&ldd,W,&ldd,&a,&lwork,NULL,&info));
  else PetscCallBLAS("LAPACKggev",LAPACKggev_("V","V",&nd,A,&ldd,B,&ldd,wr,NULL,U,&ldd,W,&ldd,&a,&lwork,NULL,&info));
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork));
  PetscCall(DSAllocateWork_Private(ds,lwork+nd,8*nd,0));
  beta  = ds->work;
  work  = ds->work + nd;
  if (useggev3) PetscCallBLAS("LAPACKggev3",LAPACKggev3_("V","V",&nd,A,&ldd,B,&ldd,wr,beta,U,&ldd,W,&ldd,work,&lwork,ds->rwork,&info));
  else PetscCallBLAS("LAPACKggev",LAPACKggev_("V","V",&nd,A,&ldd,B,&ldd,wr,beta,U,&ldd,W,&ldd,work,&lwork,ds->rwork,&info));
#else
  if (useggev3) PetscCallBLAS("LAPACKggev3",LAPACKggev3_("V","V",&nd,A,&ldd,B,&ldd,wr,wi,NULL,U,&ldd,W,&ldd,&a,&lwork,&info));
  else PetscCallBLAS("LAPACKggev",LAPACKggev_("V","V",&nd,A,&ldd,B,&ldd,wr,wi,NULL,U,&ldd,W,&ldd,&a,&lwork,&info));
  PetscCall(PetscBLASIntCast((PetscInt)a,&lwork));
  PetscCall(DSAllocateWork_Private(ds,lwork+nd,0,0));
  beta = ds->work;
  work = ds->work + nd;
  if (useggev3) PetscCallBLAS("LAPACKggev3",LAPACKggev3_("V","V",&nd,A,&ldd,B,&ldd,wr,wi,beta,U,&ldd,W,&ldd,work,&lwork,&info));
  else PetscCallBLAS("LAPACKggev",LAPACKggev_("V","V",&nd,A,&ldd,B,&ldd,wr,wi,beta,U,&ldd,W,&ldd,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo(useggev3?"ggev3":"ggev",info);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));

  /* copy eigenvalues */
  for (i=0;i<nd;i++) {
    if (beta[i]==0.0) wr[i] = (PetscRealPart(wr[i])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wr[i] /= beta[i];
#if !defined(PETSC_USE_COMPLEX)
    if (beta[i]==0.0) wi[i] = 0.0;
    else wi[i] /= beta[i];
#else
    if (wi) wi[i] = 0.0;
#endif
  }

  /* copy and normalize eigenvectors */
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Y],&Y));
  for (j=0;j<nd;j++) {
    PetscCall(PetscArraycpy(X+j*ds->ld,W+j*ldd,ds->n));
    PetscCall(PetscArraycpy(Y+j*ds->ld,U+ds->n*(ctx->d-1)+j*ldd,ds->n));
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_W],&W));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_U],&U));
  for (j=0;j<nd;j++) {
    cols = 1;
    norm = BLASnrm2_(&n,X+j*ds->ld,&one);
#if !defined(PETSC_USE_COMPLEX)
    if (wi[j] != 0.0) {
      norm = SlepcAbsEigenvalue(norm,BLASnrm2_(&n,X+(j+1)*ds->ld,&one));
      cols = 2;
    }
#endif
    PetscCallBLAS("LAPACKlascl",LAPACKlascl_("G",&zero,&zero,&norm,&done,&n,&cols,X+j*ds->ld,&ld,&info));
    SlepcCheckLapackInfo("lascl",info);
    norm = BLASnrm2_(&n,Y+j*ds->ld,&one);
#if !defined(PETSC_USE_COMPLEX)
    if (wi[j] != 0.0) norm = SlepcAbsEigenvalue(norm,BLASnrm2_(&n,Y+(j+1)*ds->ld,&one));
#endif
    PetscCallBLAS("LAPACKlascl",LAPACKlascl_("G",&zero,&zero,&norm,&done,&n,&cols,Y+j*ds->ld,&ld,&info));
    SlepcCheckLapackInfo("lascl",info);
#if !defined(PETSC_USE_COMPLEX)
    if (wi[j] != 0.0) j++;
#endif
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Y],&Y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_PEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       ld=ds->ld,k=0;
  PetscMPIInt    ldnd,rank,off=0,size,dn;
  PetscScalar    *X,*Y;

  PetscFunctionBegin;
  if (ds->state>=DS_STATE_CONDENSED) k += 2*ctx->d*ds->n*ld;
  if (eigr) k += ctx->d*ds->n;
  if (eigi) k += ctx->d*ds->n;
  PetscCall(DSAllocateWork_Private(ds,k,0,0));
  PetscCall(PetscMPIIntCast(k*sizeof(PetscScalar),&size));
  PetscCall(PetscMPIIntCast(ds->n*ctx->d*ld,&ldnd));
  PetscCall(PetscMPIIntCast(ctx->d*ds->n,&dn));
  if (ds->state>=DS_STATE_CONDENSED) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Y],&Y));
  }
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank));
  if (!rank) {
    if (ds->state>=DS_STATE_CONDENSED) {
      PetscCallMPI(MPI_Pack(X,ldnd,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Pack(Y,ldnd,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    }
    if (eigr) PetscCallMPI(MPI_Pack(eigr,dn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
#if !defined(PETSC_USE_COMPLEX)
    if (eigi) PetscCallMPI(MPI_Pack(eigi,dn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
#endif
  }
  PetscCallMPI(MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds)));
  if (rank) {
    if (ds->state>=DS_STATE_CONDENSED) {
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,X,ldnd,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,Y,ldnd,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    }
    if (eigr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigr,dn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
#if !defined(PETSC_USE_COMPLEX)
    if (eigi) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigi,dn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
#endif
  }
  if (ds->state>=DS_STATE_CONDENSED) {
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Y],&Y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSPEPSetDegree_PEP(DS ds,PetscInt d)
{
  DS_PEP *ctx = (DS_PEP*)ds->data;

  PetscFunctionBegin;
  PetscCheck(d>=0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The degree must be a non-negative integer");
  PetscCheck(d<DS_NUM_EXTRA,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Only implemented for polynomials of degree at most %d",DS_NUM_EXTRA-1);
  ctx->d = d;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSPEPSetDegree - Sets the polynomial degree for a DSPEP.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  d  - the degree

   Level: intermediate

.seealso: DSPEPGetDegree()
@*/
PetscErrorCode DSPEPSetDegree(DS ds,PetscInt d)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,d,2);
  PetscTryMethod(ds,"DSPEPSetDegree_C",(DS,PetscInt),(ds,d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSPEPGetDegree_PEP(DS ds,PetscInt *d)
{
  DS_PEP *ctx = (DS_PEP*)ds->data;

  PetscFunctionBegin;
  *d = ctx->d;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSPEPGetDegree - Returns the polynomial degree for a DSPEP.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  d - the degree

   Level: intermediate

.seealso: DSPEPSetDegree()
@*/
PetscErrorCode DSPEPGetDegree(DS ds,PetscInt *d)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(d,2);
  PetscUseMethod(ds,"DSPEPGetDegree_C",(DS,PetscInt*),(ds,d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSPEPSetCoefficients_PEP(DS ds,PetscReal *pbc)
{
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCheck(ctx->d,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"Must first specify the polynomial degree via DSPEPSetDegree()");
  if (ctx->pbc) PetscCall(PetscFree(ctx->pbc));
  PetscCall(PetscMalloc1(3*(ctx->d+1),&ctx->pbc));
  for (i=0;i<3*(ctx->d+1);i++) ctx->pbc[i] = pbc[i];
  ds->state = DS_STATE_RAW;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSPEPSetCoefficients - Sets the polynomial basis coefficients for a DSPEP.

   Logically Collective

   Input Parameters:
+  ds  - the direct solver context
-  pbc - the polynomial basis coefficients

   Notes:
   This function is required only in the case of a polynomial specified in a
   non-monomial basis, to provide the coefficients that will be used
   during the linearization, multiplying the identity blocks on the three main
   diagonal blocks. Depending on the polynomial basis (Chebyshev, Legendre, ...)
   the coefficients must be different.

   There must be a total of 3*(d+1) coefficients, where d is the degree of the
   polynomial. The coefficients are arranged in three groups, alpha, beta, and
   gamma, according to the definition of the three-term recurrence. In the case
   of the monomial basis, alpha=1 and beta=gamma=0, in which case it is not
   necessary to invoke this function.

   Level: advanced

.seealso: DSPEPGetCoefficients(), DSPEPSetDegree()
@*/
PetscErrorCode DSPEPSetCoefficients(DS ds,PetscReal pbc[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(pbc,2);
  PetscTryMethod(ds,"DSPEPSetCoefficients_C",(DS,PetscReal*),(ds,pbc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSPEPGetCoefficients_PEP(DS ds,PetscReal *pbc[])
{
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCheck(ctx->d,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"Must first specify the polynomial degree via DSPEPSetDegree()");
  PetscCall(PetscCalloc1(3*(ctx->d+1),pbc));
  if (ctx->pbc) for (i=0;i<3*(ctx->d+1);i++) (*pbc)[i] = ctx->pbc[i];
  else for (i=0;i<ctx->d+1;i++) (*pbc)[i] = 1.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSPEPGetCoefficients - Returns the polynomial basis coefficients for a DSPEP.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  pbc - the polynomial basis coefficients

   Note:
   The returned array has length 3*(d+1) and should be freed by the user.

   Fortran Notes:
   The calling sequence from Fortran is
.vb
   DSPEPGetCoefficients(eps,pbc,ierr)
   double precision pbc(d+1) output
.ve

   Level: advanced

.seealso: DSPEPSetCoefficients()
@*/
PetscErrorCode DSPEPGetCoefficients(DS ds,PetscReal *pbc[]) PeNS
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(pbc,2);
  PetscUseMethod(ds,"DSPEPGetCoefficients_C",(DS,PetscReal**),(ds,pbc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSDestroy_PEP(DS ds)
{
  DS_PEP         *ctx = (DS_PEP*)ds->data;

  PetscFunctionBegin;
  if (ctx->pbc) PetscCall(PetscFree(ctx->pbc));
  PetscCall(PetscFree(ds->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSPEPSetDegree_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSPEPGetDegree_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSPEPSetCoefficients_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSPEPGetCoefficients_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSMatGetSize_PEP(DS ds,DSMatType t,PetscInt *rows,PetscInt *cols)
{
  DS_PEP *ctx = (DS_PEP*)ds->data;

  PetscFunctionBegin;
  PetscCheck(ctx->d,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"DSPEP requires specifying the polynomial degree via DSPEPSetDegree()");
  *rows = ds->n;
  if (t==DS_MAT_A || t==DS_MAT_B || t==DS_MAT_W || t==DS_MAT_U) *rows *= ctx->d;
  *cols = ds->n;
  if (t==DS_MAT_A || t==DS_MAT_B || t==DS_MAT_W || t==DS_MAT_U || t==DS_MAT_X || t==DS_MAT_Y) *cols *= ctx->d;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSPEP - Dense Polynomial Eigenvalue Problem.

   Level: beginner

   Notes:
   The problem is expressed as P(lambda)*x = 0, where P(.) is a matrix
   polynomial of degree d. The eigenvalues lambda are the arguments
   returned by DSSolve().

   The degree of the polynomial, d, can be set with DSPEPSetDegree(), with
   the first d+1 extra matrices of the DS defining the matrix polynomial. By
   default, the polynomial is expressed in the monomial basis, but a
   different basis can be used by setting the corresponding coefficients
   via DSPEPSetCoefficients().

   The problem is solved via linearization, by building a pencil (A,B) of
   size p*n and solving the corresponding GNHEP.

   Used DS matrices:
+  DS_MAT_Ex - coefficients of the matrix polynomial
.  DS_MAT_X  - right eigenvectors
.  DS_MAT_Y  - left eigenvectors
.  DS_MAT_A  - (workspace) first matrix of the linearization
.  DS_MAT_B  - (workspace) second matrix of the linearization
.  DS_MAT_W  - (workspace) right eigenvectors of the linearization
-  DS_MAT_U  - (workspace) left eigenvectors of the linearization

   Implemented methods:
.  0 - QZ iteration on the linearization (_ggev)

.seealso: DSCreate(), DSSetType(), DSType, DSPEPSetDegree(), DSPEPSetCoefficients()
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_PEP(DS ds)
{
  DS_PEP         *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  ds->data = (void*)ctx;

  ds->ops->allocate      = DSAllocate_PEP;
  ds->ops->view          = DSView_PEP;
  ds->ops->vectors       = DSVectors_PEP;
  ds->ops->solve[0]      = DSSolve_PEP_QZ;
#if !defined(SLEPC_MISSING_LAPACK_GGES3)
  ds->ops->solve[1]      = DSSolve_PEP_QZ;
#endif
  ds->ops->sort          = DSSort_PEP;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize   = DSSynchronize_PEP;
#endif
  ds->ops->destroy       = DSDestroy_PEP;
  ds->ops->matgetsize    = DSMatGetSize_PEP;
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSPEPSetDegree_C",DSPEPSetDegree_PEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSPEPGetDegree_C",DSPEPGetDegree_PEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSPEPSetCoefficients_C",DSPEPSetCoefficients_PEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSPEPGetCoefficients_C",DSPEPGetCoefficients_PEP));
  PetscFunctionReturn(PETSC_SUCCESS);
}
