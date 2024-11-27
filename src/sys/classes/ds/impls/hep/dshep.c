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

static PetscErrorCode DSAllocate_HEP(DS ds,PetscInt ld)
{
  PetscFunctionBegin;
  if (!ds->compact) PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Q));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_T));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*   0       l           k                 n-1
    -----------------------------------------
    |*       .           .                  |
    |  *     .           .                  |
    |    *   .           .                  |
    |      * .           .                  |
    |. . . . o           o                  |
    |          o         o                  |
    |            o       o                  |
    |              o     o                  |
    |                o   o                  |
    |                  o o                  |
    |. . . . o o o o o o o x                |
    |                    x x x              |
    |                      x x x            |
    |                        x x x          |
    |                          x x x        |
    |                            x x x      |
    |                              x x x    |
    |                                x x x  |
    |                                  x x x|
    |                                    x x|
    -----------------------------------------
*/

static PetscErrorCode DSSwitchFormat_HEP(DS ds)
{
  PetscReal      *T;
  PetscScalar    *A;
  PetscInt       i,n=ds->n,k=ds->k,ld=ds->ld;

  PetscFunctionBegin;
  /* switch from compact (arrow) to dense storage */
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_A],&A));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  PetscCall(PetscArrayzero(A,ld*ld));
  for (i=0;i<k;i++) {
    A[i+i*ld] = T[i];
    A[k+i*ld] = T[i+ld];
    A[i+k*ld] = T[i+ld];
  }
  A[k+k*ld] = T[k];
  for (i=k+1;i<n;i++) {
    A[i+i*ld]     = T[i];
    A[i-1+i*ld]   = T[i-1+ld];
    A[i+(i-1)*ld] = T[i-1+ld];
  }
  if (ds->extrarow) A[n+(n-1)*ld] = T[n-1+ld];
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_A],&A));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSView_HEP(DS ds,PetscViewer viewer)
{
  PetscViewerFormat format;
  PetscInt          i,j,r,c,rows;
  PetscReal         *T,value;
  const char        *methodname[] = {
                     "Implicit QR method (_steqr)",
                     "Relatively Robust Representations (_stevr)",
                     "Divide and Conquer method (_stedc)",
                     "Block Divide and Conquer method (dsbtdc)"
  };
  const int         nmeth=PETSC_STATIC_ARRAY_LENGTH(methodname);

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
    if (ds->bs>1) PetscCall(PetscViewerASCIIPrintf(viewer,"block size: %" PetscInt_FMT "\n",ds->bs));
    if (ds->method<nmeth) PetscCall(PetscViewerASCIIPrintf(viewer,"solving the problem with: %s\n",methodname[ds->method]));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    rows = ds->extrarow? ds->n+1: ds->n;
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%% Size = %" PetscInt_FMT " %" PetscInt_FMT "\n",rows,ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = zeros(%" PetscInt_FMT ",3);\n",3*ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = [\n"));
      for (i=0;i<ds->n;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,i+1,(double)T[i]));
      for (i=0;i<rows-1;i++) {
        r = PetscMax(i+2,ds->k+1);
        c = i+1;
        PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",r,c,(double)T[i+ds->ld]));
        if (i<ds->n-1 && ds->k<ds->n) { /* do not print vertical arrow when k=n */
          PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",c,r,(double)T[i+ds->ld]));
        }
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(zzz);\n",DSMatName[DS_MAT_T]));
    } else {
      for (i=0;i<rows;i++) {
        for (j=0;j<ds->n;j++) {
          if (i==j) value = T[i];
          else if ((i<ds->k && j==ds->k) || (i==ds->k && j<ds->k)) value = T[PetscMin(i,j)+ds->ld];
          else if (i==j+1 && i>ds->k) value = T[i-1+ds->ld];
          else if (i+1==j && j>ds->k) value = T[j-1+ds->ld];
          else value = 0.0;
          PetscCall(PetscViewerASCIIPrintf(viewer," %18.16e ",(double)value));
        }
        PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      }
    }
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerFlush(viewer));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  } else PetscCall(DSViewMat(ds,viewer,DS_MAT_A));
  if (ds->state>DS_STATE_INTERMEDIATE) PetscCall(DSViewMat(ds,viewer,DS_MAT_Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_HEP(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscScalar       *Z;
  const PetscScalar *Q;
  PetscInt          ld = ds->ld;

  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_X:
    case DS_MAT_Y:
      if (j) {
        PetscCall(MatDenseGetArray(ds->omat[mat],&Z));
        if (ds->state>=DS_STATE_CONDENSED) {
          PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Q],&Q));
          PetscCall(PetscArraycpy(Z+(*j)*ld,Q+(*j)*ld,ld));
          if (rnorm) *rnorm = PetscAbsScalar(Q[ds->n-1+(*j)*ld]);
          PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_Q],&Q));
        } else {
          PetscCall(PetscArrayzero(Z+(*j)*ld,ld));
          Z[(*j)+(*j)*ld] = 1.0;
          if (rnorm) *rnorm = 0.0;
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

/*
  ARROWTRIDIAG reduces a symmetric arrowhead matrix of the form

                [ d 0 0 0 e ]
                [ 0 d 0 0 e ]
            A = [ 0 0 d 0 e ]
                [ 0 0 0 d e ]
                [ e e e e d ]

  to tridiagonal form

                [ d e 0 0 0 ]
                [ e d e 0 0 ]
   T = Q'*A*Q = [ 0 e d e 0 ],
                [ 0 0 e d e ]
                [ 0 0 0 e d ]

  where Q is an orthogonal matrix. Rutishauser's algorithm is used to
  perform the reduction, which requires O(n**2) flops. The accumulation
  of the orthogonal factor Q, however, requires O(n**3) flops.

  Arguments
  =========

  N       (input) INTEGER
          The order of the matrix A.  N >= 0.

  D       (input/output) DOUBLE PRECISION array, dimension (N)
          On entry, the diagonal entries of the matrix A to be
          reduced.
          On exit, the diagonal entries of the reduced matrix T.

  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
          On entry, the off-diagonal entries of the matrix A to be
          reduced.
          On exit, the subdiagonal entries of the reduced matrix T.

  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
          On exit, the orthogonal matrix Q.

  LDQ     (input) INTEGER
          The leading dimension of the array Q.

  Note
  ====
  Based on Fortran code contributed by Daniel Kressner
*/
PetscErrorCode DSArrowTridiag(PetscBLASInt n,PetscReal *d,PetscReal *e,PetscScalar *Q,PetscBLASInt ld)
{
  PetscBLASInt i,j,j2,one=1;
  PetscReal    c,s,p,off,temp;

  PetscFunctionBegin;
  if (n<=2) PetscFunctionReturn(PETSC_SUCCESS);

  for (j=0;j<n-2;j++) {

    /* Eliminate entry e(j) by a rotation in the planes (j,j+1) */
    temp = e[j+1];
    PetscCallBLAS("LAPACKlartg",LAPACKREALlartg_(&temp,&e[j],&c,&s,&e[j+1]));
    s = -s;

    /* Apply rotation to diagonal elements */
    temp   = d[j+1];
    e[j]   = c*s*(temp-d[j]);
    d[j+1] = s*s*d[j] + c*c*temp;
    d[j]   = c*c*d[j] + s*s*temp;

    /* Apply rotation to Q */
    j2 = j+2;
    PetscCallBLAS("BLASrot",BLASMIXEDrot_(&j2,Q+j*ld,&one,Q+(j+1)*ld,&one,&c,&s));

    /* Chase newly introduced off-diagonal entry to the top left corner */
    for (i=j-1;i>=0;i--) {
      off  = -s*e[i];
      e[i] = c*e[i];
      temp = e[i+1];
      PetscCallBLAS("LAPACKlartg",LAPACKREALlartg_(&temp,&off,&c,&s,&e[i+1]));
      s = -s;
      temp = (d[i]-d[i+1])*s - 2.0*c*e[i];
      p = s*temp;
      d[i+1] += p;
      d[i] -= p;
      e[i] = -e[i] - c*temp;
      PetscCallBLAS("BLASrot",BLASMIXEDrot_(&j2,Q+i*ld,&one,Q+(i+1)*ld,&one,&c,&s));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Reduce to tridiagonal form by means of DSArrowTridiag.
*/
static PetscErrorCode DSIntermediate_HEP(DS ds)
{
  PetscInt          i;
  PetscBLASInt      n1 = 0,n2,lwork,info,l = 0,n = 0,ld,off;
  PetscScalar       *Q,*work,*tau;
  const PetscScalar *A;
  PetscReal         *d,*e;
  Mat               At,Qt;  /* trailing submatrices */

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->l,&l));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(PetscMax(0,ds->k-l+1),&n1)); /* size of leading block, excl. locked */
  n2 = n-l;     /* n2 = n1 + size of trailing block */
  off = l+l*ld;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d+ld;
  PetscCall(DSSetIdentity(ds,DS_MAT_Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));

  if (ds->compact) {

    if (ds->state<DS_STATE_INTERMEDIATE) PetscCall(DSArrowTridiag(n1,d+l,e+l,Q+off,ld));

  } else {

    PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
    for (i=0;i<l;i++) { d[i] = PetscRealPart(A[i+i*ld]); e[i] = 0.0; }

    if (ds->state<DS_STATE_INTERMEDIATE) {
      PetscCall(MatDenseGetSubMatrix(ds->omat[DS_MAT_A],ds->l,ds->n,ds->l,ds->n,&At));
      PetscCall(MatDenseGetSubMatrix(ds->omat[DS_MAT_Q],ds->l,ds->n,ds->l,ds->n,&Qt));
      PetscCall(MatCopy(At,Qt,SAME_NONZERO_PATTERN));
      PetscCall(MatDenseRestoreSubMatrix(ds->omat[DS_MAT_A],&At));
      PetscCall(MatDenseRestoreSubMatrix(ds->omat[DS_MAT_Q],&Qt));
      PetscCall(DSAllocateWork_Private(ds,ld+ld*ld,0,0));
      tau  = ds->work;
      work = ds->work+ld;
      lwork = ld*ld;
      PetscCallBLAS("LAPACKsytrd",LAPACKsytrd_("L",&n2,Q+off,&ld,d+l,e+l,tau,work,&lwork,&info));
      SlepcCheckLapackInfo("sytrd",info);
      PetscCallBLAS("LAPACKorgtr",LAPACKorgtr_("L",&n2,Q+off,&ld,tau,work,&lwork,&info));
      SlepcCheckLapackInfo("orgtr",info);
    } else {
      /* copy tridiagonal to d,e */
      for (i=l;i<n;i++)   d[i] = PetscRealPart(A[i+i*ld]);
      for (i=l;i<n-1;i++) e[i] = PetscRealPart(A[(i+1)+i*ld]);
    }
    PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_HEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscInt       n,l,i,*perm,ld=ds->ld;
  PetscScalar    *A;
  PetscReal      *d;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(PETSC_SUCCESS);
  n = ds->n;
  l = ds->l;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  perm = ds->perm;
  if (!rr) PetscCall(DSSortEigenvaluesReal_Private(ds,d,perm));
  else PetscCall(DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_FALSE));
  for (i=l;i<n;i++) wr[i] = d[perm[i]];
  PetscCall(DSPermuteColumns_Private(ds,l,n,n,DS_MAT_Q,perm));
  for (i=l;i<n;i++) d[i] = PetscRealPart(wr[i]);
  if (!ds->compact) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    for (i=l;i<n;i++) A[i+i*ld] = wr[i];
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSUpdateExtraRow_HEP(DS ds)
{
  PetscInt          i;
  PetscBLASInt      n,ld,incx=1;
  PetscScalar       *A,*x,*y,one=1.0,zero=0.0;
  PetscReal         *T,*e,beta;
  const PetscScalar *Q;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_Q],&Q));
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    e = T+ld;
    beta = e[n-1];   /* in compact, we assume all entries are zero except the last one */
    for (i=0;i<n;i++) e[i] = PetscRealPart(beta*Q[n-1+i*ld]);
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
    ds->k = n;
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

static PetscErrorCode DSSolve_HEP_QR(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt       i;
  PetscBLASInt   n1,info,l = 0,n = 0,ld,off;
  PetscScalar    *Q,*A;
  PetscReal      *d,*e;

  PetscFunctionBegin;
  PetscCheck(ds->bs==1,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"This method is not prepared for bs>1");
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->l,&l));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  n1 = n-l;     /* n1 = size of leading block, excl. locked + size of trailing block */
  off = l+l*ld;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d+ld;

  /* Reduce to tridiagonal form */
  PetscCall(DSIntermediate_HEP(ds));

  /* Solve the tridiagonal eigenproblem */
  for (i=0;i<l;i++) wr[i] = d[i];

  PetscCall(DSAllocateWork_Private(ds,0,2*ld,0));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCallBLAS("LAPACKsteqr",LAPACKsteqr_("V",&n1,d+l,e+l,Q+off,&ld,ds->rwork,&info));
  SlepcCheckLapackInfo("steqr",info);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  for (i=l;i<n;i++) wr[i] = d[i];

  /* Create diagonal matrix as a result */
  if (ds->compact) PetscCall(PetscArrayzero(e,n-1));
  else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    for (i=l;i<n;i++) PetscCall(PetscArrayzero(A+l+i*ld,n-l));
    for (i=l;i<n;i++) A[i+i*ld] = d[i];
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));

  /* Set zero wi */
  if (wi) for (i=l;i<n;i++) wi[i] = 0.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_HEP_MRRR(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  Mat            At,Qt;  /* trailing submatrices */
  PetscInt       i;
  PetscBLASInt   n1 = 0,n2 = 0,n3,lrwork,liwork,info,l = 0,n = 0,m = 0,ld,off,il,iu,*isuppz;
  PetscScalar    *A,*Q,*W=NULL,one=1.0,zero=0.0;
  PetscReal      *d,*e,abstol=0.0,vl,vu;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
  PetscReal      *Qr,*ritz;
#endif

  PetscFunctionBegin;
  PetscCheck(ds->bs==1,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"This method is not prepared for bs>1");
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->l,&l));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(ds->k-l+1,&n1)); /* size of leading block, excl. locked */
  PetscCall(PetscBLASIntCast(n-ds->k-1,&n2)); /* size of trailing block */
  n3 = n1+n2;
  off = l+l*ld;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d+ld;

  /* Reduce to tridiagonal form */
  PetscCall(DSIntermediate_HEP(ds));

  /* Solve the tridiagonal eigenproblem */
  for (i=0;i<l;i++) wr[i] = d[i];

  if (ds->state<DS_STATE_INTERMEDIATE) {  /* Q contains useful info */
    PetscCall(DSAllocateMat_Private(ds,DS_MAT_W));
    PetscCall(MatCopy(ds->omat[DS_MAT_Q],ds->omat[DS_MAT_W],SAME_NONZERO_PATTERN));
  }
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  lrwork = 20*ld;
  liwork = 10*ld;
#if defined(PETSC_USE_COMPLEX)
  PetscCall(DSAllocateWork_Private(ds,0,lrwork+ld+ld*ld,liwork+2*ld));
#else
  PetscCall(DSAllocateWork_Private(ds,0,lrwork+ld,liwork+2*ld));
#endif
  isuppz = ds->iwork+liwork;
#if defined(PETSC_USE_COMPLEX)
  ritz = ds->rwork+lrwork;
  Qr   = ds->rwork+lrwork+ld;
  PetscCallBLAS("LAPACKstevr",LAPACKstevr_("V","A",&n3,d+l,e+l,&vl,&vu,&il,&iu,&abstol,&m,ritz+l,Qr+off,&ld,isuppz,ds->rwork,&lrwork,ds->iwork,&liwork,&info));
  for (i=l;i<n;i++) wr[i] = ritz[i];
#else
  PetscCallBLAS("LAPACKstevr",LAPACKstevr_("V","A",&n3,d+l,e+l,&vl,&vu,&il,&iu,&abstol,&m,wr+l,Q+off,&ld,isuppz,ds->rwork,&lrwork,ds->iwork,&liwork,&info));
#endif
  SlepcCheckLapackInfo("stevr",info);
#if defined(PETSC_USE_COMPLEX)
  for (i=l;i<n;i++)
    for (j=l;j<n;j++)
      Q[i+j*ld] = Qr[i+j*ld];
#endif
  if (ds->state<DS_STATE_INTERMEDIATE) {  /* accumulate previous Q */
    if (ds->compact) PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_W],&W));
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n3,&n3,&n3,&one,W+off,&ld,Q+off,&ld,&zero,A+off,&ld));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_W],&W));
    PetscCall(MatDenseGetSubMatrix(ds->omat[DS_MAT_A],ds->l,ds->n,ds->l,ds->n,&At));
    PetscCall(MatDenseGetSubMatrix(ds->omat[DS_MAT_Q],ds->l,ds->n,ds->l,ds->n,&Qt));
    PetscCall(MatCopy(At,Qt,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(ds->omat[DS_MAT_A],&At));
    PetscCall(MatDenseRestoreSubMatrix(ds->omat[DS_MAT_Q],&Qt));
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  for (i=l;i<n;i++) d[i] = PetscRealPart(wr[i]);

  /* Create diagonal matrix as a result */
  if (ds->compact) PetscCall(PetscArrayzero(e,n-1));
  else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    for (i=l;i<n;i++) PetscCall(PetscArrayzero(A+l+i*ld,n-l));
    for (i=l;i<n;i++) A[i+i*ld] = d[i];
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));

  /* Set zero wi */
  if (wi) for (i=l;i<n;i++) wi[i] = 0.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_HEP_DC(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt       i;
  PetscBLASInt   n1,info,l = 0,ld,off,lrwork,liwork;
  PetscScalar    *Q,*A;
  PetscReal      *d,*e;
#if defined(PETSC_USE_COMPLEX)
  PetscBLASInt   lwork;
  PetscInt       j;
#endif

  PetscFunctionBegin;
  PetscCheck(ds->bs==1,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"This method is not prepared for bs>1");
  PetscCall(PetscBLASIntCast(ds->l,&l));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(ds->n-ds->l,&n1));
  off = l+l*ld;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d+ld;

  /* Reduce to tridiagonal form */
  PetscCall(DSIntermediate_HEP(ds));

  /* Solve the tridiagonal eigenproblem */
  for (i=0;i<l;i++) wr[i] = d[i];

  lrwork = 5*n1*n1+3*n1+1;
  liwork = 5*n1*n1+6*n1+6;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(DSAllocateWork_Private(ds,0,lrwork,liwork));
  PetscCallBLAS("LAPACKstedc",LAPACKstedc_("V",&n1,d+l,e+l,Q+off,&ld,ds->rwork,&lrwork,ds->iwork,&liwork,&info));
#else
  lwork = ld*ld;
  PetscCall(DSAllocateWork_Private(ds,lwork,lrwork,liwork));
  PetscCallBLAS("LAPACKstedc",LAPACKstedc_("V",&n1,d+l,e+l,Q+off,&ld,ds->work,&lwork,ds->rwork,&lrwork,ds->iwork,&liwork,&info));
  /* Fixing Lapack bug*/
  for (j=ds->l;j<ds->n;j++)
    for (i=0;i<ds->l;i++) Q[i+j*ld] = 0.0;
#endif
  SlepcCheckLapackInfo("stedc",info);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  for (i=l;i<ds->n;i++) wr[i] = d[i];

  /* Create diagonal matrix as a result */
  if (ds->compact) PetscCall(PetscArrayzero(e,ds->n-1));
  else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    for (i=l;i<ds->n;i++) PetscCall(PetscArrayzero(A+l+i*ld,ds->n-l));
    for (i=l;i<ds->n;i++) A[i+i*ld] = d[i];
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));

  /* Set zero wi */
  if (wi) for (i=l;i<ds->n;i++) wi[i] = 0.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_USE_COMPLEX)
static PetscErrorCode DSSolve_HEP_BDC(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscBLASInt   i,j,k,m,n = 0,info,nblks,bs = 0,ld = 0,lde,lrwork,liwork,*ksizes,*iwork,mingapi;
  PetscScalar    *Q,*A;
  PetscReal      *D,*E,*d,*e,tol=PETSC_MACHINE_EPSILON/2,tau1=1e-16,tau2=1e-18,*rwork,mingap;

  PetscFunctionBegin;
  PetscCheck(ds->l==0,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"This method is not prepared for l>1");
  PetscCheck(!ds->compact,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented for compact storage");
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(ds->bs,&bs));
  PetscCall(PetscBLASIntCast(ds->n,&n));
  nblks = n/bs;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d+ld;
  lrwork = 4*n*n+60*n+1;
  liwork = 5*n+5*nblks-1;
  lde = 2*bs+1;
  PetscCall(DSAllocateWork_Private(ds,bs*n+lde*lde*(nblks-1),lrwork,nblks+liwork));
  D      = ds->work;
  E      = ds->work+bs*n;
  rwork  = ds->rwork;
  ksizes = ds->iwork;
  iwork  = ds->iwork+nblks;
  PetscCall(PetscArrayzero(iwork,liwork));

  /* Copy matrix to block tridiagonal format */
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  j=0;
  for (i=0;i<nblks;i++) {
    ksizes[i]=bs;
    for (k=0;k<bs;k++)
      for (m=0;m<bs;m++)
        D[k+m*bs+i*bs*bs] = PetscRealPart(A[j+k+(j+m)*n]);
    j = j + bs;
  }
  j=0;
  for (i=0;i<nblks-1;i++) {
    for (k=0;k<bs;k++)
      for (m=0;m<bs;m++)
        E[k+m*lde+i*lde*lde] = PetscRealPart(A[j+bs+k+(j+m)*n]);
    j = j + bs;
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));

  /* Solve the block tridiagonal eigenproblem */
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(BDC_dsbtdc_("D","A",n,nblks,ksizes,D,bs,bs,E,lde,lde,tol,tau1,tau2,d,Q,n,rwork,lrwork,iwork,liwork,&mingap,&mingapi,&info,1,1));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  for (i=0;i<ds->n;i++) wr[i] = d[i];

  /* Create diagonal matrix as a result */
  if (ds->compact) PetscCall(PetscArrayzero(e,ds->n-1));
  else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    for (i=0;i<ds->n;i++) PetscCall(PetscArrayzero(A+i*ld,ds->n));
    for (i=0;i<ds->n;i++) A[i+i*ld] = wr[i];
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));

  /* Set zero wi */
  if (wi) for (i=0;i<ds->n;i++) wi[i] = 0.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSTruncate_HEP(DS ds,PetscInt n,PetscBool trim)
{
  PetscInt    i,ld=ds->ld,l=ds->l;
  PetscScalar *A;

  PetscFunctionBegin;
  if (!ds->compact && ds->extrarow) PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  if (trim) {
    if (!ds->compact && ds->extrarow) {   /* clean extra row */
      for (i=l;i<ds->n;i++) A[ds->n+i*ld] = 0.0;
    }
    ds->l = 0;
    ds->k = 0;
    ds->n = n;
    ds->t = ds->n;   /* truncated length equal to the new dimension */
  } else {
    if (!ds->compact && ds->extrarow && ds->k==ds->n) {
      /* copy entries of extra row to the new position, then clean last row */
      for (i=l;i<n;i++) A[n+i*ld] = A[ds->n+i*ld];
      for (i=l;i<ds->n;i++) A[ds->n+i*ld] = 0.0;
    }
    ds->k = (ds->extrarow)? n: 0;
    ds->t = ds->n;   /* truncated length equal to previous dimension */
    ds->n = n;
  }
  if (!ds->compact && ds->extrarow) PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_HEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscInt       ld=ds->ld,l=ds->l,k=0,kr=0;
  PetscMPIInt    n,rank,off=0,size,ldn,ld3;
  PetscScalar    *A,*Q;
  PetscReal      *T;

  PetscFunctionBegin;
  if (ds->compact) kr = 3*ld;
  else k = (ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += (ds->n-l)*ld;
  if (eigr) k += (ds->n-l);
  PetscCall(DSAllocateWork_Private(ds,k+kr,0,0));
  PetscCall(PetscMPIIntCast(k*sizeof(PetscScalar)+kr*sizeof(PetscReal),&size));
  PetscCall(PetscMPIIntCast(ds->n-l,&n));
  PetscCall(PetscMPIIntCast(ld*(ds->n-l),&ldn));
  PetscCall(PetscMPIIntCast(ld*3,&ld3));
  if (ds->compact) PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  else PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  if (ds->state>DS_STATE_RAW) PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank));
  if (!rank) {
    if (ds->compact) PetscCallMPI(MPI_Pack(T,ld3,MPIU_REAL,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    else PetscCallMPI(MPI_Pack(A+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) PetscCallMPI(MPI_Pack(Q+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (eigr) PetscCallMPI(MPI_Pack(eigr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
  }
  PetscCallMPI(MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds)));
  if (rank) {
    if (ds->compact) PetscCallMPI(MPI_Unpack(ds->work,size,&off,T,ld3,MPIU_REAL,PetscObjectComm((PetscObject)ds)));
    else PetscCallMPI(MPI_Unpack(ds->work,size,&off,A+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) PetscCallMPI(MPI_Unpack(ds->work,size,&off,Q+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (eigr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
  }
  if (ds->compact) PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  else PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  if (ds->state>DS_STATE_RAW) PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSCond_HEP(DS ds,PetscReal *cond)
{
  PetscScalar    *work;
  PetscReal      *rwork;
  PetscBLASInt   *ipiv;
  PetscBLASInt   lwork,info,n,ld;
  PetscReal      hn,hin;
  PetscScalar    *A;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  lwork = 8*ld;
  PetscCall(DSAllocateWork_Private(ds,lwork,ld,ld));
  work  = ds->work;
  rwork = ds->rwork;
  ipiv  = ds->iwork;
  if (ds->compact) PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSSwitchFormat_HEP(ds));

  /* use workspace matrix W to avoid overwriting A */
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_W));
  PetscCall(MatCopy(ds->omat[DS_MAT_A],ds->omat[DS_MAT_W],SAME_NONZERO_PATTERN));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_W],&A));

  /* norm of A */
  hn = LAPACKlange_("I",&n,&n,A,&ld,rwork);

  /* norm of inv(A) */
  PetscCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n,&n,A,&ld,ipiv,&info));
  SlepcCheckLapackInfo("getrf",info);
  PetscCallBLAS("LAPACKgetri",LAPACKgetri_(&n,A,&ld,ipiv,work,&lwork,&info));
  SlepcCheckLapackInfo("getri",info);
  hin = LAPACKlange_("I",&n,&n,A,&ld,rwork);
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_W],&A));

  *cond = hn*hin;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSTranslateRKS_HEP(DS ds,PetscScalar alpha)
{
  PetscInt       i,j,k=ds->k;
  PetscScalar    *Q,*A,*R,*tau,*work;
  PetscBLASInt   ld,n1,n0,lwork,info;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(DSAllocateWork_Private(ds,ld*ld,0,0));
  tau = ds->work;
  work = ds->work+ld;
  PetscCall(PetscBLASIntCast(ld*(ld-1),&lwork));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_W));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_W],&R));

  /* copy I+alpha*A */
  PetscCall(PetscArrayzero(Q,ld*ld));
  PetscCall(PetscArrayzero(R,ld*ld));
  for (i=0;i<k;i++) {
    Q[i+i*ld] = 1.0 + alpha*A[i+i*ld];
    Q[k+i*ld] = alpha*A[k+i*ld];
  }

  /* compute qr */
  PetscCall(PetscBLASIntCast(k+1,&n1));
  PetscCall(PetscBLASIntCast(k,&n0));
  PetscCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&n1,&n0,Q,&ld,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("geqrf",info);

  /* copy R from Q */
  for (j=0;j<k;j++)
    for (i=0;i<=j;i++)
      R[i+j*ld] = Q[i+j*ld];

  /* compute orthogonal matrix in Q */
  PetscCallBLAS("LAPACKorgqr",LAPACKorgqr_(&n1,&n1,&n0,Q,&ld,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("orgqr",info);

  /* compute the updated matrix of projected problem */
  for (j=0;j<k;j++)
    for (i=0;i<k+1;i++)
      A[j*ld+i] = Q[i*ld+j];
  alpha = -1.0/alpha;
  PetscCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&n1,&n0,&alpha,R,&ld,A,&ld));
  for (i=0;i<k;i++)
    A[ld*i+i] -= alpha;

  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_W],&R));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSHermitian_HEP(DS ds,DSMatType m,PetscBool *flg)
{
  PetscFunctionBegin;
  if (m==DS_MAT_A && !ds->extrarow) *flg = PETSC_TRUE;
  else *flg = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSetCompact_HEP(DS ds,PetscBool comp)
{
  PetscFunctionBegin;
  if (!comp) PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSReallocate_HEP(DS ds,PetscInt ld)
{
  PetscInt i,*perm=ds->perm;

  PetscFunctionBegin;
  for (i=0;i<DS_NUM_MAT;i++) {
    if (!ds->compact && i==DS_MAT_A) continue;
    if (i!=DS_MAT_Q && i!=DS_MAT_T) PetscCall(MatDestroy(&ds->omat[i]));
  }

  if (!ds->compact) PetscCall(DSReallocateMat_Private(ds,DS_MAT_A,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_Q,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_T,ld));

  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscCall(PetscArraycpy(ds->perm,perm,ds->ld));
  PetscCall(PetscFree(perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSHEP - Dense Hermitian Eigenvalue Problem.

   Level: beginner

   Notes:
   The problem is expressed as A*X = X*Lambda, where A is real symmetric
   (or complex Hermitian). Lambda is a diagonal matrix whose diagonal
   elements are the arguments of DSSolve(). After solve, A is overwritten
   with Lambda.

   In the intermediate state A is reduced to tridiagonal form. In compact
   storage format, the symmetric tridiagonal matrix is stored in T.

   Used DS matrices:
+  DS_MAT_A - problem matrix (used only if compact=false)
.  DS_MAT_T - symmetric tridiagonal matrix
-  DS_MAT_Q - orthogonal/unitary transformation that reduces to tridiagonal form
   (intermediate step) or matrix of orthogonal eigenvectors, which is equal to X

   Implemented methods:
+  0 - Implicit QR (_steqr)
.  1 - Multiple Relatively Robust Representations (_stevr)
.  2 - Divide and Conquer (_stedc)
-  3 - Block Divide and Conquer (real scalars only)

.seealso: DSCreate(), DSSetType(), DSType
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_HEP(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate      = DSAllocate_HEP;
  ds->ops->view          = DSView_HEP;
  ds->ops->vectors       = DSVectors_HEP;
  ds->ops->solve[0]      = DSSolve_HEP_QR;
  ds->ops->solve[1]      = DSSolve_HEP_MRRR;
  ds->ops->solve[2]      = DSSolve_HEP_DC;
#if !defined(PETSC_USE_COMPLEX)
  ds->ops->solve[3]      = DSSolve_HEP_BDC;
#endif
  ds->ops->sort          = DSSort_HEP;
  ds->ops->truncate      = DSTruncate_HEP;
  ds->ops->update        = DSUpdateExtraRow_HEP;
  ds->ops->cond          = DSCond_HEP;
  ds->ops->transrks      = DSTranslateRKS_HEP;
  ds->ops->hermitian     = DSHermitian_HEP;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize   = DSSynchronize_HEP;
#endif
  ds->ops->setcompact    = DSSetCompact_HEP;
  ds->ops->reallocate    = DSReallocate_HEP;
  PetscFunctionReturn(PETSC_SUCCESS);
}
