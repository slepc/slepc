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
  PetscInt m;              /* number of columns */
  PetscInt t;              /* number of rows of V after truncating */
} DS_SVD;

static PetscErrorCode DSAllocate_SVD(DS ds,PetscInt ld)
{
  PetscFunctionBegin;
  if (!ds->compact) PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_U));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_V));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_T));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*   0       l           k                 m-1
    -----------------------------------------
    |*       .           .                  |
    |  *     .           .                  |
    |    *   .           .                  |
    |      * .           .                  |
    |        o           o                  |
    |          o         o                  |
    |            o       o                  |
    |              o     o                  |
    |                o   o                  |
    |                  o o                  |
    |                    o x                |
    |                      x x              |
    |                        x x            |
    |                          x x          |
    |                            x x        |
    |                              x x      |
    |                                x x    |
    |                                  x x  |
    |                                    x x|
n-1 |                                      x|
    -----------------------------------------
*/

static PetscErrorCode DSSwitchFormat_SVD(DS ds)
{
  DS_SVD         *ctx = (DS_SVD*)ds->data;
  PetscReal      *T;
  PetscScalar    *A;
  PetscInt       i,m=ctx->m,k=ds->k,ld=ds->ld;

  PetscFunctionBegin;
  PetscCheck(m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the number of columns with DSSVDSetDimensions()");
  PetscCheck(ds->compact,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Must have compact storage");
  /* switch from compact (arrow) to dense storage */
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_A],&A));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  PetscCall(PetscArrayzero(A,ld*ld));
  for (i=0;i<k;i++) {
    A[i+i*ld] = T[i];
    A[i+k*ld] = T[i+ld];
  }
  A[k+k*ld] = T[k];
  for (i=k+1;i<m;i++) {
    A[i+i*ld]   = T[i];
    A[i-1+i*ld] = T[i-1+ld];
  }
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_A],&A));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSView_SVD(DS ds,PetscViewer viewer)
{
  DS_SVD            *ctx = (DS_SVD*)ds->data;
  PetscViewerFormat format;
  PetscInt          i,j,r,c,m=ctx->m,rows,cols;
  PetscReal         *T,value;
  const char        *methodname[] = {
                     "Implicit zero-shift QR for bidiagonals (_bdsqr)",
                     "Divide and Conquer (_bdsdc or _gesdd)"
  };
  const int         nmeth=PETSC_STATIC_ARRAY_LENGTH(methodname);

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"number of columns: %" PetscInt_FMT "\n",m));
    if (ds->method<nmeth) PetscCall(PetscViewerASCIIPrintf(viewer,"solving the problem with: %s\n",methodname[ds->method]));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCheck(m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the number of columns with DSSVDSetDimensions()");
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    rows = ds->n;
    cols = ds->extrarow? m+1: m;
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%% Size = %" PetscInt_FMT " %" PetscInt_FMT "\n",rows,cols));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = zeros(%" PetscInt_FMT ",3);\n",2*ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = [\n"));
      for (i=0;i<PetscMin(ds->n,m);i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,i+1,(double)T[i]));
      for (i=0;i<cols-1;i++) {
        r = PetscMax(i+2,ds->k+1);
        c = i+1;
        PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",c,r,(double)T[i+ds->ld]));
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(zzz);\n",DSMatName[DS_MAT_T]));
    } else {
      for (i=0;i<rows;i++) {
        for (j=0;j<cols;j++) {
          if (i==j) value = T[i];
          else if (i<ds->k && j==ds->k) value = T[PetscMin(i,j)+ds->ld];
          else if (i+1==j && i>=ds->k) value = T[i+ds->ld];
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
  if (ds->state>DS_STATE_INTERMEDIATE) {
    PetscCall(DSViewMat(ds,viewer,DS_MAT_U));
    PetscCall(DSViewMat(ds,viewer,DS_MAT_V));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_SVD(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_U:
    case DS_MAT_V:
      if (rnorm) *rnorm = 0.0;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_SVD(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  DS_SVD         *ctx = (DS_SVD*)ds->data;
  PetscInt       n,l,i,*perm,ld=ds->ld;
  PetscScalar    *A;
  PetscReal      *d;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the number of columns with DSSVDSetDimensions()");
  l = ds->l;
  n = PetscMin(ds->n,ctx->m);
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  perm = ds->perm;
  if (!rr) PetscCall(DSSortEigenvaluesReal_Private(ds,d,perm));
  else PetscCall(DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_FALSE));
  for (i=l;i<n;i++) wr[i] = d[perm[i]];
  PetscCall(DSPermuteBoth_Private(ds,l,n,ds->n,ctx->m,DS_MAT_U,DS_MAT_V,perm));
  for (i=l;i<n;i++) d[i] = PetscRealPart(wr[i]);
  if (!ds->compact) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    for (i=l;i<n;i++) A[i+i*ld] = wr[i];
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSUpdateExtraRow_SVD(DS ds)
{
  DS_SVD            *ctx = (DS_SVD*)ds->data;
  PetscInt          i;
  PetscBLASInt      n=0,m=0,ld,incx=1;
  PetscScalar       *A,*x,*y,one=1.0,zero=0.0;
  PetscReal         *T,*e,beta;
  const PetscScalar *U;

  PetscFunctionBegin;
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the number of columns with DSSVDSetDimensions()");
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ctx->m,&m));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_U],&U));
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    e = T+ld;
    beta = e[m-1];   /* in compact, we assume all entries are zero except the last one */
    for (i=0;i<n;i++) e[i] = PetscRealPart(beta*U[n-1+i*ld]);
    ds->k = m;
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  } else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    PetscCall(DSAllocateWork_Private(ds,2*ld,0,0));
    x = ds->work;
    y = ds->work+ld;
    for (i=0;i<n;i++) x[i] = PetscConj(A[i+m*ld]);
    PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&one,U,&ld,x,&incx,&zero,y,&incx));
    for (i=0;i<n;i++) A[i+m*ld] = PetscConj(y[i]);
    ds->k = m;
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_U],&U));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSTruncate_SVD(DS ds,PetscInt n,PetscBool trim)
{
  PetscInt    i,ld=ds->ld,l=ds->l;
  PetscScalar *A;
  DS_SVD      *ctx = (DS_SVD*)ds->data;

  PetscFunctionBegin;
  if (!ds->compact && ds->extrarow) PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  if (trim) {
    if (!ds->compact && ds->extrarow) {   /* clean extra column */
      for (i=l;i<ds->n;i++) A[i+ctx->m*ld] = 0.0;
    }
    ds->l  = 0;
    ds->k  = 0;
    ds->n  = n;
    ctx->m = n;
    ds->t  = ds->n;   /* truncated length equal to the new dimension */
    ctx->t = ctx->m;  /* must also keep the previous dimension of V */
  } else {
    if (!ds->compact && ds->extrarow && ds->k==ds->n) {
      /* copy entries of extra column to the new position, then clean last row */
      for (i=l;i<n;i++) A[i+n*ld] = A[i+ctx->m*ld];
      for (i=l;i<ds->n;i++) A[i+ctx->m*ld] = 0.0;
    }
    ds->k  = (ds->extrarow)? n: 0;
    ds->t  = ds->n;   /* truncated length equal to previous dimension */
    ctx->t = ctx->m;  /* must also keep the previous dimension of V */
    ds->n  = n;
    ctx->m = n;
  }
  if (!ds->compact && ds->extrarow) PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  DSArrowBidiag reduces a real square arrowhead matrix of the form

                [ d 0 0 0 e ]
                [ 0 d 0 0 e ]
            A = [ 0 0 d 0 e ]
                [ 0 0 0 d e ]
                [ 0 0 0 0 d ]

  to upper bidiagonal form

                [ d e 0 0 0 ]
                [ 0 d e 0 0 ]
   B = Q'*A*P = [ 0 0 d e 0 ],
                [ 0 0 0 d e ]
                [ 0 0 0 0 d ]

  where P,Q are orthogonal matrices. Uses plane rotations with a bulge chasing scheme.
  On input, P and Q must be initialized to the identity matrix.
*/
static PetscErrorCode DSArrowBidiag(PetscBLASInt n,PetscReal *d,PetscReal *e,PetscScalar *Q,PetscBLASInt ldq,PetscScalar *P,PetscBLASInt ldp)
{
  PetscBLASInt i,j,j2,one=1;
  PetscReal    c,s,ct,st,off,temp0,temp1,temp2;

  PetscFunctionBegin;
  if (n<=2) PetscFunctionReturn(PETSC_SUCCESS);

  for (j=0;j<n-2;j++) {

    /* Eliminate entry e(j) by a rotation in the planes (j,j+1) */
    temp0 = e[j+1];
    PetscCallBLAS("LAPACKlartg",LAPACKREALlartg_(&temp0,&e[j],&c,&s,&e[j+1]));
    s = -s;

    /* Apply rotation to Q */
    j2 = j+2;
    PetscCallBLAS("BLASrot",BLASMIXEDrot_(&j2,Q+j*ldq,&one,Q+(j+1)*ldq,&one,&c,&s));

    /* Apply rotation to diagonal elements, eliminate newly introduced entry A(j+1,j) */
    temp0 = d[j+1];
    temp1 = c*temp0;
    temp2 = -s*d[j];
    PetscCallBLAS("LAPACKlartg",LAPACKREALlartg_(&temp1,&temp2,&ct,&st,&d[j+1]));
    st = -st;
    e[j] = -c*st*d[j] + s*ct*temp0;
    d[j] = c*ct*d[j] + s*st*temp0;

    /* Apply rotation to P */
    PetscCallBLAS("BLASrot",BLASMIXEDrot_(&j2,P+j*ldp,&one,P+(j+1)*ldp,&one,&ct,&st));

    /* Chase newly introduced off-diagonal entry to the top left corner */
    for (i=j-1;i>=0;i--) {

      /* Upper bulge */
      off   = -st*e[i];
      e[i]  = ct*e[i];
      temp0 = e[i+1];
      PetscCallBLAS("LAPACKlartg",LAPACKREALlartg_(&temp0,&off,&c,&s,&e[i+1]));
      s = -s;
      PetscCallBLAS("BLASrot",BLASMIXEDrot_(&j2,Q+i*ldq,&one,Q+(i+1)*ldq,&one,&c,&s));

      /* Lower bulge */
      temp0 = d[i+1];
      temp1 = -s*e[i] + c*temp0;
      temp2 = c*e[i] + s*temp0;
      off   = -s*d[i];
      PetscCallBLAS("LAPACKlartg",LAPACKREALlartg_(&temp1,&off,&ct,&st,&d[i+1]));
      st = -st;
      e[i] = -c*st*d[i] + ct*temp2;
      d[i] = c*ct*d[i] + st*temp2;
      PetscCallBLAS("BLASrot",BLASMIXEDrot_(&j2,P+i*ldp,&one,P+(i+1)*ldp,&one,&ct,&st));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Reduce to bidiagonal form by means of DSArrowBidiag.
*/
static PetscErrorCode DSIntermediate_SVD(DS ds)
{
  DS_SVD        *ctx = (DS_SVD*)ds->data;
  PetscInt      i,j;
  PetscBLASInt  n1 = 0,n2,m2,lwork,info,l = 0,n = 0,m = 0,nm,ld,off;
  PetscScalar   *A,*U,*V,*W,*work,*tauq,*taup;
  PetscReal     *d,*e;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ctx->m,&m));
  PetscCall(PetscBLASIntCast(ds->l,&l));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(PetscMax(0,ds->k-l+1),&n1)); /* size of leading block, excl. locked */
  n2 = n-l;     /* n2 = n1 + size of trailing block */
  m2 = m-l;
  off = l+l*ld;
  nm = PetscMin(n,m);
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d+ld;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_V],&V));
  PetscCall(PetscArrayzero(U,ld*ld));
  for (i=0;i<n;i++) U[i+i*ld] = 1.0;
  PetscCall(PetscArrayzero(V,ld*ld));
  for (i=0;i<m;i++) V[i+i*ld] = 1.0;

  if (ds->compact) {

    if (ds->state<DS_STATE_INTERMEDIATE) PetscCall(DSArrowBidiag(n1,d+l,e+l,U+off,ld,V+off,ld));

  } else {

    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    for (i=0;i<l;i++) { d[i] = PetscRealPart(A[i+i*ld]); e[i] = 0.0; }

    if (ds->state<DS_STATE_INTERMEDIATE) {
      lwork = (m+n)*16;
      PetscCall(DSAllocateWork_Private(ds,2*nm+ld*ld+lwork,0,0));
      tauq = ds->work;
      taup = ds->work+nm;
      W    = ds->work+2*nm;
      work = ds->work+2*nm+ld*ld;
      for (j=0;j<m;j++) PetscCall(PetscArraycpy(W+j*ld,A+j*ld,n));
      PetscCallBLAS("LAPACKgebrd",LAPACKgebrd_(&n2,&m2,W+off,&ld,d+l,e+l,tauq,taup,work,&lwork,&info));
      SlepcCheckLapackInfo("gebrd",info);
      PetscCallBLAS("LAPACKormbr",LAPACKormbr_("Q","L","N",&n2,&n2,&m2,W+off,&ld,tauq,U+off,&ld,work,&lwork,&info));
      SlepcCheckLapackInfo("ormbr",info);
      PetscCallBLAS("LAPACKormbr",LAPACKormbr_("P","R","N",&m2,&m2,&n2,W+off,&ld,taup,V+off,&ld,work,&lwork,&info));
      SlepcCheckLapackInfo("ormbr",info);
    } else {
      /* copy bidiagonal to d,e */
      for (i=l;i<nm;i++)   d[i] = PetscRealPart(A[i+i*ld]);
      for (i=l;i<nm-1;i++) e[i] = PetscRealPart(A[i+(i+1)*ld]);
    }
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_V],&V));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_SVD_QR(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  DS_SVD         *ctx = (DS_SVD*)ds->data;
  PetscInt       i,j;
  PetscBLASInt   n1,m1,info,l = 0,n = 0,m = 0,nm,ld,off,zero=0;
  PetscScalar    *A,*U,*V,*Vt;
  PetscReal      *d,*e;

  PetscFunctionBegin;
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the number of columns with DSSVDSetDimensions()");
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ctx->m,&m));
  PetscCall(PetscBLASIntCast(ds->l,&l));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  n1 = n-l;     /* n1 = size of leading block, excl. locked + size of trailing block */
  m1 = m-l;
  nm = PetscMin(n1,m1);
  off = l+l*ld;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d+ld;

  /* Reduce to bidiagonal form */
  PetscCall(DSIntermediate_SVD(ds));

  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_V],&V));

  /* solve bidiagonal SVD problem */
  for (i=0;i<l;i++) wr[i] = d[i];
  PetscCall(DSAllocateWork_Private(ds,ld*ld,4*n1,0));
  Vt = ds->work;
  for (i=l;i<m;i++) {
    for (j=l;j<m;j++) {
      Vt[i+j*ld] = PetscConj(V[j+i*ld]);  /* Lapack expects transposed VT */
    }
  }
  PetscCallBLAS("LAPACKbdsqr",LAPACKbdsqr_(n>=m?"U":"L",&nm,&m1,&n1,&zero,d+l,e+l,Vt+off,&ld,U+off,&ld,NULL,&ld,ds->rwork,&info));
  SlepcCheckLapackInfo("bdsqr",info);
  for (i=l;i<m;i++) {
    for (j=l;j<m;j++) {
      V[i+j*ld] = PetscConj(Vt[j+i*ld]);  /* transpose VT returned by Lapack */
    }
  }
  for (i=l;i<PetscMin(ds->n,ctx->m);i++) wr[i] = d[i];

  /* create diagonal matrix as a result */
  if (ds->compact) PetscCall(PetscArrayzero(e,n-1));
  else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    for (i=l;i<m;i++) PetscCall(PetscArrayzero(A+l+i*ld,n-l));
    for (i=l;i<PetscMin(ds->n,ctx->m);i++) A[i+i*ld] = d[i];
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_V],&V));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSolve_SVD_DC(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  DS_SVD         *ctx = (DS_SVD*)ds->data;
  PetscInt       i,j;
  PetscBLASInt   n1,m1,info,l = 0,n = 0,m = 0,nm,ld,off,lwork;
  PetscScalar    *A,*U,*V,*W,qwork;
  PetscReal      *d,*e,*Ur,*Vr;

  PetscFunctionBegin;
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the number of columns with DSSVDSetDimensions()");
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ctx->m,&m));
  PetscCall(PetscBLASIntCast(ds->l,&l));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  n1 = n-l;     /* n1 = size of leading block, excl. locked + size of trailing block */
  m1 = m-l;
  off = l+l*ld;
  if (ds->compact) PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_V],&V));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  e = d+ld;
  PetscCall(PetscArrayzero(U,ld*ld));
  for (i=0;i<l;i++) U[i+i*ld] = 1.0;
  PetscCall(PetscArrayzero(V,ld*ld));
  for (i=0;i<l;i++) V[i+i*ld] = 1.0;

  if (ds->state>DS_STATE_RAW) {
    /* solve bidiagonal SVD problem */
    for (i=0;i<l;i++) wr[i] = d[i];
#if defined(PETSC_USE_COMPLEX)
    PetscCall(DSAllocateWork_Private(ds,0,3*n1*n1+4*n1+2*ld*ld,8*n1));
    Ur = ds->rwork+3*n1*n1+4*n1;
    Vr = ds->rwork+3*n1*n1+4*n1+ld*ld;
#else
    PetscCall(DSAllocateWork_Private(ds,0,3*n1*n1+4*n1+ld*ld,8*n1));
    Ur = U;
    Vr = ds->rwork+3*n1*n1+4*n1;
#endif
    PetscCallBLAS("LAPACKbdsdc",LAPACKbdsdc_("U","I",&n1,d+l,e+l,Ur+off,&ld,Vr+off,&ld,NULL,NULL,ds->rwork,ds->iwork,&info));
    SlepcCheckLapackInfo("bdsdc",info);
    for (i=l;i<n;i++) {
      for (j=l;j<n;j++) {
#if defined(PETSC_USE_COMPLEX)
        U[i+j*ld] = Ur[i+j*ld];
#endif
        V[i+j*ld] = PetscConj(Vr[j+i*ld]);  /* transpose VT returned by Lapack */
      }
    }
  } else {
    /* solve general rectangular SVD problem */
    PetscCall(DSAllocateMat_Private(ds,DS_MAT_W));
    PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_W],&W));
    if (ds->compact) PetscCall(DSSwitchFormat_SVD(ds));
    for (i=0;i<l;i++) wr[i] = d[i];
    nm = PetscMin(n,m);
    PetscCall(DSAllocateWork_Private(ds,0,0,8*nm));
    lwork = -1;
#if defined(PETSC_USE_COMPLEX)
    PetscCall(DSAllocateWork_Private(ds,0,5*nm*nm+7*nm,0));
    PetscCallBLAS("LAPACKgesdd",LAPACKgesdd_("A",&n1,&m1,A+off,&ld,d+l,U+off,&ld,W+off,&ld,&qwork,&lwork,ds->rwork,ds->iwork,&info));
#else
    PetscCallBLAS("LAPACKgesdd",LAPACKgesdd_("A",&n1,&m1,A+off,&ld,d+l,U+off,&ld,W+off,&ld,&qwork,&lwork,ds->iwork,&info));
#endif
    SlepcCheckLapackInfo("gesdd",info);
    PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(qwork),&lwork));
    PetscCall(DSAllocateWork_Private(ds,lwork,0,0));
#if defined(PETSC_USE_COMPLEX)
    PetscCallBLAS("LAPACKgesdd",LAPACKgesdd_("A",&n1,&m1,A+off,&ld,d+l,U+off,&ld,W+off,&ld,ds->work,&lwork,ds->rwork,ds->iwork,&info));
#else
    PetscCallBLAS("LAPACKgesdd",LAPACKgesdd_("A",&n1,&m1,A+off,&ld,d+l,U+off,&ld,W+off,&ld,ds->work,&lwork,ds->iwork,&info));
#endif
    SlepcCheckLapackInfo("gesdd",info);
    for (i=l;i<m;i++) {
      for (j=l;j<m;j++) V[i+j*ld] = PetscConj(W[j+i*ld]);  /* transpose VT returned by Lapack */
    }
    PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_W],&W));
  }
  for (i=l;i<PetscMin(ds->n,ctx->m);i++) wr[i] = d[i];

  /* create diagonal matrix as a result */
  if (ds->compact) PetscCall(PetscArrayzero(e,n-1));
  else {
    for (i=l;i<m;i++) PetscCall(PetscArrayzero(A+l+i*ld,n-l));
    for (i=l;i<n;i++) A[i+i*ld] = d[i];
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_V],&V));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_SVD(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscInt       ld=ds->ld,l=ds->l,k=0,kr=0;
  PetscMPIInt    n,rank,off=0,size,ldn,ld3;
  PetscScalar    *A,*U,*V;
  PetscReal      *T;

  PetscFunctionBegin;
  if (ds->compact) kr = 3*ld;
  else k = (ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += 2*(ds->n-l)*ld;
  if (eigr) k += ds->n-l;
  PetscCall(DSAllocateWork_Private(ds,k+kr,0,0));
  PetscCall(PetscMPIIntCast(k*sizeof(PetscScalar)+kr*sizeof(PetscReal),&size));
  PetscCall(PetscMPIIntCast(ds->n-l,&n));
  PetscCall(PetscMPIIntCast(ld*(ds->n-l),&ldn));
  PetscCall(PetscMPIIntCast(3*ld,&ld3));
  if (ds->compact) PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  else PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  if (ds->state>DS_STATE_RAW) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_U],&U));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_V],&V));
  }
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank));
  if (!rank) {
    if (ds->compact) PetscCallMPI(MPI_Pack(T,ld3,MPIU_REAL,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    else PetscCallMPI(MPI_Pack(A+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) {
      PetscCallMPI(MPI_Pack(U+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Pack(V+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    }
    if (eigr) PetscCallMPI(MPI_Pack(eigr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
  }
  PetscCallMPI(MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds)));
  if (rank) {
    if (ds->compact) PetscCallMPI(MPI_Unpack(ds->work,size,&off,T,ld3,MPIU_REAL,PetscObjectComm((PetscObject)ds)));
    else PetscCallMPI(MPI_Unpack(ds->work,size,&off,A+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) {
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,U+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,V+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    }
    if (eigr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
  }
  if (ds->compact) PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  else PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  if (ds->state>DS_STATE_RAW) {
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_U],&U));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_V],&V));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSMatGetSize_SVD(DS ds,DSMatType t,PetscInt *rows,PetscInt *cols)
{
  DS_SVD *ctx = (DS_SVD*)ds->data;

  PetscFunctionBegin;
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the number of columns with DSSVDSetDimensions()");
  switch (t) {
    case DS_MAT_A:
      *rows = ds->n;
      *cols = ds->extrarow? ctx->m+1: ctx->m;
      break;
    case DS_MAT_T:
      *rows = ds->n;
      *cols = PetscDefined(USE_COMPLEX)? 2: 3;
      break;
    case DS_MAT_U:
      *rows = ds->state==DS_STATE_TRUNCATED? ds->t: ds->n;
      *cols = ds->n;
      break;
    case DS_MAT_V:
      *rows = ds->state==DS_STATE_TRUNCATED? ctx->t: ctx->m;
      *cols = ctx->m;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid t parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSVDSetDimensions_SVD(DS ds,PetscInt m)
{
  DS_SVD *ctx = (DS_SVD*)ds->data;

  PetscFunctionBegin;
  DSCheckAlloc(ds,1);
  if (m==PETSC_DECIDE || m==PETSC_DEFAULT) {
    ctx->m = ds->ld;
  } else {
    PetscCheck(m>0 && m<=ds->ld,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of m. Must be between 1 and ld");
    ctx->m = m;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSVDSetDimensions - Sets the number of columns for a DSSVD.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  m  - the number of columns

   Notes:
   This call is complementary to DSSetDimensions(), to provide a dimension
   that is specific to this DS type.

   Level: intermediate

.seealso: DSSVDGetDimensions(), DSSetDimensions()
@*/
PetscErrorCode DSSVDSetDimensions(DS ds,PetscInt m)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,m,2);
  PetscTryMethod(ds,"DSSVDSetDimensions_C",(DS,PetscInt),(ds,m));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSVDGetDimensions_SVD(DS ds,PetscInt *m)
{
  DS_SVD *ctx = (DS_SVD*)ds->data;

  PetscFunctionBegin;
  *m = ctx->m;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSVDGetDimensions - Returns the number of columns for a DSSVD.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  m - the number of columns

   Level: intermediate

.seealso: DSSVDSetDimensions()
@*/
PetscErrorCode DSSVDGetDimensions(DS ds,PetscInt *m)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(m,2);
  PetscUseMethod(ds,"DSSVDGetDimensions_C",(DS,PetscInt*),(ds,m));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSDestroy_SVD(DS ds)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(ds->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSSVDSetDimensions_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSSVDGetDimensions_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSetCompact_SVD(DS ds,PetscBool comp)
{
  PetscFunctionBegin;
  if (!comp) PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSReallocate_SVD(DS ds,PetscInt ld)
{
  PetscInt i,*perm=ds->perm;

  PetscFunctionBegin;
  for (i=0;i<DS_NUM_MAT;i++) {
    if (!ds->compact && i==DS_MAT_A) continue;
    if (i!=DS_MAT_U && i!=DS_MAT_V && i!=DS_MAT_T) PetscCall(MatDestroy(&ds->omat[i]));
  }

  if (!ds->compact) PetscCall(DSReallocateMat_Private(ds,DS_MAT_A,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_U,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_V,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_T,ld));

  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscCall(PetscArraycpy(ds->perm,perm,ds->ld));
  PetscCall(PetscFree(perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSSVD - Dense Singular Value Decomposition.

   Level: beginner

   Notes:
   The problem is expressed as A = U*Sigma*V', where A is rectangular in
   general, with n rows and m columns. Sigma is a diagonal matrix whose diagonal
   elements are the arguments of DSSolve(). After solve, A is overwritten
   with Sigma.

   The orthogonal (or unitary) matrices of left and right singular vectors, U
   and V, have size n and m, respectively. The number of columns m must
   be specified via DSSVDSetDimensions().

   If the DS object is in the intermediate state, A is assumed to be in upper
   bidiagonal form (possibly with an arrow) and is stored in compact format
   on matrix T. Otherwise, no particular structure is assumed. The compact
   storage is implemented for the square case only, m=n. The extra row should
   be interpreted in this case as an extra column.

   Used DS matrices:
+  DS_MAT_A - problem matrix (used only if compact=false)
.  DS_MAT_T - upper bidiagonal matrix
.  DS_MAT_U - left singular vectors
-  DS_MAT_V - right singular vectors

   Implemented methods:
+  0 - Implicit zero-shift QR for bidiagonals (_bdsqr)
-  1 - Divide and Conquer (_bdsdc or _gesdd)

.seealso: DSCreate(), DSSetType(), DSType, DSSVDSetDimensions()
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_SVD(DS ds)
{
  DS_SVD         *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  ds->data = (void*)ctx;

  ds->ops->allocate      = DSAllocate_SVD;
  ds->ops->view          = DSView_SVD;
  ds->ops->vectors       = DSVectors_SVD;
  ds->ops->solve[0]      = DSSolve_SVD_QR;
  ds->ops->solve[1]      = DSSolve_SVD_DC;
  ds->ops->sort          = DSSort_SVD;
  ds->ops->truncate      = DSTruncate_SVD;
  ds->ops->update        = DSUpdateExtraRow_SVD;
  ds->ops->destroy       = DSDestroy_SVD;
  ds->ops->matgetsize    = DSMatGetSize_SVD;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize   = DSSynchronize_SVD;
#endif
  ds->ops->setcompact    = DSSetCompact_SVD;
  ds->ops->reallocate    = DSReallocate_SVD;
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSSVDSetDimensions_C",DSSVDSetDimensions_SVD));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSSVDGetDimensions_C",DSSVDGetDimensions_SVD));
  PetscFunctionReturn(PETSC_SUCCESS);
}
