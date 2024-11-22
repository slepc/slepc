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
  PetscInt p;              /* number of rows of B */
  PetscInt tm;             /* number of rows of X after truncating */
  PetscInt tp;             /* number of rows of V after truncating */
} DS_GSVD;

static PetscErrorCode DSAllocate_GSVD(DS ds,PetscInt ld)
{
  PetscFunctionBegin;
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_B));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_X));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_U));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_V));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_T));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_D));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  In compact form, A is either in form (a) or (b):

                     (a)                                            (b)
    lower bidiagonal with upper arrow (n=m+1)         square upper bidiagonal with upper arrow (n=m)
     0       l           k                 m-1
    -----------------------------------------         0     l           k                   m-1
    |*                   .                  |        -----------------------------------------
    |  *                 .                  |        |*                 .                    |
    |    *               .                  |        |  *               .                    |
    |      *             .                  |        |    *             .                    |
  l |. . . . o           o                  |      l |. . . o           o                    |
    |          o         o                  |        |        o         o                    |
    |            o       o                  |        |          o       o                    |
    |              o     o                  |        |            o     o                    |
    |                o   o                  |        |              o   o                    |
    |                  o o                  |        |                o o                    |
  k |. . . . . . . . . . o                  |      k |. . . . . . . . . o x                  |
    |                    x x                |        |                    x x                |
    |                      x x              |        |                      x x              |
    |                        x x            |        |                        x x            |
    |                          x x          |        |                          x x          |
    |                            x x        |        |                            x x        |
    |                              x x      |        |                              x x      |
    |                                x x    |        |                                x x    |
    |                                  x x  |        |                                  x x  |
    |                                    x x|        |                                    x x|
n-1 |                                      x|    n-1 |                                      x|
    -----------------------------------------        -----------------------------------------

  and B is square bidiagonal with upper arrow (p=m)

     0       l           k                 m-1
    -----------------------------------------
    |*                   .                  |
    |  *                 .                  |
    |    *               .                  |
    |      *             .                  |
  l |. . . . o           o                  |
    |          o         o                  |
    |            o       o                  |
    |              o     o                  |
    |                o   o                  |
    |                  o o                  |
  k |. . . . . . . . . . o x                |
    |                      x x              |
    |                        x x            |
    |                          x x          |
    |                            x x        |
    |                              x x      |
    |                                x x    |
    |                                  x x  |
    |                                    x x|
p-1 |                                      x|
     ----------------------------------------
*/
static PetscErrorCode DSView_GSVD(DS ds,PetscViewer viewer)
{
  DS_GSVD           *ctx = (DS_GSVD*)ds->data;
  PetscViewerFormat format;
  PetscInt          i,j,r,k=ds->k,n=ds->n,m=ctx->m,p=ctx->p,rowsa,rowsb,colsa,colsb;
  PetscReal         *T,*S,value;

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO) PetscFunctionReturn(PETSC_SUCCESS);
  if (format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"number of columns: %" PetscInt_FMT "\n",m));
    PetscCall(PetscViewerASCIIPrintf(viewer,"number of rows of B: %" PetscInt_FMT "\n",p));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the other dimensions with DSGSVDSetDimensions()");
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSGetArrayReal(ds,DS_MAT_D,&S));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    rowsa = n;
    colsa = ds->extrarow? m+1: m;
    rowsb = p;
    colsb = ds->extrarow? m+1: m;
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%% Size = %" PetscInt_FMT " %" PetscInt_FMT "\n",rowsa,colsa));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = zeros(%" PetscInt_FMT ",3);\n",2*ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = [\n"));
      for (i=0;i<PetscMin(rowsa,colsa);i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,i+1,(double)T[i]));
      for (i=0;i<k;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,k+1,(double)T[i+ds->ld]));
      if (n>m) { /* A lower bidiagonal */
        for (i=k;i<rowsa-1;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+2,i+1,(double)T[i+ds->ld]));
      } else { /* A (square) upper bidiagonal */
        for (i=k;i<colsa-1;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,i+2,(double)T[i+ds->ld]));
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(zzz);\n",DSMatName[DS_MAT_T]));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%% Size = %" PetscInt_FMT " %" PetscInt_FMT "\n",rowsb,colsb));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = zeros(%" PetscInt_FMT ",3);\n",2*ds->n));
      PetscCall(PetscViewerASCIIPrintf(viewer,"zzz = [\n"));
      for (i=0;i<rowsb;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,i+1,(double)S[i]));
      for (i=0;i<colsb-1;i++) {
        r = PetscMax(i+2,ds->k+1);
        PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT " %" PetscInt_FMT "  %18.16e\n",i+1,r,(double)T[i+2*ds->ld]));
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(zzz);\n",DSMatName[DS_MAT_D]));
    } else {
      PetscCall(PetscViewerASCIIPrintf(viewer,"Matrix %s =\n",DSMatName[DS_MAT_T]));
      for (i=0;i<rowsa;i++) {
        for (j=0;j<colsa;j++) {
          if (i==j) value = T[i];
          else if (i<ds->k && j==ds->k) value = T[i+ds->ld];
          else if (n>m && i==j+1 && i>ds->k) value = T[j+ds->ld];
          else if (n<=m && i+1==j && i>=ds->k) value = T[i+ds->ld];
          else value = 0.0;
          PetscCall(PetscViewerASCIIPrintf(viewer," %18.16e ",(double)value));
        }
        PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"Matrix %s =\n",DSMatName[DS_MAT_D]));
      for (i=0;i<rowsb;i++) {
        for (j=0;j<colsb;j++) {
          if (i==j) value = S[i];
          else if (i<ds->k && j==ds->k) value = T[PetscMin(i,j)+2*ds->ld];
          else if (i+1==j && i>=ds->k) value = T[i+2*ds->ld];
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
  if (ds->state>DS_STATE_INTERMEDIATE) {
    PetscCall(DSViewMat(ds,viewer,DS_MAT_X));
    PetscCall(DSViewMat(ds,viewer,DS_MAT_U));
    PetscCall(DSViewMat(ds,viewer,DS_MAT_V));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_GSVD(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_U:
    case DS_MAT_V:
      if (rnorm) *rnorm = 0.0;
      break;
    case DS_MAT_X:
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_GSVD(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  DS_GSVD        *ctx = (DS_GSVD*)ds->data;
  PetscInt       t,l,ld=ds->ld,i,*perm,*perm2;
  PetscReal      *T=NULL,*D=NULL,*eig;
  PetscScalar    *A=NULL,*B=NULL;
  PetscBool      compact=ds->compact;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the other dimensions with DSGSVDSetDimensions()");
  l = ds->l;
  t = ds->t;
  perm = ds->perm;
  PetscCall(PetscMalloc2(t,&eig,t,&perm2));
  if (compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSGetArrayReal(ds,DS_MAT_D,&D));
    for (i=0;i<t;i++) eig[i] = (D[i]==0)?PETSC_INFINITY:T[i]/D[i];
  } else {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
    for (i=0;i<t;i++) eig[i] = (B[i+i*ld]==0)?PETSC_INFINITY:PetscRealPart(A[i+i*ld])/PetscRealPart(B[i*(1+ld)]);
  }
  PetscCall(DSSortEigenvaluesReal_Private(ds,eig,perm));
  PetscCall(PetscArraycpy(perm2,perm,t));
  for (i=l;i<t;i++) wr[i] = eig[perm[i]];
  if (compact) {
    PetscCall(PetscArraycpy(eig,T,t));
    for (i=l;i<t;i++) T[i] = eig[perm[i]];
    PetscCall(PetscArraycpy(eig,D,t));
    for (i=l;i<t;i++) D[i] = eig[perm[i]];
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&D));
  } else {
    for (i=l;i<t;i++) eig[i] = PetscRealPart(A[i*(1+ld)]);
    for (i=l;i<t;i++) A[i*(1+ld)] = eig[perm[i]];
    for (i=l;i<t;i++) eig[i] = PetscRealPart(B[i*(1+ld)]);
    for (i=l;i<t;i++) B[i*(1+ld)] = eig[perm[i]];
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  }
  PetscCall(DSPermuteColumns_Private(ds,l,t,ds->n,DS_MAT_U,perm2));
  PetscCall(PetscArraycpy(perm2,perm,t));
  PetscCall(DSPermuteColumns_Private(ds,l,t,ctx->m,DS_MAT_X,perm2));
  PetscCall(DSPermuteColumns_Private(ds,l,t,ctx->p,DS_MAT_V,perm));
  PetscCall(PetscFree2(eig,perm2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSUpdateExtraRow_GSVD(DS ds)
{
  DS_GSVD           *ctx = (DS_GSVD*)ds->data;
  PetscInt          i;
  PetscBLASInt      n=0,m=0,ld=0;
  const PetscScalar *U,*V;
  PetscReal         *T,*e,*f,alpha,beta,betah;

  PetscFunctionBegin;
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the other dimensions with DSGSVDSetDimensions()");
  PetscCheck(ds->compact,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented for non-compact storage");
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ctx->m,&m));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  e = T+ld;
  f = T+2*ld;
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_V],&V));
  if (n<=m) {   /* upper variant, A is square upper bidiagonal */
    beta  = e[m-1];   /* in compact, we assume all entries are zero except the last one */
    betah = f[m-1];
    for (i=0;i<m;i++) {
      e[i] = PetscRealPart(beta*U[m-1+i*ld]);
      f[i] = PetscRealPart(betah*V[m-1+i*ld]);
    }
  } else {   /* lower variant, A is (m+1)xm lower bidiagonal */
    alpha = T[m];
    betah = f[m-1];
    for (i=0;i<m;i++) {
      e[i] = PetscRealPart(alpha*U[m+i*ld]);
      f[i] = PetscRealPart(betah*V[m-1+i*ld]);
    }
    T[m] = PetscRealPart(alpha*U[m+m*ld]);
  }
  ds->k = m;
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_V],&V));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSTruncate_GSVD(DS ds,PetscInt n,PetscBool trim)
{
  DS_GSVD     *ctx = (DS_GSVD*)ds->data;
  PetscScalar *U;
  PetscReal   *T;
  PetscInt    i,m=ctx->m,ld=ds->ld;
  PetscBool   lower=(ds->n>ctx->m)?PETSC_TRUE:PETSC_FALSE;

  PetscFunctionBegin;
  PetscCheck(ds->compact,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented for non-compact storage");
  if (trim) {
    ds->l   = 0;
    ds->k   = 0;
    ds->n   = lower? n+1: n;
    ctx->m  = n;
    ctx->p  = n;
    ds->t   = ds->n;   /* truncated length equal to the new dimension */
    ctx->tm = ctx->m;  /* must also keep the previous dimension of X */
    ctx->tp = ctx->p;  /* must also keep the previous dimension of V */
  } else {
    if (lower) {
      /* move value of diagonal element of arrow (alpha) */
      PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
      T[n] = T[m];
      PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
      /* copy last column of U so that it updates the next initial vector of U1 */
      PetscCall(MatDenseGetArray(ds->omat[DS_MAT_U],&U));
      for (i=0;i<=m;i++) U[i+n*ld] = U[i+m*ld];
      PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_U],&U));
    }
    ds->k   = (ds->extrarow)? n: 0;
    ds->t   = ds->n;   /* truncated length equal to previous dimension */
    ctx->tm = ctx->m;  /* must also keep the previous dimension of X */
    ctx->tp = ctx->p;  /* must also keep the previous dimension of V */
    ds->n   = lower? n+1: n;
    ctx->m  = n;
    ctx->p  = n;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSwitchFormat_GSVD(DS ds)
{
  DS_GSVD        *ctx = (DS_GSVD*)ds->data;
  PetscReal      *T,*D;
  PetscScalar    *A,*B;
  PetscInt       i,n=ds->n,k=ds->k,ld=ds->ld,m=ctx->m;

  PetscFunctionBegin;
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the other dimensions with DSGSVDSetDimensions()");
  /* switch from compact (arrow) to dense storage */
  /* bidiagonal associated to B is stored in D and T+2*ld */
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArrayWrite(ds->omat[DS_MAT_B],&B));
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&D));
  PetscCall(PetscArrayzero(A,ld*ld));
  PetscCall(PetscArrayzero(B,ld*ld));
  for (i=0;i<k;i++) {
    A[i+i*ld] = T[i];
    A[i+k*ld] = T[i+ld];
    B[i+i*ld] = D[i];
    B[i+k*ld] = T[i+2*ld];
  }
  /* B is upper bidiagonal */
  B[k+k*ld] = D[k];
  for (i=k+1;i<m;i++) {
    B[i+i*ld]   = D[i];
    B[i-1+i*ld] = T[i-1+2*ld];
  }
  /* A can be upper (square) or lower bidiagonal */
  for (i=k;i<m;i++) A[i+i*ld] = T[i];
  if (n>m) for (i=k;i<m;i++) A[i+1+i*ld] = T[i+ld];
  else for (i=k+1;i<m;i++) A[i-1+i*ld] = T[i-1+ld];
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArrayWrite(ds->omat[DS_MAT_B],&B));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&D));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compact format is used when [A;B] has orthonormal columns.
  In this case R=I and the GSVD of (A,B) is the CS decomposition
*/
static PetscErrorCode DSSolve_GSVD(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  DS_GSVD        *ctx = (DS_GSVD*)ds->data;
  PetscInt       i,j;
  PetscBLASInt   n1,m1,info,lc = 0,n = 0,m = 0,p = 0,p1,l,k,q,ld,off,lwork,r;
  PetscScalar    *A,*B,*X,*U,*V,sone=1.0,smone=-1.0;
  PetscReal      *alpha,*beta,*T,*D;
#if !defined(SLEPC_MISSING_LAPACK_GGSVD3)
  PetscScalar    a,dummy;
  PetscReal      rdummy;
  PetscBLASInt   idummy;
#endif

  PetscFunctionBegin;
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the other dimensions with DSGSVDSetDimensions()");
  PetscCall(PetscBLASIntCast(ds->n,&m));
  PetscCall(PetscBLASIntCast(ctx->m,&n));
  PetscCall(PetscBLASIntCast(ctx->p,&p));
  PetscCall(PetscBLASIntCast(ds->l,&lc));
  PetscCheck(ds->compact || lc==0,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"DSGSVD with non-compact format does not support locking");
  /* In compact storage B is always nxn and A can be either nxn or (n+1)xn */
  PetscCheck(!ds->compact || (p==n && (m==p || m==p+1)),PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Dimensions not supported in compact format");
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  n1 = n-lc;     /* n1 = size of leading block, excl. locked + size of trailing block */
  m1 = m-lc;
  p1 = p-lc;
  off = lc+lc*ld;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_V],&V));
  PetscCall(PetscArrayzero(X,ld*ld));
  for (i=0;i<lc;i++) X[i+i*ld] = 1.0;
  PetscCall(PetscArrayzero(U,ld*ld));
  for (i=0;i<lc;i++) U[i+i*ld] = 1.0;
  PetscCall(PetscArrayzero(V,ld*ld));
  for (i=0;i<lc;i++) V[i+i*ld] = 1.0;
  if (ds->compact) PetscCall(DSSwitchFormat_GSVD(ds));

#if !defined(SLEPC_MISSING_LAPACK_GGSVD3)
  /* workspace query and memory allocation */
  lwork = -1;
#if !defined (PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKggsvd3",LAPACKggsvd3_("U","V","Q",&m1,&n1,&p1,&k,&l,&dummy,&ld,&dummy,&ld,&rdummy,&rdummy,&dummy,&ld,&dummy,&ld,&dummy,&ld,&a,&lwork,&idummy,&info));
  PetscCall(PetscBLASIntCast((PetscInt)a,&lwork));
#else
  PetscCallBLAS("LAPACKggsvd3",LAPACKggsvd3_("U","V","Q",&m1,&n1,&p1,&k,&l,&dummy,&ld,&dummy,&ld,&rdummy,&rdummy,&dummy,&ld,&dummy,&ld,&dummy,&ld,&a,&lwork,&rdummy,&idummy,&info));
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork));
#endif

#if !defined (PETSC_USE_COMPLEX)
  PetscCall(DSAllocateWork_Private(ds,lwork,2*ds->ld,ds->ld));
  alpha = ds->rwork;
  beta  = ds->rwork+ds->ld;
  PetscCallBLAS("LAPACKggsvd3",LAPACKggsvd3_("U","V","Q",&m1,&n1,&p1,&k,&l,A+off,&ld,B+off,&ld,alpha,beta,U+off,&ld,V+off,&ld,X+off,&ld,ds->work,&lwork,ds->iwork,&info));
#else
  PetscCall(DSAllocateWork_Private(ds,lwork,4*ds->ld,ds->ld));
  alpha = ds->rwork+2*ds->ld;
  beta  = ds->rwork+3*ds->ld;
  PetscCallBLAS("LAPACKggsvd3",LAPACKggsvd3_("U","V","Q",&m1,&n1,&p1,&k,&l,A+off,&ld,B+off,&ld,alpha,beta,U+off,&ld,V+off,&ld,X+off,&ld,ds->work,&lwork,ds->rwork,ds->iwork,&info));
#endif
  SlepcCheckLapackInfo("ggsvd3",info);

#else  /* defined(SLEPC_MISSING_LAPACK_GGSVD3) */

  lwork = PetscMax(PetscMax(3*n,m),p)+n;
#if !defined (PETSC_USE_COMPLEX)
  PetscCall(DSAllocateWork_Private(ds,lwork,2*ds->ld,ds->ld));
  alpha = ds->rwork;
  beta  = ds->rwork+ds->ld;
  PetscCallBLAS("LAPACKggsvd",LAPACKggsvd_("U","V","Q",&m1,&n1,&p1,&k,&l,A+off,&ld,B+off,&ld,alpha,beta,U+off,&ld,V+off,&ld,X+off,&ld,ds->work,ds->iwork,&info));
#else
  PetscCall(DSAllocateWork_Private(ds,lwork,4*ds->ld,ds->ld));
  alpha = ds->rwork+2*ds->ld;
  beta  = ds->rwork+3*ds->ld;
  PetscCallBLAS("LAPACKggsvd",LAPACKggsvd_("U","V","Q",&m1,&n1,&p1,&k,&l,A+off,&ld,B+off,&ld,alpha,beta,U+off,&ld,V+off,&ld,X+off,&ld,ds->work,ds->rwork,ds->iwork,&info));
#endif
  SlepcCheckLapackInfo("ggsvd",info);

#endif

  PetscCheck(k+l>=n1,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"The rank deficient case not supported yet");
  if (ds->compact) {
    PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
    PetscCall(DSGetArrayReal(ds,DS_MAT_D,&D));
    /* R is the identity matrix (except the sign) */
    for (i=lc;i<n;i++) {
      if (PetscRealPart(A[i+i*ld])<0.0) { /* scale column i */
        for (j=lc;j<n;j++) X[j+i*ld] = -X[j+i*ld];
      }
    }
    PetscCall(PetscArrayzero(T+ld,m-1));
    PetscCall(PetscArrayzero(T+2*ld,n-1));
    for (i=lc;i<n;i++) {
      T[i] = alpha[i-lc];
      D[i] = beta[i-lc];
      if (D[i]==0.0) wr[i] = PETSC_INFINITY;
      else wr[i] = T[i]/D[i];
    }
    ds->t = n;
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&D));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  } else {
    /* X = X*inv(R) */
    q = PetscMin(m,n);
    PetscCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&n,&q,&sone,A,&ld,X,&ld));
    if (m<n) {
      r = n-m;
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&r,&m,&sone,X,&ld,A,&ld,&smone,X+m*ld,&ld));
      PetscCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&n,&r,&sone,B+m*ld,&ld,X+m*ld,&ld));
    }
    if (k>0) {
      for (i=k;i<PetscMin(m,k+l);i++) {
        PetscCall(PetscArraycpy(X+(i-k)*ld,X+i*ld,ld));
        PetscCall(PetscArraycpy(U+(i-k)*ld,U+i*ld,ld));
      }
    }
    /* singular values */
    PetscCall(PetscArrayzero(A,ld*ld));
    PetscCall(PetscArrayzero(B,ld*ld));
    for (j=k;j<PetscMin(m,k+l);j++) {
      A[(j-k)*(1+ld)] = alpha[j];
      B[(j-k)*(1+ld)] = beta[j];
      wr[j-k] = alpha[j]/beta[j];
    }
    ds->t = PetscMin(m,k+l)-k; /* set number of computed values */
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_V],&V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSCond_GSVD(DS ds,PetscReal *cond)
{
  DS_GSVD           *ctx = (DS_GSVD*)ds->data;
  PetscBLASInt      lwork,lrwork=0,info,m,n,p,ld;
  PetscScalar       *A,*work;
  const PetscScalar *M;
  PetscReal         *sigma,conda,condb;
#if defined(PETSC_USE_COMPLEX)
  PetscReal         *rwork;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&m));
  PetscCall(PetscBLASIntCast(ctx->m,&n));
  PetscCall(PetscBLASIntCast(ctx->p,&p));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  lwork = 5*n;
#if defined(PETSC_USE_COMPLEX)
  lrwork = 5*n;
#endif
  PetscCall(DSAllocateWork_Private(ds,ld*n+lwork,n+lrwork,0));
  A     = ds->work;
  work  = ds->work+ld*n;
  sigma = ds->rwork;
#if defined(PETSC_USE_COMPLEX)
  rwork = ds->rwork+n;
#endif
  if (ds->compact) PetscCall(DSSwitchFormat_GSVD(ds));

  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&M));
  PetscCall(PetscArraycpy(A,M,ld*n));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&M));
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","N",&m,&n,A,&ld,sigma,NULL,&ld,NULL,&ld,work,&lwork,rwork,&info));
#else
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","N",&m,&n,A,&ld,sigma,NULL,&ld,NULL,&ld,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);
  conda = sigma[0]/sigma[PetscMin(m,n)-1];

  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&M));
  PetscCall(PetscArraycpy(A,M,ld*n));
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&M));
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","N",&p,&n,A,&ld,sigma,NULL,&ld,NULL,&ld,work,&lwork,rwork,&info));
#else
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","N",&p,&n,A,&ld,sigma,NULL,&ld,NULL,&ld,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);
  condb = sigma[0]/sigma[PetscMin(p,n)-1];

  *cond = PetscMax(conda,condb);
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_GSVD(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  DS_GSVD        *ctx = (DS_GSVD*)ds->data;
  PetscInt       ld=ds->ld,l=ds->l,k=0,kr=0;
  PetscMPIInt    m,rank,off=0,size,n,ldn,ld3;
  PetscScalar    *A,*U,*V,*X;
  PetscReal      *T;

  PetscFunctionBegin;
  PetscCall(PetscMPIIntCast(ctx->m,&m));
  if (ds->compact) kr = 3*ld;
  else k = 2*(m-l)*ld;
  if (ds->state>DS_STATE_RAW) k += 3*(m-l)*ld;
  if (eigr) k += m-l;
  PetscCall(DSAllocateWork_Private(ds,k+kr,0,0));
  PetscCall(PetscMPIIntCast(k*sizeof(PetscScalar)+kr*sizeof(PetscReal),&size));
  PetscCall(PetscMPIIntCast(m-l,&n));
  PetscCall(PetscMPIIntCast(ld*(m-l),&ldn));
  PetscCall(PetscMPIIntCast(3*ld,&ld3));
  if (ds->compact) PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  else PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  if (ds->state>DS_STATE_RAW) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_U],&U));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_V],&V));
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
  }
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank));
  if (!rank) {
    if (ds->compact) PetscCallMPI(MPI_Pack(T,ld3,MPIU_REAL,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    else PetscCallMPI(MPI_Pack(A+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (ds->state>DS_STATE_RAW) {
      PetscCallMPI(MPI_Pack(U+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Pack(V+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
      PetscCallMPI(MPI_Pack(X+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
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
      PetscCallMPI(MPI_Unpack(ds->work,size,&off,X+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    }
    if (eigr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
  }
  if (ds->compact) PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  else PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  if (ds->state>DS_STATE_RAW) {
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_U],&U));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_V],&V));
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSMatGetSize_GSVD(DS ds,DSMatType t,PetscInt *rows,PetscInt *cols)
{
  DS_GSVD *ctx = (DS_GSVD*)ds->data;

  PetscFunctionBegin;
  PetscCheck(ctx->m,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"You should set the other dimensions with DSGSVDSetDimensions()");
  switch (t) {
    case DS_MAT_A:
      *rows = ds->n;
      *cols = ds->extrarow? ctx->m+1: ctx->m;
      break;
    case DS_MAT_B:
      *rows = ctx->p;
      *cols = ds->extrarow? ctx->m+1: ctx->m;
      break;
    case DS_MAT_T:
      *rows = ds->n;
      *cols = PetscDefined(USE_COMPLEX)? 2: 3;
      break;
    case DS_MAT_D:
      *rows = ctx->p;
      *cols = 1;
      break;
    case DS_MAT_U:
      *rows = ds->state==DS_STATE_TRUNCATED? ds->t: ds->n;
      *cols = ds->n;
      break;
    case DS_MAT_V:
      *rows = ds->state==DS_STATE_TRUNCATED? ctx->tp: ctx->p;
      *cols = ctx->p;
      break;
    case DS_MAT_X:
      *rows = ds->state==DS_STATE_TRUNCATED? ctx->tm: ctx->m;
      *cols = ctx->m;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid t parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSGSVDSetDimensions_GSVD(DS ds,PetscInt m,PetscInt p)
{
  DS_GSVD *ctx = (DS_GSVD*)ds->data;

  PetscFunctionBegin;
  DSCheckAlloc(ds,1);
  if (m == PETSC_DETERMINE) {
    ctx->m = ds->ld;
  } else if (m != PETSC_CURRENT) {
    PetscCheck(m>0 && m<=ds->ld,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of m. Must be between 1 and ld");
    ctx->m = m;
  }
  if (p == PETSC_DETERMINE) {
    ctx->p = ds->n;
  } else if (p != PETSC_CURRENT) {
    PetscCheck(p>0 && p<=ds->ld,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of p. Must be between 1 and ld");
    ctx->p = p;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGSVDSetDimensions - Sets the number of columns and rows for a DSGSVD.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
.  m  - the number of columns
-  p  - the number of rows for the second matrix (B)

   Notes:
   This call is complementary to DSSetDimensions(), to provide two dimensions
   that are specific to this DS type. The number of rows for the first matrix (A)
   is set by DSSetDimensions().

   Use PETSC_CURRENT to leave any of the values unchanged. Use PETSC_DETERMINE
   to set m to the leading dimension and p to the number of columns of B.

   Level: intermediate

.seealso: DSGSVDGetDimensions(), DSSetDimensions()
@*/
PetscErrorCode DSGSVDSetDimensions(DS ds,PetscInt m,PetscInt p)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,m,2);
  PetscValidLogicalCollectiveInt(ds,p,3);
  PetscTryMethod(ds,"DSGSVDSetDimensions_C",(DS,PetscInt,PetscInt),(ds,m,p));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSGSVDGetDimensions_GSVD(DS ds,PetscInt *m,PetscInt *p)
{
  DS_GSVD *ctx = (DS_GSVD*)ds->data;

  PetscFunctionBegin;
  if (m) *m = ctx->m;
  if (p) *p = ctx->p;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGSVDGetDimensions - Returns the number of columns and rows for a DSGSVD.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
+  m - the number of columns
-  p - the number of rows for the second matrix (B)

   Level: intermediate

.seealso: DSGSVDSetDimensions()
@*/
PetscErrorCode DSGSVDGetDimensions(DS ds,PetscInt *m,PetscInt *p)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscUseMethod(ds,"DSGSVDGetDimensions_C",(DS,PetscInt*,PetscInt*),(ds,m,p));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSDestroy_GSVD(DS ds)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(ds->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSGSVDSetDimensions_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSGSVDGetDimensions_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSReallocate_GSVD(DS ds,PetscInt ld)
{
  PetscInt i,*perm=ds->perm;

  PetscFunctionBegin;
  for (i=0;i<DS_NUM_MAT;i++) {
    if (i!=DS_MAT_A && i!=DS_MAT_B && i!=DS_MAT_X && i!=DS_MAT_U && i!=DS_MAT_V && i!=DS_MAT_T && i!=DS_MAT_D) PetscCall(MatDestroy(&ds->omat[i]));
  }

  PetscCall(DSReallocateMat_Private(ds,DS_MAT_A,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_B,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_X,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_U,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_V,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_T,ld));
  PetscCall(DSReallocateMat_Private(ds,DS_MAT_D,ld));

  PetscCall(PetscMalloc1(ld,&ds->perm));
  PetscCall(PetscArraycpy(ds->perm,perm,ds->ld));
  PetscCall(PetscFree(perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSGSVD - Dense Generalized Singular Value Decomposition.

   Level: beginner

   Notes:
   The problem is expressed as A*X = U*C, B*X = V*S, where A and B are
   matrices with the same number of columns, m, U and V are orthogonal
   (unitary), and X is an mxm invertible matrix. The DS object does not
   expose matrices C and S, instead the singular values sigma, which are
   the ratios c_i/s_i, are returned in the arguments of DSSolve().
   Note that the number of columns of the returned X, U, V may be smaller
   in the case that some c_i or s_i are zero.

   The number of rows of A (and U) is the value n passed with DSSetDimensions().
   The number of columns m and the number of rows of B (and V) must be
   set via DSGSVDSetDimensions().

   Internally, LAPACK's representation is used, U'*A*Q = C*[0 R], V'*B*Q = S*[0 R],
   where X = Q*inv(R) is computed at the end of DSSolve().

   If the compact storage format is selected, then a simplified problem is
   solved, where A and B are bidiagonal (possibly with an arrow), and [A;B]
   is assumed to have orthonormal columns. We consider two cases: (1) A and B
   are square mxm upper bidiagonal, and (2) A is lower bidiagonal with m+1
   rows and B is square upper bidiagonal. In these cases, R=I so it
   corresponds to the CS decomposition. The first matrix is stored in two
   diagonals of DS_MAT_T, while the second matrix is stored in DS_MAT_D
   and the remaining diagonal of DS_MAT_T.

   Allowed arguments of DSVectors() are DS_MAT_U, DS_MAT_V and DS_MAT_X.

   Used DS matrices:
+  DS_MAT_A - first problem matrix
.  DS_MAT_B - second problem matrix
.  DS_MAT_T - first upper bidiagonal matrix (if compact storage is selected)
.  DS_MAT_D - second upper bidiagonal matrix (if compact storage is selected)
.  DS_MAT_U - (upper) left generalized singular vectors
.  DS_MAT_V - (lower) left generalized singular vectors
-  DS_MAT_X - right generalized singular vectors

   Implemented methods:
.  0 - Lapack (_ggsvd3 if available, or _ggsvd)

.seealso: DSCreate(), DSSetType(), DSType, DSGSVDSetDimensions()
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_GSVD(DS ds)
{
  DS_GSVD        *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  ds->data = (void*)ctx;

  ds->ops->allocate      = DSAllocate_GSVD;
  ds->ops->view          = DSView_GSVD;
  ds->ops->vectors       = DSVectors_GSVD;
  ds->ops->sort          = DSSort_GSVD;
  ds->ops->solve[0]      = DSSolve_GSVD;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize   = DSSynchronize_GSVD;
#endif
  ds->ops->truncate      = DSTruncate_GSVD;
  ds->ops->update        = DSUpdateExtraRow_GSVD;
  ds->ops->cond          = DSCond_GSVD;
  ds->ops->matgetsize    = DSMatGetSize_GSVD;
  ds->ops->destroy       = DSDestroy_GSVD;
  ds->ops->reallocate    = DSReallocate_GSVD;
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSGSVDSetDimensions_C",DSGSVDSetDimensions_GSVD));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSGSVDGetDimensions_C",DSGSVDGetDimensions_GSVD));
  PetscFunctionReturn(PETSC_SUCCESS);
}
