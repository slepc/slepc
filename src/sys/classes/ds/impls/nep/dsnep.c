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
  PetscInt       nf;                 /* number of functions in f[] */
  FN             f[DS_NUM_EXTRA];    /* functions defining the nonlinear operator */
  PetscInt       max_mid;            /* maximum minimality index */
  PetscInt       nnod;               /* number of nodes for quadrature rules */
  PetscInt       spls;               /* number of sampling columns for quadrature rules */
  PetscInt       Nit;                /* number of refinement iterations */
  PetscReal      rtol;               /* tolerance of Newton refinement */
  RG             rg;                 /* region for contour integral */
  PetscLayout    map;                /* used to distribute work among MPI processes */
  void           *computematrixctx;
  DSNEPMatrixFunctionFn *computematrix;
} DS_NEP;

/*
   DSNEPComputeMatrix - Build the matrix associated with a nonlinear operator
   T(lambda) or its derivative T'(lambda), given the parameter lambda, where
   T(lambda) = sum_i E_i*f_i(lambda). The result is written in mat.
*/
static PetscErrorCode DSNEPComputeMatrix(DS ds,PetscScalar lambda,PetscBool deriv,DSMatType mat)
{
  DS_NEP            *ctx = (DS_NEP*)ds->data;
  PetscScalar       *T,alpha;
  const PetscScalar *E;
  PetscInt          i,ld,n;
  PetscBLASInt      k,inc=1;

  PetscFunctionBegin;
  PetscCall(PetscLogEventBegin(DS_Other,ds,0,0,0));
  if (ctx->computematrix) PetscCall((*ctx->computematrix)(ds,lambda,deriv,mat,ctx->computematrixctx));
  else {
    PetscCall(DSGetDimensions(ds,&n,NULL,NULL,NULL));
    PetscCall(DSGetLeadingDimension(ds,&ld));
    PetscCall(PetscBLASIntCast(ld*n,&k));
    PetscCall(MatDenseGetArray(ds->omat[mat],&T));
    PetscCall(PetscArrayzero(T,k));
    for (i=0;i<ctx->nf;i++) {
      if (deriv) PetscCall(FNEvaluateDerivative(ctx->f[i],lambda,&alpha));
      else PetscCall(FNEvaluateFunction(ctx->f[i],lambda,&alpha));
      PetscCall(MatDenseGetArrayRead(ds->omat[DSMatExtra[i]],&E));
      PetscCallBLAS("BLASaxpy",BLASaxpy_(&k,&alpha,E,&inc,T,&inc));
      PetscCall(MatDenseRestoreArrayRead(ds->omat[DSMatExtra[i]],&E));
    }
    PetscCall(MatDenseRestoreArray(ds->omat[mat],&T));
  }
  PetscCall(PetscLogEventEnd(DS_Other,ds,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSAllocate_NEP(DS ds,PetscInt ld)
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_X));
  for (i=0;i<ctx->nf;i++) PetscCall(DSAllocateMat_Private(ds,DSMatExtra[i]));
  PetscCall(PetscFree(ds->perm));
  PetscCall(PetscMalloc1(ld*ctx->max_mid,&ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSView_NEP(DS ds,PetscViewer viewer)
{
  DS_NEP            *ctx = (DS_NEP*)ds->data;
  PetscViewerFormat format;
  PetscInt          i;
  const char        *methodname[] = {
                     "Successive Linear Problems",
                     "Contour Integral"
  };
  const int         nmeth=PETSC_STATIC_ARRAY_LENGTH(methodname);

  PetscFunctionBegin;
  PetscCall(PetscViewerGetFormat(viewer,&format));
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
    if (ds->method<nmeth) PetscCall(PetscViewerASCIIPrintf(viewer,"solving the problem with: %s\n",methodname[ds->method]));
#if defined(PETSC_USE_COMPLEX)
    if (ds->method==1) {  /* contour integral method */
      PetscCall(PetscViewerASCIIPrintf(viewer,"number of integration points: %" PetscInt_FMT "\n",ctx->nnod));
      PetscCall(PetscViewerASCIIPrintf(viewer,"maximum minimality index: %" PetscInt_FMT "\n",ctx->max_mid));
      if (ctx->spls) PetscCall(PetscViewerASCIIPrintf(viewer,"number of sampling columns for quadrature: %" PetscInt_FMT "\n",ctx->spls));
      if (ctx->Nit) PetscCall(PetscViewerASCIIPrintf(viewer,"doing iterative refinement (%" PetscInt_FMT " its, tolerance %g)\n",ctx->Nit,(double)ctx->rtol));
      PetscCall(RGView(ctx->rg,viewer));
    }
#endif
    if (format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscCall(PetscViewerASCIIPrintf(viewer,"number of functions: %" PetscInt_FMT "\n",ctx->nf));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  for (i=0;i<ctx->nf;i++) {
    PetscCall(FNView(ctx->f[i],viewer));
    PetscCall(DSViewMat(ds,viewer,DSMatExtra[i]));
  }
  if (ds->state>DS_STATE_INTERMEDIATE) PetscCall(DSViewMat(ds,viewer,DS_MAT_X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSVectors_NEP(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscFunctionBegin;
  PetscCheck(!rnorm,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
  switch (mat) {
    case DS_MAT_X:
      break;
    case DS_MAT_Y:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSort_NEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *dummy)
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;
  PetscInt       n,l,i,*perm,lds;
  PetscScalar    *Q;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(PETSC_SUCCESS);
  if (!ds->method) PetscFunctionReturn(PETSC_SUCCESS);  /* SLP computes just one eigenvalue */
  n = ds->n*ctx->max_mid;
  lds = ds->ld*ctx->max_mid;
  l = ds->l;
  perm = ds->perm;
  for (i=0;i<n;i++) perm[i] = i;
  if (rr) PetscCall(DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_FALSE));
  else PetscCall(DSSortEigenvalues_Private(ds,wr,NULL,perm,PETSC_FALSE));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  for (i=l;i<ds->t;i++) Q[i+i*lds] = wr[perm[i]];
  for (i=l;i<ds->t;i++) wr[i] = Q[i+i*lds];
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  /* n != ds->n */
  PetscCall(DSPermuteColumns_Private(ds,0,ds->t,ds->n,DS_MAT_X,perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(SLEPC_MISSING_LAPACK_GGEV3)
#define LAPGEEV "ggev"
#else
#define LAPGEEV "ggev3"
#endif

static PetscErrorCode DSSolve_NEP_SLP(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscScalar    *A,*B,*W,*X,*work,*alpha,*beta,a;
  PetscScalar    sigma,lambda,mu,re,re2,sone=1.0,szero=0.0;
  PetscBLASInt   info,n,ld,lwork,one=1,zero=0;
  PetscInt       it,pos,j,maxit=100,result;
  PetscReal      norm,tol,done=1.0;
#if !defined(PETSC_USE_COMPLEX)
  PetscReal      *alphai,im,im2;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_A));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_B));
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_W));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_W],&W));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));

  /* workspace query and memory allocation */
  lwork = -1;
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACK" LAPGEEV,LAPACKggevalt_("N","V",&n,A,&ld,B,&ld,NULL,NULL,NULL,&ld,W,&ld,&a,&lwork,NULL,&info));
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork));
  PetscCall(DSAllocateWork_Private(ds,lwork+2*ds->n,8*ds->n,0));
  alpha = ds->work;
  beta  = ds->work + ds->n;
  work  = ds->work + 2*ds->n;
#else
  PetscCallBLAS("LAPACK" LAPGEEV,LAPACKggevalt_("N","V",&n,A,&ld,B,&ld,NULL,NULL,NULL,NULL,&ld,W,&ld,&a,&lwork,&info));
  PetscCall(PetscBLASIntCast((PetscInt)a,&lwork));
  PetscCall(DSAllocateWork_Private(ds,lwork+3*ds->n,0,0));
  alpha  = ds->work;
  beta   = ds->work + ds->n;
  alphai = ds->work + 2*ds->n;
  work   = ds->work + 3*ds->n;
#endif

  sigma = 0.0;
  if (ds->sc->comparison==SlepcCompareTargetMagnitude || ds->sc->comparison==SlepcCompareTargetReal) sigma = *(PetscScalar*)ds->sc->comparisonctx;
  lambda = sigma;
  tol = n*PETSC_MACHINE_EPSILON/PetscSqrtReal(PETSC_SQRT_MACHINE_EPSILON);

  for (it=0;it<maxit;it++) {

    /* evaluate T and T' */
    PetscCall(DSNEPComputeMatrix(ds,lambda,PETSC_FALSE,DS_MAT_A));
    if (it) {
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,A,&ld,X,&one,&szero,X+ld,&one));
      norm = BLASnrm2_(&n,X+ld,&one);
      if (norm/PetscAbsScalar(lambda)<=tol) break;
    }
    PetscCall(DSNEPComputeMatrix(ds,lambda,PETSC_TRUE,DS_MAT_B));

    /* compute eigenvalue correction mu and eigenvector u */
#if defined(PETSC_USE_COMPLEX)
    PetscCallBLAS("LAPACK" LAPGEEV,LAPACKggevalt_("N","V",&n,A,&ld,B,&ld,alpha,beta,NULL,&ld,W,&ld,work,&lwork,ds->rwork,&info));
#else
    PetscCallBLAS("LAPACK" LAPGEEV,LAPACKggevalt_("N","V",&n,A,&ld,B,&ld,alpha,alphai,beta,NULL,&ld,W,&ld,work,&lwork,&info));
#endif
    SlepcCheckLapackInfo(LAPGEEV,info);

    /* find smallest eigenvalue */
    j = 0;
    if (beta[j]==0.0) re = (PetscRealPart(alpha[j])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else re = alpha[j]/beta[j];
#if !defined(PETSC_USE_COMPLEX)
    if (beta[j]==0.0) im = (alphai[j]>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else im = alphai[j]/beta[j];
#endif
    pos = 0;
    for (j=1;j<n;j++) {
      if (beta[j]==0.0) re2 = (PetscRealPart(alpha[j])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
      else re2 = alpha[j]/beta[j];
#if !defined(PETSC_USE_COMPLEX)
      if (beta[j]==0.0) im2 = (alphai[j]>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
      else im2 = alphai[j]/beta[j];
      PetscCall(SlepcCompareSmallestMagnitude(re,im,re2,im2,&result,NULL));
#else
      PetscCall(SlepcCompareSmallestMagnitude(re,0.0,re2,0.0,&result,NULL));
#endif
      if (result > 0) {
        re = re2;
#if !defined(PETSC_USE_COMPLEX)
        im = im2;
#endif
        pos = j;
      }
    }

#if !defined(PETSC_USE_COMPLEX)
    PetscCheck(im==0.0,PETSC_COMM_SELF,PETSC_ERR_SUP,"DSNEP found a complex eigenvalue; try rerunning with complex scalars");
#endif
    mu = alpha[pos]/beta[pos];
    PetscCall(PetscArraycpy(X,W+pos*ld,n));
    norm = BLASnrm2_(&n,X,&one);
    PetscCallBLAS("LAPACKlascl",LAPACKlascl_("G",&zero,&zero,&norm,&done,&n,&one,X,&n,&info));
    SlepcCheckLapackInfo("lascl",info);

    /* correct eigenvalue approximation */
    lambda = lambda - mu;
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_A],&A));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_B],&B));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_W],&W));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));

  PetscCheck(it<maxit,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"DSNEP did not converge");
  ds->t = 1;
  wr[0] = lambda;
  if (wi) wi[0] = 0.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_COMPLEX)
/*
  Newton refinement for eigenpairs computed with contour integral.
  k  - number of eigenpairs to refine
  wr - eigenvalues (eigenvectors are stored in DS_MAT_X)
*/
static PetscErrorCode DSNEPNewtonRefine(DS ds,PetscInt k,PetscScalar *wr)
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;
  PetscScalar    *X,*W,*U,*R,sone=1.0,szero=0.0;
  PetscReal      norm;
  PetscInt       i,j,ii,nwu=0,*p,jstart=0,jend=k;
  const PetscInt *range;
  PetscBLASInt   n,*perm,info,ld,one=1,n1;
  PetscMPIInt    len,size,root;
  PetscLayout    map;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  n1 = n+1;
  p  = ds->perm;
  PetscCall(PetscArrayzero(p,k));
  PetscCall(DSAllocateWork_Private(ds,(n+2)*(n+1),0,n+1));
  U    = ds->work+nwu;    nwu += (n+1)*(n+1);
  R    = ds->work+nwu;    /*nwu += n+1;*/
  perm = ds->iwork;
  if (ds->pmode==DS_PARALLEL_DISTRIBUTED) {
    PetscCall(PetscLayoutCreateFromSizes(PetscObjectComm((PetscObject)ds),PETSC_DECIDE,k,1,&map));
    PetscCall(PetscLayoutGetRange(map,&jstart,&jend));
  }
  for (ii=0;ii<ctx->Nit;ii++) {
    for (j=jstart;j<jend;j++) {
      if (p[j]<2) {
        PetscCall(DSNEPComputeMatrix(ds,wr[j],PETSC_FALSE,DS_MAT_W));
        PetscCall(MatDenseGetArray(ds->omat[DS_MAT_W],&W));
        PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,W,&ld,X+ld*j,&one,&szero,R,&one));
        PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_W],&W));
        norm = BLASnrm2_(&n,R,&one);
        if (norm/PetscAbsScalar(wr[j]) > ctx->rtol) {
          PetscCall(PetscInfo(NULL,"Refining eigenpair %" PetscInt_FMT ", residual=%g\n",j,(double)(norm/PetscAbsScalar(wr[j]))));
          p[j] = 1;
          R[n] = 0.0;
          PetscCall(MatDenseGetArray(ds->omat[DS_MAT_W],&W));
          for (i=0;i<n;i++) {
            PetscCall(PetscArraycpy(U+i*n1,W+i*ld,n));
            U[n+i*n1] = PetscConj(X[j*ld+i]);
          }
          PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_W],&W));
          U[n+n*n1] = 0.0;
          PetscCall(DSNEPComputeMatrix(ds,wr[j],PETSC_TRUE,DS_MAT_W));
          PetscCall(MatDenseGetArray(ds->omat[DS_MAT_W],&W));
          PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,W,&ld,X+ld*j,&one,&szero,U+n*(n+1),&one));
          PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_W],&W));
          /* solve system  */
          PetscCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n1,&n1,U,&n1,perm,&info));
          SlepcCheckLapackInfo("getrf",info);
          PetscCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&n1,&one,U,&n1,perm,R,&n1,&info));
          SlepcCheckLapackInfo("getrs",info);
          wr[j] -= R[n];
          for (i=0;i<n;i++) X[j*ld+i] -= R[i];
          /* normalization */
          norm = BLASnrm2_(&n,X+ld*j,&one);
          for (i=0;i<n;i++) X[ld*j+i] /= norm;
        } else p[j] = 2;
      }
    }
  }
  if (ds->pmode==DS_PARALLEL_DISTRIBUTED) {  /* communicate results */
    PetscCall(PetscMPIIntCast(k,&len));
    PetscCallMPI(MPIU_Allreduce(MPI_IN_PLACE,p,len,MPIU_INT,MPI_SUM,PetscObjectComm((PetscObject)ds)));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ds),&size));
    PetscCall(PetscLayoutGetRanges(map,&range));
    for (j=0;j<k;j++) {
      if (p[j]) {  /* j-th eigenpair has been refined */
        for (root=0;root<size;root++) if (range[root+1]>j) break;
        PetscCall(PetscMPIIntCast(1,&len));
        PetscCallMPI(MPI_Bcast(wr+j,len,MPIU_SCALAR,root,PetscObjectComm((PetscObject)ds)));
        PetscCall(PetscMPIIntCast(n,&len));
        PetscCallMPI(MPI_Bcast(X+ld*j,len,MPIU_SCALAR,root,PetscObjectComm((PetscObject)ds)));
      }
    }
    PetscCall(PetscLayoutDestroy(&map));
  }
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DSSolve_NEP_Contour(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;
  PetscScalar    *alpha,*beta,*Q,*Z,*X,*U,*V,*W,*work,*Rc,*R,*w,*z,*zn,*S;
  PetscScalar    sone=1.0,szero=0.0,center,a;
  PetscReal      *rwork,norm,radius,vscale,rgscale,*sigma;
  PetscBLASInt   info,n,*perm,p,pp,ld,lwork,k_,rk_,colA,rowA,one=1;
  PetscInt       mid,lds,nnod=ctx->nnod,k,i,ii,jj,j,s,off,rk,nwu=0,nw,lrwork,*inside,kstart=0,kend=nnod;
  PetscMPIInt    len;
  PetscBool      isellipse;
  PetscRandom    rand;

  PetscFunctionBegin;
  PetscCheck(ctx->rg,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"The contour solver requires a region passed with DSNEPSetRG()");
  /* Contour parameters */
  PetscCall(PetscObjectTypeCompare((PetscObject)ctx->rg,RGELLIPSE,&isellipse));
  PetscCheck(isellipse,PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Region must be Ellipse");
  PetscCall(RGEllipseGetParameters(ctx->rg,&center,&radius,&vscale));
  PetscCall(RGGetScale(ctx->rg,&rgscale));
  if (ds->pmode==DS_PARALLEL_DISTRIBUTED) {
    if (!ctx->map) PetscCall(PetscLayoutCreateFromSizes(PetscObjectComm((PetscObject)ds),PETSC_DECIDE,ctx->nnod,1,&ctx->map));
    PetscCall(PetscLayoutGetRange(ctx->map,&kstart,&kend));
  }

  PetscCall(DSAllocateMat_Private(ds,DS_MAT_W)); /* size n */
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Q)); /* size mid*n */
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_Z)); /* size mid*n */
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_U)); /* size mid*n */
  PetscCall(DSAllocateMat_Private(ds,DS_MAT_V)); /* size mid*n */
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Z],&Z));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_V],&V));
  mid  = ctx->max_mid;
  PetscCall(PetscBLASIntCast(ds->n,&n));
  p    = n;   /* maximum number of columns for the probing matrix */
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(mid*n,&rowA));
  nw     = 2*n*(p+mid)+3*nnod+2*mid*n*p;
  lrwork = 9*mid*n;

  /* workspace query and memory allocation */
  lwork = -1;
  PetscCallBLAS("LAPACK" LAPGEEV,LAPACKggevalt_("N","V",&rowA,Q,&rowA,Z,&rowA,NULL,NULL,NULL,&ld,V,&rowA,&a,&lwork,NULL,&info));
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork));
  PetscCall(DSAllocateWork_Private(ds,lwork+nw,lrwork,n+1));

  sigma = ds->rwork;
  rwork = ds->rwork+mid*n;
  perm  = ds->iwork;
  z     = ds->work+nwu;    nwu += nnod;         /* quadrature points */
  zn    = ds->work+nwu;    nwu += nnod;         /* normalized quadrature points */
  w     = ds->work+nwu;    nwu += nnod;         /* quadrature weights */
  Rc    = ds->work+nwu;    nwu += n*p;
  R     = ds->work+nwu;    nwu += n*p;
  alpha = ds->work+nwu;    nwu += mid*n;
  beta  = ds->work+nwu;    nwu += mid*n;
  S     = ds->work+nwu;    nwu += 2*mid*n*p;
  work  = ds->work+nwu;

  /* Compute quadrature parameters */
  PetscCall(RGComputeQuadrature(ctx->rg,RG_QUADRULE_TRAPEZOIDAL,nnod,z,zn,w));

  /* Set random matrix */
  PetscCall(PetscRandomCreate(PetscObjectComm((PetscObject)ds),&rand));
  PetscCall(PetscRandomSetSeed(rand,0x12345678));
  PetscCall(PetscRandomSeed(rand));
  for (j=0;j<p;j++)
    for (i=0;i<n;i++) PetscCall(PetscRandomGetValue(rand,Rc+i+j*n));
  PetscCall(PetscArrayzero(S,2*mid*n*p));
  /* Loop of integration points */
  for (k=kstart;k<kend;k++) {
    PetscCall(PetscInfo(NULL,"Solving integration point %" PetscInt_FMT "\n",k));
    PetscCall(PetscArraycpy(R,Rc,p*n));
    PetscCall(DSNEPComputeMatrix(ds,z[k],PETSC_FALSE,DS_MAT_W));

    /* LU factorization */
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_W],&W));
    PetscCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n,&n,W,&ld,perm,&info));
    SlepcCheckLapackInfo("getrf",info);
    PetscCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&n,&p,W,&ld,perm,R,&n,&info));
    SlepcCheckLapackInfo("getrs",info);
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_W],&W));

    /* Moments computation */
    for (s=0;s<2*ctx->max_mid;s++) {
      off = s*n*p;
      for (j=0;j<p;j++)
        for (i=0;i<n;i++) S[off+i+j*n] += w[k]*R[j*n+i];
      w[k] *= zn[k];
    }
  }

  if (ds->pmode==DS_PARALLEL_DISTRIBUTED) {  /* compute final S via reduction */
    PetscCall(PetscMPIIntCast(2*mid*n*p,&len));
    PetscCallMPI(MPIU_Allreduce(MPI_IN_PLACE,S,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)ds)));
  }
  PetscCall(PetscBLASIntCast(ctx->spls?PetscMin(ctx->spls,n):n,&p));
  pp = p;
  do {
    p = pp;
    PetscCall(PetscBLASIntCast(mid*p,&colA));

    PetscCall(PetscInfo(ds,"Computing SVD of size %" PetscBLASInt_FMT "x%" PetscBLASInt_FMT "\n",rowA,colA));
    for (jj=0;jj<mid;jj++) {
      for (ii=0;ii<mid;ii++) {
        off = jj*p*rowA+ii*n;
        for (j=0;j<p;j++)
          for (i=0;i<n;i++) Q[off+j*rowA+i] = S[((jj+ii)*n+j)*n+i];
      }
    }
    PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&rowA,&colA,Q,&rowA,sigma,U,&rowA,V,&colA,work,&lwork,rwork,&info));
    SlepcCheckLapackInfo("gesvd",info);

    rk = colA;
    for (i=1;i<colA;i++) if (sigma[i]/sigma[0]<PETSC_MACHINE_EPSILON*1e4) {rk = i; break;}
    if (rk<colA || p==n) break;
    pp *= 2;
  } while (pp<=n);
  PetscCall(PetscInfo(ds,"Solving generalized eigenproblem of size %" PetscInt_FMT "\n",rk));
  for (jj=0;jj<mid;jj++) {
    for (ii=0;ii<mid;ii++) {
      off = jj*p*rowA+ii*n;
      for (j=0;j<p;j++)
        for (i=0;i<n;i++) Q[off+j*rowA+i] = S[((jj+ii+1)*n+j)*n+i];
    }
  }
  PetscCall(PetscBLASIntCast(rk,&rk_));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&rowA,&rk_,&colA,&sone,Q,&rowA,V,&colA,&szero,Z,&rowA));
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&rk_,&rk_,&rowA,&sone,U,&rowA,Z,&rowA,&szero,Q,&rk_));
  PetscCall(PetscArrayzero(Z,n*mid*n*mid));
  for (j=0;j<rk;j++) Z[j+j*rk_] = sigma[j];
  PetscCallBLAS("LAPACK" LAPGEEV,LAPACKggevalt_("N","V",&rk_,Q,&rk_,Z,&rk_,alpha,beta,NULL,&ld,V,&rk_,work,&lwork,rwork,&info));
  SlepcCheckLapackInfo(LAPGEEV,info);
  for (i=0;i<rk;i++) wr[i] = (center+radius*alpha[i]/beta[i])*rgscale;
  PetscCall(PetscMalloc1(rk,&inside));
  PetscCall(RGCheckInside(ctx->rg,rk,wr,wi,inside));
  k=0;
  for (i=0;i<rk;i++)
    if (inside[i]==1) inside[k++] = i;
  /* Discard values outside region */
  lds = ld*mid;
  PetscCall(PetscArrayzero(Q,lds*lds));
  PetscCall(PetscArrayzero(Z,lds*lds));
  for (i=0;i<k;i++) Q[i+i*lds] = (center*beta[inside[i]]+radius*alpha[inside[i]])*rgscale;
  for (i=0;i<k;i++) Z[i+i*lds] = beta[inside[i]];
  for (i=0;i<k;i++) wr[i] = Q[i+i*lds]/Z[i+i*lds];
  for (j=0;j<k;j++) for (i=0;i<rk;i++) V[j*rk+i] = sigma[i]*V[inside[j]*rk+i];
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&k_,&rk_,&sone,U,&rowA,V,&rk_,&szero,X,&ld));
  /* Normalize */
  for (j=0;j<k;j++) {
    norm = BLASnrm2_(&n,X+ld*j,&one);
    for (i=0;i<n;i++) X[ld*j+i] /= norm;
  }
  PetscCall(PetscFree(inside));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Z],&Z));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_U],&U));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_V],&V));

  /* Newton refinement */
  if (ctx->Nit) PetscCall(DSNEPNewtonRefine(ds,k,wr));
  ds->t = k;
  PetscCall(PetscRandomDestroy(&rand));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

#if !defined(PETSC_HAVE_MPIUNI)
static PetscErrorCode DSSynchronize_NEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;
  PetscInt       ld=ds->ld,k=0;
  PetscMPIInt    n,n2,rank,size,off=0;
  PetscScalar    *X;

  PetscFunctionBegin;
  if (!ds->method) { /* SLP */
    if (ds->state>=DS_STATE_CONDENSED) k += ds->n;
    if (eigr) k += 1;
    if (eigi) k += 1;
    PetscCall(PetscMPIIntCast(1,&n));
    PetscCall(PetscMPIIntCast(ds->n,&n2));
  } else { /* Contour */
    if (ds->state>=DS_STATE_CONDENSED) k += ctx->max_mid*ds->n*ld;
    if (eigr) k += ctx->max_mid*ds->n;
    if (eigi) k += ctx->max_mid*ds->n;
    PetscCall(PetscMPIIntCast(ctx->max_mid*ds->n,&n));
    PetscCall(PetscMPIIntCast(ctx->max_mid*ds->n*ld,&n2));
  }
  PetscCall(DSAllocateWork_Private(ds,k,0,0));
  PetscCall(PetscMPIIntCast(k*sizeof(PetscScalar),&size));
  if (ds->state>=DS_STATE_CONDENSED) PetscCall(MatDenseGetArray(ds->omat[DS_MAT_X],&X));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank));
  if (!rank) {
    if (ds->state>=DS_STATE_CONDENSED) PetscCallMPI(MPI_Pack(X,n2,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
    if (eigr) PetscCallMPI(MPI_Pack(eigr,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
#if !defined(PETSC_USE_COMPLEX)
    if (eigi) PetscCallMPI(MPI_Pack(eigi,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds)));
#endif
  }
  PetscCallMPI(MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds)));
  if (rank) {
    if (ds->state>=DS_STATE_CONDENSED) PetscCallMPI(MPI_Unpack(ds->work,size,&off,X,n2,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
    if (eigr) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigr,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
#if !defined(PETSC_USE_COMPLEX)
    if (eigi) PetscCallMPI(MPI_Unpack(ds->work,size,&off,eigi,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds)));
#endif
  }
  if (ds->state>=DS_STATE_CONDENSED) PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_X],&X));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode DSNEPSetFN_NEP(DS ds,PetscInt n,FN fn[])
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCheck(n>0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Must have one or more functions, you have %" PetscInt_FMT,n);
  PetscCheck(n<=DS_NUM_EXTRA,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Too many functions, you specified %" PetscInt_FMT " but the limit is %d",n,DS_NUM_EXTRA);
  if (ds->ld) PetscCall(PetscInfo(ds,"DSNEPSetFN() called after DSAllocate()\n"));
  for (i=0;i<n;i++) PetscCall(PetscObjectReference((PetscObject)fn[i]));
  for (i=0;i<ctx->nf;i++) PetscCall(FNDestroy(&ctx->f[i]));
  for (i=0;i<n;i++) ctx->f[i] = fn[i];
  ctx->nf = n;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPSetFN - Sets a number of functions that define the nonlinear
   eigenproblem.

   Collective

   Input Parameters:
+  ds - the direct solver context
.  n  - number of functions
-  fn - array of functions

   Notes:
   The nonlinear eigenproblem is defined in terms of the split nonlinear
   operator T(lambda) = sum_i A_i*f_i(lambda).

   This function must be called before DSAllocate(). Then DSAllocate()
   will allocate an extra matrix A_i per each function, that can be
   filled in the usual way.

   Level: advanced

.seealso: DSNEPGetFN(), DSAllocate()
@*/
PetscErrorCode DSNEPSetFN(DS ds,PetscInt n,FN fn[])
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,n,2);
  PetscAssertPointer(fn,3);
  for (i=0;i<n;i++) {
    PetscValidHeaderSpecific(fn[i],FN_CLASSID,3);
    PetscCheckSameComm(ds,1,fn[i],3);
  }
  PetscTryMethod(ds,"DSNEPSetFN_C",(DS,PetscInt,FN[]),(ds,n,fn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPGetFN_NEP(DS ds,PetscInt k,FN *fn)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  PetscCheck(k>=0 && k<ctx->nf,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,ctx->nf-1);
  *fn = ctx->f[k];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPGetFN - Gets the functions associated with the nonlinear DS.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
-  k  - the index of the requested function (starting in 0)

   Output Parameter:
.  fn - the function

   Level: advanced

.seealso: DSNEPSetFN()
@*/
PetscErrorCode DSNEPGetFN(DS ds,PetscInt k,FN *fn)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(fn,3);
  PetscUseMethod(ds,"DSNEPGetFN_C",(DS,PetscInt,FN*),(ds,k,fn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPGetNumFN_NEP(DS ds,PetscInt *n)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  *n = ctx->nf;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPGetNumFN - Returns the number of functions stored internally by
   the DS.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  n - the number of functions passed in DSNEPSetFN()

   Level: advanced

.seealso: DSNEPSetFN()
@*/
PetscErrorCode DSNEPGetNumFN(DS ds,PetscInt *n)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(n,2);
  PetscUseMethod(ds,"DSNEPGetNumFN_C",(DS,PetscInt*),(ds,n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPSetMinimality_NEP(DS ds,PetscInt n)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  if (n == PETSC_DECIDE || n == PETSC_DEFAULT) ctx->max_mid = 4;
  else {
    PetscCheck(n>0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The minimality value must be > 0");
    ctx->max_mid = n;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPSetMinimality - Sets the maximum minimality index used internally by
   the DSNEP.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  n  - the maximum minimality index

   Options Database Key:
.  -ds_nep_minimality <n> - sets the maximum minimality index

   Notes:
   The maximum minimality index is used only in the contour integral method,
   and is related to the highest momemts used in the method. The default
   value is 1, an larger value might give better accuracy in some cases, but
   at a higher cost.

   Level: advanced

.seealso: DSNEPGetMinimality()
@*/
PetscErrorCode DSNEPSetMinimality(DS ds,PetscInt n)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,n,2);
  PetscTryMethod(ds,"DSNEPSetMinimality_C",(DS,PetscInt),(ds,n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPGetMinimality_NEP(DS ds,PetscInt *n)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  *n = ctx->max_mid;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPGetMinimality - Returns the maximum minimality index used internally by
   the DSNEP.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  n - the maximum minimality index passed in DSNEPSetMinimality()

   Level: advanced

.seealso: DSNEPSetMinimality()
@*/
PetscErrorCode DSNEPGetMinimality(DS ds,PetscInt *n)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(n,2);
  PetscUseMethod(ds,"DSNEPGetMinimality_C",(DS,PetscInt*),(ds,n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPSetRefine_NEP(DS ds,PetscReal tol,PetscInt its)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  if (tol == (PetscReal)PETSC_DETERMINE) {
    ctx->rtol = PETSC_MACHINE_EPSILON/PetscSqrtReal(PETSC_SQRT_MACHINE_EPSILON);
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The tolerance must be > 0");
    ctx->rtol = tol;
  }
  if (its == PETSC_DETERMINE) {
    ctx->Nit = 3;
  } else if (its != PETSC_CURRENT) {
    PetscCheck(its>=0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The number of iterations must be >= 0");
    ctx->Nit = its;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPSetRefine - Sets the tolerance and the number of iterations of Newton iterative
   refinement for eigenpairs.

   Logically Collective

   Input Parameters:
+  ds  - the direct solver context
.  tol - the tolerance
-  its - the number of iterations

   Options Database Key:
+  -ds_nep_refine_tol <tol> - sets the tolerance
-  -ds_nep_refine_its <its> - sets the number of Newton iterations

   Notes:
   Iterative refinement of eigenpairs is currently used only in the contour
   integral method.

   Use PETSC_CURRENT to retain the current value of any of the parameters.
   Use PETSC_DETERMINE for either argument to assign a default value computed
   internally.

   Level: advanced

.seealso: DSNEPGetRefine()
@*/
PetscErrorCode DSNEPSetRefine(DS ds,PetscReal tol,PetscInt its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveReal(ds,tol,2);
  PetscValidLogicalCollectiveInt(ds,its,3);
  PetscTryMethod(ds,"DSNEPSetRefine_C",(DS,PetscReal,PetscInt),(ds,tol,its));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPGetRefine_NEP(DS ds,PetscReal *tol,PetscInt *its)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  if (tol) *tol = ctx->rtol;
  if (its) *its = ctx->Nit;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPGetRefine - Returns the tolerance and the number of iterations of Newton iterative
   refinement for eigenpairs.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
+  tol - the tolerance
-  its - the number of iterations

   Level: advanced

.seealso: DSNEPSetRefine()
@*/
PetscErrorCode DSNEPGetRefine(DS ds,PetscReal *tol,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscUseMethod(ds,"DSNEPGetRefine_C",(DS,PetscReal*,PetscInt*),(ds,tol,its));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPSetIntegrationPoints_NEP(DS ds,PetscInt ip)
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  if (ip == PETSC_DECIDE || ip == PETSC_DEFAULT) ctx->nnod = 64;
  else {
    PetscCheck(ip>0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The number of integration points must be > 0");
    ctx->nnod = ip;
  }
  PetscCall(PetscLayoutDestroy(&ctx->map));  /* need to redistribute at next solve */
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPSetIntegrationPoints - Sets the number of integration points to be
   used in the contour integral method.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  ip - the number of integration points

   Options Database Key:
.  -ds_nep_integration_points <ip> - sets the number of integration points

   Notes:
   This parameter is relevant only in the contour integral method.

   Level: advanced

.seealso: DSNEPGetIntegrationPoints()
@*/
PetscErrorCode DSNEPSetIntegrationPoints(DS ds,PetscInt ip)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,ip,2);
  PetscTryMethod(ds,"DSNEPSetIntegrationPoints_C",(DS,PetscInt),(ds,ip));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPGetIntegrationPoints_NEP(DS ds,PetscInt *ip)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  *ip = ctx->nnod;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPGetIntegrationPoints - Returns the number of integration points used
   in the contour integral method.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  ip - the number of integration points

   Level: advanced

.seealso: DSNEPSetIntegrationPoints()
@*/
PetscErrorCode DSNEPGetIntegrationPoints(DS ds,PetscInt *ip)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(ip,2);
  PetscUseMethod(ds,"DSNEPGetIntegrationPoints_C",(DS,PetscInt*),(ds,ip));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPSetSamplingSize_NEP(DS ds,PetscInt p)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  if (p == PETSC_DECIDE || p == PETSC_DEFAULT) ctx->spls = 0;
  else {
    PetscCheck(p>0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The sample size must be > 0");
    PetscCheck(p>=20,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The sample size cannot be smaller than 20");
    ctx->spls = p;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPSetSamplingSize - Sets the number of sampling columns to be
   used in the contour integral method.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  p - the number of columns for the sampling matrix

   Options Database Key:
.  -ds_nep_sampling_size <p> - sets the number of sampling columns

   Notes:
   This parameter is relevant only in the contour integral method.

   Level: advanced

.seealso: DSNEPGetSamplingSize()
@*/
PetscErrorCode DSNEPSetSamplingSize(DS ds,PetscInt p)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,p,2);
  PetscTryMethod(ds,"DSNEPSetSamplingSize_C",(DS,PetscInt),(ds,p));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPGetSamplingSize_NEP(DS ds,PetscInt *p)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  *p = ctx->spls;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPGetSamplingSize - Returns the number of sampling columns used
   in the contour integral method.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  p -  the number of columns for the sampling matrix

   Level: advanced

.seealso: DSNEPSetSamplingSize()
@*/
PetscErrorCode DSNEPGetSamplingSize(DS ds,PetscInt *p)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(p,2);
  PetscUseMethod(ds,"DSNEPGetSamplingSize_C",(DS,PetscInt*),(ds,p));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPSetComputeMatrixFunction_NEP(DS ds,DSNEPMatrixFunctionFn *fun,void *ctx)
{
  DS_NEP *dsctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  dsctx->computematrix    = fun;
  dsctx->computematrixctx = ctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSNEPSetComputeMatrixFunction - Sets a user-provided subroutine to compute
   the matrices T(lambda) or T'(lambda).

   Logically Collective

   Input Parameters:
+  ds  - the direct solver context
.  fun - matrix function evaluation routine, see DSNEPMatrixFunctionFn for the calling sequence
-  ctx - a context pointer (the last parameter to the user function)

   Note:
   The result is computed as T(lambda) = sum_i E_i*f_i(lambda), and similarly
   for the derivative.

   Level: developer

.seealso: DSNEPGetComputeMatrixFunction()
@*/
PetscErrorCode DSNEPSetComputeMatrixFunction(DS ds,DSNEPMatrixFunctionFn *fun,void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscTryMethod(ds,"DSNEPSetComputeMatrixFunction_C",(DS,DSNEPMatrixFunctionFn*,void*),(ds,fun,ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPGetComputeMatrixFunction_NEP(DS ds,DSNEPMatrixFunctionFn **fun,void **ctx)
{
  DS_NEP *dsctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  if (fun) *fun = dsctx->computematrix;
  if (ctx) *ctx = dsctx->computematrixctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSNEPGetComputeMatrixFunction - Returns the user-provided callback function
   set in DSNEPSetComputeMatrixFunction().

   Not Collective

   Input Parameter:
.  ds  - the direct solver context

   Output Parameters:
+  fun - the pointer to the user function
-  ctx - the context pointer

   Level: developer

.seealso: DSNEPSetComputeMatrixFunction()
@*/
PetscErrorCode DSNEPGetComputeMatrixFunction(DS ds,DSNEPMatrixFunctionFn **fun,void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscUseMethod(ds,"DSNEPGetComputeMatrixFunction_C",(DS,DSNEPMatrixFunctionFn**,void**),(ds,fun,ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPSetRG_NEP(DS ds,RG rg)
{
  DS_NEP         *dsctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)rg));
  PetscCall(RGDestroy(&dsctx->rg));
  dsctx->rg = rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPSetRG - Associates a region object to the DSNEP solver.

   Collective

   Input Parameters:
+  ds  - the direct solver context
-  rg  - the region context

   Notes:
   The region is used only in the contour integral method, and
   should enclose the wanted eigenvalues.

   Level: developer

.seealso: DSNEPGetRG()
@*/
PetscErrorCode DSNEPSetRG(DS ds,RG rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  if (rg) {
    PetscValidHeaderSpecific(rg,RG_CLASSID,2);
    PetscCheckSameComm(ds,1,rg,2);
  }
  PetscTryMethod(ds,"DSNEPSetRG_C",(DS,RG),(ds,rg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSNEPGetRG_NEP(DS ds,RG *rg)
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  if (!ctx->rg) {
    PetscCall(RGCreate(PetscObjectComm((PetscObject)ds),&ctx->rg));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->rg,(PetscObject)ds,1));
    PetscCall(RGSetOptionsPrefix(ctx->rg,((PetscObject)ds)->prefix));
    PetscCall(RGAppendOptionsPrefix(ctx->rg,"ds_nep_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->rg,((PetscObject)ds)->options));
  }
  *rg = ctx->rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSNEPGetRG - Obtain the region object associated to the DSNEP solver.

   Collective

   Input Parameter:
.  ds  - the direct solver context

   Output Parameter:
.  rg  - the region context

   Level: developer

.seealso: DSNEPSetRG()
@*/
PetscErrorCode DSNEPGetRG(DS ds,RG *rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(rg,2);
  PetscUseMethod(ds,"DSNEPGetRG_C",(DS,RG*),(ds,rg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSSetFromOptions_NEP(DS ds,PetscOptionItems *PetscOptionsObject)
{
  PetscInt       k;
  PetscBool      flg;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      r;
  PetscBool      flg1;
  DS_NEP         *ctx = (DS_NEP*)ds->data;
#endif

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"DS NEP Options");

    PetscCall(PetscOptionsInt("-ds_nep_minimality","Maximum minimality index","DSNEPSetMinimality",4,&k,&flg));
    if (flg) PetscCall(DSNEPSetMinimality(ds,k));

    PetscCall(PetscOptionsInt("-ds_nep_integration_points","Number of integration points","DSNEPSetIntegrationPoints",64,&k,&flg));
    if (flg) PetscCall(DSNEPSetIntegrationPoints(ds,k));

    PetscCall(PetscOptionsInt("-ds_nep_sampling_size","Number of sampling columns","DSNEPSetSamplingSize",0,&k,&flg));
    if (flg) PetscCall(DSNEPSetSamplingSize(ds,k));

#if defined(PETSC_USE_COMPLEX)
    r = ctx->rtol;
    PetscCall(PetscOptionsReal("-ds_nep_refine_tol","Refinement tolerance","DSNEPSetRefine",ctx->rtol,&r,&flg1));
    k = ctx->Nit;
    PetscCall(PetscOptionsInt("-ds_nep_refine_its","Number of iterative refinement iterations","DSNEPSetRefine",ctx->Nit,&k,&flg));
    if (flg1||flg) PetscCall(DSNEPSetRefine(ds,r,k));

    if (ds->method==1) {
      if (!ctx->rg) PetscCall(DSNEPGetRG(ds,&ctx->rg));
      PetscCall(RGSetFromOptions(ctx->rg));
    }
#endif

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSDestroy_NEP(DS ds)
{
  DS_NEP         *ctx = (DS_NEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  for (i=0;i<ctx->nf;i++) PetscCall(FNDestroy(&ctx->f[i]));
  PetscCall(RGDestroy(&ctx->rg));
  PetscCall(PetscLayoutDestroy(&ctx->map));
  PetscCall(PetscFree(ds->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetFN_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetFN_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetNumFN_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetMinimality_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetMinimality_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetRefine_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetRefine_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetIntegrationPoints_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetIntegrationPoints_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetSamplingSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetSamplingSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetRG_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetRG_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetComputeMatrixFunction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetComputeMatrixFunction_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode DSMatGetSize_NEP(DS ds,DSMatType t,PetscInt *rows,PetscInt *cols)
{
  DS_NEP *ctx = (DS_NEP*)ds->data;

  PetscFunctionBegin;
  *rows = ds->n;
  if (t==DS_MAT_Q || t==DS_MAT_Z || t==DS_MAT_U || t==DS_MAT_V) *rows *= ctx->max_mid;
  *cols = ds->n;
  if (t==DS_MAT_Q || t==DS_MAT_Z || t==DS_MAT_U || t==DS_MAT_V || t==DS_MAT_X || t==DS_MAT_Y) *cols *= ctx->max_mid;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   DSNEP - Dense Nonlinear Eigenvalue Problem.

   Level: beginner

   Notes:
   The problem is expressed as T(lambda)*x = 0, where T(lambda) is a
   parameter-dependent matrix written as T(lambda) = sum_i E_i*f_i(lambda).
   The eigenvalues lambda are the arguments returned by DSSolve()..

   The coefficient matrices E_i are the extra matrices of the DS, and
   the scalar functions f_i are passed via DSNEPSetFN(). Optionally, a
   callback function to fill the E_i matrices can be set with
   DSNEPSetComputeMatrixFunction().

   Used DS matrices:
+  DS_MAT_Ex - coefficient matrices of the split form of T(lambda)
.  DS_MAT_X  - eigenvectors
.  DS_MAT_A  - (workspace) T(lambda) evaluated at a given lambda (SLP only)
.  DS_MAT_B  - (workspace) T'(lambda) evaluated at a given lambda (SLP only)
.  DS_MAT_Q  - (workspace) left Hankel matrix (contour only)
.  DS_MAT_Z  - (workspace) right Hankel matrix (contour only)
.  DS_MAT_U  - (workspace) left singular vectors (contour only)
.  DS_MAT_V  - (workspace) right singular vectors (contour only)
-  DS_MAT_W  - (workspace) auxiliary matrix of size nxn

   Implemented methods:
+  0 - Successive Linear Problems (SLP), computes just one eigenpair
-  1 - Contour integral, computes all eigenvalues inside a region

.seealso: DSCreate(), DSSetType(), DSType, DSNEPSetFN(), DSNEPSetComputeMatrixFunction()
M*/
SLEPC_EXTERN PetscErrorCode DSCreate_NEP(DS ds)
{
  DS_NEP         *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  ds->data = (void*)ctx;
  ctx->max_mid = 4;
  ctx->nnod    = 64;
  ctx->Nit     = 3;
  ctx->rtol    = PETSC_MACHINE_EPSILON/PetscSqrtReal(PETSC_SQRT_MACHINE_EPSILON);

  ds->ops->allocate       = DSAllocate_NEP;
  ds->ops->setfromoptions = DSSetFromOptions_NEP;
  ds->ops->view           = DSView_NEP;
  ds->ops->vectors        = DSVectors_NEP;
  ds->ops->solve[0]       = DSSolve_NEP_SLP;
#if defined(PETSC_USE_COMPLEX)
  ds->ops->solve[1]       = DSSolve_NEP_Contour;
#endif
  ds->ops->sort           = DSSort_NEP;
#if !defined(PETSC_HAVE_MPIUNI)
  ds->ops->synchronize    = DSSynchronize_NEP;
#endif
  ds->ops->destroy        = DSDestroy_NEP;
  ds->ops->matgetsize     = DSMatGetSize_NEP;

  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetFN_C",DSNEPSetFN_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetFN_C",DSNEPGetFN_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetNumFN_C",DSNEPGetNumFN_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetMinimality_C",DSNEPGetMinimality_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetMinimality_C",DSNEPSetMinimality_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetRefine_C",DSNEPGetRefine_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetRefine_C",DSNEPSetRefine_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetIntegrationPoints_C",DSNEPGetIntegrationPoints_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetIntegrationPoints_C",DSNEPSetIntegrationPoints_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetSamplingSize_C",DSNEPGetSamplingSize_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetSamplingSize_C",DSNEPSetSamplingSize_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetRG_C",DSNEPSetRG_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetRG_C",DSNEPGetRG_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPSetComputeMatrixFunction_C",DSNEPSetComputeMatrixFunction_NEP));
  PetscCall(PetscObjectComposeFunction((PetscObject)ds,"DSNEPGetComputeMatrixFunction_C",DSNEPGetComputeMatrixFunction_NEP));
  PetscFunctionReturn(PETSC_SUCCESS);
}
