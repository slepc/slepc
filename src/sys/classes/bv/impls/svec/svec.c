/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented as a single Vec
*/

#include <slepc/private/bvimpl.h>
#include "svec.h"

static PetscErrorCode BVMult_Svec(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  BV_SVEC           *y = (BV_SVEC*)Y->data,*x = (BV_SVEC*)X->data;
  const PetscScalar *px,*q;
  PetscScalar       *py;
  PetscInt          ldq;

  PetscFunctionBegin;
  PetscCall(VecGetArrayRead(x->v,&px));
  PetscCall(VecGetArray(y->v,&py));
  if (Q) {
    PetscCall(MatDenseGetLDA(Q,&ldq));
    PetscCall(MatDenseGetArrayRead(Q,&q));
    PetscCall(BVMult_BLAS_Private(Y,Y->n,Y->k-Y->l,X->k-X->l,alpha,px+(X->nc+X->l)*X->ld,X->ld,q+Y->l*ldq+X->l,ldq,beta,py+(Y->nc+Y->l)*Y->ld,Y->ld));
    PetscCall(MatDenseRestoreArrayRead(Q,&q));
  } else PetscCall(BVAXPY_BLAS_Private(Y,Y->n,Y->k-Y->l,alpha,px+(X->nc+X->l)*X->ld,X->ld,beta,py+(Y->nc+Y->l)*Y->ld,Y->ld));
  PetscCall(VecRestoreArrayRead(x->v,&px));
  PetscCall(VecRestoreArray(y->v,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultVec_Svec(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  BV_SVEC        *x = (BV_SVEC*)X->data;
  PetscScalar    *px,*py,*qq=q;

  PetscFunctionBegin;
  PetscCall(VecGetArray(x->v,&px));
  PetscCall(VecGetArray(y,&py));
  if (!q) PetscCall(VecGetArray(X->buffer,&qq));
  PetscCall(BVMultVec_BLAS_Private(X,X->n,X->k-X->l,alpha,px+(X->nc+X->l)*X->ld,X->ld,qq,beta,py));
  if (!q) PetscCall(VecRestoreArray(X->buffer,&qq));
  PetscCall(VecRestoreArray(x->v,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultInPlace_Svec(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_SVEC           *ctx = (BV_SVEC*)V->data;
  PetscScalar       *pv;
  const PetscScalar *q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(VecGetArray(ctx->v,&pv));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,s-V->l,e-V->l,pv+(V->nc+V->l)*V->ld,V->ld,q+V->l*ldq+V->l,ldq,PETSC_FALSE));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscCall(VecRestoreArray(ctx->v,&pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultInPlaceHermitianTranspose_Svec(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_SVEC           *ctx = (BV_SVEC*)V->data;
  PetscScalar       *pv;
  const PetscScalar *q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(VecGetArray(ctx->v,&pv));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,s-V->l,e-V->l,pv+(V->nc+V->l)*V->ld,V->ld,q+V->l*ldq+V->l,ldq,PETSC_TRUE));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscCall(VecRestoreArray(ctx->v,&pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDot_Svec(BV X,BV Y,Mat M)
{
  BV_SVEC           *x = (BV_SVEC*)X->data,*y = (BV_SVEC*)Y->data;
  const PetscScalar *px,*py;
  PetscScalar       *m;
  PetscInt          ldm;

  PetscFunctionBegin;
  PetscCall(MatDenseGetLDA(M,&ldm));
  PetscCall(VecGetArrayRead(x->v,&px));
  PetscCall(VecGetArrayRead(y->v,&py));
  PetscCall(MatDenseGetArray(M,&m));
  PetscCall(BVDot_BLAS_Private(X,Y->k-Y->l,X->k-X->l,X->n,py+(Y->nc+Y->l)*Y->ld,Y->ld,px+(X->nc+X->l)*X->ld,X->ld,m+X->l*ldm+Y->l,ldm,x->mpi));
  PetscCall(MatDenseRestoreArray(M,&m));
  PetscCall(VecRestoreArrayRead(x->v,&px));
  PetscCall(VecRestoreArrayRead(y->v,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDotVec_Svec(BV X,Vec y,PetscScalar *q)
{
  BV_SVEC           *x = (BV_SVEC*)X->data;
  const PetscScalar *px,*py;
  PetscScalar       *qq=q;
  Vec               z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(VecGetArrayRead(x->v,&px));
  PetscCall(VecGetArrayRead(z,&py));
  if (!q) PetscCall(VecGetArray(X->buffer,&qq));
  PetscCall(BVDotVec_BLAS_Private(X,X->n,X->k-X->l,px+(X->nc+X->l)*X->ld,X->ld,py,qq,x->mpi));
  if (!q) PetscCall(VecRestoreArray(X->buffer,&qq));
  PetscCall(VecRestoreArrayRead(z,&py));
  PetscCall(VecRestoreArrayRead(x->v,&px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDotVec_Local_Svec(BV X,Vec y,PetscScalar *m)
{
  BV_SVEC        *x = (BV_SVEC*)X->data;
  PetscScalar    *px,*py;
  Vec            z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(VecGetArray(x->v,&px));
  PetscCall(VecGetArray(z,&py));
  PetscCall(BVDotVec_BLAS_Private(X,X->n,X->k-X->l,px+(X->nc+X->l)*X->ld,X->ld,py,m,PETSC_FALSE));
  PetscCall(VecRestoreArray(z,&py));
  PetscCall(VecRestoreArray(x->v,&px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVScale_Svec(BV bv,PetscInt j,PetscScalar alpha)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *array;

  PetscFunctionBegin;
  if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(VecGetArray(ctx->v,&array));
  if (PetscUnlikely(j<0)) PetscCall(BVScale_BLAS_Private(bv,(bv->k-bv->l)*bv->ld,array+(bv->nc+bv->l)*bv->ld,alpha));
  else PetscCall(BVScale_BLAS_Private(bv,bv->n,array+(bv->nc+j)*bv->ld,alpha));
  PetscCall(VecRestoreArray(ctx->v,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNorm_Svec(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_SVEC           *ctx = (BV_SVEC*)bv->data;
  const PetscScalar *array;

  PetscFunctionBegin;
  PetscCall(VecGetArrayRead(ctx->v,&array));
  if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,ctx->mpi));
  else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->ld,bv->ld,type,val,ctx->mpi));
  PetscCall(VecRestoreArrayRead(ctx->v,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNorm_Local_Svec(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_SVEC           *ctx = (BV_SVEC*)bv->data;
  const PetscScalar *array;

  PetscFunctionBegin;
  PetscCall(VecGetArrayRead(ctx->v,&array));
  if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,PETSC_FALSE));
  else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->ld,bv->ld,type,val,PETSC_FALSE));
  PetscCall(VecRestoreArrayRead(ctx->v,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNormalize_Svec(BV bv,PetscScalar *eigi)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *array,*wi=NULL;

  PetscFunctionBegin;
  PetscCall(VecGetArray(ctx->v,&array));
  if (eigi) wi = eigi+bv->l;
  PetscCall(BVNormalize_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,wi,ctx->mpi));
  PetscCall(VecRestoreArray(ctx->v,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMatMult_Svec(BV V,Mat A,BV W)
{
  PetscInt       j;
  Mat            Vmat,Wmat;
  Vec            vv,ww;

  PetscFunctionBegin;
  if (V->vmm) {
    PetscCall(BVGetMat(V,&Vmat));
    PetscCall(BVGetMat(W,&Wmat));
    PetscCall(MatProductCreateWithMat(A,Vmat,NULL,Wmat));
    PetscCall(MatProductSetType(Wmat,MATPRODUCT_AB));
    PetscCall(MatProductSetFromOptions(Wmat));
    PetscCall(MatProductSymbolic(Wmat));
    PetscCall(MatProductNumeric(Wmat));
    PetscCall(MatProductClear(Wmat));
    PetscCall(BVRestoreMat(V,&Vmat));
    PetscCall(BVRestoreMat(W,&Wmat));
  } else {
    for (j=0;j<V->k-V->l;j++) {
      PetscCall(BVGetColumn(V,V->l+j,&vv));
      PetscCall(BVGetColumn(W,W->l+j,&ww));
      PetscCall(MatMult(A,vv,ww));
      PetscCall(BVRestoreColumn(V,V->l+j,&vv));
      PetscCall(BVRestoreColumn(W,W->l+j,&ww));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVCopy_Svec(BV V,BV W)
{
  BV_SVEC           *v = (BV_SVEC*)V->data,*w = (BV_SVEC*)W->data;
  const PetscScalar *pv;
  PetscScalar       *pw;
  PetscInt          j;

  PetscFunctionBegin;
  PetscCall(VecGetArrayRead(v->v,&pv));
  PetscCall(VecGetArray(w->v,&pw));
  for (j=0;j<V->k-V->l;j++) PetscCall(PetscArraycpy(pw+(W->nc+W->l+j)*W->ld,pv+(V->nc+V->l+j)*V->ld,V->n));
  PetscCall(VecRestoreArrayRead(v->v,&pv));
  PetscCall(VecRestoreArray(w->v,&pw));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVCopyColumn_Svec(BV V,PetscInt j,PetscInt i)
{
  BV_SVEC        *v = (BV_SVEC*)V->data;
  PetscScalar    *pv;

  PetscFunctionBegin;
  PetscCall(VecGetArray(v->v,&pv));
  PetscCall(PetscArraycpy(pv+(V->nc+i)*V->ld,pv+(V->nc+j)*V->ld,V->n));
  PetscCall(VecRestoreArray(v->v,&pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVResize_Svec(BV bv,PetscInt m,PetscBool copy)
{
  BV_SVEC           *ctx = (BV_SVEC*)bv->data;
  PetscScalar       *pnew;
  const PetscScalar *pv;
  PetscInt          bs;
  Vec               vnew;
  char              str[50];

  PetscFunctionBegin;
  PetscCall(PetscLayoutGetBlockSize(bv->map,&bs));
  PetscCall(VecCreate(PetscObjectComm((PetscObject)bv),&vnew));
  PetscCall(VecSetType(vnew,bv->vtype));
  PetscCall(VecSetSizes(vnew,m*bv->ld,PETSC_DECIDE));
  PetscCall(VecSetBlockSize(vnew,bs));
  if (((PetscObject)bv)->name) {
    PetscCall(PetscSNPrintf(str,sizeof(str),"%s_0",((PetscObject)bv)->name));
    PetscCall(PetscObjectSetName((PetscObject)vnew,str));
  }
  if (copy) {
    PetscCall(VecGetArrayRead(ctx->v,&pv));
    PetscCall(VecGetArray(vnew,&pnew));
    PetscCall(PetscArraycpy(pnew,pv,PetscMin(m,bv->m)*bv->ld));
    PetscCall(VecRestoreArrayRead(ctx->v,&pv));
    PetscCall(VecRestoreArray(vnew,&pnew));
  }
  PetscCall(VecDestroy(&ctx->v));
  ctx->v = vnew;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetColumn_Svec(BV bv,PetscInt j,Vec *v)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *pv;
  PetscInt       l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  PetscCall(VecGetArray(ctx->v,&pv));
  PetscCall(VecPlaceArray(bv->cv[l],pv+(bv->nc+j)*bv->ld));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVRestoreColumn_Svec(BV bv,PetscInt j,Vec *v)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscInt       l;

  PetscFunctionBegin;
  l = (j==bv->ci[0])? 0: 1;
  PetscCall(VecResetArray(bv->cv[l]));
  PetscCall(VecRestoreArray(ctx->v,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetArray_Svec(BV bv,PetscScalar **a)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  PetscCall(VecGetArray(ctx->v,a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVRestoreArray_Svec(BV bv,PetscScalar **a)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  PetscCall(VecRestoreArray(ctx->v,a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetArrayRead_Svec(BV bv,const PetscScalar **a)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  PetscCall(VecGetArrayRead(ctx->v,a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVRestoreArrayRead_Svec(BV bv,const PetscScalar **a)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  PetscCall(VecRestoreArrayRead(ctx->v,a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVView_Svec(BV bv,PetscViewer viewer)
{
  PetscInt          j;
  Vec               v;
  PetscViewerFormat format;
  PetscBool         isascii,ismatlab=PETSC_FALSE;
  const char        *bvname,*name;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(PETSC_SUCCESS);
    if (format == PETSC_VIEWER_ASCII_MATLAB) ismatlab = PETSC_TRUE;
  }
  if (ismatlab) {
    PetscCall(PetscObjectGetName((PetscObject)bv,&bvname));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%s=[];\n",bvname));
  }
  for (j=0;j<bv->m;j++) {
    PetscCall(BVGetColumn(bv,j,&v));
    PetscCall(VecView(v,viewer));
    if (ismatlab) {
      PetscCall(PetscObjectGetName((PetscObject)v,&name));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%s=[%s,%s];clear %s\n",bvname,bvname,name,name));
    }
    PetscCall(BVRestoreColumn(bv,j,&v));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDestroy_Svec(BV bv)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  PetscCall(VecDestroy(&ctx->v));
  PetscCall(VecDestroy(&bv->cv[0]));
  PetscCall(VecDestroy(&bv->cv[1]));
  PetscCall(PetscFree(bv->data));
  bv->cuda = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   BVSVEC - BVSVEC = "svec" - A basis vectors object represented by a single
   `Vec` object. The local array of the `Vec` stores the local part of the
   basis vectors column-wise.

   Level: beginner

.seealso: [](sec:bv), `BV`, `BVType`, `BVSetType()`
M*/

SLEPC_EXTERN PetscErrorCode BVCreate_Svec(BV bv)
{
  BV_SVEC           *ctx;
  PetscInt          nloc,N,bs,tglobal=0,tlocal,lsplit,j,lda;
  PetscBool         seq,isdense;
  PetscScalar       *vv;
  const PetscScalar *aa,*array,*ptr;
  char              str[50];
  BV                parent;
  Vec               vpar;
  MatType           mtype;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  bv->data = (void*)ctx;

  PetscCall(PetscStrcmpAny(bv->vtype,&bv->cuda,VECSEQCUDA,VECMPICUDA,""));
  PetscCall(PetscStrcmpAny(bv->vtype,&bv->hip,VECSEQHIP,VECMPIHIP,""));
  PetscCall(PetscStrcmpAny(bv->vtype,&ctx->mpi,VECMPI,VECMPICUDA,VECMPIHIP,""));

  PetscCall(PetscStrcmp(bv->vtype,VECSEQ,&seq));
  PetscCheck(seq || ctx->mpi || bv->cuda || bv->hip,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"BVSVEC does not support the requested vector type: %s",bv->vtype);

  PetscCall(PetscLayoutGetLocalSize(bv->map,&nloc));
  PetscCall(PetscLayoutGetSize(bv->map,&N));
  PetscCall(PetscLayoutGetBlockSize(bv->map,&bs));
  PetscCall(BV_SetDefaultLD(bv,nloc));
  tlocal = bv->m*bv->ld;
  PetscCall(PetscIntMultError(bv->m,N,&tglobal));  /* just to check integer overflow */

  if (PetscUnlikely(bv->issplit)) {
    /* split BV: create Vec sharing the memory of the parent BV */
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    vpar = ((BV_SVEC*)parent->data)->v;
    if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
      PetscCall(VecCUDAGetArrayRead(vpar,&array));
      if (bv->issplit>0) ptr = (bv->issplit==1)? array: array+lsplit*bv->ld;
      else ptr = (bv->issplit==1)? array: array-lsplit;
      PetscCall(VecCUDARestoreArrayRead(vpar,&array));
      if (ctx->mpi) PetscCall(VecCreateMPICUDAWithArray(PetscObjectComm((PetscObject)bv),bs,tlocal,PETSC_DECIDE,NULL,&ctx->v));
      else PetscCall(VecCreateSeqCUDAWithArray(PetscObjectComm((PetscObject)bv),bs,tlocal,NULL,&ctx->v));
      PetscCall(VecCUDAPlaceArray(ctx->v,ptr));
#endif
    } else if (bv->hip) {
#if defined(PETSC_HAVE_HIP)
      PetscCall(VecHIPGetArrayRead(vpar,&array));
      if (bv->issplit>0) ptr = (bv->issplit==1)? array: array+lsplit*bv->ld;
      else ptr = (bv->issplit==1)? array: array-lsplit;
      PetscCall(VecHIPRestoreArrayRead(vpar,&array));
      if (ctx->mpi) PetscCall(VecCreateMPIHIPWithArray(PetscObjectComm((PetscObject)bv),bs,tlocal,PETSC_DECIDE,NULL,&ctx->v));
      else PetscCall(VecCreateSeqHIPWithArray(PetscObjectComm((PetscObject)bv),bs,tlocal,NULL,&ctx->v));
      PetscCall(VecHIPPlaceArray(ctx->v,ptr));
#endif
    } else {
      PetscCall(VecGetArrayRead(vpar,&array));
      if (bv->issplit>0) ptr = (bv->issplit==1)? array: array+lsplit*bv->ld;
      else ptr = (bv->issplit==1)? array: array-lsplit;
      PetscCall(VecRestoreArrayRead(vpar,&array));
      if (ctx->mpi) PetscCall(VecCreateMPIWithArray(PetscObjectComm((PetscObject)bv),bs,tlocal,PETSC_DECIDE,NULL,&ctx->v));
      else PetscCall(VecCreateSeqWithArray(PetscObjectComm((PetscObject)bv),bs,tlocal,NULL,&ctx->v));
      PetscCall(VecPlaceArray(ctx->v,ptr));
    }
  } else {
    /* regular BV: create Vec to store the BV entries */
    PetscCall(VecCreate(PetscObjectComm((PetscObject)bv),&ctx->v));
    PetscCall(VecSetType(ctx->v,bv->vtype));
    PetscCall(VecSetSizes(ctx->v,tlocal,PETSC_DECIDE));
    PetscCall(VecSetBlockSize(ctx->v,bs));
  }
  if (((PetscObject)bv)->name) {
    PetscCall(PetscSNPrintf(str,sizeof(str),"%s_0",((PetscObject)bv)->name));
    PetscCall(PetscObjectSetName((PetscObject)ctx->v,str));
  }

  if (PetscUnlikely(bv->Acreate)) {
    PetscCall(MatGetType(bv->Acreate,&mtype));
    PetscCall(PetscStrcmpAny(mtype,&isdense,MATSEQDENSE,MATMPIDENSE,""));
    PetscCheck(isdense,PetscObjectComm((PetscObject)bv->Acreate),PETSC_ERR_SUP,"BVSVEC requires a dense matrix in BVCreateFromMat()");
    PetscCall(MatDenseGetArrayRead(bv->Acreate,&aa));
    PetscCall(MatDenseGetLDA(bv->Acreate,&lda));
    PetscCall(VecGetArray(ctx->v,&vv));
    for (j=0;j<bv->m;j++) PetscCall(PetscArraycpy(vv+j*bv->ld,aa+j*lda,bv->n));
    PetscCall(VecRestoreArray(ctx->v,&vv));
    PetscCall(MatDenseRestoreArrayRead(bv->Acreate,&aa));
    PetscCall(MatDestroy(&bv->Acreate));
  }

  PetscCall(BVCreateVecEmpty(bv,&bv->cv[0]));
  PetscCall(BVCreateVecEmpty(bv,&bv->cv[1]));

  if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
    bv->ops->mult             = BVMult_Svec_CUDA;
    bv->ops->multvec          = BVMultVec_Svec_CUDA;
    bv->ops->multinplace      = BVMultInPlace_Svec_CUDA;
    bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Svec_CUDA;
    bv->ops->dot              = BVDot_Svec_CUDA;
    bv->ops->dotvec           = BVDotVec_Svec_CUDA;
    bv->ops->dotvec_local     = BVDotVec_Local_Svec_CUDA;
    bv->ops->scale            = BVScale_Svec_CUDA;
    bv->ops->norm             = BVNorm_Svec_CUDA;
    bv->ops->norm_local       = BVNorm_Local_Svec_CUDA;
    bv->ops->normalize        = BVNormalize_Svec_CUDA;
    bv->ops->matmult          = BVMatMult_Svec_CUDA;
    bv->ops->copy             = BVCopy_Svec_CUDA;
    bv->ops->copycolumn       = BVCopyColumn_Svec_CUDA;
    bv->ops->resize           = BVResize_Svec_CUDA;
    bv->ops->getcolumn        = BVGetColumn_Svec_CUDA;
    bv->ops->restorecolumn    = BVRestoreColumn_Svec_CUDA;
    bv->ops->restoresplit     = BVRestoreSplit_Svec_CUDA;
    bv->ops->restoresplitrows = BVRestoreSplitRows_Svec_CUDA;
    bv->ops->getmat           = BVGetMat_Svec_CUDA;
    bv->ops->restoremat       = BVRestoreMat_Svec_CUDA;
#endif
  } else if (bv->hip) {
#if defined(PETSC_HAVE_HIP)
    bv->ops->mult             = BVMult_Svec_HIP;
    bv->ops->multvec          = BVMultVec_Svec_HIP;
    bv->ops->multinplace      = BVMultInPlace_Svec_HIP;
    bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Svec_HIP;
    bv->ops->dot              = BVDot_Svec_HIP;
    bv->ops->dotvec           = BVDotVec_Svec_HIP;
    bv->ops->dotvec_local     = BVDotVec_Local_Svec_HIP;
    bv->ops->scale            = BVScale_Svec_HIP;
    bv->ops->norm             = BVNorm_Svec_HIP;
    bv->ops->norm_local       = BVNorm_Local_Svec_HIP;
    bv->ops->normalize        = BVNormalize_Svec_HIP;
    bv->ops->matmult          = BVMatMult_Svec_HIP;
    bv->ops->copy             = BVCopy_Svec_HIP;
    bv->ops->copycolumn       = BVCopyColumn_Svec_HIP;
    bv->ops->resize           = BVResize_Svec_HIP;
    bv->ops->getcolumn        = BVGetColumn_Svec_HIP;
    bv->ops->restorecolumn    = BVRestoreColumn_Svec_HIP;
    bv->ops->restoresplit     = BVRestoreSplit_Svec_HIP;
    bv->ops->restoresplitrows = BVRestoreSplitRows_Svec_HIP;
    bv->ops->getmat           = BVGetMat_Svec_HIP;
    bv->ops->restoremat       = BVRestoreMat_Svec_HIP;
#endif
  } else {
    bv->ops->mult             = BVMult_Svec;
    bv->ops->multvec          = BVMultVec_Svec;
    bv->ops->multinplace      = BVMultInPlace_Svec;
    bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Svec;
    bv->ops->dot              = BVDot_Svec;
    bv->ops->dotvec           = BVDotVec_Svec;
    bv->ops->dotvec_local     = BVDotVec_Local_Svec;
    bv->ops->scale            = BVScale_Svec;
    bv->ops->norm             = BVNorm_Svec;
    bv->ops->norm_local       = BVNorm_Local_Svec;
    bv->ops->normalize        = BVNormalize_Svec;
    bv->ops->matmult          = BVMatMult_Svec;
    bv->ops->copy             = BVCopy_Svec;
    bv->ops->copycolumn       = BVCopyColumn_Svec;
    bv->ops->resize           = BVResize_Svec;
    bv->ops->getcolumn        = BVGetColumn_Svec;
    bv->ops->restorecolumn    = BVRestoreColumn_Svec;
    bv->ops->getmat           = BVGetMat_Default;
    bv->ops->restoremat       = BVRestoreMat_Default;
  }
  bv->ops->getarray         = BVGetArray_Svec;
  bv->ops->restorearray     = BVRestoreArray_Svec;
  bv->ops->getarrayread     = BVGetArrayRead_Svec;
  bv->ops->restorearrayread = BVRestoreArrayRead_Svec;
  bv->ops->destroy          = BVDestroy_Svec;
  if (!ctx->mpi) bv->ops->view = BVView_Svec;
  PetscFunctionReturn(PETSC_SUCCESS);
}
