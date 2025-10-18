/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented with a dense Mat
*/

#include <slepc/private/bvimpl.h>
#include "bvmat.h"

static PetscErrorCode BVMult_Mat(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  BV_MAT            *y = (BV_MAT*)Y->data,*x = (BV_MAT*)X->data;
  PetscScalar       *py;
  const PetscScalar *px,*q;
  PetscInt          ldq;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(x->A,&px));
  PetscCall(MatDenseGetArray(y->A,&py));
  if (Q) {
    PetscCall(MatDenseGetLDA(Q,&ldq));
    PetscCall(MatDenseGetArrayRead(Q,&q));
    PetscCall(BVMult_BLAS_Private(Y,Y->n,Y->k-Y->l,X->k-X->l,alpha,px+(X->nc+X->l)*X->ld,X->ld,q+Y->l*ldq+X->l,ldq,beta,py+(Y->nc+Y->l)*Y->ld,Y->ld));
    PetscCall(MatDenseRestoreArrayRead(Q,&q));
  } else PetscCall(BVAXPY_BLAS_Private(Y,Y->n,Y->k-Y->l,alpha,px+(X->nc+X->l)*X->ld,X->ld,beta,py+(Y->nc+Y->l)*Y->ld,Y->ld));
  PetscCall(MatDenseRestoreArrayRead(x->A,&px));
  PetscCall(MatDenseRestoreArray(y->A,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultVec_Mat(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  BV_MAT            *x = (BV_MAT*)X->data;
  PetscScalar       *py,*qq=q;
  const PetscScalar *px;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(x->A,&px));
  PetscCall(VecGetArray(y,&py));
  if (!q) PetscCall(VecGetArray(X->buffer,&qq));
  PetscCall(BVMultVec_BLAS_Private(X,X->n,X->k-X->l,alpha,px+(X->nc+X->l)*X->ld,X->ld,qq,beta,py));
  if (!q) PetscCall(VecRestoreArray(X->buffer,&qq));
  PetscCall(MatDenseRestoreArrayRead(x->A,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultInPlace_Mat(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_MAT            *ctx = (BV_MAT*)V->data;
  PetscScalar       *pv;
  const PetscScalar *q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(MatDenseGetArray(ctx->A,&pv));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,s-V->l,e-V->l,pv+(V->nc+V->l)*V->ld,V->ld,q+V->l*ldq+V->l,ldq,PETSC_FALSE));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscCall(MatDenseRestoreArray(ctx->A,&pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultInPlaceHermitianTranspose_Mat(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_MAT            *ctx = (BV_MAT*)V->data;
  PetscScalar       *pv;
  const PetscScalar *q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(MatDenseGetArray(ctx->A,&pv));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,s-V->l,e-V->l,pv+(V->nc+V->l)*V->ld,V->ld,q+V->l*ldq+V->l,ldq,PETSC_TRUE));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscCall(MatDenseRestoreArray(ctx->A,&pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDot_Mat(BV X,BV Y,Mat M)
{
  BV_MAT            *x = (BV_MAT*)X->data,*y = (BV_MAT*)Y->data;
  PetscScalar       *m;
  const PetscScalar *px,*py;
  PetscInt          ldm;

  PetscFunctionBegin;
  PetscCall(MatDenseGetLDA(M,&ldm));
  PetscCall(MatDenseGetArrayRead(x->A,&px));
  PetscCall(MatDenseGetArrayRead(y->A,&py));
  PetscCall(MatDenseGetArray(M,&m));
  PetscCall(BVDot_BLAS_Private(X,Y->k-Y->l,X->k-X->l,X->n,py+(Y->nc+Y->l)*Y->ld,Y->ld,px+(X->nc+X->l)*X->ld,X->ld,m+X->l*ldm+Y->l,ldm,x->mpi));
  PetscCall(MatDenseRestoreArray(M,&m));
  PetscCall(MatDenseRestoreArrayRead(x->A,&px));
  PetscCall(MatDenseRestoreArrayRead(y->A,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDotVec_Mat(BV X,Vec y,PetscScalar *q)
{
  BV_MAT            *x = (BV_MAT*)X->data;
  PetscScalar       *qq=q;
  const PetscScalar *px,*py;
  Vec               z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(MatDenseGetArrayRead(x->A,&px));
  PetscCall(VecGetArrayRead(z,&py));
  if (!q) PetscCall(VecGetArray(X->buffer,&qq));
  PetscCall(BVDotVec_BLAS_Private(X,X->n,X->k-X->l,px+(X->nc+X->l)*X->ld,X->ld,py,qq,x->mpi));
  if (!q) PetscCall(VecRestoreArray(X->buffer,&qq));
  PetscCall(VecRestoreArrayRead(z,&py));
  PetscCall(MatDenseRestoreArrayRead(x->A,&px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDotVec_Local_Mat(BV X,Vec y,PetscScalar *m)
{
  BV_MAT            *x = (BV_MAT*)X->data;
  const PetscScalar *px,*py;
  Vec               z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(MatDenseGetArrayRead(x->A,&px));
  PetscCall(VecGetArrayRead(z,&py));
  PetscCall(BVDotVec_BLAS_Private(X,X->n,X->k-X->l,px+(X->nc+X->l)*X->ld,X->ld,py,m,PETSC_FALSE));
  PetscCall(VecRestoreArrayRead(z,&py));
  PetscCall(MatDenseRestoreArrayRead(x->A,&px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVScale_Mat(BV bv,PetscInt j,PetscScalar alpha)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *array;

  PetscFunctionBegin;
  if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetArray(ctx->A,&array));
  if (PetscUnlikely(j<0)) PetscCall(BVScale_BLAS_Private(bv,(bv->k-bv->l)*bv->ld,array+(bv->nc+bv->l)*bv->ld,alpha));
  else PetscCall(BVScale_BLAS_Private(bv,bv->n,array+(bv->nc+j)*bv->ld,alpha));
  PetscCall(MatDenseRestoreArray(ctx->A,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNorm_Mat(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_MAT            *ctx = (BV_MAT*)bv->data;
  const PetscScalar *array;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(ctx->A,&array));
  if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,ctx->mpi));
  else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->ld,bv->ld,type,val,ctx->mpi));
  PetscCall(MatDenseRestoreArrayRead(ctx->A,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNorm_Local_Mat(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_MAT            *ctx = (BV_MAT*)bv->data;
  const PetscScalar *array;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(ctx->A,&array));
  if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,PETSC_FALSE));
  else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->ld,bv->ld,type,val,PETSC_FALSE));
  PetscCall(MatDenseRestoreArrayRead(ctx->A,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNormalize_Mat(BV bv,PetscScalar *eigi)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *array,*wi=NULL;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ctx->A,&array));
  if (eigi) wi = eigi+bv->l;
  PetscCall(BVNormalize_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,wi,ctx->mpi));
  PetscCall(MatDenseRestoreArray(ctx->A,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMatMult_Mat(BV V,Mat A,BV W)
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

static PetscErrorCode BVCopy_Mat(BV V,BV W)
{
  BV_MAT            *v = (BV_MAT*)V->data,*w = (BV_MAT*)W->data;
  const PetscScalar *pv;
  PetscScalar       *pw;
  PetscInt          j;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(v->A,&pv));
  PetscCall(MatDenseGetArray(w->A,&pw));
  for (j=0;j<V->k-V->l;j++) PetscCall(PetscArraycpy(pw+(W->nc+W->l+j)*W->ld,pv+(V->nc+V->l+j)*V->ld,V->n));
  PetscCall(MatDenseRestoreArrayRead(v->A,&pv));
  PetscCall(MatDenseRestoreArray(w->A,&pw));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVCopyColumn_Mat(BV V,PetscInt j,PetscInt i)
{
  BV_MAT         *v = (BV_MAT*)V->data;
  PetscScalar    *pv;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(v->A,&pv));
  PetscCall(PetscArraycpy(pv+(V->nc+i)*V->ld,pv+(V->nc+j)*V->ld,V->n));
  PetscCall(MatDenseRestoreArray(v->A,&pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVResize_Mat(BV bv,PetscInt m,PetscBool copy)
{
  BV_MAT            *ctx = (BV_MAT*)bv->data;
  Mat               A,Msrc,Mdst;
  char              str[50];

  PetscFunctionBegin;
  PetscCall(MatCreateDenseFromVecType(PetscObjectComm((PetscObject)bv),bv->vtype,bv->n,PETSC_DECIDE,bv->N,m,bv->ld,NULL,&A));
  if (((PetscObject)bv)->name) {
    PetscCall(PetscSNPrintf(str,sizeof(str),"%s_0",((PetscObject)bv)->name));
    PetscCall(PetscObjectSetName((PetscObject)A,str));
  }
  if (copy) {
    PetscCall(MatDenseGetSubMatrix(ctx->A,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PetscMin(m,bv->m),&Msrc));
    PetscCall(MatDenseGetSubMatrix(A,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PetscMin(m,bv->m),&Mdst));
    PetscCall(MatCopy(Msrc,Mdst,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(ctx->A,&Msrc));
    PetscCall(MatDenseRestoreSubMatrix(A,&Mdst));
  }
  PetscCall(MatDestroy(&ctx->A));
  ctx->A = A;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetColumn_Mat(BV bv,PetscInt j,Vec *v)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *pA;
  PetscInt       l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  PetscCall(MatDenseGetArray(ctx->A,&pA));
  PetscCall(VecPlaceArray(bv->cv[l],pA+(bv->nc+j)*bv->ld));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVRestoreColumn_Mat(BV bv,PetscInt j,Vec *v)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *pA;
  PetscInt       l;

  PetscFunctionBegin;
  l = (j==bv->ci[0])? 0: 1;
  PetscCall(VecResetArray(bv->cv[l]));
  PetscCall(MatDenseRestoreArray(ctx->A,&pA));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetArray_Mat(BV bv,PetscScalar **a)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ctx->A,a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVRestoreArray_Mat(BV bv,PetscScalar **a)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  if (a) PetscCall(MatDenseRestoreArray(ctx->A,a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetArrayRead_Mat(BV bv,const PetscScalar **a)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(ctx->A,a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVRestoreArrayRead_Mat(BV bv,const PetscScalar **a)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  if (a) PetscCall(MatDenseRestoreArrayRead(ctx->A,a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVView_Mat(BV bv,PetscViewer viewer)
{
  Mat               A;
  PetscViewerFormat format;
  PetscBool         isascii;
  const char        *bvname,*name;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCall(BVGetMat(bv,&A));
  PetscCall(MatView(A,viewer));
  if (format == PETSC_VIEWER_ASCII_MATLAB) {
    PetscCall(PetscObjectGetName((PetscObject)A,&name));
    PetscCall(PetscObjectGetName((PetscObject)bv,&bvname));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%s=%s;clear %s\n",bvname,name,name));
    if (bv->nc) PetscCall(PetscViewerASCIIPrintf(viewer,"%s=%s(:,%" PetscInt_FMT ":end);\n",bvname,bvname,bv->nc+1));
  }
  PetscCall(BVRestoreMat(bv,&A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDestroy_Mat(BV bv)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->A));
  PetscCall(VecDestroy(&bv->cv[0]));
  PetscCall(VecDestroy(&bv->cv[1]));
  PetscCall(PetscFree(bv->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   BVMAT - BVMAT = "mat" - A basis vectors object represented by a dense `Mat`.

   Level: beginner

.seealso: [](sec:bv), `BV`, `BVType`, `BVSetType()`
M*/

SLEPC_EXTERN PetscErrorCode BVCreate_Mat(BV bv)
{
  BV_MAT         *ctx;
  PetscInt       nloc,lsplit;
  PetscBool      seq;
  char           str[50];
  PetscScalar    *array,*ptr=NULL;
  BV             parent;
  Mat            Apar;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  bv->data = (void*)ctx;

  PetscCall(PetscStrcmpAny(bv->vtype,&bv->cuda,VECSEQCUDA,VECMPICUDA,""));
  PetscCall(PetscStrcmpAny(bv->vtype,&bv->hip,VECSEQHIP,VECMPIHIP,""));
  PetscCall(PetscStrcmpAny(bv->vtype,&ctx->mpi,VECMPI,VECMPICUDA,VECMPIHIP,""));

  PetscCall(PetscStrcmp(bv->vtype,VECSEQ,&seq));
  PetscCheck(seq || ctx->mpi || bv->cuda || bv->hip,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"BVMAT does not support the requested vector type: %s",bv->vtype);

  PetscCall(PetscLayoutGetLocalSize(bv->map,&nloc));
  PetscCall(BV_SetDefaultLD(bv,nloc));

  if (PetscUnlikely(bv->issplit)) {
    /* split BV: share the memory of the parent BV */
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    Apar = ((BV_MAT*)parent->data)->A;
    if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
      PetscCall(MatDenseCUDAGetArray(Apar,&array));
      if (bv->issplit>0) ptr = (bv->issplit==1)? array: array+lsplit*bv->ld;
      else ptr = (bv->issplit==1)? array: array-lsplit;
      PetscCall(MatDenseCUDARestoreArray(Apar,&array));
#endif
    } else if (bv->hip) {
#if defined(PETSC_HAVE_HIP)
      PetscCall(MatDenseHIPGetArray(Apar,&array));
      if (bv->issplit>0) ptr = (bv->issplit==1)? array: array+lsplit*bv->ld;
      else ptr = (bv->issplit==1)? array: array-lsplit;
      PetscCall(MatDenseHIPRestoreArray(Apar,&array));
#endif
    } else {
      PetscCall(MatDenseGetArray(Apar,&array));
      if (bv->issplit>0) ptr = (bv->issplit==1)? array: array+lsplit*bv->ld;
      else ptr = (bv->issplit==-1)? array: array-lsplit;
      PetscCall(MatDenseRestoreArray(Apar,&array));
    }
  }

  PetscCall(MatCreateDenseFromVecType(PetscObjectComm((PetscObject)bv),bv->vtype,nloc,PETSC_DECIDE,bv->N,bv->m,bv->ld,ptr,&ctx->A));
  if (((PetscObject)bv)->name) {
    PetscCall(PetscSNPrintf(str,sizeof(str),"%s_0",((PetscObject)bv)->name));
    PetscCall(PetscObjectSetName((PetscObject)ctx->A,str));
  }

  if (PetscUnlikely(bv->Acreate)) {
    PetscCall(MatConvert(bv->Acreate,bv->cuda?MATDENSECUDA:bv->hip?MATDENSEHIP:MATDENSE,MAT_REUSE_MATRIX,&ctx->A));
    PetscCall(MatDestroy(&bv->Acreate));
  }

  PetscCall(BVCreateVecEmpty(bv,&bv->cv[0]));
  PetscCall(BVCreateVecEmpty(bv,&bv->cv[1]));

  if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
    bv->ops->mult             = BVMult_Mat_CUDA;
    bv->ops->multvec          = BVMultVec_Mat_CUDA;
    bv->ops->multinplace      = BVMultInPlace_Mat_CUDA;
    bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Mat_CUDA;
    bv->ops->dot              = BVDot_Mat_CUDA;
    bv->ops->dotvec           = BVDotVec_Mat_CUDA;
    bv->ops->dotvec_local     = BVDotVec_Local_Mat_CUDA;
    bv->ops->scale            = BVScale_Mat_CUDA;
    bv->ops->norm             = BVNorm_Mat_CUDA;
    bv->ops->norm_local       = BVNorm_Local_Mat_CUDA;
    bv->ops->normalize        = BVNormalize_Mat_CUDA;
    bv->ops->matmult          = BVMatMult_Mat_CUDA;
    bv->ops->copy             = BVCopy_Mat_CUDA;
    bv->ops->copycolumn       = BVCopyColumn_Mat_CUDA;
    bv->ops->getcolumn        = BVGetColumn_Mat_CUDA;
    bv->ops->restorecolumn    = BVRestoreColumn_Mat_CUDA;
    bv->ops->restoresplit     = BVRestoreSplit_Mat_CUDA;
    bv->ops->restoresplitrows = BVRestoreSplitRows_Mat_CUDA;
    bv->ops->getmat           = BVGetMat_Mat_CUDA;
    bv->ops->restoremat       = BVRestoreMat_Mat_CUDA;
#endif
  } else if (bv->hip) {
#if defined(PETSC_HAVE_HIP)
    bv->ops->mult             = BVMult_Mat_HIP;
    bv->ops->multvec          = BVMultVec_Mat_HIP;
    bv->ops->multinplace      = BVMultInPlace_Mat_HIP;
    bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Mat_HIP;
    bv->ops->dot              = BVDot_Mat_HIP;
    bv->ops->dotvec           = BVDotVec_Mat_HIP;
    bv->ops->dotvec_local     = BVDotVec_Local_Mat_HIP;
    bv->ops->scale            = BVScale_Mat_HIP;
    bv->ops->norm             = BVNorm_Mat_HIP;
    bv->ops->norm_local       = BVNorm_Local_Mat_HIP;
    bv->ops->normalize        = BVNormalize_Mat_HIP;
    bv->ops->matmult          = BVMatMult_Mat_HIP;
    bv->ops->copy             = BVCopy_Mat_HIP;
    bv->ops->copycolumn       = BVCopyColumn_Mat_HIP;
    bv->ops->getcolumn        = BVGetColumn_Mat_HIP;
    bv->ops->restorecolumn    = BVRestoreColumn_Mat_HIP;
    bv->ops->restoresplit     = BVRestoreSplit_Mat_HIP;
    bv->ops->restoresplitrows = BVRestoreSplitRows_Mat_HIP;
    bv->ops->getmat           = BVGetMat_Mat_HIP;
    bv->ops->restoremat       = BVRestoreMat_Mat_HIP;
#endif
  } else {
    bv->ops->mult             = BVMult_Mat;
    bv->ops->multvec          = BVMultVec_Mat;
    bv->ops->multinplace      = BVMultInPlace_Mat;
    bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Mat;
    bv->ops->dot              = BVDot_Mat;
    bv->ops->dotvec           = BVDotVec_Mat;
    bv->ops->dotvec_local     = BVDotVec_Local_Mat;
    bv->ops->scale            = BVScale_Mat;
    bv->ops->norm             = BVNorm_Mat;
    bv->ops->norm_local       = BVNorm_Local_Mat;
    bv->ops->normalize        = BVNormalize_Mat;
    bv->ops->matmult          = BVMatMult_Mat;
    bv->ops->copy             = BVCopy_Mat;
    bv->ops->copycolumn       = BVCopyColumn_Mat;
    bv->ops->getcolumn        = BVGetColumn_Mat;
    bv->ops->restorecolumn    = BVRestoreColumn_Mat;
    bv->ops->getmat           = BVGetMat_Default;
    bv->ops->restoremat       = BVRestoreMat_Default;
  }
  bv->ops->resize           = BVResize_Mat;
  bv->ops->getarray         = BVGetArray_Mat;
  bv->ops->restorearray     = BVRestoreArray_Mat;
  bv->ops->getarrayread     = BVGetArrayRead_Mat;
  bv->ops->restorearrayread = BVRestoreArrayRead_Mat;
  bv->ops->destroy          = BVDestroy_Mat;
  bv->ops->view             = BVView_Mat;
  PetscFunctionReturn(PETSC_SUCCESS);
}
