/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented as an array of Vecs sharing a contiguous array for elements
*/

#include <slepc/private/bvimpl.h>

typedef struct {
  Vec         *V;
  PetscScalar *array;
  PetscBool   mpi;
} BV_CONTIGUOUS;

static PetscErrorCode BVMult_Contiguous(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  BV_CONTIGUOUS     *y = (BV_CONTIGUOUS*)Y->data,*x = (BV_CONTIGUOUS*)X->data;
  const PetscScalar *q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (Q) {
    PetscCall(MatDenseGetLDA(Q,&ldq));
    PetscCall(MatDenseGetArrayRead(Q,&q));
    PetscCall(BVMult_BLAS_Private(Y,Y->n,Y->k-Y->l,X->k-X->l,alpha,x->array+(X->nc+X->l)*X->ld,X->ld,q+Y->l*ldq+X->l,ldq,beta,y->array+(Y->nc+Y->l)*Y->ld,Y->ld));
    PetscCall(MatDenseRestoreArrayRead(Q,&q));
  } else PetscCall(BVAXPY_BLAS_Private(Y,Y->n,Y->k-Y->l,alpha,x->array+(X->nc+X->l)*X->ld,X->ld,beta,y->array+(Y->nc+Y->l)*Y->ld,Y->ld));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultVec_Contiguous(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  BV_CONTIGUOUS  *x = (BV_CONTIGUOUS*)X->data;
  PetscScalar    *py,*qq=q;

  PetscFunctionBegin;
  PetscCall(VecGetArray(y,&py));
  if (!q) PetscCall(VecGetArray(X->buffer,&qq));
  PetscCall(BVMultVec_BLAS_Private(X,X->n,X->k-X->l,alpha,x->array+(X->nc+X->l)*X->ld,X->ld,qq,beta,py));
  if (!q) PetscCall(VecRestoreArray(X->buffer,&qq));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultInPlace_Contiguous(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_CONTIGUOUS     *ctx = (BV_CONTIGUOUS*)V->data;
  const PetscScalar *q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,s-V->l,e-V->l,ctx->array+(V->nc+V->l)*V->ld,V->ld,q+V->l*ldq+V->l,ldq,PETSC_FALSE));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultInPlaceHermitianTranspose_Contiguous(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_CONTIGUOUS     *ctx = (BV_CONTIGUOUS*)V->data;
  const PetscScalar *q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,s-V->l,e-V->l,ctx->array+(V->nc+V->l)*V->ld,V->ld,q+V->l*ldq+V->l,ldq,PETSC_TRUE));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDot_Contiguous(BV X,BV Y,Mat M)
{
  BV_CONTIGUOUS  *x = (BV_CONTIGUOUS*)X->data,*y = (BV_CONTIGUOUS*)Y->data;
  PetscScalar    *m;
  PetscInt       ldm;

  PetscFunctionBegin;
  PetscCall(MatDenseGetLDA(M,&ldm));
  PetscCall(MatDenseGetArray(M,&m));
  PetscCall(BVDot_BLAS_Private(X,Y->k-Y->l,X->k-X->l,X->n,y->array+(Y->nc+Y->l)*Y->ld,Y->ld,x->array+(X->nc+X->l)*X->ld,X->ld,m+X->l*ldm+Y->l,ldm,x->mpi));
  PetscCall(MatDenseRestoreArray(M,&m));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDotVec_Contiguous(BV X,Vec y,PetscScalar *q)
{
  BV_CONTIGUOUS     *x = (BV_CONTIGUOUS*)X->data;
  const PetscScalar *py;
  PetscScalar       *qq=q;
  Vec               z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(VecGetArrayRead(z,&py));
  if (!q) PetscCall(VecGetArray(X->buffer,&qq));
  PetscCall(BVDotVec_BLAS_Private(X,X->n,X->k-X->l,x->array+(X->nc+X->l)*X->ld,X->ld,py,qq,x->mpi));
  if (!q) PetscCall(VecRestoreArray(X->buffer,&qq));
  PetscCall(VecRestoreArrayRead(z,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDotVec_Local_Contiguous(BV X,Vec y,PetscScalar *m)
{
  BV_CONTIGUOUS  *x = (BV_CONTIGUOUS*)X->data;
  PetscScalar    *py;
  Vec            z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(VecGetArray(z,&py));
  PetscCall(BVDotVec_BLAS_Private(X,X->n,X->k-X->l,x->array+(X->nc+X->l)*X->ld,X->ld,py,m,PETSC_FALSE));
  PetscCall(VecRestoreArray(z,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVScale_Contiguous(BV bv,PetscInt j,PetscScalar alpha)
{
  BV_CONTIGUOUS  *ctx = (BV_CONTIGUOUS*)bv->data;

  PetscFunctionBegin;
  if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
  if (PetscUnlikely(j<0)) PetscCall(BVScale_BLAS_Private(bv,(bv->k-bv->l)*bv->ld,ctx->array+(bv->nc+bv->l)*bv->ld,alpha));
  else PetscCall(BVScale_BLAS_Private(bv,bv->n,ctx->array+(bv->nc+j)*bv->ld,alpha));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNorm_Contiguous(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_CONTIGUOUS  *ctx = (BV_CONTIGUOUS*)bv->data;

  PetscFunctionBegin;
  if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,ctx->array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,ctx->mpi));
  else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,ctx->array+(bv->nc+j)*bv->ld,bv->ld,type,val,ctx->mpi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNorm_Local_Contiguous(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_CONTIGUOUS  *ctx = (BV_CONTIGUOUS*)bv->data;

  PetscFunctionBegin;
  if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,ctx->array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,PETSC_FALSE));
  else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,ctx->array+(bv->nc+j)*bv->ld,bv->ld,type,val,PETSC_FALSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNormalize_Contiguous(BV bv,PetscScalar *eigi)
{
  BV_CONTIGUOUS  *ctx = (BV_CONTIGUOUS*)bv->data;
  PetscScalar    *wi=NULL;

  PetscFunctionBegin;
  if (eigi) wi = eigi+bv->l;
  PetscCall(BVNormalize_LAPACK_Private(bv,bv->n,bv->k-bv->l,ctx->array+(bv->nc+bv->l)*bv->ld,bv->ld,wi,ctx->mpi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMatMult_Contiguous(BV V,Mat A,BV W)
{
  BV_CONTIGUOUS  *v = (BV_CONTIGUOUS*)V->data,*w = (BV_CONTIGUOUS*)W->data;
  PetscInt       j;
  Mat            Vmat,Wmat;

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
    for (j=0;j<V->k-V->l;j++) PetscCall(MatMult(A,v->V[V->nc+V->l+j],w->V[W->nc+W->l+j]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVCopy_Contiguous(BV V,BV W)
{
  BV_CONTIGUOUS  *v = (BV_CONTIGUOUS*)V->data,*w = (BV_CONTIGUOUS*)W->data;
  PetscInt       j;

  PetscFunctionBegin;
  for (j=0;j<V->k-V->l;j++) PetscCall(PetscArraycpy(w->array+(W->nc+W->l+j)*W->ld,v->array+(V->nc+V->l+j)*V->ld,V->n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVCopyColumn_Contiguous(BV V,PetscInt j,PetscInt i)
{
  BV_CONTIGUOUS  *v = (BV_CONTIGUOUS*)V->data;

  PetscFunctionBegin;
  PetscCall(PetscArraycpy(v->array+(V->nc+i)*V->ld,v->array+(V->nc+j)*V->ld,V->n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVResize_Contiguous(BV bv,PetscInt m,PetscBool copy)
{
  BV_CONTIGUOUS  *ctx = (BV_CONTIGUOUS*)bv->data;
  PetscInt       j,bs;
  PetscScalar    *newarray;
  Vec            *newV;
  char           str[50];

  PetscFunctionBegin;
  PetscCall(PetscLayoutGetBlockSize(bv->map,&bs));
  PetscCall(PetscCalloc1(m*bv->ld,&newarray));
  PetscCall(PetscMalloc1(m,&newV));
  for (j=0;j<m;j++) {
    if (ctx->mpi) PetscCall(VecCreateMPIWithArray(PetscObjectComm((PetscObject)bv),bs,bv->n,PETSC_DECIDE,newarray+j*bv->ld,newV+j));
    else PetscCall(VecCreateSeqWithArray(PetscObjectComm((PetscObject)bv),bs,bv->n,newarray+j*bv->ld,newV+j));
  }
  if (((PetscObject)bv)->name) {
    for (j=0;j<m;j++) {
      PetscCall(PetscSNPrintf(str,sizeof(str),"%s_%" PetscInt_FMT,((PetscObject)bv)->name,j));
      PetscCall(PetscObjectSetName((PetscObject)newV[j],str));
    }
  }
  if (copy) PetscCall(PetscArraycpy(newarray,ctx->array,PetscMin(m,bv->m)*bv->ld));
  PetscCall(VecDestroyVecs(bv->m,&ctx->V));
  ctx->V = newV;
  PetscCall(PetscFree(ctx->array));
  ctx->array = newarray;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetColumn_Contiguous(BV bv,PetscInt j,Vec *v)
{
  BV_CONTIGUOUS *ctx = (BV_CONTIGUOUS*)bv->data;
  PetscInt      l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  bv->cv[l] = ctx->V[bv->nc+j];
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVRestoreColumn_Contiguous(BV bv,PetscInt j,Vec *v)
{
  PetscInt l;

  PetscFunctionBegin;
  l = (j==bv->ci[0])? 0: 1;
  bv->cv[l] = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetArray_Contiguous(BV bv,PetscScalar **a)
{
  BV_CONTIGUOUS *ctx = (BV_CONTIGUOUS*)bv->data;

  PetscFunctionBegin;
  *a = ctx->array;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetArrayRead_Contiguous(BV bv,const PetscScalar **a)
{
  BV_CONTIGUOUS *ctx = (BV_CONTIGUOUS*)bv->data;

  PetscFunctionBegin;
  *a = ctx->array;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDestroy_Contiguous(BV bv)
{
  BV_CONTIGUOUS  *ctx = (BV_CONTIGUOUS*)bv->data;

  PetscFunctionBegin;
  if (!bv->issplit) {
    PetscCall(VecDestroyVecs(bv->nc+bv->m,&ctx->V));
    PetscCall(PetscFree(ctx->array));
  }
  PetscCall(PetscFree(bv->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVView_Contiguous(BV bv,PetscViewer viewer)
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

/*MC
   BVCONTIGUOUS - BVCONTIGUOUS = "contiguous" - A basis vectors object represented
   by an array of `Vec`s sharing a contiguous array of `PetscScalar` storing the
   elements column-wise.

   Level: beginner

.seealso: [](sec:bv), `BV`, `BVType`, `BVSetType()`
M*/

SLEPC_EXTERN PetscErrorCode BVCreate_Contiguous(BV bv)
{
  BV_CONTIGUOUS  *ctx;
  PetscInt       j,nloc,bs,lsplit,lda;
  PetscBool      seq,isdense;
  PetscScalar    *aa;
  char           str[50];
  PetscScalar    *array;
  BV             parent;
  Vec            *Vpar;
  MatType        mtype;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  bv->data = (void*)ctx;

  PetscCall(PetscStrcmp(bv->vtype,VECMPI,&ctx->mpi));
  if (!ctx->mpi) {
    PetscCall(PetscStrcmp(bv->vtype,VECSEQ,&seq));
    PetscCheck(seq,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Cannot create a contiguous BV from a non-standard vector type: %s",bv->vtype);
  }

  PetscCall(PetscLayoutGetLocalSize(bv->map,&nloc));
  PetscCall(PetscLayoutGetBlockSize(bv->map,&bs));
  PetscCall(BV_SetDefaultLD(bv,nloc));

  if (PetscUnlikely(bv->issplit)) {
    PetscCheck(bv->issplit>0,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"BVCONTIGUOUS does not support BVGetSplitRows()");
    /* split BV: share memory and Vecs of the parent BV */
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    Vpar   = ((BV_CONTIGUOUS*)parent->data)->V;
    ctx->V = (bv->issplit==1)? Vpar: Vpar+lsplit;
    array  = ((BV_CONTIGUOUS*)parent->data)->array;
    ctx->array = (bv->issplit==1)? array: array+lsplit*bv->ld;
  } else {
    /* regular BV: allocate memory and Vecs for the BV entries */
    PetscCall(PetscCalloc1(bv->m*bv->ld,&ctx->array));
    PetscCall(PetscMalloc1(bv->m,&ctx->V));
    for (j=0;j<bv->m;j++) {
      if (ctx->mpi) PetscCall(VecCreateMPIWithArray(PetscObjectComm((PetscObject)bv),bs,nloc,PETSC_DECIDE,ctx->array+j*bv->ld,ctx->V+j));
      else PetscCall(VecCreateSeqWithArray(PetscObjectComm((PetscObject)bv),bs,nloc,ctx->array+j*bv->ld,ctx->V+j));
    }
  }
  if (((PetscObject)bv)->name) {
    for (j=0;j<bv->m;j++) {
      PetscCall(PetscSNPrintf(str,sizeof(str),"%s_%" PetscInt_FMT,((PetscObject)bv)->name,j));
      PetscCall(PetscObjectSetName((PetscObject)ctx->V[j],str));
    }
  }

  if (PetscUnlikely(bv->Acreate)) {
    PetscCall(MatGetType(bv->Acreate,&mtype));
    PetscCall(PetscStrcmpAny(mtype,&isdense,MATSEQDENSE,MATMPIDENSE,""));
    PetscCheck(isdense,PetscObjectComm((PetscObject)bv->Acreate),PETSC_ERR_SUP,"BVCONTIGUOUS requires a dense matrix in BVCreateFromMat()");
    PetscCall(MatDenseGetArray(bv->Acreate,&aa));
    PetscCall(MatDenseGetLDA(bv->Acreate,&lda));
    for (j=0;j<bv->m;j++) PetscCall(PetscArraycpy(ctx->array+j*bv->ld,aa+j*lda,bv->n));
    PetscCall(MatDenseRestoreArray(bv->Acreate,&aa));
    PetscCall(MatDestroy(&bv->Acreate));
  }

  bv->ops->mult             = BVMult_Contiguous;
  bv->ops->multvec          = BVMultVec_Contiguous;
  bv->ops->multinplace      = BVMultInPlace_Contiguous;
  bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Contiguous;
  bv->ops->dot              = BVDot_Contiguous;
  bv->ops->dotvec           = BVDotVec_Contiguous;
  bv->ops->dotvec_local     = BVDotVec_Local_Contiguous;
  bv->ops->scale            = BVScale_Contiguous;
  bv->ops->norm             = BVNorm_Contiguous;
  bv->ops->norm_local       = BVNorm_Local_Contiguous;
  bv->ops->normalize        = BVNormalize_Contiguous;
  bv->ops->matmult          = BVMatMult_Contiguous;
  bv->ops->copy             = BVCopy_Contiguous;
  bv->ops->copycolumn       = BVCopyColumn_Contiguous;
  bv->ops->resize           = BVResize_Contiguous;
  bv->ops->getcolumn        = BVGetColumn_Contiguous;
  bv->ops->restorecolumn    = BVRestoreColumn_Contiguous;
  bv->ops->getarray         = BVGetArray_Contiguous;
  bv->ops->getarrayread     = BVGetArrayRead_Contiguous;
  bv->ops->getmat           = BVGetMat_Default;
  bv->ops->restoremat       = BVRestoreMat_Default;
  bv->ops->destroy          = BVDestroy_Contiguous;
  bv->ops->view             = BVView_Contiguous;
  PetscFunctionReturn(PETSC_SUCCESS);
}
