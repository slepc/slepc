/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2021, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented as an array of independent Vecs
*/

#include <slepc/private/bvimpl.h>

typedef struct {
  Vec      *V;
  PetscInt vmip;   /* Version of BVMultInPlace:
       0: memory-efficient version, uses VecGetArray (default in CPU)
       1: version that allocates (e-s) work vectors in every call (default in GPU) */
} BV_VECS;

PetscErrorCode BVMult_Vecs(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  BV_VECS           *y = (BV_VECS*)Y->data,*x = (BV_VECS*)X->data;
  PetscScalar       *s=NULL;
  const PetscScalar *q;
  PetscInt          i,j,ldq;
  PetscBool         trivial=(alpha==1.0)?PETSC_TRUE:PETSC_FALSE;

  PetscFunctionBegin;
  if (Q) {
    PetscCall(MatGetSize(Q,&ldq,NULL));
    if (!trivial) {
      PetscCall(BVAllocateWork_Private(Y,X->k-X->l));
      s = Y->work;
    }
    PetscCall(MatDenseGetArrayRead(Q,&q));
    for (j=Y->l;j<Y->k;j++) {
      PetscCall(VecScale(y->V[Y->nc+j],beta));
      if (!trivial) {
        for (i=X->l;i<X->k;i++) s[i-X->l] = alpha*q[i+j*ldq];
      } else s = (PetscScalar*)(q+j*ldq+X->l);
      PetscCall(VecMAXPY(y->V[Y->nc+j],X->k-X->l,s,x->V+X->nc+X->l));
    }
    PetscCall(MatDenseRestoreArrayRead(Q,&q));
  } else {
    for (j=0;j<Y->k-Y->l;j++) {
      PetscCall(VecScale(y->V[Y->nc+Y->l+j],beta));
      PetscCall(VecAXPY(y->V[Y->nc+Y->l+j],alpha,x->V[X->nc+X->l+j]));
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultVec_Vecs(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  BV_VECS        *x = (BV_VECS*)X->data;
  PetscScalar    *s=NULL,*qq=q;
  PetscInt       i;
  PetscBool      trivial=(alpha==1.0)?PETSC_TRUE:PETSC_FALSE;

  PetscFunctionBegin;
  if (!trivial) {
    PetscCall(BVAllocateWork_Private(X,X->k-X->l));
    s = X->work;
  }
  if (!q) PetscCall(VecGetArray(X->buffer,&qq));
  PetscCall(VecScale(y,beta));
  if (!trivial) {
    for (i=0;i<X->k-X->l;i++) s[i] = alpha*qq[i];
  } else s = qq;
  PetscCall(VecMAXPY(y,X->k-X->l,s,x->V+X->nc+X->l));
  if (!q) PetscCall(VecRestoreArray(X->buffer,&qq));
  PetscFunctionReturn(0);
}

/*
   BVMultInPlace_Vecs_ME - V(:,s:e-1) = V*Q(:,s:e-1) for regular vectors.

   Memory-efficient version, uses VecGetArray (default in CPU)

   Writing V = [ V1 V2 V3 ] and Q(:,s:e-1) = [ Q1 Q2 Q3 ]', where V2
   corresponds to the columns s:e-1, the computation is done as
                  V2 := V2*Q2 + V1*Q1 + V3*Q3
*/
PetscErrorCode BVMultInPlace_Vecs_ME(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_VECS           *ctx = (BV_VECS*)V->data;
  const PetscScalar *q;
  PetscInt          i,ldq;

  PetscFunctionBegin;
  PetscCall(MatGetSize(Q,&ldq,NULL));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  /* V2 := V2*Q2 */
  PetscCall(BVMultInPlace_Vecs_Private(V,V->n,e-s,ldq,ctx->V+V->nc+s,q+s*ldq+s,PETSC_FALSE));
  /* V2 += V1*Q1 + V3*Q3 */
  for (i=s;i<e;i++) {
    if (PetscUnlikely(s>V->l)) PetscCall(VecMAXPY(ctx->V[V->nc+i],s-V->l,q+i*ldq+V->l,ctx->V+V->nc+V->l));
    if (V->k>e) PetscCall(VecMAXPY(ctx->V[V->nc+i],V->k-e,q+i*ldq+e,ctx->V+V->nc+e));
  }
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscFunctionReturn(0);
}

/*
   BVMultInPlace_Vecs_Alloc - V(:,s:e-1) = V*Q(:,s:e-1) for regular vectors.

   Version that allocates (e-s) work vectors in every call (default in GPU)
*/
PetscErrorCode BVMultInPlace_Vecs_Alloc(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_VECS           *ctx = (BV_VECS*)V->data;
  const PetscScalar *q;
  PetscInt          i,ldq;
  Vec               *W;

  PetscFunctionBegin;
  PetscCall(MatGetSize(Q,&ldq,NULL));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(VecDuplicateVecs(V->t,e-s,&W));
  for (i=s;i<e;i++) PetscCall(VecMAXPY(W[i-s],V->k-V->l,q+i*ldq+V->l,ctx->V+V->nc+V->l));
  for (i=s;i<e;i++) PetscCall(VecCopy(W[i-s],ctx->V[V->nc+i]));
  PetscCall(VecDestroyVecs(e-s,&W));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscFunctionReturn(0);
}

/*
   BVMultInPlaceHermitianTranspose_Vecs - V(:,s:e-1) = V*Q'(:,s:e-1) for regular vectors.
*/
PetscErrorCode BVMultInPlaceHermitianTranspose_Vecs(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_VECS           *ctx = (BV_VECS*)V->data;
  const PetscScalar *q;
  PetscInt          i,j,ldq,n;

  PetscFunctionBegin;
  PetscCall(MatGetSize(Q,&ldq,&n));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  /* V2 := V2*Q2' */
  PetscCall(BVMultInPlace_Vecs_Private(V,V->n,e-s,ldq,ctx->V+V->nc+s,q+s*ldq+s,PETSC_TRUE));
  /* V2 += V1*Q1' + V3*Q3' */
  for (i=s;i<e;i++) {
    for (j=V->l;j<s;j++) PetscCall(VecAXPY(ctx->V[V->nc+i],q[i+j*ldq],ctx->V[V->nc+j]));
    for (j=e;j<n;j++) PetscCall(VecAXPY(ctx->V[V->nc+i],q[i+j*ldq],ctx->V[V->nc+j]));
  }
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscFunctionReturn(0);
}

PetscErrorCode BVDot_Vecs(BV X,BV Y,Mat M)
{
  BV_VECS        *x = (BV_VECS*)X->data,*y = (BV_VECS*)Y->data;
  PetscScalar    *m;
  PetscInt       j,ldm;

  PetscFunctionBegin;
  PetscCall(MatGetSize(M,&ldm,NULL));
  PetscCall(MatDenseGetArray(M,&m));
  for (j=X->l;j<X->k;j++) PetscCall(VecMDot(x->V[X->nc+j],Y->k-Y->l,y->V+Y->nc+Y->l,m+j*ldm+Y->l));
  PetscCall(MatDenseRestoreArray(M,&m));
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_Vecs(BV X,Vec y,PetscScalar *q)
{
  BV_VECS        *x = (BV_VECS*)X->data;
  Vec            z = y;
  PetscScalar    *qq=q;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  if (!q) PetscCall(VecGetArray(X->buffer,&qq));
  PetscCall(VecMDot(z,X->k-X->l,x->V+X->nc+X->l,qq));
  if (!q) PetscCall(VecRestoreArray(X->buffer,&qq));
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_Begin_Vecs(BV X,Vec y,PetscScalar *m)
{
  BV_VECS        *x = (BV_VECS*)X->data;
  Vec            z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(VecMDotBegin(z,X->k-X->l,x->V+X->nc+X->l,m));
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_End_Vecs(BV X,Vec y,PetscScalar *m)
{
  BV_VECS        *x = (BV_VECS*)X->data;

  PetscFunctionBegin;
  PetscCall(VecMDotEnd(y,X->k-X->l,x->V+X->nc+X->l,m));
  PetscFunctionReturn(0);
}

PetscErrorCode BVScale_Vecs(BV bv,PetscInt j,PetscScalar alpha)
{
  PetscInt       i;
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  if (PetscUnlikely(j<0)) {
    for (i=bv->l;i<bv->k;i++) PetscCall(VecScale(ctx->V[bv->nc+i],alpha));
  } else PetscCall(VecScale(ctx->V[bv->nc+j],alpha));
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Vecs(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscInt       i;
  PetscReal      nrm;
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  if (PetscUnlikely(j<0)) {
    PetscCheck(type==NORM_FROBENIUS,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not implemented in BVVECS");
    *val = 0.0;
    for (i=bv->l;i<bv->k;i++) {
      PetscCall(VecNorm(ctx->V[bv->nc+i],NORM_2,&nrm));
      *val += nrm*nrm;
    }
    *val = PetscSqrtReal(*val);
  } else PetscCall(VecNorm(ctx->V[bv->nc+j],type,val));
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Begin_Vecs(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  PetscCheck(j>=0,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not implemented in BVVECS");
  PetscCall(VecNormBegin(ctx->V[bv->nc+j],type,val));
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_End_Vecs(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  PetscCheck(j>=0,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not implemented in BVVECS");
  PetscCall(VecNormEnd(ctx->V[bv->nc+j],type,val));
  PetscFunctionReturn(0);
}

PetscErrorCode BVMatMult_Vecs(BV V,Mat A,BV W)
{
  BV_VECS        *v = (BV_VECS*)V->data,*w = (BV_VECS*)W->data;
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
  PetscFunctionReturn(0);
}

PetscErrorCode BVCopy_Vecs(BV V,BV W)
{
  BV_VECS        *v = (BV_VECS*)V->data,*w = (BV_VECS*)W->data;
  PetscInt       j;

  PetscFunctionBegin;
  for (j=0;j<V->k-V->l;j++) PetscCall(VecCopy(v->V[V->nc+V->l+j],w->V[W->nc+W->l+j]));
  PetscFunctionReturn(0);
}

PetscErrorCode BVCopyColumn_Vecs(BV V,PetscInt j,PetscInt i)
{
  BV_VECS        *v = (BV_VECS*)V->data;

  PetscFunctionBegin;
  PetscCall(VecCopy(v->V[V->nc+j],v->V[V->nc+i]));
  PetscFunctionReturn(0);
}

PetscErrorCode BVResize_Vecs(BV bv,PetscInt m,PetscBool copy)
{
  BV_VECS        *ctx = (BV_VECS*)bv->data;
  Vec            *newV;
  PetscInt       j;
  char           str[50];

  PetscFunctionBegin;
  PetscCall(VecDuplicateVecs(bv->t,m,&newV));
  PetscCall(PetscLogObjectParents(bv,m,newV));
  if (((PetscObject)bv)->name) {
    for (j=0;j<m;j++) {
      PetscCall(PetscSNPrintf(str,sizeof(str),"%s_%" PetscInt_FMT,((PetscObject)bv)->name,j));
      PetscCall(PetscObjectSetName((PetscObject)newV[j],str));
    }
  }
  if (copy) {
    for (j=0;j<PetscMin(m,bv->m);j++) PetscCall(VecCopy(ctx->V[j],newV[j]));
  }
  PetscCall(VecDestroyVecs(bv->m,&ctx->V));
  ctx->V = newV;
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetColumn_Vecs(BV bv,PetscInt j,Vec *v)
{
  BV_VECS  *ctx = (BV_VECS*)bv->data;
  PetscInt l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  bv->cv[l] = ctx->V[bv->nc+j];
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetArray_Vecs(BV bv,PetscScalar **a)
{
  BV_VECS           *ctx = (BV_VECS*)bv->data;
  PetscInt          j;
  const PetscScalar *p;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1((bv->nc+bv->m)*bv->n,a));
  for (j=0;j<bv->nc+bv->m;j++) {
    PetscCall(VecGetArrayRead(ctx->V[j],&p));
    PetscCall(PetscArraycpy(*a+j*bv->n,p,bv->n));
    PetscCall(VecRestoreArrayRead(ctx->V[j],&p));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreArray_Vecs(BV bv,PetscScalar **a)
{
  BV_VECS        *ctx = (BV_VECS*)bv->data;
  PetscInt       j;
  PetscScalar    *p;

  PetscFunctionBegin;
  for (j=0;j<bv->nc+bv->m;j++) {
    PetscCall(VecGetArray(ctx->V[j],&p));
    PetscCall(PetscArraycpy(p,*a+j*bv->n,bv->n));
    PetscCall(VecRestoreArray(ctx->V[j],&p));
  }
  PetscCall(PetscFree(*a));
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetArrayRead_Vecs(BV bv,const PetscScalar **a)
{
  BV_VECS           *ctx = (BV_VECS*)bv->data;
  PetscInt          j;
  const PetscScalar *p;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1((bv->nc+bv->m)*bv->n,(PetscScalar**)a));
  for (j=0;j<bv->nc+bv->m;j++) {
    PetscCall(VecGetArrayRead(ctx->V[j],&p));
    PetscCall(PetscArraycpy((PetscScalar*)*a+j*bv->n,p,bv->n));
    PetscCall(VecRestoreArrayRead(ctx->V[j],&p));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreArrayRead_Vecs(BV bv,const PetscScalar **a)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(*a));
  PetscFunctionReturn(0);
}

/*
   Sets the value of vmip flag and resets ops->multinplace accordingly
 */
static inline PetscErrorCode BVVecsSetVmip(BV bv,PetscInt vmip)
{
  typedef PetscErrorCode (*fmultinplace)(BV,Mat,PetscInt,PetscInt);
  fmultinplace multinplace[2] = {BVMultInPlace_Vecs_ME, BVMultInPlace_Vecs_Alloc};
  BV_VECS      *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  ctx->vmip            = vmip;
  bv->ops->multinplace = multinplace[vmip];
  PetscFunctionReturn(0);
}

PetscErrorCode BVSetFromOptions_Vecs(PetscOptionItems *PetscOptionsObject,BV bv)
{
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  PetscCall(PetscOptionsHead(PetscOptionsObject,"BV Vecs Options"));

    PetscCall(PetscOptionsRangeInt("-bv_vecs_vmip","Version of BVMultInPlace operation","",ctx->vmip,&ctx->vmip,NULL,0,1));
    PetscCall(BVVecsSetVmip(bv,ctx->vmip));

  PetscCall(PetscOptionsTail());
  PetscFunctionReturn(0);
}

PetscErrorCode BVView_Vecs(BV bv,PetscViewer viewer)
{
  BV_VECS           *ctx = (BV_VECS*)bv->data;
  PetscInt          j;
  PetscViewerFormat format;
  PetscBool         isascii,ismatlab=PETSC_FALSE;
  const char        *bvname,*name;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
    if (format == PETSC_VIEWER_ASCII_MATLAB) ismatlab = PETSC_TRUE;
  }
  if (ismatlab) {
    PetscCall(PetscObjectGetName((PetscObject)bv,&bvname));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%s=[];\n",bvname));
  }
  for (j=bv->nc;j<bv->nc+bv->m;j++) {
    PetscCall(VecView(ctx->V[j],viewer));
    if (ismatlab) {
      PetscCall(PetscObjectGetName((PetscObject)ctx->V[j],&name));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%s=[%s,%s];clear %s\n",bvname,bvname,name,name));
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVDestroy_Vecs(BV bv)
{
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  if (!bv->issplit) PetscCall(VecDestroyVecs(bv->nc+bv->m,&ctx->V));
  PetscCall(PetscFree(bv->data));
  PetscFunctionReturn(0);
}

PetscErrorCode BVDuplicate_Vecs(BV V,BV W)
{
  BV_VECS        *ctx = (BV_VECS*)V->data;

  PetscFunctionBegin;
  PetscCall(BVVecsSetVmip(W,ctx->vmip));
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode BVCreate_Vecs(BV bv)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx;
  PetscInt       j,lsplit;
  PetscBool      isgpu;
  char           str[50];
  BV             parent;
  Vec            *Vpar;

  PetscFunctionBegin;
  PetscCall(PetscNewLog(bv,&ctx));
  bv->data = (void*)ctx;

  if (PetscUnlikely(bv->issplit)) {
    /* split BV: share the Vecs of the parent BV */
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    Vpar   = ((BV_VECS*)parent->data)->V;
    ctx->V = (bv->issplit==1)? Vpar: Vpar+lsplit;
  } else {
    /* regular BV: create array of Vecs to store the BV columns */
    PetscCall(VecDuplicateVecs(bv->t,bv->m,&ctx->V));
    PetscCall(PetscLogObjectParents(bv,bv->m,ctx->V));
    if (((PetscObject)bv)->name) {
      for (j=0;j<bv->m;j++) {
        PetscCall(PetscSNPrintf(str,sizeof(str),"%s_%" PetscInt_FMT,((PetscObject)bv)->name,j));
        PetscCall(PetscObjectSetName((PetscObject)ctx->V[j],str));
      }
    }
  }

  if (PetscUnlikely(bv->Acreate)) {
    for (j=0;j<bv->m;j++) PetscCall(MatGetColumnVector(bv->Acreate,ctx->V[j],j));
    PetscCall(MatDestroy(&bv->Acreate));
  }

  /* Default version of BVMultInPlace */
  PetscCall(PetscObjectTypeCompareAny((PetscObject)bv->t,&isgpu,VECSEQCUDA,VECMPICUDA,""));
  ctx->vmip = isgpu? 1: 0;

  /* Default BVMatMult method */
  bv->vmm = BV_MATMULT_VECS;

  /* Deferred call to setfromoptions */
  if (bv->defersfo) {
    ierr = PetscObjectOptionsBegin((PetscObject)bv);PetscCall(ierr);
    PetscCall(BVSetFromOptions_Vecs(PetscOptionsObject,bv));
    ierr = PetscOptionsEnd();PetscCall(ierr);
  }
  PetscCall(BVVecsSetVmip(bv,ctx->vmip));

  bv->ops->mult             = BVMult_Vecs;
  bv->ops->multvec          = BVMultVec_Vecs;
  bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Vecs;
  bv->ops->dot              = BVDot_Vecs;
  bv->ops->dotvec           = BVDotVec_Vecs;
  bv->ops->dotvec_begin     = BVDotVec_Begin_Vecs;
  bv->ops->dotvec_end       = BVDotVec_End_Vecs;
  bv->ops->scale            = BVScale_Vecs;
  bv->ops->norm             = BVNorm_Vecs;
  bv->ops->norm_begin       = BVNorm_Begin_Vecs;
  bv->ops->norm_end         = BVNorm_End_Vecs;
  bv->ops->matmult          = BVMatMult_Vecs;
  bv->ops->copy             = BVCopy_Vecs;
  bv->ops->copycolumn       = BVCopyColumn_Vecs;
  bv->ops->resize           = BVResize_Vecs;
  bv->ops->getcolumn        = BVGetColumn_Vecs;
  bv->ops->getarray         = BVGetArray_Vecs;
  bv->ops->restorearray     = BVRestoreArray_Vecs;
  bv->ops->getarrayread     = BVGetArrayRead_Vecs;
  bv->ops->restorearrayread = BVRestoreArrayRead_Vecs;
  bv->ops->destroy          = BVDestroy_Vecs;
  bv->ops->duplicate        = BVDuplicate_Vecs;
  bv->ops->setfromoptions   = BVSetFromOptions_Vecs;
  bv->ops->view             = BVView_Vecs;
  PetscFunctionReturn(0);
}
