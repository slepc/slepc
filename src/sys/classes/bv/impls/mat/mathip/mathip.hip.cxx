/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented with a dense Mat (HIP version)
*/

#include <slepc/private/bvimpl.h>
#include <slepccupmblas.h>
#include "../src/sys/classes/bv/impls/mat/bvmat.h"

PetscErrorCode BVMult_Mat_HIP(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  BV_MAT            *y = (BV_MAT*)Y->data,*x = (BV_MAT*)X->data;
  const PetscScalar *d_px,*d_A,*d_B,*d_q;
  PetscScalar       *d_py,*d_C;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (!Y->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseHIPGetArrayRead(x->A,&d_px));
  if (beta==(PetscScalar)0.0) PetscCall(MatDenseHIPGetArrayWrite(y->A,&d_py));
  else PetscCall(MatDenseHIPGetArray(y->A,&d_py));
  d_A = d_px+(X->nc+X->l)*X->ld;
  d_C = d_py+(Y->nc+Y->l)*Y->ld;
  if (Q) {
    PetscCall(MatDenseGetLDA(Q,&ldq));
    PetscCall(BV_MatDenseHIPGetArrayRead(Y,Q,&d_q));
    d_B = d_q+Y->l*ldq+X->l;
    PetscCall(BVMult_BLAS_HIP(Y,Y->n,Y->k-Y->l,X->k-X->l,alpha,d_A,X->ld,d_B,ldq,beta,d_C,Y->ld));
    PetscCall(BV_MatDenseHIPRestoreArrayRead(Y,Q,&d_q));
  } else PetscCall(BVAXPY_BLAS_HIP(Y,Y->n,Y->k-Y->l,alpha,d_A,X->ld,beta,d_C,Y->ld));
  PetscCall(MatDenseHIPRestoreArrayRead(x->A,&d_px));
  if (beta==(PetscScalar)0.0) PetscCall(MatDenseHIPRestoreArrayWrite(y->A,&d_py));
  else PetscCall(MatDenseHIPRestoreArray(y->A,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVMultVec_Mat_HIP(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  BV_MAT            *x = (BV_MAT*)X->data;
  PetscScalar       *d_py,*d_q;
  const PetscScalar *d_px;

  PetscFunctionBegin;
  PetscCall(MatDenseHIPGetArrayRead(x->A,&d_px));
  if (beta==(PetscScalar)0.0) PetscCall(VecHIPGetArrayWrite(y,&d_py));
  else PetscCall(VecHIPGetArray(y,&d_py));
  if (!q) PetscCall(VecHIPGetArray(X->buffer,&d_q));
  else {
    PetscInt k=X->k-X->l;
    PetscCallHIP(hipMalloc((void**)&d_q,k*sizeof(PetscScalar)));
    PetscCallHIP(hipMemcpy(d_q,q,k*sizeof(PetscScalar),hipMemcpyHostToDevice));
    PetscCall(PetscLogCpuToGpu(k*sizeof(PetscScalar)));
  }
  PetscCall(BVMultVec_BLAS_HIP(X,X->n,X->k-X->l,alpha,d_px+(X->nc+X->l)*X->ld,X->ld,d_q,beta,d_py));
  PetscCall(MatDenseHIPRestoreArrayRead(x->A,&d_px));
  if (beta==(PetscScalar)0.0) PetscCall(VecHIPRestoreArrayWrite(y,&d_py));
  else PetscCall(VecHIPRestoreArray(y,&d_py));
  if (!q) PetscCall(VecHIPRestoreArray(X->buffer,&d_q));
  else PetscCallHIP(hipFree(d_q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVMultInPlace_Mat_HIP(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_MAT            *ctx = (BV_MAT*)V->data;
  PetscScalar       *d_pv;
  const PetscScalar *d_q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(MatDenseHIPGetArray(ctx->A,&d_pv));
  PetscCall(BV_MatDenseHIPGetArrayRead(V,Q,&d_q));
  PetscCall(BVMultInPlace_BLAS_HIP(V,V->n,V->k-V->l,s-V->l,e-V->l,d_pv+(V->nc+V->l)*V->ld,V->ld,d_q+V->l*ldq+V->l,ldq,PETSC_FALSE));
  PetscCall(BV_MatDenseHIPRestoreArrayRead(V,Q,&d_q));
  PetscCall(MatDenseHIPRestoreArray(ctx->A,&d_pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVMultInPlaceHermitianTranspose_Mat_HIP(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_MAT            *ctx = (BV_MAT*)V->data;
  PetscScalar       *d_pv;
  const PetscScalar *d_q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(MatDenseHIPGetArray(ctx->A,&d_pv));
  PetscCall(BV_MatDenseHIPGetArrayRead(V,Q,&d_q));
  PetscCall(BVMultInPlace_BLAS_HIP(V,V->n,V->k-V->l,s-V->l,e-V->l,d_pv+(V->nc+V->l)*V->ld,V->ld,d_q+V->l*ldq+V->l,ldq,PETSC_TRUE));
  PetscCall(BV_MatDenseHIPRestoreArrayRead(V,Q,&d_q));
  PetscCall(MatDenseHIPRestoreArray(ctx->A,&d_pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVDot_Mat_HIP(BV X,BV Y,Mat M)
{
  BV_MAT            *x = (BV_MAT*)X->data,*y = (BV_MAT*)Y->data;
  const PetscScalar *d_px,*d_py;
  PetscScalar       *pm;
  PetscInt          ldm;

  PetscFunctionBegin;
  PetscCall(MatDenseGetLDA(M,&ldm));
  PetscCall(MatDenseHIPGetArrayRead(x->A,&d_px));
  PetscCall(MatDenseHIPGetArrayRead(y->A,&d_py));
  PetscCall(MatDenseGetArrayWrite(M,&pm));
  PetscCall(BVDot_BLAS_HIP(X,Y->k-Y->l,X->k-X->l,X->n,d_py+(Y->nc+Y->l)*Y->ld,Y->ld,d_px+(X->nc+X->l)*X->ld,X->ld,pm+X->l*ldm+Y->l,ldm,x->mpi));
  PetscCall(MatDenseRestoreArrayWrite(M,&pm));
  PetscCall(MatDenseHIPRestoreArrayRead(x->A,&d_px));
  PetscCall(MatDenseHIPRestoreArrayRead(y->A,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVDotVec_Mat_HIP(BV X,Vec y,PetscScalar *q)
{
  BV_MAT            *x = (BV_MAT*)X->data;
  const PetscScalar *d_px,*d_py;
  Vec               z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(MatDenseHIPGetArrayRead(x->A,&d_px));
  PetscCall(VecHIPGetArrayRead(z,&d_py));
  PetscCall(BVDotVec_BLAS_HIP(X,X->n,X->k-X->l,d_px+(X->nc+X->l)*X->ld,X->ld,d_py,q,x->mpi));
  PetscCall(VecHIPRestoreArrayRead(z,&d_py));
  PetscCall(MatDenseHIPRestoreArrayRead(x->A,&d_px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVDotVec_Local_Mat_HIP(BV X,Vec y,PetscScalar *m)
{
  BV_MAT            *x = (BV_MAT*)X->data;
  const PetscScalar *d_px,*d_py;
  Vec               z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(MatDenseHIPGetArrayRead(x->A,&d_px));
  PetscCall(VecHIPGetArrayRead(z,&d_py));
  PetscCall(BVDotVec_BLAS_HIP(X,X->n,X->k-X->l,d_px+(X->nc+X->l)*X->ld,X->ld,d_py,m,PETSC_FALSE));
  PetscCall(VecHIPRestoreArrayRead(z,&d_py));
  PetscCall(MatDenseHIPRestoreArrayRead(x->A,&d_px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVScale_Mat_HIP(BV bv,PetscInt j,PetscScalar alpha)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *d_array,*d_A;
  PetscInt       n=0;

  PetscFunctionBegin;
  if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseHIPGetArray(ctx->A,&d_array));
  if (PetscUnlikely(j<0)) {
    d_A = d_array+(bv->nc+bv->l)*bv->ld;
    n = (bv->k-bv->l)*bv->ld;
  } else {
    d_A = d_array+(bv->nc+j)*bv->ld;
    n = bv->n;
  }
  PetscCall(BVScale_BLAS_HIP(bv,n,d_A,alpha));
  PetscCall(MatDenseHIPRestoreArray(ctx->A,&d_array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVNorm_Mat_HIP(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_MAT            *ctx = (BV_MAT*)bv->data;
  const PetscScalar *array,*d_array,*d_A;
  PetscInt          n=0;

  PetscFunctionBegin;
  if (!ctx->mpi && ((j<0 && type==NORM_FROBENIUS && bv->ld==bv->n) || (j>=0 && type==NORM_2))) {
    /* compute on GPU with hipBLAS - TODO: include the MPI case here */
    *val = 0.0;
    if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(MatDenseHIPGetArrayRead(ctx->A,&d_array));
    if (PetscUnlikely(j<0)) {
      d_A = d_array+(bv->nc+bv->l)*bv->ld;
      n = (bv->k-bv->l)*bv->ld;
    } else {
      d_A = d_array+(bv->nc+j)*bv->ld;
      n = bv->n;
    }
    PetscCall(BVNorm_BLAS_HIP(bv,n,d_A,val));
    PetscCall(MatDenseHIPRestoreArrayRead(ctx->A,&d_array));
  } else {
    /* compute on CPU */
    PetscCall(MatDenseGetArrayRead(ctx->A,&array));
    if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,ctx->mpi));
    else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->ld,bv->ld,type,val,ctx->mpi));
    PetscCall(MatDenseRestoreArrayRead(ctx->A,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVNorm_Local_Mat_HIP(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_MAT            *ctx = (BV_MAT*)bv->data;
  const PetscScalar *array,*d_array,*d_A;
  PetscInt          n=0;

  PetscFunctionBegin;
  if ((j<0 && type==NORM_FROBENIUS && bv->ld==bv->n) || (j>=0 && type==NORM_2)) {
    /* compute on GPU with hipBLAS */
    *val = 0.0;
    if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(MatDenseHIPGetArrayRead(ctx->A,&d_array));
    if (PetscUnlikely(j<0)) {
      d_A = d_array+(bv->nc+bv->l)*bv->ld;
      n = (bv->k-bv->l)*bv->ld;
    } else {
      d_A = d_array+(bv->nc+j)*bv->ld;
      n = bv->n;
    }
    PetscCall(BVNorm_BLAS_HIP(bv,n,d_A,val));
    PetscCall(MatDenseHIPRestoreArrayRead(ctx->A,&d_array));
  } else {
    /* compute on CPU */
    PetscCall(MatDenseGetArrayRead(ctx->A,&array));
    if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,PETSC_FALSE));
    else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->ld,bv->ld,type,val,PETSC_FALSE));
    PetscCall(MatDenseRestoreArrayRead(ctx->A,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVNormalize_Mat_HIP(BV bv,PetscScalar *eigi)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *array,*d_array,*wi=NULL;

  PetscFunctionBegin;
  if (eigi) wi = eigi+bv->l;
  if (!ctx->mpi) {
    /* compute on GPU with hipBLAS - TODO: include the MPI case here */
    if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(MatDenseHIPGetArray(ctx->A,&d_array));
    PetscCall(BVNormalize_BLAS_HIP(bv,bv->n,bv->k-bv->l,d_array+(bv->nc+bv->l)*bv->ld,bv->ld,wi));
    PetscCall(MatDenseHIPRestoreArray(ctx->A,&d_array));
  } else {
    /* compute on CPU */
    PetscCall(MatDenseGetArray(ctx->A,&array));
    PetscCall(BVNormalize_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,wi,ctx->mpi));
    PetscCall(MatDenseRestoreArray(ctx->A,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVMatMult_Mat_HIP(BV V,Mat A,BV W)
{
  BV_MAT            *v = (BV_MAT*)V->data,*w = (BV_MAT*)W->data;
  Mat               Vmat,Wmat;
  const PetscScalar *d_pv;
  PetscScalar       *d_pw;
  PetscInt          j;

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
    PetscCall(MatDenseHIPGetArrayRead(v->A,&d_pv));
    PetscCall(MatDenseHIPGetArrayWrite(w->A,&d_pw));
    for (j=0;j<V->k-V->l;j++) {
      PetscCall(VecHIPPlaceArray(V->cv[1],(PetscScalar *)d_pv+(V->nc+V->l+j)*V->ld));
      PetscCall(VecHIPPlaceArray(W->cv[1],d_pw+(W->nc+W->l+j)*W->ld));
      PetscCall(MatMult(A,V->cv[1],W->cv[1]));
      PetscCall(VecHIPResetArray(V->cv[1]));
      PetscCall(VecHIPResetArray(W->cv[1]));
    }
    PetscCall(MatDenseHIPRestoreArrayRead(v->A,&d_pv));
    PetscCall(MatDenseHIPRestoreArrayWrite(w->A,&d_pw));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVCopy_Mat_HIP(BV V,BV W)
{
  BV_MAT            *v = (BV_MAT*)V->data,*w = (BV_MAT*)W->data;
  const PetscScalar *d_pv;
  PetscScalar       *d_pw;

  PetscFunctionBegin;
  PetscCall(MatDenseHIPGetArrayRead(v->A,&d_pv));
  PetscCall(MatDenseHIPGetArray(w->A,&d_pw));
  PetscCallHIP(hipMemcpy2D(d_pw+(W->nc+W->l)*W->ld,W->ld*sizeof(PetscScalar),d_pv+(V->nc+V->l)*V->ld,V->ld*sizeof(PetscScalar),V->n*sizeof(PetscScalar),V->k-V->l,hipMemcpyDeviceToDevice));
  PetscCall(MatDenseHIPRestoreArrayRead(v->A,&d_pv));
  PetscCall(MatDenseHIPRestoreArray(w->A,&d_pw));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVCopyColumn_Mat_HIP(BV V,PetscInt j,PetscInt i)
{
  BV_MAT         *v = (BV_MAT*)V->data;
  PetscScalar    *d_pv;

  PetscFunctionBegin;
  PetscCall(MatDenseHIPGetArray(v->A,&d_pv));
  PetscCallHIP(hipMemcpy(d_pv+(V->nc+i)*V->ld,d_pv+(V->nc+j)*V->ld,V->n*sizeof(PetscScalar),hipMemcpyDeviceToDevice));
  PetscCall(MatDenseHIPRestoreArray(v->A,&d_pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVGetColumn_Mat_HIP(BV bv,PetscInt j,Vec*)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *d_pv;
  PetscInt       l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  PetscCall(MatDenseHIPGetArray(ctx->A,&d_pv));
  PetscCall(VecHIPPlaceArray(bv->cv[l],d_pv+(bv->nc+j)*bv->ld));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreColumn_Mat_HIP(BV bv,PetscInt j,Vec*)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscInt       l;

  PetscFunctionBegin;
  l = (j==bv->ci[0])? 0: 1;
  PetscCall(VecHIPResetArray(bv->cv[l]));
  PetscCall(MatDenseHIPRestoreArray(ctx->A,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreSplit_Mat_HIP(BV bv,BV *L,BV *R)
{
  Mat               A;
  const PetscScalar *d_pv;
  PetscObjectState  lstate,rstate;
  PetscBool         change=PETSC_FALSE;

  PetscFunctionBegin;
  /* force sync flag to PETSC_OFFLOAD_BOTH */
  if (L) {
    PetscCall(PetscObjectStateGet((PetscObject)*L,&lstate));
    if (lstate != bv->lstate) {
      A = ((BV_MAT*)bv->L->data)->A;
      PetscCall(MatDenseHIPGetArrayRead(A,&d_pv));
      PetscCall(MatDenseHIPRestoreArrayRead(A,&d_pv));
      change = PETSC_TRUE;
    }
  }
  if (R) {
    PetscCall(PetscObjectStateGet((PetscObject)*R,&rstate));
    if (rstate != bv->rstate) {
      A = ((BV_MAT*)bv->R->data)->A;
      PetscCall(MatDenseHIPGetArrayRead(A,&d_pv));
      PetscCall(MatDenseHIPRestoreArrayRead(A,&d_pv));
      change = PETSC_TRUE;
    }
  }
  if (change) {
    A = ((BV_MAT*)bv->data)->A;
    PetscCall(MatDenseHIPGetArray(A,(PetscScalar **)&d_pv));
    PetscCall(MatDenseHIPRestoreArray(A,(PetscScalar **)&d_pv));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreSplitRows_Mat_HIP(BV bv,IS,IS,BV *U,BV *L)
{
  Mat               A;
  const PetscScalar *d_pv;
  PetscObjectState  lstate,rstate;
  PetscBool         change=PETSC_FALSE;

  PetscFunctionBegin;
  /* force sync flag to PETSC_OFFLOAD_BOTH */
  if (U) {
    PetscCall(PetscObjectStateGet((PetscObject)*U,&rstate));
    if (rstate != bv->rstate) {
      A = ((BV_MAT*)bv->R->data)->A;
      PetscCall(MatDenseHIPGetArrayRead(A,&d_pv));
      PetscCall(MatDenseHIPRestoreArrayRead(A,&d_pv));
      change = PETSC_TRUE;
    }
  }
  if (L) {
    PetscCall(PetscObjectStateGet((PetscObject)*L,&lstate));
    if (lstate != bv->lstate) {
      A = ((BV_MAT*)bv->L->data)->A;
      PetscCall(MatDenseHIPGetArrayRead(A,&d_pv));
      PetscCall(MatDenseHIPRestoreArrayRead(A,&d_pv));
      change = PETSC_TRUE;
    }
  }
  if (change) {
    A = ((BV_MAT*)bv->data)->A;
    PetscCall(MatDenseHIPGetArray(A,(PetscScalar **)&d_pv));
    PetscCall(MatDenseHIPRestoreArray(A,(PetscScalar **)&d_pv));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVGetMat_Mat_HIP(BV bv,Mat *A)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *vv,*aa;
  PetscBool      create=PETSC_FALSE;
  PetscInt       m,cols;

  PetscFunctionBegin;
  m = bv->k-bv->l;
  if (!bv->Aget) create=PETSC_TRUE;
  else {
    PetscCall(MatDenseHIPGetArray(bv->Aget,&aa));
    PetscCheck(!aa,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"BVGetMat already called on this BV");
    PetscCall(MatGetSize(bv->Aget,NULL,&cols));
    if (cols!=m) {
      PetscCall(MatDestroy(&bv->Aget));
      create=PETSC_TRUE;
    }
  }
  PetscCall(MatDenseHIPGetArray(ctx->A,&vv));
  if (create) {
    PetscCall(MatCreateDenseFromVecType(PetscObjectComm((PetscObject)bv),bv->vtype,bv->n,PETSC_DECIDE,bv->N,m,bv->ld,vv,&bv->Aget)); /* pass a pointer to avoid allocation of storage */
    PetscCall(MatDenseHIPReplaceArray(bv->Aget,NULL));  /* replace with a null pointer, the value after BVRestoreMat */
  }
  PetscCall(MatDenseHIPPlaceArray(bv->Aget,vv+(bv->nc+bv->l)*bv->ld));  /* set the actual pointer */
  *A = bv->Aget;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreMat_Mat_HIP(BV bv,Mat *A)
{
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *vv,*aa;

  PetscFunctionBegin;
  PetscCall(MatDenseHIPGetArray(bv->Aget,&aa));
  vv = aa-(bv->nc+bv->l)*bv->ld;
  PetscCall(MatDenseHIPResetArray(bv->Aget));
  PetscCall(MatDenseHIPRestoreArray(ctx->A,&vv));
  *A = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}
