/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented as a single Vec (HIP version)
*/

#include <slepc/private/bvimpl.h>
#include <slepccupmblas.h>
#include "../src/sys/classes/bv/impls/svec/svec.h"

PetscErrorCode BVMult_Svec_HIP(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  BV_SVEC           *y = (BV_SVEC*)Y->data,*x = (BV_SVEC*)X->data;
  const PetscScalar *d_px,*d_A,*d_B,*d_q;
  PetscScalar       *d_py,*d_C;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (!Y->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(VecHIPGetArrayRead(x->v,&d_px));
  if (beta==(PetscScalar)0.0) PetscCall(VecHIPGetArrayWrite(y->v,&d_py));
  else PetscCall(VecHIPGetArray(y->v,&d_py));
  d_A = d_px+(X->nc+X->l)*X->ld;
  d_C = d_py+(Y->nc+Y->l)*Y->ld;
  if (Q) {
    PetscCall(MatDenseGetLDA(Q,&ldq));
    PetscCall(BV_MatDenseHIPGetArrayRead(Y,Q,&d_q));
    d_B = d_q+Y->l*ldq+X->l;
    PetscCall(BVMult_BLAS_HIP(Y,Y->n,Y->k-Y->l,X->k-X->l,alpha,d_A,X->ld,d_B,ldq,beta,d_C,Y->ld));
    PetscCall(BV_MatDenseHIPRestoreArrayRead(Y,Q,&d_q));
  } else PetscCall(BVAXPY_BLAS_HIP(Y,Y->n,Y->k-Y->l,alpha,d_A,X->ld,beta,d_C,Y->ld));
  PetscCall(VecHIPRestoreArrayRead(x->v,&d_px));
  if (beta==(PetscScalar)0.0) PetscCall(VecHIPRestoreArrayWrite(y->v,&d_py));
  else PetscCall(VecHIPRestoreArray(y->v,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVMultVec_Svec_HIP(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  BV_SVEC           *x = (BV_SVEC*)X->data;
  PetscScalar       *d_py,*d_q;
  const PetscScalar *d_px;

  PetscFunctionBegin;
  PetscCall(VecHIPGetArrayRead(x->v,&d_px));
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
  PetscCall(VecHIPRestoreArrayRead(x->v,&d_px));
  if (beta==(PetscScalar)0.0) PetscCall(VecHIPRestoreArrayWrite(y,&d_py));
  else PetscCall(VecHIPRestoreArray(y,&d_py));
  if (!q) PetscCall(VecHIPRestoreArray(X->buffer,&d_q));
  else PetscCallHIP(hipFree(d_q));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVMultInPlace_Svec_HIP(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_SVEC           *ctx = (BV_SVEC*)V->data;
  PetscScalar       *d_pv;
  const PetscScalar *d_q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(VecHIPGetArray(ctx->v,&d_pv));
  PetscCall(BV_MatDenseHIPGetArrayRead(V,Q,&d_q));
  PetscCall(BVMultInPlace_BLAS_HIP(V,V->n,V->k-V->l,s-V->l,e-V->l,d_pv+(V->nc+V->l)*V->ld,V->ld,d_q+V->l*ldq+V->l,ldq,PETSC_FALSE));
  PetscCall(BV_MatDenseHIPRestoreArrayRead(V,Q,&d_q));
  PetscCall(VecHIPRestoreArray(ctx->v,&d_pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVMultInPlaceHermitianTranspose_Svec_HIP(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_SVEC           *ctx = (BV_SVEC*)V->data;
  PetscScalar       *d_pv;
  const PetscScalar *d_q;
  PetscInt          ldq;

  PetscFunctionBegin;
  if (s>=e || !V->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(VecHIPGetArray(ctx->v,&d_pv));
  PetscCall(BV_MatDenseHIPGetArrayRead(V,Q,&d_q));
  PetscCall(BVMultInPlace_BLAS_HIP(V,V->n,V->k-V->l,s-V->l,e-V->l,d_pv+(V->nc+V->l)*V->ld,V->ld,d_q+V->l*ldq+V->l,ldq,PETSC_TRUE));
  PetscCall(BV_MatDenseHIPRestoreArrayRead(V,Q,&d_q));
  PetscCall(VecHIPRestoreArray(ctx->v,&d_pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVDot_Svec_HIP(BV X,BV Y,Mat M)
{
  BV_SVEC           *x = (BV_SVEC*)X->data,*y = (BV_SVEC*)Y->data;
  const PetscScalar *d_px,*d_py;
  PetscScalar       *pm;
  PetscInt          ldm;

  PetscFunctionBegin;
  PetscCall(MatDenseGetLDA(M,&ldm));
  PetscCall(VecHIPGetArrayRead(x->v,&d_px));
  PetscCall(VecHIPGetArrayRead(y->v,&d_py));
  PetscCall(MatDenseGetArrayWrite(M,&pm));
  PetscCall(BVDot_BLAS_HIP(X,Y->k-Y->l,X->k-X->l,X->n,d_py+(Y->nc+Y->l)*Y->ld,Y->ld,d_px+(X->nc+X->l)*X->ld,X->ld,pm+X->l*ldm+Y->l,ldm,x->mpi));
  PetscCall(MatDenseRestoreArrayWrite(M,&pm));
  PetscCall(VecHIPRestoreArrayRead(x->v,&d_px));
  PetscCall(VecHIPRestoreArrayRead(y->v,&d_py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVDotVec_Svec_HIP(BV X,Vec y,PetscScalar *q)
{
  BV_SVEC           *x = (BV_SVEC*)X->data;
  const PetscScalar *d_px,*d_py;
  Vec               z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(VecHIPGetArrayRead(x->v,&d_px));
  PetscCall(VecHIPGetArrayRead(z,&d_py));
  PetscCall(BVDotVec_BLAS_HIP(X,X->n,X->k-X->l,d_px+(X->nc+X->l)*X->ld,X->ld,d_py,q,x->mpi));
  PetscCall(VecHIPRestoreArrayRead(z,&d_py));
  PetscCall(VecHIPRestoreArrayRead(x->v,&d_px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVDotVec_Local_Svec_HIP(BV X,Vec y,PetscScalar *m)
{
  BV_SVEC           *x = (BV_SVEC*)X->data;
  const PetscScalar *d_px,*d_py;
  Vec               z = y;

  PetscFunctionBegin;
  if (PetscUnlikely(X->matrix)) {
    PetscCall(BV_IPMatMult(X,y));
    z = X->Bx;
  }
  PetscCall(VecHIPGetArrayRead(x->v,&d_px));
  PetscCall(VecHIPGetArrayRead(z,&d_py));
  PetscCall(BVDotVec_BLAS_HIP(X,X->n,X->k-X->l,d_px+(X->nc+X->l)*X->ld,X->ld,d_py,m,PETSC_FALSE));
  PetscCall(VecHIPRestoreArrayRead(z,&d_py));
  PetscCall(VecHIPRestoreArrayRead(x->v,&d_px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVScale_Svec_HIP(BV bv,PetscInt j,PetscScalar alpha)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *d_array,*d_A;
  PetscInt       n=0;

  PetscFunctionBegin;
  if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(VecHIPGetArray(ctx->v,&d_array));
  if (PetscUnlikely(j<0)) {
    d_A = d_array+(bv->nc+bv->l)*bv->ld;
    n = (bv->k-bv->l)*bv->ld;
  } else {
    d_A = d_array+(bv->nc+j)*bv->ld;
    n = bv->n;
  }
  PetscCall(BVScale_BLAS_HIP(bv,n,d_A,alpha));
  PetscCall(VecHIPRestoreArray(ctx->v,&d_array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVNorm_Svec_HIP(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_SVEC           *ctx = (BV_SVEC*)bv->data;
  const PetscScalar *array,*d_array,*d_A;
  PetscInt          n=0;

  PetscFunctionBegin;
  if (!ctx->mpi && ((j<0 && type==NORM_FROBENIUS && bv->ld==bv->n) || (j>=0 && type==NORM_2))) {
    /* compute on GPU with hipBLAS - TODO: include the MPI case here */
    *val = 0.0;
    if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(VecHIPGetArrayRead(ctx->v,&d_array));
    if (PetscUnlikely(j<0)) {
      d_A = d_array+(bv->nc+bv->l)*bv->ld;
      n = (bv->k-bv->l)*bv->ld;
    } else {
      d_A = d_array+(bv->nc+j)*bv->ld;
      n = bv->n;
    }
    PetscCall(BVNorm_BLAS_HIP(bv,n,d_A,val));
    PetscCall(VecHIPRestoreArrayRead(ctx->v,&d_array));
  } else {
    /* compute on CPU */
    PetscCall(VecGetArrayRead(ctx->v,&array));
    if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,ctx->mpi));
    else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->ld,bv->ld,type,val,ctx->mpi));
    PetscCall(VecRestoreArrayRead(ctx->v,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVNorm_Local_Svec_HIP(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_SVEC           *ctx = (BV_SVEC*)bv->data;
  const PetscScalar *array,*d_array,*d_A;
  PetscInt          n=0;

  PetscFunctionBegin;
  if ((j<0 && type==NORM_FROBENIUS && bv->ld==bv->n) || (j>=0 && type==NORM_2)) {
    /* compute on GPU with hipBLAS */
    *val = 0.0;
    if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(VecHIPGetArrayRead(ctx->v,&d_array));
    if (PetscUnlikely(j<0)) {
      d_A = d_array+(bv->nc+bv->l)*bv->ld;
      n = (bv->k-bv->l)*bv->ld;
    } else {
      d_A = d_array+(bv->nc+j)*bv->ld;
      n = bv->n;
    }
    PetscCall(BVNorm_BLAS_HIP(bv,n,d_A,val));
    PetscCall(VecHIPRestoreArrayRead(ctx->v,&d_array));
  } else {
    /* compute on CPU */
    PetscCall(VecGetArrayRead(ctx->v,&array));
    if (PetscUnlikely(j<0)) PetscCall(BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,type,val,PETSC_FALSE));
    else PetscCall(BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->ld,bv->ld,type,val,PETSC_FALSE));
    PetscCall(VecRestoreArrayRead(ctx->v,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVNormalize_Svec_HIP(BV bv,PetscScalar *eigi)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *array,*d_array,*wi=NULL;

  PetscFunctionBegin;
  if (eigi) wi = eigi+bv->l;
  if (!ctx->mpi) {
    /* compute on GPU with hipBLAS - TODO: include the MPI case here */
    if (!bv->n) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(VecHIPGetArray(ctx->v,&d_array));
    PetscCall(BVNormalize_BLAS_HIP(bv,bv->n,bv->k-bv->l,d_array+(bv->nc+bv->l)*bv->ld,bv->ld,wi));
    PetscCall(VecHIPRestoreArray(ctx->v,&d_array));
  } else {
    /* compute on CPU */
    PetscCall(VecGetArray(ctx->v,&array));
    PetscCall(BVNormalize_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->ld,bv->ld,wi,ctx->mpi));
    PetscCall(VecRestoreArray(ctx->v,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVMatMult_Svec_HIP(BV V,Mat A,BV W)
{
  BV_SVEC           *v = (BV_SVEC*)V->data,*w = (BV_SVEC*)W->data;
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
    PetscCall(VecHIPGetArrayRead(v->v,&d_pv));
    PetscCall(VecHIPGetArrayWrite(w->v,&d_pw));
    for (j=0;j<V->k-V->l;j++) {
      PetscCall(VecHIPPlaceArray(V->cv[1],(PetscScalar *)d_pv+(V->nc+V->l+j)*V->ld));
      PetscCall(VecHIPPlaceArray(W->cv[1],d_pw+(W->nc+W->l+j)*W->ld));
      PetscCall(MatMult(A,V->cv[1],W->cv[1]));
      PetscCall(VecHIPResetArray(V->cv[1]));
      PetscCall(VecHIPResetArray(W->cv[1]));
    }
    PetscCall(VecHIPRestoreArrayRead(v->v,&d_pv));
    PetscCall(VecHIPRestoreArrayWrite(w->v,&d_pw));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVCopy_Svec_HIP(BV V,BV W)
{
  BV_SVEC           *v = (BV_SVEC*)V->data,*w = (BV_SVEC*)W->data;
  const PetscScalar *d_pv;
  PetscScalar       *d_pw;

  PetscFunctionBegin;
  PetscCall(VecHIPGetArrayRead(v->v,&d_pv));
  PetscCall(VecHIPGetArray(w->v,&d_pw));
  PetscCallHIP(hipMemcpy2D(d_pw+(W->nc+W->l)*W->ld,W->ld*sizeof(PetscScalar),d_pv+(V->nc+V->l)*V->ld,V->ld*sizeof(PetscScalar),V->n*sizeof(PetscScalar),V->k-V->l,hipMemcpyDeviceToDevice));
  PetscCall(VecHIPRestoreArrayRead(v->v,&d_pv));
  PetscCall(VecHIPRestoreArray(w->v,&d_pw));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVCopyColumn_Svec_HIP(BV V,PetscInt j,PetscInt i)
{
  BV_SVEC        *v = (BV_SVEC*)V->data;
  PetscScalar    *d_pv;

  PetscFunctionBegin;
  PetscCall(VecHIPGetArray(v->v,&d_pv));
  PetscCallHIP(hipMemcpy(d_pv+(V->nc+i)*V->ld,d_pv+(V->nc+j)*V->ld,V->n*sizeof(PetscScalar),hipMemcpyDeviceToDevice));
  PetscCall(VecHIPRestoreArray(v->v,&d_pv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVResize_Svec_HIP(BV bv,PetscInt m,PetscBool copy)
{
  BV_SVEC           *ctx = (BV_SVEC*)bv->data;
  const PetscScalar *d_pv;
  PetscScalar       *d_pnew;
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
    PetscCall(VecHIPGetArrayRead(ctx->v,&d_pv));
    PetscCall(VecHIPGetArrayWrite(vnew,&d_pnew));
    PetscCallHIP(hipMemcpy(d_pnew,d_pv,PetscMin(m,bv->m)*bv->ld*sizeof(PetscScalar),hipMemcpyDeviceToDevice));
    PetscCall(VecHIPRestoreArrayRead(ctx->v,&d_pv));
    PetscCall(VecHIPRestoreArrayWrite(vnew,&d_pnew));
  }
  PetscCall(VecDestroy(&ctx->v));
  ctx->v = vnew;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVGetColumn_Svec_HIP(BV bv,PetscInt j,Vec*)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *d_pv;
  PetscInt       l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  PetscCall(VecHIPGetArray(ctx->v,&d_pv));
  PetscCall(VecHIPPlaceArray(bv->cv[l],d_pv+(bv->nc+j)*bv->ld));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreColumn_Svec_HIP(BV bv,PetscInt j,Vec*)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscInt       l;

  PetscFunctionBegin;
  l = (j==bv->ci[0])? 0: 1;
  PetscCall(VecHIPResetArray(bv->cv[l]));
  PetscCall(VecHIPRestoreArray(ctx->v,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreSplit_Svec_HIP(BV bv,BV *L,BV *R)
{
  Vec               v;
  const PetscScalar *d_pv;
  PetscObjectState  lstate,rstate;
  PetscBool         change=PETSC_FALSE;

  PetscFunctionBegin;
  /* force sync flag to PETSC_OFFLOAD_BOTH */
  if (L) {
    PetscCall(PetscObjectStateGet((PetscObject)*L,&lstate));
    if (lstate != bv->lstate) {
      v = ((BV_SVEC*)bv->L->data)->v;
      PetscCall(VecHIPGetArrayRead(v,&d_pv));
      PetscCall(VecHIPRestoreArrayRead(v,&d_pv));
      change = PETSC_TRUE;
    }
  }
  if (R) {
    PetscCall(PetscObjectStateGet((PetscObject)*R,&rstate));
    if (rstate != bv->rstate) {
      v = ((BV_SVEC*)bv->R->data)->v;
      PetscCall(VecHIPGetArrayRead(v,&d_pv));
      PetscCall(VecHIPRestoreArrayRead(v,&d_pv));
      change = PETSC_TRUE;
    }
  }
  if (change) {
    v = ((BV_SVEC*)bv->data)->v;
    PetscCall(VecHIPGetArray(v,(PetscScalar **)&d_pv));
    PetscCall(VecHIPRestoreArray(v,(PetscScalar **)&d_pv));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreSplitRows_Svec_HIP(BV bv,IS,IS,BV *U,BV *L)
{
  Vec               v;
  const PetscScalar *d_pv;
  PetscObjectState  lstate,rstate;
  PetscBool         change=PETSC_FALSE;

  PetscFunctionBegin;
  /* force sync flag to PETSC_OFFLOAD_BOTH */
  if (U) {
    PetscCall(PetscObjectStateGet((PetscObject)*U,&rstate));
    if (rstate != bv->rstate) {
      v = ((BV_SVEC*)bv->R->data)->v;
      PetscCall(VecHIPGetArrayRead(v,&d_pv));
      PetscCall(VecHIPRestoreArrayRead(v,&d_pv));
      change = PETSC_TRUE;
    }
  }
  if (L) {
    PetscCall(PetscObjectStateGet((PetscObject)*L,&lstate));
    if (lstate != bv->lstate) {
      v = ((BV_SVEC*)bv->L->data)->v;
      PetscCall(VecHIPGetArrayRead(v,&d_pv));
      PetscCall(VecHIPRestoreArrayRead(v,&d_pv));
      change = PETSC_TRUE;
    }
  }
  if (change) {
    v = ((BV_SVEC*)bv->data)->v;
    PetscCall(VecHIPGetArray(v,(PetscScalar **)&d_pv));
    PetscCall(VecHIPRestoreArray(v,(PetscScalar **)&d_pv));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVGetMat_Svec_HIP(BV bv,Mat *A)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
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
  PetscCall(VecHIPGetArray(ctx->v,&vv));
  if (create) {
    PetscCall(MatCreateDenseFromVecType(PetscObjectComm((PetscObject)bv),bv->vtype,bv->n,PETSC_DECIDE,bv->N,m,bv->ld,vv,&bv->Aget)); /* pass a pointer to avoid allocation of storage */
    PetscCall(MatDenseHIPReplaceArray(bv->Aget,NULL));  /* replace with a null pointer, the value after BVRestoreMat */
  }
  PetscCall(MatDenseHIPPlaceArray(bv->Aget,vv+(bv->nc+bv->l)*bv->ld));  /* set the actual pointer */
  *A = bv->Aget;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreMat_Svec_HIP(BV bv,Mat *A)
{
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *vv,*aa;

  PetscFunctionBegin;
  PetscCall(MatDenseHIPGetArray(bv->Aget,&aa));
  vv = aa-(bv->nc+bv->l)*bv->ld;
  PetscCall(MatDenseHIPResetArray(bv->Aget));
  PetscCall(VecHIPRestoreArray(ctx->v,&vv));
  *A = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}
