/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Tensor BV that is represented in compact form as V = (I otimes U) S
*/

#include <slepc/private/bvimpl.h>      /*I "slepcbv.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  BV          U;        /* first factor */
  Mat         S;        /* second factor */
  PetscScalar *qB;      /* auxiliary matrix used in non-standard inner products */
  PetscScalar *sw;      /* work space */
  PetscInt    d;        /* degree of the tensor BV */
  PetscInt    ld;       /* leading dimension of a single block in S */
  PetscInt    puk;      /* copy of the k value */
  Vec         u;        /* auxiliary work vector */
} BV_TENSOR;

static PetscErrorCode BVMultInPlace_Tensor(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_TENSOR         *ctx = (BV_TENSOR*)V->data;
  PetscScalar       *pS;
  const PetscScalar *q;
  PetscInt          ldq,lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(MatDenseGetArray(ctx->S,&pS));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(BVMultInPlace_BLAS_Private(V,lds,V->k-V->l,s-V->l,e-V->l,pS+(V->nc+V->l)*lds,lds,q+V->l*ldq+V->l,ldq,PETSC_FALSE));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscCall(MatDenseRestoreArray(ctx->S,&pS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVMultInPlaceHermitianTranspose_Tensor(BV V,Mat Q,PetscInt s,PetscInt e)
{
  BV_TENSOR         *ctx = (BV_TENSOR*)V->data;
  PetscScalar       *pS;
  const PetscScalar *q;
  PetscInt          ldq,lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(MatDenseGetArray(ctx->S,&pS));
  PetscCall(MatDenseGetArrayRead(Q,&q));
  PetscCall(BVMultInPlace_BLAS_Private(V,lds,V->k-V->l,s-V->l,e-V->l,pS+(V->nc+V->l)*lds,lds,q+V->l*ldq+V->l,ldq,PETSC_TRUE));
  PetscCall(MatDenseRestoreArrayRead(Q,&q));
  PetscCall(MatDenseRestoreArray(ctx->S,&pS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDot_Tensor(BV X,BV Y,Mat M)
{
  BV_TENSOR         *x = (BV_TENSOR*)X->data,*y = (BV_TENSOR*)Y->data;
  PetscScalar       *m;
  const PetscScalar *px,*py;
  PetscInt          ldm,lds = x->ld*x->d;

  PetscFunctionBegin;
  PetscCheck(x->U==y->U,PetscObjectComm((PetscObject)X),PETSC_ERR_SUP,"BVDot() in BVTENSOR requires that both operands have the same U factor");
  PetscCheck(lds==y->ld*y->d,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Mismatching dimensions ld*d %" PetscInt_FMT " %" PetscInt_FMT,lds,y->ld*y->d);
  PetscCall(MatDenseGetLDA(M,&ldm));
  PetscCall(MatDenseGetArrayRead(x->S,&px));
  PetscCall(MatDenseGetArrayRead(y->S,&py));
  PetscCall(MatDenseGetArray(M,&m));
  PetscCall(BVDot_BLAS_Private(X,Y->k-Y->l,X->k-X->l,lds,py+(Y->nc+Y->l)*lds,lds,px+(X->nc+X->l)*lds,lds,m+X->l*ldm+Y->l,ldm,PETSC_FALSE));
  PetscCall(MatDenseRestoreArray(M,&m));
  PetscCall(MatDenseRestoreArrayRead(x->S,&px));
  PetscCall(MatDenseRestoreArrayRead(y->S,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVScale_Tensor(BV bv,PetscInt j,PetscScalar alpha)
{
  BV_TENSOR      *ctx = (BV_TENSOR*)bv->data;
  PetscScalar    *pS;
  PetscInt       lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ctx->S,&pS));
  if (PetscUnlikely(j<0)) PetscCall(BVScale_BLAS_Private(bv,(bv->k-bv->l)*lds,pS+(bv->nc+bv->l)*lds,alpha));
  else PetscCall(BVScale_BLAS_Private(bv,lds,pS+(bv->nc+j)*lds,alpha));
  PetscCall(MatDenseRestoreArray(ctx->S,&pS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVNorm_Tensor(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  BV_TENSOR         *ctx = (BV_TENSOR*)bv->data;
  const PetscScalar *pS;
  PetscInt          lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(ctx->S,&pS));
  if (j<0) PetscCall(BVNorm_LAPACK_Private(bv,lds,bv->k-bv->l,pS+(bv->nc+bv->l)*lds,lds,type,val,PETSC_FALSE));
  else PetscCall(BVNorm_LAPACK_Private(bv,lds,1,pS+(bv->nc+j)*lds,lds,type,val,PETSC_FALSE));
  PetscCall(MatDenseRestoreArrayRead(ctx->S,&pS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVCopyColumn_Tensor(BV V,PetscInt j,PetscInt i)
{
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscScalar    *pS;
  PetscInt       lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ctx->S,&pS));
  PetscCall(PetscArraycpy(pS+(V->nc+i)*lds,pS+(V->nc+j)*lds,lds));
  PetscCall(MatDenseRestoreArray(ctx->S,&pS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVTensorNormColumn(BV bv,PetscInt j,PetscReal *norm)
{
  BV_TENSOR         *ctx = (BV_TENSOR*)bv->data;
  PetscBLASInt      one=1,lds_;
  PetscScalar       sone=1.0,szero=0.0,*x,dot;
  const PetscScalar *S;
  PetscReal         alpha=1.0,scale=0.0,aval;
  PetscInt          i,lds,ld=ctx->ld;

  PetscFunctionBegin;
  lds = ld*ctx->d;
  PetscCall(MatDenseGetArrayRead(ctx->S,&S));
  PetscCall(PetscBLASIntCast(lds,&lds_));
  if (PetscUnlikely(ctx->qB)) {
    x = ctx->sw;
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&lds_,&sone,ctx->qB,&lds_,S+j*lds,&one,&szero,x,&one));
    dot = PetscRealPart(BLASdot_(&lds_,S+j*lds,&one,x,&one));
    PetscCall(BV_SafeSqrt(bv,dot,norm));
  } else {
    /* Compute *norm = BLASnrm2_(&lds_,S+j*lds,&one); */
    if (lds==1) *norm = PetscAbsScalar(S[j*lds]);
    else {
      for (i=0;i<lds;i++) {
        aval = PetscAbsScalar(S[i+j*lds]);
        if (aval!=0.0) {
          if (PetscUnlikely(scale<aval)) {
            alpha = 1.0 + alpha*PetscSqr(scale/aval);
            scale = aval;
          } else alpha += PetscSqr(aval/scale);
        }
      }
      *norm = scale*PetscSqrtReal(alpha);
    }
  }
  PetscCall(MatDenseRestoreArrayRead(ctx->S,&S));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVOrthogonalizeGS1_Tensor(BV bv,PetscInt k,Vec v,PetscBool *which,PetscScalar *h,PetscScalar *c,PetscReal *onorm,PetscReal *norm)
{
  BV_TENSOR         *ctx = (BV_TENSOR*)bv->data;
  PetscScalar       *pS,*cc,*x,dot,sonem=-1.0,sone=1.0,szero=0.0;
  PetscInt          i,lds = ctx->ld*ctx->d;
  PetscBLASInt      lds_,k_,one=1;
  const PetscScalar *omega;

  PetscFunctionBegin;
  PetscCheck(!v,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Orthogonalization against an external vector is not allowed in BVTENSOR");
  PetscCall(MatDenseGetArray(ctx->S,&pS));
  if (!c) PetscCall(VecGetArray(bv->buffer,&cc));
  else cc = c;
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(k,&k_));

  if (onorm) PetscCall(BVTensorNormColumn(bv,k,onorm));

  if (ctx->qB) x = ctx->sw;
  else x = pS+k*lds;

  if (PetscUnlikely(bv->orthog_type==BV_ORTHOG_MGS)) {  /* modified Gram-Schmidt */

    if (PetscUnlikely(bv->indef)) { /* signature */
      PetscCall(VecGetArrayRead(bv->omega,&omega));
    }
    for (i=-bv->nc;i<k;i++) {
      if (which && i>=0 && !which[i]) continue;
      if (ctx->qB) PetscCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&lds_,&sone,ctx->qB,&lds_,pS+k*lds,&one,&szero,x,&one));
      /* c_i = (s_k, s_i) */
      dot = PetscRealPart(BLASdot_(&lds_,pS+i*lds,&one,x,&one));
      if (bv->indef) dot /= PetscRealPart(omega[i]);
      PetscCall(BV_SetValue(bv,i,0,cc,dot));
      /* s_k = s_k - c_i s_i */
      dot = -dot;
      PetscCallBLAS("BLASaxpy",BLASaxpy_(&lds_,&dot,pS+i*lds,&one,pS+k*lds,&one));
    }
    if (PetscUnlikely(bv->indef)) PetscCall(VecRestoreArrayRead(bv->omega,&omega));

  } else {  /* classical Gram-Schmidt */
    if (ctx->qB) PetscCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&lds_,&sone,ctx->qB,&lds_,pS+k*lds,&one,&szero,x,&one));

    /* cc = S_{0:k-1}^* s_k */
    PetscCallBLAS("BLASgemv",BLASgemv_("C",&lds_,&k_,&sone,pS,&lds_,x,&one,&szero,cc,&one));

    /* s_k = s_k - S_{0:k-1} cc */
    if (PetscUnlikely(bv->indef)) PetscCall(BV_ApplySignature(bv,k,cc,PETSC_TRUE));
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&k_,&sonem,pS,&lds_,cc,&one,&sone,pS+k*lds,&one));
    if (PetscUnlikely(bv->indef)) PetscCall(BV_ApplySignature(bv,k,cc,PETSC_FALSE));
  }

  if (norm) PetscCall(BVTensorNormColumn(bv,k,norm));
  PetscCall(BV_AddCoefficients(bv,k,h,cc));
  PetscCall(MatDenseRestoreArray(ctx->S,&pS));
  PetscCall(VecRestoreArray(bv->buffer,&cc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVView_Tensor(BV bv,PetscViewer viewer)
{
  BV_TENSOR         *ctx = (BV_TENSOR*)bv->data;
  PetscViewerFormat format;
  PetscBool         isascii;
  const char        *bvname,*uname,*sname;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"number of tensor blocks (degree): %" PetscInt_FMT "\n",ctx->d));
      PetscCall(PetscViewerASCIIPrintf(viewer,"number of columns of U factor: %" PetscInt_FMT "\n",ctx->ld));
      PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscCall(BVView(ctx->U,viewer));
    PetscCall(MatView(ctx->S,viewer));
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      PetscCall(PetscObjectGetName((PetscObject)bv,&bvname));
      PetscCall(PetscObjectGetName((PetscObject)ctx->U,&uname));
      PetscCall(PetscObjectGetName((PetscObject)ctx->S,&sname));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%s=kron(eye(%" PetscInt_FMT "),%s)*%s(:,1:%" PetscInt_FMT ");\n",bvname,ctx->d,uname,sname,bv->k));
    }
  } else {
    PetscCall(BVView(ctx->U,viewer));
    PetscCall(MatView(ctx->S,viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVTensorUpdateMatrix(BV V,PetscInt ini,PetscInt end)
{
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscInt       i,j,r,c,l,k,ld=ctx->ld,lds=ctx->d*ctx->ld;
  PetscScalar    *qB,*sqB;
  Vec            u;
  Mat            A;

  PetscFunctionBegin;
  if (!V->matrix) PetscFunctionReturn(PETSC_SUCCESS);
  l = ctx->U->l; k = ctx->U->k;
  /* update inner product matrix */
  if (!ctx->qB) {
    PetscCall(PetscCalloc2(lds*lds,&ctx->qB,lds,&ctx->sw));
    PetscCall(BVCreateVec(ctx->U,&ctx->u));
  }
  ctx->U->l = 0;
  for (r=0;r<ctx->d;r++) {
    for (c=0;c<=r;c++) {
      PetscCall(MatNestGetSubMat(V->matrix,r,c,&A));
      if (A) {
        qB = ctx->qB+c*ld*lds+r*ld;
        for (i=ini;i<end;i++) {
          PetscCall(BVGetColumn(ctx->U,i,&u));
          PetscCall(MatMult(A,u,ctx->u));
          ctx->U->k = i+1;
          PetscCall(BVDotVec(ctx->U,ctx->u,qB+i*lds));
          PetscCall(BVRestoreColumn(ctx->U,i,&u));
          for (j=0;j<i;j++) qB[i+j*lds] = PetscConj(qB[j+i*lds]);
          qB[i*lds+i] = PetscRealPart(qB[i+i*lds]);
        }
        if (PetscUnlikely(c!=r)) {
          sqB = ctx->qB+r*ld*lds+c*ld;
          for (i=ini;i<end;i++) for (j=0;j<=i;j++) {
            sqB[i+j*lds] = PetscConj(qB[j+i*lds]);
            sqB[j+i*lds] = qB[j+i*lds];
          }
        }
      }
    }
  }
  ctx->U->l = l; ctx->U->k = k;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVTensorBuildFirstColumn_Tensor(BV V,PetscInt k)
{
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscInt       i,nq=0;
  PetscScalar    *pS,*omega;
  PetscReal      norm;
  PetscBool      breakdown=PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(ctx->S,&pS));
  for (i=0;i<ctx->d;i++) {
    if (i>=k) PetscCall(BVSetRandomColumn(ctx->U,nq));
    else PetscCall(BVCopyColumn(ctx->U,i,nq));
    PetscCall(BVOrthogonalizeColumn(ctx->U,nq,pS+i*ctx->ld,&norm,&breakdown));
    if (!breakdown) {
      PetscCall(BVScaleColumn(ctx->U,nq,1.0/norm));
      pS[nq+i*ctx->ld] = norm;
      nq++;
    }
  }
  PetscCall(MatDenseRestoreArray(ctx->S,&pS));
  PetscCheck(nq,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Cannot build first column of tensor BV; U should contain k=%" PetscInt_FMT " nonzero columns",k);
  PetscCall(BVTensorUpdateMatrix(V,0,nq));
  PetscCall(BVTensorNormColumn(V,0,&norm));
  PetscCall(BVScale_Tensor(V,0,1.0/norm));
  if (V->indef) {
    PetscCall(BV_AllocateSignature(V));
    PetscCall(VecGetArray(V->omega,&omega));
    omega[0] = (norm<0.0)? -1.0: 1.0;
    PetscCall(VecRestoreArray(V->omega,&omega));
  }
  /* set active columns */
  ctx->U->l = 0;
  ctx->U->k = nq;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVTensorBuildFirstColumn - Builds the first column of the tensor basis vectors
   V from the data contained in the first k columns of U.

   Collective

   Input Parameters:
+  V - the basis vectors context
-  k - the number of columns of U with relevant information

   Notes:
   At most d columns are considered, where d is the degree of the tensor BV.
   Given V = (I otimes U) S, this function computes the first column of V, that
   is, it computes the coefficients of the first column of S by orthogonalizing
   the first d columns of U. If k is less than d (or linearly dependent columns
   are found) then additional random columns are used.

   The computed column has unit norm.

   Level: advanced

.seealso: BVCreateTensor()
@*/
PetscErrorCode BVTensorBuildFirstColumn(BV V,PetscInt k)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(V,k,2);
  PetscUseMethod(V,"BVTensorBuildFirstColumn_C",(BV,PetscInt),(V,k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVTensorCompress_Tensor(BV V,PetscInt newc)
{
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscInt       nwu=0,nnc,nrow,lwa,r,c;
  PetscInt       i,j,k,n,lds=ctx->ld*ctx->d,deg=ctx->d,lock,cs1=V->k,rs1=ctx->U->k,rk=0,offu;
  PetscScalar    *S,*M,*Z,*pQ,*SS,*SS2,t,sone=1.0,zero=0.0,mone=-1.0,*p,*tau,*work,*qB,*sqB;
  PetscReal      *sg,tol,*rwork;
  PetscBLASInt   ld_,cs1_,rs1_,cs1tdeg,n_,info,lw_,newc_,newctdeg,nnc_,nrow_,nnctdeg,lds_,rk_;
  Mat            Q,A;

  PetscFunctionBegin;
  if (!cs1) PetscFunctionReturn(PETSC_SUCCESS);
  lwa = 6*ctx->ld*lds+2*cs1;
  n = PetscMin(rs1,deg*cs1);
  lock = ctx->U->l;
  nnc = cs1-lock-newc;
  nrow = rs1-lock;
  PetscCall(PetscCalloc6(deg*newc*nnc,&SS,newc*nnc,&SS2,(rs1+lock+newc)*n,&pQ,deg*rs1,&tau,lwa,&work,6*n,&rwork));
  offu = lock*(rs1+1);
  M = work+nwu;
  nwu += rs1*cs1*deg;
  sg = rwork;
  Z = work+nwu;
  nwu += deg*cs1*n;
  PetscCall(PetscBLASIntCast(n,&n_));
  PetscCall(PetscBLASIntCast(nnc,&nnc_));
  PetscCall(PetscBLASIntCast(cs1,&cs1_));
  PetscCall(PetscBLASIntCast(rs1,&rs1_));
  PetscCall(PetscBLASIntCast(newc,&newc_));
  PetscCall(PetscBLASIntCast(newc*deg,&newctdeg));
  PetscCall(PetscBLASIntCast(nnc*deg,&nnctdeg));
  PetscCall(PetscBLASIntCast(cs1*deg,&cs1tdeg));
  PetscCall(PetscBLASIntCast(lwa-nwu,&lw_));
  PetscCall(PetscBLASIntCast(nrow,&nrow_));
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(MatDenseGetArray(ctx->S,&S));

  if (newc>0) {
    /* truncate columns associated with new converged eigenpairs */
    for (j=0;j<deg;j++) {
      for (i=lock;i<lock+newc;i++) PetscCall(PetscArraycpy(M+(i-lock+j*newc)*nrow,S+i*lds+j*ctx->ld+lock,nrow));
    }
    PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined (PETSC_USE_COMPLEX)
    PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&nrow_,&newctdeg,M,&nrow_,sg,pQ+offu,&rs1_,Z,&n_,work+nwu,&lw_,&info));
#else
    PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&nrow_,&newctdeg,M,&nrow_,sg,pQ+offu,&rs1_,Z,&n_,work+nwu,&lw_,rwork+n,&info));
#endif
    SlepcCheckLapackInfo("gesvd",info);
    PetscCall(PetscFPTrapPop());
    /* SVD has rank min(newc,nrow) */
    rk = PetscMin(newc,nrow);
    for (i=0;i<rk;i++) {
      t = sg[i];
      PetscCallBLAS("BLASscal",BLASscal_(&newctdeg,&t,Z+i,&n_));
    }
    for (i=0;i<deg;i++) {
      for (j=lock;j<lock+newc;j++) {
        PetscCall(PetscArraycpy(S+j*lds+i*ctx->ld+lock,Z+(newc*i+j-lock)*n,rk));
        PetscCall(PetscArrayzero(S+j*lds+i*ctx->ld+lock+rk,(ctx->ld-lock-rk)));
      }
    }
    /*
      update columns associated with non-converged vectors, orthogonalize
      against pQ so that next M has rank nnc+d-1 instead of nrow+d-1
    */
    for (i=0;i<deg;i++) {
      PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&newc_,&nnc_,&nrow_,&sone,pQ+offu,&rs1_,S+(lock+newc)*lds+i*ctx->ld+lock,&lds_,&zero,PetscSafePointerPlusOffset(SS,i*nnc*newc),&newc_));
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&nrow_,&nnc_,&newc_,&mone,pQ+offu,&rs1_,PetscSafePointerPlusOffset(SS,i*nnc*newc),&newc_,&sone,S+(lock+newc)*lds+i*ctx->ld+lock,&lds_));
      /* repeat orthogonalization step */
      PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&newc_,&nnc_,&nrow_,&sone,pQ+offu,&rs1_,S+(lock+newc)*lds+i*ctx->ld+lock,&lds_,&zero,SS2,&newc_));
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&nrow_,&nnc_,&newc_,&mone,pQ+offu,&rs1_,SS2,&newc_,&sone,S+(lock+newc)*lds+i*ctx->ld+lock,&lds_));
      for (j=0;j<newc*nnc;j++) *(SS+i*newc*nnc+j) += SS2[j];
    }
  }

  /* truncate columns associated with non-converged eigenpairs */
  for (j=0;j<deg;j++) {
    for (i=lock+newc;i<cs1;i++) PetscCall(PetscArraycpy(M+(i-lock-newc+j*nnc)*nrow,S+i*lds+j*ctx->ld+lock,nrow));
  }
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined (PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&nrow_,&nnctdeg,M,&nrow_,sg,pQ+offu+newc*rs1,&rs1_,Z,&n_,work+nwu,&lw_,&info));
#else
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&nrow_,&nnctdeg,M,&nrow_,sg,pQ+offu+newc*rs1,&rs1_,Z,&n_,work+nwu,&lw_,rwork+n,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);
  PetscCall(PetscFPTrapPop());
  tol = PetscMax(rs1,deg*cs1)*PETSC_MACHINE_EPSILON*sg[0];
  rk = 0;
  for (i=0;i<PetscMin(nrow,nnctdeg);i++) if (sg[i]>tol) rk++;
  rk = PetscMin(nnc+deg-1,rk);
  /* the SVD has rank (at most) nnc+deg-1 */
  for (i=0;i<rk;i++) {
    t = sg[i];
    PetscCallBLAS("BLASscal",BLASscal_(&nnctdeg,&t,Z+i,&n_));
  }
  /* update S */
  PetscCall(PetscArrayzero(S+cs1*lds,(V->m-cs1)*lds));
  k = ctx->ld-lock-newc-rk;
  for (i=0;i<deg;i++) {
    for (j=lock+newc;j<cs1;j++) {
      PetscCall(PetscArraycpy(S+j*lds+i*ctx->ld+lock+newc,Z+(nnc*i+j-lock-newc)*n,rk));
      PetscCall(PetscArrayzero(S+j*lds+i*ctx->ld+lock+newc+rk,k));
    }
  }
  if (newc>0) {
    for (i=0;i<deg;i++) {
      p = PetscSafePointerPlusOffset(SS,i*nnc*newc);
      for (j=lock+newc;j<cs1;j++) {
        for (k=0;k<newc;k++) S[j*lds+i*ctx->ld+lock+k] = *(p++);
      }
    }
  }

  /* orthogonalize pQ */
  rk = rk+newc;
  PetscCall(PetscBLASIntCast(rk,&rk_));
  PetscCall(PetscBLASIntCast(cs1-lock,&nnc_));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&nrow_,&rk_,pQ+offu,&rs1_,tau,work+nwu,&lw_,&info));
  SlepcCheckLapackInfo("geqrf",info);
  for (i=0;i<deg;i++) {
    PetscCallBLAS("BLAStrmm",BLAStrmm_("L","U","N","N",&rk_,&nnc_,&sone,pQ+offu,&rs1_,S+lock*lds+lock+i*ctx->ld,&lds_));
  }
  PetscCallBLAS("LAPACKorgqr",LAPACKorgqr_(&nrow_,&rk_,&rk_,pQ+offu,&rs1_,tau,work+nwu,&lw_,&info));
  SlepcCheckLapackInfo("orgqr",info);
  PetscCall(PetscFPTrapPop());

  /* update vectors U(:,idx) = U*Q(:,idx) */
  rk = rk+lock;
  for (i=0;i<lock;i++) pQ[i*(1+rs1)] = 1.0;
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,rs1,rk,pQ,&Q));
  ctx->U->k = rs1;
  PetscCall(BVMultInPlace(ctx->U,Q,lock,rk));
  PetscCall(MatDestroy(&Q));

  if (ctx->qB) {
   /* update matrix qB */
    PetscCall(PetscBLASIntCast(ctx->ld,&ld_));
    PetscCall(PetscBLASIntCast(rk,&rk_));
    for (r=0;r<ctx->d;r++) {
      for (c=0;c<=r;c++) {
        PetscCall(MatNestGetSubMat(V->matrix,r,c,&A));
        if (A) {
          qB = ctx->qB+r*ctx->ld+c*ctx->ld*lds;
          PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&rs1_,&rk_,&rs1_,&sone,qB,&lds_,pQ,&rs1_,&zero,work+nwu,&rs1_));
          PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&rk_,&rk_,&rs1_,&sone,pQ,&rs1_,work+nwu,&rs1_,&zero,qB,&lds_));
          for (i=0;i<rk;i++) {
            for (j=0;j<i;j++) qB[i+j*lds] = PetscConj(qB[j+i*lds]);
            qB[i+i*lds] = PetscRealPart(qB[i+i*lds]);
          }
          for (i=rk;i<ctx->ld;i++) PetscCall(PetscArrayzero(qB+i*lds,ctx->ld));
          for (i=0;i<rk;i++) PetscCall(PetscArrayzero(qB+i*lds+rk,(ctx->ld-rk)));
          if (c!=r) {
            sqB = ctx->qB+r*ctx->ld*lds+c*ctx->ld;
            for (i=0;i<ctx->ld;i++) for (j=0;j<ctx->ld;j++) sqB[i+j*lds] = PetscConj(qB[j+i*lds]);
          }
        }
      }
    }
  }

  /* free work space */
  PetscCall(PetscFree6(SS,SS2,pQ,tau,work,rwork));
  PetscCall(MatDenseRestoreArray(ctx->S,&S));

  /* set active columns */
  if (newc) ctx->U->l += newc;
  ctx->U->k = rk;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVTensorCompress - Updates the U and S factors of the tensor basis vectors
   object V by means of an SVD, removing redundant information.

   Collective

   Input Parameters:
+  V - the tensor basis vectors context
-  newc - additional columns to be locked

   Notes:
   This function is typically used when restarting Krylov solvers. Truncating a
   tensor BV V = (I otimes U) S to its leading columns amounts to keeping the
   leading columns of S. However, to effectively reduce the size of the
   decomposition, it is necessary to compress it in a way that fewer columns of
   U are employed. This can be achieved by means of an update that involves the
   SVD of the low-rank matrix [S_0 S_1 ... S_{d-1}], where S_i are the pieces of S.

   If newc is nonzero, then newc columns are added to the leading columns of V.
   This means that the corresponding columns of the U and S factors will remain
   invariant in subsequent operations.

   Level: advanced

.seealso: BVCreateTensor(), BVSetActiveColumns()
@*/
PetscErrorCode BVTensorCompress(BV V,PetscInt newc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(V,newc,2);
  PetscUseMethod(V,"BVTensorCompress_C",(BV,PetscInt),(V,newc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVTensorGetDegree_Tensor(BV bv,PetscInt *d)
{
  BV_TENSOR *ctx = (BV_TENSOR*)bv->data;

  PetscFunctionBegin;
  *d = ctx->d;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVTensorGetDegree - Returns the number of blocks (degree) of the tensor BV.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  d - the degree

   Level: advanced

.seealso: BVCreateTensor()
@*/
PetscErrorCode BVTensorGetDegree(BV bv,PetscInt *d)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(d,2);
  PetscUseMethod(bv,"BVTensorGetDegree_C",(BV,PetscInt*),(bv,d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVTensorGetFactors_Tensor(BV V,BV *U,Mat *S)
{
  BV_TENSOR *ctx = (BV_TENSOR*)V->data;

  PetscFunctionBegin;
  PetscCheck(ctx->puk==-1,PetscObjectComm((PetscObject)V),PETSC_ERR_ORDER,"Previous call to BVTensonGetFactors without a BVTensorRestoreFactors call");
  ctx->puk = ctx->U->k;
  if (U) *U = ctx->U;
  if (S) *S = ctx->S;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVTensorGetFactors - Returns the two factors involved in the definition of the
   tensor basis vectors object, V = (I otimes U) S.

   Logically Collective

   Input Parameter:
.  V - the basis vectors context

   Output Parameters:
+  U - the BV factor
-  S - the Mat factor

   Notes:
   The returned factors are references (not copies) of the internal factors,
   so modifying them will change the tensor BV as well. Some operations of the
   tensor BV assume that U has orthonormal columns, so if the user modifies U
   this restriction must be taken into account.

   The returned factors must not be destroyed. BVTensorRestoreFactors() must
   be called when they are no longer needed.

   Pass a NULL vector for any of the arguments that is not needed.

   Level: advanced

.seealso: BVTensorRestoreFactors()
@*/
PetscErrorCode BVTensorGetFactors(BV V,BV *U,Mat *S)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscUseMethod(V,"BVTensorGetFactors_C",(BV,BV*,Mat*),(V,U,S));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVTensorRestoreFactors_Tensor(BV V,BV *U,Mat *S)
{
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  if (U) *U = NULL;
  if (S) *S = NULL;
  PetscCall(BVTensorUpdateMatrix(V,ctx->puk,ctx->U->k));
  ctx->puk = -1;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVTensorRestoreFactors - Restore the two factors that were obtained with
   BVTensorGetFactors().

   Logically Collective

   Input Parameters:
+  V - the basis vectors context
.  U - the BV factor (or NULL)
-  S - the Mat factor (or NULL)

   Notes:
   The arguments must match the corresponding call to BVTensorGetFactors().

   Level: advanced

.seealso: BVTensorGetFactors()
@*/
PetscErrorCode BVTensorRestoreFactors(BV V,BV *U,Mat *S)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  if (U) PetscValidHeaderSpecific(*U,BV_CLASSID,2);
  if (S) PetscValidHeaderSpecific(*S,MAT_CLASSID,3);
  PetscUseMethod(V,"BVTensorRestoreFactors_C",(BV,BV*,Mat*),(V,U,S));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVDestroy_Tensor(BV bv)
{
  BV_TENSOR      *ctx = (BV_TENSOR*)bv->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&ctx->U));
  PetscCall(MatDestroy(&ctx->S));
  if (ctx->u) {
    PetscCall(PetscFree2(ctx->qB,ctx->sw));
    PetscCall(VecDestroy(&ctx->u));
  }
  PetscCall(PetscFree(bv->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorBuildFirstColumn_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorCompress_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorGetDegree_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorGetFactors_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorRestoreFactors_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode BVCreate_Tensor(BV bv)
{
  BV_TENSOR      *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  bv->data = (void*)ctx;
  ctx->puk = -1;

  bv->ops->multinplace      = BVMultInPlace_Tensor;
  bv->ops->multinplacetrans = BVMultInPlaceHermitianTranspose_Tensor;
  bv->ops->dot              = BVDot_Tensor;
  bv->ops->scale            = BVScale_Tensor;
  bv->ops->norm             = BVNorm_Tensor;
  bv->ops->copycolumn       = BVCopyColumn_Tensor;
  bv->ops->gramschmidt      = BVOrthogonalizeGS1_Tensor;
  bv->ops->destroy          = BVDestroy_Tensor;
  bv->ops->view             = BVView_Tensor;

  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorBuildFirstColumn_C",BVTensorBuildFirstColumn_Tensor));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorCompress_C",BVTensorCompress_Tensor));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorGetDegree_C",BVTensorGetDegree_Tensor));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorGetFactors_C",BVTensorGetFactors_Tensor));
  PetscCall(PetscObjectComposeFunction((PetscObject)bv,"BVTensorRestoreFactors_C",BVTensorRestoreFactors_Tensor));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCreateTensor - Creates a tensor BV that is represented in compact form
   as V = (I otimes U) S, where U has orthonormal columns.

   Collective

   Input Parameters:
+  U - a basis vectors object
-  d - the number of blocks (degree) of the tensor BV

   Output Parameter:
.  V - the new basis vectors context

   Notes:
   The new basis vectors object is V = (I otimes U) S, where otimes denotes
   the Kronecker product, I is the identity matrix of order d, and S is a
   sequential matrix allocated internally. This compact representation is
   used e.g. to represent the Krylov basis generated with the linearization
   of a matrix polynomial of degree d.

   The size of V (number of rows) is equal to d times n, where n is the size
   of U. The dimensions of S are d times m rows and m-d+1 columns, where m is
   the number of columns of U, so m should be at least d.

   The communicator of V will be the same as U.

   On input, the content of U is irrelevant. Alternatively, it may contain
   some nonzero columns that will be used by BVTensorBuildFirstColumn().

   Level: advanced

.seealso: BVTensorGetDegree(), BVTensorGetFactors(), BVTensorBuildFirstColumn()
@*/
PetscErrorCode BVCreateTensor(BV U,PetscInt d,BV *V)
{
  PetscBool      match;
  PetscInt       n,N,m;
  VecType        vtype;
  BV_TENSOR      *ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(U,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(U,d,2);
  PetscCall(PetscObjectTypeCompare((PetscObject)U,BVTENSOR,&match));
  PetscCheck(!match,PetscObjectComm((PetscObject)U),PETSC_ERR_SUP,"U cannot be of type tensor");

  PetscCall(BVCreate(PetscObjectComm((PetscObject)U),V));
  PetscCall(BVGetSizes(U,&n,&N,&m));
  PetscCheck(m>=d,PetscObjectComm((PetscObject)U),PETSC_ERR_ARG_SIZ,"U has %" PetscInt_FMT " columns, it should have at least d=%" PetscInt_FMT,m,d);
  PetscCall(BVSetSizes(*V,d*n,d*N,m-d+1));
  PetscCall(BVGetVecType(U,&vtype));
  PetscCall(BVSetVecType(*V,vtype));
  PetscCall(PetscObjectChangeTypeName((PetscObject)*V,BVTENSOR));
  PetscCall(PetscLogEventBegin(BV_Create,*V,0,0,0));
  PetscCall(BVCreate_Tensor(*V));
  PetscCall(PetscLogEventEnd(BV_Create,*V,0,0,0));

  ctx = (BV_TENSOR*)(*V)->data;
  ctx->U  = U;
  ctx->d  = d;
  ctx->ld = m;
  PetscCall(PetscObjectReference((PetscObject)U));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,d*m,m-d+1,NULL,&ctx->S));
  PetscCall(PetscObjectSetName((PetscObject)ctx->S,"S"));

  /* Copy user-provided attributes of U */
  (*V)->orthog_type  = U->orthog_type;
  (*V)->orthog_ref   = U->orthog_ref;
  (*V)->orthog_eta   = U->orthog_eta;
  (*V)->orthog_block = U->orthog_block;
  (*V)->vmm          = U->vmm;
  (*V)->rrandom      = U->rrandom;
  PetscFunctionReturn(PETSC_SUCCESS);
}
