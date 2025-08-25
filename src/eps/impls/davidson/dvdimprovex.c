/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "davidson"

   Step: improve the eigenvectors X
*/

#include "davidson.h"
#include <slepcblaslapack.h>

/**** JD update step (I - Kfg'/(g'Kf)) K(A - sB) (I - Kfg'/(g'Kf)) t = (I - Kfg'/(g'Kf))r  ****/

typedef struct {
  PetscInt     size_X;
  KSP          ksp;                /* correction equation solver */
  Vec          friends;            /* reference vector for composite vectors */
  PetscScalar  theta[4],thetai[2]; /* the shifts used in the correction eq. */
  PetscInt     maxits;             /* maximum number of iterations */
  PetscInt     r_s,r_e;            /* the selected eigenpairs to improve */
  PetscInt     ksp_max_size;       /* the ksp maximum subvectors size */
  PetscReal    tol;                /* the maximum solution tolerance */
  PetscReal    lastTol;            /* last tol for dynamic stopping criterion */
  PetscReal    fix;                /* tolerance for using the approx. eigenvalue */
  PetscBool    dynamic;            /* if the dynamic stopping criterion is applied */
  dvdDashboard *d;                 /* the current dvdDashboard reference */
  PC           old_pc;             /* old pc in ksp */
  BV           KZ;                 /* KZ vecs for the projector KZ*inv(X'*KZ)*X' */
  BV           U;                  /* new X vectors */
  PetscScalar  *XKZ;               /* X'*KZ */
  PetscScalar  *iXKZ;              /* inverse of XKZ */
  PetscInt     ldXKZ;              /* leading dimension of XKZ */
  PetscInt     size_iXKZ;          /* size of iXKZ */
  PetscInt     ldiXKZ;             /* leading dimension of iXKZ */
  PetscInt     size_cX;            /* last value of d->size_cX */
  PetscInt     old_size_X;         /* last number of improved vectors */
  PetscBLASInt *iXKZPivots;        /* array of pivots */
} dvdImprovex_jd;

/*
   Compute (I - KZ*iXKZ*X')*V where,
   V, the vectors to apply the projector,
   cV, the number of vectors in V,
*/
static PetscErrorCode dvd_improvex_apply_proj(dvdDashboard *d,Vec *V,PetscInt cV)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscInt       i,ldh,k,l;
  PetscScalar    *h;
  PetscBLASInt   cV_,n,info,ld;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
#endif

  PetscFunctionBegin;
  PetscAssert(cV<=2,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");

  /* h <- X'*V */
  PetscCall(PetscMalloc1(data->size_iXKZ*cV,&h));
  ldh = data->size_iXKZ;
  PetscCall(BVGetActiveColumns(data->U,&l,&k));
  PetscAssert(ldh==k,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");
  PetscCall(BVSetActiveColumns(data->U,0,k));
  for (i=0;i<cV;i++) {
    PetscCall(BVDotVec(data->U,V[i],&h[ldh*i]));
#if defined(PETSC_USE_COMPLEX)
    for (j=0; j<k; j++) h[ldh*i+j] = PetscConj(h[ldh*i+j]);
#endif
  }
  PetscCall(BVSetActiveColumns(data->U,l,k));

  /* h <- iXKZ\h */
  PetscCall(PetscBLASIntCast(cV,&cV_));
  PetscCall(PetscBLASIntCast(data->size_iXKZ,&n));
  PetscCall(PetscBLASIntCast(data->ldiXKZ,&ld));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&n,&cV_,data->iXKZ,&ld,data->iXKZPivots,h,&n,&info));
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("getrs",info);

  /* V <- V - KZ*h */
  PetscCall(BVSetActiveColumns(data->KZ,0,k));
  for (i=0;i<cV;i++) PetscCall(BVMultVec(data->KZ,-1.0,1.0,V[i],&h[ldh*i]));
  PetscCall(BVSetActiveColumns(data->KZ,l,k));
  PetscCall(PetscFree(h));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Compute (I - X*iXKZ*KZ')*V where,
   V, the vectors to apply the projector,
   cV, the number of vectors in V,
*/
static PetscErrorCode dvd_improvex_applytrans_proj(dvdDashboard *d,Vec *V,PetscInt cV)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscInt       i,ldh,k,l;
  PetscScalar    *h;
  PetscBLASInt   cV_, n, info, ld;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
#endif

  PetscFunctionBegin;
  PetscAssert(cV<=2,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");

  /* h <- KZ'*V */
  PetscCall(PetscMalloc1(data->size_iXKZ*cV,&h));
  ldh = data->size_iXKZ;
  PetscCall(BVGetActiveColumns(data->U,&l,&k));
  PetscAssert(ldh==k,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");
  PetscCall(BVSetActiveColumns(data->KZ,0,k));
  for (i=0;i<cV;i++) {
    PetscCall(BVDotVec(data->KZ,V[i],&h[ldh*i]));
#if defined(PETSC_USE_COMPLEX)
    for (j=0;j<k;j++) h[ldh*i+j] = PetscConj(h[ldh*i+j]);
#endif
  }
  PetscCall(BVSetActiveColumns(data->KZ,l,k));

  /* h <- iXKZ\h */
  PetscCall(PetscBLASIntCast(cV,&cV_));
  PetscCall(PetscBLASIntCast(data->size_iXKZ,&n));
  PetscCall(PetscBLASIntCast(data->ldiXKZ,&ld));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKgetrs",LAPACKgetrs_("C",&n,&cV_,data->iXKZ,&ld,data->iXKZPivots,h,&n,&info));
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("getrs",info);

  /* V <- V - U*h */
  PetscCall(BVSetActiveColumns(data->U,0,k));
  for (i=0;i<cV;i++) PetscCall(BVMultVec(data->U,-1.0,1.0,V[i],&h[ldh*i]));
  PetscCall(BVSetActiveColumns(data->U,l,k));
  PetscCall(PetscFree(h));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_improvex_jd_end(dvdDashboard *d)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;

  PetscFunctionBegin;
  PetscCall(VecDestroy(&data->friends));

  /* Restore the pc of ksp */
  if (data->old_pc) {
    PetscCall(KSPSetPC(data->ksp, data->old_pc));
    PetscCall(PCDestroy(&data->old_pc));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_improvex_jd_d(dvdDashboard *d)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;

  PetscFunctionBegin;
  /* Free local data and objects */
  PetscCall(PetscFree(data->XKZ));
  PetscCall(PetscFree(data->iXKZ));
  PetscCall(PetscFree(data->iXKZPivots));
  PetscCall(BVDestroy(&data->KZ));
  PetscCall(BVDestroy(&data->U));
  PetscCall(PetscFree(data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   y <- theta[1]A*x - theta[0]*B*x
   auxV, two auxiliary vectors
 */
static inline PetscErrorCode dvd_aux_matmult(dvdImprovex_jd *data,const Vec *x,const Vec *y)
{
  PetscInt       n,i;
  const Vec      *Bx;
  Vec            *auxV;

  PetscFunctionBegin;
  n = data->r_e - data->r_s;
  for (i=0;i<n;i++) PetscCall(MatMult(data->d->A,x[i],y[i]));

  PetscCall(SlepcVecPoolGetVecs(data->d->auxV,2,&auxV));
  for (i=0;i<n;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (PetscUnlikely(data->d->eigi[data->r_s+i] != 0.0)) {
      if (data->d->B) {
        PetscCall(MatMult(data->d->B,x[i],auxV[0]));
        PetscCall(MatMult(data->d->B,x[i+1],auxV[1]));
        Bx = auxV;
      } else Bx = &x[i];

      /* y_i   <- [ t_2i+1*A*x_i   - t_2i*Bx_i + ti_i*Bx_i+1;
         y_i+1      t_2i+1*A*x_i+1 - ti_i*Bx_i - t_2i*Bx_i+1  ] */
      PetscCall(VecAXPBYPCZ(y[i],-data->theta[2*i],data->thetai[i],data->theta[2*i+1],Bx[0],Bx[1]));
      PetscCall(VecAXPBYPCZ(y[i+1],-data->thetai[i],-data->theta[2*i],data->theta[2*i+1],Bx[0],Bx[1]));
      i++;
    } else
#endif
    {
      if (data->d->B) {
        PetscCall(MatMult(data->d->B,x[i],auxV[0]));
        Bx = auxV;
      } else Bx = &x[i];
      PetscCall(VecAXPBY(y[i],-data->theta[i*2],data->theta[i*2+1],Bx[0]));
    }
  }
  PetscCall(SlepcVecPoolRestoreVecs(data->d->auxV,2,&auxV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   y <- theta[1]'*A'*x - theta[0]'*B'*x
 */
static inline PetscErrorCode dvd_aux_matmulttrans(dvdImprovex_jd *data,const Vec *x,const Vec *y)
{
  PetscInt       n,i;
  const Vec      *Bx;
  Vec            *auxV;

  PetscFunctionBegin;
  n = data->r_e - data->r_s;
  for (i=0;i<n;i++) PetscCall(MatMultTranspose(data->d->A,x[i],y[i]));

  PetscCall(SlepcVecPoolGetVecs(data->d->auxV,2,&auxV));
  for (i=0;i<n;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (data->d->eigi[data->r_s+i] != 0.0) {
      if (data->d->B) {
        PetscCall(MatMultTranspose(data->d->B,x[i],auxV[0]));
        PetscCall(MatMultTranspose(data->d->B,x[i+1],auxV[1]));
        Bx = auxV;
      } else Bx = &x[i];

      /* y_i   <- [ t_2i+1*A*x_i   - t_2i*Bx_i - ti_i*Bx_i+1;
         y_i+1      t_2i+1*A*x_i+1 + ti_i*Bx_i - t_2i*Bx_i+1  ] */
      PetscCall(VecAXPBYPCZ(y[i],-data->theta[2*i],-data->thetai[i],data->theta[2*i+1],Bx[0],Bx[1]));
      PetscCall(VecAXPBYPCZ(y[i+1],data->thetai[i],-data->theta[2*i],data->theta[2*i+1],Bx[0],Bx[1]));
      i++;
    } else
#endif
    {
      if (data->d->B) {
        PetscCall(MatMultTranspose(data->d->B,x[i],auxV[0]));
        Bx = auxV;
      } else Bx = &x[i];
      PetscCall(VecAXPBY(y[i],PetscConj(-data->theta[i*2]),PetscConj(data->theta[i*2+1]),Bx[0]));
    }
  }
  PetscCall(SlepcVecPoolRestoreVecs(data->d->auxV,2,&auxV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCApplyBA_dvd(PC pc,PCSide side,Vec in,Vec out,Vec w)
{
  dvdImprovex_jd *data;
  PetscInt       n,i;
  const Vec      *inx,*outx,*wx;
  Vec            *auxV;
  Mat            A;

  PetscFunctionBegin;
  PetscCall(PCGetOperators(pc,&A,NULL));
  PetscCall(MatShellGetContext(A,&data));
  PetscCall(VecCompGetSubVecs(in,NULL,&inx));
  PetscCall(VecCompGetSubVecs(out,NULL,&outx));
  PetscCall(VecCompGetSubVecs(w,NULL,&wx));
  n = data->r_e - data->r_s;
  PetscCall(SlepcVecPoolGetVecs(data->d->auxV,n,&auxV));
  switch (side) {
  case PC_LEFT:
    /* aux <- theta[1]A*in - theta[0]*B*in */
    PetscCall(dvd_aux_matmult(data,inx,auxV));

    /* out <- K * aux */
    for (i=0;i<n;i++) PetscCall(data->d->improvex_precond(data->d,data->r_s+i,auxV[i],outx[i]));
    break;
  case PC_RIGHT:
    /* aux <- K * in */
    for (i=0;i<n;i++) PetscCall(data->d->improvex_precond(data->d,data->r_s+i,inx[i],auxV[i]));

    /* out <- theta[1]A*auxV - theta[0]*B*auxV */
    PetscCall(dvd_aux_matmult(data,auxV,outx));
    break;
  case PC_SYMMETRIC:
    /* aux <- K^{1/2} * in */
    for (i=0;i<n;i++) PetscCall(PCApplySymmetricRight(data->old_pc,inx[i],auxV[i]));

    /* wx <- theta[1]A*auxV - theta[0]*B*auxV */
    PetscCall(dvd_aux_matmult(data,auxV,wx));

    /* aux <- K^{1/2} * in */
    for (i=0;i<n;i++) PetscCall(PCApplySymmetricLeft(data->old_pc,wx[i],outx[i]));
    break;
  default:
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unsupported KSP side");
  }
  /* out <- out - v*(u'*out) */
  PetscCall(dvd_improvex_apply_proj(data->d,(Vec*)outx,n));
  PetscCall(SlepcVecPoolRestoreVecs(data->d->auxV,n,&auxV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCApply_dvd(PC pc,Vec in,Vec out)
{
  dvdImprovex_jd *data;
  PetscInt       n,i;
  const Vec      *inx, *outx;
  Mat            A;

  PetscFunctionBegin;
  PetscCall(PCGetOperators(pc,&A,NULL));
  PetscCall(MatShellGetContext(A,&data));
  PetscCall(VecCompGetSubVecs(in,NULL,&inx));
  PetscCall(VecCompGetSubVecs(out,NULL,&outx));
  n = data->r_e - data->r_s;
  /* out <- K * in */
  for (i=0;i<n;i++) PetscCall(data->d->improvex_precond(data->d,data->r_s+i,inx[i],outx[i]));
  /* out <- out - v*(u'*out) */
  PetscCall(dvd_improvex_apply_proj(data->d,(Vec*)outx,n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCApplyTranspose_dvd(PC pc,Vec in,Vec out)
{
  dvdImprovex_jd *data;
  PetscInt       n,i;
  const Vec      *inx, *outx;
  Vec            *auxV;
  Mat            A;

  PetscFunctionBegin;
  PetscCall(PCGetOperators(pc,&A,NULL));
  PetscCall(MatShellGetContext(A,&data));
  PetscCall(VecCompGetSubVecs(in,NULL,&inx));
  PetscCall(VecCompGetSubVecs(out,NULL,&outx));
  n = data->r_e - data->r_s;
  PetscCall(SlepcVecPoolGetVecs(data->d->auxV,n,&auxV));
  /* auxV <- in */
  for (i=0;i<n;i++) PetscCall(VecCopy(inx[i],auxV[i]));
  /* auxV <- auxV - u*(v'*auxV) */
  PetscCall(dvd_improvex_applytrans_proj(data->d,auxV,n));
  /* out <- K' * aux */
  for (i=0;i<n;i++) PetscCall(PCApplyTranspose(data->old_pc,auxV[i],outx[i]));
  PetscCall(SlepcVecPoolRestoreVecs(data->d->auxV,n,&auxV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_dvd_jd(Mat A,Vec in,Vec out)
{
  dvdImprovex_jd *data;
  PetscInt       n;
  const Vec      *inx, *outx;
  PCSide         side;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&data));
  PetscCall(VecCompGetSubVecs(in,NULL,&inx));
  PetscCall(VecCompGetSubVecs(out,NULL,&outx));
  n = data->r_e - data->r_s;
  /* out <- theta[1]A*in - theta[0]*B*in */
  PetscCall(dvd_aux_matmult(data,inx,outx));
  PetscCall(KSPGetPCSide(data->ksp,&side));
  if (side == PC_RIGHT) {
    /* out <- out - v*(u'*out) */
    PetscCall(dvd_improvex_apply_proj(data->d,(Vec*)outx,n));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultTranspose_dvd_jd(Mat A,Vec in,Vec out)
{
  dvdImprovex_jd *data;
  PetscInt       n,i;
  const Vec      *inx,*outx,*r;
  Vec            *auxV;
  PCSide         side;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&data));
  PetscCall(VecCompGetSubVecs(in,NULL,&inx));
  PetscCall(VecCompGetSubVecs(out,NULL,&outx));
  n = data->r_e - data->r_s;
  PetscCall(KSPGetPCSide(data->ksp,&side));
  if (side == PC_RIGHT) {
    /* auxV <- in */
    PetscCall(SlepcVecPoolGetVecs(data->d->auxV,n,&auxV));
    for (i=0;i<n;i++) PetscCall(VecCopy(inx[i],auxV[i]));
    /* auxV <- auxV - v*(u'*auxV) */
    PetscCall(dvd_improvex_applytrans_proj(data->d,auxV,n));
    r = auxV;
  } else r = inx;
  /* out <- theta[1]A*r - theta[0]*B*r */
  PetscCall(dvd_aux_matmulttrans(data,r,outx));
  if (side == PC_RIGHT) PetscCall(SlepcVecPoolRestoreVecs(data->d->auxV,n,&auxV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatCreateVecs_dvd_jd(Mat A,Vec *right,Vec *left)
{
  Vec            *r,*l;
  dvdImprovex_jd *data;
  PetscInt       n,i;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&data));
  n = data->ksp_max_size;
  if (right) PetscCall(PetscMalloc1(n,&r));
  if (left) PetscCall(PetscMalloc1(n,&l));
  for (i=0;i<n;i++) PetscCall(MatCreateVecs(data->d->A,right?&r[i]:NULL,left?&l[i]:NULL));
  if (right) {
    PetscCall(VecCreateCompWithVecs(r,n,data->friends,right));
    for (i=0;i<n;i++) PetscCall(VecDestroy(&r[i]));
  }
  if (left) {
    PetscCall(VecCreateCompWithVecs(l,n,data->friends,left));
    for (i=0;i<n;i++) PetscCall(VecDestroy(&l[i]));
  }

  if (right) PetscCall(PetscFree(r));
  if (left) PetscCall(PetscFree(l));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_improvex_jd_start(dvdDashboard *d)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscInt       rA, cA, rlA, clA;
  Mat            A;
  PetscBool      t;
  PC             pc;
  Vec            v0[2];

  PetscFunctionBegin;
  data->size_cX = data->old_size_X = 0;
  data->lastTol = data->dynamic?0.5:0.0;

  /* Setup the ksp */
  if (data->ksp) {
    /* Create the reference vector */
    PetscCall(BVGetColumn(d->eps->V,0,&v0[0]));
    v0[1] = v0[0];
    PetscCall(VecCreateCompWithVecs(v0,data->ksp_max_size,NULL,&data->friends));
    PetscCall(BVRestoreColumn(d->eps->V,0,&v0[0]));

    /* Save the current pc and set a PCNONE */
    PetscCall(KSPGetPC(data->ksp, &data->old_pc));
    PetscCall(PetscObjectTypeCompare((PetscObject)data->old_pc,PCNONE,&t));
    data->lastTol = 0.5;
    if (t) data->old_pc = NULL;
    else {
      PetscCall(PetscObjectReference((PetscObject)data->old_pc));
      PetscCall(PCCreate(PetscObjectComm((PetscObject)d->eps),&pc));
      PetscCall(PCSetType(pc,PCSHELL));
      PetscCall(PCSetOperators(pc,d->A,d->A));
      PetscCall(PCSetReusePreconditioner(pc,PETSC_TRUE));
      PetscCall(PCShellSetApply(pc,PCApply_dvd));
      PetscCall(PCShellSetApplyBA(pc,PCApplyBA_dvd));
      PetscCall(PCShellSetApplyTranspose(pc,PCApplyTranspose_dvd));
      PetscCall(KSPSetPC(data->ksp,pc));
      PetscCall(PCDestroy(&pc));
    }

    /* Create the (I-v*u')*K*(A-s*B) matrix */
    PetscCall(MatGetSize(d->A,&rA,&cA));
    PetscCall(MatGetLocalSize(d->A,&rlA,&clA));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)d->A),rlA*data->ksp_max_size,clA*data->ksp_max_size,rA*data->ksp_max_size,cA*data->ksp_max_size,data,&A));
    PetscCall(MatShellSetOperation(A,MATOP_MULT,(PetscErrorCodeFn*)MatMult_dvd_jd));
    PetscCall(MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_dvd_jd));
    PetscCall(MatShellSetOperation(A,MATOP_CREATE_VECS,(PetscErrorCodeFn*)MatCreateVecs_dvd_jd));

    /* Try to avoid KSPReset */
    PetscCall(KSPGetOperatorsSet(data->ksp,&t,NULL));
    if (t) {
      Mat      M;
      PetscInt rM;
      PetscCall(KSPGetOperators(data->ksp,&M,NULL));
      PetscCall(MatGetSize(M,&rM,NULL));
      if (rM != rA*data->ksp_max_size) PetscCall(KSPReset(data->ksp));
    }
    PetscCall(EPS_KSPSetOperators(data->ksp,A,A));
    PetscCall(KSPSetReusePreconditioner(data->ksp,PETSC_TRUE));
    PetscCall(KSPSetUp(data->ksp));
    PetscCall(MatDestroy(&A));
  } else {
    data->old_pc = NULL;
    data->friends = NULL;
  }
  PetscCall(BVSetActiveColumns(data->KZ,0,0));
  PetscCall(BVSetActiveColumns(data->U,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute: u <- X, v <- K*(theta[0]*A+theta[1]*B)*X,
  kr <- K^{-1}*(A-eig*B)*X, being X <- V*pX[i_s..i_e-1], Y <- W*pY[i_s..i_e-1]
  where
  pX,pY, the right and left eigenvectors of the projected system
  ld, the leading dimension of pX and pY
*/
static PetscErrorCode dvd_improvex_jd_proj_cuv(dvdDashboard *d,PetscInt i_s,PetscInt i_e,Vec *kr,PetscScalar *theta,PetscScalar *thetai,PetscScalar *pX,PetscScalar *pY,PetscInt ld)
{
  PetscInt          n=i_e-i_s,size_KZ,V_new,rm,i,lv,kv,lKZ,kKZ;
  dvdImprovex_jd    *data = (dvdImprovex_jd*)d->improveX_data;
  const PetscScalar *array;
  Mat               M;
  Vec               u[2],v[2];
  PetscBLASInt      s,ldXKZ,info;

  PetscFunctionBegin;
  /* Check consistency */
  PetscCall(BVGetActiveColumns(d->eps->V,&lv,&kv));
  V_new = lv - data->size_cX;
  PetscAssert(V_new<=data->old_size_X,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");
  data->old_size_X = n;
  data->size_cX = lv;

  /* KZ <- KZ(rm:rm+max_cX-1) */
  PetscCall(BVGetActiveColumns(data->KZ,&lKZ,&kKZ));
  rm = PetscMax(V_new+lKZ,0);
  if (rm > 0) {
    for (i=0;i<lKZ;i++) {
      PetscCall(BVCopyColumn(data->KZ,i+rm,i));
      PetscCall(BVCopyColumn(data->U,i+rm,i));
    }
  }

  /* XKZ <- XKZ(rm:rm+max_cX-1,rm:rm+max_cX-1) */
  if (rm > 0) {
    for (i=0;i<lKZ;i++) PetscCall(PetscArraycpy(&data->XKZ[i*data->ldXKZ+i],&data->XKZ[(i+rm)*data->ldXKZ+i+rm],lKZ));
  }
  lKZ = PetscMin(0,lKZ+V_new);
  PetscCall(BVSetActiveColumns(data->KZ,lKZ,lKZ+n));
  PetscCall(BVSetActiveColumns(data->U,lKZ,lKZ+n));

  /* Compute X, KZ and KR */
  PetscCall(BVGetColumn(data->U,lKZ,u));
  if (n>1) PetscCall(BVGetColumn(data->U,lKZ+1,&u[1]));
  PetscCall(BVGetColumn(data->KZ,lKZ,v));
  if (n>1) PetscCall(BVGetColumn(data->KZ,lKZ+1,&v[1]));
  PetscCall(d->improvex_jd_proj_uv(d,i_s,i_e,u,v,kr,theta,thetai,pX,pY,ld));
  PetscCall(BVRestoreColumn(data->U,lKZ,u));
  if (n>1) PetscCall(BVRestoreColumn(data->U,lKZ+1,&u[1]));
  PetscCall(BVRestoreColumn(data->KZ,lKZ,v));
  if (n>1) PetscCall(BVRestoreColumn(data->KZ,lKZ+1,&v[1]));

  /* XKZ <- U'*KZ */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,lKZ+n,lKZ+n,NULL,&M));
  PetscCall(BVMatProject(data->KZ,NULL,data->U,M));
  PetscCall(MatDenseGetArrayRead(M,&array));
  for (i=lKZ;i<lKZ+n;i++) { /* upper part */
    PetscCall(PetscArraycpy(&data->XKZ[data->ldXKZ*i],&array[i*(lKZ+n)],lKZ));
  }
  for (i=0;i<lKZ+n;i++) { /* lower part */
    PetscCall(PetscArraycpy(&data->XKZ[data->ldXKZ*i+lKZ],&array[i*(lKZ+n)+lKZ],n));
  }
  PetscCall(MatDenseRestoreArrayRead(M,&array));
  PetscCall(MatDestroy(&M));

  /* iXKZ <- inv(XKZ) */
  size_KZ = lKZ+n;
  PetscCall(PetscBLASIntCast(lKZ+n,&s));
  data->ldiXKZ = data->size_iXKZ = size_KZ;
  for (i=0;i<size_KZ;i++) PetscCall(PetscArraycpy(&data->iXKZ[data->ldiXKZ*i],&data->XKZ[data->ldXKZ*i],size_KZ));
  PetscCall(PetscBLASIntCast(data->ldiXKZ,&ldXKZ));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKgetrf",LAPACKgetrf_(&s,&s,data->iXKZ,&ldXKZ,data->iXKZPivots,&info));
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("getrf",info);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_improvex_jd_gen(dvdDashboard *d,PetscInt r_s,PetscInt r_e,PetscInt *size_D)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscInt       i,j,n,maxits,maxits0,lits,s,ld,k,max_size_D,lV,kV;
  PetscScalar    *pX,*pY;
  PetscReal      tol,tol0;
  Vec            *kr,kr_comp,D_comp,D[2],kr0[2];
  PetscBool      odd_situation = PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&lV,&kV));
  max_size_D = d->eps->ncv-kV;
  /* Quick exit */
  if ((max_size_D == 0) || r_e-r_s <= 0) {
   *size_D = 0;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  n = PetscMin(PetscMin(data->size_X, max_size_D), r_e-r_s);
  PetscAssert(n>0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"n == 0");
  PetscAssert(data->size_X>=r_e-r_s,PETSC_COMM_SELF,PETSC_ERR_PLIB,"size_X < r_e-r_s");

  PetscCall(DSGetLeadingDimension(d->eps->ds,&ld));

  /* Restart lastTol if a new pair converged */
  if (data->dynamic && data->size_cX < lV)
    data->lastTol = 0.5;

  for (i=0;i<n;i+=s) {
    /* If the selected eigenvalue is complex, but the arithmetic is real... */
#if !defined(PETSC_USE_COMPLEX)
    if (d->eigi[r_s+i] != 0.0) {
      if (i+2 <= max_size_D) s=2;
      else break;
    } else
#endif
      s=1;

    data->r_s = r_s+i;
    data->r_e = r_s+i+s;
    PetscCall(SlepcVecPoolGetVecs(d->auxV,s,&kr));

    /* Compute theta, maximum iterations and tolerance */
    maxits = 0;
    tol = 1;
    for (j=0;j<s;j++) {
      PetscCall(d->improvex_jd_lit(d,r_s+i+j,&data->theta[2*j],&data->thetai[j],&maxits0,&tol0));
      maxits += maxits0;
      tol *= tol0;
    }
    maxits/= s;
    tol = data->dynamic?data->lastTol:PetscExpReal(PetscLogReal(tol)/s);

    /* Compute u, v and kr */
    k = r_s+i;
    PetscCall(DSVectors(d->eps->ds,DS_MAT_X,&k,NULL));
    k = r_s+i;
    PetscCall(DSVectors(d->eps->ds,DS_MAT_Y,&k,NULL));
    PetscCall(DSGetArray(d->eps->ds,DS_MAT_X,&pX));
    PetscCall(DSGetArray(d->eps->ds,DS_MAT_Y,&pY));
    PetscCall(dvd_improvex_jd_proj_cuv(d,r_s+i,r_s+i+s,kr,data->theta,data->thetai,pX,pY,ld));
    PetscCall(DSRestoreArray(d->eps->ds,DS_MAT_X,&pX));
    PetscCall(DSRestoreArray(d->eps->ds,DS_MAT_Y,&pY));

    /* Check if the first eigenpairs are converged */
    if (i == 0) {
      PetscInt oldnpreconv = d->npreconv;
      PetscCall(d->preTestConv(d,0,r_s+s,r_s+s,&d->npreconv));
      if (d->npreconv > oldnpreconv) break;
    }

    /* Test the odd situation of solving Ax=b with A=I */
#if !defined(PETSC_USE_COMPLEX)
    odd_situation = (data->ksp && data->theta[0] == 1. && data->theta[1] == 0. && data->thetai[0] == 0. && d->B == NULL)? PETSC_TRUE: PETSC_FALSE;
#else
    odd_situation = (data->ksp && data->theta[0] == 1. && data->theta[1] == 0. && d->B == NULL)? PETSC_TRUE: PETSC_FALSE;
#endif
    /* If JD */
    if (data->ksp && !odd_situation) {
      /* kr <- -kr */
      for (j=0;j<s;j++) PetscCall(VecScale(kr[j],-1.0));

      /* Compose kr and D */
      kr0[0] = kr[0];
      kr0[1] = (s==2 ? kr[1] : NULL);
      PetscCall(VecCreateCompWithVecs(kr0,data->ksp_max_size,data->friends,&kr_comp));
      PetscCall(BVGetColumn(d->eps->V,kV+i,&D[0]));
      if (s==2) PetscCall(BVGetColumn(d->eps->V,kV+i+1,&D[1]));
      else D[1] = NULL;
      PetscCall(VecCreateCompWithVecs(D,data->ksp_max_size,data->friends,&D_comp));
      PetscCall(VecCompSetSubVecs(data->friends,s,NULL));

      /* Solve the correction equation */
      PetscCall(KSPSetTolerances(data->ksp,tol,PETSC_CURRENT,PETSC_CURRENT,maxits));
      PetscCall(KSPSolve(data->ksp,kr_comp,D_comp));
      PetscCall(KSPGetIterationNumber(data->ksp,&lits));

      /* Destroy the composed ks and D */
      PetscCall(VecDestroy(&kr_comp));
      PetscCall(VecDestroy(&D_comp));
      PetscCall(BVRestoreColumn(d->eps->V,kV+i,&D[0]));
      if (s==2) PetscCall(BVRestoreColumn(d->eps->V,kV+i+1,&D[1]));

    /* If GD */
    } else {
      PetscCall(BVGetColumn(d->eps->V,kV+i,&D[0]));
      if (s==2) PetscCall(BVGetColumn(d->eps->V,kV+i+1,&D[1]));
      for (j=0;j<s;j++) PetscCall(d->improvex_precond(d,r_s+i+j,kr[j],D[j]));
      PetscCall(dvd_improvex_apply_proj(d,D,s));
      PetscCall(BVRestoreColumn(d->eps->V,kV+i,&D[0]));
      if (s==2) PetscCall(BVRestoreColumn(d->eps->V,kV+i+1,&D[1]));
    }
    /* Prevent that short vectors are discarded in the orthogonalization */
    if (i == 0 && d->eps->errest[d->nconv+r_s] > PETSC_MACHINE_EPSILON && d->eps->errest[d->nconv+r_s] < PETSC_MAX_REAL) {
      for (j=0;j<s;j++) PetscCall(BVScaleColumn(d->eps->V,kV+i+j,1.0/d->eps->errest[d->nconv+r_s]));
    }
    PetscCall(SlepcVecPoolRestoreVecs(d->auxV,s,&kr));
  }
  *size_D = i;
  if (data->dynamic) data->lastTol = PetscMax(data->lastTol/2.0,PETSC_MACHINE_EPSILON*10.0);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvd_improvex_jd(dvdDashboard *d,dvdBlackboard *b,KSP ksp,PetscInt max_bs,PetscBool dynamic)
{
  dvdImprovex_jd *data;
  PetscBool      useGD;
  PC             pc;
  PetscInt       size_P;

  PetscFunctionBegin;
  /* Setting configuration constrains */
  PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&useGD));

  /* If the arithmetic is real and the problem is not Hermitian, then
     the block size is incremented in one */
#if !defined(PETSC_USE_COMPLEX)
  if (!DVD_IS(d->sEP,DVD_EP_HERMITIAN)) {
    max_bs++;
    b->max_size_P = PetscMax(b->max_size_P,2);
  } else
#endif
  {
    b->max_size_P = PetscMax(b->max_size_P,1);
  }
  b->max_size_X = PetscMax(b->max_size_X,max_bs);
  size_P = b->max_size_P;

  /* Setup the preconditioner */
  if (ksp) {
    PetscCall(KSPGetPC(ksp,&pc));
    PetscCall(dvd_static_precond_PC(d,b,pc));
  } else PetscCall(dvd_static_precond_PC(d,b,NULL));

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    PetscCall(PetscNew(&data));
    data->dynamic = dynamic;
    PetscCall(PetscMalloc1(size_P*size_P,&data->XKZ));
    PetscCall(PetscMalloc1(size_P*size_P,&data->iXKZ));
    PetscCall(PetscMalloc1(size_P,&data->iXKZPivots));
    data->ldXKZ = size_P;
    data->size_X = b->max_size_X;
    d->improveX_data = data;
    data->ksp = useGD? NULL: ksp;
    data->d = d;
    d->improveX = dvd_improvex_jd_gen;
#if !defined(PETSC_USE_COMPLEX)
    if (!DVD_IS(d->sEP,DVD_EP_HERMITIAN)) data->ksp_max_size = 2;
    else
#endif
      data->ksp_max_size = 1;
    /* Create various vector basis */
    PetscCall(BVDuplicateResize(d->eps->V,size_P,&data->KZ));
    PetscCall(BVSetMatrix(data->KZ,NULL,PETSC_FALSE));
    PetscCall(BVDuplicate(data->KZ,&data->U));

    PetscCall(EPSDavidsonFLAdd(&d->startList,dvd_improvex_jd_start));
    PetscCall(EPSDavidsonFLAdd(&d->endList,dvd_improvex_jd_end));
    PetscCall(EPSDavidsonFLAdd(&d->destroyList,dvd_improvex_jd_d));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_USE_COMPLEX)
static inline PetscErrorCode dvd_complex_rayleigh_quotient(Vec ur,Vec ui,Vec Axr,Vec Axi,Vec Bxr,Vec Bxi,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscScalar    rAr,iAr,rAi,iAi,rBr,iBr,rBi,iBi,b0,b2,b4,b6,b7;

  PetscFunctionBegin;
  /* eigr = [(rAr+iAi)*(rBr+iBi) + (rAi-iAr)*(rBi-iBr)]/k
     eigi = [(rAi-iAr)*(rBr+iBi) - (rAr+iAi)*(rBi-iBr)]/k
     k    =  (rBr+iBi)*(rBr+iBi) + (rBi-iBr)*(rBi-iBr)    */
  PetscCall(VecDotBegin(Axr,ur,&rAr)); /* r*A*r */
  PetscCall(VecDotBegin(Axr,ui,&iAr)); /* i*A*r */
  PetscCall(VecDotBegin(Axi,ur,&rAi)); /* r*A*i */
  PetscCall(VecDotBegin(Axi,ui,&iAi)); /* i*A*i */
  PetscCall(VecDotBegin(Bxr,ur,&rBr)); /* r*B*r */
  PetscCall(VecDotBegin(Bxr,ui,&iBr)); /* i*B*r */
  PetscCall(VecDotBegin(Bxi,ur,&rBi)); /* r*B*i */
  PetscCall(VecDotBegin(Bxi,ui,&iBi)); /* i*B*i */
  PetscCall(VecDotEnd(Axr,ur,&rAr)); /* r*A*r */
  PetscCall(VecDotEnd(Axr,ui,&iAr)); /* i*A*r */
  PetscCall(VecDotEnd(Axi,ur,&rAi)); /* r*A*i */
  PetscCall(VecDotEnd(Axi,ui,&iAi)); /* i*A*i */
  PetscCall(VecDotEnd(Bxr,ur,&rBr)); /* r*B*r */
  PetscCall(VecDotEnd(Bxr,ui,&iBr)); /* i*B*r */
  PetscCall(VecDotEnd(Bxi,ur,&rBi)); /* r*B*i */
  PetscCall(VecDotEnd(Bxi,ui,&iBi)); /* i*B*i */
  b0 = rAr+iAi; /* rAr+iAi */
  b2 = rAi-iAr; /* rAi-iAr */
  b4 = rBr+iBi; /* rBr+iBi */
  b6 = rBi-iBr; /* rBi-iBr */
  b7 = b4*b4 + b6*b6; /* k */
  *eigr = (b0*b4 + b2*b6) / b7; /* eig_r */
  *eigi = (b2*b4 - b0*b6) / b7; /* eig_i */
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static inline PetscErrorCode dvd_compute_n_rr(PetscInt i_s,PetscInt n,PetscScalar *eigr,PetscScalar *eigi,Vec *u,Vec *Ax,Vec *Bx)
{
  PetscInt       i;
  PetscScalar    b0,b1;

  PetscFunctionBegin;
  for (i=0; i<n; i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (eigi[i_s+i] != 0.0) {
      PetscScalar eigr0=0.0,eigi0=0.0;
      PetscCall(dvd_complex_rayleigh_quotient(u[i],u[i+1],Ax[i],Ax[i+1],Bx[i],Bx[i+1],&eigr0,&eigi0));
      if (PetscAbsScalar(eigr[i_s+i]-eigr0)/PetscAbsScalar(eigr[i_s+i]) > 1e-10 || PetscAbsScalar(eigi[i_s+i]-eigi0)/PetscAbsScalar(eigi[i_s+i]) > 1e-10) PetscCall(PetscInfo(u[0],"The eigenvalue %g%+gi is far from its Rayleigh quotient value %g%+gi\n",(double)eigr[i_s+i],(double)eigi[i_s+i],(double)eigr0,(double)eigi0));
      i++;
    } else
#endif
    {
      PetscCall(VecDotBegin(Ax[i],u[i],&b0));
      PetscCall(VecDotBegin(Bx[i],u[i],&b1));
      PetscCall(VecDotEnd(Ax[i],u[i],&b0));
      PetscCall(VecDotEnd(Bx[i],u[i],&b1));
      b0 = b0/b1;
      if (PetscAbsScalar(eigr[i_s+i]-b0)/PetscAbsScalar(eigr[i_s+i]) > 1e-10) PetscCall(PetscInfo(u[0],"The eigenvalue %g+%g is far from its Rayleigh quotient value %g+%g\n",(double)PetscRealPart(eigr[i_s+i]),(double)PetscImaginaryPart(eigr[i_s+i]),(double)PetscRealPart(b0),(double)PetscImaginaryPart(b0)));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute: u <- X, v <- K*(theta[0]*A+theta[1]*B)*X,
  kr <- K^{-1}*(A-eig*B)*X, being X <- V*pX[i_s..i_e-1], Y <- W*pY[i_s..i_e-1]
  where
  pX,pY, the right and left eigenvectors of the projected system
  ld, the leading dimension of pX and pY
*/
static PetscErrorCode dvd_improvex_jd_proj_uv_KZX(dvdDashboard *d,PetscInt i_s,PetscInt i_e,Vec *u,Vec *v,Vec *kr,PetscScalar *theta,PetscScalar *thetai,PetscScalar *pX,PetscScalar *pY,PetscInt ld)
{
  PetscInt       n = i_e-i_s,i;
  PetscScalar    *b;
  Vec            *Ax,*Bx,*r;
  Mat            M;
  BV             X;

  PetscFunctionBegin;
  PetscCall(BVDuplicateResize(d->eps->V,4,&X));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,4,4,NULL,&M));
  /* u <- X(i) */
  PetscCall(dvd_improvex_compute_X(d,i_s,i_e,u,pX,ld));

  /* v <- theta[0]A*u + theta[1]*B*u */

  /* Bx <- B*X(i) */
  Bx = kr;
  if (d->BX) {
    for (i=i_s; i<i_e; ++i) PetscCall(BVMultVec(d->BX,1.0,0.0,Bx[i-i_s],&pX[ld*i]));
  } else {
    for (i=0;i<n;i++) {
      if (d->B) PetscCall(MatMult(d->B, u[i], Bx[i]));
      else PetscCall(VecCopy(u[i], Bx[i]));
    }
  }

  /* Ax <- A*X(i) */
  PetscCall(SlepcVecPoolGetVecs(d->auxV,i_e-i_s,&r));
  Ax = r;
  for (i=i_s; i<i_e; ++i) PetscCall(BVMultVec(d->AX,1.0,0.0,Ax[i-i_s],&pX[ld*i]));

  /* v <- Y(i) */
  for (i=i_s; i<i_e; ++i) PetscCall(BVMultVec(d->W?d->W:d->eps->V,1.0,0.0,v[i-i_s],&pY[ld*i]));

  /* Recompute the eigenvalue */
  PetscCall(dvd_compute_n_rr(i_s,n,d->eigr,d->eigi,v,Ax,Bx));

  for (i=0;i<n;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (d->eigi[i_s+i] != 0.0) {
      /* [r_i r_i+1 kr_i kr_i+1]*= [ theta_2i'    0            1        0
                                       0         theta_2i'     0        1
                                     theta_2i+1 -thetai_i   -eigr_i -eigi_i
                                     thetai_i    theta_2i+1  eigi_i -eigr_i ] */
      PetscCall(MatDenseGetArrayWrite(M,&b));
      b[0] = b[5] = PetscConj(theta[2*i]);
      b[2] = b[7] = -theta[2*i+1];
      b[6] = -(b[3] = thetai[i]);
      b[1] = b[4] = 0.0;
      b[8] = b[13] = 1.0/d->nX[i_s+i];
      b[10] = b[15] = -d->eigr[i_s+i]/d->nX[i_s+i];
      b[14] = -(b[11] = d->eigi[i_s+i]/d->nX[i_s+i]);
      b[9] = b[12] = 0.0;
      PetscCall(MatDenseRestoreArrayWrite(M,&b));
      PetscCall(BVInsertVec(X,0,Ax[i]));
      PetscCall(BVInsertVec(X,1,Ax[i+1]));
      PetscCall(BVInsertVec(X,2,Bx[i]));
      PetscCall(BVInsertVec(X,3,Bx[i+1]));
      PetscCall(BVSetActiveColumns(X,0,4));
      PetscCall(BVMultInPlace(X,M,0,4));
      PetscCall(BVCopyVec(X,0,Ax[i]));
      PetscCall(BVCopyVec(X,1,Ax[i+1]));
      PetscCall(BVCopyVec(X,2,Bx[i]));
      PetscCall(BVCopyVec(X,3,Bx[i+1]));
      i++;
    } else
#endif
    {
      /* [Ax_i Bx_i]*= [ theta_2i'    1/nX_i
                        theta_2i+1  -eig_i/nX_i ] */
      PetscCall(MatDenseGetArrayWrite(M,&b));
      b[0] = PetscConj(theta[i*2]);
      b[1] = theta[i*2+1];
      b[4] = 1.0/d->nX[i_s+i];
      b[5] = -d->eigr[i_s+i]/d->nX[i_s+i];
      PetscCall(MatDenseRestoreArrayWrite(M,&b));
      PetscCall(BVInsertVec(X,0,Ax[i]));
      PetscCall(BVInsertVec(X,1,Bx[i]));
      PetscCall(BVSetActiveColumns(X,0,2));
      PetscCall(BVMultInPlace(X,M,0,2));
      PetscCall(BVCopyVec(X,0,Ax[i]));
      PetscCall(BVCopyVec(X,1,Bx[i]));
    }
  }
  for (i=0; i<n; i++) d->nX[i_s+i] = 1.0;

  /* v <- K^{-1} r = K^{-1}(theta_2i'*Ax + theta_2i+1*Bx) */
  for (i=0;i<n;i++) PetscCall(d->improvex_precond(d,i_s+i,r[i],v[i]));
  PetscCall(SlepcVecPoolRestoreVecs(d->auxV,i_e-i_s,&r));

  /* kr <- P*(Ax - eig_i*Bx) */
  PetscCall(d->calcpairs_proj_res(d,i_s,i_e,kr));
  PetscCall(BVDestroy(&X));
  PetscCall(MatDestroy(&M));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_improvex_jd_lit_const_0(dvdDashboard *d,PetscInt i,PetscScalar* theta,PetscScalar* thetai,PetscInt *maxits,PetscReal *tol)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscReal      a;

  PetscFunctionBegin;
  a = SlepcAbsEigenvalue(d->eigr[i],d->eigi[i]);

  if (d->nR[i] < data->fix*a) {
    theta[0] = d->eigr[i];
    theta[1] = 1.0;
#if !defined(PETSC_USE_COMPLEX)
    *thetai = d->eigi[i];
#endif
  } else {
    theta[0] = d->target[0];
    theta[1] = d->target[1];
#if !defined(PETSC_USE_COMPLEX)
    *thetai = 0.0;
#endif
}

#if defined(PETSC_USE_COMPLEX)
  if (thetai) *thetai = 0.0;
#endif
  *maxits = data->maxits;
  *tol = data->tol;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvd_improvex_jd_lit_const(dvdDashboard *d,dvdBlackboard *b,PetscInt maxits,PetscReal tol,PetscReal fix)
{
  dvdImprovex_jd  *data = (dvdImprovex_jd*)d->improveX_data;

  PetscFunctionBegin;
  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    data->maxits = maxits;
    data->tol = tol;
    data->fix = fix;
    d->improvex_jd_lit = dvd_improvex_jd_lit_const_0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvd_improvex_jd_proj_uv(dvdDashboard *d,dvdBlackboard *b)
{
  PetscFunctionBegin;
  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    d->improvex_jd_proj_uv = dvd_improvex_jd_proj_uv_KZX;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvd_improvex_compute_X(dvdDashboard *d,PetscInt i_s,PetscInt i_e,Vec *u_,PetscScalar *pX,PetscInt ld)
{
  PetscInt       n = i_e - i_s,i;
  Vec            *u;

  PetscFunctionBegin;
  if (u_) u = u_;
  else if (d->correctXnorm) PetscCall(SlepcVecPoolGetVecs(d->auxV,i_e-i_s,&u));
  if (u_ || d->correctXnorm) {
    for (i=0;i<n;i++) PetscCall(BVMultVec(d->eps->V,1.0,0.0,u[i],&pX[ld*(i+i_s)]));
  }
  /* nX(i) <- ||X(i)|| */
  if (d->correctXnorm) {
    for (i=0;i<n;i++) PetscCall(VecNormBegin(u[i],NORM_2,&d->nX[i_s+i]));
    for (i=0;i<n;i++) PetscCall(VecNormEnd(u[i],NORM_2,&d->nX[i_s+i]));
#if !defined(PETSC_USE_COMPLEX)
    for (i=0;i<n;i++) {
      if (d->eigi[i_s+i] != 0.0) {
        d->nX[i_s+i] = d->nX[i_s+i+1] = PetscSqrtScalar(d->nX[i_s+i]*d->nX[i_s+i]+d->nX[i_s+i+1]*d->nX[i_s+i+1]);
        i++;
      }
    }
#endif
  } else {
    for (i=0;i<n;i++) d->nX[i_s+i] = 1.0;
  }
  if (d->correctXnorm && !u_) PetscCall(SlepcVecPoolRestoreVecs(d->auxV,i_e-i_s,&u));
  PetscFunctionReturn(PETSC_SUCCESS);
}
