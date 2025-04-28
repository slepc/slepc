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

   Step: calculate the best eigenpairs in the subspace V

   For that, performs these steps:
     1) Update W <- A * V
     2) Update H <- V' * W
     3) Obtain eigenpairs of H
     4) Select some eigenpairs
     5) Compute the Ritz pairs of the selected ones
*/

#include "davidson.h"
#include <slepcblaslapack.h>

static PetscErrorCode dvd_calcpairs_qz_start(dvdDashboard *d)
{
  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(d->eps->V,0,0));
  if (d->W) PetscCall(BVSetActiveColumns(d->W,0,0));
  PetscCall(BVSetActiveColumns(d->AX,0,0));
  if (d->BX) PetscCall(BVSetActiveColumns(d->BX,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_calcpairs_qz_d(dvdDashboard *d)
{
  PetscFunctionBegin;
  PetscCall(BVDestroy(&d->W));
  PetscCall(BVDestroy(&d->AX));
  PetscCall(BVDestroy(&d->BX));
  PetscCall(BVDestroy(&d->auxBV));
  PetscCall(MatDestroy(&d->H));
  PetscCall(MatDestroy(&d->G));
  PetscCall(MatDestroy(&d->auxM));
  PetscCall(SlepcVecPoolDestroy(&d->auxV));
  PetscCall(PetscFree(d->nBds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* in complex, d->size_H real auxiliary values are needed */
static PetscErrorCode dvd_calcpairs_projeig_solve(dvdDashboard *d)
{
  Vec               v;
  Mat               A,B,H0,G0;
  PetscScalar       *pA;
  const PetscScalar *pv;
  PetscInt          i,lV,kV,n,ld;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&lV,&kV));
  n = kV-lV;
  if (!n) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DSSetDimensions(d->eps->ds,n,0,0));
  PetscCall(DSGetMat(d->eps->ds,DS_MAT_A,&A));
  PetscCall(MatDenseGetSubMatrix(d->H,lV,lV+n,lV,lV+n,&H0));
  PetscCall(MatCopy(H0,A,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseRestoreSubMatrix(d->H,&H0));
  PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_A,&A));
  if (d->G) {
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_B,&B));
    PetscCall(MatDenseGetSubMatrix(d->G,lV,lV+n,lV,lV+n,&G0));
    PetscCall(MatCopy(G0,B,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(d->G,&G0));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_B,&B));
  }
  /* Set the signature on projected matrix B */
  if (DVD_IS(d->sEP,DVD_EP_INDEFINITE)) {
    PetscCall(DSGetLeadingDimension(d->eps->ds,&ld));
    PetscCall(DSGetArray(d->eps->ds,DS_MAT_B,&pA));
    PetscCall(PetscArrayzero(pA,n*ld));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF,kV,&v));
    PetscCall(BVGetSignature(d->eps->V,v));
    PetscCall(VecGetArrayRead(v,&pv));
    for (i=0;i<n;i++) {
      pA[i+ld*i] = d->nBds[i] = PetscRealPart(pv[lV+i]);
    }
    PetscCall(VecRestoreArrayRead(v,&pv));
    PetscCall(VecDestroy(&v));
    PetscCall(DSRestoreArray(d->eps->ds,DS_MAT_B,&pA));
  }
  PetscCall(DSSetState(d->eps->ds,DS_STATE_RAW));
  PetscCall(DSSolve(d->eps->ds,d->eigr,d->eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   A(lA:kA-1,lA:kA-1) <- Z(l:k-1)'*A(l:k-1,l:k-1)*Q(l,k-1), where k=l+kA-lA
 */
static PetscErrorCode EPSXDUpdateProj(Mat Q,Mat Z,PetscInt l,Mat A,PetscInt lA,PetscInt kA,Mat aux)
{
  PetscScalar       one=1.0,zero=0.0;
  PetscInt          i,j,dA_=kA-lA,m0,n0,ldA_,ldQ_,ldZ_,nQ_;
  PetscBLASInt      dA,nQ,ldA,ldQ,ldZ;
  PetscScalar       *pA,*pW;
  const PetscScalar *pQ,*pZ;
  PetscBool         symm=PETSC_FALSE,set,flg;

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,&m0,&n0));
  PetscCall(MatDenseGetLDA(A,&ldA_));
  PetscAssert(m0==n0,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"A should be square");
  PetscAssert(lA>=0 && lA<=m0,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Invalid initial row, column in A");
  PetscAssert(kA>=0 && kA>=lA && kA<=m0,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Invalid final row, column in A");
  PetscCall(MatIsHermitianKnown(A,&set,&flg));
  symm = set? flg: PETSC_FALSE;
  PetscCall(MatGetSize(Q,&m0,&n0)); nQ_=m0;
  PetscCall(MatDenseGetLDA(Q,&ldQ_));
  PetscAssert(l>=0 && l<=n0 && l+dA_<=n0,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Invalid initial column in Q");
  PetscCall(MatGetSize(Z,&m0,&n0));
  PetscCall(MatDenseGetLDA(Z,&ldZ_));
  PetscAssert(l>=0 && l<=n0 && l+dA_<=n0,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Invalid initial column in Z");
  PetscCall(MatGetSize(aux,&m0,&n0));
  PetscAssert(m0*n0>=nQ_*dA_,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"aux should be larger");
  PetscCall(PetscBLASIntCast(dA_,&dA));
  PetscCall(PetscBLASIntCast(nQ_,&nQ));
  PetscCall(PetscBLASIntCast(ldA_,&ldA));
  PetscCall(PetscBLASIntCast(ldQ_,&ldQ));
  PetscCall(PetscBLASIntCast(ldZ_,&ldZ));
  PetscCall(MatDenseGetArray(A,&pA));
  PetscCall(MatDenseGetArrayRead(Q,&pQ));
  if (Q!=Z) PetscCall(MatDenseGetArrayRead(Z,&pZ));
  else pZ = pQ;
  PetscCall(MatDenseGetArrayWrite(aux,&pW));
  /* W = A*Q */
  if (symm) {
    /* symmetrize before multiplying */
    for (i=lA+1;i<lA+nQ;i++) {
      for (j=lA;j<i;j++) pA[i+j*ldA] = PetscConj(pA[j+i*ldA]);
    }
  }
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&nQ,&dA,&nQ,&one,&pA[ldA*lA+lA],&ldA,&pQ[ldQ*l+l],&ldQ,&zero,pW,&nQ));
  /* A = Q'*W */
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&dA,&dA,&nQ,&one,&pZ[ldZ*l+l],&ldZ,pW,&nQ,&zero,&pA[ldA*lA+lA],&ldA));
  PetscCall(MatDenseRestoreArray(A,&pA));
  PetscCall(MatDenseRestoreArrayRead(Q,&pQ));
  if (Q!=Z) PetscCall(MatDenseRestoreArrayRead(Z,&pZ));
  PetscCall(MatDenseRestoreArrayWrite(aux,&pW));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_calcpairs_updateproj(dvdDashboard *d)
{
  Mat            Q,Z;
  PetscInt       lV,kV;
  PetscBool      symm;

  PetscFunctionBegin;
  PetscCall(DSGetMat(d->eps->ds,DS_MAT_Q,&Q));
  if (d->W) PetscCall(DSGetMat(d->eps->ds,DS_MAT_Z,&Z));
  else Z = Q;
  PetscCall(BVGetActiveColumns(d->eps->V,&lV,&kV));
  PetscCall(EPSXDUpdateProj(Q,Z,0,d->H,lV,lV+d->V_tra_e,d->auxM));
  if (d->G) PetscCall(EPSXDUpdateProj(Q,Z,0,d->G,lV,lV+d->V_tra_e,d->auxM));
  PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Q,&Q));
  if (d->W) PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Z,&Z));

  PetscCall(PetscObjectTypeCompareAny((PetscObject)d->eps->ds,&symm,DSHEP,DSGHIEP,DSGHEP,""));
  if (d->V_tra_s==0 || symm) PetscFunctionReturn(PETSC_SUCCESS);
  /* Compute upper part of H (and G): H(0:l-1,l:k-1) <- W(0:l-1)' * AV(l:k-1), where
     k=l+d->V_tra_s */
  PetscCall(BVSetActiveColumns(d->W?d->W:d->eps->V,0,lV));
  PetscCall(BVSetActiveColumns(d->AX,lV,lV+d->V_tra_s));
  PetscCall(BVDot(d->AX,d->W?d->W:d->eps->V,d->H));
  if (d->G) {
    PetscCall(BVSetActiveColumns(d->BX?d->BX:d->eps->V,lV,lV+d->V_tra_s));
    PetscCall(BVDot(d->BX?d->BX:d->eps->V,d->W?d->W:d->eps->V,d->G));
  }
  PetscCall(PetscObjectTypeCompare((PetscObject)d->eps->ds,DSGHEP,&symm));
  if (!symm) {
    PetscCall(BVSetActiveColumns(d->W?d->W:d->eps->V,lV,lV+d->V_tra_s));
    PetscCall(BVSetActiveColumns(d->AX,0,lV));
    PetscCall(BVDot(d->AX,d->W?d->W:d->eps->V,d->H));
    if (d->G) {
      PetscCall(BVSetActiveColumns(d->BX?d->BX:d->eps->V,0,lV));
      PetscCall(BVDot(d->BX?d->BX:d->eps->V,d->W?d->W:d->eps->V,d->G));
    }
  }
  PetscCall(BVSetActiveColumns(d->eps->V,lV,kV));
  PetscCall(BVSetActiveColumns(d->AX,lV,kV));
  if (d->BX) PetscCall(BVSetActiveColumns(d->BX,lV,kV));
  if (d->W) PetscCall(BVSetActiveColumns(d->W,lV,kV));
  if (d->W) PetscCall(dvd_harm_updateproj(d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV <- BV*MT
 */
static inline PetscErrorCode dvd_calcpairs_updateBV0_gen(dvdDashboard *d,BV bv,DSMatType mat)
{
  PetscInt       l,k,n;
  Mat            M,M0,auxM,auxM0;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&l,&k));
  PetscCall(DSGetDimensions(d->eps->ds,&n,NULL,NULL,NULL));
  PetscAssert(k-l==n,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");
  PetscCall(DSGetMat(d->eps->ds,mat,&M));
  PetscCall(MatDenseGetSubMatrix(M,0,n,0,d->V_tra_e,&M0));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&auxM));
  PetscCall(MatDenseGetSubMatrix(auxM,l,l+n,l,l+d->V_tra_e,&auxM0));
  PetscCall(MatCopy(M0,auxM0,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseRestoreSubMatrix(auxM,&auxM0));
  PetscCall(MatDenseRestoreSubMatrix(M,&M0));
  PetscCall(DSRestoreMat(d->eps->ds,mat,&M));
  PetscCall(BVMultInPlace(bv,auxM,l,l+d->V_tra_e));
  PetscCall(MatDestroy(&auxM));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_calcpairs_proj(dvdDashboard *d)
{
  PetscInt       i,l,k;
  Vec            v1,v2;
  PetscScalar    *pv;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&l,&k));
  /* Update AV, BV, W and the projected matrices */
  /* 1. S <- S*MT */
  if (d->V_tra_s != d->V_tra_e || d->V_tra_e > 0) {
    PetscCall(dvd_calcpairs_updateBV0_gen(d,d->eps->V,DS_MAT_Q));
    if (d->W) PetscCall(dvd_calcpairs_updateBV0_gen(d,d->W,DS_MAT_Z));
    PetscCall(dvd_calcpairs_updateBV0_gen(d,d->AX,DS_MAT_Q));
    if (d->BX) PetscCall(dvd_calcpairs_updateBV0_gen(d,d->BX,DS_MAT_Q));
    PetscCall(dvd_calcpairs_updateproj(d));
    /* Update signature */
    if (d->nBds) {
      PetscCall(VecCreateSeq(PETSC_COMM_SELF,l+d->V_tra_e,&v1));
      PetscCall(BVSetActiveColumns(d->eps->V,0,l+d->V_tra_e));
      PetscCall(BVGetSignature(d->eps->V,v1));
      PetscCall(VecGetArray(v1,&pv));
      for (i=0;i<d->V_tra_e;i++) pv[l+i] = d->nBds[i];
      PetscCall(VecRestoreArray(v1,&pv));
      PetscCall(BVSetSignature(d->eps->V,v1));
      PetscCall(BVSetActiveColumns(d->eps->V,l,k));
      PetscCall(VecDestroy(&v1));
    }
    k = l+d->V_tra_e;
    l+= d->V_tra_s;
  } else {
    /* 2. V <- orth(V, V_new) */
    PetscCall(dvd_orthV(d->eps->V,l+d->V_new_s,l+d->V_new_e));
    /* 3. AV <- [AV A * V(V_new_s:V_new_e-1)] */
    /* Check consistency */
    PetscAssert(k-l==d->V_new_s,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");
    for (i=l+d->V_new_s;i<l+d->V_new_e;i++) {
      PetscCall(BVGetColumn(d->eps->V,i,&v1));
      PetscCall(BVGetColumn(d->AX,i,&v2));
      PetscCall(MatMult(d->A,v1,v2));
      PetscCall(BVRestoreColumn(d->eps->V,i,&v1));
      PetscCall(BVRestoreColumn(d->AX,i,&v2));
    }
    /* 4. BV <- [BV B * V(V_new_s:V_new_e-1)] */
    if (d->BX) {
      /* Check consistency */
      PetscAssert(k-l==d->V_new_s,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");
      for (i=l+d->V_new_s;i<l+d->V_new_e;i++) {
        PetscCall(BVGetColumn(d->eps->V,i,&v1));
        PetscCall(BVGetColumn(d->BX,i,&v2));
        PetscCall(MatMult(d->B,v1,v2));
        PetscCall(BVRestoreColumn(d->eps->V,i,&v1));
        PetscCall(BVRestoreColumn(d->BX,i,&v2));
      }
    }
    /* 5. W <- [W f(AV,BV)] */
    if (d->W) {
      PetscCall(d->calcpairs_W(d));
      PetscCall(dvd_orthV(d->W,l+d->V_new_s,l+d->V_new_e));
    }
    /* 6. H <- W' * AX; G <- W' * BX */
    PetscCall(BVSetActiveColumns(d->eps->V,l+d->V_new_s,l+d->V_new_e));
    PetscCall(BVSetActiveColumns(d->AX,l+d->V_new_s,l+d->V_new_e));
    if (d->BX) PetscCall(BVSetActiveColumns(d->BX,l+d->V_new_s,l+d->V_new_e));
    if (d->W) PetscCall(BVSetActiveColumns(d->W,l+d->V_new_s,l+d->V_new_e));
    PetscCall(BVMatProject(d->AX,NULL,d->W?d->W:d->eps->V,d->H));
    if (d->G) PetscCall(BVMatProject(d->BX?d->BX:d->eps->V,NULL,d->W?d->W:d->eps->V,d->G));
    PetscCall(BVSetActiveColumns(d->eps->V,l,k));
    PetscCall(BVSetActiveColumns(d->AX,l,k));
    if (d->BX) PetscCall(BVSetActiveColumns(d->BX,l,k));
    if (d->W) PetscCall(BVSetActiveColumns(d->W,l,k));

    /* Perform the transformation on the projected problem */
    if (d->W) PetscCall(d->calcpairs_proj_trans(d));
    k = l+d->V_new_e;
  }
  PetscCall(BVSetActiveColumns(d->eps->V,l,k));
  PetscCall(BVSetActiveColumns(d->AX,l,k));
  if (d->BX) PetscCall(BVSetActiveColumns(d->BX,l,k));
  if (d->W) PetscCall(BVSetActiveColumns(d->W,l,k));

  /* Solve the projected problem */
  PetscCall(dvd_calcpairs_projeig_solve(d));

  d->V_tra_s = d->V_tra_e = 0;
  d->V_new_s = d->V_new_e;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_calcpairs_apply_arbitrary(dvdDashboard *d,PetscInt r_s,PetscInt r_e,PetscScalar *rr,PetscScalar *ri)
{
  PetscInt       i,k,ld;
  PetscScalar    *pX;
  Vec            *X,xr,xi;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       N=1;
#else
  PetscInt       N=2,j;
#endif

  PetscFunctionBegin;
  /* Quick exit without neither arbitrary selection nor harmonic extraction */
  if (!d->eps->arbitrary && !d->calcpairs_eig_backtrans) PetscFunctionReturn(PETSC_SUCCESS);

  /* Quick exit without arbitrary selection, but with harmonic extraction */
  if (d->calcpairs_eig_backtrans) {
    for (i=r_s; i<r_e; i++) PetscCall(d->calcpairs_eig_backtrans(d,d->eigr[i],d->eigi[i],&rr[i-r_s],&ri[i-r_s]));
  }
  if (!d->eps->arbitrary) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(SlepcVecPoolGetVecs(d->auxV,N,&X));
  PetscCall(DSGetLeadingDimension(d->eps->ds,&ld));
  for (i=r_s;i<r_e;i++) {
    k = i;
    PetscCall(DSVectors(d->eps->ds,DS_MAT_X,&k,NULL));
    PetscCall(DSGetArray(d->eps->ds,DS_MAT_X,&pX));
    PetscCall(dvd_improvex_compute_X(d,i,k+1,X,pX,ld));
    PetscCall(DSRestoreArray(d->eps->ds,DS_MAT_X,&pX));
#if !defined(PETSC_USE_COMPLEX)
    if (d->nX[i] != 1.0) {
      for (j=i;j<k+1;j++) PetscCall(VecScale(X[j-i],1.0/d->nX[i]));
    }
    xr = X[0];
    xi = X[1];
    if (i == k) PetscCall(VecSet(xi,0.0));
#else
    xr = X[0];
    xi = NULL;
    if (d->nX[i] != 1.0) PetscCall(VecScale(xr,1.0/d->nX[i]));
#endif
    PetscCall(d->eps->arbitrary(rr[i-r_s],ri[i-r_s],xr,xi,&rr[i-r_s],&ri[i-r_s],d->eps->arbitraryctx));
#if !defined(PETSC_USE_COMPLEX)
    if (i != k) {
      rr[i+1-r_s] = rr[i-r_s];
      ri[i+1-r_s] = ri[i-r_s];
      i++;
    }
#endif
  }
  PetscCall(SlepcVecPoolRestoreVecs(d->auxV,N,&X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_calcpairs_selectPairs(dvdDashboard *d,PetscInt n)
{
  PetscInt       k,lV,kV,nV;
  PetscScalar    *rr,*ri;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&lV,&kV));
  nV = kV - lV;
  n = PetscMin(n,nV);
  if (n <= 0) PetscFunctionReturn(PETSC_SUCCESS);
  /* Put the best n pairs at the beginning. Useful for restarting */
  if (d->eps->arbitrary || d->calcpairs_eig_backtrans) {
    PetscCall(PetscMalloc1(nV,&rr));
    PetscCall(PetscMalloc1(nV,&ri));
    PetscCall(dvd_calcpairs_apply_arbitrary(d,0,nV,rr,ri));
  } else {
    rr = d->eigr;
    ri = d->eigi;
  }
  k = n;
  PetscCall(DSSort(d->eps->ds,d->eigr,d->eigi,rr,ri,&k));
  /* Put the best pair at the beginning. Useful to check its residual */
#if !defined(PETSC_USE_COMPLEX)
  if (n != 1 && (n != 2 || d->eigi[0] == 0.0))
#else
  if (n != 1)
#endif
  {
    PetscCall(dvd_calcpairs_apply_arbitrary(d,0,nV,rr,ri));
    k = 1;
    PetscCall(DSSort(d->eps->ds,d->eigr,d->eigi,rr,ri,&k));
  }
  PetscCall(DSSynchronize(d->eps->ds,d->eigr,d->eigi));

  if (d->calcpairs_eigs_trans) PetscCall(d->calcpairs_eigs_trans(d));
  if (d->eps->arbitrary || d->calcpairs_eig_backtrans) {
    PetscCall(PetscFree(rr));
    PetscCall(PetscFree(ri));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSXDComputeDSConv(dvdDashboard *d)
{
  PetscInt          i,ld;
  Vec               v;
  Mat               A,B,H0,G0;
  PetscScalar       *pA;
  const PetscScalar *pv;
  PetscBool         symm;

  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(d->eps->V,0,d->eps->nconv));
  PetscCall(PetscObjectTypeCompare((PetscObject)d->eps->ds,DSHEP,&symm));
  if (symm || !d->eps->nconv) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DSSetDimensions(d->eps->ds,d->eps->nconv,0,0));
  PetscCall(DSGetMat(d->eps->ds,DS_MAT_A,&A));
  PetscCall(MatDenseGetSubMatrix(d->H,0,d->eps->nconv,0,d->eps->nconv,&H0));
  PetscCall(MatCopy(H0,A,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseRestoreSubMatrix(d->H,&H0));
  PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_A,&A));
  if (d->G) {
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_B,&B));
    PetscCall(MatDenseGetSubMatrix(d->G,0,d->eps->nconv,0,d->eps->nconv,&G0));
    PetscCall(MatCopy(G0,B,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(d->G,&G0));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_B,&B));
  }
  /* Set the signature on projected matrix B */
  if (DVD_IS(d->sEP,DVD_EP_INDEFINITE)) {
    PetscCall(DSGetLeadingDimension(d->eps->ds,&ld));
    PetscCall(DSGetArray(d->eps->ds,DS_MAT_B,&pA));
    PetscCall(PetscArrayzero(pA,d->eps->nconv*ld));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF,d->eps->nconv,&v));
    PetscCall(BVGetSignature(d->eps->V,v));
    PetscCall(VecGetArrayRead(v,&pv));
    for (i=0;i<d->eps->nconv;i++) pA[i+ld*i] = pv[i];
    PetscCall(VecRestoreArrayRead(v,&pv));
    PetscCall(VecDestroy(&v));
    PetscCall(DSRestoreArray(d->eps->ds,DS_MAT_B,&pA));
  }
  PetscCall(DSSetState(d->eps->ds,DS_STATE_RAW));
  PetscCall(DSSolve(d->eps->ds,d->eps->eigr,d->eps->eigi));
  PetscCall(DSSynchronize(d->eps->ds,d->eps->eigr,d->eps->eigi));
  if (d->W) {
    for (i=0;i<d->eps->nconv;i++) PetscCall(d->calcpairs_eig_backtrans(d,d->eps->eigr[i],d->eps->eigi[i],&d->eps->eigr[i],&d->eps->eigi[i]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Compute the residual vectors R(i) <- (AV - BV*eigr(i))*pX(i), and also
   the norm associated to the Schur pair, where i = r_s..r_e
*/
static PetscErrorCode dvd_calcpairs_res_0(dvdDashboard *d,PetscInt r_s,PetscInt r_e)
{
  PetscInt       i,ldpX;
  PetscScalar    *pX;
  BV             BX = d->BX?d->BX:d->eps->V;
  Vec            *R;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(d->eps->ds,&ldpX));
  PetscCall(DSGetArray(d->eps->ds,DS_MAT_Q,&pX));
  /* nX(i) <- ||X(i)|| */
  PetscCall(dvd_improvex_compute_X(d,r_s,r_e,NULL,pX,ldpX));
  PetscCall(SlepcVecPoolGetVecs(d->auxV,r_e-r_s,&R));
  for (i=r_s;i<r_e;i++) {
    /* R(i-r_s) <- AV*pX(i) */
    PetscCall(BVMultVec(d->AX,1.0,0.0,R[i-r_s],&pX[ldpX*i]));
    /* R(i-r_s) <- R(i-r_s) - eigr(i)*BX*pX(i) */
    PetscCall(BVMultVec(BX,-d->eigr[i],1.0,R[i-r_s],&pX[ldpX*i]));
  }
  PetscCall(DSRestoreArray(d->eps->ds,DS_MAT_Q,&pX));
  PetscCall(d->calcpairs_proj_res(d,r_s,r_e,R));
  PetscCall(SlepcVecPoolRestoreVecs(d->auxV,r_e-r_s,&R));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_calcpairs_proj_res(dvdDashboard *d,PetscInt r_s,PetscInt r_e,Vec *R)
{
  PetscInt       i,l,k;
  PetscBool      lindep=PETSC_FALSE;
  BV             cX;

  PetscFunctionBegin;
  if (d->W) cX = d->W; /* If left subspace exists, R <- orth(cY, R), nR[i] <- ||R[i]|| */
  else if (!(DVD_IS(d->sEP, DVD_EP_STD) && DVD_IS(d->sEP, DVD_EP_HERMITIAN))) cX = d->eps->V; /* If not HEP, R <- orth(cX, R), nR[i] <- ||R[i]|| */
  else cX = NULL; /* Otherwise, nR[i] <- ||R[i]|| */

  if (cX) {
    PetscCall(BVGetActiveColumns(cX,&l,&k));
    PetscCall(BVSetActiveColumns(cX,0,l));
    for (i=0;i<r_e-r_s;i++) PetscCall(BVOrthogonalizeVec(cX,R[i],NULL,&d->nR[r_s+i],&lindep));
    PetscCall(BVSetActiveColumns(cX,l,k));
    if (lindep || (PetscAbs(d->nR[r_s+i]) < PETSC_MACHINE_EPSILON)) PetscCall(PetscInfo(d->eps,"The computed eigenvector residual %" PetscInt_FMT " is too low, %g!\n",r_s+i,(double)d->nR[r_s+i]));
  } else {
    for (i=0;i<r_e-r_s;i++) PetscCall(VecNormBegin(R[i],NORM_2,&d->nR[r_s+i]));
    for (i=0;i<r_e-r_s;i++) PetscCall(VecNormEnd(R[i],NORM_2,&d->nR[r_s+i]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvd_calcpairs_qz(dvdDashboard *d,dvdBlackboard *b,PetscBool borth,PetscBool harm)
{
  PetscBool      std_probl,her_probl,ind_probl;
  DSType         dstype;
  Vec            v1;

  PetscFunctionBegin;
  std_probl = DVD_IS(d->sEP,DVD_EP_STD)? PETSC_TRUE: PETSC_FALSE;
  her_probl = DVD_IS(d->sEP,DVD_EP_HERMITIAN)? PETSC_TRUE: PETSC_FALSE;
  ind_probl = DVD_IS(d->sEP,DVD_EP_INDEFINITE)? PETSC_TRUE: PETSC_FALSE;

  /* Setting configuration constrains */
  b->max_size_proj = PetscMax(b->max_size_proj,b->max_size_V);
  d->W_shift = d->B? PETSC_TRUE: PETSC_FALSE;

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    d->max_size_P = b->max_size_P;
    d->max_size_proj = b->max_size_proj;
    /* Create a DS if the method works with Schur decompositions */
    d->calcPairs = dvd_calcpairs_proj;
    d->calcpairs_residual = dvd_calcpairs_res_0;
    d->calcpairs_proj_res = dvd_calcpairs_proj_res;
    d->calcpairs_selectPairs = dvd_calcpairs_selectPairs;
    /* Create and configure a DS for solving the projected problems */
    if (d->W) dstype = DSGNHEP;    /* If we use harmonics */
    else {
      if (ind_probl) dstype = DSGHIEP;
      else if (std_probl) dstype = her_probl? DSHEP : DSNHEP;
      else dstype = her_probl? DSGHEP : DSGNHEP;
    }
    PetscCall(DSSetType(d->eps->ds,dstype));
    PetscCall(DSAllocate(d->eps->ds,d->eps->ncv));
    /* Create various vector basis */
    if (harm) {
      PetscCall(BVDuplicateResize(d->eps->V,d->eps->ncv,&d->W));
      PetscCall(BVSetMatrix(d->W,NULL,PETSC_FALSE));
    } else d->W = NULL;
    PetscCall(BVDuplicateResize(d->eps->V,d->eps->ncv,&d->AX));
    PetscCall(BVSetMatrix(d->AX,NULL,PETSC_FALSE));
    PetscCall(BVDuplicateResize(d->eps->V,d->eps->ncv,&d->auxBV));
    PetscCall(BVSetMatrix(d->auxBV,NULL,PETSC_FALSE));
    if (d->B) {
      PetscCall(BVDuplicateResize(d->eps->V,d->eps->ncv,&d->BX));
      PetscCall(BVSetMatrix(d->BX,NULL,PETSC_FALSE));
    } else d->BX = NULL;
    PetscCall(MatCreateVecsEmpty(d->A,&v1,NULL));
    PetscCall(SlepcVecPoolCreate(v1,0,&d->auxV));
    PetscCall(VecDestroy(&v1));
    /* Create projected problem matrices */
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,d->eps->ncv,d->eps->ncv,NULL,&d->H));
    if (!std_probl) PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,d->eps->ncv,d->eps->ncv,NULL,&d->G));
    else d->G = NULL;
    if (her_probl) {
      PetscCall(MatSetOption(d->H,MAT_HERMITIAN,PETSC_TRUE));
      if (d->G) PetscCall(MatSetOption(d->G,MAT_HERMITIAN,PETSC_TRUE));
    }

    if (ind_probl) PetscCall(PetscMalloc1(d->eps->ncv,&d->nBds));
    else d->nBds = NULL;
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,d->eps->ncv,d->eps->ncv,NULL,&d->auxM));

    PetscCall(EPSDavidsonFLAdd(&d->startList,dvd_calcpairs_qz_start));
    PetscCall(EPSDavidsonFLAdd(&d->endList,EPSXDComputeDSConv));
    PetscCall(EPSDavidsonFLAdd(&d->destroyList,dvd_calcpairs_qz_d));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
