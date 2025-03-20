/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc polynomial eigensolver: "toar"

   Method: TOAR

   Algorithm:

       Two-Level Orthogonal Arnoldi.

   References:

       [1] Y. Su, J. Zhang and Z. Bai, "A compact Arnoldi algorithm for
           polynomial eigenvalue problems", talk presented at RANMEP 2008.

       [2] C. Campos and J.E. Roman, "Parallel Krylov solvers for the
           polynomial eigenvalue problem in SLEPc", SIAM J. Sci. Comput.
           38(5):S385-S411, 2016.

       [3] D. Lu, Y. Su and Z. Bai, "Stability analysis of the two-level
           orthogonal Arnoldi procedure", SIAM J. Matrix Anal. App.
           37(1):195-214, 2016.
*/

#include <slepc/private/pepimpl.h>    /*I "slepcpep.h" I*/
#include "../src/pep/impls/krylov/pepkrylov.h"
#include <slepcblaslapack.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-pep,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"Parallel {Krylov} solvers for the polynomial eigenvalue problem in {SLEPc}\",\n"
  "   journal = \"{SIAM} J. Sci. Comput.\",\n"
  "   volume = \"38\",\n"
  "   number = \"5\",\n"
  "   pages = \"S385--S411\",\n"
  "   year = \"2016,\"\n"
  "   doi = \"https://doi.org/10.1137/15M1022458\"\n"
  "}\n";

static PetscErrorCode PEPSetUp_TOAR(PEP pep)
{
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  PetscBool      sinv,flg;
  PetscInt       i;

  PetscFunctionBegin;
  PEPCheckShiftSinvert(pep);
  PetscCall(PEPSetDimensions_Default(pep,pep->nev,&pep->ncv,&pep->mpd));
  PetscCheck(ctx->lock || pep->mpd>=pep->ncv,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Should not use mpd parameter in non-locking variant");
  if (pep->max_it==PETSC_DETERMINE) pep->max_it = PetscMax(100,2*(pep->nmat-1)*pep->n/pep->ncv);
  if (!pep->which) PetscCall(PEPSetWhichEigenpairs_Default(pep));
  PetscCheck(pep->which!=PEP_ALL,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues");
  if (pep->problem_type!=PEP_GENERAL) PetscCall(PetscInfo(pep,"Problem type ignored, performing a non-symmetric linearization\n"));

  if (!ctx->keep) ctx->keep = 0.5;

  PetscCall(PEPAllocateSolution(pep,pep->nmat-1));
  PetscCall(PEPSetWorkVecs(pep,3));
  PetscCall(DSSetType(pep->ds,DSNHEP));
  PetscCall(DSSetExtraRow(pep->ds,PETSC_TRUE));
  PetscCall(DSAllocate(pep->ds,pep->ncv+1));

  PetscCall(PEPBasisCoefficients(pep,pep->pbc));
  PetscCall(STGetTransform(pep->st,&flg));
  if (!flg) {
    PetscCall(PetscFree(pep->solvematcoeffs));
    PetscCall(PetscMalloc1(pep->nmat,&pep->solvematcoeffs));
    PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv));
    if (sinv) PetscCall(PEPEvaluateBasis(pep,pep->target,0,pep->solvematcoeffs,NULL));
    else {
      for (i=0;i<pep->nmat-1;i++) pep->solvematcoeffs[i] = 0.0;
      pep->solvematcoeffs[pep->nmat-1] = 1.0;
    }
  }
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVCreateTensor(pep->V,pep->nmat-1,&ctx->V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Extend the TOAR basis by applying the matrix operator
  over a vector which is decomposed in the TOAR way
  Input:
    - pbc: array containing the polynomial basis coefficients
    - S,V: define the latest Arnoldi vector (nv vectors in V)
  Output:
    - t: new vector extending the TOAR basis
    - r: temporary coefficients to compute the TOAR coefficients
         for the new Arnoldi vector
  Workspace: t_ (two vectors)
*/
static PetscErrorCode PEPTOARExtendBasis(PEP pep,PetscBool sinvert,PetscScalar sigma,PetscScalar *S,PetscInt ls,PetscInt nv,BV V,Vec t,PetscScalar *r,PetscInt lr,Vec *t_)
{
  PetscInt       nmat=pep->nmat,deg=nmat-1,k,j,off=0,lss;
  Vec            v=t_[0],ve=t_[1],q=t_[2];
  PetscScalar    alpha=1.0,*ss,a;
  PetscReal      *ca=pep->pbc,*cb=pep->pbc+nmat,*cg=pep->pbc+2*nmat;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(pep->V,0,nv));
  PetscCall(STGetTransform(pep->st,&flg));
  if (sinvert) {
    for (j=0;j<nv;j++) {
      if (deg>1) r[lr+j] = S[j]/ca[0];
      if (deg>2) r[2*lr+j] = (S[ls+j]+(sigma-cb[1])*r[lr+j])/ca[1];
    }
    for (k=2;k<deg-1;k++) {
      for (j=0;j<nv;j++) r[(k+1)*lr+j] = (S[k*ls+j]+(sigma-cb[k])*r[k*lr+j]-cg[k]*r[(k-1)*lr+j])/ca[k];
    }
    k = deg-1;
    for (j=0;j<nv;j++) r[j] = (S[k*ls+j]+(sigma-cb[k])*r[k*lr+j]-cg[k]*r[(k-1)*lr+j])/ca[k];
    ss = r; lss = lr; off = 1; alpha = -1.0; a = pep->sfactor;
  } else {
    ss = S; lss = ls; off = 0; alpha = -ca[deg-1]; a = 1.0;
  }
  PetscCall(BVMultVec(V,1.0,0.0,v,ss+off*lss));
  if (PetscUnlikely(pep->Dr)) { /* balancing */
    PetscCall(VecPointwiseMult(v,v,pep->Dr));
  }
  PetscCall(STMatMult(pep->st,off,v,q));
  PetscCall(VecScale(q,a));
  for (j=1+off;j<deg+off-1;j++) {
    PetscCall(BVMultVec(V,1.0,0.0,v,ss+j*lss));
    if (PetscUnlikely(pep->Dr)) PetscCall(VecPointwiseMult(v,v,pep->Dr));
    PetscCall(STMatMult(pep->st,j,v,t));
    a *= pep->sfactor;
    PetscCall(VecAXPY(q,a,t));
  }
  if (sinvert) {
    PetscCall(BVMultVec(V,1.0,0.0,v,ss));
    if (PetscUnlikely(pep->Dr)) PetscCall(VecPointwiseMult(v,v,pep->Dr));
    PetscCall(STMatMult(pep->st,deg,v,t));
    a *= pep->sfactor;
    PetscCall(VecAXPY(q,a,t));
  } else {
    PetscCall(BVMultVec(V,1.0,0.0,ve,ss+(deg-1)*lss));
    if (PetscUnlikely(pep->Dr)) PetscCall(VecPointwiseMult(ve,ve,pep->Dr));
    a *= pep->sfactor;
    PetscCall(STMatMult(pep->st,deg-1,ve,t));
    PetscCall(VecAXPY(q,a,t));
    a *= pep->sfactor;
  }
  if (flg || !sinvert) alpha /= a;
  PetscCall(STMatSolve(pep->st,q,t));
  PetscCall(VecScale(t,alpha));
  if (!sinvert) {
    PetscCall(VecAXPY(t,cg[deg-1],v));
    PetscCall(VecAXPY(t,cb[deg-1],ve));
  }
  if (PetscUnlikely(pep->Dr)) PetscCall(VecPointwiseDivide(t,t,pep->Dr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute TOAR coefficients of the blocks of the new Arnoldi vector computed
*/
static PetscErrorCode PEPTOARCoefficients(PEP pep,PetscBool sinvert,PetscScalar sigma,PetscInt nv,PetscScalar *S,PetscInt ls,PetscScalar *r,PetscInt lr,PetscScalar *x)
{
  PetscInt    k,j,nmat=pep->nmat,d=nmat-1;
  PetscReal   *ca=pep->pbc,*cb=pep->pbc+nmat,*cg=pep->pbc+2*nmat;
  PetscScalar t=1.0,tp=0.0,tt;

  PetscFunctionBegin;
  if (sinvert) {
    for (k=1;k<d;k++) {
      tt = t;
      t = ((sigma-cb[k-1])*t-cg[k-1]*tp)/ca[k-1]; /* k-th basis polynomial */
      tp = tt;
      for (j=0;j<=nv;j++) r[k*lr+j] += t*x[j];
    }
  } else {
    for (j=0;j<=nv;j++) r[j] = (cb[0]-sigma)*S[j]+ca[0]*S[ls+j];
    for (k=1;k<d-1;k++) {
      for (j=0;j<=nv;j++) r[k*lr+j] = (cb[k]-sigma)*S[k*ls+j]+ca[k]*S[(k+1)*ls+j]+cg[k]*S[(k-1)*ls+j];
    }
    if (sigma!=0.0) for (j=0;j<=nv;j++) r[(d-1)*lr+j] -= sigma*S[(d-1)*ls+j];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute a run of Arnoldi iterations dim(work)=ld
*/
static PetscErrorCode PEPTOARrun(PEP pep,PetscScalar sigma,Mat A,PetscInt k,PetscInt *M,PetscReal *beta,PetscBool *breakdown,Vec *t_)
{
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  PetscInt       j,m=*M,deg=pep->nmat-1,ld;
  PetscInt       ldh,lds,nqt,l;
  Vec            t;
  PetscReal      norm=0.0;
  PetscBool      flg,sinvert=PETSC_FALSE,lindep;
  PetscScalar    *H,*x,*S;
  Mat            MS;

  PetscFunctionBegin;
  *beta = 0.0;
  PetscCall(MatDenseGetArray(A,&H));
  PetscCall(MatDenseGetLDA(A,&ldh));
  PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
  PetscCall(MatDenseGetArray(MS,&S));
  PetscCall(BVGetSizes(pep->V,NULL,NULL,&ld));
  lds = ld*deg;
  PetscCall(BVGetActiveColumns(pep->V,&l,&nqt));
  PetscCall(STGetTransform(pep->st,&flg));
  if (!flg) {
    /* spectral transformation handled by the solver */
    PetscCall(PetscObjectTypeCompareAny((PetscObject)pep->st,&flg,STSINVERT,STSHIFT,""));
    PetscCheck(flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"ST type not supported for TOAR without transforming matrices");
    PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinvert));
  }
  PetscCall(BVSetActiveColumns(ctx->V,0,m));
  for (j=k;j<m;j++) {
    /* apply operator */
    PetscCall(BVGetColumn(pep->V,nqt,&t));
    PetscCall(PEPTOARExtendBasis(pep,sinvert,sigma,S+j*lds,ld,nqt,pep->V,t,S+(j+1)*lds,ld,t_));
    PetscCall(BVRestoreColumn(pep->V,nqt,&t));

    /* orthogonalize */
    if (sinvert) x = S+(j+1)*lds;
    else x = S+(deg-1)*ld+(j+1)*lds;
    PetscCall(BVOrthogonalizeColumn(pep->V,nqt,x,&norm,&lindep));
    if (!lindep) {
      x[nqt] = norm;
      PetscCall(BVScaleColumn(pep->V,nqt,1.0/norm));
      nqt++;
    }

    PetscCall(PEPTOARCoefficients(pep,sinvert,sigma,nqt-1,S+j*lds,ld,S+(j+1)*lds,ld,x));

    /* level-2 orthogonalization */
    PetscCall(BVOrthogonalizeColumn(ctx->V,j+1,H+j*ldh,&norm,breakdown));
    H[j+1+ldh*j] = norm;
    if (PetscUnlikely(*breakdown)) {
      *M = j+1;
      break;
    }
    PetscCall(BVScaleColumn(ctx->V,j+1,1.0/norm));
    PetscCall(BVSetActiveColumns(pep->V,l,nqt));
  }
  *beta = norm;
  PetscCall(BVSetActiveColumns(ctx->V,0,*M));
  PetscCall(MatDenseRestoreArray(MS,&S));
  PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
  PetscCall(MatDenseRestoreArray(A,&H));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Computes T_j = phi_idx(T). In T_j and T_p are phi_{idx-1}(T)
   and phi_{idx-2}(T) respectively or null if idx=0,1.
   Tp and Tj are input/output arguments
*/
static PetscErrorCode PEPEvaluateBasisM(PEP pep,PetscInt k,PetscScalar *T,PetscInt ldt,PetscInt idx,PetscScalar **Tp,PetscScalar **Tj)
{
  PetscInt       i;
  PetscReal      *ca,*cb,*cg;
  PetscScalar    *pt,g,a;
  PetscBLASInt   k_,ldt_;

  PetscFunctionBegin;
  if (idx==0) {
    PetscCall(PetscArrayzero(*Tj,k*k));
    PetscCall(PetscArrayzero(*Tp,k*k));
    for (i=0;i<k;i++) (*Tj)[i+i*k] = 1.0;
  } else {
    PetscCall(PetscBLASIntCast(ldt,&ldt_));
    PetscCall(PetscBLASIntCast(k,&k_));
    ca = pep->pbc; cb = pep->pbc+pep->nmat; cg = pep->pbc+2*pep->nmat;
    for (i=0;i<k;i++) T[i*ldt+i] -= cb[idx-1];
    a = 1/ca[idx-1];
    g = (idx==1)?0.0:-cg[idx-1]/ca[idx-1];
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&a,T,&ldt_,*Tj,&k_,&g,*Tp,&k_));
    pt = *Tj; *Tj = *Tp; *Tp = pt;
    for (i=0;i<k;i++) T[i*ldt+i] += cb[idx-1];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPExtractInvariantPair(PEP pep,PetscScalar sigma,PetscInt sr,PetscInt k,PetscScalar *S,PetscInt ld,PetscInt deg,Mat HH)
{
  PetscInt       i,j,jj,ldh,lds,ldt,d=pep->nmat-1,idxcpy=0;
  PetscScalar    *H,*At,*Bt,*Hj,*Hp,*T,sone=1.0,g,a,*pM,*work;
  PetscBLASInt   k_,sr_,lds_,ldh_,info,*p,lwork,ldt_;
  PetscBool      transf=PETSC_FALSE,flg;
  PetscReal      norm,maxnrm,*rwork;
  BV             *R,Y;
  Mat            M,*A;

  PetscFunctionBegin;
  if (k==0) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatDenseGetArray(HH,&H));
  PetscCall(MatDenseGetLDA(HH,&ldh));
  lds = deg*ld;
  PetscCall(PetscCalloc6(k,&p,sr*k,&At,k*k,&Bt,k*k,&Hj,k*k,&Hp,sr*k,&work));
  PetscCall(PetscBLASIntCast(sr,&sr_));
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(ldh,&ldh_));
  PetscCall(STGetTransform(pep->st,&flg));
  if (!flg) {
    PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&flg));
    if (flg || sigma!=0.0) transf=PETSC_TRUE;
  }
  if (transf) {
    PetscCall(PetscMalloc1(k*k,&T));
    ldt = k;
    for (i=0;i<k;i++) PetscCall(PetscArraycpy(T+k*i,H+i*ldh,k));
    if (flg) {
      PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
      PetscCallBLAS("LAPACKgetrf",LAPACKgetrf_(&k_,&k_,T,&k_,p,&info));
      SlepcCheckLapackInfo("getrf",info);
      PetscCall(PetscBLASIntCast(sr*k,&lwork));
      PetscCallBLAS("LAPACKgetri",LAPACKgetri_(&k_,T,&k_,p,work,&lwork,&info));
      SlepcCheckLapackInfo("getri",info);
      PetscCall(PetscFPTrapPop());
    }
    if (sigma!=0.0) for (i=0;i<k;i++) T[i+k*i] += sigma;
  } else {
    T = H; ldt = ldh;
  }
  PetscCall(PetscBLASIntCast(ldt,&ldt_));
  switch (pep->extract) {
  case PEP_EXTRACT_NONE:
    break;
  case PEP_EXTRACT_NORM:
    if (pep->basis == PEP_BASIS_MONOMIAL) {
      PetscCall(PetscBLASIntCast(ldt,&ldt_));
      PetscCall(PetscMalloc1(k,&rwork));
      norm = LAPACKlange_("F",&k_,&k_,T,&ldt_,rwork);
      PetscCall(PetscFree(rwork));
      if (norm>1.0) idxcpy = d-1;
    } else {
      PetscCall(PetscBLASIntCast(ldt,&ldt_));
      PetscCall(PetscMalloc1(k,&rwork));
      maxnrm = 0.0;
      for (i=0;i<pep->nmat-1;i++) {
        PetscCall(PEPEvaluateBasisM(pep,k,T,ldt,i,&Hp,&Hj));
        norm = LAPACKlange_("F",&k_,&k_,Hj,&k_,rwork);
        if (norm > maxnrm) {
          idxcpy = i;
          maxnrm = norm;
        }
      }
      PetscCall(PetscFree(rwork));
    }
    if (idxcpy>0) {
      /* copy block idxcpy of S to the first one */
      for (j=0;j<k;j++) PetscCall(PetscArraycpy(S+j*lds,S+idxcpy*ld+j*lds,sr));
    }
    break;
  case PEP_EXTRACT_RESIDUAL:
    PetscCall(STGetTransform(pep->st,&flg));
    if (flg) {
      PetscCall(PetscMalloc1(pep->nmat,&A));
      for (i=0;i<pep->nmat;i++) PetscCall(STGetMatrixTransformed(pep->st,i,A+i));
    } else A = pep->A;
    PetscCall(PetscMalloc1(pep->nmat-1,&R));
    for (i=0;i<pep->nmat-1;i++) PetscCall(BVDuplicateResize(pep->V,k,R+i));
    PetscCall(BVDuplicateResize(pep->V,sr,&Y));
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,sr,k,NULL,&M));
    g = 0.0; a = 1.0;
    PetscCall(BVSetActiveColumns(pep->V,0,sr));
    for (j=0;j<pep->nmat;j++) {
      PetscCall(BVMatMult(pep->V,A[j],Y));
      PetscCall(PEPEvaluateBasisM(pep,k,T,ldt,i,&Hp,&Hj));
      for (i=0;i<pep->nmat-1;i++) {
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&sr_,&k_,&k_,&a,S+i*ld,&lds_,Hj,&k_,&g,At,&sr_));
        PetscCall(MatDenseGetArray(M,&pM));
        for (jj=0;jj<k;jj++) PetscCall(PetscArraycpy(pM+jj*sr,At+jj*sr,sr));
        PetscCall(MatDenseRestoreArray(M,&pM));
        PetscCall(BVMult(R[i],1.0,(i==0)?0.0:1.0,Y,M));
      }
    }

    /* frobenius norm */
    maxnrm = 0.0;
    for (i=0;i<pep->nmat-1;i++) {
      PetscCall(BVNorm(R[i],NORM_FROBENIUS,&norm));
      if (maxnrm > norm) {
        maxnrm = norm;
        idxcpy = i;
      }
    }
    if (idxcpy>0) {
      /* copy block idxcpy of S to the first one */
      for (j=0;j<k;j++) PetscCall(PetscArraycpy(S+j*lds,S+idxcpy*ld+j*lds,sr));
    }
    if (flg) PetscCall(PetscFree(A));
    for (i=0;i<pep->nmat-1;i++) PetscCall(BVDestroy(&R[i]));
    PetscCall(PetscFree(R));
    PetscCall(BVDestroy(&Y));
    PetscCall(MatDestroy(&M));
    break;
  case PEP_EXTRACT_STRUCTURED:
    for (j=0;j<k;j++) Bt[j+j*k] = 1.0;
    for (j=0;j<sr;j++) {
      for (i=0;i<k;i++) At[j*k+i] = PetscConj(S[i*lds+j]);
    }
    PetscCall(PEPEvaluateBasisM(pep,k,T,ldt,0,&Hp,&Hj));
    for (i=1;i<deg;i++) {
      PetscCall(PEPEvaluateBasisM(pep,k,T,ldt,i,&Hp,&Hj));
      PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&k_,&sr_,&k_,&sone,Hj,&k_,S+i*ld,&lds_,&sone,At,&k_));
      PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&k_,&k_,&k_,&sone,Hj,&k_,Hj,&k_,&sone,Bt,&k_));
    }
    PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
    PetscCallBLAS("LAPACKgesv",LAPACKgesv_(&k_,&sr_,Bt,&k_,p,At,&k_,&info));
    PetscCall(PetscFPTrapPop());
    SlepcCheckLapackInfo("gesv",info);
    for (j=0;j<sr;j++) {
      for (i=0;i<k;i++) S[i*lds+j] = PetscConj(At[j*k+i]);
    }
    break;
  }
  if (transf) PetscCall(PetscFree(T));
  PetscCall(PetscFree6(p,At,Bt,Hj,Hp,work));
  PetscCall(MatDenseRestoreArray(HH,&H));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSolve_TOAR(PEP pep)
{
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  PetscInt       i,j,k,l,nv=0,ld,lds,nq=0,nconv=0;
  PetscInt       nmat=pep->nmat,deg=nmat-1;
  PetscScalar    *S,sigma;
  PetscReal      beta;
  PetscBool      breakdown=PETSC_FALSE,flg,falselock=PETSC_FALSE,sinv=PETSC_FALSE;
  Mat            H,MS,MQ;

  PetscFunctionBegin;
  PetscCall(PetscCitationsRegister(citation,&cited));
  if (ctx->lock) {
    /* undocumented option to use a cheaper locking instead of the true locking */
    PetscCall(PetscOptionsGetBool(NULL,NULL,"-pep_toar_falselocking",&falselock,NULL));
  }
  PetscCall(STGetShift(pep->st,&sigma));

  /* update polynomial basis coefficients */
  PetscCall(STGetTransform(pep->st,&flg));
  if (pep->sfactor!=1.0) {
    for (i=0;i<nmat;i++) {
      pep->pbc[nmat+i] /= pep->sfactor;
      pep->pbc[2*nmat+i] /= pep->sfactor*pep->sfactor;
    }
    if (!flg) {
      pep->target /= pep->sfactor;
      PetscCall(RGPushScale(pep->rg,1.0/pep->sfactor));
      PetscCall(STScaleShift(pep->st,1.0/pep->sfactor));
      sigma /= pep->sfactor;
    } else {
      PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv));
      pep->target = sinv?pep->target*pep->sfactor:pep->target/pep->sfactor;
      PetscCall(RGPushScale(pep->rg,sinv?pep->sfactor:1.0/pep->sfactor));
      PetscCall(STScaleShift(pep->st,sinv?pep->sfactor:1.0/pep->sfactor));
    }
  }

  if (flg) sigma = 0.0;

  /* clean projected matrix (including the extra-arrow) */
  PetscCall(DSSetDimensions(pep->ds,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(DSGetMat(pep->ds,DS_MAT_A,&H));
  PetscCall(MatZeroEntries(H));
  PetscCall(DSRestoreMat(pep->ds,DS_MAT_A,&H));

  /* Get the starting Arnoldi vector */
  PetscCall(BVTensorBuildFirstColumn(ctx->V,pep->nini));

  /* restart loop */
  l = 0;
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;

    /* compute an nv-step Lanczos factorization */
    nv = PetscMax(PetscMin(nconv+pep->mpd,pep->ncv),nv);
    PetscCall(DSGetMat(pep->ds,DS_MAT_A,&H));
    PetscCall(PEPTOARrun(pep,sigma,H,pep->nconv+l,&nv,&beta,&breakdown,pep->work));
    PetscCall(DSRestoreMat(pep->ds,DS_MAT_A,&H));
    PetscCall(DSSetDimensions(pep->ds,nv,pep->nconv,pep->nconv+l));
    PetscCall(DSSetState(pep->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(ctx->V,pep->nconv,nv));

    /* solve projected problem */
    PetscCall(DSSolve(pep->ds,pep->eigr,pep->eigi));
    PetscCall(DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(pep->ds));
    PetscCall(DSSynchronize(pep->ds,pep->eigr,pep->eigi));

    /* check convergence */
    PetscCall(PEPKrylovConvergence(pep,PETSC_FALSE,pep->nconv,nv-pep->nconv,beta,&k));
    PetscCall((*pep->stopping)(pep,pep->its,pep->max_it,k,pep->nev,&pep->reason,pep->stoppingctx));

    /* update l */
    if (pep->reason != PEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = (nv==k)?0:PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
      PetscCall(DSGetTruncateSize(pep->ds,k,nv,&l));
      if (!breakdown) {
        /* prepare the Rayleigh quotient for restart */
        PetscCall(DSTruncate(pep->ds,k+l,PETSC_FALSE));
      }
    }
    nconv = k;
    if (!ctx->lock && pep->reason == PEP_CONVERGED_ITERATING && !breakdown) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(pep,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    /* update S */
    PetscCall(DSGetMat(pep->ds,DS_MAT_Q,&MQ));
    PetscCall(BVMultInPlace(ctx->V,MQ,pep->nconv,k+l));
    PetscCall(DSRestoreMat(pep->ds,DS_MAT_Q,&MQ));

    /* copy last column of S */
    PetscCall(BVCopyColumn(ctx->V,nv,k+l));

    if (PetscUnlikely(breakdown && pep->reason == PEP_CONVERGED_ITERATING)) {
      /* stop if breakdown */
      PetscCall(PetscInfo(pep,"Breakdown TOAR method (it=%" PetscInt_FMT " norm=%g)\n",pep->its,(double)beta));
      pep->reason = PEP_DIVERGED_BREAKDOWN;
    }
    if (pep->reason != PEP_CONVERGED_ITERATING) l--;
    /* truncate S */
    PetscCall(BVGetActiveColumns(pep->V,NULL,&nq));
    if (k+l+deg<=nq) {
      PetscCall(BVSetActiveColumns(ctx->V,pep->nconv,k+l+1));
      if (!falselock && ctx->lock) PetscCall(BVTensorCompress(ctx->V,k-pep->nconv));
      else PetscCall(BVTensorCompress(ctx->V,0));
    }
    pep->nconv = k;
    PetscCall(PEPMonitor(pep,pep->its,nconv,pep->eigr,pep->eigi,pep->errest,nv));
  }
  if (pep->nconv>0) {
    /* {V*S_nconv^i}_{i=0}^{d-1} has rank nconv instead of nconv+d-1. Force zeros in each S_nconv^i block */
    PetscCall(BVSetActiveColumns(ctx->V,0,pep->nconv));
    PetscCall(BVGetActiveColumns(pep->V,NULL,&nq));
    PetscCall(BVSetActiveColumns(pep->V,0,nq));
    if (nq>pep->nconv) {
      PetscCall(BVTensorCompress(ctx->V,pep->nconv));
      PetscCall(BVSetActiveColumns(pep->V,0,pep->nconv));
      nq = pep->nconv;
    }

    /* perform Newton refinement if required */
    if (pep->refine==PEP_REFINE_MULTIPLE && pep->rits>0) {
      /* extract invariant pair */
      PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
      PetscCall(MatDenseGetArray(MS,&S));
      PetscCall(DSGetMat(pep->ds,DS_MAT_A,&H));
      PetscCall(BVGetSizes(pep->V,NULL,NULL,&ld));
      lds = deg*ld;
      PetscCall(PEPExtractInvariantPair(pep,sigma,nq,pep->nconv,S,ld,deg,H));
      PetscCall(DSRestoreMat(pep->ds,DS_MAT_A,&H));
      PetscCall(DSSetDimensions(pep->ds,pep->nconv,0,0));
      PetscCall(DSSetState(pep->ds,DS_STATE_RAW));
      PetscCall(PEPNewtonRefinement_TOAR(pep,sigma,&pep->rits,NULL,pep->nconv,S,lds));
      PetscCall(DSSolve(pep->ds,pep->eigr,pep->eigi));
      PetscCall(DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL));
      PetscCall(DSSynchronize(pep->ds,pep->eigr,pep->eigi));
      PetscCall(DSGetMat(pep->ds,DS_MAT_Q,&MQ));
      PetscCall(BVMultInPlace(ctx->V,MQ,0,pep->nconv));
      PetscCall(DSRestoreMat(pep->ds,DS_MAT_Q,&MQ));
      PetscCall(MatDenseRestoreArray(MS,&S));
      PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
    }
  }
  PetscCall(STGetTransform(pep->st,&flg));
  if (pep->refine!=PEP_REFINE_MULTIPLE || pep->rits==0) {
    if (!flg) PetscTryTypeMethod(pep,backtransform);
    if (pep->sfactor!=1.0) {
      for (j=0;j<pep->nconv;j++) {
        pep->eigr[j] *= pep->sfactor;
        pep->eigi[j] *= pep->sfactor;
      }
      /* restore original values */
      for (i=0;i<pep->nmat;i++) {
        pep->pbc[pep->nmat+i] *= pep->sfactor;
        pep->pbc[2*pep->nmat+i] *= pep->sfactor*pep->sfactor;
      }
    }
  }
  /* restore original values */
  if (!flg) {
    pep->target *= pep->sfactor;
    PetscCall(STScaleShift(pep->st,pep->sfactor));
  } else {
    PetscCall(STScaleShift(pep->st,sinv?1.0/pep->sfactor:pep->sfactor));
    pep->target = (sinv)?pep->target/pep->sfactor:pep->target*pep->sfactor;
  }
  if (pep->sfactor!=1.0) PetscCall(RGPopScale(pep->rg));

  /* change the state to raw so that DSVectors() computes eigenvectors from scratch */
  PetscCall(DSSetDimensions(pep->ds,pep->nconv,0,0));
  PetscCall(DSSetState(pep->ds,DS_STATE_RAW));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPTOARSetRestart_TOAR(PEP pep,PetscReal keep)
{
  PEP_TOAR *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  if (keep==(PetscReal)PETSC_DEFAULT || keep==(PetscReal)PETSC_DECIDE) ctx->keep = 0.5;
  else {
    PetscCheck(keep>=0.1 && keep<=0.9,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument must be in the range [0.1,0.9]");
    ctx->keep = keep;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPTOARSetRestart - Sets the restart parameter for the TOAR
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective

   Input Parameters:
+  pep  - the eigenproblem solver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -pep_toar_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: PEPTOARGetRestart()
@*/
PetscErrorCode PEPTOARSetRestart(PEP pep,PetscReal keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,keep,2);
  PetscTryMethod(pep,"PEPTOARSetRestart_C",(PEP,PetscReal),(pep,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPTOARGetRestart_TOAR(PEP pep,PetscReal *keep)
{
  PEP_TOAR *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  *keep = ctx->keep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPTOARGetRestart - Gets the restart parameter used in the TOAR method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: PEPTOARSetRestart()
@*/
PetscErrorCode PEPTOARGetRestart(PEP pep,PetscReal *keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(keep,2);
  PetscUseMethod(pep,"PEPTOARGetRestart_C",(PEP,PetscReal*),(pep,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPTOARSetLocking_TOAR(PEP pep,PetscBool lock)
{
  PEP_TOAR *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPTOARSetLocking - Choose between locking and non-locking variants of
   the TOAR method.

   Logically Collective

   Input Parameters:
+  pep  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -pep_toar_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: PEPTOARGetLocking()
@*/
PetscErrorCode PEPTOARSetLocking(PEP pep,PetscBool lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,lock,2);
  PetscTryMethod(pep,"PEPTOARSetLocking_C",(PEP,PetscBool),(pep,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPTOARGetLocking_TOAR(PEP pep,PetscBool *lock)
{
  PEP_TOAR *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPTOARGetLocking - Gets the locking flag used in the TOAR method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: PEPTOARSetLocking()
@*/
PetscErrorCode PEPTOARGetLocking(PEP pep,PetscBool *lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(lock,2);
  PetscUseMethod(pep,"PEPTOARGetLocking_C",(PEP,PetscBool*),(pep,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetFromOptions_TOAR(PEP pep,PetscOptionItems PetscOptionsObject)
{
  PetscBool      flg,lock;
  PetscReal      keep;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"PEP TOAR Options");

    PetscCall(PetscOptionsReal("-pep_toar_restart","Proportion of vectors kept after restart","PEPTOARSetRestart",0.5,&keep,&flg));
    if (flg) PetscCall(PEPTOARSetRestart(pep,keep));

    PetscCall(PetscOptionsBool("-pep_toar_locking","Choose between locking and non-locking variants","PEPTOARSetLocking",PETSC_FALSE,&lock,&flg));
    if (flg) PetscCall(PEPTOARSetLocking(pep,lock));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPView_TOAR(PEP pep,PetscViewer viewer)
{
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*ctx->keep)));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPDestroy_TOAR(PEP pep)
{
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(PetscFree(pep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPTOARSetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPTOARGetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPTOARSetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPTOARGetLocking_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_TOAR(PEP pep)
{
  PEP_TOAR       *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  pep->data = (void*)ctx;

  pep->lineariz = PETSC_TRUE;
  ctx->lock     = PETSC_TRUE;

  pep->ops->solve          = PEPSolve_TOAR;
  pep->ops->setup          = PEPSetUp_TOAR;
  pep->ops->setfromoptions = PEPSetFromOptions_TOAR;
  pep->ops->destroy        = PEPDestroy_TOAR;
  pep->ops->view           = PEPView_TOAR;
  pep->ops->backtransform  = PEPBackTransform_Default;
  pep->ops->computevectors = PEPComputeVectors_Default;
  pep->ops->extractvectors = PEPExtractVectors_TOAR;

  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPTOARSetRestart_C",PEPTOARSetRestart_TOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPTOARGetRestart_C",PEPTOARGetRestart_TOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPTOARSetLocking_C",PEPTOARSetLocking_TOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPTOARGetLocking_C",PEPTOARGetLocking_TOAR));
  PetscFunctionReturn(PETSC_SUCCESS);
}
