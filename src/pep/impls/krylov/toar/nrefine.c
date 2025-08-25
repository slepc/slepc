/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Newton refinement for polynomial eigenproblems.

   References:

       [1] T. Betcke and D. Kressner, "Perturbation, extraction and refinement
           of invariant pairs for matrix polynomials", Linear Algebra Appl.
           435(3):514-536, 2011.

       [2] C. Campos and J.E. Roman, "Parallel iterative refinement in
           polynomial eigenvalue problems", Numer. Linear Algebra Appl. 23(4):
           730-745, 2016.
*/

#include <slepc/private/pepimpl.h>
#include <slepcblaslapack.h>

typedef struct {
  Mat          *A,M1;
  BV           V,M2,M3,W;
  PetscInt     k,nmat;
  PetscScalar  *fih,*work,*M4;
  PetscBLASInt *pM4;
  PetscBool    compM1;
  Vec          t;
} PEP_REFINE_MATSHELL;

typedef struct {
  Mat          E[2],M1;
  Vec          tN,ttN,t1,vseq;
  VecScatter   scatterctx;
  PetscBool    compM1;
  PetscInt     *map0,*map1,*idxg,*idxp;
  PetscSubcomm subc;
  VecScatter   scatter_sub;
  VecScatter   *scatter_id,*scatterp_id;
  Mat          *A;
  BV           V,W,M2,M3,Wt;
  PetscScalar  *M4,*w,*wt,*d,*dt;
  Vec          t,tg,Rv,Vi,tp,tpg;
  PetscInt     idx,*cols;
} PEP_REFINE_EXPLICIT;

static PetscErrorCode MatMult_FS(Mat M ,Vec x,Vec y)
{
  PEP_REFINE_MATSHELL *ctx;
  PetscInt            k,i;
  PetscScalar         *c;
  PetscBLASInt        k_,one=1,info;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(VecCopy(x,ctx->t));
  k    = ctx->k;
  c    = ctx->work;
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(MatMult(ctx->M1,x,y));
  PetscCall(VecConjugate(ctx->t));
  PetscCall(BVDotVec(ctx->M3,ctx->t,c));
  for (i=0;i<k;i++) c[i] = PetscConj(c[i]);
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&k_,&one,ctx->M4,&k_,ctx->pM4,c,&k_,&info));
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("getrs",info);
  PetscCall(BVMultVec(ctx->M2,-1.0,1.0,y,c));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Evaluates the first d elements of the polynomial basis
  on a given matrix H which is considered to be triangular
*/
static PetscErrorCode PEPEvaluateBasisforMatrix(PEP pep,PetscInt nm,PetscInt k,PetscScalar *H,PetscInt ldh,PetscScalar *fH)
{
  PetscInt       i,j,ldfh=nm*k,off,nmat=pep->nmat;
  PetscReal      *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat,t;
  PetscScalar    corr=0.0,alpha,beta;
  PetscBLASInt   k_,ldh_,ldfh_;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ldh,&ldh_));
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(ldfh,&ldfh_));
  PetscCall(PetscArrayzero(fH,nm*k*k));
  if (nm>0) for (j=0;j<k;j++) fH[j+j*ldfh] = 1.0;
  if (nm>1) {
    t = b[0]/a[0];
    off = k;
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) fH[off+i+j*ldfh] = H[i+j*ldh]/a[0];
      fH[j+j*ldfh] -= t;
    }
  }
  for (i=2;i<nm;i++) {
    off = i*k;
    if (i==2) {
      for (j=0;j<k;j++) {
        fH[off+j+j*ldfh] = 1.0;
        H[j+j*ldh] -= b[1];
      }
    } else {
      for (j=0;j<k;j++) {
        PetscCall(PetscArraycpy(fH+off+j*ldfh,fH+(i-2)*k+j*ldfh,k));
        H[j+j*ldh] += corr-b[i-1];
      }
    }
    corr  = b[i-1];
    beta  = -g[i-1]/a[i-1];
    alpha = 1/a[i-1];
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&alpha,H,&ldh_,fH+(i-1)*k,&ldfh_,&beta,fH+off,&ldfh_));
  }
  for (j=0;j<k;j++) H[j+j*ldh] += corr;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSysSetup_shell(PEP pep,PetscInt k,PetscScalar *fH,PetscScalar *S,PetscInt lds,PetscScalar *fh,PetscScalar h,PEP_REFINE_MATSHELL *ctx)
{
  PetscScalar       *DHii,*T12,*Tr,*Ts,*array,s,ss,sone=1.0,zero=0.0,*M4=ctx->M4,t,*v,*T;
  const PetscScalar *m3,*m2;
  PetscInt          i,d,j,nmat=pep->nmat,lda=nmat*k,deg=nmat-1,nloc,ld2,ld3;
  PetscReal         *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat;
  PetscBLASInt      k_,lda_,lds_,nloc_,ld2_,one=1,info;
  Mat               *A=ctx->A,Mk,M1=ctx->M1,P;
  BV                V=ctx->V,M2=ctx->M2,M3=ctx->M3,W=ctx->W;
  MatStructure      str;
  Vec               vc;

  PetscFunctionBegin;
  PetscCall(STGetMatStructure(pep->st,&str));
  PetscCall(PetscMalloc3(nmat*k*k,&T12,k*k,&Tr,PetscMax(k*k,nmat),&Ts));
  DHii = T12;
  PetscCall(PetscArrayzero(DHii,k*k*nmat));
  for (i=0;i<k;i++) DHii[k+i+i*lda] = 1.0/a[0];
  for (d=2;d<nmat;d++) {
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) {
        DHii[d*k+i+j*lda] = ((h-b[d-1])*DHii[(d-1)*k+i+j*lda]+fH[(d-1)*k+i+j*lda]-g[d-1]*DHii[(d-2)*k+i+j*lda])/(a[d-1]);
      }
    }
  }
  /* T11 */
  if (!ctx->compM1) {
    PetscCall(MatCopy(A[0],M1,DIFFERENT_NONZERO_PATTERN));
    PetscCall(PEPEvaluateBasis(pep,h,0,Ts,NULL));
    for (j=1;j<nmat;j++) PetscCall(MatAXPY(M1,Ts[j],A[j],str));
  }

  /* T22 */
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(lda,&lda_));
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,S,&lds_,S,&lds_,&zero,Tr,&k_));
  for (i=1;i<deg;i++) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,Tr,&k_,DHii+i*k,&lda_,&zero,Ts,&k_));
    s = (i==1)?0.0:1.0;
    PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,fH+i*k,&lda_,Ts,&k_,&s,M4,&k_));
  }
  for (i=0;i<k;i++) for (j=0;j<i;j++) { t=M4[i+j*k];M4[i+j*k]=M4[j+i*k];M4[j+i*k]=t; }

  /* T12 */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&Mk));
  for (i=1;i<nmat;i++) {
    PetscCall(MatDenseGetArrayWrite(Mk,&array));
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,S,&lds_,DHii+i*k,&lda_,&zero,array,&k_));
    PetscCall(MatDenseRestoreArrayWrite(Mk,&array));
    PetscCall(BVSetActiveColumns(W,0,k));
    PetscCall(BVMult(W,1.0,0.0,V,Mk));
    if (i==1) PetscCall(BVMatMult(W,A[i],M2));
    else {
      PetscCall(BVMatMult(W,A[i],M3)); /* using M3 as work space */
      PetscCall(BVMult(M2,1.0,1.0,M3,NULL));
    }
  }

  /* T21 */
  PetscCall(MatDenseGetArrayWrite(Mk,&array));
  for (i=1;i<deg;i++) {
    s = (i==1)?0.0:1.0;
    ss = PetscConj(fh[i]);
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&ss,S,&lds_,fH+i*k,&lda_,&s,array,&k_));
  }
  PetscCall(MatDenseRestoreArrayWrite(Mk,&array));
  PetscCall(BVSetActiveColumns(M3,0,k));
  PetscCall(BVMult(M3,1.0,0.0,V,Mk));
  for (i=0;i<k;i++) {
    PetscCall(BVGetColumn(M3,i,&vc));
    PetscCall(VecConjugate(vc));
    PetscCall(BVRestoreColumn(M3,i,&vc));
  }
  PetscCall(MatDestroy(&Mk));
  PetscCall(PetscFree3(T12,Tr,Ts));

  PetscCall(VecGetLocalSize(ctx->t,&nloc));
  PetscCall(PetscBLASIntCast(nloc,&nloc_));
  PetscCall(PetscMalloc1(nloc*k,&T));
  PetscCall(KSPGetOperators(pep->refineksp,NULL,&P));
  if (!ctx->compM1) PetscCall(MatCopy(ctx->M1,P,SAME_NONZERO_PATTERN));
  PetscCall(BVGetArrayRead(ctx->M2,&m2));
  PetscCall(BVGetLeadingDimension(ctx->M2,&ld2));
  PetscCall(PetscBLASIntCast(ld2,&ld2_));
  PetscCall(BVGetArrayRead(ctx->M3,&m3));
  PetscCall(BVGetLeadingDimension(ctx->M3,&ld3));
  PetscCall(VecGetArray(ctx->t,&v));
  for (i=0;i<nloc;i++) for (j=0;j<k;j++) T[j+i*k] = m3[i+j*ld3];
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKgesv",LAPACKgesv_(&k_,&nloc_,ctx->M4,&k_,ctx->pM4,T,&k_,&info));
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("gesv",info);
  for (i=0;i<nloc;i++) v[i] = BLASdot_(&k_,m2+i,&ld2_,T+i*k,&one);
  PetscCall(VecRestoreArray(ctx->t,&v));
  PetscCall(BVRestoreArrayRead(ctx->M2,&m2));
  PetscCall(BVRestoreArrayRead(ctx->M3,&m3));
  PetscCall(MatDiagonalSet(P,ctx->t,ADD_VALUES));
  PetscCall(PetscFree(T));
  PetscCall(KSPSetUp(pep->refineksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSysSolve_shell(KSP ksp,PetscInt nmat,Vec Rv,PetscScalar *Rh,PetscInt k,Vec dVi,PetscScalar *dHi)
{
  PetscScalar         *t0;
  PetscBLASInt        k_,one=1,info,lda_;
  PetscInt            i,lda=nmat*k;
  Mat                 M;
  PEP_REFINE_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(KSPGetOperators(ksp,&M,NULL));
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(PetscCalloc1(k,&t0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCall(PetscBLASIntCast(lda,&lda_));
  PetscCall(PetscBLASIntCast(k,&k_));
  for (i=0;i<k;i++) t0[i] = Rh[i];
  PetscCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&k_,&one,ctx->M4,&k_,ctx->pM4,t0,&k_,&info));
  SlepcCheckLapackInfo("getrs",info);
  PetscCall(BVMultVec(ctx->M2,-1.0,1.0,Rv,t0));
  PetscCall(KSPSolve(ksp,Rv,dVi));
  PetscCall(VecConjugate(dVi));
  PetscCall(BVDotVec(ctx->M3,dVi,dHi));
  PetscCall(VecConjugate(dVi));
  for (i=0;i<k;i++) dHi[i] = Rh[i]-PetscConj(dHi[i]);
  PetscCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&k_,&one,ctx->M4,&k_,ctx->pM4,dHi,&k_,&info));
  SlepcCheckLapackInfo("getrs",info);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscFree(t0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Computes the residual P(H,V*S)*e_j for the polynomial
*/
static PetscErrorCode NRefRightSide(PetscInt nmat,PetscReal *pcf,Mat *A,PetscInt k,BV V,PetscScalar *S,PetscInt lds,PetscInt j,PetscScalar *H,PetscInt ldh,PetscScalar *fH,PetscScalar *DfH,PetscScalar *dH,BV dV,PetscScalar *dVS,PetscInt rds,Vec Rv,PetscScalar *Rh,BV W,Vec t)
{
  PetscScalar    *DS0,*DS1,*F,beta=0.0,sone=1.0,none=-1.0,tt=0.0,*h,zero=0.0,*Z,*c0;
  PetscReal      *a=pcf,*b=pcf+nmat,*g=b+nmat;
  PetscInt       i,ii,jj,lda;
  PetscBLASInt   lda_,k_,ldh_,lds_,nmat_,k2_,krds_,j_,one=1;
  Mat            M0;
  Vec            w;

  PetscFunctionBegin;
  PetscCall(PetscMalloc4(k*nmat,&h,k*k,&DS0,k*k,&DS1,k*k,&Z));
  lda = k*nmat;
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(lda,&lda_));
  PetscCall(PetscBLASIntCast(nmat,&nmat_));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&nmat_,&k_,&sone,S,&lds_,fH+j*lda,&k_,&zero,h,&k_));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,nmat,h,&M0));
  PetscCall(BVSetActiveColumns(W,0,nmat));
  PetscCall(BVMult(W,1.0,0.0,V,M0));
  PetscCall(MatDestroy(&M0));

  PetscCall(BVGetColumn(W,0,&w));
  PetscCall(MatMult(A[0],w,Rv));
  PetscCall(BVRestoreColumn(W,0,&w));
  for (i=1;i<nmat;i++) {
    PetscCall(BVGetColumn(W,i,&w));
    PetscCall(MatMult(A[i],w,t));
    PetscCall(BVRestoreColumn(W,i,&w));
    PetscCall(VecAXPY(Rv,1.0,t));
  }
  /* Update right-hand side */
  if (j) {
    PetscCall(PetscBLASIntCast(ldh,&ldh_));
    PetscCall(PetscArrayzero(Z,k*k));
    PetscCall(PetscArrayzero(DS0,k*k));
    PetscCall(PetscArraycpy(Z+(j-1)*k,dH+(j-1)*k,k));
    /* Update DfH */
    for (i=1;i<nmat;i++) {
      if (i>1) {
        beta = -g[i-1];
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,fH+(i-1)*k,&lda_,Z,&k_,&beta,DS0,&k_));
        tt += -b[i-1];
        for (ii=0;ii<k;ii++) H[ii+ii*ldh] += tt;
        tt = b[i-1];
        beta = 1.0/a[i-1];
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&beta,DS1,&k_,H,&ldh_,&beta,DS0,&k_));
        F = DS0; DS0 = DS1; DS1 = F;
      } else {
        PetscCall(PetscArrayzero(DS1,k*k));
        for (ii=0;ii<k;ii++) DS1[ii+(j-1)*k] = Z[ii+(j-1)*k]/a[0];
      }
      for (jj=j;jj<k;jj++) {
        for (ii=0;ii<k;ii++) DfH[k*i+ii+jj*lda] += DS1[ii+jj*k];
      }
    }
    for (ii=0;ii<k;ii++) H[ii+ii*ldh] += tt;
    /* Update right-hand side */
    PetscCall(PetscBLASIntCast(2*k,&k2_));
    PetscCall(PetscBLASIntCast(j,&j_));
    PetscCall(PetscBLASIntCast(k+rds,&krds_));
    c0 = DS0;
    PetscCall(PetscArrayzero(Rh,k));
    for (i=0;i<nmat;i++) {
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&krds_,&j_,&sone,dVS,&k2_,fH+j*lda+i*k,&one,&zero,h,&one));
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&k_,&k_,&sone,S,&lds_,DfH+i*k+j*lda,&one,&sone,h,&one));
      PetscCall(BVMultVec(V,1.0,0.0,t,h));
      PetscCall(BVSetActiveColumns(dV,0,rds));
      PetscCall(BVMultVec(dV,1.0,1.0,t,h+k));
      PetscCall(BVGetColumn(W,i,&w));
      PetscCall(MatMult(A[i],t,w));
      PetscCall(BVRestoreColumn(W,i,&w));
      if (i>0 && i<nmat-1) {
        PetscCallBLAS("BLASgemv",BLASgemv_("C",&k_,&k_,&sone,S,&lds_,h,&one,&zero,c0,&one));
        PetscCallBLAS("BLASgemv",BLASgemv_("C",&k_,&k_,&none,fH+i*k,&lda_,c0,&one,&sone,Rh,&one));
      }
    }

    for (i=0;i<nmat;i++) h[i] = -1.0;
    PetscCall(BVMultVec(W,1.0,1.0,Rv,h));
  }
  PetscCall(PetscFree4(h,DS0,DS1,Z));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSysSolve_mbe(PetscInt k,PetscInt sz,BV W,PetscScalar *w,BV Wt,PetscScalar *wt,PetscScalar *d,PetscScalar *dt,KSP ksp,BV T2,BV T3 ,PetscScalar *T4,PetscBool trans,Vec x1,PetscScalar *x2,Vec sol1,PetscScalar *sol2,Vec vw)
{
  PetscInt       i,j,incf,incc;
  PetscScalar    *y,*g,*xx2,*ww,y2,*dd;
  Vec            v,t,xx1;
  BV             WW,T;

  PetscFunctionBegin;
  PetscCall(PetscMalloc3(sz,&y,sz,&g,k,&xx2));
  if (trans) {
    WW = W; ww = w; dd = d; T = T3; incf = 0; incc = 1;
  } else {
    WW = Wt; ww = wt; dd = dt; T = T2; incf = 1; incc = 0;
  }
  xx1 = vw;
  PetscCall(VecCopy(x1,xx1));
  PetscCall(PetscArraycpy(xx2,x2,sz));
  PetscCall(PetscArrayzero(sol2,k));
  for (i=sz-1;i>=0;i--) {
    PetscCall(BVGetColumn(WW,i,&v));
    PetscCall(VecConjugate(v));
    PetscCall(VecDot(xx1,v,y+i));
    PetscCall(VecConjugate(v));
    PetscCall(BVRestoreColumn(WW,i,&v));
    for (j=0;j<i;j++) y[i] += ww[j+i*k]*xx2[j];
    y[i] = -(y[i]-xx2[i])/dd[i];
    PetscCall(BVGetColumn(T,i,&t));
    PetscCall(VecAXPY(xx1,-y[i],t));
    PetscCall(BVRestoreColumn(T,i,&t));
    for (j=0;j<=i;j++) xx2[j] -= y[i]*T4[j*incf+incc*i+(i*incf+incc*j)*k];
    g[i] = xx2[i];
  }
  if (trans) PetscCall(KSPSolveTranspose(ksp,xx1,sol1));
  else PetscCall(KSPSolve(ksp,xx1,sol1));
  if (trans) {
    WW = Wt; ww = wt; dd = dt; T = T2; incf = 1; incc = 0;
  } else {
    WW = W; ww = w; dd = d; T = T3; incf = 0; incc = 1;
  }
  for (i=0;i<sz;i++) {
    PetscCall(BVGetColumn(T,i,&t));
    PetscCall(VecConjugate(t));
    PetscCall(VecDot(sol1,t,&y2));
    PetscCall(VecConjugate(t));
    PetscCall(BVRestoreColumn(T,i,&t));
    for (j=0;j<i;j++) y2 += sol2[j]*T4[j*incf+incc*i+(i*incf+incc*j)*k];
    y2 = (g[i]-y2)/dd[i];
    PetscCall(BVGetColumn(WW,i,&v));
    PetscCall(VecAXPY(sol1,-y2,v));
    for (j=0;j<i;j++) sol2[j] -= ww[j+i*k]*y2;
    sol2[i] = y[i]+y2;
    PetscCall(BVRestoreColumn(WW,i,&v));
  }
  PetscCall(PetscFree3(y,g,xx2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSysSetup_mbe(PEP pep,PetscInt k,KSP ksp,PetscScalar *fH,PetscScalar *S,PetscInt lds,PetscScalar *fh,PetscScalar h,BV V,PEP_REFINE_EXPLICIT *matctx)
{
  PetscInt       i,j,l,nmat=pep->nmat,lda=nmat*k,deg=nmat-1;
  Mat            M1=matctx->M1,*A,*At,Mk;
  PetscReal      *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat;
  PetscScalar    s,ss,*DHii,*T12,*array,*Ts,*Tr,*M4=matctx->M4,sone=1.0,zero=0.0;
  PetscScalar    *w=matctx->w,*wt=matctx->wt,*d=matctx->d,*dt=matctx->dt;
  PetscBLASInt   lds_,lda_,k_;
  MatStructure   str;
  PetscBool      flg;
  BV             M2=matctx->M2,M3=matctx->M3,W=matctx->W,Wt=matctx->Wt;
  Vec            vc,vc2;

  PetscFunctionBegin;
  PetscCall(PetscMalloc3(nmat*k*k,&T12,k*k,&Tr,PetscMax(k*k,nmat),&Ts));
  PetscCall(STGetMatStructure(pep->st,&str));
  PetscCall(STGetTransform(pep->st,&flg));
  if (flg) {
    PetscCall(PetscMalloc1(pep->nmat,&At));
    for (i=0;i<pep->nmat;i++) PetscCall(STGetMatrixTransformed(pep->st,i,&At[i]));
  } else At = pep->A;
  if (matctx->subc) A = matctx->A;
  else A = At;
  /* Form the explicit system matrix */
  DHii = T12;
  PetscCall(PetscArrayzero(DHii,k*k*nmat));
  for (i=0;i<k;i++) DHii[k+i+i*lda] = 1.0/a[0];
  for (l=2;l<nmat;l++) {
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) {
        DHii[l*k+i+j*lda] = ((h-b[l-1])*DHii[(l-1)*k+i+j*lda]+fH[(l-1)*k+i+j*lda]-g[l-1]*DHii[(l-2)*k+i+j*lda])/a[l-1];
      }
    }
  }

  /* T11 */
  if (!matctx->compM1) {
    PetscCall(MatCopy(A[0],M1,DIFFERENT_NONZERO_PATTERN));
    PetscCall(PEPEvaluateBasis(pep,h,0,Ts,NULL));
    for (j=1;j<nmat;j++) PetscCall(MatAXPY(M1,Ts[j],A[j],str));
  }

  /* T22 */
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(lda,&lda_));
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,S,&lds_,S,&lds_,&zero,Tr,&k_));
  for (i=1;i<deg;i++) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,Tr,&k_,DHii+i*k,&lda_,&zero,Ts,&k_));
    s = (i==1)?0.0:1.0;
    PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,fH+i*k,&lda_,Ts,&k_,&s,M4,&k_));
  }

  /* T12 */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&Mk));
  for (i=1;i<nmat;i++) {
    PetscCall(MatDenseGetArrayWrite(Mk,&array));
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,S,&lds_,DHii+i*k,&lda_,&zero,array,&k_));
    PetscCall(MatDenseRestoreArrayWrite(Mk,&array));
    PetscCall(BVSetActiveColumns(W,0,k));
    PetscCall(BVMult(W,1.0,0.0,V,Mk));
    if (i==1) PetscCall(BVMatMult(W,A[i],M2));
    else {
      PetscCall(BVMatMult(W,A[i],M3)); /* using M3 as work space */
      PetscCall(BVMult(M2,1.0,1.0,M3,NULL));
    }
  }

  /* T21 */
  PetscCall(MatDenseGetArrayWrite(Mk,&array));
  for (i=1;i<deg;i++) {
    s = (i==1)?0.0:1.0;
    ss = PetscConj(fh[i]);
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&ss,S,&lds_,fH+i*k,&lda_,&s,array,&k_));
  }
  PetscCall(MatDenseRestoreArrayWrite(Mk,&array));
  PetscCall(BVSetActiveColumns(M3,0,k));
  PetscCall(BVMult(M3,1.0,0.0,V,Mk));
  for (i=0;i<k;i++) {
    PetscCall(BVGetColumn(M3,i,&vc));
    PetscCall(VecConjugate(vc));
    PetscCall(BVRestoreColumn(M3,i,&vc));
  }

  PetscCall(PEP_KSPSetOperators(ksp,M1,M1));
  PetscCall(KSPSetUp(ksp));
  PetscCall(MatDestroy(&Mk));

  /* Set up for BEMW */
  for (i=0;i<k;i++) {
    PetscCall(BVGetColumn(M2,i,&vc));
    PetscCall(BVGetColumn(W,i,&vc2));
    PetscCall(NRefSysSolve_mbe(k,i,W,w,Wt,wt,d,dt,ksp,M2,M3,M4,PETSC_FALSE,vc,M4+i*k,vc2,w+i*k,matctx->t));
    PetscCall(BVRestoreColumn(M2,i,&vc));
    PetscCall(BVGetColumn(M3,i,&vc));
    PetscCall(VecConjugate(vc));
    PetscCall(VecDot(vc2,vc,&d[i]));
    PetscCall(VecConjugate(vc));
    PetscCall(BVRestoreColumn(M3,i,&vc));
    for (j=0;j<i;j++) d[i] += M4[i+j*k]*w[j+i*k];
    d[i] = M4[i+i*k]-d[i];
    PetscCall(BVRestoreColumn(W,i,&vc2));

    PetscCall(BVGetColumn(M3,i,&vc));
    PetscCall(BVGetColumn(Wt,i,&vc2));
    for (j=0;j<=i;j++) Ts[j] = M4[i+j*k];
    PetscCall(NRefSysSolve_mbe(k,i,W,w,Wt,wt,d,dt,ksp,M2,M3,M4,PETSC_TRUE,vc,Ts,vc2,wt+i*k,matctx->t));
    PetscCall(BVRestoreColumn(M3,i,&vc));
    PetscCall(BVGetColumn(M2,i,&vc));
    PetscCall(VecConjugate(vc2));
    PetscCall(VecDot(vc,vc2,&dt[i]));
    PetscCall(VecConjugate(vc2));
    PetscCall(BVRestoreColumn(M2,i,&vc));
    for (j=0;j<i;j++) dt[i] += M4[j+i*k]*wt[j+i*k];
    dt[i] = M4[i+i*k]-dt[i];
    PetscCall(BVRestoreColumn(Wt,i,&vc2));
  }

  if (flg) PetscCall(PetscFree(At));
  PetscCall(PetscFree3(T12,Tr,Ts));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSysSetup_explicit(PEP pep,PetscInt k,KSP ksp,PetscScalar *fH,PetscScalar *S,PetscInt lds,PetscScalar *fh,PetscScalar h,BV V,PEP_REFINE_EXPLICIT *matctx,BV W)
{
  PetscInt          i,j,d,n,n0,m0,n1,m1,nmat=pep->nmat,lda=nmat*k,deg=nmat-1;
  PetscInt          *idxg=matctx->idxg,*idxp=matctx->idxp,idx,ncols;
  Mat               M,*E=matctx->E,*A,*At,Mk,Md;
  PetscReal         *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat;
  PetscScalar       s,ss,*DHii,*T22,*T21,*T12,*Ts,*Tr,*array,*ts,sone=1.0,zero=0.0;
  PetscBLASInt      lds_,lda_,k_;
  const PetscInt    *idxmc;
  const PetscScalar *valsc,*carray;
  MatStructure      str;
  Vec               vc,vc0;
  PetscBool         flg;

  PetscFunctionBegin;
  PetscCall(PetscMalloc5(k*k,&T22,k*k,&T21,nmat*k*k,&T12,k*k,&Tr,k*k,&Ts));
  PetscCall(STGetMatStructure(pep->st,&str));
  PetscCall(KSPGetOperators(ksp,&M,NULL));
  PetscCall(MatGetOwnershipRange(E[1],&n1,&m1));
  PetscCall(MatGetOwnershipRange(E[0],&n0,&m0));
  PetscCall(MatGetOwnershipRange(M,&n,NULL));
  PetscCall(PetscMalloc1(nmat,&ts));
  PetscCall(STGetTransform(pep->st,&flg));
  if (flg) {
    PetscCall(PetscMalloc1(pep->nmat,&At));
    for (i=0;i<pep->nmat;i++) PetscCall(STGetMatrixTransformed(pep->st,i,&At[i]));
  } else At = pep->A;
  if (matctx->subc) A = matctx->A;
  else A = At;
  /* Form the explicit system matrix */
  DHii = T12;
  PetscCall(PetscArrayzero(DHii,k*k*nmat));
  for (i=0;i<k;i++) DHii[k+i+i*lda] = 1.0/a[0];
  for (d=2;d<nmat;d++) {
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) {
        DHii[d*k+i+j*lda] = ((h-b[d-1])*DHii[(d-1)*k+i+j*lda]+fH[(d-1)*k+i+j*lda]-g[d-1]*DHii[(d-2)*k+i+j*lda])/a[d-1];
      }
    }
  }

  /* T11 */
  if (!matctx->compM1) {
    PetscCall(MatCopy(A[0],E[0],DIFFERENT_NONZERO_PATTERN));
    PetscCall(PEPEvaluateBasis(pep,h,0,Ts,NULL));
    for (j=1;j<nmat;j++) PetscCall(MatAXPY(E[0],Ts[j],A[j],str));
  }
  for (i=n0;i<m0;i++) {
    PetscCall(MatGetRow(E[0],i,&ncols,&idxmc,&valsc));
    idx = n+i-n0;
    for (j=0;j<ncols;j++) {
      idxg[j] = matctx->map0[idxmc[j]];
    }
    PetscCall(MatSetValues(M,1,&idx,ncols,idxg,valsc,INSERT_VALUES));
    PetscCall(MatRestoreRow(E[0],i,&ncols,&idxmc,&valsc));
  }

  /* T22 */
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(lda,&lda_));
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,S,&lds_,S,&lds_,&zero,Tr,&k_));
  for (i=1;i<deg;i++) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,Tr,&k_,DHii+i*k,&lda_,&zero,Ts,&k_));
    s = (i==1)?0.0:1.0;
    PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,fH+i*k,&lda_,Ts,&k_,&s,T22,&k_));
  }
  for (j=0;j<k;j++) idxp[j] = matctx->map1[j];
  for (i=0;i<m1-n1;i++) {
    idx = n+m0-n0+i;
    for (j=0;j<k;j++) {
      Tr[j] = T22[n1+i+j*k];
    }
    PetscCall(MatSetValues(M,1,&idx,k,idxp,Tr,INSERT_VALUES));
  }

  /* T21 */
  for (i=1;i<deg;i++) {
    s = (i==1)?0.0:1.0;
    ss = PetscConj(fh[i]);
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&ss,S,&lds_,fH+i*k,&lda_,&s,T21,&k_));
  }
  PetscCall(BVSetActiveColumns(W,0,k));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,T21,&Mk));
  PetscCall(BVMult(W,1.0,0.0,V,Mk));
  for (i=0;i<k;i++) {
    PetscCall(BVGetColumn(W,i,&vc));
    PetscCall(VecConjugate(vc));
    PetscCall(VecGetArrayRead(vc,&carray));
    idx = matctx->map1[i];
    PetscCall(MatSetValues(M,1,&idx,m0-n0,matctx->map0+n0,carray,INSERT_VALUES));
    PetscCall(VecRestoreArrayRead(vc,&carray));
    PetscCall(BVRestoreColumn(W,i,&vc));
  }

  /* T12 */
  for (i=1;i<nmat;i++) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,S,&lds_,DHii+i*k,&lda_,&zero,Ts,&k_));
    for (j=0;j<k;j++) PetscCall(PetscArraycpy(T12+i*k+j*lda,Ts+j*k,k));
  }
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,nmat-1,NULL,&Md));
  for (i=0;i<nmat;i++) ts[i] = 1.0;
  for (j=0;j<k;j++) {
    PetscCall(MatDenseGetArrayWrite(Md,&array));
    PetscCall(PetscArraycpy(array,T12+k+j*lda,(nmat-1)*k));
    PetscCall(MatDenseRestoreArrayWrite(Md,&array));
    PetscCall(BVSetActiveColumns(W,0,nmat-1));
    PetscCall(BVMult(W,1.0,0.0,V,Md));
    for (i=nmat-1;i>0;i--) {
      PetscCall(BVGetColumn(W,i-1,&vc0));
      PetscCall(BVGetColumn(W,i,&vc));
      PetscCall(MatMult(A[i],vc0,vc));
      PetscCall(BVRestoreColumn(W,i-1,&vc0));
      PetscCall(BVRestoreColumn(W,i,&vc));
    }
    PetscCall(BVSetActiveColumns(W,1,nmat));
    PetscCall(BVGetColumn(W,0,&vc0));
    PetscCall(BVMultVec(W,1.0,0.0,vc0,ts));
    PetscCall(VecGetArrayRead(vc0,&carray));
    idx = matctx->map1[j];
    PetscCall(MatSetValues(M,m0-n0,matctx->map0+n0,1,&idx,carray,INSERT_VALUES));
    PetscCall(VecRestoreArrayRead(vc0,&carray));
    PetscCall(BVRestoreColumn(W,0,&vc0));
  }
  PetscCall(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY));
  PetscCall(PEP_KSPSetOperators(ksp,M,M));
  PetscCall(KSPSetUp(ksp));
  PetscCall(PetscFree(ts));
  PetscCall(MatDestroy(&Mk));
  PetscCall(MatDestroy(&Md));
  if (flg) PetscCall(PetscFree(At));
  PetscCall(PetscFree5(T22,T21,T12,Tr,Ts));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSysSolve_explicit(PetscInt k,KSP ksp,Vec Rv,PetscScalar *Rh,Vec dVi,PetscScalar *dHi,PEP_REFINE_EXPLICIT *matctx)
{
  PetscInt          n0,m0,n1,m1,i;
  PetscScalar       *arrayV;
  const PetscScalar *array;

  PetscFunctionBegin;
  PetscCall(MatGetOwnershipRange(matctx->E[1],&n1,&m1));
  PetscCall(MatGetOwnershipRange(matctx->E[0],&n0,&m0));

  /* Right side */
  PetscCall(VecGetArrayRead(Rv,&array));
  PetscCall(VecSetValues(matctx->tN,m0-n0,matctx->map0+n0,array,INSERT_VALUES));
  PetscCall(VecRestoreArrayRead(Rv,&array));
  PetscCall(VecSetValues(matctx->tN,m1-n1,matctx->map1+n1,Rh+n1,INSERT_VALUES));
  PetscCall(VecAssemblyBegin(matctx->tN));
  PetscCall(VecAssemblyEnd(matctx->tN));

  /* Solve */
  PetscCall(KSPSolve(ksp,matctx->tN,matctx->ttN));

  /* Retrieve solution */
  PetscCall(VecGetArray(dVi,&arrayV));
  PetscCall(VecGetArrayRead(matctx->ttN,&array));
  PetscCall(PetscArraycpy(arrayV,array,m0-n0));
  PetscCall(VecRestoreArray(dVi,&arrayV));
  if (!matctx->subc) {
    PetscCall(VecGetArray(matctx->t1,&arrayV));
    for (i=0;i<m1-n1;i++) arrayV[i] =  array[m0-n0+i];
    PetscCall(VecRestoreArray(matctx->t1,&arrayV));
    PetscCall(VecRestoreArrayRead(matctx->ttN,&array));
    PetscCall(VecScatterBegin(matctx->scatterctx,matctx->t1,matctx->vseq,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(matctx->scatterctx,matctx->t1,matctx->vseq,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGetArrayRead(matctx->vseq,&array));
    for (i=0;i<k;i++) dHi[i] = array[i];
    PetscCall(VecRestoreArrayRead(matctx->vseq,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSysIter(PetscInt i,PEP pep,PetscInt k,KSP ksp,PetscScalar *fH,PetscScalar *S,PetscInt lds,PetscScalar *fh,PetscScalar *H,PetscInt ldh,Vec Rv,PetscScalar *Rh,BV V,Vec dVi,PetscScalar *dHi,PEP_REFINE_EXPLICIT *matctx,BV W)
{
  PetscInt            j,m,lda=pep->nmat*k,n0,m0,idx;
  PetscMPIInt         root,len;
  PetscScalar         *array2,h;
  const PetscScalar   *array;
  Vec                 R,Vi;
  PEP_REFINE_MATSHELL *ctx;
  Mat                 M;

  PetscFunctionBegin;
  if (!matctx || !matctx->subc) {
    for (j=0;j<pep->nmat;j++) fh[j] = fH[j*k+i+i*lda];
    h   = H[i+i*ldh];
    idx = i;
    R   = Rv;
    Vi  = dVi;
    switch (pep->scheme) {
    case PEP_REFINE_SCHEME_EXPLICIT:
      PetscCall(NRefSysSetup_explicit(pep,k,ksp,fH,S,lds,fh,h,V,matctx,W));
      matctx->compM1 = PETSC_FALSE;
      break;
    case PEP_REFINE_SCHEME_MBE:
      PetscCall(NRefSysSetup_mbe(pep,k,ksp,fH,S,lds,fh,h,V,matctx));
      matctx->compM1 = PETSC_FALSE;
      break;
    case PEP_REFINE_SCHEME_SCHUR:
      PetscCall(KSPGetOperators(ksp,&M,NULL));
      PetscCall(MatShellGetContext(M,&ctx));
      PetscCall(NRefSysSetup_shell(pep,k,fH,S,lds,fh,h,ctx));
      ctx->compM1 = PETSC_FALSE;
      break;
    }
  } else {
    if (i%matctx->subc->n==0 && (idx=i+matctx->subc->color)<k) {
      for (j=0;j<pep->nmat;j++) fh[j] = fH[j*k+idx+idx*lda];
      h = H[idx+idx*ldh];
      matctx->idx = idx;
      switch (pep->scheme) {
      case PEP_REFINE_SCHEME_EXPLICIT:
        PetscCall(NRefSysSetup_explicit(pep,k,ksp,fH,S,lds,fh,h,matctx->V,matctx,matctx->W));
        matctx->compM1 = PETSC_FALSE;
        break;
      case PEP_REFINE_SCHEME_MBE:
        PetscCall(NRefSysSetup_mbe(pep,k,ksp,fH,S,lds,fh,h,matctx->V,matctx));
        matctx->compM1 = PETSC_FALSE;
        break;
      case PEP_REFINE_SCHEME_SCHUR:
        break;
      }
    } else idx = matctx->idx;
    PetscCall(VecScatterBegin(matctx->scatter_id[i%matctx->subc->n],Rv,matctx->tg,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(matctx->scatter_id[i%matctx->subc->n],Rv,matctx->tg,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGetArrayRead(matctx->tg,&array));
    PetscCall(VecPlaceArray(matctx->t,array));
    PetscCall(VecCopy(matctx->t,matctx->Rv));
    PetscCall(VecResetArray(matctx->t));
    PetscCall(VecRestoreArrayRead(matctx->tg,&array));
    R  = matctx->Rv;
    Vi = matctx->Vi;
  }
  if (idx==i && idx<k) {
    switch (pep->scheme) {
      case PEP_REFINE_SCHEME_EXPLICIT:
        PetscCall(NRefSysSolve_explicit(k,ksp,R,Rh,Vi,dHi,matctx));
        break;
      case PEP_REFINE_SCHEME_MBE:
        PetscCall(NRefSysSolve_mbe(k,k,matctx->W,matctx->w,matctx->Wt,matctx->wt,matctx->d,matctx->dt,ksp,matctx->M2,matctx->M3 ,matctx->M4,PETSC_FALSE,R,Rh,Vi,dHi,matctx->t));
        break;
      case PEP_REFINE_SCHEME_SCHUR:
        PetscCall(NRefSysSolve_shell(ksp,pep->nmat,R,Rh,k,Vi,dHi));
        break;
    }
  }
  if (matctx && matctx->subc) {
    PetscCall(VecGetLocalSize(Vi,&m));
    PetscCall(VecGetArrayRead(Vi,&array));
    PetscCall(VecGetArray(matctx->tg,&array2));
    PetscCall(PetscArraycpy(array2,array,m));
    PetscCall(VecRestoreArray(matctx->tg,&array2));
    PetscCall(VecRestoreArrayRead(Vi,&array));
    PetscCall(VecScatterBegin(matctx->scatter_id[i%matctx->subc->n],matctx->tg,dVi,INSERT_VALUES,SCATTER_REVERSE));
    PetscCall(VecScatterEnd(matctx->scatter_id[i%matctx->subc->n],matctx->tg,dVi,INSERT_VALUES,SCATTER_REVERSE));
    switch (pep->scheme) {
    case PEP_REFINE_SCHEME_EXPLICIT:
      PetscCall(MatGetOwnershipRange(matctx->E[0],&n0,&m0));
      PetscCall(VecGetArrayRead(matctx->ttN,&array));
      PetscCall(VecPlaceArray(matctx->tp,array+m0-n0));
      PetscCall(VecScatterBegin(matctx->scatterp_id[i%matctx->subc->n],matctx->tp,matctx->tpg,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecScatterEnd(matctx->scatterp_id[i%matctx->subc->n],matctx->tp,matctx->tpg,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecResetArray(matctx->tp));
      PetscCall(VecRestoreArrayRead(matctx->ttN,&array));
      PetscCall(VecGetArrayRead(matctx->tpg,&array));
      for (j=0;j<k;j++) dHi[j] = array[j];
      PetscCall(VecRestoreArrayRead(matctx->tpg,&array));
      break;
     case PEP_REFINE_SCHEME_MBE:
      root = 0;
      for (j=0;j<i%matctx->subc->n;j++) root += matctx->subc->subsize[j];
      PetscCall(PetscMPIIntCast(k,&len));
      PetscCallMPI(MPI_Bcast(dHi,len,MPIU_SCALAR,root,matctx->subc->dupparent));
      break;
    case PEP_REFINE_SCHEME_SCHUR:
      break;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPNRefForwardSubstitution(PEP pep,PetscInt k,PetscScalar *S,PetscInt lds,PetscScalar *H,PetscInt ldh,PetscScalar *fH,BV dV,PetscScalar *dVS,PetscInt *rds,PetscScalar *dH,PetscInt lddh,KSP ksp,PEP_REFINE_EXPLICIT *matctx)
{
  PetscInt            i,nmat=pep->nmat,lda=nmat*k;
  PetscScalar         *fh,*Rh,*DfH;
  PetscReal           norm;
  BV                  W;
  Vec                 Rv,t,dvi;
  PEP_REFINE_MATSHELL *ctx;
  Mat                 M,*At;
  PetscBool           flg,lindep;

  PetscFunctionBegin;
  PetscCall(PetscMalloc2(nmat*k*k,&DfH,k,&Rh));
  *rds = 0;
  PetscCall(BVCreateVec(pep->V,&Rv));
  switch (pep->scheme) {
  case PEP_REFINE_SCHEME_EXPLICIT:
    PetscCall(BVCreateVec(pep->V,&t));
    PetscCall(BVDuplicateResize(pep->V,PetscMax(k,nmat),&W));
    PetscCall(PetscMalloc1(nmat,&fh));
    break;
  case PEP_REFINE_SCHEME_MBE:
    if (matctx->subc) {
      PetscCall(BVCreateVec(pep->V,&t));
      PetscCall(BVDuplicateResize(pep->V,PetscMax(k,nmat),&W));
    } else {
      W = matctx->W;
      PetscCall(PetscObjectReference((PetscObject)W));
      t = matctx->t;
      PetscCall(PetscObjectReference((PetscObject)t));
    }
    PetscCall(BVScale(matctx->W,0.0));
    PetscCall(BVScale(matctx->Wt,0.0));
    PetscCall(BVScale(matctx->M2,0.0));
    PetscCall(BVScale(matctx->M3,0.0));
    PetscCall(PetscMalloc1(nmat,&fh));
    break;
  case PEP_REFINE_SCHEME_SCHUR:
    PetscCall(KSPGetOperators(ksp,&M,NULL));
    PetscCall(MatShellGetContext(M,&ctx));
    PetscCall(BVCreateVec(pep->V,&t));
    PetscCall(BVDuplicateResize(pep->V,PetscMax(k,nmat),&W));
    fh = ctx->fih;
    break;
  }
  PetscCall(PetscArrayzero(dVS,2*k*k));
  PetscCall(PetscArrayzero(DfH,lda*k));
  PetscCall(STGetTransform(pep->st,&flg));
  if (flg) {
    PetscCall(PetscMalloc1(pep->nmat,&At));
    for (i=0;i<pep->nmat;i++) PetscCall(STGetMatrixTransformed(pep->st,i,&At[i]));
  } else At = pep->A;

  /* Main loop for computing the i-th columns of dX and dS */
  for (i=0;i<k;i++) {
    /* Compute and update i-th column of the right hand side */
    PetscCall(PetscArrayzero(Rh,k));
    PetscCall(NRefRightSide(nmat,pep->pbc,At,k,pep->V,S,lds,i,H,ldh,fH,DfH,dH,dV,dVS,*rds,Rv,Rh,W,t));

    /* Update and solve system */
    PetscCall(BVGetColumn(dV,i,&dvi));
    PetscCall(NRefSysIter(i,pep,k,ksp,fH,S,lds,fh,H,ldh,Rv,Rh,pep->V,dvi,dH+i*k,matctx,W));
    /* Orthogonalize computed solution */
    PetscCall(BVOrthogonalizeVec(pep->V,dvi,dVS+i*2*k,&norm,&lindep));
    PetscCall(BVRestoreColumn(dV,i,&dvi));
    if (!lindep) {
      PetscCall(BVOrthogonalizeColumn(dV,i,dVS+k+i*2*k,&norm,&lindep));
      if (!lindep) {
        dVS[k+i+i*2*k] = norm;
        PetscCall(BVScaleColumn(dV,i,1.0/norm));
        (*rds)++;
      }
    }
  }
  PetscCall(BVSetActiveColumns(dV,0,*rds));
  PetscCall(VecDestroy(&t));
  PetscCall(VecDestroy(&Rv));
  PetscCall(BVDestroy(&W));
  if (flg) PetscCall(PetscFree(At));
  PetscCall(PetscFree2(DfH,Rh));
  if (pep->scheme!=PEP_REFINE_SCHEME_SCHUR) PetscCall(PetscFree(fh));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefOrthogStep(PEP pep,PetscInt k,PetscScalar *H,PetscInt ldh,PetscScalar *fH,PetscScalar *S,PetscInt lds)
{
  PetscInt       j,nmat=pep->nmat,deg=nmat-1,lda=nmat*k,ldg;
  PetscScalar    *G,*tau,sone=1.0,zero=0.0,*work;
  PetscBLASInt   lds_,k_,ldh_,info,ldg_,lda_;

  PetscFunctionBegin;
  PetscCall(PetscMalloc3(k,&tau,k,&work,deg*k*k,&G));
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(lda,&lda_));
  PetscCall(PetscBLASIntCast(k,&k_));

  /* Form auxiliary matrix for the orthogonalization step */
  ldg = deg*k;
  PetscCall(PEPEvaluateBasisforMatrix(pep,nmat,k,H,ldh,fH));
  PetscCall(PetscBLASIntCast(ldg,&ldg_));
  PetscCall(PetscBLASIntCast(ldh,&ldh_));
  for (j=0;j<deg;j++) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,S,&lds_,fH+j*k,&lda_,&zero,G+j*k,&ldg_));
  }
  /* Orthogonalize and update S */
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&ldg_,&k_,G,&ldg_,tau,work,&k_,&info));
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("geqrf",info);
  PetscCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&k_,&k_,&sone,G,&ldg_,S,&lds_));

  /* Update H */
  PetscCallBLAS("BLAStrmm",BLAStrmm_("L","U","N","N",&k_,&k_,&sone,G,&ldg_,H,&ldh_));
  PetscCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&k_,&k_,&sone,G,&ldg_,H,&ldh_));
  PetscCall(PetscFree3(tau,work,G));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPNRefUpdateInvPair(PEP pep,PetscInt k,PetscScalar *H,PetscInt ldh,PetscScalar *fH,PetscScalar *dH,PetscScalar *S,PetscInt lds,BV dV,PetscScalar *dVS,PetscInt rds)
{
  PetscInt       i,j,nmat=pep->nmat,lda=nmat*k;
  PetscScalar    *tau,*array,*work;
  PetscBLASInt   lds_,k_,lda_,ldh_,kdrs_,info,k2_;
  Mat            M0;

  PetscFunctionBegin;
  PetscCall(PetscMalloc2(k,&tau,k,&work));
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(lda,&lda_));
  PetscCall(PetscBLASIntCast(ldh,&ldh_));
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(2*k,&k2_));
  PetscCall(PetscBLASIntCast((k+rds),&kdrs_));
  /* Update H */
  for (j=0;j<k;j++) {
    for (i=0;i<k;i++) H[i+j*ldh] -= dH[i+j*k];
  }
  /* Update V */
  for (j=0;j<k;j++) {
    for (i=0;i<k;i++) dVS[i+j*2*k] = -dVS[i+j*2*k]+S[i+j*lds];
    for (i=k;i<2*k;i++) dVS[i+j*2*k] = -dVS[i+j*2*k];
  }
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&kdrs_,&k_,dVS,&k2_,tau,work,&k_,&info));
  SlepcCheckLapackInfo("geqrf",info);
  /* Copy triangular matrix in S */
  for (j=0;j<k;j++) {
    for (i=0;i<=j;i++) S[i+j*lds] = dVS[i+j*2*k];
    for (i=j+1;i<k;i++) S[i+j*lds] = 0.0;
  }
  PetscCallBLAS("LAPACKorgqr",LAPACKorgqr_(&k2_,&k_,&k_,dVS,&k2_,tau,work,&k_,&info));
  SlepcCheckLapackInfo("orgqr",info);
  PetscCall(PetscFPTrapPop());
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M0));
  PetscCall(MatDenseGetArrayWrite(M0,&array));
  for (j=0;j<k;j++) PetscCall(PetscArraycpy(array+j*k,dVS+j*2*k,k));
  PetscCall(MatDenseRestoreArrayWrite(M0,&array));
  PetscCall(BVMultInPlace(pep->V,M0,0,k));
  if (rds) {
    PetscCall(MatDenseGetArrayWrite(M0,&array));
    for (j=0;j<k;j++) PetscCall(PetscArraycpy(array+j*k,dVS+k+j*2*k,rds));
    PetscCall(MatDenseRestoreArrayWrite(M0,&array));
    PetscCall(BVMultInPlace(dV,M0,0,k));
    PetscCall(BVMult(pep->V,1.0,1.0,dV,NULL));
  }
  PetscCall(MatDestroy(&M0));
  PetscCall(NRefOrthogStep(pep,k,H,ldh,fH,S,lds));
  PetscCall(PetscFree2(tau,work));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPNRefSetUp(PEP pep,PetscInt k,PetscScalar *H,PetscInt ldh,PEP_REFINE_EXPLICIT *matctx,PetscBool ini)
{
  PEP_REFINE_MATSHELL *ctx;
  PetscScalar         t,*coef;
  const PetscScalar   *array;
  MatStructure        str;
  PetscInt            j,nmat=pep->nmat,n0,m0,n1,m1,n0_,m0_,n1_,m1_,N0,N1,p,*idx1,*idx2,count,si,i,l0;
  MPI_Comm            comm;
  PetscMPIInt         np;
  const PetscInt      *rgs0,*rgs1;
  Mat                 B,C,*E,*A,*At;
  IS                  is1,is2;
  Vec                 v;
  PetscBool           flg;
  Mat                 M,P;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(nmat,&coef));
  PetscCall(STGetTransform(pep->st,&flg));
  if (flg) {
    PetscCall(PetscMalloc1(pep->nmat,&At));
    for (i=0;i<pep->nmat;i++) PetscCall(STGetMatrixTransformed(pep->st,i,&At[i]));
  } else At = pep->A;
  switch (pep->scheme) {
  case PEP_REFINE_SCHEME_EXPLICIT:
    if (ini) {
      if (matctx->subc) {
        A = matctx->A;
        PetscCall(PetscSubcommGetChild(matctx->subc,&comm));
      } else {
        A = At;
        PetscCall(PetscObjectGetComm((PetscObject)pep,&comm));
      }
      E = matctx->E;
      PetscCall(STGetMatStructure(pep->st,&str));
      PetscCall(MatDuplicate(A[0],MAT_COPY_VALUES,&E[0]));
      j = (matctx->subc)?matctx->subc->color:0;
      PetscCall(PEPEvaluateBasis(pep,H[j+j*ldh],0,coef,NULL));
      for (j=1;j<nmat;j++) PetscCall(MatAXPY(E[0],coef[j],A[j],str));
      PetscCall(MatCreateDense(comm,PETSC_DECIDE,PETSC_DECIDE,k,k,NULL,&E[1]));
      PetscCall(MatGetOwnershipRange(E[0],&n0,&m0));
      PetscCall(MatGetOwnershipRange(E[1],&n1,&m1));
      PetscCall(MatGetOwnershipRangeColumn(E[0],&n0_,&m0_));
      PetscCall(MatGetOwnershipRangeColumn(E[1],&n1_,&m1_));
      /* T12 and T21 are computed from V and V*, so,
         they must have the same column and row ranges */
      PetscCheck(m0_-n0_==m0-n0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Inconsistent dimensions");
      PetscCall(MatCreateDense(comm,m0-n0,m1_-n1_,PETSC_DECIDE,PETSC_DECIDE,NULL,&B));
      PetscCall(MatCreateDense(comm,m1-n1,m0_-n0_,PETSC_DECIDE,PETSC_DECIDE,NULL,&C));
      PetscCall(MatCreateTile(1.0,E[0],1.0,B,1.0,C,1.0,E[1],&M));
      PetscCall(MatDestroy(&B));
      PetscCall(MatDestroy(&C));
      matctx->compM1 = PETSC_TRUE;
      PetscCall(MatGetSize(E[0],NULL,&N0));
      PetscCall(MatGetSize(E[1],NULL,&N1));
      PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)M),&np));
      PetscCall(MatGetOwnershipRanges(E[0],&rgs0));
      PetscCall(MatGetOwnershipRanges(E[1],&rgs1));
      PetscCall(PetscMalloc4(PetscMax(k,N1),&matctx->idxp,N0,&matctx->idxg,N0,&matctx->map0,N1,&matctx->map1));
      /* Create column (and row) mapping */
      for (p=0;p<np;p++) {
        for (j=rgs0[p];j<rgs0[p+1];j++) matctx->map0[j] = j+rgs1[p];
        for (j=rgs1[p];j<rgs1[p+1];j++) matctx->map1[j] = j+rgs0[p+1];
      }
      PetscCall(MatCreateVecs(M,NULL,&matctx->tN));
      PetscCall(MatCreateVecs(matctx->E[1],NULL,&matctx->t1));
      PetscCall(VecDuplicate(matctx->tN,&matctx->ttN));
      if (matctx->subc) {
        PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np));
        count = np*k;
        PetscCall(PetscMalloc2(count,&idx1,count,&idx2));
        PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)pep),m1-n1,PETSC_DECIDE,&matctx->tp));
        PetscCall(VecGetOwnershipRange(matctx->tp,&l0,NULL));
        PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)pep),k,PETSC_DECIDE,&matctx->tpg));
        for (si=0;si<matctx->subc->n;si++) {
          if (matctx->subc->color==si) {
            j=0;
            if (matctx->subc->color==si) {
              for (p=0;p<np;p++) {
                for (i=n1;i<m1;i++) {
                  idx1[j] = l0+i-n1;
                  idx2[j++] =p*k+i;
                }
              }
            }
            count = np*(m1-n1);
          } else count =0;
          PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pep),count,idx1,PETSC_COPY_VALUES,&is1));
          PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pep),count,idx2,PETSC_COPY_VALUES,&is2));
          PetscCall(VecScatterCreate(matctx->tp,is1,matctx->tpg,is2,&matctx->scatterp_id[si]));
          PetscCall(ISDestroy(&is1));
          PetscCall(ISDestroy(&is2));
        }
        PetscCall(PetscFree2(idx1,idx2));
      } else PetscCall(VecScatterCreateToAll(matctx->t1,&matctx->scatterctx,&matctx->vseq));
      P = M;
    } else {
      if (matctx->subc) {
        /* Scatter vectors pep->V */
        for (i=0;i<k;i++) {
          PetscCall(BVGetColumn(pep->V,i,&v));
          PetscCall(VecScatterBegin(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD));
          PetscCall(VecScatterEnd(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD));
          PetscCall(BVRestoreColumn(pep->V,i,&v));
          PetscCall(VecGetArrayRead(matctx->tg,&array));
          PetscCall(VecPlaceArray(matctx->t,(const PetscScalar*)array));
          PetscCall(BVInsertVec(matctx->V,i,matctx->t));
          PetscCall(VecResetArray(matctx->t));
          PetscCall(VecRestoreArrayRead(matctx->tg,&array));
        }
      }
    }
    break;
  case PEP_REFINE_SCHEME_MBE:
    if (ini) {
      if (matctx->subc) {
        A = matctx->A;
        PetscCall(PetscSubcommGetChild(matctx->subc,&comm));
      } else {
        matctx->V = pep->V;
        A = At;
        PetscCall(PetscObjectGetComm((PetscObject)pep,&comm));
        PetscCall(MatCreateVecs(pep->A[0],&matctx->t,NULL));
      }
      PetscCall(STGetMatStructure(pep->st,&str));
      PetscCall(MatDuplicate(A[0],MAT_COPY_VALUES,&matctx->M1));
      j = (matctx->subc)?matctx->subc->color:0;
      PetscCall(PEPEvaluateBasis(pep,H[j+j*ldh],0,coef,NULL));
      for (j=1;j<nmat;j++) PetscCall(MatAXPY(matctx->M1,coef[j],A[j],str));
      PetscCall(BVDuplicateResize(matctx->V,PetscMax(k,pep->nmat),&matctx->W));
      PetscCall(BVDuplicateResize(matctx->V,k,&matctx->M2));
      PetscCall(BVDuplicate(matctx->M2,&matctx->M3));
      PetscCall(BVDuplicate(matctx->M2,&matctx->Wt));
      PetscCall(PetscMalloc5(k*k,&matctx->M4,k*k,&matctx->w,k*k,&matctx->wt,k,&matctx->d,k,&matctx->dt));
      matctx->compM1 = PETSC_TRUE;
      M = matctx->M1;
      P = M;
    }
    break;
  case PEP_REFINE_SCHEME_SCHUR:
    if (ini) {
      PetscCall(PetscObjectGetComm((PetscObject)pep,&comm));
      PetscCall(MatGetSize(At[0],&m0,&n0));
      PetscCall(PetscMalloc1(1,&ctx));
      PetscCall(STGetMatStructure(pep->st,&str));
      /* Create a shell matrix to solve the linear system */
      ctx->V = pep->V;
      ctx->k = k; ctx->nmat = nmat;
      PetscCall(PetscMalloc5(nmat,&ctx->A,k*k,&ctx->M4,k,&ctx->pM4,2*k*k,&ctx->work,nmat,&ctx->fih));
      for (i=0;i<nmat;i++) ctx->A[i] = At[i];
      PetscCall(PetscArrayzero(ctx->M4,k*k));
      PetscCall(MatCreateShell(comm,PETSC_DECIDE,PETSC_DECIDE,m0,n0,ctx,&M));
      PetscCall(MatShellSetOperation(M,MATOP_MULT,(PetscErrorCodeFn*)MatMult_FS));
      PetscCall(BVDuplicateResize(ctx->V,PetscMax(k,pep->nmat),&ctx->W));
      PetscCall(BVDuplicateResize(ctx->V,k,&ctx->M2));
      PetscCall(BVDuplicate(ctx->M2,&ctx->M3));
      PetscCall(BVCreateVec(pep->V,&ctx->t));
      PetscCall(MatDuplicate(At[0],MAT_COPY_VALUES,&ctx->M1));
      PetscCall(PEPEvaluateBasis(pep,H[0],0,coef,NULL));
      for (j=1;j<nmat;j++) PetscCall(MatAXPY(ctx->M1,coef[j],At[j],str));
      PetscCall(MatDuplicate(At[0],MAT_COPY_VALUES,&P));
      /* Compute a precond matrix for the system */
      t = H[0];
      PetscCall(PEPEvaluateBasis(pep,t,0,coef,NULL));
      for (j=1;j<nmat;j++) PetscCall(MatAXPY(P,coef[j],At[j],str));
      ctx->compM1 = PETSC_TRUE;
    }
    break;
  }
  if (ini) {
    PetscCall(PEPRefineGetKSP(pep,&pep->refineksp));
    PetscCall(KSPSetErrorIfNotConverged(pep->refineksp,PETSC_TRUE));
    PetscCall(PEP_KSPSetOperators(pep->refineksp,M,P));
    PetscCall(KSPSetFromOptions(pep->refineksp));
  }

  if (!ini && matctx && matctx->subc) {
     /* Scatter vectors pep->V */
    for (i=0;i<k;i++) {
      PetscCall(BVGetColumn(pep->V,i,&v));
      PetscCall(VecScatterBegin(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecScatterEnd(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(BVRestoreColumn(pep->V,i,&v));
      PetscCall(VecGetArrayRead(matctx->tg,&array));
      PetscCall(VecPlaceArray(matctx->t,(const PetscScalar*)array));
      PetscCall(BVInsertVec(matctx->V,i,matctx->t));
      PetscCall(VecResetArray(matctx->t));
      PetscCall(VecRestoreArrayRead(matctx->tg,&array));
    }
   }
  PetscCall(PetscFree(coef));
  if (flg) PetscCall(PetscFree(At));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSubcommSetup(PEP pep,PetscInt k,PEP_REFINE_EXPLICIT *matctx,PetscInt nsubc)
{
  PetscInt          i,si,j,m0,n0,nloc0,nloc_sub,*idx1,*idx2;
  IS                is1,is2;
  BVType            type;
  Vec               v;
  const PetscScalar *array;
  Mat               *A;
  PetscBool         flg;
  MPI_Comm          contpar,child;

  PetscFunctionBegin;
  PetscCall(STGetTransform(pep->st,&flg));
  if (flg) {
    PetscCall(PetscMalloc1(pep->nmat,&A));
    for (i=0;i<pep->nmat;i++) PetscCall(STGetMatrixTransformed(pep->st,i,&A[i]));
  } else A = pep->A;
  PetscCall(PetscSubcommGetChild(matctx->subc,&child));
  PetscCall(PetscSubcommGetContiguousParent(matctx->subc,&contpar));

  /* Duplicate pep matrices */
  PetscCall(PetscMalloc3(pep->nmat,&matctx->A,nsubc,&matctx->scatter_id,nsubc,&matctx->scatterp_id));
  for (i=0;i<pep->nmat;i++) PetscCall(MatCreateRedundantMatrix(A[i],0,child,MAT_INITIAL_MATRIX,&matctx->A[i]));

  /* Create Scatter */
  PetscCall(MatCreateVecs(matctx->A[0],&matctx->t,NULL));
  PetscCall(MatGetLocalSize(matctx->A[0],&nloc_sub,NULL));
  PetscCall(VecCreateMPI(contpar,nloc_sub,PETSC_DECIDE,&matctx->tg));
  PetscCall(BVGetColumn(pep->V,0,&v));
  PetscCall(VecGetOwnershipRange(v,&n0,&m0));
  nloc0 = m0-n0;
  PetscCall(PetscMalloc2(matctx->subc->n*nloc0,&idx1,matctx->subc->n*nloc0,&idx2));
  j = 0;
  for (si=0;si<matctx->subc->n;si++) {
    for (i=n0;i<m0;i++) {
      idx1[j]   = i;
      idx2[j++] = i+pep->n*si;
    }
  }
  PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pep),matctx->subc->n*nloc0,idx1,PETSC_COPY_VALUES,&is1));
  PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pep),matctx->subc->n*nloc0,idx2,PETSC_COPY_VALUES,&is2));
  PetscCall(VecScatterCreate(v,is1,matctx->tg,is2,&matctx->scatter_sub));
  PetscCall(ISDestroy(&is1));
  PetscCall(ISDestroy(&is2));
  for (si=0;si<matctx->subc->n;si++) {
    j=0;
    for (i=n0;i<m0;i++) {
      idx1[j] = i;
      idx2[j++] = i+pep->n*si;
    }
    PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pep),nloc0,idx1,PETSC_COPY_VALUES,&is1));
    PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pep),nloc0,idx2,PETSC_COPY_VALUES,&is2));
    PetscCall(VecScatterCreate(v,is1,matctx->tg,is2,&matctx->scatter_id[si]));
    PetscCall(ISDestroy(&is1));
    PetscCall(ISDestroy(&is2));
  }
  PetscCall(BVRestoreColumn(pep->V,0,&v));
  PetscCall(PetscFree2(idx1,idx2));

  /* Duplicate pep->V vecs */
  PetscCall(BVGetType(pep->V,&type));
  PetscCall(BVCreate(child,&matctx->V));
  PetscCall(BVSetType(matctx->V,type));
  PetscCall(BVSetSizesFromVec(matctx->V,matctx->t,k));
  if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) PetscCall(BVDuplicateResize(matctx->V,PetscMax(k,pep->nmat),&matctx->W));
  for (i=0;i<k;i++) {
    PetscCall(BVGetColumn(pep->V,i,&v));
    PetscCall(VecScatterBegin(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(BVRestoreColumn(pep->V,i,&v));
    PetscCall(VecGetArrayRead(matctx->tg,&array));
    PetscCall(VecPlaceArray(matctx->t,(const PetscScalar*)array));
    PetscCall(BVInsertVec(matctx->V,i,matctx->t));
    PetscCall(VecResetArray(matctx->t));
    PetscCall(VecRestoreArrayRead(matctx->tg,&array));
  }

  PetscCall(VecDuplicate(matctx->t,&matctx->Rv));
  PetscCall(VecDuplicate(matctx->t,&matctx->Vi));
  if (flg) PetscCall(PetscFree(A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NRefSubcommDestroy(PEP pep,PEP_REFINE_EXPLICIT *matctx)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscCall(VecScatterDestroy(&matctx->scatter_sub));
  for (i=0;i<matctx->subc->n;i++) PetscCall(VecScatterDestroy(&matctx->scatter_id[i]));
  for (i=0;i<pep->nmat;i++) PetscCall(MatDestroy(&matctx->A[i]));
  if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) {
    for (i=0;i<matctx->subc->n;i++) PetscCall(VecScatterDestroy(&matctx->scatterp_id[i]));
    PetscCall(VecDestroy(&matctx->tp));
    PetscCall(VecDestroy(&matctx->tpg));
    PetscCall(BVDestroy(&matctx->W));
  }
  PetscCall(PetscFree3(matctx->A,matctx->scatter_id,matctx->scatterp_id));
  PetscCall(BVDestroy(&matctx->V));
  PetscCall(VecDestroy(&matctx->t));
  PetscCall(VecDestroy(&matctx->tg));
  PetscCall(VecDestroy(&matctx->Rv));
  PetscCall(VecDestroy(&matctx->Vi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPNewtonRefinement_TOAR(PEP pep,PetscScalar sigma,PetscInt *maxits,PetscReal *tol,PetscInt k,PetscScalar *S,PetscInt lds)
{
  PetscScalar         *H,*work,*dH,*fH,*dVS;
  PetscInt            ldh,i,j,its=1,nmat=pep->nmat,nsubc=pep->npart,rds;
  PetscBLASInt        k_,ld_,*p,info;
  BV                  dV;
  PetscBool           sinvert,flg;
  PEP_REFINE_EXPLICIT *matctx=NULL;
  Vec                 v;
  Mat                 M,P;
  PEP_REFINE_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(PetscLogEventBegin(PEP_Refine,pep,0,0,0));
  PetscCheck(k<=pep->n,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Multiple Refinement available only for invariant pairs of dimension smaller than n=%" PetscInt_FMT,pep->n);
  /* the input tolerance is not being taken into account (by the moment) */
  its = *maxits;
  PetscCall(PetscMalloc3(k*k,&dH,nmat*k*k,&fH,k,&work));
  PetscCall(DSGetLeadingDimension(pep->ds,&ldh));
  PetscCall(PetscMalloc1(2*k*k,&dVS));
  PetscCall(STGetTransform(pep->st,&flg));
  if (!flg && pep->st && pep->ops->backtransform) { /* BackTransform */
    PetscCall(PetscBLASIntCast(k,&k_));
    PetscCall(PetscBLASIntCast(ldh,&ld_));
    PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinvert));
    if (sinvert) {
      PetscCall(DSGetArray(pep->ds,DS_MAT_A,&H));
      PetscCall(PetscMalloc1(k,&p));
      PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
      PetscCallBLAS("LAPACKgetrf",LAPACKgetrf_(&k_,&k_,H,&ld_,p,&info));
      SlepcCheckLapackInfo("getrf",info);
      PetscCallBLAS("LAPACKgetri",LAPACKgetri_(&k_,H,&ld_,p,work,&k_,&info));
      SlepcCheckLapackInfo("getri",info);
      PetscCall(PetscFPTrapPop());
      PetscCall(DSRestoreArray(pep->ds,DS_MAT_A,&H));
      pep->ops->backtransform = NULL;
    }
    if (sigma!=0.0) {
      PetscCall(DSGetArray(pep->ds,DS_MAT_A,&H));
      for (i=0;i<k;i++) H[i+ldh*i] += sigma;
      PetscCall(DSRestoreArray(pep->ds,DS_MAT_A,&H));
      pep->ops->backtransform = NULL;
    }
  }
  if ((pep->scale==PEP_SCALE_BOTH || pep->scale==PEP_SCALE_SCALAR) && pep->sfactor!=1.0) {
    PetscCall(DSGetArray(pep->ds,DS_MAT_A,&H));
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) H[i+j*ldh] *= pep->sfactor;
    }
    PetscCall(DSRestoreArray(pep->ds,DS_MAT_A,&H));
    if (!flg) {
      /* Restore original values */
      for (i=0;i<pep->nmat;i++) {
        pep->pbc[pep->nmat+i] *= pep->sfactor;
        pep->pbc[2*pep->nmat+i] *= pep->sfactor*pep->sfactor;
      }
    }
  }
  if ((pep->scale==PEP_SCALE_DIAGONAL || pep->scale==PEP_SCALE_BOTH) && pep->Dr) {
    for (i=0;i<k;i++) {
      PetscCall(BVGetColumn(pep->V,i,&v));
      PetscCall(VecPointwiseMult(v,v,pep->Dr));
      PetscCall(BVRestoreColumn(pep->V,i,&v));
    }
  }
  PetscCall(DSGetArray(pep->ds,DS_MAT_A,&H));

  PetscCall(NRefOrthogStep(pep,k,H,ldh,fH,S,lds));
  /* check if H is in Schur form */
  for (i=0;i<k-1;i++) {
#if !defined(PETSC_USE_COMPLEX)
    PetscCheck(H[i+1+i*ldh]==0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Iterative Refinement requires the complex Schur form of the projected matrix");
#else
    PetscCheck(H[i+1+i*ldh]==0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Iterative Refinement requires an upper triangular projected matrix");
#endif
  }
  PetscCheck(nsubc<=k,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Number of subcommunicators should not be larger than the invariant pair dimension");
  PetscCall(BVSetActiveColumns(pep->V,0,k));
  PetscCall(BVDuplicateResize(pep->V,k,&dV));
  if (pep->scheme!=PEP_REFINE_SCHEME_SCHUR) {
    PetscCall(PetscMalloc1(1,&matctx));
    if (nsubc>1) { /* splitting in subcommunicators */
      matctx->subc = pep->refinesubc;
      PetscCall(NRefSubcommSetup(pep,k,matctx,nsubc));
    } else matctx->subc=NULL;
  }

  /* Loop performing iterative refinements */
  for (i=0;i<its;i++) {
    /* Pre-compute the polynomial basis evaluated in H */
    PetscCall(PEPEvaluateBasisforMatrix(pep,nmat,k,H,ldh,fH));
    PetscCall(PEPNRefSetUp(pep,k,H,ldh,matctx,PetscNot(i)));
    /* Solve the linear system */
    PetscCall(PEPNRefForwardSubstitution(pep,k,S,lds,H,ldh,fH,dV,dVS,&rds,dH,k,pep->refineksp,matctx));
    /* Update X (=V*S) and H, and orthogonalize [X;X*fH1;...;XfH(deg-1)] */
    PetscCall(PEPNRefUpdateInvPair(pep,k,H,ldh,fH,dH,S,lds,dV,dVS,rds));
  }
  PetscCall(DSRestoreArray(pep->ds,DS_MAT_A,&H));
  if (!flg && sinvert) PetscCall(PetscFree(p));
  PetscCall(PetscFree3(dH,fH,work));
  PetscCall(PetscFree(dVS));
  PetscCall(BVDestroy(&dV));
  switch (pep->scheme) {
  case PEP_REFINE_SCHEME_EXPLICIT:
    for (i=0;i<2;i++) PetscCall(MatDestroy(&matctx->E[i]));
    PetscCall(PetscFree4(matctx->idxp,matctx->idxg,matctx->map0,matctx->map1));
    PetscCall(VecDestroy(&matctx->tN));
    PetscCall(VecDestroy(&matctx->ttN));
    PetscCall(VecDestroy(&matctx->t1));
    if (nsubc>1) PetscCall(NRefSubcommDestroy(pep,matctx));
    else {
      PetscCall(VecDestroy(&matctx->vseq));
      PetscCall(VecScatterDestroy(&matctx->scatterctx));
    }
    PetscCall(PetscFree(matctx));
    PetscCall(KSPGetOperators(pep->refineksp,&M,NULL));
    PetscCall(MatDestroy(&M));
    break;
  case PEP_REFINE_SCHEME_MBE:
    PetscCall(BVDestroy(&matctx->W));
    PetscCall(BVDestroy(&matctx->Wt));
    PetscCall(BVDestroy(&matctx->M2));
    PetscCall(BVDestroy(&matctx->M3));
    PetscCall(MatDestroy(&matctx->M1));
    PetscCall(VecDestroy(&matctx->t));
    PetscCall(PetscFree5(matctx->M4,matctx->w,matctx->wt,matctx->d,matctx->dt));
    if (nsubc>1) PetscCall(NRefSubcommDestroy(pep,matctx));
    PetscCall(PetscFree(matctx));
    break;
  case PEP_REFINE_SCHEME_SCHUR:
    PetscCall(KSPGetOperators(pep->refineksp,&M,&P));
    PetscCall(MatShellGetContext(M,&ctx));
    PetscCall(PetscFree5(ctx->A,ctx->M4,ctx->pM4,ctx->work,ctx->fih));
    PetscCall(MatDestroy(&ctx->M1));
    PetscCall(BVDestroy(&ctx->M2));
    PetscCall(BVDestroy(&ctx->M3));
    PetscCall(BVDestroy(&ctx->W));
    PetscCall(VecDestroy(&ctx->t));
    PetscCall(PetscFree(ctx));
    PetscCall(MatDestroy(&M));
    PetscCall(MatDestroy(&P));
    break;
  }
  PetscCall(PetscLogEventEnd(PEP_Refine,pep,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}
