/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "nleigs"

   Method: NLEIGS

   Algorithm:

       Fully rational Krylov method for nonlinear eigenvalue problems.

   References:

       [1] S. Guttel et al., "NLEIGS: A class of robust fully rational Krylov
           method for nonlinear eigenvalue problems", SIAM J. Sci. Comput.
           36(6):A2842-A2864, 2014.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <slepcblaslapack.h>
#include "nleigs.h"

PetscErrorCode NEPNLEIGSBackTransform(PetscObject ob,PetscInt n,PetscScalar *valr,PetscScalar *vali)
{
  NEP         nep;
  PetscInt    j;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar t;
#endif

  PetscFunctionBegin;
  nep = (NEP)ob;
#if !defined(PETSC_USE_COMPLEX)
  for (j=0;j<n;j++) {
    if (vali[j] == 0) valr[j] = 1.0 / valr[j] + nep->target;
    else {
      t = valr[j] * valr[j] + vali[j] * vali[j];
      valr[j] = valr[j] / t + nep->target;
      vali[j] = - vali[j] / t;
    }
  }
#else
  for (j=0;j<n;j++) {
    valr[j] = 1.0 / valr[j] + nep->target;
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Computes the roots of a polynomial */
static PetscErrorCode NEPNLEIGSAuxiliarPRootFinder(PetscInt deg,PetscScalar *polcoeffs,PetscScalar *wr,PetscScalar *wi,PetscBool *avail)
{
  PetscScalar    *C;
  PetscBLASInt   n_,lwork;
  PetscInt       i;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork=NULL;
#endif
  PetscScalar    *work;
  PetscBLASInt   info;

  PetscFunctionBegin;
  *avail = PETSC_TRUE;
  if (deg>0) {
    PetscCall(PetscCalloc1(deg*deg,&C));
    PetscCall(PetscBLASIntCast(deg,&n_));
    for (i=0;i<deg-1;i++) {
      C[(deg+1)*i+1]   = 1.0;
      C[(deg-1)*deg+i] = -polcoeffs[deg-i]/polcoeffs[0];
    }
    C[deg*deg+-1] = -polcoeffs[1]/polcoeffs[0];
    PetscCall(PetscBLASIntCast(3*deg,&lwork));

    PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(PetscMalloc1(lwork,&work));
    PetscCallBLAS("LAPACKgeev",LAPACKgeev_("N","N",&n_,C,&n_,wr,wi,NULL,&n_,NULL,&n_,work,&lwork,&info));
    if (info) *avail = PETSC_FALSE;
    PetscCall(PetscFree(work));
#else
    PetscCall(PetscMalloc2(2*deg,&rwork,lwork,&work));
    PetscCallBLAS("LAPACKgeev",LAPACKgeev_("N","N",&n_,C,&n_,wr,NULL,&n_,NULL,&n_,work,&lwork,rwork,&info));
    if (info) *avail = PETSC_FALSE;
    PetscCall(PetscFree2(rwork,work));
#endif
    PetscCall(PetscFPTrapPop());
    PetscCall(PetscFree(C));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSAuxiliarRmDuplicates(PetscInt nin,PetscScalar *pin,PetscInt *nout,PetscScalar *pout,PetscInt max)
{
  PetscInt i,j;

  PetscFunctionBegin;
  for (i=0;i<nin;i++) {
    if (max && *nout>=max) break;
    pout[(*nout)++] = pin[i];
    for (j=0;j<*nout-1;j++)
      if (PetscAbsScalar(pin[i]-pout[j])<PETSC_MACHINE_EPSILON*100) {
        (*nout)--;
        break;
      }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSFNSingularities(FN f,PetscInt *nisol,PetscScalar **isol,PetscBool *rational)
{
  FNCombineType  ctype;
  FN             f1,f2;
  PetscInt       i,nq,nisol1,nisol2;
  PetscScalar    *qcoeff,*wr,*wi,*isol1,*isol2;
  PetscBool      flg,avail,rat1,rat2;

  PetscFunctionBegin;
  *rational = PETSC_FALSE;
  PetscCall(PetscObjectTypeCompare((PetscObject)f,FNRATIONAL,&flg));
  if (flg) {
    *rational = PETSC_TRUE;
    PetscCall(FNRationalGetDenominator(f,&nq,&qcoeff));
    if (nq>1) {
      PetscCall(PetscMalloc2(nq-1,&wr,nq-1,&wi));
      PetscCall(NEPNLEIGSAuxiliarPRootFinder(nq-1,qcoeff,wr,wi,&avail));
      if (avail) {
        PetscCall(PetscCalloc1(nq-1,isol));
        *nisol = 0;
        for (i=0;i<nq-1;i++)
#if !defined(PETSC_USE_COMPLEX)
          if (wi[i]==0)
#endif
            (*isol)[(*nisol)++] = wr[i];
        nq = *nisol; *nisol = 0;
        for (i=0;i<nq;i++) wr[i] = (*isol)[i];
        PetscCall(NEPNLEIGSAuxiliarRmDuplicates(nq,wr,nisol,*isol,0));
        PetscCall(PetscFree2(wr,wi));
      } else { *nisol=0; *isol = NULL; }
    } else { *nisol = 0; *isol = NULL; }
    PetscCall(PetscFree(qcoeff));
  }
  PetscCall(PetscObjectTypeCompare((PetscObject)f,FNCOMBINE,&flg));
  if (flg) {
    PetscCall(FNCombineGetChildren(f,&ctype,&f1,&f2));
    if (ctype != FN_COMBINE_COMPOSE && ctype != FN_COMBINE_DIVIDE) {
      PetscCall(NEPNLEIGSFNSingularities(f1,&nisol1,&isol1,&rat1));
      PetscCall(NEPNLEIGSFNSingularities(f2,&nisol2,&isol2,&rat2));
      if (nisol1+nisol2>0) {
        PetscCall(PetscCalloc1(nisol1+nisol2,isol));
        *nisol = 0;
        PetscCall(NEPNLEIGSAuxiliarRmDuplicates(nisol1,isol1,nisol,*isol,0));
        PetscCall(NEPNLEIGSAuxiliarRmDuplicates(nisol2,isol2,nisol,*isol,0));
      }
      *rational = (rat1&&rat2)?PETSC_TRUE:PETSC_FALSE;
      PetscCall(PetscFree(isol1));
      PetscCall(PetscFree(isol2));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSRationalSingularities(NEP nep,PetscInt *ndptx,PetscScalar *dxi,PetscBool *rational)
{
  PetscInt       nt,i,nisol;
  FN             f;
  PetscScalar    *isol;
  PetscBool      rat;

  PetscFunctionBegin;
  *rational = PETSC_TRUE;
  *ndptx = 0;
  PetscCall(NEPGetSplitOperatorInfo(nep,&nt,NULL));
  for (i=0;i<nt;i++) {
    PetscCall(NEPGetSplitOperatorTerm(nep,i,NULL,&f));
    PetscCall(NEPNLEIGSFNSingularities(f,&nisol,&isol,&rat));
    if (nisol) {
      PetscCall(NEPNLEIGSAuxiliarRmDuplicates(nisol,isol,ndptx,dxi,0));
      PetscCall(PetscFree(isol));
    }
    *rational = ((*rational)&&rat)?PETSC_TRUE:PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(SLEPC_MISSING_LAPACK_GGEV3)
#define LAPGEEV "ggev"
#else
#define LAPGEEV "ggev3"
#endif

/*  Adaptive Anderson-Antoulas algorithm */
static PetscErrorCode NEPNLEIGSAAAComputation(NEP nep,PetscInt ndpt,PetscScalar *ds,PetscScalar *F,PetscInt *ndptx,PetscScalar *dxi)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscScalar    mean=0.0,*z,*f,*C,*A,*VT,*work,*ww,szero=0.0,sone=1.0;
  PetscScalar    *N,*D;
  PetscReal      *S,norm,err,*R;
  PetscInt       i,k,j,idx=0,cont;
  PetscBLASInt   n_,m_,lda_,lwork,info,one=1;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(8*ndpt,&lwork));
  PetscCall(PetscMalloc5(ndpt,&R,ndpt,&z,ndpt,&f,ndpt*ndpt,&C,ndpt,&ww));
  PetscCall(PetscMalloc6(ndpt*ndpt,&A,ndpt,&S,ndpt*ndpt,&VT,lwork,&work,ndpt,&D,ndpt,&N));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscMalloc1(8*ndpt,&rwork));
#endif
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  norm = 0.0;
  for (i=0;i<ndpt;i++) {
    mean += F[i];
    norm = PetscMax(PetscAbsScalar(F[i]),norm);
  }
  mean /= ndpt;
  PetscCall(PetscBLASIntCast(ndpt,&lda_));
  for (i=0;i<ndpt;i++) R[i] = PetscAbsScalar(F[i]-mean);
  /* next support point */
  err = 0.0;
  for (i=0;i<ndpt;i++) if (R[i]>=err) {idx = i; err = R[i];}
  for (k=0;k<ndpt-1;k++) {
    z[k] = ds[idx]; f[k] = F[idx]; R[idx] = -1.0;
    /* next column of Cauchy matrix */
    for (i=0;i<ndpt;i++) {
      C[i+k*ndpt] = 1.0/(ds[i]-ds[idx]);
    }

    PetscCall(PetscArrayzero(A,ndpt*ndpt));
    cont = 0;
    for (i=0;i<ndpt;i++) {
      if (R[i]!=-1.0) {
        for (j=0;j<=k;j++)A[cont+j*ndpt] = C[i+j*ndpt]*F[i]-C[i+j*ndpt]*f[j];
        cont++;
      }
    }
    PetscCall(PetscBLASIntCast(cont,&m_));
    PetscCall(PetscBLASIntCast(k+1,&n_));
#if defined(PETSC_USE_COMPLEX)
    PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","A",&m_,&n_,A,&lda_,S,NULL,&lda_,VT,&lda_,work,&lwork,rwork,&info));
#else
    PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","A",&m_,&n_,A,&lda_,S,NULL,&lda_,VT,&lda_,work,&lwork,&info));
#endif
    SlepcCheckLapackInfo("gesvd",info);
    for (i=0;i<=k;i++) {
      ww[i] = PetscConj(VT[i*ndpt+k]);
      D[i] = ww[i]*f[i];
    }
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&lda_,&n_,&sone,C,&lda_,D,&one,&szero,N,&one));
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&lda_,&n_,&sone,C,&lda_,ww,&one,&szero,D,&one));
    for (i=0;i<ndpt;i++) if (R[i]>=0) R[i] = PetscAbsScalar(F[i]-N[i]/D[i]);
    /* next support point */
    err = 0.0;
    for (i=0;i<ndpt;i++) if (R[i]>=err) {idx = i; err = R[i];}
    if (err <= ctx->ddtol*norm) break;
  }

  PetscCheck(k<ndpt-1,PetscObjectComm((PetscObject)nep),PETSC_ERR_CONV_FAILED,"Failed to determine singularities automatically in general problem");
  /* poles */
  PetscCall(PetscArrayzero(C,ndpt*ndpt));
  PetscCall(PetscArrayzero(A,ndpt*ndpt));
  for (i=0;i<=k;i++) {
    C[i+ndpt*i] = 1.0;
    A[(i+1)*ndpt] = ww[i];
    A[i+1] = 1.0;
    A[i+1+(i+1)*ndpt] = z[i];
  }
  C[0] = 0.0; C[k+1+(k+1)*ndpt] = 1.0;
  n_++;
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACK" LAPGEEV,LAPACKggevalt_("N","N",&n_,A,&lda_,C,&lda_,D,N,NULL,&lda_,NULL,&lda_,work,&lwork,rwork,&info));
#else
  PetscCallBLAS("LAPACK" LAPGEEV,LAPACKggevalt_("N","N",&n_,A,&lda_,C,&lda_,D,VT,N,NULL,&lda_,NULL,&lda_,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo(LAPGEEV,info);
  cont = 0.0;
  for (i=0;i<n_;i++) if (N[i]!=0.0) {
    dxi[cont++] = D[i]/N[i];
  }
  *ndptx = cont;
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscFree5(R,z,f,C,ww));
  PetscCall(PetscFree6(A,S,VT,work,D,N));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree(rwork));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*  Singularities using Adaptive Anderson-Antoulas algorithm */
static PetscErrorCode NEPNLEIGSAAASingularities(NEP nep,PetscInt ndpt,PetscScalar *ds,PetscInt *ndptx,PetscScalar *dxi)
{
  Vec            u,v,w;
  PetscRandom    rand=NULL;
  PetscScalar    *F,*isol;
  PetscInt       i,k,nisol,nt;
  Mat            T;
  FN             f;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(ndpt,&F));
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    PetscCall(PetscMalloc1(ndpt,&isol));
    *ndptx = 0;
    PetscCall(NEPGetSplitOperatorInfo(nep,&nt,NULL));
    nisol = *ndptx;
    for (k=0;k<nt;k++) {
      PetscCall(NEPGetSplitOperatorTerm(nep,k,NULL,&f));
      for (i=0;i<ndpt;i++) PetscCall(FNEvaluateFunction(f,ds[i],&F[i]));
      PetscCall(NEPNLEIGSAAAComputation(nep,ndpt,ds,F,&nisol,isol));
      if (nisol) PetscCall(NEPNLEIGSAuxiliarRmDuplicates(nisol,isol,ndptx,dxi,ndpt));
    }
    PetscCall(PetscFree(isol));
  } else {
    PetscCall(MatCreateVecs(nep->function,&u,NULL));
    PetscCall(VecDuplicate(u,&v));
    PetscCall(VecDuplicate(u,&w));
    if (nep->V) PetscCall(BVGetRandomContext(nep->V,&rand));
    PetscCall(VecSetRandom(u,rand));
    PetscCall(VecNormalize(u,NULL));
    PetscCall(VecSetRandom(v,rand));
    PetscCall(VecNormalize(v,NULL));
    T = nep->function;
    for (i=0;i<ndpt;i++) {
      PetscCall(NEPComputeFunction(nep,ds[i],T,T));
      PetscCall(MatMult(T,v,w));
      PetscCall(VecDot(w,u,&F[i]));
    }
    PetscCall(NEPNLEIGSAAAComputation(nep,ndpt,ds,F,ndptx,dxi));
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&v));
    PetscCall(VecDestroy(&w));
  }
  PetscCall(PetscFree(F));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSLejaBagbyPoints(NEP nep)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       i,k,ndpt=NDPOINTS,ndptx=NDPOINTS;
  PetscScalar    *ds,*dsi,*dxi,*nrs,*nrxi,*s=ctx->s,*xi=ctx->xi,*beta=ctx->beta;
  PetscReal      maxnrs,minnrxi;
  PetscBool      rational;
#if !defined(PETSC_USE_COMPLEX)
  PetscReal      a,b,h;
#endif

  PetscFunctionBegin;
  if (!ctx->computesingularities && nep->problem_type!=NEP_RATIONAL) ndpt = ndptx = LBPOINTS;
  PetscCall(PetscMalloc5(ndpt+1,&ds,ndpt+1,&dsi,ndpt,&dxi,ndpt+1,&nrs,ndpt,&nrxi));

  /* Discretize the target region boundary */
  PetscCall(RGComputeContour(nep->rg,ndpt,ds,dsi));
#if !defined(PETSC_USE_COMPLEX)
  for (i=0;i<ndpt;i++) if (dsi[i]!=0.0) break;
  if (i<ndpt) {
    PetscCheck(nep->problem_type==NEP_RATIONAL,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NLEIGS with real arithmetic requires the target set to be included in the real axis");
    /* Select a segment in the real axis */
    PetscCall(RGComputeBoundingBox(nep->rg,&a,&b,NULL,NULL));
    PetscCheck(a>-PETSC_MAX_REAL && b<PETSC_MAX_REAL,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"NLEIGS requires a bounded target set");
    h = (b-a)/ndpt;
    for (i=0;i<ndpt;i++) {ds[i] = a+h*i; dsi[i] = 0.0;}
  }
#endif
  /* Discretize the singularity region */
  if (ctx->computesingularities) PetscCall(ctx->computesingularities(nep,&ndptx,dxi,ctx->singularitiesctx));
  else {
    if (nep->problem_type==NEP_RATIONAL) {
      PetscCall(NEPNLEIGSRationalSingularities(nep,&ndptx,dxi,&rational));
      PetscCheck(rational,PetscObjectComm((PetscObject)nep),PETSC_ERR_CONV_FAILED,"Failed to determine singularities automatically in rational problem; consider solving the problem as general");
    } else {
      /* AAA algorithm */
      PetscCall(NEPNLEIGSAAASingularities(nep,ndpt,ds,&ndptx,dxi));
    }
  }
  /* Look for Leja-Bagby points in the discretization sets */
  s[0]  = ds[0];
  xi[0] = (ndptx>0)?dxi[0]:PETSC_INFINITY;
  PetscCheck(PetscAbsScalar(xi[0])>=10*PETSC_MACHINE_EPSILON,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"Singularity point 0 is nearly zero: %g; consider removing the singularity or shifting the problem",(double)PetscAbsScalar(xi[0]));
  beta[0] = 1.0; /* scaling factors are also computed here */
  for (i=0;i<ndpt;i++) {
    nrs[i] = 1.0;
    nrxi[i] = 1.0;
  }
  for (k=1;k<ctx->ddmaxit;k++) {
    maxnrs = 0.0;
    minnrxi = PETSC_MAX_REAL;
    for (i=0;i<ndpt;i++) {
      nrs[i] *= ((ds[i]-s[k-1])/(1.0-ds[i]/xi[k-1]))/beta[k-1];
      if (PetscAbsScalar(nrs[i])>maxnrs) {maxnrs = PetscAbsScalar(nrs[i]); s[k] = ds[i];}
    }
    if (ndptx>k) {
      for (i=1;i<ndptx;i++) {
        nrxi[i] *= ((dxi[i]-s[k-1])/(1.0-dxi[i]/xi[k-1]))/beta[k-1];
        if (PetscAbsScalar(nrxi[i])<minnrxi) {minnrxi = PetscAbsScalar(nrxi[i]); xi[k] = dxi[i];}
      }
      PetscCheck(PetscAbsScalar(xi[k])>=10*PETSC_MACHINE_EPSILON,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"Singularity point %" PetscInt_FMT " is nearly zero: %g; consider removing the singularity or shifting the problem",k,(double)PetscAbsScalar(xi[k]));
    } else xi[k] = PETSC_INFINITY;
    beta[k] = maxnrs;
  }
  PetscCall(PetscFree5(ds,dsi,dxi,nrs,nrxi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPNLEIGSEvalNRTFunct(NEP nep,PetscInt k,PetscScalar sigma,PetscScalar *b)
{
  NEP_NLEIGS  *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt    i;
  PetscScalar *beta=ctx->beta,*s=ctx->s,*xi=ctx->xi;

  PetscFunctionBegin;
  b[0] = 1.0/beta[0];
  for (i=0;i<k;i++) {
    b[i+1] = ((sigma-s[i])*b[i])/(beta[i+1]*(1.0-sigma/xi[i]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_Fun(Mat A,Vec x,Vec y)
{
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            i;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(MatMult(ctx->A[0],x,y));
  if (ctx->coeff[0]!=1.0) PetscCall(VecScale(y,ctx->coeff[0]));
  for (i=1;i<ctx->nmat;i++) {
    PetscCall(MatMult(ctx->A[i],x,ctx->t));
    PetscCall(VecAXPY(y,ctx->coeff[i],ctx->t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultTranspose_Fun(Mat A,Vec x,Vec y)
{
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            i;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(MatMultTranspose(ctx->A[0],x,y));
  if (ctx->coeff[0]!=1.0) PetscCall(VecScale(y,ctx->coeff[0]));
  for (i=1;i<ctx->nmat;i++) {
    PetscCall(MatMultTranspose(ctx->A[i],x,ctx->t));
    PetscCall(VecAXPY(y,ctx->coeff[i],ctx->t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatGetDiagonal_Fun(Mat A,Vec diag)
{
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            i;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(MatGetDiagonal(ctx->A[0],diag));
  if (ctx->coeff[0]!=1.0) PetscCall(VecScale(diag,ctx->coeff[0]));
  for (i=1;i<ctx->nmat;i++) {
    PetscCall(MatGetDiagonal(ctx->A[i],ctx->t));
    PetscCall(VecAXPY(diag,ctx->coeff[i],ctx->t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDuplicate_Fun(Mat A,MatDuplicateOption op,Mat *B)
{
  PetscInt            m,n,M,N,i;
  NEP_NLEIGS_MATSHELL *ctxnew,*ctx;
  void                (*fun)(void);

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(PetscNew(&ctxnew));
  ctxnew->nmat = ctx->nmat;
  ctxnew->maxnmat = ctx->maxnmat;
  PetscCall(PetscMalloc2(ctxnew->maxnmat,&ctxnew->A,ctxnew->maxnmat,&ctxnew->coeff));
  for (i=0;i<ctx->nmat;i++) {
    PetscCall(PetscObjectReference((PetscObject)ctx->A[i]));
    ctxnew->A[i] = ctx->A[i];
    ctxnew->coeff[i] = ctx->coeff[i];
  }
  PetscCall(MatGetSize(ctx->A[0],&M,&N));
  PetscCall(MatGetLocalSize(ctx->A[0],&m,&n));
  PetscCall(VecDuplicate(ctx->t,&ctxnew->t));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)A),m,n,M,N,(void*)ctxnew,B));
  PetscCall(MatShellSetManageScalingShifts(*B));
  PetscCall(MatShellGetOperation(A,MATOP_MULT,&fun));
  PetscCall(MatShellSetOperation(*B,MATOP_MULT,fun));
  PetscCall(MatShellGetOperation(A,MATOP_MULT_TRANSPOSE,&fun));
  PetscCall(MatShellSetOperation(*B,MATOP_MULT_TRANSPOSE,fun));
  PetscCall(MatShellGetOperation(A,MATOP_GET_DIAGONAL,&fun));
  PetscCall(MatShellSetOperation(*B,MATOP_GET_DIAGONAL,fun));
  PetscCall(MatShellGetOperation(A,MATOP_DUPLICATE,&fun));
  PetscCall(MatShellSetOperation(*B,MATOP_DUPLICATE,fun));
  PetscCall(MatShellGetOperation(A,MATOP_DESTROY,&fun));
  PetscCall(MatShellSetOperation(*B,MATOP_DESTROY,fun));
  PetscCall(MatShellGetOperation(A,MATOP_AXPY,&fun));
  PetscCall(MatShellSetOperation(*B,MATOP_AXPY,fun));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_Fun(Mat A)
{
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            i;

  PetscFunctionBeginUser;
  if (A) {
    PetscCall(MatShellGetContext(A,&ctx));
    for (i=0;i<ctx->nmat;i++) PetscCall(MatDestroy(&ctx->A[i]));
    PetscCall(VecDestroy(&ctx->t));
    PetscCall(PetscFree2(ctx->A,ctx->coeff));
    PetscCall(PetscFree(ctx));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatAXPY_Fun(Mat Y,PetscScalar a,Mat X,MatStructure str)
{
  NEP_NLEIGS_MATSHELL *ctxY,*ctxX;
  PetscInt            i,j;
  PetscBool           found;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(Y,&ctxY));
  PetscCall(MatShellGetContext(X,&ctxX));
  for (i=0;i<ctxX->nmat;i++) {
    found = PETSC_FALSE;
    for (j=0;!found&&j<ctxY->nmat;j++) {
      if (ctxX->A[i]==ctxY->A[j]) {
        found = PETSC_TRUE;
        ctxY->coeff[j] += a*ctxX->coeff[i];
      }
    }
    if (!found) {
      ctxY->coeff[ctxY->nmat] = a*ctxX->coeff[i];
      ctxY->A[ctxY->nmat++] = ctxX->A[i];
      PetscCall(PetscObjectReference((PetscObject)ctxX->A[i]));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatScale_Fun(Mat M,PetscScalar a)
{
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            i;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(M,&ctx));
  for (i=0;i<ctx->nmat;i++) ctx->coeff[i] *= a;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NLEIGSMatToMatShellArray(Mat A,Mat *Ms,PetscInt maxnmat)
{
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            m,n,M,N;
  PetscBool           has;

  PetscFunctionBegin;
  PetscCall(MatHasOperation(A,MATOP_DUPLICATE,&has));
  PetscCheck(has,PetscObjectComm((PetscObject)A),PETSC_ERR_USER,"MatDuplicate operation required");
  PetscCall(PetscNew(&ctx));
  ctx->maxnmat = maxnmat;
  PetscCall(PetscMalloc2(ctx->maxnmat,&ctx->A,ctx->maxnmat,&ctx->coeff));
  PetscCall(MatDuplicate(A,MAT_COPY_VALUES,&ctx->A[0]));
  ctx->nmat = 1;
  ctx->coeff[0] = 1.0;
  PetscCall(MatCreateVecs(A,&ctx->t,NULL));
  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)A),m,n,M,N,(void*)ctx,Ms));
  PetscCall(MatShellSetManageScalingShifts(*Ms));
  PetscCall(MatShellSetOperation(*Ms,MATOP_MULT,(void(*)(void))MatMult_Fun));
  PetscCall(MatShellSetOperation(*Ms,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTranspose_Fun));
  PetscCall(MatShellSetOperation(*Ms,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Fun));
  PetscCall(MatShellSetOperation(*Ms,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_Fun));
  PetscCall(MatShellSetOperation(*Ms,MATOP_DESTROY,(void(*)(void))MatDestroy_Fun));
  PetscCall(MatShellSetOperation(*Ms,MATOP_AXPY,(void(*)(void))MatAXPY_Fun));
  PetscCall(MatShellSetOperation(*Ms,MATOP_SCALE,(void(*)(void))MatScale_Fun));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   MatIsShellAny - returns true if any of the n matrices is a shell matrix
 */
static PetscErrorCode MatIsShellAny(Mat *A,PetscInt n,PetscBool *shell)
{
  PetscInt       i;
  PetscBool      flg;

  PetscFunctionBegin;
  *shell = PETSC_FALSE;
  for (i=0;i<n;i++) {
    PetscCall(MatIsShell(A[i],&flg));
    if (flg) { *shell = PETSC_TRUE; break; }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSDividedDifferences_split(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       k,j,i,maxnmat,nmax;
  PetscReal      norm0,norm,*matnorm;
  PetscScalar    *s=ctx->s,*beta=ctx->beta,*xi=ctx->xi,*b,alpha,*coeffs,*pK,*pH,sone=1.0;
  Mat            T,P,Ts,K,H;
  PetscBool      shell,hasmnorm=PETSC_FALSE,matrix=PETSC_TRUE;
  PetscBLASInt   n_;

  PetscFunctionBegin;
  nmax = ctx->ddmaxit;
  PetscCall(PetscMalloc1(nep->nt*nmax,&ctx->coeffD));
  PetscCall(PetscMalloc3(nmax+1,&b,nmax+1,&coeffs,nep->nt,&matnorm));
  for (j=0;j<nep->nt;j++) {
    PetscCall(MatHasOperation(nep->A[j],MATOP_NORM,&hasmnorm));
    if (!hasmnorm) break;
    PetscCall(MatNorm(nep->A[j],NORM_INFINITY,matnorm+j));
  }
  /* Try matrix functions scheme */
  PetscCall(PetscCalloc2(nmax*nmax,&pK,nmax*nmax,&pH));
  for (i=0;i<nmax-1;i++) {
    pK[(nmax+1)*i]   = 1.0;
    pK[(nmax+1)*i+1] = beta[i+1]/xi[i];
    pH[(nmax+1)*i]   = s[i];
    pH[(nmax+1)*i+1] = beta[i+1];
  }
  pH[nmax*nmax-1] = s[nmax-1];
  pK[nmax*nmax-1] = 1.0;
  PetscCall(PetscBLASIntCast(nmax,&n_));
  PetscCallBLAS("BLAStrsm",BLAStrsm_("R","L","N","U",&n_,&n_,&sone,pK,&n_,pH,&n_));
  /* The matrix to be used is in H. K will be a work-space matrix */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,nmax,nmax,pH,&H));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,nmax,nmax,pK,&K));
  for (j=0;matrix&&j<nep->nt;j++) {
    PetscCall(PetscPushErrorHandler(PetscReturnErrorHandler,NULL));
    ierr = FNEvaluateFunctionMat(nep->f[j],H,K);
    PetscCall(PetscPopErrorHandler());
    if (!ierr) {
      for (i=0;i<nmax;i++) ctx->coeffD[j+i*nep->nt] = pK[i]*beta[0];
    } else {
      matrix = PETSC_FALSE;
      PetscCall(PetscFPTrapPop());
    }
  }
  PetscCall(MatDestroy(&H));
  PetscCall(MatDestroy(&K));
  if (!matrix) {
    for (j=0;j<nep->nt;j++) {
      PetscCall(FNEvaluateFunction(nep->f[j],s[0],ctx->coeffD+j));
      ctx->coeffD[j] *= beta[0];
    }
  }
  if (hasmnorm) {
    norm0 = 0.0;
    for (j=0;j<nep->nt;j++) norm0 += matnorm[j]*PetscAbsScalar(ctx->coeffD[j]);
  } else {
    norm0 = 0.0;
    for (j=0;j<nep->nt;j++) norm0 = PetscMax(PetscAbsScalar(ctx->coeffD[j]),norm0);
  }
  ctx->nmat = ctx->ddmaxit;
  for (k=1;k<ctx->ddmaxit;k++) {
    if (!matrix) {
      PetscCall(NEPNLEIGSEvalNRTFunct(nep,k,s[k],b));
      for (i=0;i<nep->nt;i++) {
        PetscCall(FNEvaluateFunction(nep->f[i],s[k],ctx->coeffD+k*nep->nt+i));
        for (j=0;j<k;j++) {
          ctx->coeffD[k*nep->nt+i] -= b[j]*ctx->coeffD[i+nep->nt*j];
        }
        ctx->coeffD[k*nep->nt+i] /= b[k];
      }
    }
    if (hasmnorm) {
      norm = 0.0;
      for (j=0;j<nep->nt;j++) norm += matnorm[j]*PetscAbsScalar(ctx->coeffD[k*nep->nt+j]);
    } else {
      norm = 0.0;
      for (j=0;j<nep->nt;j++) norm = PetscMax(PetscAbsScalar(ctx->coeffD[k*nep->nt+j]),norm);
    }
    if (k>1 && norm/norm0 < ctx->ddtol) {
      ctx->nmat = k+1;
      break;
    }
  }
  if (!ctx->ksp) PetscCall(NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp));
  PetscCall(MatIsShellAny(nep->A,nep->nt,&shell));
  maxnmat = PetscMax(ctx->ddmaxit,nep->nt);
  for (i=0;i<ctx->nshiftsw;i++) {
    PetscCall(NEPNLEIGSEvalNRTFunct(nep,ctx->nmat-1,ctx->shifts[i],coeffs));
    if (!shell) PetscCall(MatDuplicate(nep->A[0],MAT_COPY_VALUES,&T));
    else PetscCall(NLEIGSMatToMatShellArray(nep->A[0],&T,maxnmat));
    if (nep->P) { /* user-defined preconditioner */
      PetscCall(MatDuplicate(nep->P[0],MAT_COPY_VALUES,&P));
    } else P=T;
    alpha = 0.0;
    for (j=0;j<ctx->nmat;j++) alpha += coeffs[j]*ctx->coeffD[j*nep->nt];
    PetscCall(MatScale(T,alpha));
    if (nep->P) PetscCall(MatScale(P,alpha));
    for (k=1;k<nep->nt;k++) {
      alpha = 0.0;
      for (j=0;j<ctx->nmat;j++) alpha += coeffs[j]*ctx->coeffD[j*nep->nt+k];
      if (shell) PetscCall(NLEIGSMatToMatShellArray(nep->A[k],&Ts,maxnmat));
      PetscCall(MatAXPY(T,alpha,shell?Ts:nep->A[k],nep->mstr));
      if (nep->P) PetscCall(MatAXPY(P,alpha,nep->P[k],nep->mstrp));
      if (shell) PetscCall(MatDestroy(&Ts));
    }
    PetscCall(NEP_KSPSetOperators(ctx->ksp[i],T,P));
    PetscCall(KSPSetUp(ctx->ksp[i]));
    PetscCall(MatDestroy(&T));
    if (nep->P) PetscCall(MatDestroy(&P));
  }
  PetscCall(PetscFree3(b,coeffs,matnorm));
  PetscCall(PetscFree2(pK,pH));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSDividedDifferences_callback(NEP nep)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       k,j,i,maxnmat;
  PetscReal      norm0,norm;
  PetscScalar    *s=ctx->s,*beta=ctx->beta,*b,*coeffs;
  Mat            *D=ctx->D,*DP,T,P;
  PetscBool      shell,has,vec=PETSC_FALSE,precond=(nep->function_pre!=nep->function)?PETSC_TRUE:PETSC_FALSE;
  PetscRandom    rand=NULL;
  Vec            w[2];

  PetscFunctionBegin;
  PetscCall(PetscMalloc2(ctx->ddmaxit+1,&b,ctx->ddmaxit+1,&coeffs));
  if (nep->V) PetscCall(BVGetRandomContext(nep->V,&rand));
  T = nep->function;
  P = nep->function_pre;
  PetscCall(NEPComputeFunction(nep,s[0],T,P));
  PetscCall(MatIsShell(T,&shell));
  maxnmat = PetscMax(ctx->ddmaxit,nep->nt);
  if (!shell) PetscCall(MatDuplicate(T,MAT_COPY_VALUES,&D[0]));
  else PetscCall(NLEIGSMatToMatShellArray(T,&D[0],maxnmat));
  if (beta[0]!=1.0) PetscCall(MatScale(D[0],1.0/beta[0]));
  PetscCall(MatHasOperation(D[0],MATOP_NORM,&has));
  if (has) PetscCall(MatNorm(D[0],NORM_FROBENIUS,&norm0));
  else {
    PetscCall(MatCreateVecs(D[0],NULL,&w[0]));
    PetscCall(VecDuplicate(w[0],&w[1]));
    PetscCall(VecDuplicate(w[0],&ctx->vrn));
    PetscCall(VecSetRandomNormal(ctx->vrn,rand,w[0],w[1]));
    PetscCall(VecNormalize(ctx->vrn,NULL));
    vec = PETSC_TRUE;
    PetscCall(MatNormEstimate(D[0],ctx->vrn,w[0],&norm0));
  }
  if (precond) {
    PetscCall(PetscMalloc1(ctx->ddmaxit,&DP));
    PetscCall(MatDuplicate(P,MAT_COPY_VALUES,&DP[0]));
  }
  ctx->nmat = ctx->ddmaxit;
  for (k=1;k<ctx->ddmaxit;k++) {
    PetscCall(NEPNLEIGSEvalNRTFunct(nep,k,s[k],b));
    PetscCall(NEPComputeFunction(nep,s[k],T,P));
    if (!shell) PetscCall(MatDuplicate(T,MAT_COPY_VALUES,&D[k]));
    else PetscCall(NLEIGSMatToMatShellArray(T,&D[k],maxnmat));
    for (j=0;j<k;j++) PetscCall(MatAXPY(D[k],-b[j],D[j],nep->mstr));
    PetscCall(MatScale(D[k],1.0/b[k]));
    PetscCall(MatHasOperation(D[k],MATOP_NORM,&has));
    if (has) PetscCall(MatNorm(D[k],NORM_FROBENIUS,&norm));
    else {
      if (!vec) {
        PetscCall(MatCreateVecs(D[k],NULL,&w[0]));
        PetscCall(VecDuplicate(w[0],&w[1]));
        PetscCall(VecDuplicate(w[0],&ctx->vrn));
        PetscCall(VecSetRandomNormal(ctx->vrn,rand,w[0],w[1]));
        PetscCall(VecNormalize(ctx->vrn,NULL));
        vec = PETSC_TRUE;
      }
      PetscCall(MatNormEstimate(D[k],ctx->vrn,w[0],&norm));
    }
    if (precond) {
      PetscCall(MatDuplicate(P,MAT_COPY_VALUES,&DP[k]));
      for (j=0;j<k;j++) PetscCall(MatAXPY(DP[k],-b[j],DP[j],nep->mstrp));
      PetscCall(MatScale(DP[k],1.0/b[k]));
    }
    if (k>1 && norm/norm0 < ctx->ddtol && k>1) {
      ctx->nmat = k+1;
      break;
    }
  }
  if (!ctx->ksp) PetscCall(NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp));
  for (i=0;i<ctx->nshiftsw;i++) {
    PetscCall(NEPNLEIGSEvalNRTFunct(nep,ctx->nmat-1,ctx->shifts[i],coeffs));
    PetscCall(MatDuplicate(D[0],MAT_COPY_VALUES,&T));
    if (coeffs[0]!=1.0) PetscCall(MatScale(T,coeffs[0]));
    for (j=1;j<ctx->nmat;j++) PetscCall(MatAXPY(T,coeffs[j],D[j],nep->mstr));
    if (precond) {
      PetscCall(MatDuplicate(DP[0],MAT_COPY_VALUES,&P));
      if (coeffs[0]!=1.0) PetscCall(MatScale(P,coeffs[0]));
      for (j=1;j<ctx->nmat;j++) PetscCall(MatAXPY(P,coeffs[j],DP[j],nep->mstrp));
    } else P=T;
    PetscCall(NEP_KSPSetOperators(ctx->ksp[i],T,P));
    PetscCall(KSPSetUp(ctx->ksp[i]));
    PetscCall(MatDestroy(&T));
  }
  PetscCall(PetscFree2(b,coeffs));
  if (vec) {
    PetscCall(VecDestroy(&w[0]));
    PetscCall(VecDestroy(&w[1]));
  }
  if (precond) {
    PetscCall(MatDestroy(&P));
    PetscCall(MatDestroyMatrices(ctx->nmat,&DP));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   NEPKrylovConvergence - This is the analogue to EPSKrylovConvergence.
*/
static PetscErrorCode NEPNLEIGSKrylovConvergence(NEP nep,PetscBool getall,PetscInt kini,PetscInt nits,PetscReal betah,PetscScalar betak,PetscInt *kout,Vec *w)
{
  PetscInt       k,newk,marker,inside;
  PetscScalar    re,im;
  PetscReal      resnorm,tt;
  PetscBool      istrivial;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  PetscCall(RGIsTrivial(nep->rg,&istrivial));
  marker = -1;
  if (nep->trackall) getall = PETSC_TRUE;
  for (k=kini;k<kini+nits;k++) {
    /* eigenvalue */
    re = nep->eigr[k];
    im = nep->eigi[k];
    if (!istrivial) {
      if (!ctx->nshifts) PetscCall(NEPNLEIGSBackTransform((PetscObject)nep,1,&re,&im));
      PetscCall(RGCheckInside(nep->rg,1,&re,&im,&inside));
      if (marker==-1 && inside<0) marker = k;
    }
    newk = k;
    PetscCall(DSVectors(nep->ds,DS_MAT_X,&newk,&resnorm));
    tt = (ctx->nshifts)?SlepcAbsEigenvalue(betak-nep->eigr[k]*betah,nep->eigi[k]*betah):betah;
    resnorm *=  PetscAbsReal(tt);
    /* error estimate */
    PetscCall((*nep->converged)(nep,nep->eigr[k],nep->eigi[k],resnorm,&nep->errest[k],nep->convergedctx));
    if (marker==-1 && nep->errest[k] >= nep->tol) marker = k;
    if (newk==k+1) {
      nep->errest[k+1] = nep->errest[k];
      k++;
    }
    if (marker!=-1 && !getall) break;
  }
  if (marker!=-1) k = marker;
  *kout = k;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSetUp_NLEIGS(NEP nep)
{
  PetscInt       k,in;
  PetscScalar    zero=0.0;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  SlepcSC        sc;
  PetscBool      istrivial;

  PetscFunctionBegin;
  PetscCall(NEPSetDimensions_Default(nep,nep->nev,&nep->ncv,&nep->mpd));
  PetscCheck(nep->ncv<=nep->nev+nep->mpd,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (nep->max_it==PETSC_DETERMINE) nep->max_it = PetscMax(5000,2*nep->n/nep->ncv);
  if (!ctx->ddmaxit) ctx->ddmaxit = LBPOINTS;
  PetscCall(RGIsTrivial(nep->rg,&istrivial));
  PetscCheck(!istrivial,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NEPNLEIGS requires a nontrivial region defining the target set");
  if (!nep->which) nep->which = NEP_TARGET_MAGNITUDE;
  PetscCheck(nep->which==NEP_TARGET_MAGNITUDE || nep->which==NEP_TARGET_REAL || nep->which==NEP_TARGET_IMAGINARY || nep->which==NEP_WHICH_USER,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver supports only target selection of eigenvalues");

  /* Initialize the NLEIGS context structure */
  k = ctx->ddmaxit;
  PetscCall(PetscMalloc4(k,&ctx->s,k,&ctx->xi,k,&ctx->beta,k,&ctx->D));
  nep->data = ctx;
  if (nep->tol==(PetscReal)PETSC_DETERMINE) nep->tol = SLEPC_DEFAULT_TOL;
  if (ctx->ddtol==(PetscReal)PETSC_DETERMINE) ctx->ddtol = nep->tol/10.0;
  if (!ctx->keep) ctx->keep = 0.5;

  /* Compute Leja-Bagby points and scaling values */
  PetscCall(NEPNLEIGSLejaBagbyPoints(nep));
  if (nep->problem_type!=NEP_RATIONAL) {
    PetscCall(RGCheckInside(nep->rg,1,&nep->target,&zero,&in));
    PetscCheck(in>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The target is not inside the target set");
  }

  /* Compute the divided difference matrices */
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) PetscCall(NEPNLEIGSDividedDifferences_split(nep));
  else PetscCall(NEPNLEIGSDividedDifferences_callback(nep));
  PetscCall(NEPAllocateSolution(nep,ctx->nmat-1));
  PetscCall(NEPSetWorkVecs(nep,4));
  if (!ctx->fullbasis) {
    PetscCheck(!nep->twosided,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Two-sided variant requires the full-basis option, rerun with -nep_nleigs_full_basis");
    /* set-up DS and transfer split operator functions */
    PetscCall(DSSetType(nep->ds,ctx->nshifts?DSGNHEP:DSNHEP));
    PetscCall(DSAllocate(nep->ds,nep->ncv+1));
    PetscCall(DSGetSlepcSC(nep->ds,&sc));
    if (!ctx->nshifts) sc->map = NEPNLEIGSBackTransform;
    PetscCall(DSSetExtraRow(nep->ds,PETSC_TRUE));
    sc->mapobj        = (PetscObject)nep;
    sc->rg            = nep->rg;
    sc->comparison    = nep->sc->comparison;
    sc->comparisonctx = nep->sc->comparisonctx;
    PetscCall(BVDestroy(&ctx->V));
    PetscCall(BVCreateTensor(nep->V,ctx->nmat-1,&ctx->V));
    nep->ops->solve          = NEPSolve_NLEIGS;
    nep->ops->computevectors = NEPComputeVectors_Schur;
  } else {
    PetscCall(NEPSetUp_NLEIGS_FullBasis(nep));
    nep->ops->solve          = NEPSolve_NLEIGS_FullBasis;
    nep->ops->computevectors = NULL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Extend the TOAR basis by applying the matrix operator
  over a vector which is decomposed on the TOAR way
  Input:
    - S,V: define the latest Arnoldi vector (nv vectors in V)
  Output:
    - t: new vector extending the TOAR basis
    - r: temporally coefficients to compute the TOAR coefficients
         for the new Arnoldi vector
  Workspace: t_ (two vectors)
*/
static PetscErrorCode NEPTOARExtendBasis(NEP nep,PetscInt idxrktg,PetscScalar *S,PetscInt ls,PetscInt nv,BV W,BV V,Vec t,PetscScalar *r,PetscInt lr,Vec *t_)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       deg=ctx->nmat-1,k,j;
  Vec            v=t_[0],q=t_[1],w;
  PetscScalar    *beta=ctx->beta,*s=ctx->s,*xi=ctx->xi,*coeffs,sigma;

  PetscFunctionBegin;
  if (!ctx->ksp) PetscCall(NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp));
  sigma = ctx->shifts[idxrktg];
  PetscCall(BVSetActiveColumns(nep->V,0,nv));
  PetscCall(PetscMalloc1(ctx->nmat,&coeffs));
  PetscCheck(PetscAbsScalar(s[deg-2]-sigma)>100*PETSC_MACHINE_EPSILON,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"Breakdown in NLEIGS");
  /* i-part stored in (i-1) position */
  for (j=0;j<nv;j++) {
    r[(deg-2)*lr+j] = (S[(deg-2)*ls+j]+(beta[deg-1]/xi[deg-2])*S[(deg-1)*ls+j])/(s[deg-2]-sigma);
  }
  PetscCall(BVSetActiveColumns(W,0,deg));
  PetscCall(BVGetColumn(W,deg-1,&w));
  PetscCall(BVMultVec(V,1.0/beta[deg],0,w,S+(deg-1)*ls));
  PetscCall(BVRestoreColumn(W,deg-1,&w));
  PetscCall(BVGetColumn(W,deg-2,&w));
  PetscCall(BVMultVec(V,1.0,0.0,w,r+(deg-2)*lr));
  PetscCall(BVRestoreColumn(W,deg-2,&w));
  for (k=deg-2;k>0;k--) {
    PetscCheck(PetscAbsScalar(s[k-1]-sigma)>100*PETSC_MACHINE_EPSILON,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"Breakdown in NLEIGS");
    for (j=0;j<nv;j++) r[(k-1)*lr+j] = (S[(k-1)*ls+j]+(beta[k]/xi[k-1])*S[k*ls+j]-beta[k]*(1.0-sigma/xi[k-1])*r[(k)*lr+j])/(s[k-1]-sigma);
    PetscCall(BVGetColumn(W,k-1,&w));
    PetscCall(BVMultVec(V,1.0,0.0,w,r+(k-1)*lr));
    PetscCall(BVRestoreColumn(W,k-1,&w));
  }
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    for (j=0;j<ctx->nmat-2;j++) coeffs[j] = ctx->coeffD[nep->nt*j];
    coeffs[ctx->nmat-2] = ctx->coeffD[nep->nt*(ctx->nmat-1)];
    PetscCall(BVMultVec(W,1.0,0.0,v,coeffs));
    PetscCall(MatMult(nep->A[0],v,q));
    for (k=1;k<nep->nt;k++) {
      for (j=0;j<ctx->nmat-2;j++) coeffs[j] = ctx->coeffD[nep->nt*j+k];
      coeffs[ctx->nmat-2] = ctx->coeffD[nep->nt*(ctx->nmat-1)+k];
      PetscCall(BVMultVec(W,1.0,0,v,coeffs));
      PetscCall(MatMult(nep->A[k],v,t));
      PetscCall(VecAXPY(q,1.0,t));
    }
    PetscCall(KSPSolve(ctx->ksp[idxrktg],q,t));
    PetscCall(VecScale(t,-1.0));
  } else {
    for (k=0;k<deg-1;k++) {
      PetscCall(BVGetColumn(W,k,&w));
      PetscCall(MatMult(ctx->D[k],w,q));
      PetscCall(BVRestoreColumn(W,k,&w));
      PetscCall(BVInsertVec(W,k,q));
    }
    PetscCall(BVGetColumn(W,deg-1,&w));
    PetscCall(MatMult(ctx->D[deg],w,q));
    PetscCall(BVRestoreColumn(W,k,&w));
    PetscCall(BVInsertVec(W,k,q));
    for (j=0;j<ctx->nmat-1;j++) coeffs[j] = 1.0;
    PetscCall(BVMultVec(W,1.0,0.0,q,coeffs));
    PetscCall(KSPSolve(ctx->ksp[idxrktg],q,t));
    PetscCall(VecScale(t,-1.0));
  }
  PetscCall(PetscFree(coeffs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute TOAR coefficients of the blocks of the new Arnoldi vector computed
*/
static PetscErrorCode NEPTOARCoefficients(NEP nep,PetscScalar sigma,PetscInt nv,PetscScalar *S,PetscInt ls,PetscScalar *r,PetscInt lr,PetscScalar *x,PetscScalar *work)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       k,j,d=ctx->nmat-1;
  PetscScalar    *t=work;

  PetscFunctionBegin;
  PetscCall(NEPNLEIGSEvalNRTFunct(nep,d-1,sigma,t));
  for (k=0;k<d-1;k++) {
    for (j=0;j<=nv;j++) r[k*lr+j] += t[k]*x[j];
  }
  for (j=0;j<=nv;j++) r[(d-1)*lr+j] = t[d-1]*x[j];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute continuation vector coefficients for the Rational-Krylov run.
  dim(work) >= (end-ini)*(end-ini+1) + end+1 + 2*(end-ini+1), dim(t) = end.
*/
static PetscErrorCode NEPNLEIGS_RKcontinuation(NEP nep,PetscInt ini,PetscInt end,PetscScalar *K,PetscScalar *H,PetscInt ld,PetscScalar sigma,PetscScalar *S,PetscInt lds,PetscScalar *cont,PetscScalar *t,PetscScalar *work)
{
  PetscScalar    *x,*W,*tau,sone=1.0,szero=0.0;
  PetscInt       i,j,n1,n,nwu=0;
  PetscBLASInt   info,n_,n1_,one=1,dim,lds_;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (!ctx->nshifts || !end) {
    t[0] = 1;
    PetscCall(PetscArraycpy(cont,S+end*lds,lds));
  } else {
    n   = end-ini;
    n1  = n+1;
    x   = work+nwu;
    nwu += end+1;
    tau = work+nwu;
    nwu += n;
    W   = work+nwu;
    nwu += n1*n;
    for (j=ini;j<end;j++) {
      for (i=ini;i<=end;i++) W[(j-ini)*n1+i-ini] = K[j*ld+i] -H[j*ld+i]*sigma;
    }
    PetscCall(PetscBLASIntCast(n,&n_));
    PetscCall(PetscBLASIntCast(n1,&n1_));
    PetscCall(PetscBLASIntCast(end+1,&dim));
    PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
    PetscCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&n1_,&n_,W,&n1_,tau,work+nwu,&n1_,&info));
    SlepcCheckLapackInfo("geqrf",info);
    for (i=0;i<end;i++) t[i] = 0.0;
    t[end] = 1.0;
    for (j=n-1;j>=0;j--) {
      for (i=0;i<ini+j;i++) x[i] = 0.0;
      x[ini+j] = 1.0;
      for (i=j+1;i<n1;i++) x[i+ini] = W[i+n1*j];
      tau[j] = PetscConj(tau[j]);
      PetscCallBLAS("LAPACKlarf",LAPACKlarf_("L",&dim,&one,x,&one,tau+j,t,&dim,work+nwu));
    }
    PetscCall(PetscBLASIntCast(lds,&lds_));
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&n1_,&sone,S,&lds_,t,&one,&szero,cont,&one));
    PetscCall(PetscFPTrapPop());
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute a run of Arnoldi iterations
*/
static PetscErrorCode NEPNLEIGSTOARrun(NEP nep,Mat MK,Mat MH,BV W,PetscInt k,PetscInt *M,PetscReal *betah,PetscScalar *betak,PetscBool *breakdown,Vec *t_)
{
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt       i,j,m=*M,lwa,deg=ctx->nmat-1,lds,nqt,ld,l,ldh;
  Vec            t;
  PetscReal      norm=0.0;
  PetscScalar    *x,*work,*tt,sigma=1.0,*cont,*S,*K=NULL,*H;
  PetscBool      lindep;
  Mat            MS;

  PetscFunctionBegin;
  *betah = 0.0; *betak = 0.0;
  PetscCall(MatDenseGetArray(MH,&H));
  if (MK) PetscCall(MatDenseGetArray(MK,&K));
  PetscCall(MatDenseGetLDA(MH,&ldh));
  PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
  PetscCall(MatDenseGetArray(MS,&S));
  PetscCall(BVGetSizes(nep->V,NULL,NULL,&ld));
  lds = ld*deg;
  PetscCall(BVGetActiveColumns(nep->V,&l,&nqt));
  lwa = PetscMax(ld,deg)+(m+1)*(m+1)+4*(m+1);
  PetscCall(PetscMalloc4(ld,&x,lwa,&work,m+1,&tt,lds,&cont));
  PetscCall(BVSetActiveColumns(ctx->V,0,m));
  for (j=k;j<m;j++) {
    sigma = ctx->shifts[(++(ctx->idxrk))%ctx->nshiftsw];

    /* Continuation vector */
    PetscCall(NEPNLEIGS_RKcontinuation(nep,0,j,K,H,ldh,sigma,S,lds,cont,tt,work));

    /* apply operator */
    PetscCall(BVGetColumn(nep->V,nqt,&t));
    PetscCall(NEPTOARExtendBasis(nep,(ctx->idxrk)%ctx->nshiftsw,cont,ld,nqt,W,nep->V,t,S+(j+1)*lds,ld,t_));
    PetscCall(BVRestoreColumn(nep->V,nqt,&t));

    /* orthogonalize */
    PetscCall(BVOrthogonalizeColumn(nep->V,nqt,x,&norm,&lindep));
    if (!lindep) {
      x[nqt] = norm;
      PetscCall(BVScaleColumn(nep->V,nqt,1.0/norm));
      nqt++;
    } else x[nqt] = 0.0;

    PetscCall(NEPTOARCoefficients(nep,sigma,nqt-1,cont,ld,S+(j+1)*lds,ld,x,work));

    /* Level-2 orthogonalization */
    PetscCall(BVOrthogonalizeColumn(ctx->V,j+1,H+j*ldh,&norm,breakdown));
    H[j+1+ldh*j] = norm;
    if (ctx->nshifts && MK) {
      for (i=0;i<=j;i++) K[i+ldh*j] = sigma*H[i+ldh*j] + tt[i];
      K[j+1+ldh*j] = sigma*H[j+1+ldh*j];
    }
    if (*breakdown) {
      *M = j+1;
      break;
    }
    PetscCall(BVScaleColumn(ctx->V,j+1,1.0/norm));
    PetscCall(BVSetActiveColumns(nep->V,l,nqt));
  }
  *betah = norm;
  if (ctx->nshifts) *betak = norm*sigma;
  PetscCall(PetscFree4(x,work,tt,cont));
  PetscCall(MatDenseRestoreArray(MS,&S));
  PetscCall(MatDenseRestoreArray(MH,&H));
  if (MK) PetscCall(MatDenseRestoreArray(MK,&K));
  PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPSolve_NLEIGS(NEP nep)
{
  NEP_NLEIGS        *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt          i,k=0,l,nv=0,ld,lds,nq;
  PetscInt          deg=ctx->nmat-1,nconv=0,dsn,dsk;
  PetscScalar       *pU,betak=0,*eigr,*eigi;
  const PetscScalar *S;
  PetscReal         betah;
  PetscBool         falselock=PETSC_FALSE,breakdown=PETSC_FALSE;
  BV                W;
  Mat               H,K=NULL,MS,MQ,U;

  PetscFunctionBegin;
  if (ctx->lock) {
    /* undocumented option to use a cheaper locking instead of the true locking */
    PetscCall(PetscOptionsGetBool(NULL,NULL,"-nep_nleigs_falselocking",&falselock,NULL));
  }

  PetscCall(BVGetSizes(nep->V,NULL,NULL,&ld));
  lds = deg*ld;
  if (!ctx->nshifts) PetscCall(PetscMalloc2(nep->ncv,&eigr,nep->ncv,&eigi));
  else { eigr = nep->eigr; eigi = nep->eigi; }
  PetscCall(BVDuplicateResize(nep->V,PetscMax(nep->nt-1,ctx->nmat-1),&W));

  /* clean projected matrix (including the extra-arrow) */
  PetscCall(DSSetDimensions(nep->ds,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(DSGetMat(nep->ds,DS_MAT_A,&H));
  PetscCall(MatZeroEntries(H));
  PetscCall(DSRestoreMat(nep->ds,DS_MAT_A,&H));
  if (ctx->nshifts) {
    PetscCall(DSGetMat(nep->ds,DS_MAT_B,&H));
    PetscCall(MatZeroEntries(H));
    PetscCall(DSRestoreMat(nep->ds,DS_MAT_B,&H));
  }

  /* Get the starting Arnoldi vector */
  PetscCall(BVTensorBuildFirstColumn(ctx->V,nep->nini));

  /* Restart loop */
  l = 0;
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;

    /* Compute an nv-step Krylov relation */
    nv = PetscMin(nep->nconv+nep->mpd,nep->ncv);
    if (ctx->nshifts) PetscCall(DSGetMat(nep->ds,DS_MAT_A,&K));
    PetscCall(DSGetMat(nep->ds,ctx->nshifts?DS_MAT_B:DS_MAT_A,&H));
    PetscCall(NEPNLEIGSTOARrun(nep,K,H,W,nep->nconv+l,&nv,&betah,&betak,&breakdown,nep->work));
    PetscCall(DSRestoreMat(nep->ds,ctx->nshifts?DS_MAT_B:DS_MAT_A,&H));
    if (ctx->nshifts) PetscCall(DSRestoreMat(nep->ds,DS_MAT_A,&K));
    PetscCall(DSSetDimensions(nep->ds,nv,nep->nconv,nep->nconv+l));
    if (l==0) PetscCall(DSSetState(nep->ds,DS_STATE_INTERMEDIATE));
    else PetscCall(DSSetState(nep->ds,DS_STATE_RAW));

    /* Solve projected problem */
    PetscCall(DSSolve(nep->ds,nep->eigr,nep->eigi));
    PetscCall(DSSort(nep->ds,nep->eigr,nep->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(nep->ds));
    PetscCall(DSSynchronize(nep->ds,nep->eigr,nep->eigi));

    /* Check convergence */
    PetscCall(NEPNLEIGSKrylovConvergence(nep,PETSC_FALSE,nep->nconv,nv-nep->nconv,betah,betak,&k,nep->work));
    PetscCall((*nep->stopping)(nep,nep->its,nep->max_it,k,nep->nev,&nep->reason,nep->stoppingctx));

    /* Update l */
    if (nep->reason != NEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
      PetscCall(DSGetTruncateSize(nep->ds,k,nv,&l));
      if (!breakdown) {
        /* Prepare the Rayleigh quotient for restart */
        PetscCall(DSGetDimensions(nep->ds,&dsn,NULL,&dsk,NULL));
        PetscCall(DSSetDimensions(nep->ds,dsn,k,dsk));
        PetscCall(DSTruncate(nep->ds,k+l,PETSC_FALSE));
      }
    }
    nconv = k;
    if (!ctx->lock && nep->reason == NEP_CONVERGED_ITERATING && !breakdown) { l += k; k = 0; }
    if (l) PetscCall(PetscInfo(nep,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    /* Update S */
    PetscCall(DSGetMat(nep->ds,ctx->nshifts?DS_MAT_Z:DS_MAT_Q,&MQ));
    PetscCall(BVMultInPlace(ctx->V,MQ,nep->nconv,k+l));
    PetscCall(DSRestoreMat(nep->ds,ctx->nshifts?DS_MAT_Z:DS_MAT_Q,&MQ));

    /* Copy last column of S */
    PetscCall(BVCopyColumn(ctx->V,nv,k+l));

    if (breakdown && nep->reason == NEP_CONVERGED_ITERATING) {
      /* Stop if breakdown */
      PetscCall(PetscInfo(nep,"Breakdown (it=%" PetscInt_FMT " norm=%g)\n",nep->its,(double)betah));
      nep->reason = NEP_DIVERGED_BREAKDOWN;
    }
    if (nep->reason != NEP_CONVERGED_ITERATING) l--;
    /* truncate S */
    PetscCall(BVGetActiveColumns(nep->V,NULL,&nq));
    if (k+l+deg<=nq) {
      PetscCall(BVSetActiveColumns(ctx->V,nep->nconv,k+l+1));
      if (!falselock && ctx->lock) PetscCall(BVTensorCompress(ctx->V,k-nep->nconv));
      else PetscCall(BVTensorCompress(ctx->V,0));
    }
    nep->nconv = k;
    if (!ctx->nshifts) {
      for (i=0;i<nv;i++) { eigr[i] = nep->eigr[i]; eigi[i] = nep->eigi[i]; }
      PetscCall(NEPNLEIGSBackTransform((PetscObject)nep,nv,eigr,eigi));
    }
    PetscCall(NEPMonitor(nep,nep->its,nconv,eigr,eigi,nep->errest,nv));
  }
  nep->nconv = nconv;
  if (nep->nconv>0) {
    PetscCall(BVSetActiveColumns(ctx->V,0,nep->nconv));
    PetscCall(BVGetActiveColumns(nep->V,NULL,&nq));
    PetscCall(BVSetActiveColumns(nep->V,0,nq));
    if (nq>nep->nconv) {
      PetscCall(BVTensorCompress(ctx->V,nep->nconv));
      PetscCall(BVSetActiveColumns(nep->V,0,nep->nconv));
      nq = nep->nconv;
    }
    if (ctx->nshifts) {
      PetscCall(DSGetMat(nep->ds,DS_MAT_B,&MQ));
      PetscCall(BVMultInPlace(ctx->V,MQ,0,nep->nconv));
      PetscCall(DSRestoreMat(nep->ds,DS_MAT_B,&MQ));
    }
    PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
    PetscCall(MatDenseGetArrayRead(MS,&S));
    PetscCall(PetscMalloc1(nq*nep->nconv,&pU));
    for (i=0;i<nep->nconv;i++) PetscCall(PetscArraycpy(pU+i*nq,S+i*lds,nq));
    PetscCall(MatDenseRestoreArrayRead(MS,&S));
    PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,nq,nep->nconv,pU,&U));
    PetscCall(BVSetActiveColumns(nep->V,0,nq));
    PetscCall(BVMultInPlace(nep->V,U,0,nep->nconv));
    PetscCall(BVSetActiveColumns(nep->V,0,nep->nconv));
    PetscCall(MatDestroy(&U));
    PetscCall(PetscFree(pU));
    PetscCall(DSTruncate(nep->ds,nep->nconv,PETSC_TRUE));
  }

  /* Map eigenvalues back to the original problem */
  if (!ctx->nshifts) {
    PetscCall(NEPNLEIGSBackTransform((PetscObject)nep,nep->nconv,nep->eigr,nep->eigi));
    PetscCall(PetscFree2(eigr,eigi));
  }
  PetscCall(BVDestroy(&W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSSetSingularitiesFunction_NLEIGS(NEP nep,PetscErrorCode (*fun)(NEP,PetscInt*,PetscScalar*,void*),void *ctx)
{
  NEP_NLEIGS *nepctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (fun) nepctx->computesingularities = fun;
  if (ctx) nepctx->singularitiesctx     = ctx;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPNLEIGSSetSingularitiesFunction - Sets a user function to compute a discretization
   of the singularity set (where T(.) is not analytic).

   Logically Collective

   Input Parameters:
+  nep - the NEP context
.  fun - user function (if NULL then NEP retains any previously set value)
-  ctx - [optional] user-defined context for private data for the function
         (may be NULL, in which case NEP retains any previously set value)

   Calling sequence of fun:
$  PetscErrorCode fun(NEP nep,PetscInt *maxnp,PetscScalar *xi,void *ctx)
+   nep   - the NEP context
.   maxnp - on input number of requested points in the discretization (can be set)
.   xi    - computed values of the discretization
-   ctx   - optional context, as set by NEPNLEIGSSetSingularitiesFunction()

   Notes:
   The user-defined function can set a smaller value of maxnp if necessary.
   It is wrong to return a larger value.

   If the problem type has been set to rational with NEPSetProblemType(),
   then it is not necessary to set the singularities explicitly since the
   solver will try to determine them automatically.

   Level: intermediate

.seealso: NEPNLEIGSGetSingularitiesFunction(), NEPSetProblemType()
@*/
PetscErrorCode NEPNLEIGSSetSingularitiesFunction(NEP nep,PetscErrorCode (*fun)(NEP nep,PetscInt *maxnp,PetscScalar *xi,void *ctx),void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscTryMethod(nep,"NEPNLEIGSSetSingularitiesFunction_C",(NEP,PetscErrorCode(*)(NEP,PetscInt*,PetscScalar*,void*),void*),(nep,fun,ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSGetSingularitiesFunction_NLEIGS(NEP nep,PetscErrorCode (**fun)(NEP,PetscInt*,PetscScalar*,void*),void **ctx)
{
  NEP_NLEIGS *nepctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (fun) *fun = nepctx->computesingularities;
  if (ctx) *ctx = nepctx->singularitiesctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPNLEIGSGetSingularitiesFunction - Returns the Function and optionally the user
   provided context for computing a discretization of the singularity set.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  fun - location to put the function (or NULL)
-  ctx - location to stash the function context (or NULL)

   Level: advanced

.seealso: NEPNLEIGSSetSingularitiesFunction()
@*/
PetscErrorCode NEPNLEIGSGetSingularitiesFunction(NEP nep,PetscErrorCode (**fun)(NEP,PetscInt*,PetscScalar*,void*),void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscUseMethod(nep,"NEPNLEIGSGetSingularitiesFunction_C",(NEP,PetscErrorCode(**)(NEP,PetscInt*,PetscScalar*,void*),void**),(nep,fun,ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSSetRestart_NLEIGS(NEP nep,PetscReal keep)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (keep==(PetscReal)PETSC_DEFAULT || keep==(PetscReal)PETSC_DECIDE) ctx->keep = 0.5;
  else {
    PetscCheck(keep>=0.1 && keep<=0.9,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument must be in the range [0.1,0.9]");
    ctx->keep = keep;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSSetRestart - Sets the restart parameter for the NLEIGS
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective

   Input Parameters:
+  nep  - the nonlinear eigensolver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -nep_nleigs_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: NEPNLEIGSGetRestart()
@*/
PetscErrorCode NEPNLEIGSSetRestart(NEP nep,PetscReal keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,keep,2);
  PetscTryMethod(nep,"NEPNLEIGSSetRestart_C",(NEP,PetscReal),(nep,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSGetRestart_NLEIGS(NEP nep,PetscReal *keep)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  *keep = ctx->keep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSGetRestart - Gets the restart parameter used in the NLEIGS method.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: NEPNLEIGSSetRestart()
@*/
PetscErrorCode NEPNLEIGSGetRestart(NEP nep,PetscReal *keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(keep,2);
  PetscUseMethod(nep,"NEPNLEIGSGetRestart_C",(NEP,PetscReal*),(nep,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSSetLocking_NLEIGS(NEP nep,PetscBool lock)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSSetLocking - Choose between locking and non-locking variants of
   the NLEIGS method.

   Logically Collective

   Input Parameters:
+  nep  - the nonlinear eigensolver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -nep_nleigs_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: NEPNLEIGSGetLocking()
@*/
PetscErrorCode NEPNLEIGSSetLocking(NEP nep,PetscBool lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,lock,2);
  PetscTryMethod(nep,"NEPNLEIGSSetLocking_C",(NEP,PetscBool),(nep,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSGetLocking_NLEIGS(NEP nep,PetscBool *lock)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSGetLocking - Gets the locking flag used in the NLEIGS method.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: NEPNLEIGSSetLocking()
@*/
PetscErrorCode NEPNLEIGSGetLocking(NEP nep,PetscBool *lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(lock,2);
  PetscUseMethod(nep,"NEPNLEIGSGetLocking_C",(NEP,PetscBool*),(nep,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSSetInterpolation_NLEIGS(NEP nep,PetscReal tol,PetscInt degree)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (tol == (PetscReal)PETSC_DETERMINE) {
    ctx->ddtol = PETSC_DETERMINE;
    nep->state = NEP_STATE_INITIAL;
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    ctx->ddtol = tol;
  }
  if (degree == PETSC_DETERMINE) {
    ctx->ddmaxit = 0;
    if (nep->state) PetscCall(NEPReset(nep));
    nep->state = NEP_STATE_INITIAL;
  } else if (degree != PETSC_CURRENT) {
    PetscCheck(degree>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of degree. Must be > 0");
    if (ctx->ddmaxit != degree) {
      ctx->ddmaxit = degree;
      if (nep->state) PetscCall(NEPReset(nep));
      nep->state = NEP_STATE_INITIAL;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSSetInterpolation - Sets the tolerance and maximum degree
   when building the interpolation via divided differences.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  tol    - tolerance to stop computing divided differences
-  degree - maximum degree of interpolation

   Options Database Key:
+  -nep_nleigs_interpolation_tol <tol> - Sets the tolerance to stop computing divided differences
-  -nep_nleigs_interpolation_degree <degree> - Sets the maximum degree of interpolation

   Notes:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: NEPNLEIGSGetInterpolation()
@*/
PetscErrorCode NEPNLEIGSSetInterpolation(NEP nep,PetscReal tol,PetscInt degree)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,tol,2);
  PetscValidLogicalCollectiveInt(nep,degree,3);
  PetscTryMethod(nep,"NEPNLEIGSSetInterpolation_C",(NEP,PetscReal,PetscInt),(nep,tol,degree));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSGetInterpolation_NLEIGS(NEP nep,PetscReal *tol,PetscInt *degree)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (tol)    *tol    = ctx->ddtol;
  if (degree) *degree = ctx->ddmaxit;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSGetInterpolation - Gets the tolerance and maximum degree
   when building the interpolation via divided differences.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  tol    - tolerance to stop computing divided differences
-  degree - maximum degree of interpolation

   Level: advanced

.seealso: NEPNLEIGSSetInterpolation()
@*/
PetscErrorCode NEPNLEIGSGetInterpolation(NEP nep,PetscReal *tol,PetscInt *degree)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscTryMethod(nep,"NEPNLEIGSGetInterpolation_C",(NEP,PetscReal*,PetscInt*),(nep,tol,degree));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSSetRKShifts_NLEIGS(NEP nep,PetscInt ns,PetscScalar *shifts)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCheck(ns>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"Number of shifts must be non-negative");
  if (ctx->nshifts) PetscCall(PetscFree(ctx->shifts));
  for (i=0;i<ctx->nshiftsw;i++) PetscCall(KSPDestroy(&ctx->ksp[i]));
  PetscCall(PetscFree(ctx->ksp));
  ctx->ksp = NULL;
  if (ns) {
    PetscCall(PetscMalloc1(ns,&ctx->shifts));
    for (i=0;i<ns;i++) ctx->shifts[i] = shifts[i];
  }
  ctx->nshifts = ns;
  nep->state   = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSSetRKShifts - Sets a list of shifts to be used in the Rational
   Krylov method.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  ns     - number of shifts
-  shifts - array of scalar values specifying the shifts

   Options Database Key:
.  -nep_nleigs_rk_shifts - Sets the list of shifts

   Notes:
   If only one shift is provided, the built subspace built is equivalent to
   shift-and-invert Krylov-Schur (provided that the absolute convergence
   criterion is used).

   In the case of real scalars, complex shifts are not allowed. In the
   command line, a comma-separated list of complex values can be provided with
   the format [+/-][realnumber][+/-]realnumberi with no spaces, e.g.
   -nep_nleigs_rk_shifts 1.0+2.0i,1.5+2.0i,1.0+1.5i

   Use ns=0 to remove previously set shifts.

   Level: advanced

.seealso: NEPNLEIGSGetRKShifts()
@*/
PetscErrorCode NEPNLEIGSSetRKShifts(NEP nep,PetscInt ns,PetscScalar shifts[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,ns,2);
  if (ns) PetscAssertPointer(shifts,3);
  PetscTryMethod(nep,"NEPNLEIGSSetRKShifts_C",(NEP,PetscInt,PetscScalar*),(nep,ns,shifts));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSGetRKShifts_NLEIGS(NEP nep,PetscInt *ns,PetscScalar **shifts)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       i;

  PetscFunctionBegin;
  *ns = ctx->nshifts;
  if (ctx->nshifts) {
    PetscCall(PetscMalloc1(ctx->nshifts,shifts));
    for (i=0;i<ctx->nshifts;i++) (*shifts)[i] = ctx->shifts[i];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPNLEIGSGetRKShifts - Gets the list of shifts used in the Rational
   Krylov method.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  ns     - number of shifts
-  shifts - array of shifts

   Note:
   The user is responsible for deallocating the returned array.

   Level: advanced

.seealso: NEPNLEIGSSetRKShifts()
@*/
PetscErrorCode NEPNLEIGSGetRKShifts(NEP nep,PetscInt *ns,PetscScalar *shifts[]) PeNS
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(ns,2);
  PetscAssertPointer(shifts,3);
  PetscTryMethod(nep,"NEPNLEIGSGetRKShifts_C",(NEP,PetscInt*,PetscScalar**),(nep,ns,shifts));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSGetKSPs_NLEIGS(NEP nep,PetscInt *nsolve,KSP **ksp)
{
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt       i;
  PC             pc;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    PetscCall(NEPNLEIGSSetShifts(nep,&ctx->nshiftsw));
    PetscCall(PetscMalloc1(ctx->nshiftsw,&ctx->ksp));
    for (i=0;i<ctx->nshiftsw;i++) {
      PetscCall(KSPCreate(PetscObjectComm((PetscObject)nep),&ctx->ksp[i]));
      PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->ksp[i],(PetscObject)nep,1));
      PetscCall(KSPSetOptionsPrefix(ctx->ksp[i],((PetscObject)nep)->prefix));
      PetscCall(KSPAppendOptionsPrefix(ctx->ksp[i],"nep_nleigs_"));
      PetscCall(PetscObjectSetOptions((PetscObject)ctx->ksp[i],((PetscObject)nep)->options));
      PetscCall(KSPSetErrorIfNotConverged(ctx->ksp[i],PETSC_TRUE));
      PetscCall(KSPSetTolerances(ctx->ksp[i],1e-3*SlepcDefaultTol(nep->tol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
      PetscCall(KSPGetPC(ctx->ksp[i],&pc));
      if ((nep->fui==NEP_USER_INTERFACE_SPLIT && nep->P) || (nep->fui==NEP_USER_INTERFACE_CALLBACK && nep->function_pre!=nep->function)) {
        PetscCall(KSPSetType(ctx->ksp[i],KSPBCGS));
        PetscCall(PCSetType(pc,PCBJACOBI));
      } else {
        PetscCall(KSPSetType(ctx->ksp[i],KSPPREONLY));
        PetscCall(PCSetType(pc,PCLU));
      }
    }
  }
  if (nsolve) *nsolve = ctx->nshiftsw;
  if (ksp)    *ksp    = ctx->ksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPNLEIGSGetKSPs - Retrieve the array of linear solver objects associated with
   the nonlinear eigenvalue solver.

   Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameters:
+  nsolve - number of returned KSP objects
-  ksp - array of linear solver object

   Notes:
   The number of KSP objects is equal to the number of shifts provided by the user,
   or 1 if the user did not provide shifts.

   Level: advanced

.seealso: NEPNLEIGSSetRKShifts()
@*/
PetscErrorCode NEPNLEIGSGetKSPs(NEP nep,PetscInt *nsolve,KSP **ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscUseMethod(nep,"NEPNLEIGSGetKSPs_C",(NEP,PetscInt*,KSP**),(nep,nsolve,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSSetFullBasis_NLEIGS(NEP nep,PetscBool fullbasis)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (fullbasis!=ctx->fullbasis) {
    ctx->fullbasis = fullbasis;
    nep->state     = NEP_STATE_INITIAL;
    nep->useds     = PetscNot(fullbasis);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSSetFullBasis - Choose between TOAR-basis (default) and full-basis
   variants of the NLEIGS method.

   Logically Collective

   Input Parameters:
+  nep  - the nonlinear eigensolver context
-  fullbasis - true if the full-basis variant must be selected

   Options Database Key:
.  -nep_nleigs_full_basis - Sets the full-basis flag

   Notes:
   The default is to use a compact representation of the Krylov basis, that is,
   V = (I otimes U) S, with a tensor BV. This behaviour can be changed so that
   the full basis V is explicitly stored and operated with. This variant is more
   expensive in terms of memory and computation, but is necessary in some cases,
   particularly for two-sided computations, see NEPSetTwoSided().

   In the full-basis variant, the NLEIGS solver uses an EPS object to explicitly
   solve the linearized eigenproblem, see NEPNLEIGSGetEPS().

   Level: advanced

.seealso: NEPNLEIGSGetFullBasis(), NEPNLEIGSGetEPS(), NEPSetTwoSided(), BVCreateTensor()
@*/
PetscErrorCode NEPNLEIGSSetFullBasis(NEP nep,PetscBool fullbasis)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,fullbasis,2);
  PetscTryMethod(nep,"NEPNLEIGSSetFullBasis_C",(NEP,PetscBool),(nep,fullbasis));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPNLEIGSGetFullBasis_NLEIGS(NEP nep,PetscBool *fullbasis)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  *fullbasis = ctx->fullbasis;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPNLEIGSGetFullBasis - Gets the flag that indicates if NLEIGS is using the
   full-basis variant.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  fullbasis - the flag

   Level: advanced

.seealso: NEPNLEIGSSetFullBasis()
@*/
PetscErrorCode NEPNLEIGSGetFullBasis(NEP nep,PetscBool *fullbasis)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(fullbasis,2);
  PetscUseMethod(nep,"NEPNLEIGSGetFullBasis_C",(NEP,PetscBool*),(nep,fullbasis));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#define SHIFTMAX 30

static PetscErrorCode NEPSetFromOptions_NLEIGS(NEP nep,PetscOptionItems PetscOptionsObject)
{
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt       i=0,k;
  PetscBool      flg1,flg2,b;
  PetscReal      r;
  PetscScalar    array[SHIFTMAX];

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"NEP NLEIGS Options");

    PetscCall(PetscOptionsReal("-nep_nleigs_restart","Proportion of vectors kept after restart","NEPNLEIGSSetRestart",0.5,&r,&flg1));
    if (flg1) PetscCall(NEPNLEIGSSetRestart(nep,r));

    PetscCall(PetscOptionsBool("-nep_nleigs_locking","Choose between locking and non-locking variants","NEPNLEIGSSetLocking",PETSC_FALSE,&b,&flg1));
    if (flg1) PetscCall(NEPNLEIGSSetLocking(nep,b));

    PetscCall(PetscOptionsBool("-nep_nleigs_full_basis","Choose between TOAR and full-basis variants","NEPNLEIGSSetFullBasis",PETSC_FALSE,&b,&flg1));
    if (flg1) PetscCall(NEPNLEIGSSetFullBasis(nep,b));

    PetscCall(NEPNLEIGSGetInterpolation(nep,&r,&i));
    if (!i) i = PETSC_DETERMINE;
    PetscCall(PetscOptionsInt("-nep_nleigs_interpolation_degree","Maximum number of terms for interpolation via divided differences","NEPNLEIGSSetInterpolation",i,&i,&flg1));
    PetscCall(PetscOptionsReal("-nep_nleigs_interpolation_tol","Tolerance for interpolation via divided differences","NEPNLEIGSSetInterpolation",r,&r,&flg2));
    if (flg1 || flg2) PetscCall(NEPNLEIGSSetInterpolation(nep,r,i));

    k = SHIFTMAX;
    for (i=0;i<k;i++) array[i] = 0;
    PetscCall(PetscOptionsScalarArray("-nep_nleigs_rk_shifts","Shifts for Rational Krylov","NEPNLEIGSSetRKShifts",array,&k,&flg1));
    if (flg1) PetscCall(NEPNLEIGSSetRKShifts(nep,k,array));

  PetscOptionsHeadEnd();

  if (!ctx->ksp) PetscCall(NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp));
  for (i=0;i<ctx->nshiftsw;i++) PetscCall(KSPSetFromOptions(ctx->ksp[i]));

  if (ctx->fullbasis) {
    if (!ctx->eps) PetscCall(NEPNLEIGSGetEPS(nep,&ctx->eps));
    PetscCall(EPSSetFromOptions(ctx->eps));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPView_NLEIGS(NEP nep,PetscViewer viewer)
{
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscBool      isascii;
  PetscInt       i;
  char           str[50];

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*ctx->keep)));
    if (ctx->fullbasis) PetscCall(PetscViewerASCIIPrintf(viewer,"  using the full-basis variant\n"));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-"));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  divided difference terms: used=%" PetscInt_FMT ", max=%" PetscInt_FMT "\n",ctx->nmat,ctx->ddmaxit));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  tolerance for divided difference convergence: %g\n",(double)ctx->ddtol));
    if (ctx->nshifts) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"  RK shifts: "));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
      for (i=0;i<ctx->nshifts;i++) {
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),ctx->shifts[i],PETSC_FALSE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"%s%s",str,(i<ctx->nshifts-1)?",":""));
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    }
    if (!ctx->ksp) PetscCall(NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(KSPView(ctx->ksp[0],viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
    if (ctx->fullbasis) {
      if (!ctx->eps) PetscCall(NEPNLEIGSGetEPS(nep,&ctx->eps));
      PetscCall(PetscViewerASCIIPushTab(viewer));
      PetscCall(EPSView(ctx->eps,viewer));
      PetscCall(PetscViewerASCIIPopTab(viewer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPReset_NLEIGS(NEP nep)
{
  PetscInt       k;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) PetscCall(PetscFree(ctx->coeffD));
  else {
    for (k=0;k<ctx->nmat;k++) PetscCall(MatDestroy(&ctx->D[k]));
  }
  PetscCall(PetscFree4(ctx->s,ctx->xi,ctx->beta,ctx->D));
  for (k=0;k<ctx->nshiftsw;k++) PetscCall(KSPReset(ctx->ksp[k]));
  PetscCall(VecDestroy(&ctx->vrn));
  if (ctx->fullbasis) {
    PetscCall(MatDestroy(&ctx->A));
    PetscCall(EPSReset(ctx->eps));
    for (k=0;k<4;k++) PetscCall(VecDestroy(&ctx->w[k]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDestroy_NLEIGS(NEP nep)
{
  PetscInt       k;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&ctx->V));
  for (k=0;k<ctx->nshiftsw;k++) PetscCall(KSPDestroy(&ctx->ksp[k]));
  PetscCall(PetscFree(ctx->ksp));
  if (ctx->nshifts) PetscCall(PetscFree(ctx->shifts));
  if (ctx->fullbasis) PetscCall(EPSDestroy(&ctx->eps));
  PetscCall(PetscFree(nep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetSingularitiesFunction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetSingularitiesFunction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetInterpolation_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetInterpolation_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetRKShifts_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetRKShifts_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetKSPs_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetFullBasis_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetFullBasis_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetEPS_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_NLEIGS(NEP nep)
{
  NEP_NLEIGS     *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  nep->data  = (void*)ctx;
  ctx->lock  = PETSC_TRUE;
  ctx->ddtol = PETSC_DETERMINE;

  nep->useds = PETSC_TRUE;

  nep->ops->setup          = NEPSetUp_NLEIGS;
  nep->ops->setfromoptions = NEPSetFromOptions_NLEIGS;
  nep->ops->view           = NEPView_NLEIGS;
  nep->ops->destroy        = NEPDestroy_NLEIGS;
  nep->ops->reset          = NEPReset_NLEIGS;

  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetSingularitiesFunction_C",NEPNLEIGSSetSingularitiesFunction_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetSingularitiesFunction_C",NEPNLEIGSGetSingularitiesFunction_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetRestart_C",NEPNLEIGSSetRestart_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetRestart_C",NEPNLEIGSGetRestart_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetLocking_C",NEPNLEIGSSetLocking_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetLocking_C",NEPNLEIGSGetLocking_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetInterpolation_C",NEPNLEIGSSetInterpolation_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetInterpolation_C",NEPNLEIGSGetInterpolation_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetRKShifts_C",NEPNLEIGSSetRKShifts_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetRKShifts_C",NEPNLEIGSGetRKShifts_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetKSPs_C",NEPNLEIGSGetKSPs_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetFullBasis_C",NEPNLEIGSSetFullBasis_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetFullBasis_C",NEPNLEIGSGetFullBasis_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetEPS_C",NEPNLEIGSSetEPS_NLEIGS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetEPS_C",NEPNLEIGSGetEPS_NLEIGS));
  PetscFunctionReturn(PETSC_SUCCESS);
}
