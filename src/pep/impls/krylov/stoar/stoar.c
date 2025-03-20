/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc polynomial eigensolver: "stoar"

   Method: S-TOAR

   Algorithm:

       Symmetric Two-Level Orthogonal Arnoldi.

   References:

       [1] C. Campos and J.E. Roman, "Restarted Q-Arnoldi-type methods
           exploiting symmetry in quadratic eigenvalue problems", BIT
           Numer. Math. 56(4):1213-1236, 2016.
*/

#include <slepc/private/pepimpl.h>         /*I "slepcpep.h" I*/
#include "../src/pep/impls/krylov/pepkrylov.h"
#include <slepcblaslapack.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-stoar,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"Restarted {Q-Arnoldi-type} methods exploiting symmetry in quadratic eigenvalue problems\",\n"
  "   journal = \"{BIT} Numer. Math.\",\n"
  "   volume = \"56\",\n"
  "   number = \"4\",\n"
  "   pages = \"1213--1236\",\n"
  "   year = \"2016,\"\n"
  "   doi = \"https://doi.org/10.1007/s10543-016-0601-5\"\n"
  "}\n";

typedef struct {
  PetscReal   scal[2];
  Mat         A[2];
  Vec         t;
} PEP_STOAR_MATSHELL;

static PetscErrorCode MatMult_STOAR(Mat A,Vec x,Vec y)
{
  PEP_STOAR_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(MatMult(ctx->A[0],x,y));
  PetscCall(VecScale(y,ctx->scal[0]));
  if (ctx->scal[1]) {
    PetscCall(MatMult(ctx->A[1],x,ctx->t));
    PetscCall(VecAXPY(y,ctx->scal[1],ctx->t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_STOAR(Mat A)
{
  PEP_STOAR_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(VecDestroy(&ctx->t));
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPSTOARSetUpInnerMatrix(PEP pep,Mat *B)
{
  Mat                pB[4],Bs[3],D[3];
  PetscInt           i,j,n,m;
  PEP_STOAR_MATSHELL *ctxMat[3];
  PEP_STOAR          *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  for (i=0;i<3;i++) {
    PetscCall(STGetMatrixTransformed(pep->st,i,&D[i])); /* D[2] = M */
  }
  PetscCall(MatGetLocalSize(D[2],&m,&n));

  for (j=0;j<3;j++) {
    PetscCall(PetscNew(ctxMat+j));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)pep),m,n,PETSC_DETERMINE,PETSC_DETERMINE,ctxMat[j],&Bs[j]));
    PetscCall(MatShellSetOperation(Bs[j],MATOP_MULT,(void(*)(void))MatMult_STOAR));
    PetscCall(MatShellSetOperation(Bs[j],MATOP_DESTROY,(void(*)(void))MatDestroy_STOAR));
  }
  for (i=0;i<4;i++) pB[i] = NULL;
  if (ctx->alpha) {
    ctxMat[0]->A[0] = D[0]; ctxMat[0]->scal[0] = ctx->alpha; ctxMat[0]->scal[1] = 0.0;
    ctxMat[2]->A[0] = D[2]; ctxMat[2]->scal[0] = -ctx->alpha*pep->sfactor*pep->sfactor; ctxMat[2]->scal[1] = 0.0;
    pB[0] = Bs[0]; pB[3] = Bs[2];
  }
  if (ctx->beta) {
    i = (ctx->alpha)?1:0;
    ctxMat[0]->scal[1] = 0.0;
    ctxMat[0]->A[i] = D[1]; ctxMat[0]->scal[i] = -ctx->beta*pep->sfactor;
    ctxMat[1]->A[0] = D[2]; ctxMat[1]->scal[0] = -ctx->beta*pep->sfactor*pep->sfactor; ctxMat[1]->scal[1] = 0.0;
    pB[0] = Bs[0]; pB[1] = pB[2] = Bs[1];
  }
  PetscCall(BVCreateVec(pep->V,&ctxMat[0]->t));
  PetscCall(MatCreateNest(PetscObjectComm((PetscObject)pep),2,NULL,2,NULL,pB,B));
  for (j=0;j<3;j++) PetscCall(MatDestroy(&Bs[j]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetUp_STOAR(PEP pep)
{
  PetscBool         sinv,flg;
  PEP_STOAR         *ctx = (PEP_STOAR*)pep->data;
  PetscInt          ld,i;
  PetscReal         eta;
  BVOrthogType      otype;
  BVOrthogBlockType obtype;

  PetscFunctionBegin;
  PEPCheckHermitian(pep);
  PEPCheckQuadratic(pep);
  PEPCheckShiftSinvert(pep);
  /* spectrum slicing requires special treatment of default values */
  if (pep->which==PEP_ALL) {
    pep->ops->solve = PEPSolve_STOAR_QSlice;
    pep->ops->extractvectors = NULL;
    pep->ops->setdefaultst   = NULL;
    PetscCall(PEPSetUp_STOAR_QSlice(pep));
  } else {
    PetscCall(PEPSetDimensions_Default(pep,pep->nev,&pep->ncv,&pep->mpd));
    PetscCheck(ctx->lock || pep->mpd>=pep->ncv,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Should not use mpd parameter in non-locking variant");
    if (pep->max_it==PETSC_DETERMINE) pep->max_it = PetscMax(100,2*(pep->nmat-1)*pep->n/pep->ncv);
    pep->ops->solve = PEPSolve_STOAR;
    ld   = pep->ncv+2;
    PetscCall(DSSetType(pep->ds,DSGHIEP));
    PetscCall(DSSetCompact(pep->ds,PETSC_TRUE));
    PetscCall(DSSetExtraRow(pep->ds,PETSC_TRUE));
    PetscCall(DSAllocate(pep->ds,ld));
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
  }
  if (!pep->which) PetscCall(PEPSetWhichEigenpairs_Default(pep));
  PEPCheckUnsupported(pep,PEP_FEATURE_NONMONOMIAL | PEP_FEATURE_REGION);

  PetscCall(PEPAllocateSolution(pep,2));
  PetscCall(PEPSetWorkVecs(pep,4));
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVCreateTensor(pep->V,pep->nmat-1,&ctx->V));
  PetscCall(BVGetOrthogonalization(pep->V,&otype,NULL,&eta,&obtype));
  PetscCall(BVSetOrthogonalization(ctx->V,otype,BV_ORTHOG_REFINE_ALWAYS,eta,obtype));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute a run of Lanczos iterations. dim(work)=(ctx->ld)*4
*/
static PetscErrorCode PEPSTOARrun(PEP pep,PetscReal *a,PetscReal *b,PetscReal *omega,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscBool *symmlost,Vec *t_)
{
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscInt       i,j,m=*M,l,lock;
  PetscInt       lds,d,ld,offq,nqt,ldds;
  Vec            v=t_[0],t=t_[1],q=t_[2];
  PetscReal      norm,sym=0.0,fro=0.0,*f;
  PetscScalar    *y,*S,*x,sigma;
  PetscBLASInt   j_,one=1;
  PetscBool      lindep,flg,sinvert=PETSC_FALSE;
  Mat            MS;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(*M,&y));
  PetscCall(BVGetSizes(pep->V,NULL,NULL,&ld));
  PetscCall(BVTensorGetDegree(ctx->V,&d));
  PetscCall(BVGetActiveColumns(pep->V,&lock,&nqt));
  lds = d*ld;
  offq = ld;
  PetscCall(DSGetLeadingDimension(pep->ds,&ldds));
  *breakdown = PETSC_FALSE; /* ----- */
  PetscCall(DSGetDimensions(pep->ds,NULL,&l,NULL,NULL));
  PetscCall(BVSetActiveColumns(ctx->V,0,m));
  PetscCall(BVSetActiveColumns(pep->V,0,nqt));
  PetscCall(STGetTransform(pep->st,&flg));
  if (!flg) {
    /* spectral transformation handled by the solver */
    PetscCall(PetscObjectTypeCompareAny((PetscObject)pep->st,&flg,STSINVERT,STSHIFT,""));
    PetscCheck(flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"ST type not supported for TOAR without transforming matrices");
    PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinvert));
    PetscCall(STGetShift(pep->st,&sigma));
  }
  for (j=k;j<m;j++) {
    /* apply operator */
    PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
    PetscCall(MatDenseGetArray(MS,&S));
    PetscCall(BVGetColumn(pep->V,nqt,&t));
    PetscCall(BVMultVec(pep->V,1.0,0.0,v,S+j*lds));
    if (!sinvert) {
      PetscCall(STMatMult(pep->st,0,v,q));
      PetscCall(BVMultVec(pep->V,1.0,0.0,v,S+offq+j*lds));
      PetscCall(STMatMult(pep->st,1,v,t));
      PetscCall(VecAXPY(q,pep->sfactor,t));
      if (ctx->beta && ctx->alpha) {
        PetscCall(STMatMult(pep->st,2,v,t));
        PetscCall(VecAXPY(q,-pep->sfactor*pep->sfactor*ctx->beta/ctx->alpha,t));
      }
      PetscCall(STMatSolve(pep->st,q,t));
      PetscCall(VecScale(t,-1.0/(pep->sfactor*pep->sfactor)));
    } else {
      PetscCall(STMatMult(pep->st,1,v,q));
      PetscCall(STMatMult(pep->st,2,v,t));
      PetscCall(VecAXPY(q,sigma*pep->sfactor,t));
      PetscCall(VecScale(q,pep->sfactor));
      PetscCall(BVMultVec(pep->V,1.0,0.0,v,S+offq+j*lds));
      PetscCall(STMatMult(pep->st,2,v,t));
      PetscCall(VecAXPY(q,pep->sfactor*pep->sfactor,t));
      PetscCall(STMatSolve(pep->st,q,t));
      PetscCall(VecScale(t,-1.0));
    }
    PetscCall(BVRestoreColumn(pep->V,nqt,&t));

    /* orthogonalize */
    if (!sinvert) x = S+offq+(j+1)*lds;
    else x = S+(j+1)*lds;
    PetscCall(BVOrthogonalizeColumn(pep->V,nqt,x,&norm,&lindep));

    if (!lindep) {
      if (!sinvert) *(S+offq+(j+1)*lds+nqt) = norm;
      else *(S+(j+1)*lds+nqt) = norm;
      PetscCall(BVScaleColumn(pep->V,nqt,1.0/norm));
      nqt++;
    }
    if (!sinvert) {
      for (i=0;i<=nqt-1;i++) *(S+(j+1)*lds+i) = *(S+offq+j*lds+i);
      if (ctx->beta && ctx->alpha) {
        for (i=0;i<=nqt-1;i++) *(S+(j+1)*lds+offq+i) -= *(S+(j+1)*lds+i)*ctx->beta/ctx->alpha;
      }
    } else for (i=0;i<nqt;i++) *(S+(j+1)*lds+offq+i) = *(S+j*lds+i)+sigma*(*(S+(j+1)*lds+i));
    PetscCall(BVSetActiveColumns(pep->V,0,nqt));
    PetscCall(MatDenseRestoreArray(MS,&S));
    PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));

    /* level-2 orthogonalization */
    PetscCall(BVOrthogonalizeColumn(ctx->V,j+1,y,&norm,&lindep));
    a[j] = PetscRealPart(y[j]);
    omega[j+1] = (norm > 0)?1.0:-1.0;
    PetscCall(BVScaleColumn(ctx->V,j+1,1.0/norm));
    b[j] = PetscAbsReal(norm);

    /* check symmetry */
    PetscCall(DSGetArrayReal(pep->ds,DS_MAT_T,&f));
    if (j==k) {
      for (i=l;i<j-1;i++) y[i] = PetscAbsScalar(y[i])-PetscAbsReal(f[2*ldds+i]);
      for (i=0;i<l;i++) y[i] = 0.0;
    }
    PetscCall(DSRestoreArrayReal(pep->ds,DS_MAT_T,&f));
    if (j>0) y[j-1] = PetscAbsScalar(y[j-1])-PetscAbsReal(b[j-1]);
    PetscCall(PetscBLASIntCast(j,&j_));
    sym = SlepcAbs(BLASnrm2_(&j_,y,&one),sym);
    fro = SlepcAbs(fro,SlepcAbs(a[j],b[j]));
    if (j>0) fro = SlepcAbs(fro,b[j-1]);
    if (sym/fro>PetscMax(PETSC_SQRT_MACHINE_EPSILON,10*pep->tol)) {
      *symmlost = PETSC_TRUE;
      *M=j;
      break;
    }
  }
  PetscCall(BVSetActiveColumns(pep->V,lock,nqt));
  PetscCall(BVSetActiveColumns(ctx->V,0,*M));
  PetscCall(PetscFree(y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if 0
static PetscErrorCode PEPSTOARpreKConvergence(PEP pep,PetscInt nv,PetscReal *norm,Vec *w)
{
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscBLASInt   n_,one=1;
  PetscInt       lds=2*ctx->ld;
  PetscReal      t1,t2;
  PetscScalar    *S=ctx->S;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(nv+2,&n_));
  t1 = BLASnrm2_(&n_,S+nv*2*ctx->ld,&one);
  t2 = BLASnrm2_(&n_,S+(nv*2+1)*ctx->ld,&one);
  *norm = SlepcAbs(t1,t2);
  PetscCall(BVSetActiveColumns(pep->V,0,nv+2));
  PetscCall(BVMultVec(pep->V,1.0,0.0,w[1],S+nv*lds));
  PetscCall(STMatMult(pep->st,0,w[1],w[2]));
  PetscCall(VecNorm(w[2],NORM_2,&t1));
  PetscCall(BVMultVec(pep->V,1.0,0.0,w[1],S+ctx->ld+nv*lds));
  PetscCall(STMatMult(pep->st,2,w[1],w[2]));
  PetscCall(VecNorm(w[2],NORM_2,&t2));
  t2 *= pep->sfactor*pep->sfactor;
  *norm = PetscMax(*norm,SlepcAbs(t1,t2));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

PetscErrorCode PEPSolve_STOAR(PEP pep)
{
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscInt       j,k,l,nv=0,ld,ldds,t,nq=0;
  PetscInt       nconv=0,deg=pep->nmat-1;
  PetscScalar    sigma;
  PetscReal      beta,norm=1.0,*omega,*a,*b;
  PetscBool      breakdown,symmlost=PETSC_FALSE,sinv=PETSC_FALSE,falselock=PETSC_TRUE,flg;
  Mat            MQ,A,D;
  Vec            vomega;

  PetscFunctionBegin;
  PetscCall(PetscCitationsRegister(citation,&cited));
  PetscCall(PEPSTOARSetUpInnerMatrix(pep,&A));
  PetscCall(BVSetMatrix(ctx->V,A,PETSC_TRUE));
  PetscCall(MatDestroy(&A));
  if (ctx->lock) {
    /* undocumented option to use a cheaper locking instead of the true locking */
    PetscCall(PetscOptionsGetBool(NULL,NULL,"-pep_stoar_falselocking",&falselock,NULL));
  }
  PetscCall(BVGetSizes(pep->V,NULL,NULL,&ld));
  PetscCall(STGetShift(pep->st,&sigma));
  PetscCall(STGetTransform(pep->st,&flg));
  if (pep->sfactor!=1.0) {
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

  /* Get the starting Arnoldi vector */
  PetscCall(BVTensorBuildFirstColumn(ctx->V,pep->nini));
  PetscCall(DSSetDimensions(pep->ds,1,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(BVSetActiveColumns(ctx->V,0,1));
  PetscCall(DSGetMatAndColumn(pep->ds,DS_MAT_D,0,&D,&vomega));
  PetscCall(BVGetSignature(ctx->V,vomega));
  PetscCall(DSRestoreMatAndColumn(pep->ds,DS_MAT_D,0,&D,&vomega));

  /* Restart loop */
  l = 0;
  PetscCall(DSGetLeadingDimension(pep->ds,&ldds));
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;
    PetscCall(DSGetArrayReal(pep->ds,DS_MAT_T,&a));
    b = a+ldds;
    PetscCall(DSGetArrayReal(pep->ds,DS_MAT_D,&omega));

    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(pep->nconv+pep->mpd,pep->ncv);
    PetscCall(DSSetDimensions(pep->ds,nv,pep->nconv,pep->nconv+l));
    PetscCall(PEPSTOARrun(pep,a,b,omega,pep->nconv+l,&nv,&breakdown,&symmlost,pep->work));
    beta = b[nv-1];
    if (symmlost && nv==pep->nconv+l) {
      pep->reason = PEP_DIVERGED_SYMMETRY_LOST;
      pep->nconv = nconv;
      if (falselock || !ctx->lock) {
       PetscCall(BVSetActiveColumns(ctx->V,0,pep->nconv));
       PetscCall(BVTensorCompress(ctx->V,0));
      }
      break;
    }
    PetscCall(DSRestoreArrayReal(pep->ds,DS_MAT_T,&a));
    PetscCall(DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega));
    PetscCall(DSSetDimensions(pep->ds,nv,pep->nconv,pep->nconv+l));
    if (l==0) PetscCall(DSSetState(pep->ds,DS_STATE_INTERMEDIATE));
    else PetscCall(DSSetState(pep->ds,DS_STATE_RAW));

    /* Solve projected problem */
    PetscCall(DSSolve(pep->ds,pep->eigr,pep->eigi));
    PetscCall(DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(pep->ds));
    PetscCall(DSSynchronize(pep->ds,pep->eigr,pep->eigi));

    /* Check convergence */
    /* PetscCall(PEPSTOARpreKConvergence(pep,nv,&norm,pep->work));*/
    norm = 1.0;
    PetscCall(DSGetDimensions(pep->ds,NULL,NULL,NULL,&t));
    PetscCall(PEPKrylovConvergence(pep,PETSC_FALSE,pep->nconv,t-pep->nconv,PetscAbsReal(beta)*norm,&k));
    PetscCall((*pep->stopping)(pep,pep->its,pep->max_it,k,pep->nev,&pep->reason,pep->stoppingctx));

    /* Update l */
    if (pep->reason != PEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)/2));
      l = PetscMin(l,t);
      PetscCall(DSGetTruncateSize(pep->ds,k,t,&l));
      if (!breakdown) {
        /* Prepare the Rayleigh quotient for restart */
        PetscCall(DSTruncate(pep->ds,k+l,PETSC_FALSE));
      }
    }
    nconv = k;
    if (!ctx->lock && pep->reason == PEP_CONVERGED_ITERATING && !breakdown) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(pep,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    /* Update S */
    PetscCall(DSGetMat(pep->ds,DS_MAT_Q,&MQ));
    PetscCall(BVMultInPlace(ctx->V,MQ,pep->nconv,k+l));
    PetscCall(DSRestoreMat(pep->ds,DS_MAT_Q,&MQ));

    /* Copy last column of S */
    PetscCall(BVCopyColumn(ctx->V,nv,k+l));
    PetscCall(BVSetActiveColumns(ctx->V,0,k+l));
    PetscCall(DSSetDimensions(pep->ds,k+l,PETSC_DETERMINE,PETSC_DETERMINE));
    PetscCall(DSGetMatAndColumn(pep->ds,DS_MAT_D,0,&D,&vomega));
    PetscCall(BVSetSignature(ctx->V,vomega));
    PetscCall(DSRestoreMatAndColumn(pep->ds,DS_MAT_D,0,&D,&vomega));

    if (breakdown && pep->reason == PEP_CONVERGED_ITERATING) {
      /* stop if breakdown */
      PetscCall(PetscInfo(pep,"Breakdown TOAR method (it=%" PetscInt_FMT " norm=%g)\n",pep->its,(double)beta));
      pep->reason = PEP_DIVERGED_BREAKDOWN;
    }
    if (pep->reason != PEP_CONVERGED_ITERATING) l--;
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
    PetscCall(BVSetActiveColumns(ctx->V,0,pep->nconv));
    PetscCall(BVGetActiveColumns(pep->V,NULL,&nq));
    PetscCall(BVSetActiveColumns(pep->V,0,nq));
    if (nq>pep->nconv) {
      PetscCall(BVTensorCompress(ctx->V,pep->nconv));
      PetscCall(BVSetActiveColumns(pep->V,0,pep->nconv));
    }
  }
  PetscCall(STGetTransform(pep->st,&flg));
  if (!flg) PetscTryTypeMethod(pep,backtransform);
  if (pep->sfactor!=1.0) {
    for (j=0;j<pep->nconv;j++) {
      pep->eigr[j] *= pep->sfactor;
      pep->eigi[j] *= pep->sfactor;
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

  PetscCall(DSTruncate(pep->ds,pep->nconv,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetFromOptions_STOAR(PEP pep,PetscOptionItems PetscOptionsObject)
{
  PetscBool      flg,lock,b,f1,f2,f3;
  PetscInt       i,j,k;
  PetscReal      array[2]={0,0};
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"PEP STOAR Options");

    PetscCall(PetscOptionsBool("-pep_stoar_locking","Choose between locking and non-locking variants","PEPSTOARSetLocking",PETSC_FALSE,&lock,&flg));
    if (flg) PetscCall(PEPSTOARSetLocking(pep,lock));

    b = ctx->detect;
    PetscCall(PetscOptionsBool("-pep_stoar_detect_zeros","Check zeros during factorizations at interval boundaries","PEPSTOARSetDetectZeros",ctx->detect,&b,&flg));
    if (flg) PetscCall(PEPSTOARSetDetectZeros(pep,b));

    i = 1;
    j = k = PETSC_DECIDE;
    PetscCall(PetscOptionsInt("-pep_stoar_nev","Number of eigenvalues to compute in each subsolve (only for spectrum slicing)","PEPSTOARSetDimensions",20,&i,&f1));
    PetscCall(PetscOptionsInt("-pep_stoar_ncv","Number of basis vectors in each subsolve (only for spectrum slicing)","PEPSTOARSetDimensions",40,&j,&f2));
    PetscCall(PetscOptionsInt("-pep_stoar_mpd","Maximum dimension of projected problem in each subsolve (only for spectrum slicing)","PEPSTOARSetDimensions",40,&k,&f3));
    if (f1 || f2 || f3) PetscCall(PEPSTOARSetDimensions(pep,i,j,k));

    k = 2;
    PetscCall(PetscOptionsRealArray("-pep_stoar_linearization","Parameters of the linearization","PEPSTOARSetLinearization",array,&k,&flg));
    if (flg) PetscCall(PEPSTOARSetLinearization(pep,array[0],array[1]));

    b = ctx->checket;
    PetscCall(PetscOptionsBool("-pep_stoar_check_eigenvalue_type","Check eigenvalue type during spectrum slicing","PEPSTOARSetCheckEigenvalueType",ctx->checket,&b,&flg));
    if (flg) PetscCall(PEPSTOARSetCheckEigenvalueType(pep,b));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARSetLocking_STOAR(PEP pep,PetscBool lock)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARSetLocking - Choose between locking and non-locking variants of
   the STOAR method.

   Logically Collective

   Input Parameters:
+  pep  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -pep_stoar_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: PEPSTOARGetLocking()
@*/
PetscErrorCode PEPSTOARSetLocking(PEP pep,PetscBool lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,lock,2);
  PetscTryMethod(pep,"PEPSTOARSetLocking_C",(PEP,PetscBool),(pep,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARGetLocking_STOAR(PEP pep,PetscBool *lock)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARGetLocking - Gets the locking flag used in the STOAR method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: PEPSTOARSetLocking()
@*/
PetscErrorCode PEPSTOARGetLocking(PEP pep,PetscBool *lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(lock,2);
  PetscUseMethod(pep,"PEPSTOARGetLocking_C",(PEP,PetscBool*),(pep,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARGetInertias_STOAR(PEP pep,PetscInt *n,PetscReal *shifts[],PetscInt *inertias[])
{
  PetscInt       i,numsh;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PEP_SR         sr = ctx->sr;

  PetscFunctionBegin;
  PetscCheck(pep->state,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"Must call PEPSetUp() first");
  PetscCheck(ctx->sr,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see PEPSetInterval()");
  switch (pep->state) {
  case PEP_STATE_INITIAL:
    break;
  case PEP_STATE_SETUP:
    if (n) *n = 2;
    if (shifts) {
      PetscCall(PetscMalloc1(2,shifts));
      (*shifts)[0] = pep->inta;
      (*shifts)[1] = pep->intb;
    }
    if (inertias) {
      PetscCall(PetscMalloc1(2,inertias));
      (*inertias)[0] = (sr->dir==1)?sr->inertia0:sr->inertia1;
      (*inertias)[1] = (sr->dir==1)?sr->inertia1:sr->inertia0;
    }
    break;
  case PEP_STATE_SOLVED:
  case PEP_STATE_EIGENVECTORS:
    numsh = ctx->nshifts;
    if (n) *n = numsh;
    if (shifts) {
      PetscCall(PetscMalloc1(numsh,shifts));
      for (i=0;i<numsh;i++) (*shifts)[i] = ctx->shifts[i];
    }
    if (inertias) {
      PetscCall(PetscMalloc1(numsh,inertias));
      for (i=0;i<numsh;i++) (*inertias)[i] = ctx->inertias[i];
    }
    break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPSTOARGetInertias - Gets the values of the shifts and their
   corresponding inertias in case of doing spectrum slicing for a
   computational interval.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameters:
+  n        - number of shifts, including the endpoints of the interval
.  shifts   - the values of the shifts used internally in the solver
-  inertias - the values of the inertia in each shift

   Notes:
   If called after PEPSolve(), all shifts used internally by the solver are
   returned (including both endpoints and any intermediate ones). If called
   before PEPSolve() and after PEPSetUp() then only the information of the
   endpoints of subintervals is available.

   This function is only available for spectrum slicing runs.

   The returned arrays should be freed by the user. Can pass NULL in any of
   the two arrays if not required.

   Fortran Notes:
   The calling sequence from Fortran is
.vb
   PEPSTOARGetInertias(pep,n,shifts,inertias,ierr)
   integer n
   double precision shifts(*)
   integer inertias(*)
.ve
   The arrays should be at least of length n. The value of n can be determined
   by an initial call
.vb
   PEPSTOARGetInertias(pep,n,PETSC_NULL_REAL,PETSC_NULL_INTEGER,ierr)
.ve

   Level: advanced

.seealso: PEPSetInterval()
@*/
PetscErrorCode PEPSTOARGetInertias(PEP pep,PetscInt *n,PetscReal *shifts[],PetscInt *inertias[]) PeNS
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(n,2);
  PetscUseMethod(pep,"PEPSTOARGetInertias_C",(PEP,PetscInt*,PetscReal**,PetscInt**),(pep,n,shifts,inertias));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARSetDetectZeros_STOAR(PEP pep,PetscBool detect)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ctx->detect = detect;
  pep->state  = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARSetDetectZeros - Sets a flag to enforce detection of
   zeros during the factorizations throughout the spectrum slicing computation.

   Logically Collective

   Input Parameters:
+  pep    - the eigenproblem solver context
-  detect - check for zeros

   Options Database Key:
.  -pep_stoar_detect_zeros - Check for zeros; this takes an optional
   bool value (0/1/no/yes/true/false)

   Notes:
   A zero in the factorization indicates that a shift coincides with an eigenvalue.

   This flag is turned off by default, and may be necessary in some cases.
   This feature currently requires an external package for factorizations
   with support for zero detection, e.g. MUMPS.

   Level: advanced

.seealso: PEPSetInterval()
@*/
PetscErrorCode PEPSTOARSetDetectZeros(PEP pep,PetscBool detect)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,detect,2);
  PetscTryMethod(pep,"PEPSTOARSetDetectZeros_C",(PEP,PetscBool),(pep,detect));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARGetDetectZeros_STOAR(PEP pep,PetscBool *detect)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  *detect = ctx->detect;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARGetDetectZeros - Gets the flag that enforces zero detection
   in spectrum slicing.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  detect - whether zeros detection is enforced during factorizations

   Level: advanced

.seealso: PEPSTOARSetDetectZeros()
@*/
PetscErrorCode PEPSTOARGetDetectZeros(PEP pep,PetscBool *detect)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(detect,2);
  PetscUseMethod(pep,"PEPSTOARGetDetectZeros_C",(PEP,PetscBool*),(pep,detect));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARSetLinearization_STOAR(PEP pep,PetscReal alpha,PetscReal beta)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  PetscCheck(beta!=0.0 || alpha!=0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"Parameters alpha and beta cannot be zero simultaneously");
  ctx->alpha = alpha;
  ctx->beta  = beta;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARSetLinearization - Set the coefficients that define
   the linearization of a quadratic eigenproblem.

   Logically Collective

   Input Parameters:
+  pep   - polynomial eigenvalue solver
.  alpha - first parameter of the linearization
-  beta  - second parameter of the linearization

   Options Database Key:
.  -pep_stoar_linearization <alpha,beta> - Sets the coefficients

   Notes:
   Cannot pass zero for both alpha and beta. The default values are
   alpha=1 and beta=0.

   Level: advanced

.seealso: PEPSTOARGetLinearization()
@*/
PetscErrorCode PEPSTOARSetLinearization(PEP pep,PetscReal alpha,PetscReal beta)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,alpha,2);
  PetscValidLogicalCollectiveReal(pep,beta,3);
  PetscTryMethod(pep,"PEPSTOARSetLinearization_C",(PEP,PetscReal,PetscReal),(pep,alpha,beta));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARGetLinearization_STOAR(PEP pep,PetscReal *alpha,PetscReal *beta)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  if (alpha) *alpha = ctx->alpha;
  if (beta)  *beta  = ctx->beta;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARGetLinearization - Returns the coefficients that define
   the linearization of a quadratic eigenproblem.

   Not Collective

   Input Parameter:
.  pep  - polynomial eigenvalue solver

   Output Parameters:
+  alpha - the first parameter of the linearization
-  beta  - the second parameter of the linearization

   Level: advanced

.seealso: PEPSTOARSetLinearization()
@*/
PetscErrorCode PEPSTOARGetLinearization(PEP pep,PetscReal *alpha,PetscReal *beta)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscUseMethod(pep,"PEPSTOARGetLinearization_C",(PEP,PetscReal*,PetscReal*),(pep,alpha,beta));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARSetDimensions_STOAR(PEP pep,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  if (nev != PETSC_CURRENT) {
    PetscCheck(nev>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
    ctx->nev = nev;
  }
  if (ncv == PETSC_DETERMINE) {
    ctx->ncv = PETSC_DETERMINE;
  } else if (ncv != PETSC_CURRENT) {
    PetscCheck(ncv>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    ctx->ncv = ncv;
  }
  if (mpd == PETSC_DETERMINE) {
    ctx->mpd = PETSC_DETERMINE;
  } else if (mpd != PETSC_CURRENT) {
    PetscCheck(mpd>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    ctx->mpd = mpd;
  }
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARSetDimensions - Sets the dimensions used for each subsolve
   step in case of doing spectrum slicing for a computational interval.
   The meaning of the parameters is the same as in PEPSetDimensions().

   Logically Collective

   Input Parameters:
+  pep - the eigenproblem solver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the subsolve
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Key:
+  -pep_stoar_nev <nev> - Sets the number of eigenvalues
.  -pep_stoar_ncv <ncv> - Sets the dimension of the subspace
-  -pep_stoar_mpd <mpd> - Sets the maximum projected dimension

   Level: advanced

.seealso: PEPSTOARGetDimensions(), PEPSetDimensions(), PEPSetInterval()
@*/
PetscErrorCode PEPSTOARSetDimensions(PEP pep,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,nev,2);
  PetscValidLogicalCollectiveInt(pep,ncv,3);
  PetscValidLogicalCollectiveInt(pep,mpd,4);
  PetscTryMethod(pep,"PEPSTOARSetDimensions_C",(PEP,PetscInt,PetscInt,PetscInt),(pep,nev,ncv,mpd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARGetDimensions_STOAR(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  if (nev) *nev = ctx->nev;
  if (ncv) *ncv = ctx->ncv;
  if (mpd) *mpd = ctx->mpd;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARGetDimensions - Gets the dimensions used for each subsolve
   step in case of doing spectrum slicing for a computational interval.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the subsolve
-  mpd - the maximum dimension allowed for the projected problem

   Level: advanced

.seealso: PEPSTOARSetDimensions()
@*/
PetscErrorCode PEPSTOARGetDimensions(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscUseMethod(pep,"PEPSTOARGetDimensions_C",(PEP,PetscInt*,PetscInt*,PetscInt*),(pep,nev,ncv,mpd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARSetCheckEigenvalueType_STOAR(PEP pep,PetscBool checket)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ctx->checket = checket;
  pep->state   = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARSetCheckEigenvalueType - Sets a flag to check that all the eigenvalues
   obtained throughout the spectrum slicing computation have the same definite type.

   Logically Collective

   Input Parameters:
+  pep     - the eigenproblem solver context
-  checket - check eigenvalue type

   Options Database Key:
.  -pep_stoar_check_eigenvalue_type - Check eigenvalue type; this takes an optional
   bool value (0/1/no/yes/true/false)

   Notes:
   This option is relevant only for spectrum slicing computations, but it is
   ignored if the problem type is PEP_HYPERBOLIC.

   This flag is turned on by default, to guarantee that the computed eigenvalues
   have the same type (otherwise the computed solution might be wrong). But since
   the check is computationally quite expensive, the check may be turned off if
   the user knows for sure that all eigenvalues in the requested interval have
   the same type.

   Level: advanced

.seealso: PEPSetProblemType(), PEPSetInterval()
@*/
PetscErrorCode PEPSTOARSetCheckEigenvalueType(PEP pep,PetscBool checket)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,checket,2);
  PetscTryMethod(pep,"PEPSTOARSetCheckEigenvalueType_C",(PEP,PetscBool),(pep,checket));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOARGetCheckEigenvalueType_STOAR(PEP pep,PetscBool *checket)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  *checket = ctx->checket;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSTOARGetCheckEigenvalueType - Gets the flag for the eigenvalue type
   check in spectrum slicing.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  checket - whether eigenvalue type must be checked during spectrum slcing

   Level: advanced

.seealso: PEPSTOARSetCheckEigenvalueType()
@*/
PetscErrorCode PEPSTOARGetCheckEigenvalueType(PEP pep,PetscBool *checket)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(checket,2);
  PetscUseMethod(pep,"PEPSTOARGetCheckEigenvalueType_C",(PEP,PetscBool*),(pep,checket));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPView_STOAR(PEP pep,PetscViewer viewer)
{
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-"));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  linearization parameters: alpha=%g beta=%g\n",(double)ctx->alpha,(double)ctx->beta));
    if (pep->which==PEP_ALL && !ctx->hyperbolic) PetscCall(PetscViewerASCIIPrintf(viewer,"  checking eigenvalue type: %s\n",ctx->checket?"enabled":"disabled"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPReset_STOAR(PEP pep)
{
  PetscFunctionBegin;
  if (pep->which==PEP_ALL) PetscCall(PEPReset_STOAR_QSlice(pep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPDestroy_STOAR(PEP pep)
{
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(PetscFree(pep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetDetectZeros_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetDetectZeros_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetInertias_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetDimensions_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetDimensions_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetLinearization_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetLinearization_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetCheckEigenvalueType_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetCheckEigenvalueType_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_STOAR(PEP pep)
{
  PEP_STOAR      *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  pep->data = (void*)ctx;

  pep->lineariz = PETSC_TRUE;
  ctx->lock     = PETSC_TRUE;
  ctx->nev      = 1;
  ctx->ncv      = PETSC_DETERMINE;
  ctx->mpd      = PETSC_DETERMINE;
  ctx->alpha    = 1.0;
  ctx->beta     = 0.0;
  ctx->checket  = PETSC_TRUE;

  pep->ops->setup          = PEPSetUp_STOAR;
  pep->ops->setfromoptions = PEPSetFromOptions_STOAR;
  pep->ops->destroy        = PEPDestroy_STOAR;
  pep->ops->view           = PEPView_STOAR;
  pep->ops->backtransform  = PEPBackTransform_Default;
  pep->ops->computevectors = PEPComputeVectors_Default;
  pep->ops->extractvectors = PEPExtractVectors_TOAR;
  pep->ops->reset          = PEPReset_STOAR;

  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetLocking_C",PEPSTOARSetLocking_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetLocking_C",PEPSTOARGetLocking_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetDetectZeros_C",PEPSTOARSetDetectZeros_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetDetectZeros_C",PEPSTOARGetDetectZeros_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetInertias_C",PEPSTOARGetInertias_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetDimensions_C",PEPSTOARGetDimensions_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetDimensions_C",PEPSTOARSetDimensions_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetLinearization_C",PEPSTOARSetLinearization_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetLinearization_C",PEPSTOARGetLinearization_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetCheckEigenvalueType_C",PEPSTOARSetCheckEigenvalueType_STOAR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetCheckEigenvalueType_C",PEPSTOARGetCheckEigenvalueType_STOAR));
  PetscFunctionReturn(PETSC_SUCCESS);
}
