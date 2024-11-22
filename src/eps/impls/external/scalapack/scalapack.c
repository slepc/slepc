/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to eigensolvers in ScaLAPACK.
*/

#include <slepc/private/epsimpl.h>    /*I "slepceps.h" I*/
#include <slepc/private/slepcscalapack.h>

typedef struct {
  Mat As,Bs;        /* converted matrices */
} EPS_ScaLAPACK;

static PetscErrorCode EPSSetUp_ScaLAPACK(EPS eps)
{
  EPS_ScaLAPACK  *ctx = (EPS_ScaLAPACK*)eps->data;
  Mat            A,B;
  PetscInt       nmat;
  PetscBool      isshift;
  PetscScalar    shift;

  PetscFunctionBegin;
  EPSCheckHermitianDefinite(eps);
  EPSCheckNotStructured(eps);
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isshift));
  PetscCheck(isshift,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support spectral transformations");
  if (eps->nev==0) eps->nev = 1;
  eps->ncv = eps->n;
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = 1;
  if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  PetscCheck(eps->which!=EPS_ALL || eps->inta==eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support interval computation");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_STOPPING);
  PetscCall(EPSAllocateSolution(eps,0));

  /* convert matrices */
  PetscCall(MatDestroy(&ctx->As));
  PetscCall(MatDestroy(&ctx->Bs));
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  PetscCall(MatConvert(A,MATSCALAPACK,MAT_INITIAL_MATRIX,&ctx->As));
  if (nmat>1) {
    PetscCall(STGetMatrix(eps->st,1,&B));
    PetscCall(MatConvert(B,MATSCALAPACK,MAT_INITIAL_MATRIX,&ctx->Bs));
  }
  PetscCall(STGetShift(eps->st,&shift));
  if (shift != 0.0) {
    if (nmat>1) PetscCall(MatAXPY(ctx->As,-shift,ctx->Bs,SAME_NONZERO_PATTERN));
    else PetscCall(MatShift(ctx->As,-shift));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_ScaLAPACK(EPS eps)
{
  EPS_ScaLAPACK  *ctx = (EPS_ScaLAPACK*)eps->data;
  Mat            A = ctx->As,B = ctx->Bs,Q,V;
  Mat_ScaLAPACK  *a = (Mat_ScaLAPACK*)A->data,*b,*q;
  PetscReal      rdummy=0.0,abstol=0.0,*gap=NULL,orfac=-1.0,*w = eps->errest;  /* used to store real eigenvalues */
  PetscScalar    *work,minlwork[3];
  PetscBLASInt   i,m,info,idummy=0,lwork=-1,liwork=-1,minliwork,*iwork,*ifail=NULL,*iclustr=NULL,one=1;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork,minlrwork[3];
  PetscBLASInt   lrwork=-1;
#endif

  PetscFunctionBegin;
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&Q));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  q = (Mat_ScaLAPACK*)Q->data;

  if (B) {

    b = (Mat_ScaLAPACK*)B->data;
    PetscCall(PetscMalloc3(a->grid->nprow*a->grid->npcol,&gap,a->N,&ifail,2*a->grid->nprow*a->grid->npcol,&iclustr));
#if !defined(PETSC_USE_COMPLEX)
    /* allocate workspace */
    PetscCallBLAS("SCALAPACKsygvx",SCALAPACKsygvx_(&one,"V","A","L",&a->N,a->loc,&one,&one,a->desc,b->loc,&one,&one,b->desc,&rdummy,&rdummy,&idummy,&idummy,&abstol,&m,&idummy,w,&orfac,q->loc,&one,&one,q->desc,minlwork,&lwork,&minliwork,&liwork,ifail,iclustr,gap,&info));
    PetscCheckScaLapackInfo("sygvx",info);
    PetscCall(PetscBLASIntCast((PetscInt)minlwork[0],&lwork));
    liwork = minliwork;
    /* call computational routine */
    PetscCall(PetscMalloc2(lwork,&work,liwork,&iwork));
    PetscCallBLAS("SCALAPACKsygvx",SCALAPACKsygvx_(&one,"V","A","L",&a->N,a->loc,&one,&one,a->desc,b->loc,&one,&one,b->desc,&rdummy,&rdummy,&idummy,&idummy,&abstol,&m,&idummy,w,&orfac,q->loc,&one,&one,q->desc,work,&lwork,iwork,&liwork,ifail,iclustr,gap,&info));
    PetscCheckScaLapackInfo("sygvx",info);
    PetscCall(PetscFree2(work,iwork));
#else
    /* allocate workspace */
    PetscCallBLAS("SCALAPACKsygvx",SCALAPACKsygvx_(&one,"V","A","L",&a->N,a->loc,&one,&one,a->desc,b->loc,&one,&one,b->desc,&rdummy,&rdummy,&idummy,&idummy,&abstol,&m,&idummy,w,&orfac,q->loc,&one,&one,q->desc,minlwork,&lwork,minlrwork,&lrwork,&minliwork,&liwork,ifail,iclustr,gap,&info));
    PetscCheckScaLapackInfo("sygvx",info);
    PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(minlwork[0]),&lwork));
    PetscCall(PetscBLASIntCast((PetscInt)minlrwork[0],&lrwork));
    lrwork += a->N*a->N;
    liwork = minliwork;
    /* call computational routine */
    PetscCall(PetscMalloc3(lwork,&work,lrwork,&rwork,liwork,&iwork));
    PetscCallBLAS("SCALAPACKsygvx",SCALAPACKsygvx_(&one,"V","A","L",&a->N,a->loc,&one,&one,a->desc,b->loc,&one,&one,b->desc,&rdummy,&rdummy,&idummy,&idummy,&abstol,&m,&idummy,w,&orfac,q->loc,&one,&one,q->desc,work,&lwork,rwork,&lrwork,iwork,&liwork,ifail,iclustr,gap,&info));
    PetscCheckScaLapackInfo("sygvx",info);
    PetscCall(PetscFree3(work,rwork,iwork));
#endif
    PetscCall(PetscFree3(gap,ifail,iclustr));

  } else {

#if !defined(PETSC_USE_COMPLEX)
    /* allocate workspace */
    PetscCallBLAS("SCALAPACKsyev",SCALAPACKsyev_("V","L",&a->N,a->loc,&one,&one,a->desc,w,q->loc,&one,&one,q->desc,minlwork,&lwork,&info));
    PetscCheckScaLapackInfo("syev",info);
    PetscCall(PetscBLASIntCast((PetscInt)minlwork[0],&lwork));
    PetscCall(PetscMalloc1(lwork,&work));
    /* call computational routine */
    PetscCallBLAS("SCALAPACKsyev",SCALAPACKsyev_("V","L",&a->N,a->loc,&one,&one,a->desc,w,q->loc,&one,&one,q->desc,work,&lwork,&info));
    PetscCheckScaLapackInfo("syev",info);
    PetscCall(PetscFree(work));
#else
    /* allocate workspace */
    PetscCallBLAS("SCALAPACKsyev",SCALAPACKsyev_("V","L",&a->N,a->loc,&one,&one,a->desc,w,q->loc,&one,&one,q->desc,minlwork,&lwork,minlrwork,&lrwork,&info));
    PetscCheckScaLapackInfo("syev",info);
    PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(minlwork[0]),&lwork));
    lrwork = 4*a->N;  /* PetscCall(PetscBLASIntCast((PetscInt)minlrwork[0],&lrwork)); */
    PetscCall(PetscMalloc2(lwork,&work,lrwork,&rwork));
    /* call computational routine */
    PetscCallBLAS("SCALAPACKsyev",SCALAPACKsyev_("V","L",&a->N,a->loc,&one,&one,a->desc,w,q->loc,&one,&one,q->desc,work,&lwork,rwork,&lrwork,&info));
    PetscCheckScaLapackInfo("syev",info);
    PetscCall(PetscFree2(work,rwork));
#endif

  }
  PetscCall(PetscFPTrapPop());

  for (i=0;i<eps->ncv;i++) {
    eps->eigr[i]   = eps->errest[i];
    eps->errest[i] = PETSC_MACHINE_EPSILON;
  }

  PetscCall(BVGetMat(eps->V,&V));
  PetscCall(MatConvert(Q,MATDENSE,MAT_REUSE_MATRIX,&V));
  PetscCall(BVRestoreMat(eps->V,&V));
  PetscCall(MatDestroy(&Q));

  eps->nconv  = eps->ncv;
  eps->its    = 1;
  eps->reason = EPS_CONVERGED_TOL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_ScaLAPACK(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_ScaLAPACK(EPS eps)
{
  EPS_ScaLAPACK  *ctx = (EPS_ScaLAPACK*)eps->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->As));
  PetscCall(MatDestroy(&ctx->Bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_ScaLAPACK(EPS eps)
{
  EPS_ScaLAPACK  *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->categ = EPS_CATEGORY_OTHER;

  eps->ops->solve          = EPSSolve_ScaLAPACK;
  eps->ops->setup          = EPSSetUp_ScaLAPACK;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->destroy        = EPSDestroy_ScaLAPACK;
  eps->ops->reset          = EPSReset_ScaLAPACK;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_NoFactor;
  PetscFunctionReturn(PETSC_SUCCESS);
}
