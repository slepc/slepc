/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "ciss"

   Method: Contour Integral Spectral Slicing

   Algorithm:

       Contour integral based on Sakurai-Sugiura method to construct a
       subspace, with various eigenpair extractions (Rayleigh-Ritz,
       explicit moment).

   Based on code contributed by Y. Maeda, T. Sakurai.

   References:

       [1] J. Asakura, T. Sakurai, H. Tadano, T. Ikegami, K. Kimura, "A
           numerical method for nonlinear eigenvalue problems using contour
           integrals", JSIAM Lett. 1:52-55, 2009.

       [2] S. Yokota and T. Sakurai, "A projection method for nonlinear
           eigenvalue problems using contour integrals", JSIAM Lett.
           5:41-44, 2013.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <slepc/private/slepccontour.h>

typedef struct _n_nep_ciss_project *NEP_CISS_PROJECT;
typedef struct {
  /* parameters */
  PetscInt          N;             /* number of integration points (32) */
  PetscInt          L;             /* block size (16) */
  PetscInt          M;             /* moment degree (N/4 = 4) */
  PetscReal         delta;         /* threshold of singular value (1e-12) */
  PetscInt          L_max;         /* maximum number of columns of the source matrix V */
  PetscReal         spurious_threshold; /* discard spurious eigenpairs */
  PetscBool         isreal;        /* T(z) is real for real z */
  PetscInt          npart;         /* number of partitions */
  PetscInt          refine_inner;
  PetscInt          refine_blocksize;
  NEPCISSExtraction extraction;
  /* private data */
  SlepcContourData  contour;
  PetscReal         *sigma;        /* threshold for numerical rank */
  PetscScalar       *weight;
  PetscScalar       *omega;
  PetscScalar       *pp;
  BV                V;
  BV                S;
  BV                Y;
  PetscBool         useconj;
  Mat               J;             /* auxiliary matrix when using subcomm */
  BV                pV;
  NEP_CISS_PROJECT  dsctxf;
  PetscObjectId     rgid;
  PetscObjectState  rgstate;
} NEP_CISS;

struct _n_nep_ciss_project {
  NEP  nep;
  BV   Q;
};

static PetscErrorCode NEPContourDSComputeMatrix(DS ds,PetscScalar lambda,PetscBool deriv,DSMatType mat,void *ctx)
{
  NEP_CISS_PROJECT proj = (NEP_CISS_PROJECT)ctx;
  Mat              M,fun;

  PetscFunctionBegin;
  if (!deriv) {
    PetscCall(NEPComputeFunction(proj->nep,lambda,proj->nep->function,proj->nep->function));
    fun = proj->nep->function;
  } else {
    PetscCall(NEPComputeJacobian(proj->nep,lambda,proj->nep->jacobian));
    fun = proj->nep->jacobian;
  }
  PetscCall(DSGetMat(ds,mat,&M));
  PetscCall(BVMatProject(proj->Q,fun,proj->Q,M));
  PetscCall(DSRestoreMat(ds,mat,&M));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPComputeFunctionSubcomm(NEP nep,PetscScalar lambda,Mat T,Mat P,PetscBool deriv)
{
  PetscInt       i;
  PetscScalar    alpha;
  NEP_CISS       *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  PetscAssert(nep->fui!=NEP_USER_INTERFACE_CALLBACK,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONGSTATE,"Should not arrive here with callbacks");
  PetscCall(MatZeroEntries(T));
  if (!deriv && T != P) PetscCall(MatZeroEntries(P));
  for (i=0;i<nep->nt;i++) {
    if (!deriv) PetscCall(FNEvaluateFunction(nep->f[i],lambda,&alpha));
    else PetscCall(FNEvaluateDerivative(nep->f[i],lambda,&alpha));
    PetscCall(MatAXPY(T,alpha,ctx->contour->pA[i],nep->mstr));
    if (!deriv && T != P) PetscCall(MatAXPY(P,alpha,ctx->contour->pP[i],nep->mstrp));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Set up KSP solvers for every integration point
*/
static PetscErrorCode NEPCISSSetUp(NEP nep,Mat T,Mat P)
{
  NEP_CISS         *ctx = (NEP_CISS*)nep->data;
  SlepcContourData contour;
  PetscInt         i,p_id;
  Mat              Amat,Pmat;

  PetscFunctionBegin;
  if (!ctx->contour || !ctx->contour->ksp) PetscCall(NEPCISSGetKSPs(nep,NULL,NULL));
  contour = ctx->contour;
  PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)nep),PETSC_ERR_PLIB,"Something went wrong with NEPCISSGetKSPs()");
  for (i=0;i<contour->npoints;i++) {
    p_id = i*contour->subcomm->n + contour->subcomm->color;
    PetscCall(MatDuplicate(T,MAT_DO_NOT_COPY_VALUES,&Amat));
    if (T != P) PetscCall(MatDuplicate(P,MAT_DO_NOT_COPY_VALUES,&Pmat)); else Pmat = Amat;
    if (contour->subcomm->n == 1 || nep->fui==NEP_USER_INTERFACE_CALLBACK) PetscCall(NEPComputeFunction(nep,ctx->omega[p_id],Amat,Pmat));
    else PetscCall(NEPComputeFunctionSubcomm(nep,ctx->omega[p_id],Amat,Pmat,PETSC_FALSE));
    PetscCall(NEP_KSPSetOperators(contour->ksp[i],Amat,Pmat));
    PetscCall(MatDestroy(&Amat));
    if (T != P) PetscCall(MatDestroy(&Pmat));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Y_i = F(z_i)^{-1}Fp(z_i)V for every integration point, Y=[Y_i] is in the context
*/
static PetscErrorCode NEPCISSSolve(NEP nep,Mat dT,BV V,PetscInt L_start,PetscInt L_end)
{
  NEP_CISS         *ctx = (NEP_CISS*)nep->data;
  SlepcContourData contour = ctx->contour;
  PetscInt         i,p_id;
  Mat              MV,BMV=NULL,MC;

  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(V,L_start,L_end));
  PetscCall(BVGetMat(V,&MV));
  for (i=0;i<contour->npoints;i++) {
    p_id = i*contour->subcomm->n + contour->subcomm->color;
    if (contour->subcomm->n==1 || nep->fui==NEP_USER_INTERFACE_CALLBACK) PetscCall(NEPComputeJacobian(nep,ctx->omega[p_id],dT));
    else PetscCall(NEPComputeFunctionSubcomm(nep,ctx->omega[p_id],dT,NULL,PETSC_TRUE));
    PetscCall(BVSetActiveColumns(ctx->Y,i*ctx->L+L_start,i*ctx->L+L_end));
    PetscCall(BVGetMat(ctx->Y,&MC));
    if (!i) {
      PetscCall(MatProductCreate(dT,MV,NULL,&BMV));
      PetscCall(MatProductSetType(BMV,MATPRODUCT_AB));
      PetscCall(MatProductSetFromOptions(BMV));
      PetscCall(MatProductSymbolic(BMV));
    }
    PetscCall(MatProductNumeric(BMV));
    PetscCall(KSPMatSolve(contour->ksp[i],BMV,MC));
    PetscCall(BVRestoreMat(ctx->Y,&MC));
  }
  PetscCall(MatDestroy(&BMV));
  PetscCall(BVRestoreMat(V,&MV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSetUp_CISS(NEP nep)
{
  NEP_CISS         *ctx = (NEP_CISS*)nep->data;
  SlepcContourData contour;
  PetscInt         nwork;
  PetscBool        istrivial,isellipse,flg;
  NEP_CISS_PROJECT dsctxf;
  PetscObjectId    id;
  PetscObjectState state;
  Vec              v0;

  PetscFunctionBegin;
  if (nep->ncv==PETSC_DETERMINE) nep->ncv = ctx->L_max*ctx->M;
  else {
    ctx->L_max = nep->ncv/ctx->M;
    if (!ctx->L_max) {
      ctx->L_max = 1;
      nep->ncv = ctx->L_max*ctx->M;
    }
  }
  ctx->L = PetscMin(ctx->L,ctx->L_max);
  if (nep->max_it==PETSC_DETERMINE) nep->max_it = 5;
  if (nep->mpd==PETSC_DETERMINE) nep->mpd = nep->ncv;
  if (!nep->which) nep->which = NEP_ALL;
  PetscCheck(nep->which==NEP_ALL,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver supports only computing all eigenvalues");
  NEPCheckUnsupported(nep,NEP_FEATURE_STOPPING | NEP_FEATURE_TWOSIDED);

  /* check region */
  PetscCall(RGIsTrivial(nep->rg,&istrivial));
  PetscCheck(!istrivial,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"CISS requires a nontrivial region, e.g. -rg_type ellipse ...");
  PetscCall(RGGetComplement(nep->rg,&flg));
  PetscCheck(!flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"A region with complement flag set is not allowed");
  PetscCall(PetscObjectTypeCompare((PetscObject)nep->rg,RGELLIPSE,&isellipse));
  PetscCheck(isellipse,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Currently only implemented for elliptic regions");

  /* if the region has changed, then reset contour data */
  PetscCall(PetscObjectGetId((PetscObject)nep->rg,&id));
  PetscCall(PetscObjectStateGet((PetscObject)nep->rg,&state));
  if (ctx->rgid && (id != ctx->rgid || state != ctx->rgstate)) {
    PetscCall(SlepcContourDataDestroy(&ctx->contour));
    PetscCall(PetscInfo(nep,"Resetting the contour data structure due to a change of region\n"));
    ctx->rgid = id; ctx->rgstate = state;
  }

  /* create contour data structure */
  if (!ctx->contour) {
    PetscCall(RGCanUseConjugates(nep->rg,ctx->isreal,&ctx->useconj));
    PetscCall(SlepcContourDataCreate(ctx->useconj?ctx->N/2:ctx->N,ctx->npart,(PetscObject)nep,&ctx->contour));
  }

  PetscCall(NEPAllocateSolution(nep,0));
  if (ctx->weight) PetscCall(PetscFree4(ctx->weight,ctx->omega,ctx->pp,ctx->sigma));
  PetscCall(PetscMalloc4(ctx->N,&ctx->weight,ctx->N,&ctx->omega,ctx->N,&ctx->pp,ctx->L_max*ctx->M,&ctx->sigma));

  /* allocate basis vectors */
  PetscCall(BVDestroy(&ctx->S));
  PetscCall(BVDuplicateResize(nep->V,ctx->L*ctx->M,&ctx->S));
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVDuplicateResize(nep->V,ctx->L,&ctx->V));

  contour = ctx->contour;
  if (contour->subcomm && contour->subcomm->n != 1 && nep->fui==NEP_USER_INTERFACE_CALLBACK) {
    PetscCall(NEPComputeFunction(nep,0,nep->function,nep->function_pre));
    PetscCall(SlepcContourRedundantMat(contour,1,&nep->function,(nep->function!=nep->function_pre)?&nep->function_pre:NULL));
  } else PetscCall(SlepcContourRedundantMat(contour,nep->nt,nep->A,nep->P));
  if (contour->pA) {
    if (!ctx->J) PetscCall(MatDuplicate(contour->pA[0],MAT_DO_NOT_COPY_VALUES,&ctx->J));
    PetscCall(BVGetColumn(ctx->V,0,&v0));
    PetscCall(SlepcContourScatterCreate(contour,v0));
    PetscCall(BVRestoreColumn(ctx->V,0,&v0));
    PetscCall(BVDestroy(&ctx->pV));
    PetscCall(BVCreate(PetscObjectComm((PetscObject)contour->xsub),&ctx->pV));
    PetscCall(BVSetSizesFromVec(ctx->pV,contour->xsub,nep->n));
    PetscCall(BVSetFromOptions(ctx->pV));
    PetscCall(BVResize(ctx->pV,ctx->L,PETSC_FALSE));
  }

  PetscCall(BVDestroy(&ctx->Y));
  if (contour->pA) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)contour->xsub),&ctx->Y));
    PetscCall(BVSetSizesFromVec(ctx->Y,contour->xsub,nep->n));
    PetscCall(BVSetFromOptions(ctx->Y));
    PetscCall(BVResize(ctx->Y,contour->npoints*ctx->L,PETSC_FALSE));
  } else PetscCall(BVDuplicateResize(nep->V,contour->npoints*ctx->L,&ctx->Y));

  if (ctx->extraction == NEP_CISS_EXTRACTION_RITZ)  {
    PetscCall(DSSetMethod(nep->ds,1));
    PetscCall(DSNEPSetRG(nep->ds,nep->rg));
    if (nep->fui==NEP_USER_INTERFACE_SPLIT) PetscCall(DSNEPSetFN(nep->ds,nep->nt,nep->f));
    else {
      PetscCall(PetscNew(&dsctxf));
      PetscCall(DSNEPSetComputeMatrixFunction(nep->ds,NEPContourDSComputeMatrix,dsctxf));
      dsctxf->nep = nep;
      ctx->dsctxf = dsctxf;
    }
  }
  PetscCall(DSAllocate(nep->ds,nep->ncv));
  nwork = (nep->fui==NEP_USER_INTERFACE_SPLIT)? 2: 1;
  PetscCall(NEPSetWorkVecs(nep,nwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSolve_CISS(NEP nep)
{
  NEP_CISS         *ctx = (NEP_CISS*)nep->data;
  SlepcContourData contour = ctx->contour;
  Mat              X,M,E,T,P,J;
  BV               V;
  PetscInt         i,j,ld,L_add=0,nv=0,L_base=ctx->L,inner,*inside;
  PetscScalar      *Mu,*H0,*H1,*rr,*temp,center;
  PetscReal        error,max_error,radius,rgscale,est_eig,eta;
  PetscBool        isellipse,*fl1;
  Vec              si;
  SlepcSC          sc;
  PetscRandom      rand;

  PetscFunctionBegin;
  PetscCall(DSGetSlepcSC(nep->ds,&sc));
  sc->comparison    = SlepcCompareLargestMagnitude;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSGetLeadingDimension(nep->ds,&ld));
  PetscCall(RGComputeQuadrature(nep->rg,RG_QUADRULE_TRAPEZOIDAL,ctx->N,ctx->omega,ctx->pp,ctx->weight));
  if (contour->pA) {
    T = contour->pA[0];
    P = ((nep->fui==NEP_USER_INTERFACE_SPLIT && nep->P) || (nep->fui==NEP_USER_INTERFACE_CALLBACK && contour->pP))? contour->pP[0]: T;
  } else {
    T = nep->function;
    P = nep->function_pre? nep->function_pre: nep->function;
  }
  PetscCall(NEPCISSSetUp(nep,T,P));
  PetscCall(BVSetActiveColumns(ctx->V,0,ctx->L));
  PetscCall(BVSetRandomSign(ctx->V));
  PetscCall(BVGetRandomContext(ctx->V,&rand));
  if (contour->pA) {
    J = ctx->J;
    V = ctx->pV;
    PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
  } else {
    J = nep->jacobian;
    V = ctx->V;
  }
  PetscCall(NEPCISSSolve(nep,J,V,0,ctx->L));
  PetscCall(PetscObjectTypeCompare((PetscObject)nep->rg,RGELLIPSE,&isellipse));
  if (isellipse) {
    PetscCall(BVTraceQuadrature(ctx->Y,ctx->V,ctx->L,ctx->L,ctx->weight,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj,&est_eig));
    PetscCall(PetscInfo(nep,"Estimated eigenvalue count: %f\n",(double)est_eig));
    eta = PetscPowReal(10.0,-PetscLog10Real(nep->tol)/ctx->N);
    L_add = PetscMax(0,(PetscInt)PetscCeilReal((est_eig*eta)/ctx->M)-ctx->L);
    if (L_add>ctx->L_max-ctx->L) {
      PetscCall(PetscInfo(nep,"Number of eigenvalues inside the contour path may be too large\n"));
      L_add = ctx->L_max-ctx->L;
    }
  }
  /* Updates L after estimate the number of eigenvalue */
  if (L_add>0) {
    PetscCall(PetscInfo(nep,"Changing L %" PetscInt_FMT " -> %" PetscInt_FMT " by Estimate #Eig\n",ctx->L,ctx->L+L_add));
    PetscCall(BVCISSResizeBases(ctx->S,contour->pA?ctx->pV:ctx->V,ctx->Y,ctx->L,ctx->L+L_add,ctx->M,contour->npoints));
    PetscCall(BVSetActiveColumns(ctx->V,ctx->L,ctx->L+L_add));
    PetscCall(BVSetRandomSign(ctx->V));
    if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
    ctx->L += L_add;
    PetscCall(NEPCISSSolve(nep,J,V,ctx->L-L_add,ctx->L));
  }

  PetscCall(PetscMalloc2(ctx->L*ctx->L*ctx->M*2,&Mu,ctx->L*ctx->M*ctx->L*ctx->M,&H0));
  for (i=0;i<ctx->refine_blocksize;i++) {
    PetscCall(BVDotQuadrature(ctx->Y,V,Mu,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->subcomm,contour->npoints,ctx->useconj));
    PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
    PetscCall(PetscLogEventBegin(NEP_CISS_SVD,nep,0,0,0));
    PetscCall(SlepcCISS_BH_SVD(H0,ctx->L*ctx->M,ctx->delta,ctx->sigma,&nv));
    PetscCall(PetscLogEventEnd(NEP_CISS_SVD,nep,0,0,0));
    if (ctx->sigma[0]<=ctx->delta || nv < ctx->L*ctx->M || ctx->L == ctx->L_max) break;
    L_add = L_base;
    if (ctx->L+L_add>ctx->L_max) L_add = ctx->L_max-ctx->L;
    PetscCall(PetscInfo(nep,"Changing L %" PetscInt_FMT " -> %" PetscInt_FMT " by SVD(H0)\n",ctx->L,ctx->L+L_add));
    PetscCall(BVCISSResizeBases(ctx->S,contour->pA?ctx->pV:ctx->V,ctx->Y,ctx->L,ctx->L+L_add,ctx->M,contour->npoints));
    PetscCall(BVSetActiveColumns(ctx->V,ctx->L,ctx->L+L_add));
    PetscCall(BVSetRandomSign(ctx->V));
    if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
    ctx->L += L_add;
    PetscCall(NEPCISSSolve(nep,J,V,ctx->L-L_add,ctx->L));
    if (L_add) {
      PetscCall(PetscFree2(Mu,H0));
      PetscCall(PetscMalloc2(ctx->L*ctx->L*ctx->M*2,&Mu,ctx->L*ctx->M*ctx->L*ctx->M,&H0));
    }
  }

  PetscCall(RGGetScale(nep->rg,&rgscale));
  PetscCall(RGEllipseGetParameters(nep->rg,&center,&radius,NULL));

  if (ctx->extraction == NEP_CISS_EXTRACTION_HANKEL) PetscCall(PetscMalloc1(ctx->L*ctx->M*ctx->L*ctx->M,&H1));

  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;
    for (inner=0;inner<=ctx->refine_inner;inner++) {
      if (ctx->extraction == NEP_CISS_EXTRACTION_HANKEL) {
        PetscCall(BVDotQuadrature(ctx->Y,V,Mu,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
        PetscCall(PetscLogEventBegin(NEP_CISS_SVD,nep,0,0,0));
        PetscCall(SlepcCISS_BH_SVD(H0,ctx->L*ctx->M,ctx->delta,ctx->sigma,&nv));
        PetscCall(PetscLogEventEnd(NEP_CISS_SVD,nep,0,0,0));
      } else {
        PetscCall(BVSumQuadrature(ctx->S,ctx->Y,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj));
        /* compute SVD of S */
        PetscCall(BVSVDAndRank(ctx->S,ctx->M,ctx->L,ctx->delta,(ctx->extraction==NEP_CISS_EXTRACTION_CAA)?BV_SVD_METHOD_QR_CAA:BV_SVD_METHOD_QR,H0,ctx->sigma,&nv));
      }
      PetscCall(PetscInfo(nep,"Estimated rank: nv = %" PetscInt_FMT "\n",nv));
      if (ctx->sigma[0]>ctx->delta && nv==ctx->L*ctx->M && inner!=ctx->refine_inner) {
        PetscCall(BVSumQuadrature(ctx->S,ctx->Y,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(BVSetActiveColumns(ctx->S,0,ctx->L));
        PetscCall(BVSetActiveColumns(ctx->V,0,ctx->L));
        PetscCall(BVCopy(ctx->S,ctx->V));
        if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
        PetscCall(NEPCISSSolve(nep,J,V,0,ctx->L));
      } else break;
    }
    nep->nconv = 0;
    if (nv == 0) { nep->reason = NEP_CONVERGED_TOL; break; }
    else {
      /* Extracting eigenpairs */
      PetscCall(DSSetDimensions(nep->ds,nv,0,0));
      PetscCall(DSSetState(nep->ds,DS_STATE_RAW));
      if (ctx->extraction == NEP_CISS_EXTRACTION_HANKEL) {
        PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
        PetscCall(CISS_BlockHankel(Mu,1,ctx->L,ctx->M,H1));
        PetscCall(DSGetArray(nep->ds,DS_MAT_A,&temp));
        for (j=0;j<nv;j++)
          for (i=0;i<nv;i++)
            temp[i+j*ld] = H1[i+j*ctx->L*ctx->M];
        PetscCall(DSRestoreArray(nep->ds,DS_MAT_A,&temp));
        PetscCall(DSGetArray(nep->ds,DS_MAT_B,&temp));
        for (j=0;j<nv;j++)
          for (i=0;i<nv;i++)
            temp[i+j*ld] = H0[i+j*ctx->L*ctx->M];
        PetscCall(DSRestoreArray(nep->ds,DS_MAT_B,&temp));
      } else if (ctx->extraction == NEP_CISS_EXTRACTION_CAA) {
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
        PetscCall(DSGetArray(nep->ds,DS_MAT_A,&temp));
        for (i=0;i<nv;i++) PetscCall(PetscArraycpy(temp+i*ld,H0+i*nv,nv));
        PetscCall(DSRestoreArray(nep->ds,DS_MAT_A,&temp));
      } else {
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
        if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
          for (i=0;i<nep->nt;i++) {
            PetscCall(DSGetMat(nep->ds,DSMatExtra[i],&E));
            PetscCall(BVMatProject(ctx->S,nep->A[i],ctx->S,E));
            PetscCall(DSRestoreMat(nep->ds,DSMatExtra[i],&E));
          }
        } else { ctx->dsctxf->Q = ctx->S; }
      }
      PetscCall(DSSolve(nep->ds,nep->eigr,nep->eigi));
      PetscCall(DSSynchronize(nep->ds,nep->eigr,nep->eigi));
      PetscCall(DSGetDimensions(nep->ds,NULL,NULL,NULL,&nv));
      if (ctx->extraction == NEP_CISS_EXTRACTION_CAA || ctx->extraction == NEP_CISS_EXTRACTION_HANKEL) {
        for (i=0;i<nv;i++) {
          nep->eigr[i] = (nep->eigr[i]*radius+center)*rgscale;
        }
      }
      PetscCall(PetscMalloc3(nv,&fl1,nv,&inside,nv,&rr));
      PetscCall(DSVectors(nep->ds,DS_MAT_X,NULL,NULL));
      PetscCall(DSGetMat(nep->ds,DS_MAT_X,&X));
      PetscCall(SlepcCISS_isGhost(X,nv,ctx->sigma,ctx->spurious_threshold,fl1));
      PetscCall(DSRestoreMat(nep->ds,DS_MAT_X,&X));
      PetscCall(RGCheckInside(nep->rg,nv,nep->eigr,nep->eigi,inside));
      for (i=0;i<nv;i++) {
        if (fl1[i] && inside[i]>=0) {
          rr[i] = 1.0;
          nep->nconv++;
        } else rr[i] = 0.0;
      }
      PetscCall(DSSort(nep->ds,nep->eigr,nep->eigi,rr,NULL,&nep->nconv));
      PetscCall(DSSynchronize(nep->ds,nep->eigr,nep->eigi));
      if (ctx->extraction == NEP_CISS_EXTRACTION_CAA || ctx->extraction == NEP_CISS_EXTRACTION_HANKEL) {
        for (i=0;i<nv;i++) nep->eigr[i] = (nep->eigr[i]*radius+center)*rgscale;
      }
      PetscCall(PetscFree3(fl1,inside,rr));
      PetscCall(BVSetActiveColumns(nep->V,0,nv));
      PetscCall(DSVectors(nep->ds,DS_MAT_X,NULL,NULL));
      if (ctx->extraction == NEP_CISS_EXTRACTION_HANKEL) {
        PetscCall(BVSumQuadrature(ctx->S,ctx->Y,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
        PetscCall(BVCopy(ctx->S,nep->V));
        PetscCall(DSGetMat(nep->ds,DS_MAT_X,&X));
        PetscCall(BVMultInPlace(ctx->S,X,0,nep->nconv));
        PetscCall(BVMultInPlace(nep->V,X,0,nep->nconv));
        PetscCall(DSRestoreMat(nep->ds,DS_MAT_X,&X));
      } else {
        PetscCall(DSGetMat(nep->ds,DS_MAT_X,&X));
        PetscCall(BVMultInPlace(ctx->S,X,0,nep->nconv));
        PetscCall(DSRestoreMat(nep->ds,DS_MAT_X,&X));
        PetscCall(BVSetActiveColumns(ctx->S,0,nep->nconv));
        PetscCall(BVCopy(ctx->S,nep->V));
      }
      max_error = 0.0;
      for (i=0;i<nep->nconv;i++) {
        PetscCall(BVGetColumn(nep->V,i,&si));
        PetscCall(VecNormalize(si,NULL));
        PetscCall(NEPComputeResidualNorm_Private(nep,PETSC_FALSE,nep->eigr[i],si,nep->work,&error));
        PetscCall((*nep->converged)(nep,nep->eigr[i],0,error,&error,nep->convergedctx));
        PetscCall(BVRestoreColumn(nep->V,i,&si));
        max_error = PetscMax(max_error,error);
      }
      if (max_error <= nep->tol) nep->reason = NEP_CONVERGED_TOL;
      else if (nep->its > nep->max_it) nep->reason = NEP_DIVERGED_ITS;
      else {
        if (nep->nconv > ctx->L) nv = nep->nconv;
        else if (ctx->L > nv) nv = ctx->L;
        nv = PetscMin(nv,ctx->L*ctx->M);
        PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,nv,ctx->L,NULL,&M));
        PetscCall(MatSetRandom(M,rand));
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
        PetscCall(BVMultInPlace(ctx->S,M,0,ctx->L));
        PetscCall(MatDestroy(&M));
        PetscCall(BVSetActiveColumns(ctx->S,0,ctx->L));
        PetscCall(BVSetActiveColumns(ctx->V,0,ctx->L));
        PetscCall(BVCopy(ctx->S,ctx->V));
        if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
        PetscCall(NEPCISSSolve(nep,J,V,0,ctx->L));
      }
    }
  }
  PetscCall(PetscFree2(Mu,H0));
  if (ctx->extraction == NEP_CISS_EXTRACTION_HANKEL) PetscCall(PetscFree(H1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSSetSizes_CISS(NEP nep,PetscInt ip,PetscInt bs,PetscInt ms,PetscInt npart,PetscInt bsmax,PetscBool realmats)
{
  NEP_CISS       *ctx = (NEP_CISS*)nep->data;
  PetscInt       oN,oL,oM,oLmax,onpart;
  PetscMPIInt    size;

  PetscFunctionBegin;
  oN = ctx->N;
  if (ip == PETSC_DETERMINE) {
    if (ctx->N!=32) { ctx->N =32; ctx->M = ctx->N/4; }
  } else if (ip != PETSC_CURRENT) {
    PetscCheck(ip>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The ip argument must be > 0");
    PetscCheck(ip%2==0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The ip argument must be an even number");
    if (ctx->N!=ip) { ctx->N = ip; ctx->M = ctx->N/4; }
  }
  oL = ctx->L;
  if (bs == PETSC_DETERMINE) {
    ctx->L = 16;
  } else if (bs != PETSC_CURRENT) {
    PetscCheck(bs>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The bs argument must be > 0");
    ctx->L = bs;
  }
  oM = ctx->M;
  if (ms == PETSC_DETERMINE) {
    ctx->M = ctx->N/4;
  } else if (ms != PETSC_CURRENT) {
    PetscCheck(ms>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The ms argument must be > 0");
    PetscCheck(ms<=ctx->N,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The ms argument must be less than or equal to the number of integration points");
    ctx->M = ms;
  }
  onpart = ctx->npart;
  if (npart == PETSC_DETERMINE) {
    ctx->npart = 1;
  } else if (npart != PETSC_CURRENT) {
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)nep),&size));
    PetscCheck(npart>0 && npart<=size,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of npart");
    ctx->npart = npart;
  }
  oLmax = ctx->L_max;
  if (bsmax == PETSC_DETERMINE) {
    ctx->L_max = 64;
  } else if (bsmax != PETSC_CURRENT) {
    PetscCheck(bsmax>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The bsmax argument must be > 0");
    ctx->L_max = PetscMax(bsmax,ctx->L);
  }
  if (onpart != ctx->npart || oN != ctx->N || realmats != ctx->isreal) {
    PetscCall(SlepcContourDataDestroy(&ctx->contour));
    PetscCall(PetscInfo(nep,"Resetting the contour data structure due to a change of parameters\n"));
    nep->state = NEP_STATE_INITIAL;
  }
  ctx->isreal = realmats;
  if (oL != ctx->L || oM != ctx->M || oLmax != ctx->L_max) nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPCISSSetSizes - Sets the values of various size parameters in the CISS solver.

   Logically Collective

   Input Parameters:
+  nep   - the nonlinear eigensolver context
.  ip    - number of integration points
.  bs    - block size
.  ms    - moment size
.  npart - number of partitions when splitting the communicator
.  bsmax - max block size
-  realmats - T(z) is real for real z

   Options Database Keys:
+  -nep_ciss_integration_points - Sets the number of integration points
.  -nep_ciss_blocksize - Sets the block size
.  -nep_ciss_moments - Sets the moment size
.  -nep_ciss_partitions - Sets the number of partitions
.  -nep_ciss_maxblocksize - Sets the maximum block size
-  -nep_ciss_realmats - T(z) is real for real z

   Notes:
   For all integer arguments, you can use PETSC_CURRENT to keep the current value, and
   PETSC_DETERMINE to set them to a default value.

   The default number of partitions is 1. This means the internal KSP object is shared
   among all processes of the NEP communicator. Otherwise, the communicator is split
   into npart communicators, so that npart KSP solves proceed simultaneously.

   The realmats flag can be set to true when T(.) is guaranteed to be real
   when the argument is a real value, for example, when all matrices in
   the split form are real. When set to true, the solver avoids some computations.

   Level: advanced

.seealso: NEPCISSGetSizes()
@*/
PetscErrorCode NEPCISSSetSizes(NEP nep,PetscInt ip,PetscInt bs,PetscInt ms,PetscInt npart,PetscInt bsmax,PetscBool realmats)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,ip,2);
  PetscValidLogicalCollectiveInt(nep,bs,3);
  PetscValidLogicalCollectiveInt(nep,ms,4);
  PetscValidLogicalCollectiveInt(nep,npart,5);
  PetscValidLogicalCollectiveInt(nep,bsmax,6);
  PetscValidLogicalCollectiveBool(nep,realmats,7);
  PetscTryMethod(nep,"NEPCISSSetSizes_C",(NEP,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool),(nep,ip,bs,ms,npart,bsmax,realmats));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSGetSizes_CISS(NEP nep,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats)
{
  NEP_CISS *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  if (ip) *ip = ctx->N;
  if (bs) *bs = ctx->L;
  if (ms) *ms = ctx->M;
  if (npart) *npart = ctx->npart;
  if (bsmax) *bsmax = ctx->L_max;
  if (realmats) *realmats = ctx->isreal;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPCISSGetSizes - Gets the values of various size parameters in the CISS solver.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  ip    - number of integration points
.  bs    - block size
.  ms    - moment size
.  npart - number of partitions when splitting the communicator
.  bsmax - max block size
-  realmats - T(z) is real for real z

   Level: advanced

.seealso: NEPCISSSetSizes()
@*/
PetscErrorCode NEPCISSGetSizes(NEP nep,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscUseMethod(nep,"NEPCISSGetSizes_C",(NEP,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*),(nep,ip,bs,ms,npart,bsmax,realmats));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSSetThreshold_CISS(NEP nep,PetscReal delta,PetscReal spur)
{
  NEP_CISS *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  if (delta == (PetscReal)PETSC_DETERMINE) {
    ctx->delta = SLEPC_DEFAULT_TOL*1e-4;
  } else if (delta != (PetscReal)PETSC_CURRENT) {
    PetscCheck(delta>0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The delta argument must be > 0.0");
    ctx->delta = delta;
  }
  if (spur == (PetscReal)PETSC_DETERMINE) {
    ctx->spurious_threshold = PetscSqrtReal(SLEPC_DEFAULT_TOL);
  } else if (spur != (PetscReal)PETSC_CURRENT) {
    PetscCheck(spur>0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The spurious threshold argument must be > 0.0");
    ctx->spurious_threshold = spur;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPCISSSetThreshold - Sets the values of various threshold parameters in
   the CISS solver.

   Logically Collective

   Input Parameters:
+  nep   - the nonlinear eigensolver context
.  delta - threshold for numerical rank
-  spur  - spurious threshold (to discard spurious eigenpairs)

   Options Database Keys:
+  -nep_ciss_delta - Sets the delta
-  -nep_ciss_spurious_threshold - Sets the spurious threshold

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: NEPCISSGetThreshold()
@*/
PetscErrorCode NEPCISSSetThreshold(NEP nep,PetscReal delta,PetscReal spur)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,delta,2);
  PetscValidLogicalCollectiveReal(nep,spur,3);
  PetscTryMethod(nep,"NEPCISSSetThreshold_C",(NEP,PetscReal,PetscReal),(nep,delta,spur));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSGetThreshold_CISS(NEP nep,PetscReal *delta,PetscReal *spur)
{
  NEP_CISS *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  if (delta) *delta = ctx->delta;
  if (spur)  *spur = ctx->spurious_threshold;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPCISSGetThreshold - Gets the values of various threshold parameters in
   the CISS solver.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  delta - threshold for numerical rank
-  spur  - spurious threshold (to discard spurious eigenpairs)

   Level: advanced

.seealso: NEPCISSSetThreshold()
@*/
PetscErrorCode NEPCISSGetThreshold(NEP nep,PetscReal *delta,PetscReal *spur)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscUseMethod(nep,"NEPCISSGetThreshold_C",(NEP,PetscReal*,PetscReal*),(nep,delta,spur));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSSetRefinement_CISS(NEP nep,PetscInt inner,PetscInt blsize)
{
  NEP_CISS *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  if (inner == PETSC_DETERMINE) {
    ctx->refine_inner = 0;
  } else if (inner != PETSC_CURRENT) {
    PetscCheck(inner>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The refine inner argument must be >= 0");
    ctx->refine_inner = inner;
  }
  if (blsize == PETSC_DETERMINE) {
    ctx->refine_blocksize = 0;
  } else if (blsize != PETSC_CURRENT) {
    PetscCheck(blsize>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The refine blocksize argument must be >= 0");
    ctx->refine_blocksize = blsize;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPCISSSetRefinement - Sets the values of various refinement parameters
   in the CISS solver.

   Logically Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  inner  - number of iterative refinement iterations (inner loop)
-  blsize - number of iterative refinement iterations (blocksize loop)

   Options Database Keys:
+  -nep_ciss_refine_inner - Sets number of inner iterations
-  -nep_ciss_refine_blocksize - Sets number of blocksize iterations

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: NEPCISSGetRefinement()
@*/
PetscErrorCode NEPCISSSetRefinement(NEP nep,PetscInt inner,PetscInt blsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,inner,2);
  PetscValidLogicalCollectiveInt(nep,blsize,3);
  PetscTryMethod(nep,"NEPCISSSetRefinement_C",(NEP,PetscInt,PetscInt),(nep,inner,blsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSGetRefinement_CISS(NEP nep,PetscInt *inner,PetscInt *blsize)
{
  NEP_CISS *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  if (inner)  *inner = ctx->refine_inner;
  if (blsize) *blsize = ctx->refine_blocksize;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPCISSGetRefinement - Gets the values of various refinement parameters
   in the CISS solver.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  inner  - number of iterative refinement iterations (inner loop)
-  blsize - number of iterative refinement iterations (blocksize loop)

   Level: advanced

.seealso: NEPCISSSetRefinement()
@*/
PetscErrorCode NEPCISSGetRefinement(NEP nep, PetscInt *inner, PetscInt *blsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscUseMethod(nep,"NEPCISSGetRefinement_C",(NEP,PetscInt*,PetscInt*),(nep,inner,blsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSSetExtraction_CISS(NEP nep,NEPCISSExtraction extraction)
{
  NEP_CISS *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  if (ctx->extraction != extraction) {
    ctx->extraction = extraction;
    nep->state      = NEP_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPCISSSetExtraction - Sets the extraction technique used in the CISS solver.

   Logically Collective

   Input Parameters:
+  nep        - the nonlinear eigensolver context
-  extraction - the extraction technique

   Options Database Key:
.  -nep_ciss_extraction - Sets the extraction technique (either 'ritz', 'hankel' or 'caa')

   Notes:
   By default, the Rayleigh-Ritz extraction is used (NEP_CISS_EXTRACTION_RITZ).

   If the 'hankel' or the 'caa' option is specified (NEP_CISS_EXTRACTION_HANKEL or
   NEP_CISS_EXTRACTION_CAA), then the Block Hankel method, or the Communication-avoiding
   Arnoldi method, respectively, is used for extracting eigenpairs.

   Level: advanced

.seealso: NEPCISSGetExtraction(), NEPCISSExtraction
@*/
PetscErrorCode NEPCISSSetExtraction(NEP nep,NEPCISSExtraction extraction)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(nep,extraction,2);
  PetscTryMethod(nep,"NEPCISSSetExtraction_C",(NEP,NEPCISSExtraction),(nep,extraction));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSGetExtraction_CISS(NEP nep,NEPCISSExtraction *extraction)
{
  NEP_CISS *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  *extraction = ctx->extraction;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPCISSGetExtraction - Gets the extraction technique used in the CISS solver.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
.  extraction - extraction technique

   Level: advanced

.seealso: NEPCISSSetExtraction() NEPCISSExtraction
@*/
PetscErrorCode NEPCISSGetExtraction(NEP nep,NEPCISSExtraction *extraction)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(extraction,2);
  PetscUseMethod(nep,"NEPCISSGetExtraction_C",(NEP,NEPCISSExtraction*),(nep,extraction));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPCISSGetKSPs_CISS(NEP nep,PetscInt *nsolve,KSP **ksp)
{
  NEP_CISS         *ctx = (NEP_CISS*)nep->data;
  SlepcContourData contour;
  PetscInt         i;
  PC               pc;
  MPI_Comm         child;

  PetscFunctionBegin;
  if (!ctx->contour) {  /* initialize contour data structure first */
    PetscCall(RGCanUseConjugates(nep->rg,ctx->isreal,&ctx->useconj));
    PetscCall(SlepcContourDataCreate(ctx->useconj?ctx->N/2:ctx->N,ctx->npart,(PetscObject)nep,&ctx->contour));
  }
  contour = ctx->contour;
  if (!contour->ksp) {
    PetscCall(PetscMalloc1(contour->npoints,&contour->ksp));
    PetscCall(PetscSubcommGetChild(contour->subcomm,&child));
    for (i=0;i<contour->npoints;i++) {
      PetscCall(KSPCreate(child,&contour->ksp[i]));
      PetscCall(PetscObjectIncrementTabLevel((PetscObject)contour->ksp[i],(PetscObject)nep,1));
      PetscCall(KSPSetOptionsPrefix(contour->ksp[i],((PetscObject)nep)->prefix));
      PetscCall(KSPAppendOptionsPrefix(contour->ksp[i],"nep_ciss_"));
      PetscCall(PetscObjectSetOptions((PetscObject)contour->ksp[i],((PetscObject)nep)->options));
      PetscCall(KSPSetErrorIfNotConverged(contour->ksp[i],PETSC_TRUE));
      PetscCall(KSPSetTolerances(contour->ksp[i],SlepcDefaultTol(nep->tol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
      PetscCall(KSPGetPC(contour->ksp[i],&pc));
      if ((nep->fui==NEP_USER_INTERFACE_SPLIT && nep->P) || (nep->fui==NEP_USER_INTERFACE_CALLBACK && nep->function_pre!=nep->function)) {
        PetscCall(KSPSetType(contour->ksp[i],KSPBCGS));
        PetscCall(PCSetType(pc,PCBJACOBI));
      } else {
        PetscCall(KSPSetType(contour->ksp[i],KSPPREONLY));
        PetscCall(PCSetType(pc,PCLU));
      }
    }
  }
  if (nsolve) *nsolve = contour->npoints;
  if (ksp)    *ksp    = contour->ksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPCISSGetKSPs - Retrieve the array of linear solver objects associated with
   the CISS solver.

   Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameters:
+  nsolve - number of solver objects
-  ksp - array of linear solver object

   Notes:
   The number of KSP solvers is equal to the number of integration points divided by
   the number of partitions. This value is halved in the case of real matrices with
   a region centered at the real axis.

   Level: advanced

.seealso: NEPCISSSetSizes()
@*/
PetscErrorCode NEPCISSGetKSPs(NEP nep,PetscInt *nsolve,KSP **ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscUseMethod(nep,"NEPCISSGetKSPs_C",(NEP,PetscInt*,KSP**),(nep,nsolve,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPReset_CISS(NEP nep)
{
  NEP_CISS       *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&ctx->S));
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVDestroy(&ctx->Y));
  PetscCall(SlepcContourDataReset(ctx->contour));
  PetscCall(MatDestroy(&ctx->J));
  PetscCall(BVDestroy(&ctx->pV));
  if (ctx->extraction == NEP_CISS_EXTRACTION_RITZ && nep->fui==NEP_USER_INTERFACE_CALLBACK) PetscCall(PetscFree(ctx->dsctxf));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSetFromOptions_CISS(NEP nep,PetscOptionItems PetscOptionsObject)
{
  NEP_CISS          *ctx = (NEP_CISS*)nep->data;
  PetscReal         r1,r2;
  PetscInt          i,i1,i2,i3,i4,i5,i6,i7;
  PetscBool         b1,flg,flg2,flg3,flg4,flg5,flg6;
  NEPCISSExtraction extraction;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"NEP CISS Options");

    PetscCall(NEPCISSGetSizes(nep,&i1,&i2,&i3,&i4,&i5,&b1));
    PetscCall(PetscOptionsInt("-nep_ciss_integration_points","Number of integration points","NEPCISSSetSizes",i1,&i1,&flg));
    PetscCall(PetscOptionsInt("-nep_ciss_blocksize","Block size","NEPCISSSetSizes",i2,&i2,&flg2));
    PetscCall(PetscOptionsInt("-nep_ciss_moments","Moment size","NEPCISSSetSizes",i3,&i3,&flg3));
    PetscCall(PetscOptionsInt("-nep_ciss_partitions","Number of partitions","NEPCISSSetSizes",i4,&i4,&flg4));
    PetscCall(PetscOptionsInt("-nep_ciss_maxblocksize","Maximum block size","NEPCISSSetSizes",i5,&i5,&flg5));
    PetscCall(PetscOptionsBool("-nep_ciss_realmats","True if T(z) is real for real z","NEPCISSSetSizes",b1,&b1,&flg6));
    if (flg || flg2 || flg3 || flg4 || flg5 || flg6) PetscCall(NEPCISSSetSizes(nep,i1,i2,i3,i4,i5,b1));

    PetscCall(NEPCISSGetThreshold(nep,&r1,&r2));
    PetscCall(PetscOptionsReal("-nep_ciss_delta","Threshold for numerical rank","NEPCISSSetThreshold",r1,&r1,&flg));
    PetscCall(PetscOptionsReal("-nep_ciss_spurious_threshold","Threshold for the spurious eigenpairs","NEPCISSSetThreshold",r2,&r2,&flg2));
    if (flg || flg2) PetscCall(NEPCISSSetThreshold(nep,r1,r2));

    PetscCall(NEPCISSGetRefinement(nep,&i6,&i7));
    PetscCall(PetscOptionsInt("-nep_ciss_refine_inner","Number of inner iterative refinement iterations","NEPCISSSetRefinement",i6,&i6,&flg));
    PetscCall(PetscOptionsInt("-nep_ciss_refine_blocksize","Number of blocksize iterative refinement iterations","NEPCISSSetRefinement",i7,&i7,&flg2));
    if (flg || flg2) PetscCall(NEPCISSSetRefinement(nep,i6,i7));

    PetscCall(PetscOptionsEnum("-nep_ciss_extraction","Extraction technique","NEPCISSSetExtraction",NEPCISSExtractions,(PetscEnum)ctx->extraction,(PetscEnum*)&extraction,&flg));
    if (flg) PetscCall(NEPCISSSetExtraction(nep,extraction));

  PetscOptionsHeadEnd();

  if (!nep->rg) PetscCall(NEPGetRG(nep,&nep->rg));
  PetscCall(RGSetFromOptions(nep->rg)); /* this is necessary here to set useconj */
  if (!ctx->contour || !ctx->contour->ksp) PetscCall(NEPCISSGetKSPs(nep,NULL,NULL));
  PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)nep),PETSC_ERR_PLIB,"Something went wrong with NEPCISSGetKSPs()");
  for (i=0;i<ctx->contour->npoints;i++) PetscCall(KSPSetFromOptions(ctx->contour->ksp[i]));
  PetscCall(PetscSubcommSetFromOptions(ctx->contour->subcomm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDestroy_CISS(NEP nep)
{
  NEP_CISS       *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  PetscCall(SlepcContourDataDestroy(&ctx->contour));
  PetscCall(PetscFree4(ctx->weight,ctx->omega,ctx->pp,ctx->sigma));
  PetscCall(PetscFree(nep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSSetSizes_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetSizes_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSSetThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSSetRefinement_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetRefinement_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSSetExtraction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetExtraction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetKSPs_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPView_CISS(NEP nep,PetscViewer viewer)
{
  NEP_CISS       *ctx = (NEP_CISS*)nep->data;
  PetscBool      isascii;
  PetscViewer    sviewer;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  sizes { integration points: %" PetscInt_FMT ", block size: %" PetscInt_FMT ", moment size: %" PetscInt_FMT ", partitions: %" PetscInt_FMT ", maximum block size: %" PetscInt_FMT " }\n",ctx->N,ctx->L,ctx->M,ctx->npart,ctx->L_max));
    if (ctx->isreal) PetscCall(PetscViewerASCIIPrintf(viewer,"  exploiting symmetry of integration points\n"));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  threshold { delta: %g, spurious threshold: %g }\n",(double)ctx->delta,(double)ctx->spurious_threshold));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  iterative refinement  { inner: %" PetscInt_FMT ", blocksize: %" PetscInt_FMT " }\n",ctx->refine_inner, ctx->refine_blocksize));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  extraction: %s\n",NEPCISSExtractions[ctx->extraction]));
    if (!ctx->contour || !ctx->contour->ksp) PetscCall(NEPCISSGetKSPs(nep,NULL,NULL));
    PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)nep),PETSC_ERR_PLIB,"Something went wrong with NEPCISSGetKSPs()");
    PetscCall(PetscViewerASCIIPushTab(viewer));
    if (ctx->npart>1 && ctx->contour->subcomm) {
      PetscCall(PetscViewerGetSubViewer(viewer,ctx->contour->subcomm->child,&sviewer));
      if (!ctx->contour->subcomm->color) PetscCall(KSPView(ctx->contour->ksp[0],sviewer));
      PetscCall(PetscViewerFlush(sviewer));
      PetscCall(PetscViewerRestoreSubViewer(viewer,ctx->contour->subcomm->child,&sviewer));
      /* extra call needed because of the two calls to PetscViewerASCIIPushSynchronized() in PetscViewerGetSubViewer() */
      PetscCall(PetscViewerASCIIPopSynchronized(viewer));
    } else PetscCall(KSPView(ctx->contour->ksp[0],viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSetDSType_CISS(NEP nep)
{
  NEP_CISS       *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  switch (ctx->extraction) {
    case NEP_CISS_EXTRACTION_RITZ:   PetscCall(DSSetType(nep->ds,DSNEP)); break;
    case NEP_CISS_EXTRACTION_HANKEL: PetscCall(DSSetType(nep->ds,DSGNHEP)); break;
    case NEP_CISS_EXTRACTION_CAA:    PetscCall(DSSetType(nep->ds,DSNHEP)); break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_CISS(NEP nep)
{
  NEP_CISS       *ctx = (NEP_CISS*)nep->data;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  nep->data = ctx;
  /* set default values of parameters */
  ctx->N                  = 32;
  ctx->L                  = 16;
  ctx->M                  = ctx->N/4;
  ctx->delta              = SLEPC_DEFAULT_TOL*1e-4;
  ctx->L_max              = 64;
  ctx->spurious_threshold = PetscSqrtReal(SLEPC_DEFAULT_TOL);
  ctx->isreal             = PETSC_FALSE;
  ctx->npart              = 1;

  nep->useds = PETSC_TRUE;

  nep->ops->solve          = NEPSolve_CISS;
  nep->ops->setup          = NEPSetUp_CISS;
  nep->ops->setfromoptions = NEPSetFromOptions_CISS;
  nep->ops->reset          = NEPReset_CISS;
  nep->ops->destroy        = NEPDestroy_CISS;
  nep->ops->view           = NEPView_CISS;
  nep->ops->setdstype      = NEPSetDSType_CISS;

  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSSetSizes_C",NEPCISSSetSizes_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetSizes_C",NEPCISSGetSizes_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSSetThreshold_C",NEPCISSSetThreshold_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetThreshold_C",NEPCISSGetThreshold_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSSetRefinement_C",NEPCISSSetRefinement_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetRefinement_C",NEPCISSGetRefinement_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSSetExtraction_C",NEPCISSSetExtraction_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetExtraction_C",NEPCISSGetExtraction_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPCISSGetKSPs_C",NEPCISSGetKSPs_CISS));
  PetscFunctionReturn(PETSC_SUCCESS);
}
