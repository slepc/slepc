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
           numerical method for polynomial eigenvalue problems using contour
           integral", Japan J. Indust. Appl. Math. 27:73-90, 2010.
*/

#include <slepc/private/pepimpl.h>         /*I "slepcpep.h" I*/
#include <slepc/private/slepccontour.h>

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
  PEPCISSExtraction extraction;
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
  Mat               J,*Psplit;     /* auxiliary matrices */
  BV                pV;
  PetscObjectId     rgid;
  PetscObjectState  rgstate;
} PEP_CISS;

static PetscErrorCode PEPComputeFunction(PEP pep,PetscScalar lambda,Mat T,Mat P,PetscBool deriv)
{
  PetscInt         i;
  PetscScalar      *coeff;
  Mat              *A,*K;
  MatStructure     str,strp;
  PEP_CISS         *ctx = (PEP_CISS*)pep->data;
  SlepcContourData contour = ctx->contour;

  PetscFunctionBegin;
  A = (contour->pA)?contour->pA:pep->A;
  K = (contour->pP)?contour->pP:ctx->Psplit;
  PetscCall(PetscMalloc1(pep->nmat,&coeff));
  if (deriv) PetscCall(PEPEvaluateBasisDerivative(pep,lambda,0,coeff,NULL));
  else PetscCall(PEPEvaluateBasis(pep,lambda,0,coeff,NULL));
  PetscCall(STGetMatStructure(pep->st,&str));
  PetscCall(MatZeroEntries(T));
  if (!deriv && T != P) {
    PetscCall(STGetSplitPreconditionerInfo(pep->st,NULL,&strp));
    PetscCall(MatZeroEntries(P));
  }
  i = deriv?1:0;
  for (;i<pep->nmat;i++) {
    PetscCall(MatAXPY(T,coeff[i],A[i],str));
    if (!deriv && T != P) PetscCall(MatAXPY(P,coeff[i],K[i],strp));
  }
  PetscCall(PetscFree(coeff));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Set up KSP solvers for every integration point
*/
static PetscErrorCode PEPCISSSetUp(PEP pep,Mat T,Mat P)
{
  PEP_CISS         *ctx = (PEP_CISS*)pep->data;
  SlepcContourData contour;
  PetscInt         i,p_id;
  Mat              Amat,Pmat;

  PetscFunctionBegin;
  if (!ctx->contour || !ctx->contour->ksp) PetscCall(PEPCISSGetKSPs(pep,NULL,NULL));
  contour = ctx->contour;
  PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Something went wrong with PEPCISSGetKSPs()");
  for (i=0;i<contour->npoints;i++) {
    p_id = i*contour->subcomm->n + contour->subcomm->color;
    PetscCall(MatDuplicate(T,MAT_DO_NOT_COPY_VALUES,&Amat));
    if (T != P) PetscCall(MatDuplicate(P,MAT_DO_NOT_COPY_VALUES,&Pmat)); else Pmat = Amat;
    PetscCall(PEPComputeFunction(pep,ctx->omega[p_id],Amat,Pmat,PETSC_FALSE));
    PetscCall(PEP_KSPSetOperators(contour->ksp[i],Amat,Pmat));
    PetscCall(MatDestroy(&Amat));
    if (T != P) PetscCall(MatDestroy(&Pmat));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Y_i = F(z_i)^{-1}Fp(z_i)V for every integration point, Y=[Y_i] is in the context
*/
static PetscErrorCode PEPCISSSolve(PEP pep,Mat dT,BV V,PetscInt L_start,PetscInt L_end)
{
  PEP_CISS         *ctx = (PEP_CISS*)pep->data;
  SlepcContourData contour;
  PetscInt         i,p_id;
  Mat              MV,BMV=NULL,MC;

  PetscFunctionBegin;
  contour = ctx->contour;
  PetscCall(BVSetActiveColumns(V,L_start,L_end));
  PetscCall(BVGetMat(V,&MV));
  for (i=0;i<contour->npoints;i++) {
    p_id = i*contour->subcomm->n + contour->subcomm->color;
    PetscCall(PEPComputeFunction(pep,ctx->omega[p_id],dT,NULL,PETSC_TRUE));
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

static PetscErrorCode PEPSetUp_CISS(PEP pep)
{
  PEP_CISS         *ctx = (PEP_CISS*)pep->data;
  SlepcContourData contour;
  PetscInt         i,nwork,nsplit;
  PetscBool        istrivial,isellipse,flg;
  PetscObjectId    id;
  PetscObjectState state;
  Vec              v0;

  PetscFunctionBegin;
  if (pep->ncv==PETSC_DETERMINE) pep->ncv = ctx->L_max*ctx->M;
  else {
    ctx->L_max = pep->ncv/ctx->M;
    if (!ctx->L_max) {
      ctx->L_max = 1;
      pep->ncv = ctx->L_max*ctx->M;
    }
  }
  ctx->L = PetscMin(ctx->L,ctx->L_max);
  if (pep->max_it==PETSC_DETERMINE) pep->max_it = 5;
  if (pep->mpd==PETSC_DETERMINE) pep->mpd = pep->ncv;
  if (!pep->which) pep->which = PEP_ALL;
  PetscCheck(pep->which==PEP_ALL,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"This solver supports only computing all eigenvalues");
  PEPCheckUnsupported(pep,PEP_FEATURE_STOPPING);
  PEPCheckIgnored(pep,PEP_FEATURE_SCALE);

  /* check region */
  PetscCall(RGIsTrivial(pep->rg,&istrivial));
  PetscCheck(!istrivial,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"CISS requires a nontrivial region, e.g. -rg_type ellipse ...");
  PetscCall(RGGetComplement(pep->rg,&flg));
  PetscCheck(!flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"A region with complement flag set is not allowed");
  PetscCall(PetscObjectTypeCompare((PetscObject)pep->rg,RGELLIPSE,&isellipse));
  PetscCheck(isellipse,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Currently only implemented for elliptic regions");

  /* if the region has changed, then reset contour data */
  PetscCall(PetscObjectGetId((PetscObject)pep->rg,&id));
  PetscCall(PetscObjectStateGet((PetscObject)pep->rg,&state));
  if (ctx->rgid && (id != ctx->rgid || state != ctx->rgstate)) {
    PetscCall(SlepcContourDataDestroy(&ctx->contour));
    PetscCall(PetscInfo(pep,"Resetting the contour data structure due to a change of region\n"));
    ctx->rgid = id; ctx->rgstate = state;
  }

  /* create contour data structure */
  if (!ctx->contour) {
    PetscCall(RGCanUseConjugates(pep->rg,ctx->isreal,&ctx->useconj));
    PetscCall(SlepcContourDataCreate(ctx->useconj?ctx->N/2:ctx->N,ctx->npart,(PetscObject)pep,&ctx->contour));
  }

  PetscCall(PEPAllocateSolution(pep,0));
  if (ctx->weight) PetscCall(PetscFree4(ctx->weight,ctx->omega,ctx->pp,ctx->sigma));
  PetscCall(PetscMalloc4(ctx->N,&ctx->weight,ctx->N,&ctx->omega,ctx->N,&ctx->pp,ctx->L_max*ctx->M,&ctx->sigma));

  /* allocate basis vectors */
  PetscCall(BVDestroy(&ctx->S));
  PetscCall(BVDuplicateResize(pep->V,ctx->L*ctx->M,&ctx->S));
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVDuplicateResize(pep->V,ctx->L,&ctx->V));

  /* check if a user-defined split preconditioner has been set */
  PetscCall(STGetSplitPreconditionerInfo(pep->st,&nsplit,NULL));
  if (nsplit) {
    PetscCall(PetscFree(ctx->Psplit));
    PetscCall(PetscMalloc1(nsplit,&ctx->Psplit));
    for (i=0;i<nsplit;i++) PetscCall(STGetSplitPreconditionerTerm(pep->st,i,&ctx->Psplit[i]));
  }

  contour = ctx->contour;
  PetscCall(SlepcContourRedundantMat(contour,pep->nmat,pep->A,ctx->Psplit));
  if (!ctx->J) PetscCall(MatDuplicate(contour->pA?contour->pA[0]:pep->A[0],MAT_DO_NOT_COPY_VALUES,&ctx->J));
  if (contour->pA) {
    PetscCall(BVGetColumn(ctx->V,0,&v0));
    PetscCall(SlepcContourScatterCreate(contour,v0));
    PetscCall(BVRestoreColumn(ctx->V,0,&v0));
    PetscCall(BVDestroy(&ctx->pV));
    PetscCall(BVCreate(PetscObjectComm((PetscObject)contour->xsub),&ctx->pV));
    PetscCall(BVSetSizesFromVec(ctx->pV,contour->xsub,pep->n));
    PetscCall(BVSetFromOptions(ctx->pV));
    PetscCall(BVResize(ctx->pV,ctx->L,PETSC_FALSE));
  }

  PetscCall(BVDestroy(&ctx->Y));
  if (contour->pA) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)contour->xsub),&ctx->Y));
    PetscCall(BVSetSizesFromVec(ctx->Y,contour->xsub,pep->n));
    PetscCall(BVSetFromOptions(ctx->Y));
    PetscCall(BVResize(ctx->Y,contour->npoints*ctx->L,PETSC_FALSE));
  } else PetscCall(BVDuplicateResize(pep->V,contour->npoints*ctx->L,&ctx->Y));

  if (ctx->extraction == PEP_CISS_EXTRACTION_RITZ) {
    PetscCall(DSPEPSetDegree(pep->ds,pep->nmat-1));
    PetscCall(DSPEPSetCoefficients(pep->ds,pep->pbc));
  }
  PetscCall(DSAllocate(pep->ds,pep->ncv));
  nwork = 2;
  PetscCall(PEPSetWorkVecs(pep,nwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSolve_CISS(PEP pep)
{
  PEP_CISS         *ctx = (PEP_CISS*)pep->data;
  SlepcContourData contour = ctx->contour;
  Mat              X,M,E,T,P;
  PetscInt         i,j,ld,L_add=0,nv=0,L_base=ctx->L,inner,*inside,nsplit;
  PetscScalar      *Mu,*H0,*H1,*rr,*temp,center;
  PetscReal        error,max_error,radius,rgscale,est_eig,eta;
  PetscBool        isellipse,*fl1;
  Vec              si;
  SlepcSC          sc;
  PetscRandom      rand;

  PetscFunctionBegin;
  PetscCall(DSGetSlepcSC(pep->ds,&sc));
  sc->comparison    = SlepcCompareLargestMagnitude;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSGetLeadingDimension(pep->ds,&ld));
  PetscCall(RGComputeQuadrature(pep->rg,RG_QUADRULE_TRAPEZOIDAL,ctx->N,ctx->omega,ctx->pp,ctx->weight));
  PetscCall(STGetSplitPreconditionerInfo(pep->st,&nsplit,NULL));
  if (contour->pA) {
    T = contour->pA[0];
    P = nsplit? contour->pP[0]: T;
  } else {
    T = pep->A[0];
    P = nsplit? ctx->Psplit[0]: T;
  }
  PetscCall(PEPCISSSetUp(pep,T,P));
  PetscCall(BVSetActiveColumns(ctx->V,0,ctx->L));
  PetscCall(BVSetRandomSign(ctx->V));
  PetscCall(BVGetRandomContext(ctx->V,&rand));
  if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
  PetscCall(PEPCISSSolve(pep,ctx->J,(contour->pA)?ctx->pV:ctx->V,0,ctx->L));
  PetscCall(PetscObjectTypeCompare((PetscObject)pep->rg,RGELLIPSE,&isellipse));
  if (isellipse) {
    PetscCall(BVTraceQuadrature(ctx->Y,ctx->V,ctx->L,ctx->L,ctx->weight,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj,&est_eig));
    PetscCall(PetscInfo(pep,"Estimated eigenvalue count: %f\n",(double)est_eig));
    eta = PetscPowReal(10.0,-PetscLog10Real(pep->tol)/ctx->N);
    L_add = PetscMax(0,(PetscInt)PetscCeilReal((est_eig*eta)/ctx->M)-ctx->L);
    if (L_add>ctx->L_max-ctx->L) {
      PetscCall(PetscInfo(pep,"Number of eigenvalues inside the contour path may be too large\n"));
      L_add = ctx->L_max-ctx->L;
    }
  }
  /* Updates L after estimate the number of eigenvalue */
  if (L_add>0) {
    PetscCall(PetscInfo(pep,"Changing L %" PetscInt_FMT " -> %" PetscInt_FMT " by Estimate #Eig\n",ctx->L,ctx->L+L_add));
    PetscCall(BVCISSResizeBases(ctx->S,contour->pA?ctx->pV:ctx->V,ctx->Y,ctx->L,ctx->L+L_add,ctx->M,contour->npoints));
    PetscCall(BVSetActiveColumns(ctx->V,ctx->L,ctx->L+L_add));
    PetscCall(BVSetRandomSign(ctx->V));
    if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
    ctx->L += L_add;
    PetscCall(PEPCISSSolve(pep,ctx->J,(contour->pA)?ctx->pV:ctx->V,ctx->L-L_add,ctx->L));
  }

  PetscCall(PetscMalloc2(ctx->L*ctx->L*ctx->M*2,&Mu,ctx->L*ctx->M*ctx->L*ctx->M,&H0));
  for (i=0;i<ctx->refine_blocksize;i++) {
    PetscCall(BVDotQuadrature(ctx->Y,(contour->pA)?ctx->pV:ctx->V,Mu,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->subcomm,contour->npoints,ctx->useconj));
    PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
    PetscCall(PetscLogEventBegin(PEP_CISS_SVD,pep,0,0,0));
    PetscCall(SlepcCISS_BH_SVD(H0,ctx->L*ctx->M,ctx->delta,ctx->sigma,&nv));
    PetscCall(PetscLogEventEnd(PEP_CISS_SVD,pep,0,0,0));
    if (ctx->sigma[0]<=ctx->delta || nv < ctx->L*ctx->M || ctx->L == ctx->L_max) break;
    L_add = L_base;
    if (ctx->L+L_add>ctx->L_max) L_add = ctx->L_max-ctx->L;
    PetscCall(PetscInfo(pep,"Changing L %" PetscInt_FMT " -> %" PetscInt_FMT " by SVD(H0)\n",ctx->L,ctx->L+L_add));
    PetscCall(BVCISSResizeBases(ctx->S,contour->pA?ctx->pV:ctx->V,ctx->Y,ctx->L,ctx->L+L_add,ctx->M,contour->npoints));
    PetscCall(BVSetActiveColumns(ctx->V,ctx->L,ctx->L+L_add));
    PetscCall(BVSetRandomSign(ctx->V));
    if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
    ctx->L += L_add;
    PetscCall(PEPCISSSolve(pep,ctx->J,(contour->pA)?ctx->pV:ctx->V,ctx->L-L_add,ctx->L));
    if (L_add) {
      PetscCall(PetscFree2(Mu,H0));
      PetscCall(PetscMalloc2(ctx->L*ctx->L*ctx->M*2,&Mu,ctx->L*ctx->M*ctx->L*ctx->M,&H0));
    }
  }

  PetscCall(RGGetScale(pep->rg,&rgscale));
  PetscCall(RGEllipseGetParameters(pep->rg,&center,&radius,NULL));

  if (ctx->extraction == PEP_CISS_EXTRACTION_HANKEL) PetscCall(PetscMalloc1(ctx->L*ctx->M*ctx->L*ctx->M,&H1));

  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;
    for (inner=0;inner<=ctx->refine_inner;inner++) {
      if (ctx->extraction == PEP_CISS_EXTRACTION_HANKEL) {
        PetscCall(BVDotQuadrature(ctx->Y,(contour->pA)?ctx->pV:ctx->V,Mu,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
        PetscCall(PetscLogEventBegin(PEP_CISS_SVD,pep,0,0,0));
        PetscCall(SlepcCISS_BH_SVD(H0,ctx->L*ctx->M,ctx->delta,ctx->sigma,&nv));
        PetscCall(PetscLogEventEnd(PEP_CISS_SVD,pep,0,0,0));
      } else {
        PetscCall(BVSumQuadrature(ctx->S,ctx->Y,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj));
        /* compute SVD of S */
        PetscCall(BVSVDAndRank(ctx->S,ctx->M,ctx->L,ctx->delta,(ctx->extraction==PEP_CISS_EXTRACTION_CAA)?BV_SVD_METHOD_QR_CAA:BV_SVD_METHOD_QR,H0,ctx->sigma,&nv));
      }
      PetscCall(PetscInfo(pep,"Estimated rank: nv = %" PetscInt_FMT "\n",nv));
      if (ctx->sigma[0]>ctx->delta && nv==ctx->L*ctx->M && inner!=ctx->refine_inner) {
        PetscCall(BVSumQuadrature(ctx->S,ctx->Y,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(BVSetActiveColumns(ctx->S,0,ctx->L));
        PetscCall(BVSetActiveColumns(ctx->V,0,ctx->L));
        PetscCall(BVCopy(ctx->S,ctx->V));
        if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
        PetscCall(PEPCISSSolve(pep,ctx->J,(contour->pA)?ctx->pV:ctx->V,0,ctx->L));
      } else break;
    }
    pep->nconv = 0;
    if (nv == 0) { pep->reason = PEP_CONVERGED_TOL; break; }
    else {
      /* Extracting eigenpairs */
      PetscCall(DSSetDimensions(pep->ds,nv,0,0));
      PetscCall(DSSetState(pep->ds,DS_STATE_RAW));
      if (ctx->extraction == PEP_CISS_EXTRACTION_HANKEL) {
        PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
        PetscCall(CISS_BlockHankel(Mu,1,ctx->L,ctx->M,H1));
        PetscCall(DSGetArray(pep->ds,DS_MAT_A,&temp));
        for (j=0;j<nv;j++)
          for (i=0;i<nv;i++)
            temp[i+j*ld] = H1[i+j*ctx->L*ctx->M];
        PetscCall(DSRestoreArray(pep->ds,DS_MAT_A,&temp));
        PetscCall(DSGetArray(pep->ds,DS_MAT_B,&temp));
        for (j=0;j<nv;j++)
          for (i=0;i<nv;i++)
            temp[i+j*ld] = H0[i+j*ctx->L*ctx->M];
        PetscCall(DSRestoreArray(pep->ds,DS_MAT_B,&temp));
      } else if (ctx->extraction == PEP_CISS_EXTRACTION_CAA) {
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
        PetscCall(DSGetArray(pep->ds,DS_MAT_A,&temp));
        for (i=0;i<nv;i++) PetscCall(PetscArraycpy(temp+i*ld,H0+i*nv,nv));
        PetscCall(DSRestoreArray(pep->ds,DS_MAT_A,&temp));
      } else {
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
        for (i=0;i<pep->nmat;i++) {
          PetscCall(DSGetMat(pep->ds,DSMatExtra[i],&E));
          PetscCall(BVMatProject(ctx->S,pep->A[i],ctx->S,E));
          PetscCall(DSRestoreMat(pep->ds,DSMatExtra[i],&E));
        }
        nv = (pep->nmat-1)*nv;
      }
      PetscCall(DSSolve(pep->ds,pep->eigr,pep->eigi));
      PetscCall(DSSynchronize(pep->ds,pep->eigr,pep->eigi));
      if (ctx->extraction == PEP_CISS_EXTRACTION_CAA || ctx->extraction == PEP_CISS_EXTRACTION_HANKEL) {
        for (i=0;i<nv;i++) {
          pep->eigr[i] = (pep->eigr[i]*radius+center)*rgscale;
        }
      }
      PetscCall(PetscMalloc3(nv,&fl1,nv,&inside,nv,&rr));
      PetscCall(DSVectors(pep->ds,DS_MAT_X,NULL,NULL));
      PetscCall(DSGetMat(pep->ds,DS_MAT_X,&X));
      PetscCall(SlepcCISS_isGhost(X,nv,ctx->sigma,ctx->spurious_threshold,fl1));
      PetscCall(DSRestoreMat(pep->ds,DS_MAT_X,&X));
      PetscCall(RGCheckInside(pep->rg,nv,pep->eigr,pep->eigi,inside));
      for (i=0;i<nv;i++) {
        if (fl1[i] && inside[i]>=0) {
          rr[i] = 1.0;
          pep->nconv++;
        } else rr[i] = 0.0;
      }
      PetscCall(DSSort(pep->ds,pep->eigr,pep->eigi,rr,NULL,&pep->nconv));
      PetscCall(DSSynchronize(pep->ds,pep->eigr,pep->eigi));
      if (ctx->extraction == PEP_CISS_EXTRACTION_CAA || ctx->extraction == PEP_CISS_EXTRACTION_HANKEL) {
        for (i=0;i<nv;i++) pep->eigr[i] = (pep->eigr[i]*radius+center)*rgscale;
      }
      PetscCall(PetscFree3(fl1,inside,rr));
      PetscCall(BVSetActiveColumns(pep->V,0,nv));
      PetscCall(DSVectors(pep->ds,DS_MAT_X,NULL,NULL));
      if (ctx->extraction == PEP_CISS_EXTRACTION_HANKEL) {
        PetscCall(BVSumQuadrature(ctx->S,ctx->Y,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
        PetscCall(BVCopy(ctx->S,pep->V));
        PetscCall(DSGetMat(pep->ds,DS_MAT_X,&X));
        PetscCall(BVMultInPlace(ctx->S,X,0,pep->nconv));
        PetscCall(BVMultInPlace(pep->V,X,0,pep->nconv));
        PetscCall(DSRestoreMat(pep->ds,DS_MAT_X,&X));
      } else {
        PetscCall(DSGetMat(pep->ds,DS_MAT_X,&X));
        PetscCall(BVMultInPlace(ctx->S,X,0,pep->nconv));
        PetscCall(DSRestoreMat(pep->ds,DS_MAT_X,&X));
        PetscCall(BVSetActiveColumns(ctx->S,0,pep->nconv));
        PetscCall(BVCopy(ctx->S,pep->V));
      }
      max_error = 0.0;
      for (i=0;i<pep->nconv;i++) {
        PetscCall(BVGetColumn(pep->V,i,&si));
        PetscCall(VecNormalize(si,NULL));
        PetscCall(PEPComputeResidualNorm_Private(pep,pep->eigr[i],0,si,NULL,pep->work,&error));
        PetscCall((*pep->converged)(pep,pep->eigr[i],0,error,&error,pep->convergedctx));
        PetscCall(BVRestoreColumn(pep->V,i,&si));
        max_error = PetscMax(max_error,error);
      }
      if (max_error <= pep->tol) pep->reason = PEP_CONVERGED_TOL;
      else if (pep->its > pep->max_it) pep->reason = PEP_DIVERGED_ITS;
      else {
        if (pep->nconv > ctx->L) nv = pep->nconv;
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
        PetscCall(PEPCISSSolve(pep,ctx->J,(contour->pA)?ctx->pV:ctx->V,0,ctx->L));
      }
    }
  }
  PetscCall(PetscFree2(Mu,H0));
  if (ctx->extraction == PEP_CISS_EXTRACTION_HANKEL) PetscCall(PetscFree(H1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSSetSizes_CISS(PEP pep,PetscInt ip,PetscInt bs,PetscInt ms,PetscInt npart,PetscInt bsmax,PetscBool realmats)
{
  PEP_CISS       *ctx = (PEP_CISS*)pep->data;
  PetscInt       oN,oL,oM,oLmax,onpart;
  PetscMPIInt    size;

  PetscFunctionBegin;
  oN = ctx->N;
  if (ip == PETSC_DETERMINE) {
    if (ctx->N!=32) { ctx->N =32; ctx->M = ctx->N/4; }
  } else if (ip != PETSC_CURRENT) {
    PetscCheck(ip>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The ip argument must be > 0");
    PetscCheck(ip%2==0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The ip argument must be an even number");
    if (ctx->N!=ip) { ctx->N = ip; ctx->M = ctx->N/4; }
  }
  oL = ctx->L;
  if (bs == PETSC_DETERMINE) {
    ctx->L = 16;
  } else if (bs != PETSC_CURRENT) {
    PetscCheck(bs>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The bs argument must be > 0");
    ctx->L = bs;
  }
  oM = ctx->M;
  if (ms == PETSC_DETERMINE) {
    ctx->M = ctx->N/4;
  } else if (ms != PETSC_CURRENT) {
    PetscCheck(ms>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The ms argument must be > 0");
    PetscCheck(ms<=ctx->N,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The ms argument must be less than or equal to the number of integration points");
    ctx->M = ms;
  }
  onpart = ctx->npart;
  if (npart == PETSC_DETERMINE) {
    ctx->npart = 1;
  } else if (npart != PETSC_CURRENT) {
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pep),&size));
    PetscCheck(npart>0 && npart<=size,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of npart");
    ctx->npart = npart;
  }
  oLmax = ctx->L_max;
  if (bsmax == PETSC_DETERMINE) {
    ctx->L_max = 64;
  } else if (bsmax != PETSC_CURRENT) {
    PetscCheck(bsmax>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The bsmax argument must be > 0");
    ctx->L_max = PetscMax(bsmax,ctx->L);
  }
  if (onpart != ctx->npart || oN != ctx->N || realmats != ctx->isreal) {
    PetscCall(SlepcContourDataDestroy(&ctx->contour));
    PetscCall(PetscInfo(pep,"Resetting the contour data structure due to a change of parameters\n"));
    pep->state = PEP_STATE_INITIAL;
  }
  ctx->isreal = realmats;
  if (oL != ctx->L || oM != ctx->M || oLmax != ctx->L_max) pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPCISSSetSizes - Sets the values of various size parameters in the CISS solver.

   Logically Collective

   Input Parameters:
+  pep   - the polynomial eigensolver context
.  ip    - number of integration points
.  bs    - block size
.  ms    - moment size
.  npart - number of partitions when splitting the communicator
.  bsmax - max block size
-  realmats - all coefficient matrices of P(.) are real

   Options Database Keys:
+  -pep_ciss_integration_points - Sets the number of integration points
.  -pep_ciss_blocksize - Sets the block size
.  -pep_ciss_moments - Sets the moment size
.  -pep_ciss_partitions - Sets the number of partitions
.  -pep_ciss_maxblocksize - Sets the maximum block size
-  -pep_ciss_realmats - all coefficient matrices of P(.) are real

   Notes:
   For all integer arguments, you can use PETSC_CURRENT to keep the current value, and
   PETSC_DETERMINE to set them to a default value.

   The default number of partitions is 1. This means the internal KSP object is shared
   among all processes of the PEP communicator. Otherwise, the communicator is split
   into npart communicators, so that npart KSP solves proceed simultaneously.

   Level: advanced

.seealso: PEPCISSGetSizes()
@*/
PetscErrorCode PEPCISSSetSizes(PEP pep,PetscInt ip,PetscInt bs,PetscInt ms,PetscInt npart,PetscInt bsmax,PetscBool realmats)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,ip,2);
  PetscValidLogicalCollectiveInt(pep,bs,3);
  PetscValidLogicalCollectiveInt(pep,ms,4);
  PetscValidLogicalCollectiveInt(pep,npart,5);
  PetscValidLogicalCollectiveInt(pep,bsmax,6);
  PetscValidLogicalCollectiveBool(pep,realmats,7);
  PetscTryMethod(pep,"PEPCISSSetSizes_C",(PEP,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool),(pep,ip,bs,ms,npart,bsmax,realmats));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSGetSizes_CISS(PEP pep,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats)
{
  PEP_CISS *ctx = (PEP_CISS*)pep->data;

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
   PEPCISSGetSizes - Gets the values of various size parameters in the CISS solver.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  ip    - number of integration points
.  bs    - block size
.  ms    - moment size
.  npart - number of partitions when splitting the communicator
.  bsmax - max block size
-  realmats - all coefficient matrices of P(.) are real

   Level: advanced

.seealso: PEPCISSSetSizes()
@*/
PetscErrorCode PEPCISSGetSizes(PEP pep,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscUseMethod(pep,"PEPCISSGetSizes_C",(PEP,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*),(pep,ip,bs,ms,npart,bsmax,realmats));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSSetThreshold_CISS(PEP pep,PetscReal delta,PetscReal spur)
{
  PEP_CISS *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  if (delta == (PetscReal)PETSC_DETERMINE) {
    ctx->delta = SLEPC_DEFAULT_TOL*1e-4;
  } else if (delta != (PetscReal)PETSC_CURRENT) {
    PetscCheck(delta>0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The delta argument must be > 0.0");
    ctx->delta = delta;
  }
  if (spur == (PetscReal)PETSC_DETERMINE) {
    ctx->spurious_threshold = PetscSqrtReal(SLEPC_DEFAULT_TOL);
  } else if (spur != (PetscReal)PETSC_CURRENT) {
    PetscCheck(spur>0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The spurious threshold argument must be > 0.0");
    ctx->spurious_threshold = spur;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPCISSSetThreshold - Sets the values of various threshold parameters in
   the CISS solver.

   Logically Collective

   Input Parameters:
+  pep   - the polynomial eigensolver context
.  delta - threshold for numerical rank
-  spur  - spurious threshold (to discard spurious eigenpairs)

   Options Database Keys:
+  -pep_ciss_delta - Sets the delta
-  -pep_ciss_spurious_threshold - Sets the spurious threshold

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: PEPCISSGetThreshold()
@*/
PetscErrorCode PEPCISSSetThreshold(PEP pep,PetscReal delta,PetscReal spur)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,delta,2);
  PetscValidLogicalCollectiveReal(pep,spur,3);
  PetscTryMethod(pep,"PEPCISSSetThreshold_C",(PEP,PetscReal,PetscReal),(pep,delta,spur));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSGetThreshold_CISS(PEP pep,PetscReal *delta,PetscReal *spur)
{
  PEP_CISS *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  if (delta) *delta = ctx->delta;
  if (spur)  *spur = ctx->spurious_threshold;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPCISSGetThreshold - Gets the values of various threshold parameters in
   the CISS solver.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  delta - threshold for numerical rank
-  spur  - spurious threshold (to discard spurious eigenpairs)

   Level: advanced

.seealso: PEPCISSSetThreshold()
@*/
PetscErrorCode PEPCISSGetThreshold(PEP pep,PetscReal *delta,PetscReal *spur)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscUseMethod(pep,"PEPCISSGetThreshold_C",(PEP,PetscReal*,PetscReal*),(pep,delta,spur));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSSetRefinement_CISS(PEP pep,PetscInt inner,PetscInt blsize)
{
  PEP_CISS *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  if (inner == PETSC_DETERMINE) {
    ctx->refine_inner = 0;
  } else if (inner != PETSC_CURRENT) {
    PetscCheck(inner>=0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The refine inner argument must be >= 0");
    ctx->refine_inner = inner;
  }
  if (blsize == PETSC_DETERMINE) {
    ctx->refine_blocksize = 0;
  } else if (blsize != PETSC_CURRENT) {
    PetscCheck(blsize>=0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The refine blocksize argument must be >= 0");
    ctx->refine_blocksize = blsize;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPCISSSetRefinement - Sets the values of various refinement parameters
   in the CISS solver.

   Logically Collective

   Input Parameters:
+  pep    - the polynomial eigensolver context
.  inner  - number of iterative refinement iterations (inner loop)
-  blsize - number of iterative refinement iterations (blocksize loop)

   Options Database Keys:
+  -pep_ciss_refine_inner - Sets number of inner iterations
-  -pep_ciss_refine_blocksize - Sets number of blocksize iterations

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: PEPCISSGetRefinement()
@*/
PetscErrorCode PEPCISSSetRefinement(PEP pep,PetscInt inner,PetscInt blsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,inner,2);
  PetscValidLogicalCollectiveInt(pep,blsize,3);
  PetscTryMethod(pep,"PEPCISSSetRefinement_C",(PEP,PetscInt,PetscInt),(pep,inner,blsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSGetRefinement_CISS(PEP pep,PetscInt *inner,PetscInt *blsize)
{
  PEP_CISS *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  if (inner)  *inner = ctx->refine_inner;
  if (blsize) *blsize = ctx->refine_blocksize;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPCISSGetRefinement - Gets the values of various refinement parameters
   in the CISS solver.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  inner  - number of iterative refinement iterations (inner loop)
-  blsize - number of iterative refinement iterations (blocksize loop)

   Level: advanced

.seealso: PEPCISSSetRefinement()
@*/
PetscErrorCode PEPCISSGetRefinement(PEP pep, PetscInt *inner, PetscInt *blsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscUseMethod(pep,"PEPCISSGetRefinement_C",(PEP,PetscInt*,PetscInt*),(pep,inner,blsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSSetExtraction_CISS(PEP pep,PEPCISSExtraction extraction)
{
  PEP_CISS *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  if (ctx->extraction != extraction) {
    ctx->extraction = extraction;
    pep->state      = PEP_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPCISSSetExtraction - Sets the extraction technique used in the CISS solver.

   Logically Collective

   Input Parameters:
+  pep        - the polynomial eigensolver context
-  extraction - the extraction technique

   Options Database Key:
.  -pep_ciss_extraction - Sets the extraction technique (either 'ritz', 'hankel' or 'caa')

   Notes:
   By default, the Rayleigh-Ritz extraction is used (PEP_CISS_EXTRACTION_RITZ).

   If the 'hankel' or the 'caa' option is specified (PEP_CISS_EXTRACTION_HANKEL or
   PEP_CISS_EXTRACTION_CAA), then the Block Hankel method, or the Communication-avoiding
   Arnoldi method, respectively, is used for extracting eigenpairs.

   Level: advanced

.seealso: PEPCISSGetExtraction(), PEPCISSExtraction
@*/
PetscErrorCode PEPCISSSetExtraction(PEP pep,PEPCISSExtraction extraction)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,extraction,2);
  PetscTryMethod(pep,"PEPCISSSetExtraction_C",(PEP,PEPCISSExtraction),(pep,extraction));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSGetExtraction_CISS(PEP pep,PEPCISSExtraction *extraction)
{
  PEP_CISS *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  *extraction = ctx->extraction;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPCISSGetExtraction - Gets the extraction technique used in the CISS solver.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
.  extraction - extraction technique

   Level: advanced

.seealso: PEPCISSSetExtraction() PEPCISSExtraction
@*/
PetscErrorCode PEPCISSGetExtraction(PEP pep,PEPCISSExtraction *extraction)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(extraction,2);
  PetscUseMethod(pep,"PEPCISSGetExtraction_C",(PEP,PEPCISSExtraction*),(pep,extraction));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPCISSGetKSPs_CISS(PEP pep,PetscInt *nsolve,KSP **ksp)
{
  PEP_CISS         *ctx = (PEP_CISS*)pep->data;
  SlepcContourData contour;
  PetscInt         i,nsplit;
  PC               pc;
  MPI_Comm         child;

  PetscFunctionBegin;
  if (!ctx->contour) {  /* initialize contour data structure first */
    PetscCall(RGCanUseConjugates(pep->rg,ctx->isreal,&ctx->useconj));
    PetscCall(SlepcContourDataCreate(ctx->useconj?ctx->N/2:ctx->N,ctx->npart,(PetscObject)pep,&ctx->contour));
  }
  contour = ctx->contour;
  if (!contour->ksp) {
    PetscCall(PetscMalloc1(contour->npoints,&contour->ksp));
    PetscCall(PEPGetST(pep,&pep->st));
    PetscCall(STGetSplitPreconditionerInfo(pep->st,&nsplit,NULL));
    PetscCall(PetscSubcommGetChild(contour->subcomm,&child));
    for (i=0;i<contour->npoints;i++) {
      PetscCall(KSPCreate(child,&contour->ksp[i]));
      PetscCall(PetscObjectIncrementTabLevel((PetscObject)contour->ksp[i],(PetscObject)pep,1));
      PetscCall(KSPSetOptionsPrefix(contour->ksp[i],((PetscObject)pep)->prefix));
      PetscCall(KSPAppendOptionsPrefix(contour->ksp[i],"pep_ciss_"));
      PetscCall(PetscObjectSetOptions((PetscObject)contour->ksp[i],((PetscObject)pep)->options));
      PetscCall(KSPSetErrorIfNotConverged(contour->ksp[i],PETSC_TRUE));
      PetscCall(KSPSetTolerances(contour->ksp[i],SlepcDefaultTol(pep->tol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
      PetscCall(KSPGetPC(contour->ksp[i],&pc));
      if (nsplit) {
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
   PEPCISSGetKSPs - Retrieve the array of linear solver objects associated with
   the CISS solver.

   Collective

   Input Parameter:
.  pep - polynomial eigenvalue solver

   Output Parameters:
+  nsolve - number of solver objects
-  ksp - array of linear solver object

   Notes:
   The number of KSP solvers is equal to the number of integration points divided by
   the number of partitions. This value is halved in the case of real matrices with
   a region centered at the real axis.

   Level: advanced

.seealso: PEPCISSSetSizes()
@*/
PetscErrorCode PEPCISSGetKSPs(PEP pep,PetscInt *nsolve,KSP **ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscUseMethod(pep,"PEPCISSGetKSPs_C",(PEP,PetscInt*,KSP**),(pep,nsolve,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPReset_CISS(PEP pep)
{
  PEP_CISS       *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&ctx->S));
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVDestroy(&ctx->Y));
  PetscCall(SlepcContourDataReset(ctx->contour));
  PetscCall(MatDestroy(&ctx->J));
  PetscCall(BVDestroy(&ctx->pV));
  PetscCall(PetscFree(ctx->Psplit));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetFromOptions_CISS(PEP pep,PetscOptionItems PetscOptionsObject)
{
  PEP_CISS          *ctx = (PEP_CISS*)pep->data;
  PetscReal         r1,r2;
  PetscInt          i,i1,i2,i3,i4,i5,i6,i7;
  PetscBool         b1,flg,flg2,flg3,flg4,flg5,flg6;
  PEPCISSExtraction extraction;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"PEP CISS Options");

    PetscCall(PEPCISSGetSizes(pep,&i1,&i2,&i3,&i4,&i5,&b1));
    PetscCall(PetscOptionsInt("-pep_ciss_integration_points","Number of integration points","PEPCISSSetSizes",i1,&i1,&flg));
    PetscCall(PetscOptionsInt("-pep_ciss_blocksize","Block size","PEPCISSSetSizes",i2,&i2,&flg2));
    PetscCall(PetscOptionsInt("-pep_ciss_moments","Moment size","PEPCISSSetSizes",i3,&i3,&flg3));
    PetscCall(PetscOptionsInt("-pep_ciss_partitions","Number of partitions","PEPCISSSetSizes",i4,&i4,&flg4));
    PetscCall(PetscOptionsInt("-pep_ciss_maxblocksize","Maximum block size","PEPCISSSetSizes",i5,&i5,&flg5));
    PetscCall(PetscOptionsBool("-pep_ciss_realmats","True if all coefficient matrices of P(.) are real","PEPCISSSetSizes",b1,&b1,&flg6));
    if (flg || flg2 || flg3 || flg4 || flg5 || flg6) PetscCall(PEPCISSSetSizes(pep,i1,i2,i3,i4,i5,b1));

    PetscCall(PEPCISSGetThreshold(pep,&r1,&r2));
    PetscCall(PetscOptionsReal("-pep_ciss_delta","Threshold for numerical rank","PEPCISSSetThreshold",r1,&r1,&flg));
    PetscCall(PetscOptionsReal("-pep_ciss_spurious_threshold","Threshold for the spurious eigenpairs","PEPCISSSetThreshold",r2,&r2,&flg2));
    if (flg || flg2) PetscCall(PEPCISSSetThreshold(pep,r1,r2));

    PetscCall(PEPCISSGetRefinement(pep,&i6,&i7));
    PetscCall(PetscOptionsInt("-pep_ciss_refine_inner","Number of inner iterative refinement iterations","PEPCISSSetRefinement",i6,&i6,&flg));
    PetscCall(PetscOptionsInt("-pep_ciss_refine_blocksize","Number of blocksize iterative refinement iterations","PEPCISSSetRefinement",i7,&i7,&flg2));
    if (flg || flg2) PetscCall(PEPCISSSetRefinement(pep,i6,i7));

    PetscCall(PetscOptionsEnum("-pep_ciss_extraction","Extraction technique","PEPCISSSetExtraction",PEPCISSExtractions,(PetscEnum)ctx->extraction,(PetscEnum*)&extraction,&flg));
    if (flg) PetscCall(PEPCISSSetExtraction(pep,extraction));

  PetscOptionsHeadEnd();

  if (!pep->rg) PetscCall(PEPGetRG(pep,&pep->rg));
  PetscCall(RGSetFromOptions(pep->rg)); /* this is necessary here to set useconj */
  if (!ctx->contour || !ctx->contour->ksp) PetscCall(PEPCISSGetKSPs(pep,NULL,NULL));
  PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Something went wrong with PEPCISSGetKSPs()");
  for (i=0;i<ctx->contour->npoints;i++) PetscCall(KSPSetFromOptions(ctx->contour->ksp[i]));
  PetscCall(PetscSubcommSetFromOptions(ctx->contour->subcomm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPDestroy_CISS(PEP pep)
{
  PEP_CISS       *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  PetscCall(SlepcContourDataDestroy(&ctx->contour));
  PetscCall(PetscFree4(ctx->weight,ctx->omega,ctx->pp,ctx->sigma));
  PetscCall(PetscFree(pep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSSetSizes_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetSizes_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSSetThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSSetRefinement_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetRefinement_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSSetExtraction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetExtraction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetKSPs_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPView_CISS(PEP pep,PetscViewer viewer)
{
  PEP_CISS       *ctx = (PEP_CISS*)pep->data;
  PetscBool      isascii;
  PetscViewer    sviewer;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  sizes { integration points: %" PetscInt_FMT ", block size: %" PetscInt_FMT ", moment size: %" PetscInt_FMT ", partitions: %" PetscInt_FMT ", maximum block size: %" PetscInt_FMT " }\n",ctx->N,ctx->L,ctx->M,ctx->npart,ctx->L_max));
    if (ctx->isreal) PetscCall(PetscViewerASCIIPrintf(viewer,"  exploiting symmetry of integration points\n"));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  threshold { delta: %g, spurious threshold: %g }\n",(double)ctx->delta,(double)ctx->spurious_threshold));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  iterative refinement  { inner: %" PetscInt_FMT ", blocksize: %" PetscInt_FMT " }\n",ctx->refine_inner, ctx->refine_blocksize));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  extraction: %s\n",PEPCISSExtractions[ctx->extraction]));
    if (!ctx->contour || !ctx->contour->ksp) PetscCall(PEPCISSGetKSPs(pep,NULL,NULL));
    PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Something went wrong with PEPCISSGetKSPs()");
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

static PetscErrorCode PEPSetDSType_CISS(PEP pep)
{
  PEP_CISS       *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  switch (ctx->extraction) {
    case PEP_CISS_EXTRACTION_RITZ:   PetscCall(DSSetType(pep->ds,DSPEP)); break;
    case PEP_CISS_EXTRACTION_HANKEL: PetscCall(DSSetType(pep->ds,DSGNHEP)); break;
    case PEP_CISS_EXTRACTION_CAA:    PetscCall(DSSetType(pep->ds,DSNHEP)); break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_CISS(PEP pep)
{
  PEP_CISS       *ctx = (PEP_CISS*)pep->data;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  pep->data = ctx;
  /* set default values of parameters */
  ctx->N                  = 32;
  ctx->L                  = 16;
  ctx->M                  = ctx->N/4;
  ctx->delta              = SLEPC_DEFAULT_TOL*1e-4;
  ctx->L_max              = 64;
  ctx->spurious_threshold = PetscSqrtReal(SLEPC_DEFAULT_TOL);
  ctx->isreal             = PETSC_FALSE;
  ctx->npart              = 1;

  pep->ops->solve          = PEPSolve_CISS;
  pep->ops->setup          = PEPSetUp_CISS;
  pep->ops->setfromoptions = PEPSetFromOptions_CISS;
  pep->ops->reset          = PEPReset_CISS;
  pep->ops->destroy        = PEPDestroy_CISS;
  pep->ops->view           = PEPView_CISS;
  pep->ops->setdstype      = PEPSetDSType_CISS;

  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSSetSizes_C",PEPCISSSetSizes_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetSizes_C",PEPCISSGetSizes_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSSetThreshold_C",PEPCISSSetThreshold_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetThreshold_C",PEPCISSGetThreshold_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSSetRefinement_C",PEPCISSSetRefinement_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetRefinement_C",PEPCISSGetRefinement_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSSetExtraction_C",PEPCISSSetExtraction_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetExtraction_C",PEPCISSGetExtraction_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPCISSGetKSPs_C",PEPCISSGetKSPs_CISS));
  PetscFunctionReturn(PETSC_SUCCESS);
}
