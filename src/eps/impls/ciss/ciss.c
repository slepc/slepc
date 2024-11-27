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

       [1] T. Sakurai and H. Sugiura, "A projection method for generalized
           eigenvalue problems", J. Comput. Appl. Math. 159:119-128, 2003.

       [2] T. Sakurai and H. Tadano, "CIRR: a Rayleigh-Ritz type method with
           contour integral for generalized eigenvalue problems", Hokkaido
           Math. J. 36:745-757, 2007.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <slepc/private/slepccontour.h>
#include <slepcblaslapack.h>

typedef struct {
  /* user parameters */
  PetscInt          N;          /* number of integration points (32) */
  PetscInt          L;          /* block size (16) */
  PetscInt          M;          /* moment degree (N/4 = 4) */
  PetscReal         delta;      /* threshold of singular value (1e-12) */
  PetscInt          L_max;      /* maximum number of columns of the source matrix V */
  PetscReal         spurious_threshold; /* discard spurious eigenpairs */
  PetscBool         isreal;     /* A and B are real */
  PetscInt          npart;      /* number of partitions */
  PetscInt          refine_inner;
  PetscInt          refine_blocksize;
  EPSCISSQuadRule   quad;
  EPSCISSExtraction extraction;
  PetscBool         usest;
  /* private data */
  SlepcContourData  contour;
  PetscReal         *sigma;     /* threshold for numerical rank */
  PetscScalar       *weight;
  PetscScalar       *omega;
  PetscScalar       *pp;
  BV                V;
  BV                S;
  BV                pV;
  BV                Y;
  PetscBool         useconj;
  PetscBool         usest_set;  /* whether the user set the usest flag or not */
  PetscObjectId     rgid;
  PetscObjectState  rgstate;
} EPS_CISS;

/*
  Set up KSP solvers for every integration point, only called if !ctx->usest
*/
static PetscErrorCode EPSCISSSetUp(EPS eps,Mat A,Mat B,Mat Pa,Mat Pb)
{
  EPS_CISS         *ctx = (EPS_CISS*)eps->data;
  SlepcContourData contour;
  PetscInt         i,p_id,nsplit;
  Mat              Amat,Pmat;
  MatStructure     str,strp;

  PetscFunctionBegin;
  if (!ctx->contour || !ctx->contour->ksp) PetscCall(EPSCISSGetKSPs(eps,NULL,NULL));
  PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Something went wrong with EPSCISSGetKSPs()");
  contour = ctx->contour;
  PetscCall(STGetMatStructure(eps->st,&str));
  PetscCall(STGetSplitPreconditionerInfo(eps->st,&nsplit,&strp));
  for (i=0;i<contour->npoints;i++) {
    p_id = i*contour->subcomm->n + contour->subcomm->color;
    PetscCall(MatDuplicate(A,MAT_COPY_VALUES,&Amat));
    if (B) PetscCall(MatAXPY(Amat,-ctx->omega[p_id],B,str));
    else PetscCall(MatShift(Amat,-ctx->omega[p_id]));
    if (nsplit) {
      PetscCall(MatDuplicate(Pa,MAT_COPY_VALUES,&Pmat));
      if (Pb) PetscCall(MatAXPY(Pmat,-ctx->omega[p_id],Pb,strp));
      else PetscCall(MatShift(Pmat,-ctx->omega[p_id]));
    } else Pmat = Amat;
    PetscCall(EPS_KSPSetOperators(contour->ksp[i],Amat,Amat));
    PetscCall(MatDestroy(&Amat));
    if (nsplit) PetscCall(MatDestroy(&Pmat));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Y_i = (A-z_i B)^{-1}BV for every integration point, Y=[Y_i] is in the context
*/
static PetscErrorCode EPSCISSSolve(EPS eps,Mat B,BV V,PetscInt L_start,PetscInt L_end)
{
  EPS_CISS         *ctx = (EPS_CISS*)eps->data;
  SlepcContourData contour;
  PetscInt         i,p_id;
  Mat              MV,BMV=NULL,MC;
  KSP              ksp;

  PetscFunctionBegin;
  if (!ctx->contour || !ctx->contour->ksp) PetscCall(EPSCISSGetKSPs(eps,NULL,NULL));
  contour = ctx->contour;
  PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Something went wrong with EPSCISSGetKSPs()");
  PetscCall(BVSetActiveColumns(V,L_start,L_end));
  PetscCall(BVGetMat(V,&MV));
  for (i=0;i<contour->npoints;i++) {
    p_id = i*contour->subcomm->n + contour->subcomm->color;
    if (ctx->usest)  {
      PetscCall(STSetShift(eps->st,ctx->omega[p_id]));
      PetscCall(STGetKSP(eps->st,&ksp));
    } else ksp = contour->ksp[i];
    PetscCall(BVSetActiveColumns(ctx->Y,i*ctx->L+L_start,i*ctx->L+L_end));
    PetscCall(BVGetMat(ctx->Y,&MC));
    if (B) {
      if (!i) {
        PetscCall(MatProductCreate(B,MV,NULL,&BMV));
        PetscCall(MatProductSetType(BMV,MATPRODUCT_AB));
        PetscCall(MatProductSetFromOptions(BMV));
        PetscCall(MatProductSymbolic(BMV));
      }
      PetscCall(MatProductNumeric(BMV));
      PetscCall(KSPMatSolve(ksp,BMV,MC));
    } else PetscCall(KSPMatSolve(ksp,MV,MC));
    PetscCall(BVRestoreMat(ctx->Y,&MC));
    if (ctx->usest && i<contour->npoints-1) PetscCall(KSPReset(ksp));
  }
  PetscCall(MatDestroy(&BMV));
  PetscCall(BVRestoreMat(V,&MV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode rescale_eig(EPS eps,PetscInt nv)
{
  EPS_CISS       *ctx = (EPS_CISS*)eps->data;
  PetscInt       i;
  PetscScalar    center;
  PetscReal      radius,a,b,c,d,rgscale;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      start_ang,end_ang,vscale,theta;
#endif
  PetscBool      isring,isellipse,isinterval;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->rg,RGELLIPSE,&isellipse));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->rg,RGRING,&isring));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->rg,RGINTERVAL,&isinterval));
  PetscCall(RGGetScale(eps->rg,&rgscale));
  if (isinterval) {
    PetscCall(RGIntervalGetEndpoints(eps->rg,NULL,NULL,&c,&d));
    if (c==d) {
      for (i=0;i<nv;i++) {
#if defined(PETSC_USE_COMPLEX)
        eps->eigr[i] = PetscRealPart(eps->eigr[i]);
#else
        eps->eigi[i] = 0;
#endif
      }
    }
  }
  if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) {
    if (isellipse) {
      PetscCall(RGEllipseGetParameters(eps->rg,&center,&radius,NULL));
      for (i=0;i<nv;i++) eps->eigr[i] = rgscale*(center + radius*eps->eigr[i]);
    } else if (isinterval) {
      PetscCall(RGIntervalGetEndpoints(eps->rg,&a,&b,&c,&d));
      if (ctx->quad == EPS_CISS_QUADRULE_CHEBYSHEV) {
        for (i=0;i<nv;i++) {
          if (c==d) eps->eigr[i] = ((b-a)*(eps->eigr[i]+1.0)/2.0+a)*rgscale;
          if (a==b) {
#if defined(PETSC_USE_COMPLEX)
            eps->eigr[i] = ((d-c)*(eps->eigr[i]+1.0)/2.0+c)*rgscale*PETSC_i;
#else
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Integration points on a vertical line require complex arithmetic");
#endif
          }
        }
      } else {
        center = (b+a)/2.0+(d+c)/2.0*PETSC_PI;
        radius = PetscSqrtReal(PetscPowRealInt((b-a)/2.0,2)+PetscPowRealInt((d-c)/2.0,2));
        for (i=0;i<nv;i++) eps->eigr[i] = center + radius*eps->eigr[i];
      }
    } else if (isring) {  /* only supported in complex scalars */
#if defined(PETSC_USE_COMPLEX)
      PetscCall(RGRingGetParameters(eps->rg,&center,&radius,&vscale,&start_ang,&end_ang,NULL));
      if (ctx->quad == EPS_CISS_QUADRULE_CHEBYSHEV) {
        for (i=0;i<nv;i++) {
          theta = (start_ang*2.0+(end_ang-start_ang)*(PetscRealPart(eps->eigr[i])+1.0))*PETSC_PI;
          eps->eigr[i] = rgscale*center + (rgscale*radius+PetscImaginaryPart(eps->eigr[i]))*PetscCMPLX(PetscCosReal(theta),vscale*PetscSinReal(theta));
        }
      } else {
        for (i=0;i<nv;i++) eps->eigr[i] = rgscale*(center + radius*eps->eigr[i]);
      }
#endif
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUp_CISS(EPS eps)
{
  EPS_CISS         *ctx = (EPS_CISS*)eps->data;
  SlepcContourData contour;
  PetscBool        istrivial,isring,isellipse,isinterval,flg;
  PetscReal        c,d;
  PetscInt         nsplit;
  PetscRandom      rand;
  PetscObjectId    id;
  PetscObjectState state;
  Mat              A[2],Psplit[2];
  Vec              v0;

  PetscFunctionBegin;
  EPSCheckNotStructured(eps);
  if (eps->ncv==PETSC_DETERMINE) {
    if (eps->nev==0) eps->nev = 1;
    eps->ncv = ctx->L_max*ctx->M;
    if (eps->ncv>eps->n) {
      eps->ncv = eps->n;
      ctx->L_max = eps->ncv/ctx->M;
      PetscCheck(ctx->L_max,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Cannot adjust solver parameters, try setting a smaller value of M (moment size)");
    }
  } else {
    PetscCall(EPSSetDimensions_Default(eps,&eps->nev,&eps->ncv,&eps->mpd));
    ctx->L_max = eps->ncv/ctx->M;
    if (!ctx->L_max) {
      ctx->L_max = 1;
      eps->ncv = ctx->L_max*ctx->M;
    }
  }
  ctx->L = PetscMin(ctx->L,ctx->L_max);
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = 5;
  if (eps->mpd==PETSC_DETERMINE) eps->mpd = eps->ncv;
  if (!eps->which) eps->which = EPS_ALL;
  PetscCheck(eps->which==EPS_ALL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver supports only computing all eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_EXTRACTION | EPS_FEATURE_STOPPING | EPS_FEATURE_TWOSIDED);

  /* check region */
  PetscCall(RGIsTrivial(eps->rg,&istrivial));
  PetscCheck(!istrivial,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"CISS requires a nontrivial region, e.g. -rg_type ellipse ...");
  PetscCall(RGGetComplement(eps->rg,&flg));
  PetscCheck(!flg,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"A region with complement flag set is not allowed");
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->rg,RGELLIPSE,&isellipse));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->rg,RGRING,&isring));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->rg,RGINTERVAL,&isinterval));
  PetscCheck(isellipse || isring || isinterval,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Currently only implemented for interval, elliptic or ring regions");

  /* if the region has changed, then reset contour data */
  PetscCall(PetscObjectGetId((PetscObject)eps->rg,&id));
  PetscCall(PetscObjectStateGet((PetscObject)eps->rg,&state));
  if (ctx->rgid && (id != ctx->rgid || state != ctx->rgstate)) {
    PetscCall(SlepcContourDataDestroy(&ctx->contour));
    PetscCall(PetscInfo(eps,"Resetting the contour data structure due to a change of region\n"));
    ctx->rgid = id; ctx->rgstate = state;
  }

#if !defined(PETSC_USE_COMPLEX)
  PetscCheck(!isring,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Ring region only supported for complex scalars");
#endif
  if (isinterval) {
    PetscCall(RGIntervalGetEndpoints(eps->rg,NULL,NULL,&c,&d));
#if !defined(PETSC_USE_COMPLEX)
    PetscCheck(c==d && c==0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"In real scalars, endpoints of the imaginary axis must be both zero");
#endif
    if (!ctx->quad && c==d) ctx->quad = EPS_CISS_QUADRULE_CHEBYSHEV;
  }
  if (!ctx->quad) ctx->quad = EPS_CISS_QUADRULE_TRAPEZOIDAL;

  /* create contour data structure */
  if (!ctx->contour) {
    PetscCall(RGCanUseConjugates(eps->rg,ctx->isreal,&ctx->useconj));
    PetscCall(SlepcContourDataCreate(ctx->useconj?ctx->N/2:ctx->N,ctx->npart,(PetscObject)eps,&ctx->contour));
  }

  PetscCall(EPSAllocateSolution(eps,0));
  PetscCall(BVGetRandomContext(eps->V,&rand));  /* make sure the random context is available when duplicating */
  if (ctx->weight) PetscCall(PetscFree4(ctx->weight,ctx->omega,ctx->pp,ctx->sigma));
  PetscCall(PetscMalloc4(ctx->N,&ctx->weight,ctx->N+1,&ctx->omega,ctx->N,&ctx->pp,ctx->L_max*ctx->M,&ctx->sigma));

  /* allocate basis vectors */
  PetscCall(BVDestroy(&ctx->S));
  PetscCall(BVDuplicateResize(eps->V,ctx->L*ctx->M,&ctx->S));
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVDuplicateResize(eps->V,ctx->L,&ctx->V));

  PetscCall(STGetMatrix(eps->st,0,&A[0]));
  PetscCall(MatIsShell(A[0],&flg));
  PetscCheck(!flg,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Matrix type shell is not supported in this solver");
  if (eps->isgeneralized) PetscCall(STGetMatrix(eps->st,1,&A[1]));

  if (!ctx->usest_set) ctx->usest = (ctx->npart>1)? PETSC_FALSE: PETSC_TRUE;
  PetscCheck(!ctx->usest || ctx->npart==1,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The usest flag is not supported when partitions > 1");

  /* check if a user-defined split preconditioner has been set */
  PetscCall(STGetSplitPreconditionerInfo(eps->st,&nsplit,NULL));
  if (nsplit) {
    PetscCall(STGetSplitPreconditionerTerm(eps->st,0,&Psplit[0]));
    if (eps->isgeneralized) PetscCall(STGetSplitPreconditionerTerm(eps->st,1,&Psplit[1]));
  }

  contour = ctx->contour;
  PetscCall(SlepcContourRedundantMat(contour,eps->isgeneralized?2:1,A,nsplit?Psplit:NULL));
  if (contour->pA) {
    PetscCall(BVGetColumn(ctx->V,0,&v0));
    PetscCall(SlepcContourScatterCreate(contour,v0));
    PetscCall(BVRestoreColumn(ctx->V,0,&v0));
    PetscCall(BVDestroy(&ctx->pV));
    PetscCall(BVCreate(PetscObjectComm((PetscObject)contour->xsub),&ctx->pV));
    PetscCall(BVSetSizesFromVec(ctx->pV,contour->xsub,eps->n));
    PetscCall(BVSetFromOptions(ctx->pV));
    PetscCall(BVResize(ctx->pV,ctx->L,PETSC_FALSE));
  }

  EPSCheckDefinite(eps);
  EPSCheckSinvertCondition(eps,ctx->usest," (with the usest flag set)");

  PetscCall(BVDestroy(&ctx->Y));
  if (contour->pA) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)contour->xsub),&ctx->Y));
    PetscCall(BVSetSizesFromVec(ctx->Y,contour->xsub,eps->n));
    PetscCall(BVSetFromOptions(ctx->Y));
    PetscCall(BVResize(ctx->Y,contour->npoints*ctx->L,PETSC_FALSE));
  } else PetscCall(BVDuplicateResize(eps->V,contour->npoints*ctx->L,&ctx->Y));

  if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) PetscCall(DSSetType(eps->ds,DSGNHEP));
  else if (eps->isgeneralized) {
    if (eps->ishermitian && eps->ispositive) PetscCall(DSSetType(eps->ds,DSGHEP));
    else PetscCall(DSSetType(eps->ds,DSGNHEP));
  } else {
    if (eps->ishermitian) PetscCall(DSSetType(eps->ds,DSHEP));
    else PetscCall(DSSetType(eps->ds,DSNHEP));
  }
  PetscCall(DSAllocate(eps->ds,eps->ncv));

#if !defined(PETSC_USE_COMPLEX)
  PetscCall(EPSSetWorkVecs(eps,3));
  if (!eps->ishermitian) PetscCall(PetscInfo(eps,"Warning: complex eigenvalues are not calculated exactly without --with-scalar-type=complex in PETSc\n"));
#else
  PetscCall(EPSSetWorkVecs(eps,2));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUpSort_CISS(EPS eps)
{
  SlepcSC        sc;

  PetscFunctionBegin;
  /* fill sorting criterion context */
  eps->sc->comparison    = SlepcCompareSmallestReal;
  eps->sc->comparisonctx = NULL;
  eps->sc->map           = NULL;
  eps->sc->mapobj        = NULL;

  /* fill sorting criterion for DS */
  PetscCall(DSGetSlepcSC(eps->ds,&sc));
  sc->comparison    = SlepcCompareLargestMagnitude;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_CISS(EPS eps)
{
  EPS_CISS         *ctx = (EPS_CISS*)eps->data;
  SlepcContourData contour = ctx->contour;
  Mat              A,B,X,M,pA,pB,T,J,Pa=NULL,Pb=NULL;
  BV               V;
  PetscInt         i,j,ld,nmat,L_add=0,nv=0,L_base=ctx->L,inner,nlocal,*inside,nsplit;
  PetscScalar      *Mu,*H0,*H1=NULL,*rr,*temp;
  PetscReal        error,max_error,norm;
  PetscBool        *fl1;
  Vec              si,si1=NULL,w[3];
  PetscRandom      rand;
#if defined(PETSC_USE_COMPLEX)
  PetscBool        isellipse;
  PetscReal        est_eig,eta;
#else
  PetscReal        normi;
#endif

  PetscFunctionBegin;
  w[0] = eps->work[0];
#if defined(PETSC_USE_COMPLEX)
  w[1] = NULL;
#else
  w[1] = eps->work[2];
#endif
  w[2] = eps->work[1];
  PetscCall(VecGetLocalSize(w[0],&nlocal));
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(RGComputeQuadrature(eps->rg,ctx->quad==EPS_CISS_QUADRULE_CHEBYSHEV?RG_QUADRULE_CHEBYSHEV:RG_QUADRULE_TRAPEZOIDAL,ctx->N,ctx->omega,ctx->pp,ctx->weight));
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  if (nmat>1) PetscCall(STGetMatrix(eps->st,1,&B));
  else B = NULL;
  J = (contour->pA && nmat>1)? contour->pA[1]: B;
  V = contour->pA? ctx->pV: ctx->V;
  if (!ctx->usest) {
    T = contour->pA? contour->pA[0]: A;
    PetscCall(STGetSplitPreconditionerInfo(eps->st,&nsplit,NULL));
    if (nsplit) {
      if (contour->pA) {
        Pa = contour->pP[0];
        if (nsplit>1) Pb = contour->pP[1];
      } else {
        PetscCall(STGetSplitPreconditionerTerm(eps->st,0,&Pa));
        if (nsplit>1) PetscCall(STGetSplitPreconditionerTerm(eps->st,1,&Pb));
      }
    }
    PetscCall(EPSCISSSetUp(eps,T,J,Pa,Pb));
  }
  PetscCall(BVSetActiveColumns(ctx->V,0,ctx->L));
  PetscCall(BVSetRandomSign(ctx->V));
  PetscCall(BVGetRandomContext(ctx->V,&rand));

  if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
  PetscCall(EPSCISSSolve(eps,J,V,0,ctx->L));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->rg,RGELLIPSE,&isellipse));
  if (isellipse) {
    PetscCall(BVTraceQuadrature(ctx->Y,ctx->V,ctx->L,ctx->L,ctx->weight,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj,&est_eig));
    PetscCall(PetscInfo(eps,"Estimated eigenvalue count: %f\n",(double)est_eig));
    eta = PetscPowReal(10.0,-PetscLog10Real(eps->tol)/ctx->N);
    L_add = PetscMax(0,(PetscInt)PetscCeilReal((est_eig*eta)/ctx->M)-ctx->L);
    if (L_add>ctx->L_max-ctx->L) {
      PetscCall(PetscInfo(eps,"Number of eigenvalues inside the contour path may be too large\n"));
      L_add = ctx->L_max-ctx->L;
    }
  }
#endif
  if (L_add>0) {
    PetscCall(PetscInfo(eps,"Changing L %" PetscInt_FMT " -> %" PetscInt_FMT " by Estimate #Eig\n",ctx->L,ctx->L+L_add));
    PetscCall(BVCISSResizeBases(ctx->S,contour->pA?ctx->pV:ctx->V,ctx->Y,ctx->L,ctx->L+L_add,ctx->M,contour->npoints));
    PetscCall(BVSetActiveColumns(ctx->V,ctx->L,ctx->L+L_add));
    PetscCall(BVSetRandomSign(ctx->V));
    if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
    ctx->L += L_add;
    PetscCall(EPSCISSSolve(eps,J,V,ctx->L-L_add,ctx->L));
  }
  PetscCall(PetscMalloc2(ctx->L*ctx->L*ctx->M*2,&Mu,ctx->L*ctx->M*ctx->L*ctx->M,&H0));
  for (i=0;i<ctx->refine_blocksize;i++) {
    PetscCall(BVDotQuadrature(ctx->Y,(contour->pA)?ctx->pV:ctx->V,Mu,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->subcomm,contour->npoints,ctx->useconj));
    PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
    PetscCall(PetscLogEventBegin(EPS_CISS_SVD,eps,0,0,0));
    PetscCall(SlepcCISS_BH_SVD(H0,ctx->L*ctx->M,ctx->delta,ctx->sigma,&nv));
    PetscCall(PetscLogEventEnd(EPS_CISS_SVD,eps,0,0,0));
    if (ctx->sigma[0]<=ctx->delta || nv < ctx->L*ctx->M || ctx->L == ctx->L_max) break;
    L_add = L_base;
    if (ctx->L+L_add>ctx->L_max) L_add = ctx->L_max-ctx->L;
    PetscCall(PetscInfo(eps,"Changing L %" PetscInt_FMT " -> %" PetscInt_FMT " by SVD(H0)\n",ctx->L,ctx->L+L_add));
    PetscCall(BVCISSResizeBases(ctx->S,contour->pA?ctx->pV:ctx->V,ctx->Y,ctx->L,ctx->L+L_add,ctx->M,contour->npoints));
    PetscCall(BVSetActiveColumns(ctx->V,ctx->L,ctx->L+L_add));
    PetscCall(BVSetRandomSign(ctx->V));
    if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
    ctx->L += L_add;
    PetscCall(EPSCISSSolve(eps,J,V,ctx->L-L_add,ctx->L));
    if (L_add) {
      PetscCall(PetscFree2(Mu,H0));
      PetscCall(PetscMalloc2(ctx->L*ctx->L*ctx->M*2,&Mu,ctx->L*ctx->M*ctx->L*ctx->M,&H0));
    }
  }
  if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) PetscCall(PetscMalloc1(ctx->L*ctx->M*ctx->L*ctx->M,&H1));

  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;
    for (inner=0;inner<=ctx->refine_inner;inner++) {
      if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) {
        PetscCall(BVDotQuadrature(ctx->Y,(contour->pA)?ctx->pV:ctx->V,Mu,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
        PetscCall(PetscLogEventBegin(EPS_CISS_SVD,eps,0,0,0));
        PetscCall(SlepcCISS_BH_SVD(H0,ctx->L*ctx->M,ctx->delta,ctx->sigma,&nv));
        PetscCall(PetscLogEventEnd(EPS_CISS_SVD,eps,0,0,0));
        break;
      } else {
        PetscCall(BVSumQuadrature(ctx->S,ctx->Y,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(BVSetActiveColumns(ctx->S,0,ctx->L));
        PetscCall(BVSetActiveColumns(ctx->V,0,ctx->L));
        PetscCall(BVCopy(ctx->S,ctx->V));
        PetscCall(BVSVDAndRank(ctx->S,ctx->M,ctx->L,ctx->delta,BV_SVD_METHOD_REFINE,H0,ctx->sigma,&nv));
        if (ctx->sigma[0]>ctx->delta && nv==ctx->L*ctx->M && inner!=ctx->refine_inner) {
          if (contour->pA) PetscCall(BVScatter(ctx->V,ctx->pV,contour->scatterin,contour->xdup));
          PetscCall(EPSCISSSolve(eps,J,V,0,ctx->L));
        } else break;
      }
    }
    eps->nconv = 0;
    if (nv == 0) eps->reason = EPS_CONVERGED_TOL;
    else {
      PetscCall(DSSetDimensions(eps->ds,nv,0,0));
      PetscCall(DSSetState(eps->ds,DS_STATE_RAW));

      if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) {
        PetscCall(CISS_BlockHankel(Mu,0,ctx->L,ctx->M,H0));
        PetscCall(CISS_BlockHankel(Mu,1,ctx->L,ctx->M,H1));
        PetscCall(DSGetArray(eps->ds,DS_MAT_A,&temp));
        for (j=0;j<nv;j++) {
          for (i=0;i<nv;i++) {
            temp[i+j*ld] = H1[i+j*ctx->L*ctx->M];
          }
        }
        PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&temp));
        PetscCall(DSGetArray(eps->ds,DS_MAT_B,&temp));
        for (j=0;j<nv;j++) {
          for (i=0;i<nv;i++) {
            temp[i+j*ld] = H0[i+j*ctx->L*ctx->M];
          }
        }
        PetscCall(DSRestoreArray(eps->ds,DS_MAT_B,&temp));
      } else {
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
        PetscCall(DSGetMat(eps->ds,DS_MAT_A,&pA));
        PetscCall(MatZeroEntries(pA));
        PetscCall(BVMatProject(ctx->S,A,ctx->S,pA));
        PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&pA));
        if (B) {
          PetscCall(DSGetMat(eps->ds,DS_MAT_B,&pB));
          PetscCall(MatZeroEntries(pB));
          PetscCall(BVMatProject(ctx->S,B,ctx->S,pB));
          PetscCall(DSRestoreMat(eps->ds,DS_MAT_B,&pB));
        }
      }

      PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
      PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

      PetscCall(PetscMalloc3(nv,&fl1,nv,&inside,nv,&rr));
      PetscCall(rescale_eig(eps,nv));
      PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));
      PetscCall(DSGetMat(eps->ds,DS_MAT_X,&X));
      PetscCall(SlepcCISS_isGhost(X,nv,ctx->sigma,ctx->spurious_threshold,fl1));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&X));
      PetscCall(RGCheckInside(eps->rg,nv,eps->eigr,eps->eigi,inside));
      for (i=0;i<nv;i++) {
        if (fl1[i] && inside[i]>=0) {
          rr[i] = 1.0;
          eps->nconv++;
        } else rr[i] = 0.0;
      }
      PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,rr,NULL,&eps->nconv));
      PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));
      PetscCall(rescale_eig(eps,nv));
      PetscCall(PetscFree3(fl1,inside,rr));
      PetscCall(BVSetActiveColumns(eps->V,0,nv));
      if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) {
        PetscCall(BVSumQuadrature(ctx->S,ctx->Y,ctx->M,ctx->L,ctx->L,ctx->weight,ctx->pp,contour->scatterin,contour->subcomm,contour->npoints,ctx->useconj));
        PetscCall(BVSetActiveColumns(ctx->S,0,ctx->L));
        PetscCall(BVCopy(ctx->S,ctx->V));
        PetscCall(BVSetActiveColumns(ctx->S,0,nv));
      }
      PetscCall(BVCopy(ctx->S,eps->V));

      PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));
      PetscCall(DSGetMat(eps->ds,DS_MAT_X,&X));
      PetscCall(BVMultInPlace(ctx->S,X,0,eps->nconv));
      if (eps->ishermitian) PetscCall(BVMultInPlace(eps->V,X,0,eps->nconv));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&X));
      max_error = 0.0;
      for (i=0;i<eps->nconv;i++) {
        PetscCall(BVGetColumn(ctx->S,i,&si));
#if !defined(PETSC_USE_COMPLEX)
        if (eps->eigi[i]!=0.0) PetscCall(BVGetColumn(ctx->S,i+1,&si1));
#endif
        PetscCall(EPSComputeResidualNorm_Private(eps,PETSC_FALSE,eps->eigr[i],eps->eigi[i],si,si1,w,&error));
        if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) {  /* vector is not normalized */
          PetscCall(VecNorm(si,NORM_2,&norm));
#if !defined(PETSC_USE_COMPLEX)
          if (eps->eigi[i]!=0.0) {
            PetscCall(VecNorm(si1,NORM_2,&normi));
            norm = SlepcAbsEigenvalue(norm,normi);
          }
#endif
          error /= norm;
        }
        PetscCall((*eps->converged)(eps,eps->eigr[i],eps->eigi[i],error,&error,eps->convergedctx));
        PetscCall(BVRestoreColumn(ctx->S,i,&si));
#if !defined(PETSC_USE_COMPLEX)
        if (eps->eigi[i]!=0.0) {
          PetscCall(BVRestoreColumn(ctx->S,i+1,&si1));
          i++;
        }
#endif
        max_error = PetscMax(max_error,error);
      }

      if (max_error <= eps->tol) eps->reason = EPS_CONVERGED_TOL;
      else if (eps->its >= eps->max_it) eps->reason = EPS_DIVERGED_ITS;
      else {
        if (eps->nconv > ctx->L) nv = eps->nconv;
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
        PetscCall(EPSCISSSolve(eps,J,V,0,ctx->L));
      }
    }
  }
  if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) PetscCall(PetscFree(H1));
  PetscCall(PetscFree2(Mu,H0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_CISS(EPS eps)
{
  EPS_CISS       *ctx = (EPS_CISS*)eps->data;
  PetscInt       n;
  Mat            Z,B=NULL;

  PetscFunctionBegin;
  if (eps->ishermitian) {
    if (eps->isgeneralized && !eps->ispositive) PetscCall(EPSComputeVectors_Indefinite(eps));
    else PetscCall(EPSComputeVectors_Hermitian(eps));
    if (eps->isgeneralized && eps->ispositive && ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) {
      /* normalize to have unit B-norm */
      PetscCall(STGetMatrix(eps->st,1,&B));
      PetscCall(BVSetMatrix(eps->V,B,PETSC_FALSE));
      PetscCall(BVNormalize(eps->V,NULL));
      PetscCall(BVSetMatrix(eps->V,NULL,PETSC_FALSE));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCall(DSGetDimensions(eps->ds,&n,NULL,NULL,NULL));
  PetscCall(BVSetActiveColumns(eps->V,0,n));

  /* right eigenvectors */
  PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));

  /* V = V * Z */
  PetscCall(DSGetMat(eps->ds,DS_MAT_X,&Z));
  PetscCall(BVMultInPlace(eps->V,Z,0,n));
  PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&Z));
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));

  /* normalize */
  if (ctx->extraction == EPS_CISS_EXTRACTION_HANKEL) PetscCall(BVNormalize(eps->V,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSSetSizes_CISS(EPS eps,PetscInt ip,PetscInt bs,PetscInt ms,PetscInt npart,PetscInt bsmax,PetscBool realmats)
{
  EPS_CISS       *ctx = (EPS_CISS*)eps->data;
  PetscInt       oN,oL,oM,oLmax,onpart;
  PetscMPIInt    size;

  PetscFunctionBegin;
  oN = ctx->N;
  if (ip == PETSC_DETERMINE) {
    if (ctx->N!=32) { ctx->N =32; ctx->M = ctx->N/4; }
  } else if (ip != PETSC_CURRENT) {
    PetscCheck(ip>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The ip argument must be > 0");
    PetscCheck(ip%2==0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The ip argument must be an even number");
    if (ctx->N!=ip) { ctx->N = ip; ctx->M = ctx->N/4; }
  }
  oL = ctx->L;
  if (bs == PETSC_DETERMINE) {
    ctx->L = 16;
  } else if (bs != PETSC_CURRENT) {
    PetscCheck(bs>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The bs argument must be > 0");
    ctx->L = bs;
  }
  oM = ctx->M;
  if (ms == PETSC_DETERMINE) {
    ctx->M = ctx->N/4;
  } else if (ms != PETSC_CURRENT) {
    PetscCheck(ms>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The ms argument must be > 0");
    PetscCheck(ms<=ctx->N,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The ms argument must be less than or equal to the number of integration points");
    ctx->M = ms;
  }
  onpart = ctx->npart;
  if (npart == PETSC_DETERMINE) {
    ctx->npart = 1;
  } else if (npart != PETSC_CURRENT) {
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)eps),&size));
    PetscCheck(npart>0 && npart<=size,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of npart");
    ctx->npart = npart;
  }
  oLmax = ctx->L_max;
  if (bsmax == PETSC_DETERMINE) {
    ctx->L_max = 64;
  } else if (bsmax != PETSC_CURRENT) {
    PetscCheck(bsmax>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The bsmax argument must be > 0");
    ctx->L_max = PetscMax(bsmax,ctx->L);
  }
  if (onpart != ctx->npart || oN != ctx->N || realmats != ctx->isreal) {
    PetscCall(SlepcContourDataDestroy(&ctx->contour));
    PetscCall(PetscInfo(eps,"Resetting the contour data structure due to a change of parameters\n"));
    eps->state = EPS_STATE_INITIAL;
  }
  ctx->isreal = realmats;
  if (oL != ctx->L || oM != ctx->M || oLmax != ctx->L_max) eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSSetSizes - Sets the values of various size parameters in the CISS solver.

   Logically Collective

   Input Parameters:
+  eps   - the eigenproblem solver context
.  ip    - number of integration points
.  bs    - block size
.  ms    - moment size
.  npart - number of partitions when splitting the communicator
.  bsmax - max block size
-  realmats - A and B are real

   Options Database Keys:
+  -eps_ciss_integration_points - Sets the number of integration points
.  -eps_ciss_blocksize - Sets the block size
.  -eps_ciss_moments - Sets the moment size
.  -eps_ciss_partitions - Sets the number of partitions
.  -eps_ciss_maxblocksize - Sets the maximum block size
-  -eps_ciss_realmats - A and B are real

   Note:
   For all integer arguments, you can use PETSC_CURRENT to keep the current value, and
   PETSC_DETERMINE to set them to a default value.

   The default number of partitions is 1. This means the internal KSP object is shared
   among all processes of the EPS communicator. Otherwise, the communicator is split
   into npart communicators, so that npart KSP solves proceed simultaneously.

   Level: advanced

.seealso: EPSCISSGetSizes()
@*/
PetscErrorCode EPSCISSSetSizes(EPS eps,PetscInt ip,PetscInt bs,PetscInt ms,PetscInt npart,PetscInt bsmax,PetscBool realmats)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,ip,2);
  PetscValidLogicalCollectiveInt(eps,bs,3);
  PetscValidLogicalCollectiveInt(eps,ms,4);
  PetscValidLogicalCollectiveInt(eps,npart,5);
  PetscValidLogicalCollectiveInt(eps,bsmax,6);
  PetscValidLogicalCollectiveBool(eps,realmats,7);
  PetscTryMethod(eps,"EPSCISSSetSizes_C",(EPS,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool),(eps,ip,bs,ms,npart,bsmax,realmats));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSGetSizes_CISS(EPS eps,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

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
   EPSCISSGetSizes - Gets the values of various size parameters in the CISS solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  ip    - number of integration points
.  bs    - block size
.  ms    - moment size
.  npart - number of partitions when splitting the communicator
.  bsmax - max block size
-  realmats - A and B are real

   Level: advanced

.seealso: EPSCISSSetSizes()
@*/
PetscErrorCode EPSCISSGetSizes(EPS eps,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSCISSGetSizes_C",(EPS,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*),(eps,ip,bs,ms,npart,bsmax,realmats));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSSetThreshold_CISS(EPS eps,PetscReal delta,PetscReal spur)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  if (delta == (PetscReal)PETSC_DETERMINE) {
    ctx->delta = SLEPC_DEFAULT_TOL*1e-4;
  } else if (delta != (PetscReal)PETSC_CURRENT) {
    PetscCheck(delta>0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The delta argument must be > 0.0");
    ctx->delta = delta;
  }
  if (spur == (PetscReal)PETSC_DETERMINE) {
    ctx->spurious_threshold = PetscSqrtReal(SLEPC_DEFAULT_TOL);
  } else if (spur != (PetscReal)PETSC_CURRENT) {
    PetscCheck(spur>0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The spurious threshold argument must be > 0.0");
    ctx->spurious_threshold = spur;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSSetThreshold - Sets the values of various threshold parameters in
   the CISS solver.

   Logically Collective

   Input Parameters:
+  eps   - the eigenproblem solver context
.  delta - threshold for numerical rank
-  spur  - spurious threshold (to discard spurious eigenpairs)

   Options Database Keys:
+  -eps_ciss_delta - Sets the delta
-  -eps_ciss_spurious_threshold - Sets the spurious threshold

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: EPSCISSGetThreshold()
@*/
PetscErrorCode EPSCISSSetThreshold(EPS eps,PetscReal delta,PetscReal spur)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,delta,2);
  PetscValidLogicalCollectiveReal(eps,spur,3);
  PetscTryMethod(eps,"EPSCISSSetThreshold_C",(EPS,PetscReal,PetscReal),(eps,delta,spur));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSGetThreshold_CISS(EPS eps,PetscReal *delta,PetscReal *spur)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  if (delta) *delta = ctx->delta;
  if (spur)  *spur = ctx->spurious_threshold;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSGetThreshold - Gets the values of various threshold parameters
   in the CISS solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  delta - threshold for numerical rank
-  spur  - spurious threshold (to discard spurious eigenpairs)

   Level: advanced

.seealso: EPSCISSSetThreshold()
@*/
PetscErrorCode EPSCISSGetThreshold(EPS eps,PetscReal *delta,PetscReal *spur)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSCISSGetThreshold_C",(EPS,PetscReal*,PetscReal*),(eps,delta,spur));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSSetRefinement_CISS(EPS eps,PetscInt inner,PetscInt blsize)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  if (inner == PETSC_DETERMINE) {
    ctx->refine_inner = 0;
  } else if (inner != PETSC_CURRENT) {
    PetscCheck(inner>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The refine inner argument must be >= 0");
    ctx->refine_inner = inner;
  }
  if (blsize == PETSC_DETERMINE) {
    ctx->refine_blocksize = 0;
  } else if (blsize != PETSC_CURRENT) {
    PetscCheck(blsize>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The refine blocksize argument must be >= 0");
    ctx->refine_blocksize = blsize;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSSetRefinement - Sets the values of various refinement parameters
   in the CISS solver.

   Logically Collective

   Input Parameters:
+  eps    - the eigenproblem solver context
.  inner  - number of iterative refinement iterations (inner loop)
-  blsize - number of iterative refinement iterations (blocksize loop)

   Options Database Keys:
+  -eps_ciss_refine_inner - Sets number of inner iterations
-  -eps_ciss_refine_blocksize - Sets number of blocksize iterations

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default of 0.

   Level: advanced

.seealso: EPSCISSGetRefinement()
@*/
PetscErrorCode EPSCISSSetRefinement(EPS eps,PetscInt inner,PetscInt blsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,inner,2);
  PetscValidLogicalCollectiveInt(eps,blsize,3);
  PetscTryMethod(eps,"EPSCISSSetRefinement_C",(EPS,PetscInt,PetscInt),(eps,inner,blsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSGetRefinement_CISS(EPS eps,PetscInt *inner,PetscInt *blsize)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  if (inner)  *inner = ctx->refine_inner;
  if (blsize) *blsize = ctx->refine_blocksize;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSGetRefinement - Gets the values of various refinement parameters
   in the CISS solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  inner  - number of iterative refinement iterations (inner loop)
-  blsize - number of iterative refinement iterations (blocksize loop)

   Level: advanced

.seealso: EPSCISSSetRefinement()
@*/
PetscErrorCode EPSCISSGetRefinement(EPS eps, PetscInt *inner, PetscInt *blsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSCISSGetRefinement_C",(EPS,PetscInt*,PetscInt*),(eps,inner,blsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSSetUseST_CISS(EPS eps,PetscBool usest)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  ctx->usest     = usest;
  ctx->usest_set = PETSC_TRUE;
  eps->state     = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSSetUseST - Sets a flag indicating that the CISS solver will
   use the ST object for the linear solves.

   Logically Collective

   Input Parameters:
+  eps    - the eigenproblem solver context
-  usest  - boolean flag to use the ST object or not

   Options Database Keys:
.  -eps_ciss_usest <bool> - whether the ST object will be used or not

   Level: advanced

.seealso: EPSCISSGetUseST()
@*/
PetscErrorCode EPSCISSSetUseST(EPS eps,PetscBool usest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,usest,2);
  PetscTryMethod(eps,"EPSCISSSetUseST_C",(EPS,PetscBool),(eps,usest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSGetUseST_CISS(EPS eps,PetscBool *usest)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  *usest = ctx->usest;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSGetUseST - Gets the flag for using the ST object
   in the CISS solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  usest - boolean flag indicating if the ST object is being used

   Level: advanced

.seealso: EPSCISSSetUseST()
@*/
PetscErrorCode EPSCISSGetUseST(EPS eps,PetscBool *usest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(usest,2);
  PetscUseMethod(eps,"EPSCISSGetUseST_C",(EPS,PetscBool*),(eps,usest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSSetQuadRule_CISS(EPS eps,EPSCISSQuadRule quad)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  if (ctx->quad != quad) {
    ctx->quad  = quad;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSSetQuadRule - Sets the quadrature rule used in the CISS solver.

   Logically Collective

   Input Parameters:
+  eps  - the eigenproblem solver context
-  quad - the quadrature rule

   Options Database Key:
.  -eps_ciss_quadrule - Sets the quadrature rule (either 'trapezoidal' or
                           'chebyshev')

   Notes:
   By default, the trapezoidal rule is used (EPS_CISS_QUADRULE_TRAPEZOIDAL).

   If the 'chebyshev' option is specified (EPS_CISS_QUADRULE_CHEBYSHEV), then
   Chebyshev points are used as quadrature points.

   Level: advanced

.seealso: EPSCISSGetQuadRule(), EPSCISSQuadRule
@*/
PetscErrorCode EPSCISSSetQuadRule(EPS eps,EPSCISSQuadRule quad)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,quad,2);
  PetscTryMethod(eps,"EPSCISSSetQuadRule_C",(EPS,EPSCISSQuadRule),(eps,quad));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSGetQuadRule_CISS(EPS eps,EPSCISSQuadRule *quad)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  *quad = ctx->quad;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSGetQuadRule - Gets the quadrature rule used in the CISS solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  quad - quadrature rule

   Level: advanced

.seealso: EPSCISSSetQuadRule() EPSCISSQuadRule
@*/
PetscErrorCode EPSCISSGetQuadRule(EPS eps,EPSCISSQuadRule *quad)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(quad,2);
  PetscUseMethod(eps,"EPSCISSGetQuadRule_C",(EPS,EPSCISSQuadRule*),(eps,quad));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSSetExtraction_CISS(EPS eps,EPSCISSExtraction extraction)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  if (ctx->extraction != extraction) {
    ctx->extraction = extraction;
    eps->state      = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSSetExtraction - Sets the extraction technique used in the CISS solver.

   Logically Collective

   Input Parameters:
+  eps        - the eigenproblem solver context
-  extraction - the extraction technique

   Options Database Key:
.  -eps_ciss_extraction - Sets the extraction technique (either 'ritz' or
                           'hankel')

   Notes:
   By default, the Rayleigh-Ritz extraction is used (EPS_CISS_EXTRACTION_RITZ).

   If the 'hankel' option is specified (EPS_CISS_EXTRACTION_HANKEL), then
   the Block Hankel method is used for extracting eigenpairs.

   Level: advanced

.seealso: EPSCISSGetExtraction(), EPSCISSExtraction
@*/
PetscErrorCode EPSCISSSetExtraction(EPS eps,EPSCISSExtraction extraction)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,extraction,2);
  PetscTryMethod(eps,"EPSCISSSetExtraction_C",(EPS,EPSCISSExtraction),(eps,extraction));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSGetExtraction_CISS(EPS eps,EPSCISSExtraction *extraction)
{
  EPS_CISS *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  *extraction = ctx->extraction;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCISSGetExtraction - Gets the extraction technique used in the CISS solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  extraction - extraction technique

   Level: advanced

.seealso: EPSCISSSetExtraction() EPSCISSExtraction
@*/
PetscErrorCode EPSCISSGetExtraction(EPS eps,EPSCISSExtraction *extraction)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(extraction,2);
  PetscUseMethod(eps,"EPSCISSGetExtraction_C",(EPS,EPSCISSExtraction*),(eps,extraction));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCISSGetKSPs_CISS(EPS eps,PetscInt *nsolve,KSP **ksp)
{
  EPS_CISS         *ctx = (EPS_CISS*)eps->data;
  SlepcContourData contour;
  PetscInt         i,nsplit;
  PC               pc;
  MPI_Comm         child;

  PetscFunctionBegin;
  if (!ctx->contour) {  /* initialize contour data structure first */
    PetscCall(RGCanUseConjugates(eps->rg,ctx->isreal,&ctx->useconj));
    PetscCall(SlepcContourDataCreate(ctx->useconj?ctx->N/2:ctx->N,ctx->npart,(PetscObject)eps,&ctx->contour));
  }
  contour = ctx->contour;
  if (!contour->ksp) {
    PetscCall(PetscMalloc1(contour->npoints,&contour->ksp));
    PetscCall(EPSGetST(eps,&eps->st));
    PetscCall(STGetSplitPreconditionerInfo(eps->st,&nsplit,NULL));
    PetscCall(PetscSubcommGetChild(contour->subcomm,&child));
    for (i=0;i<contour->npoints;i++) {
      PetscCall(KSPCreate(child,&contour->ksp[i]));
      PetscCall(PetscObjectIncrementTabLevel((PetscObject)contour->ksp[i],(PetscObject)eps,1));
      PetscCall(KSPSetOptionsPrefix(contour->ksp[i],((PetscObject)eps)->prefix));
      PetscCall(KSPAppendOptionsPrefix(contour->ksp[i],"eps_ciss_"));
      PetscCall(PetscObjectSetOptions((PetscObject)contour->ksp[i],((PetscObject)eps)->options));
      PetscCall(KSPSetErrorIfNotConverged(contour->ksp[i],PETSC_TRUE));
      PetscCall(KSPSetTolerances(contour->ksp[i],SlepcDefaultTol(eps->tol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
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
   EPSCISSGetKSPs - Retrieve the array of linear solver objects associated with
   the CISS solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver solver

   Output Parameters:
+  nsolve - number of solver objects
-  ksp - array of linear solver object

   Notes:
   The number of KSP solvers is equal to the number of integration points divided by
   the number of partitions. This value is halved in the case of real matrices with
   a region centered at the real axis.

   Level: advanced

.seealso: EPSCISSSetSizes()
@*/
PetscErrorCode EPSCISSGetKSPs(EPS eps,PetscInt *nsolve,KSP **ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSCISSGetKSPs_C",(EPS,PetscInt*,KSP**),(eps,nsolve,ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_CISS(EPS eps)
{
  EPS_CISS       *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&ctx->S));
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVDestroy(&ctx->Y));
  if (!ctx->usest) PetscCall(SlepcContourDataReset(ctx->contour));
  PetscCall(BVDestroy(&ctx->pV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_CISS(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  PetscReal         r3,r4;
  PetscInt          i,i1,i2,i3,i4,i5,i6,i7;
  PetscBool         b1,b2,flg,flg2,flg3,flg4,flg5,flg6;
  EPS_CISS          *ctx = (EPS_CISS*)eps->data;
  EPSCISSQuadRule   quad;
  EPSCISSExtraction extraction;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS CISS Options");

    PetscCall(EPSCISSGetSizes(eps,&i1,&i2,&i3,&i4,&i5,&b1));
    PetscCall(PetscOptionsInt("-eps_ciss_integration_points","Number of integration points","EPSCISSSetSizes",i1,&i1,&flg));
    PetscCall(PetscOptionsInt("-eps_ciss_blocksize","Block size","EPSCISSSetSizes",i2,&i2,&flg2));
    PetscCall(PetscOptionsInt("-eps_ciss_moments","Moment size","EPSCISSSetSizes",i3,&i3,&flg3));
    PetscCall(PetscOptionsInt("-eps_ciss_partitions","Number of partitions","EPSCISSSetSizes",i4,&i4,&flg4));
    PetscCall(PetscOptionsInt("-eps_ciss_maxblocksize","Maximum block size","EPSCISSSetSizes",i5,&i5,&flg5));
    PetscCall(PetscOptionsBool("-eps_ciss_realmats","True if A and B are real","EPSCISSSetSizes",b1,&b1,&flg6));
    if (flg || flg2 || flg3 || flg4 || flg5 || flg6) PetscCall(EPSCISSSetSizes(eps,i1,i2,i3,i4,i5,b1));

    PetscCall(EPSCISSGetThreshold(eps,&r3,&r4));
    PetscCall(PetscOptionsReal("-eps_ciss_delta","Threshold for numerical rank","EPSCISSSetThreshold",r3,&r3,&flg));
    PetscCall(PetscOptionsReal("-eps_ciss_spurious_threshold","Threshold for the spurious eigenpairs","EPSCISSSetThreshold",r4,&r4,&flg2));
    if (flg || flg2) PetscCall(EPSCISSSetThreshold(eps,r3,r4));

    PetscCall(EPSCISSGetRefinement(eps,&i6,&i7));
    PetscCall(PetscOptionsInt("-eps_ciss_refine_inner","Number of inner iterative refinement iterations","EPSCISSSetRefinement",i6,&i6,&flg));
    PetscCall(PetscOptionsInt("-eps_ciss_refine_blocksize","Number of blocksize iterative refinement iterations","EPSCISSSetRefinement",i7,&i7,&flg2));
    if (flg || flg2) PetscCall(EPSCISSSetRefinement(eps,i6,i7));

    PetscCall(EPSCISSGetUseST(eps,&b2));
    PetscCall(PetscOptionsBool("-eps_ciss_usest","Use ST for linear solves","EPSCISSSetUseST",b2,&b2,&flg));
    if (flg) PetscCall(EPSCISSSetUseST(eps,b2));

    PetscCall(PetscOptionsEnum("-eps_ciss_quadrule","Quadrature rule","EPSCISSSetQuadRule",EPSCISSQuadRules,(PetscEnum)ctx->quad,(PetscEnum*)&quad,&flg));
    if (flg) PetscCall(EPSCISSSetQuadRule(eps,quad));

    PetscCall(PetscOptionsEnum("-eps_ciss_extraction","Extraction technique","EPSCISSSetExtraction",EPSCISSExtractions,(PetscEnum)ctx->extraction,(PetscEnum*)&extraction,&flg));
    if (flg) PetscCall(EPSCISSSetExtraction(eps,extraction));

  PetscOptionsHeadEnd();

  if (!eps->rg) PetscCall(EPSGetRG(eps,&eps->rg));
  PetscCall(RGSetFromOptions(eps->rg)); /* this is necessary here to set useconj */
  if (!ctx->contour || !ctx->contour->ksp) PetscCall(EPSCISSGetKSPs(eps,NULL,NULL));
  PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Something went wrong with EPSCISSGetKSPs()");
  for (i=0;i<ctx->contour->npoints;i++) PetscCall(KSPSetFromOptions(ctx->contour->ksp[i]));
  PetscCall(PetscSubcommSetFromOptions(ctx->contour->subcomm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_CISS(EPS eps)
{
  EPS_CISS       *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  PetscCall(SlepcContourDataDestroy(&ctx->contour));
  PetscCall(PetscFree4(ctx->weight,ctx->omega,ctx->pp,ctx->sigma));
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetSizes_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetSizes_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetRefinement_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetRefinement_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetUseST_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetUseST_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetQuadRule_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetQuadRule_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetExtraction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetExtraction_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetKSPs_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_CISS(EPS eps,PetscViewer viewer)
{
  EPS_CISS       *ctx = (EPS_CISS*)eps->data;
  PetscBool      isascii;
  PetscViewer    sviewer;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  sizes { integration points: %" PetscInt_FMT ", block size: %" PetscInt_FMT ", moment size: %" PetscInt_FMT ", partitions: %" PetscInt_FMT ", maximum block size: %" PetscInt_FMT " }\n",ctx->N,ctx->L,ctx->M,ctx->npart,ctx->L_max));
    if (ctx->isreal) PetscCall(PetscViewerASCIIPrintf(viewer,"  exploiting symmetry of integration points\n"));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  threshold { delta: %g, spurious threshold: %g }\n",(double)ctx->delta,(double)ctx->spurious_threshold));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  iterative refinement { inner: %" PetscInt_FMT ", blocksize: %" PetscInt_FMT " }\n",ctx->refine_inner, ctx->refine_blocksize));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  extraction: %s\n",EPSCISSExtractions[ctx->extraction]));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  quadrature rule: %s\n",EPSCISSQuadRules[ctx->quad]));
    if (ctx->usest) PetscCall(PetscViewerASCIIPrintf(viewer,"  using ST for linear solves\n"));
    else {
      if (!ctx->contour || !ctx->contour->ksp) PetscCall(EPSCISSGetKSPs(eps,NULL,NULL));
      PetscAssert(ctx->contour && ctx->contour->ksp,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Something went wrong with EPSCISSGetKSPs()");
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
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetDefaultST_CISS(EPS eps)
{
  EPS_CISS       *ctx = (EPS_CISS*)eps->data;
  PetscBool      usest = ctx->usest;
  KSP            ksp;
  PC             pc;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) {
    if (!ctx->usest_set) usest = (ctx->npart>1)? PETSC_FALSE: PETSC_TRUE;
    if (usest) PetscCall(STSetType(eps->st,STSINVERT));
    else {
      /* we are not going to use ST, so avoid factorizing the matrix */
      PetscCall(STSetType(eps->st,STSHIFT));
      if (eps->isgeneralized) {
        PetscCall(STGetKSP(eps->st,&ksp));
        PetscCall(KSPGetPC(ksp,&pc));
        PetscCall(PCSetType(pc,PCNONE));
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_CISS(EPS eps)
{
  EPS_CISS       *ctx = (EPS_CISS*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = ctx;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_CONTOUR;

  eps->ops->solve          = EPSSolve_CISS;
  eps->ops->setup          = EPSSetUp_CISS;
  eps->ops->setupsort      = EPSSetUpSort_CISS;
  eps->ops->setfromoptions = EPSSetFromOptions_CISS;
  eps->ops->destroy        = EPSDestroy_CISS;
  eps->ops->reset          = EPSReset_CISS;
  eps->ops->view           = EPSView_CISS;
  eps->ops->computevectors = EPSComputeVectors_CISS;
  eps->ops->setdefaultst   = EPSSetDefaultST_CISS;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetSizes_C",EPSCISSSetSizes_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetSizes_C",EPSCISSGetSizes_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetThreshold_C",EPSCISSSetThreshold_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetThreshold_C",EPSCISSGetThreshold_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetRefinement_C",EPSCISSSetRefinement_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetRefinement_C",EPSCISSGetRefinement_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetUseST_C",EPSCISSSetUseST_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetUseST_C",EPSCISSGetUseST_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetQuadRule_C",EPSCISSSetQuadRule_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetQuadRule_C",EPSCISSGetQuadRule_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSSetExtraction_C",EPSCISSSetExtraction_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetExtraction_C",EPSCISSGetExtraction_CISS));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCISSGetKSPs_C",EPSCISSGetKSPs_CISS));

  /* set default values of parameters */
  ctx->N                  = 32;
  ctx->L                  = 16;
  ctx->M                  = ctx->N/4;
  ctx->delta              = SLEPC_DEFAULT_TOL*1e-4;
  ctx->L_max              = 64;
  ctx->spurious_threshold = PetscSqrtReal(SLEPC_DEFAULT_TOL);
  ctx->usest              = PETSC_TRUE;
  ctx->usest_set          = PETSC_FALSE;
  ctx->isreal             = PETSC_FALSE;
  ctx->refine_inner       = 0;
  ctx->refine_blocksize   = 0;
  ctx->npart              = 1;
  ctx->quad               = (EPSCISSQuadRule)0;
  ctx->extraction         = EPS_CISS_EXTRACTION_RITZ;
  PetscFunctionReturn(PETSC_SUCCESS);
}
