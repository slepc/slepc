/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "interpol"

   Method: Polynomial interpolation

   Algorithm:

       Uses a PEP object to solve the interpolated NEP. Currently supports
       only Chebyshev interpolation on an interval.

   References:

       [1] C. Effenberger and D. Kresser, "Chebyshev interpolation for
           nonlinear eigenvalue problems", BIT 52:933-951, 2012.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/

typedef struct {
  PEP       pep;
  PetscReal tol;       /* tolerance for norm of polynomial coefficients */
  PetscInt  maxdeg;    /* maximum degree of interpolation polynomial */
  PetscInt  deg;       /* actual degree of interpolation polynomial */
} NEP_INTERPOL;

static PetscErrorCode NEPSetUp_Interpol(NEP nep)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;
  ST             st;
  RG             rg;
  PetscReal      a,b,c,d,s,tol;
  PetscScalar    zero=0.0;
  PetscBool      flg,istrivial,trackall;
  PetscInt       its,in;

  PetscFunctionBegin;
  PetscCall(NEPSetDimensions_Default(nep,nep->nev,&nep->ncv,&nep->mpd));
  PetscCheck(nep->ncv<=nep->nev+nep->mpd,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (nep->max_it==PETSC_DETERMINE) nep->max_it = PetscMax(5000,2*nep->n/nep->ncv);
  if (!nep->which) nep->which = NEP_TARGET_MAGNITUDE;
  PetscCheck(nep->which==NEP_TARGET_MAGNITUDE,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver supports only target magnitude eigenvalues");
  NEPCheckUnsupported(nep,NEP_FEATURE_CALLBACK | NEP_FEATURE_STOPPING | NEP_FEATURE_TWOSIDED);

  /* transfer PEP options */
  if (!ctx->pep) PetscCall(NEPInterpolGetPEP(nep,&ctx->pep));
  PetscCall(PEPSetBasis(ctx->pep,PEP_BASIS_CHEBYSHEV1));
  PetscCall(PEPSetWhichEigenpairs(ctx->pep,PEP_TARGET_MAGNITUDE));
  PetscCall(PetscObjectTypeCompare((PetscObject)ctx->pep,PEPJD,&flg));
  if (!flg) {
    PetscCall(PEPGetST(ctx->pep,&st));
    PetscCall(STSetType(st,STSINVERT));
  }
  PetscCall(PEPSetDimensions(ctx->pep,nep->nev,nep->ncv,nep->mpd));
  PetscCall(PEPGetTolerances(ctx->pep,&tol,&its));
  if (tol==(PetscReal)PETSC_DETERMINE) tol = SlepcDefaultTol(nep->tol);
  if (ctx->tol==(PetscReal)PETSC_DETERMINE) ctx->tol = tol;
  if (its==PETSC_DETERMINE) its = nep->max_it;
  PetscCall(PEPSetTolerances(ctx->pep,tol,its));
  PetscCall(NEPGetTrackAll(nep,&trackall));
  PetscCall(PEPSetTrackAll(ctx->pep,trackall));

  /* transfer region options */
  PetscCall(RGIsTrivial(nep->rg,&istrivial));
  PetscCheck(!istrivial,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NEPINTERPOL requires a nontrivial region");
  PetscCall(PetscObjectTypeCompare((PetscObject)nep->rg,RGINTERVAL,&flg));
  PetscCheck(flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Only implemented for interval regions");
  PetscCall(RGIntervalGetEndpoints(nep->rg,&a,&b,&c,&d));
  PetscCheck(a>-PETSC_MAX_REAL && b<PETSC_MAX_REAL,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Only implemented for bounded intervals");
  PetscCall(PEPGetRG(ctx->pep,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
  PetscCheck(a!=b,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Only implemented for intervals on the real axis");
  s = 2.0/(b-a);
  c = c*s;
  d = d*s;
  PetscCall(RGIntervalSetEndpoints(rg,-1.0,1.0,c,d));
  PetscCall(RGCheckInside(nep->rg,1,&nep->target,&zero,&in));
  PetscCheck(in>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The target is not inside the target set");
  PetscCall(PEPSetTarget(ctx->pep,(nep->target-(a+b)/2)*s));

  PetscCall(NEPAllocateSolution(nep,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Input:
    d, number of nodes to compute
    a,b, interval extremes
  Output:
    *x, array containing the d Chebyshev nodes of the interval [a,b]
    *dct2, coefficients to compute a discrete cosine transformation (DCT-II)
*/
static PetscErrorCode ChebyshevNodes(PetscInt d,PetscReal a,PetscReal b,PetscScalar *x,PetscReal *dct2)
{
  PetscInt  j,i;
  PetscReal t;

  PetscFunctionBegin;
  for (j=0;j<d+1;j++) {
    t = ((2*j+1)*PETSC_PI)/(2*(d+1));
    x[j] = (a+b)/2.0+((b-a)/2.0)*PetscCosReal(t);
    for (i=0;i<d+1;i++) dct2[j*(d+1)+i] = PetscCosReal(i*t);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSolve_Interpol(NEP nep)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;
  Mat            *A,*P;
  PetscScalar    *x,*fx,t;
  PetscReal      *cs,a,b,s,aprox,aprox0=1.0,*matnorm;
  PetscInt       i,j,k,deg=ctx->maxdeg;
  PetscBool      hasmnorm=PETSC_FALSE;
  Vec            vr,vi=NULL;
  ST             st;

  PetscFunctionBegin;
  PetscCall(PetscMalloc4((deg+1)*(deg+1),&cs,deg+1,&x,(deg+1)*nep->nt,&fx,nep->nt,&matnorm));
  for  (j=0;j<nep->nt;j++) {
    PetscCall(MatHasOperation(nep->A[j],MATOP_NORM,&hasmnorm));
    if (!hasmnorm) break;
    PetscCall(MatNorm(nep->A[j],NORM_INFINITY,matnorm+j));
  }
  if (!hasmnorm) for (j=0;j<nep->nt;j++) matnorm[j] = 1.0;
  PetscCall(RGIntervalGetEndpoints(nep->rg,&a,&b,NULL,NULL));
  PetscCall(ChebyshevNodes(deg,a,b,x,cs));
  for (j=0;j<nep->nt;j++) {
    for (i=0;i<=deg;i++) PetscCall(FNEvaluateFunction(nep->f[j],x[i],&fx[i+j*(deg+1)]));
  }
  /* Polynomial coefficients */
  PetscCall(PetscMalloc1(deg+1,&A));
  if (nep->P) PetscCall(PetscMalloc1(deg+1,&P));
  ctx->deg = deg;
  for (k=0;k<=deg;k++) {
    PetscCall(MatDuplicate(nep->A[0],MAT_COPY_VALUES,&A[k]));
    if (nep->P) PetscCall(MatDuplicate(nep->P[0],MAT_COPY_VALUES,&P[k]));
    t = 0.0;
    for (i=0;i<deg+1;i++) t += fx[i]*cs[i*(deg+1)+k];
    t *= 2.0/(deg+1);
    if (k==0) t /= 2.0;
    aprox = matnorm[0]*PetscAbsScalar(t);
    PetscCall(MatScale(A[k],t));
    if (nep->P) PetscCall(MatScale(P[k],t));
    for (j=1;j<nep->nt;j++) {
      t = 0.0;
      for (i=0;i<deg+1;i++) t += fx[i+j*(deg+1)]*cs[i*(deg+1)+k];
      t *= 2.0/(deg+1);
      if (k==0) t /= 2.0;
      aprox += matnorm[j]*PetscAbsScalar(t);
      PetscCall(MatAXPY(A[k],t,nep->A[j],nep->mstr));
      if (nep->P) PetscCall(MatAXPY(P[k],t,nep->P[j],nep->mstrp));
    }
    if (k==0) aprox0 = aprox;
    if (k>1 && aprox/aprox0<ctx->tol) { ctx->deg = k; deg = k; break; }
  }
  PetscCall(PEPSetOperators(ctx->pep,deg+1,A));
  PetscCall(MatDestroyMatrices(deg+1,&A));
  if (nep->P) {
    PetscCall(PEPGetST(ctx->pep,&st));
    PetscCall(STSetSplitPreconditioner(st,deg+1,P,nep->mstrp));
    PetscCall(MatDestroyMatrices(deg+1,&P));
  }
  PetscCall(PetscFree4(cs,x,fx,matnorm));

  /* Solve polynomial eigenproblem */
  PetscCall(PEPSolve(ctx->pep));
  PetscCall(PEPGetConverged(ctx->pep,&nep->nconv));
  PetscCall(PEPGetIterationNumber(ctx->pep,&nep->its));
  PetscCall(PEPGetConvergedReason(ctx->pep,(PEPConvergedReason*)&nep->reason));
  PetscCall(BVSetActiveColumns(nep->V,0,nep->nconv));
  PetscCall(BVCreateVec(nep->V,&vr));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(VecDuplicate(vr,&vi));
#endif
  s = 2.0/(b-a);
  for (i=0;i<nep->nconv;i++) {
    PetscCall(PEPGetEigenpair(ctx->pep,i,&nep->eigr[i],&nep->eigi[i],vr,vi));
    nep->eigr[i] /= s;
    nep->eigr[i] += (a+b)/2.0;
    nep->eigi[i] /= s;
    PetscCall(BVInsertVec(nep->V,i,vr));
#if !defined(PETSC_USE_COMPLEX)
    if (nep->eigi[i]!=0.0) PetscCall(BVInsertVec(nep->V,++i,vi));
#endif
  }
  PetscCall(VecDestroy(&vr));
  PetscCall(VecDestroy(&vi));

  nep->state = NEP_STATE_EIGENVECTORS;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPMonitor_Interpol(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscInt       i,n;
  NEP            nep = (NEP)ctx;
  PetscReal      a,b,s;
  ST             st;

  PetscFunctionBegin;
  n = PetscMin(nest,nep->ncv);
  for (i=0;i<n;i++) {
    nep->eigr[i]   = eigr[i];
    nep->eigi[i]   = eigi[i];
    nep->errest[i] = errest[i];
  }
  PetscCall(PEPGetST(pep,&st));
  PetscCall(STBackTransform(st,n,nep->eigr,nep->eigi));
  PetscCall(RGIntervalGetEndpoints(nep->rg,&a,&b,NULL,NULL));
  s = 2.0/(b-a);
  for (i=0;i<n;i++) {
    nep->eigr[i] /= s;
    nep->eigr[i] += (a+b)/2.0;
    nep->eigi[i] /= s;
  }
  PetscCall(NEPMonitor(nep,its,nconv,nep->eigr,nep->eigi,nep->errest,nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSetFromOptions_Interpol(NEP nep,PetscOptionItems PetscOptionsObject)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;
  PetscInt       i;
  PetscBool      flg1,flg2;
  PetscReal      r;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"NEP Interpol Options");

    PetscCall(NEPInterpolGetInterpolation(nep,&r,&i));
    if (!i) i = PETSC_DETERMINE;
    PetscCall(PetscOptionsInt("-nep_interpol_interpolation_degree","Maximum degree of polynomial interpolation","NEPInterpolSetInterpolation",i,&i,&flg1));
    PetscCall(PetscOptionsReal("-nep_interpol_interpolation_tol","Tolerance for interpolation coefficients","NEPInterpolSetInterpolation",r,&r,&flg2));
    if (flg1 || flg2) PetscCall(NEPInterpolSetInterpolation(nep,r,i));

  PetscOptionsHeadEnd();

  if (!ctx->pep) PetscCall(NEPInterpolGetPEP(nep,&ctx->pep));
  PetscCall(PEPSetFromOptions(ctx->pep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPInterpolSetInterpolation_Interpol(NEP nep,PetscReal tol,PetscInt degree)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  if (tol == (PetscReal)PETSC_DETERMINE) {
    ctx->tol   = PETSC_DETERMINE;
    nep->state = NEP_STATE_INITIAL;
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    ctx->tol = tol;
  }
  if (degree == PETSC_DETERMINE) {
    ctx->maxdeg = 0;
    if (nep->state) PetscCall(NEPReset(nep));
    nep->state = NEP_STATE_INITIAL;
  } else if (degree != PETSC_CURRENT) {
    PetscCheck(degree>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of degree. Must be > 0");
    if (ctx->maxdeg != degree) {
      ctx->maxdeg = degree;
      if (nep->state) PetscCall(NEPReset(nep));
      nep->state = NEP_STATE_INITIAL;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPInterpolSetInterpolation - Sets the tolerance and maximum degree when building
   the interpolation polynomial.

   Collective

   Input Parameters:
+  nep - nonlinear eigenvalue solver
.  tol - tolerance to stop computing polynomial coefficients
-  deg - maximum degree of interpolation

   Options Database Key:
+  -nep_interpol_interpolation_tol <tol> - Sets the tolerance to stop computing polynomial coefficients
-  -nep_interpol_interpolation_degree <degree> - Sets the maximum degree of interpolation

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: NEPInterpolGetInterpolation()
@*/
PetscErrorCode NEPInterpolSetInterpolation(NEP nep,PetscReal tol,PetscInt deg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,tol,2);
  PetscValidLogicalCollectiveInt(nep,deg,3);
  PetscTryMethod(nep,"NEPInterpolSetInterpolation_C",(NEP,PetscReal,PetscInt),(nep,tol,deg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPInterpolGetInterpolation_Interpol(NEP nep,PetscReal *tol,PetscInt *deg)
{
  NEP_INTERPOL *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  if (tol) *tol = ctx->tol;
  if (deg) *deg = ctx->maxdeg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPInterpolGetInterpolation - Gets the tolerance and maximum degree when building
   the interpolation polynomial.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameters:
+  tol - tolerance to stop computing polynomial coefficients
-  deg - maximum degree of interpolation

   Level: advanced

.seealso: NEPInterpolSetInterpolation()
@*/
PetscErrorCode NEPInterpolGetInterpolation(NEP nep,PetscReal *tol,PetscInt *deg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscUseMethod(nep,"NEPInterpolGetInterpolation_C",(NEP,PetscReal*,PetscInt*),(nep,tol,deg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPInterpolSetPEP_Interpol(NEP nep,PEP pep)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)pep));
  PetscCall(PEPDestroy(&ctx->pep));
  ctx->pep = pep;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPInterpolSetPEP - Associate a polynomial eigensolver object (PEP) to the
   nonlinear eigenvalue solver.

   Collective

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  pep - the polynomial eigensolver object

   Level: advanced

.seealso: NEPInterpolGetPEP()
@*/
PetscErrorCode NEPInterpolSetPEP(NEP nep,PEP pep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(pep,PEP_CLASSID,2);
  PetscCheckSameComm(nep,1,pep,2);
  PetscTryMethod(nep,"NEPInterpolSetPEP_C",(NEP,PEP),(nep,pep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPInterpolGetPEP_Interpol(NEP nep,PEP *pep)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  if (!ctx->pep) {
    PetscCall(PEPCreate(PetscObjectComm((PetscObject)nep),&ctx->pep));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->pep,(PetscObject)nep,1));
    PetscCall(PEPSetOptionsPrefix(ctx->pep,((PetscObject)nep)->prefix));
    PetscCall(PEPAppendOptionsPrefix(ctx->pep,"nep_interpol_"));
    PetscCall(PetscObjectSetOptions((PetscObject)ctx->pep,((PetscObject)nep)->options));
    PetscCall(PEPMonitorSet(ctx->pep,PEPMonitor_Interpol,nep,NULL));
  }
  *pep = ctx->pep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPInterpolGetPEP - Retrieve the polynomial eigensolver object (PEP)
   associated with the nonlinear eigenvalue solver.

   Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  pep - the polynomial eigensolver object

   Level: advanced

.seealso: NEPInterpolSetPEP()
@*/
PetscErrorCode NEPInterpolGetPEP(NEP nep,PEP *pep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(pep,2);
  PetscUseMethod(nep,"NEPInterpolGetPEP_C",(NEP,PEP*),(nep,pep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPView_Interpol(NEP nep,PetscViewer viewer)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (!ctx->pep) PetscCall(NEPInterpolGetPEP(nep,&ctx->pep));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  polynomial degree %" PetscInt_FMT ", max=%" PetscInt_FMT "\n",ctx->deg,ctx->maxdeg));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  tolerance for norm of polynomial coefficients %g\n",(double)ctx->tol));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(PEPView(ctx->pep,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPReset_Interpol(NEP nep)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  PetscCall(PEPReset(ctx->pep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDestroy_Interpol(NEP nep)
{
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  PetscCall(PEPDestroy(&ctx->pep));
  PetscCall(PetscFree(nep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolSetInterpolation_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolGetInterpolation_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolSetPEP_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolGetPEP_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_Interpol(NEP nep)
{
  NEP_INTERPOL   *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  nep->data   = (void*)ctx;
  ctx->maxdeg = 5;
  ctx->tol    = PETSC_DETERMINE;

  nep->ops->solve          = NEPSolve_Interpol;
  nep->ops->setup          = NEPSetUp_Interpol;
  nep->ops->setfromoptions = NEPSetFromOptions_Interpol;
  nep->ops->reset          = NEPReset_Interpol;
  nep->ops->destroy        = NEPDestroy_Interpol;
  nep->ops->view           = NEPView_Interpol;

  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolSetInterpolation_C",NEPInterpolSetInterpolation_Interpol));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolGetInterpolation_C",NEPInterpolGetInterpolation_Interpol));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolSetPEP_C",NEPInterpolSetPEP_Interpol));
  PetscCall(PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolGetPEP_C",NEPInterpolGetPEP_Interpol));
  PetscFunctionReturn(PETSC_SUCCESS);
}
