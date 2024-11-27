/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the FEAST solver in MKL
*/

#include <petscsys.h>
#if defined(PETSC_HAVE_MKL_INTEL_ILP64)
#define MKL_ILP64
#endif
#include <mkl.h>
#include <slepc/private/epsimpl.h>        /*I "slepceps.h" I*/

#if defined(PETSC_USE_COMPLEX)
#  if defined(PETSC_USE_REAL_SINGLE)
#    define FEAST_RCI cfeast_hrci
#    define SCALAR_CAST (MKL_Complex8*)
#  else
#    define FEAST_RCI zfeast_hrci
#    define SCALAR_CAST (MKL_Complex16*)
#  endif
#else
#  if defined(PETSC_USE_REAL_SINGLE)
#    define FEAST_RCI sfeast_srci
#  else
#    define FEAST_RCI dfeast_srci
#  endif
#  define SCALAR_CAST
#endif

typedef struct {
  PetscInt      npoints;          /* number of contour points */
  PetscScalar   *work1,*Aq,*Bq;   /* workspace */
#if defined(PETSC_USE_REAL_SINGLE)
  MKL_Complex8  *work2;
#else
  MKL_Complex16 *work2;
#endif
} EPS_FEAST;

static PetscErrorCode EPSSetUp_FEAST(EPS eps)
{
  PetscInt       ncv;
  EPS_FEAST      *ctx = (EPS_FEAST*)eps->data;
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)eps),&size));
  PetscCheck(size==1,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The FEAST interface is supported for sequential runs only");
  EPSCheckHermitianDefinite(eps);
  EPSCheckSinvertCayley(eps);
  if (eps->nev==0) eps->nev = 1;
  if (eps->ncv!=PETSC_DETERMINE) {
    PetscCheck(eps->ncv>=eps->nev+2,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The value of ncv must be at least nev+2");
  } else eps->ncv = PetscMin(PetscMax(20,2*eps->nev+1),eps->n); /* set default value of ncv */
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = 20;
  if (!eps->which) eps->which = EPS_ALL;
  PetscCheck(eps->which==EPS_ALL && eps->inta!=eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver must be used with a computational interval");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_STOPPING | EPS_FEATURE_TWOSIDED);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION);

  if (!ctx->npoints) ctx->npoints = 8;

  ncv = eps->ncv;
  PetscCall(PetscFree4(ctx->work1,ctx->work2,ctx->Aq,ctx->Bq));
  PetscCall(PetscMalloc4(eps->nloc*ncv,&ctx->work1,eps->nloc*ncv,&ctx->work2,ncv*ncv,&ctx->Aq,ncv*ncv,&ctx->Bq));

  PetscCall(EPSAllocateSolution(eps,0));
  PetscCall(EPSSetWorkVecs(eps,2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_FEAST(EPS eps)
{
  EPS_FEAST      *ctx = (EPS_FEAST*)eps->data;
  MKL_INT        fpm[128],ijob,n,ncv,nconv,loop,info;
  PetscReal      *evals,epsout=0.0;
  PetscInt       i,k,nmat,ld;
  PetscScalar    *pV,*pz,*X=NULL;
  Vec            x,y,w=eps->work[0],z=eps->work[1];
  Mat            A,B;
#if defined(PETSC_USE_REAL_SINGLE)
  MKL_Complex8   Ze;
#else
  MKL_Complex16  Ze;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_HAVE_64BIT_BLAS_INDICES)
  PetscCall(PetscBLASIntCast(eps->ncv,&ncv));
  PetscCall(PetscBLASIntCast(eps->nloc,&n));
#else
  ncv = eps->ncv;
  n   = eps->nloc;
#endif

  /* parameters */
  feastinit(fpm);
  fpm[0] = (eps->numbermonitors>0)? 1: 0;   /* runtime comments */
  fpm[1] = (MKL_INT)ctx->npoints;           /* contour points */
#if !defined(PETSC_USE_REAL_SINGLE)
  fpm[2] = -PetscLog10Real(eps->tol);       /* tolerance for trace */
#endif
  fpm[3] = (MKL_INT)eps->max_it;            /* refinement loops */
  fpm[5] = 1;                               /* second stopping criterion */
#if defined(PETSC_USE_REAL_SINGLE)
  fpm[6] = -PetscLog10Real(eps->tol);       /* tolerance for trace */
#endif

  PetscCall(PetscMalloc1(eps->ncv,&evals));
  PetscCall(BVGetLeadingDimension(eps->V,&ld));
  PetscCall(BVGetArray(eps->V,&pV));
  if (ld==n) X = pV;
  else PetscCall(PetscMalloc1(eps->ncv*n,&X));

  ijob = -1;           /* first call to reverse communication interface */
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  if (nmat>1) PetscCall(STGetMatrix(eps->st,1,&B));
  else B = NULL;
  PetscCall(MatCreateVecsEmpty(A,&x,&y));

  do {

    FEAST_RCI(&ijob,&n,&Ze,SCALAR_CAST ctx->work1,ctx->work2,SCALAR_CAST ctx->Aq,SCALAR_CAST ctx->Bq,fpm,&epsout,&loop,&eps->inta,&eps->intb,&ncv,evals,SCALAR_CAST X,&nconv,eps->errest,&info);

    PetscCheck(ncv==eps->ncv,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"FEAST changed value of ncv to %d",(int)ncv);
    if (ijob == 10) {
      /* set new quadrature point */
      PetscCall(STSetShift(eps->st,Ze.real));
    } else if (ijob == 20) {
      /* use same quadrature point and factorization for transpose solve */
    } else if (ijob == 11 || ijob == 21) {
      /* linear solve (A-sigma*B)\work2, overwrite work2 */
      for (k=0;k<ncv;k++) {
        PetscCall(VecGetArray(z,&pz));
#if defined(PETSC_USE_COMPLEX)
        for (i=0;i<eps->nloc;i++) pz[i] = PetscCMPLX(ctx->work2[eps->nloc*k+i].real,ctx->work2[eps->nloc*k+i].imag);
#else
        for (i=0;i<eps->nloc;i++) pz[i] = ctx->work2[eps->nloc*k+i].real;
#endif
        PetscCall(VecRestoreArray(z,&pz));
        if (ijob == 11) PetscCall(STMatSolve(eps->st,z,w));
        else PetscCall(STMatSolveHermitianTranspose(eps->st,z,w));
        PetscCall(VecGetArray(w,&pz));
#if defined(PETSC_USE_COMPLEX)
        for (i=0;i<eps->nloc;i++) {
          ctx->work2[eps->nloc*k+i].real = PetscRealPart(pz[i]);
          ctx->work2[eps->nloc*k+i].imag = PetscImaginaryPart(pz[i]);
        }
#else
        for (i=0;i<eps->nloc;i++) ctx->work2[eps->nloc*k+i].real = pz[i];
#endif
        PetscCall(VecRestoreArray(w,&pz));
      }
    } else if (ijob == 30 || ijob == 40) {
      /* multiplication A*V or B*V, result in work1 */
      for (k=fpm[23]-1;k<fpm[23]+fpm[24]-1;k++) {
        PetscCall(VecPlaceArray(x,&X[k*eps->nloc]));
        PetscCall(VecPlaceArray(y,&ctx->work1[k*eps->nloc]));
        if (ijob == 30) PetscCall(MatMult(A,x,y));
        else if (nmat>1) PetscCall(MatMult(B,x,y));
        else PetscCall(VecCopy(x,y));
        PetscCall(VecResetArray(x));
        PetscCall(VecResetArray(y));
      }
    } else PetscCheck(ijob==0 || ijob==-2,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Internal error in FEAST reverse communication interface (ijob=%d)",(int)ijob);

  } while (ijob);

  eps->reason = EPS_CONVERGED_TOL;
  eps->its    = loop;
  eps->nconv  = nconv;
  if (info) {
    switch (info) {
      case 1:  /* No eigenvalue has been found in the proposed search interval */
        eps->nconv = 0;
        break;
      case 2:   /* FEAST did not converge "yet" */
        eps->reason = EPS_DIVERGED_ITS;
        break;
      default:
        SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Error reported by FEAST (%d)",(int)info);
    }
  }

  for (i=0;i<eps->nconv;i++) eps->eigr[i] = evals[i];
  if (ld!=n) {
    for (i=0;i<eps->nconv;i++) PetscCall(PetscArraycpy(pV+i*ld,X+i*n,n));
    PetscCall(PetscFree(X));
  }
  PetscCall(BVRestoreArray(eps->V,&pV));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&y));
  PetscCall(PetscFree(evals));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_FEAST(EPS eps)
{
  EPS_FEAST      *ctx = (EPS_FEAST*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscFree4(ctx->work1,ctx->work2,ctx->Aq,ctx->Bq));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_FEAST(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSFEASTSetNumPoints_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSFEASTGetNumPoints_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_FEAST(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  EPS_FEAST      *ctx = (EPS_FEAST*)eps->data;
  PetscInt       n;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS FEAST Options");

    n = ctx->npoints;
    PetscCall(PetscOptionsInt("-eps_feast_num_points","Number of contour integration points","EPSFEASTSetNumPoints",n,&n,&flg));
    if (flg) PetscCall(EPSFEASTSetNumPoints(eps,n));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_FEAST(EPS eps,PetscViewer viewer)
{
  EPS_FEAST      *ctx = (EPS_FEAST*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) PetscCall(PetscViewerASCIIPrintf(viewer,"  number of contour integration points=%" PetscInt_FMT "\n",ctx->npoints));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetDefaultST_FEAST(EPS eps)
{
  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) PetscCall(STSetType(eps->st,STSINVERT));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSFEASTSetNumPoints_FEAST(EPS eps,PetscInt npoints)
{
  EPS_FEAST *ctx = (EPS_FEAST*)eps->data;

  PetscFunctionBegin;
  if (npoints == PETSC_DEFAULT || npoints == PETSC_DECIDE) ctx->npoints = 8;
  else ctx->npoints = npoints;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSFEASTSetNumPoints - Sets the number of contour integration points for
   the FEAST package.

   Logically Collective

   Input Parameters:
+  eps     - the eigenproblem solver context
-  npoints - number of contour integration points

   Options Database Key:
.  -eps_feast_num_points - Sets the number of points

   Level: advanced

.seealso: EPSFEASTGetNumPoints()
@*/
PetscErrorCode EPSFEASTSetNumPoints(EPS eps,PetscInt npoints)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,npoints,2);
  PetscTryMethod(eps,"EPSFEASTSetNumPoints_C",(EPS,PetscInt),(eps,npoints));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSFEASTGetNumPoints_FEAST(EPS eps,PetscInt *npoints)
{
  EPS_FEAST *ctx = (EPS_FEAST*)eps->data;

  PetscFunctionBegin;
  *npoints = ctx->npoints;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSFEASTGetNumPoints - Gets the number of contour integration points for
   the FEAST package.

   Not Collective

   Input Parameter:
.  eps     - the eigenproblem solver context

   Output Parameter:
.  npoints - number of contour integration points

   Level: advanced

.seealso: EPSFEASTSetNumPoints()
@*/
PetscErrorCode EPSFEASTGetNumPoints(EPS eps,PetscInt *npoints)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(npoints,2);
  PetscUseMethod(eps,"EPSFEASTGetNumPoints_C",(EPS,PetscInt*),(eps,npoints));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_FEAST(EPS eps)
{
  EPS_FEAST      *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->categ = EPS_CATEGORY_CONTOUR;

  eps->ops->solve          = EPSSolve_FEAST;
  eps->ops->setup          = EPSSetUp_FEAST;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->setfromoptions = EPSSetFromOptions_FEAST;
  eps->ops->destroy        = EPSDestroy_FEAST;
  eps->ops->reset          = EPSReset_FEAST;
  eps->ops->view           = EPSView_FEAST;
  eps->ops->setdefaultst   = EPSSetDefaultST_FEAST;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSFEASTSetNumPoints_C",EPSFEASTSetNumPoints_FEAST));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSFEASTGetNumPoints_C",EPSFEASTGetNumPoints_FEAST));
  PetscFunctionReturn(PETSC_SUCCESS);
}
