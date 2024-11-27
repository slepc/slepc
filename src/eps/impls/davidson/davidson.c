/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Skeleton of Davidson solver. Actual solvers are GD and JD.

   References:

       [1] E. Romero and J.E. Roman, "A parallel implementation of Davidson
           methods for large-scale eigenvalue problems in SLEPc", ACM Trans.
           Math. Software 40(2):13, 2014.
*/

#include "davidson.h"

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-davidson,\n"
  "   author = \"E. Romero and J. E. Roman\",\n"
  "   title = \"A parallel implementation of {Davidson} methods for large-scale eigenvalue problems in {SLEPc}\",\n"
  "   journal = \"{ACM} Trans. Math. Software\",\n"
  "   volume = \"40\",\n"
  "   number = \"2\",\n"
  "   pages = \"13:1--13:29\",\n"
  "   year = \"2014,\"\n"
  "   doi = \"https://doi.org/10.1145/2543696\"\n"
  "}\n";

PetscErrorCode EPSSetUp_XD(EPS eps)
{
  EPS_DAVIDSON   *data = (EPS_DAVIDSON*)eps->data;
  dvdDashboard   *dvd = &data->ddb;
  dvdBlackboard  b;
  PetscInt       min_size_V,bs,initv,nmat;
  Mat            A,B;
  KSP            ksp;
  PetscBool      ipB,ispositive;
  HarmType_t     harm;
  InitType_t     init;
  PetscScalar    target;

  PetscFunctionBegin;
  EPSCheckNotStructured(eps);
  /* Setup EPS options and get the problem specification */
  if (eps->nev==0) eps->nev = 1;
  bs = data->blocksize;
  if (bs <= 0) bs = 1;
  if (eps->ncv!=PETSC_DETERMINE) {
    PetscCheck(eps->ncv>=eps->nev,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The value of ncv must be at least nev");
  } else if (eps->mpd!=PETSC_DETERMINE) eps->ncv = eps->mpd + eps->nev + bs;
  else if (eps->n < 10) eps->ncv = eps->n+eps->nev+bs;
  else if (eps->nev < 500) eps->ncv = PetscMax(eps->nev,PetscMin(eps->n-bs,PetscMax(2*eps->nev,eps->nev+15))+bs);
  else eps->ncv = PetscMax(eps->nev,PetscMin(eps->n-bs,eps->nev+500)+bs);
  if (eps->mpd==PETSC_DETERMINE) eps->mpd = PetscMin(eps->n,eps->ncv);
  PetscCheck(eps->mpd<=eps->ncv,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The mpd parameter has to be less than or equal to ncv");
  PetscCheck(eps->mpd>=2,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The mpd parameter has to be greater than 2");
  if (eps->max_it == PETSC_DETERMINE) eps->max_it = PetscMax(100*eps->ncv,2*eps->n);
  if (!eps->which) eps->which = EPS_LARGEST_MAGNITUDE;
  PetscCheck(eps->nev+bs<=eps->ncv,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The value of ncv has to be greater than nev plus blocksize");
  PetscCheck(!eps->trueres,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"-eps_true_residual is disabled in this solver.");
  EPSCheckUnsupported(eps,EPS_FEATURE_REGION | EPS_FEATURE_TWOSIDED | EPS_FEATURE_THRESHOLD);

  if (!data->minv) data->minv = (eps->n && eps->n<10)? 1: PetscMin(PetscMax(bs,6),eps->mpd/2);
  min_size_V = data->minv;
  PetscCheck(min_size_V+bs<=eps->mpd,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The value of minv must be less than mpd minus blocksize");
  if (data->plusk == PETSC_DETERMINE) {
    if (eps->problem_type == EPS_GHIEP || eps->nev+bs>eps->ncv) data->plusk = 0;
    else data->plusk = 1;
  }
  if (!data->initialsize) data->initialsize = (eps->n && eps->n<10)? 1: 6;
  initv = data->initialsize;
  PetscCheck(eps->mpd>=initv,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The initv parameter has to be less than or equal to mpd");

  /* Change the default sigma to inf if necessary */
  if (eps->which == EPS_LARGEST_MAGNITUDE || eps->which == EPS_LARGEST_REAL || eps->which == EPS_LARGEST_IMAGINARY) PetscCall(STSetDefaultShift(eps->st,PETSC_MAX_REAL));

  /* Set up preconditioner */
  PetscCall(STSetUp(eps->st));

  /* Setup problem specification in dvd */
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  if (nmat>1) PetscCall(STGetMatrix(eps->st,1,&B));
  PetscCall(EPSReset_XD(eps));
  PetscCall(PetscMemzero(dvd,sizeof(dvdDashboard)));
  dvd->A = A; dvd->B = eps->isgeneralized? B: NULL;
  ispositive = eps->ispositive;
  dvd->sA = DVD_MAT_IMPLICIT | (eps->ishermitian? DVD_MAT_HERMITIAN: 0) | ((ispositive && !eps->isgeneralized) ? DVD_MAT_POS_DEF: 0);
  /* Assume -eps_hermitian means hermitian-definite in generalized problems */
  if (!ispositive && !eps->isgeneralized && eps->ishermitian) ispositive = PETSC_TRUE;
  if (!eps->isgeneralized) dvd->sB = DVD_MAT_IMPLICIT | DVD_MAT_HERMITIAN | DVD_MAT_IDENTITY | DVD_MAT_UNITARY | DVD_MAT_POS_DEF;
  else dvd->sB = DVD_MAT_IMPLICIT | (eps->ishermitian? DVD_MAT_HERMITIAN: 0) | (ispositive? DVD_MAT_POS_DEF: 0);
  ipB = (dvd->B && data->ipB && DVD_IS(dvd->sB,DVD_MAT_HERMITIAN))?PETSC_TRUE:PETSC_FALSE;
  dvd->sEP = ((!eps->isgeneralized || (eps->isgeneralized && ipB))? DVD_EP_STD: 0) | (ispositive? DVD_EP_HERMITIAN: 0) | ((eps->problem_type == EPS_GHIEP && ipB) ? DVD_EP_INDEFINITE : 0);
  if (data->ipB && !ipB) data->ipB = PETSC_FALSE;
  dvd->correctXnorm = (dvd->B && (DVD_IS(dvd->sB,DVD_MAT_HERMITIAN)||DVD_IS(dvd->sEP,DVD_EP_INDEFINITE)))?PETSC_TRUE:PETSC_FALSE;
  dvd->nev        = eps->nev;
  dvd->which      = eps->which;
  dvd->withTarget = PETSC_TRUE;
  switch (eps->which) {
    case EPS_TARGET_MAGNITUDE:
    case EPS_TARGET_IMAGINARY:
      dvd->target[0] = target = eps->target;
      dvd->target[1] = 1.0;
      break;
    case EPS_TARGET_REAL:
      dvd->target[0] = PetscRealPart(target = eps->target);
      dvd->target[1] = 1.0;
      break;
    case EPS_LARGEST_REAL:
    case EPS_LARGEST_MAGNITUDE:
    case EPS_LARGEST_IMAGINARY: /* TODO: think about this case */
      dvd->target[0] = 1.0;
      dvd->target[1] = target = 0.0;
      break;
    case EPS_SMALLEST_MAGNITUDE:
    case EPS_SMALLEST_REAL:
    case EPS_SMALLEST_IMAGINARY: /* TODO: think about this case */
      dvd->target[0] = target = 0.0;
      dvd->target[1] = 1.0;
      break;
    case EPS_WHICH_USER:
      PetscCall(STGetShift(eps->st,&target));
      dvd->target[0] = target;
      dvd->target[1] = 1.0;
      break;
    case EPS_ALL:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues");
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported value of option 'which'");
  }
  dvd->tol = SlepcDefaultTol(eps->tol);
  dvd->eps = eps;

  /* Setup the extraction technique */
  if (!eps->extraction) {
    if (ipB || ispositive) eps->extraction = EPS_RITZ;
    else {
      switch (eps->which) {
        case EPS_TARGET_REAL:
        case EPS_TARGET_MAGNITUDE:
        case EPS_TARGET_IMAGINARY:
        case EPS_SMALLEST_MAGNITUDE:
        case EPS_SMALLEST_REAL:
        case EPS_SMALLEST_IMAGINARY:
          eps->extraction = EPS_HARMONIC;
          break;
        case EPS_LARGEST_REAL:
        case EPS_LARGEST_MAGNITUDE:
        case EPS_LARGEST_IMAGINARY:
          eps->extraction = EPS_HARMONIC_LARGEST;
          break;
        default:
          eps->extraction = EPS_RITZ;
      }
    }
  }
  switch (eps->extraction) {
    case EPS_RITZ:              harm = DVD_HARM_NONE; break;
    case EPS_HARMONIC:          harm = DVD_HARM_RR; break;
    case EPS_HARMONIC_RELATIVE: harm = DVD_HARM_RRR; break;
    case EPS_HARMONIC_RIGHT:    harm = DVD_HARM_REIGS; break;
    case EPS_HARMONIC_LARGEST:  harm = DVD_HARM_LEIGS; break;
    default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");
  }

  /* Setup the type of starting subspace */
  init = data->krylovstart? DVD_INITV_KRYLOV: DVD_INITV_CLASSIC;

  /* Preconfigure dvd */
  PetscCall(STGetKSP(eps->st,&ksp));
  PetscCall(dvd_schm_basic_preconf(dvd,&b,eps->mpd,min_size_V,bs,initv,PetscAbs(eps->nini),data->plusk,harm,ksp,init,eps->trackall,data->ipB,data->doubleexp));

  /* Allocate memory */
  PetscCall(EPSAllocateSolution(eps,0));

  /* Setup orthogonalization */
  PetscCall(EPS_SetInnerProduct(eps));
  if (!(ipB && dvd->B)) PetscCall(BVSetMatrix(eps->V,NULL,PETSC_FALSE));

  /* Configure dvd for a basic GD */
  PetscCall(dvd_schm_basic_conf(dvd,&b,eps->mpd,min_size_V,bs,initv,PetscAbs(eps->nini),data->plusk,harm,dvd->withTarget,target,ksp,data->fix,init,eps->trackall,data->ipB,data->dynamic,data->doubleexp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_XD(EPS eps)
{
  EPS_DAVIDSON   *data = (EPS_DAVIDSON*)eps->data;
  dvdDashboard   *d = &data->ddb;
  PetscInt       l,k;

  PetscFunctionBegin;
  PetscCall(PetscCitationsRegister(citation,&cited));
  /* Call the starting routines */
  PetscCall(EPSDavidsonFLCall(d->startList,d));

  while (eps->reason == EPS_CONVERGED_ITERATING) {

    /* Initialize V, if it is needed */
    PetscCall(BVGetActiveColumns(d->eps->V,&l,&k));
    if (PetscUnlikely(l == k)) PetscCall(d->initV(d));

    /* Find the best approximated eigenpairs in V, X */
    PetscCall(d->calcPairs(d));

    /* Test for convergence */
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx));
    if (eps->reason != EPS_CONVERGED_ITERATING) break;

    /* Expand the subspace */
    PetscCall(d->updateV(d));

    /* Monitor progress */
    eps->nconv = d->nconv;
    eps->its++;
    PetscCall(BVGetActiveColumns(d->eps->V,NULL,&k));
    PetscCall(EPSMonitor(eps,eps->its,eps->nconv+d->npreconv,eps->eigr,eps->eigi,eps->errest,PetscMin(k,eps->nev)));
  }

  /* Call the ending routines */
  PetscCall(EPSDavidsonFLCall(d->endList,d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSReset_XD(EPS eps)
{
  EPS_DAVIDSON   *data = (EPS_DAVIDSON*)eps->data;
  dvdDashboard   *dvd = &data->ddb;

  PetscFunctionBegin;
  /* Call step destructors and destroys the list */
  PetscCall(EPSDavidsonFLCall(dvd->destroyList,dvd));
  PetscCall(EPSDavidsonFLDestroy(&dvd->destroyList));
  PetscCall(EPSDavidsonFLDestroy(&dvd->startList));
  PetscCall(EPSDavidsonFLDestroy(&dvd->endList));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDSetKrylovStart_XD(EPS eps,PetscBool krylovstart)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  data->krylovstart = krylovstart;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDGetKrylovStart_XD(EPS eps,PetscBool *krylovstart)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *krylovstart = data->krylovstart;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDSetBlockSize_XD(EPS eps,PetscInt blocksize)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (blocksize == PETSC_DEFAULT || blocksize == PETSC_DECIDE) blocksize = 1;
  PetscCheck(blocksize>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid blocksize value, must be >0");
  if (data->blocksize != blocksize) {
    data->blocksize = blocksize;
    eps->state      = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDGetBlockSize_XD(EPS eps,PetscInt *blocksize)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *blocksize = data->blocksize;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDSetRestart_XD(EPS eps,PetscInt minv,PetscInt plusk)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (minv == PETSC_DETERMINE) {
    if (data->minv != 0) eps->state = EPS_STATE_INITIAL;
    data->minv = 0;
  } else if (minv != PETSC_CURRENT) {
    PetscCheck(minv>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid minv value, must be >0");
    if (data->minv != minv) eps->state = EPS_STATE_INITIAL;
    data->minv = minv;
  }
  if (plusk == PETSC_DETERMINE) {
    if (data->plusk != PETSC_DETERMINE) eps->state = EPS_STATE_INITIAL;
    data->plusk = PETSC_DETERMINE;
  } else if (plusk != PETSC_CURRENT) {
    PetscCheck(plusk>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid plusk value, must be >0");
    if (data->plusk != plusk) eps->state = EPS_STATE_INITIAL;
    data->plusk = plusk;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDGetRestart_XD(EPS eps,PetscInt *minv,PetscInt *plusk)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (minv) *minv = data->minv;
  if (plusk) *plusk = data->plusk;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDGetInitialSize_XD(EPS eps,PetscInt *initialsize)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *initialsize = data->initialsize;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDSetInitialSize_XD(EPS eps,PetscInt initialsize)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (initialsize == PETSC_DEFAULT || initialsize == PETSC_DECIDE) initialsize = 0;
  else PetscCheck(initialsize>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid initial size value, must be >0");
  if (data->initialsize != initialsize) {
    data->initialsize = initialsize;
    eps->state        = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDSetBOrth_XD(EPS eps,PetscBool borth)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  data->ipB = borth;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSXDGetBOrth_XD(EPS eps,PetscBool *borth)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *borth = data->ipB;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSComputeVectors_XD - Compute eigenvectors from the vectors
  provided by the eigensolver. This version is intended for solvers
  that provide Schur vectors from the QZ decomposition. Given the partial
  Schur decomposition OP*V=V*T, the following steps are performed:
      1) compute eigenvectors of (S,T): S*Z=T*Z*D
      2) compute eigenvectors of OP: X=V*Z
 */
PetscErrorCode EPSComputeVectors_XD(EPS eps)
{
  Mat            X;
  PetscBool      symm;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->ds,DSHEP,&symm));
  if (symm) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));

  /* V <- V * X */
  PetscCall(DSGetMat(eps->ds,DS_MAT_X,&X));
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
  PetscCall(BVMultInPlace(eps->V,X,0,eps->nconv));
  PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&X));
  PetscFunctionReturn(PETSC_SUCCESS);
}
