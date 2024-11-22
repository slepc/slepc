/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the PRIMME package
*/

#include <slepc/private/epsimpl.h>    /*I "slepceps.h" I*/

#include <primme.h>

#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#define PRIMME_DRIVER cprimme
#else
#define PRIMME_DRIVER zprimme
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
#define PRIMME_DRIVER sprimme
#else
#define PRIMME_DRIVER dprimme
#endif
#endif

#if defined(PRIMME_VERSION_MAJOR) && PRIMME_VERSION_MAJOR*100+PRIMME_VERSION_MINOR >= 202
#define SLEPC_HAVE_PRIMME2p2
#endif

typedef struct {
  primme_params        primme;    /* param struct */
  PetscInt             bs;        /* block size */
  primme_preset_method method;    /* primme method */
  Mat                  A,B;       /* problem matrices */
  KSP                  ksp;       /* linear solver and preconditioner */
  Vec                  x,y;       /* auxiliary vectors */
  double               target;    /* a copy of eps's target */
} EPS_PRIMME;

static void par_GlobalSumReal(void *sendBuf,void *recvBuf,int *count,primme_params *primme,int *ierr)
{
  if (sendBuf == recvBuf) {
    *ierr = MPI_Allreduce(MPI_IN_PLACE,recvBuf,*count,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)primme->commInfo));
  } else {
    *ierr = MPI_Allreduce(sendBuf,recvBuf,*count,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)primme->commInfo));
  }
}

#if defined(SLEPC_HAVE_PRIMME3)
static void par_broadcastReal(void *buf,int *count,primme_params *primme,int *ierr)
{
  *ierr = MPI_Bcast(buf,*count,MPIU_REAL,0/*root*/,PetscObjectComm((PetscObject)primme->commInfo));
}
#endif

#if defined(SLEPC_HAVE_PRIMME2p2)
static void convTestFun(double *eval,void *evec,double *resNorm,int *isconv,primme_params *primme,int *err)
{
  PetscErrorCode ierr;
  EPS            eps = (EPS)primme->commInfo;
  PetscScalar    eigvr = eval?*eval:0.0;
  PetscReal      r = resNorm?*resNorm:HUGE_VAL,errest;

  ierr = (*eps->converged)(eps,eigvr,0.0,r,&errest,eps->convergedctx);
  if (ierr) *err = 1;
  else {
    *isconv = (errest<=eps->tol?1:0);
    *err = 0;
  }
}

static void monitorFun(void *basisEvals,int *basisSize,int *basisFlags,int *iblock,int *blockSize,void *basisNorms,int *numConverged,void *lockedEvals,int *numLocked,int *lockedFlags,void *lockedNorms,int *inner_its,void *LSRes,
#if defined(SLEPC_HAVE_PRIMME3)
                       const char *msg,double *time,
#endif
                       primme_event *event,struct primme_params *primme,int *err)
{
  PetscErrorCode ierr = PETSC_SUCCESS;
  EPS            eps = (EPS)primme->commInfo;
  PetscInt       i,k,nerrest;

  switch (*event) {
    case primme_event_outer_iteration:
      /* Update EPS */
      PetscCallVoid(PetscIntCast(primme->stats.numOuterIterations,&eps->its));
      eps->nconv = primme->initSize;
      k=0;
      if (lockedEvals && numLocked) for (i=0; i<*numLocked && k<eps->ncv; i++) eps->eigr[k++] = ((PetscReal*)lockedEvals)[i];
      nerrest = k;
      if (iblock && blockSize) {
        for (i=0; i<*blockSize && k+iblock[i]<eps->ncv; i++) eps->errest[k+iblock[i]] = ((PetscReal*)basisNorms)[i];
        nerrest = k+(*blockSize>0?1+iblock[*blockSize-1]:0);
      }
      if (basisEvals && basisSize) for (i=0; i<*basisSize && k<eps->ncv; i++) eps->eigr[k++] = ((PetscReal*)basisEvals)[i];
      /* Show progress */
      ierr = EPSMonitor(eps,eps->its,numConverged?*numConverged:0,eps->eigr,eps->eigi,eps->errest,nerrest);
      break;
#if defined(SLEPC_HAVE_PRIMME3)
    case primme_event_message:
      /* Print PRIMME information messages */
      ierr = PetscInfo(eps,"%s\n",msg);
      break;
#endif
    default:
      break;
  }
  *err = (ierr!=0)? 1: 0;
}
#endif /* SLEPC_HAVE_PRIMME2p2 */

static void matrixMatvec_PRIMME(void *xa,PRIMME_INT *ldx,void *ya,PRIMME_INT *ldy,int *blockSize,struct primme_params *primme,int *ierr)
{
  PetscInt   i;
  EPS_PRIMME *ops = (EPS_PRIMME*)primme->matrix;
  Vec        x = ops->x,y = ops->y;
  Mat        A = ops->A;

  PetscFunctionBegin;
  for (i=0;i<*blockSize;i++) {
    PetscCallAbort(PetscObjectComm((PetscObject)A),VecPlaceArray(x,(PetscScalar*)xa+(*ldx)*i));
    PetscCallAbort(PetscObjectComm((PetscObject)A),VecPlaceArray(y,(PetscScalar*)ya+(*ldy)*i));
    PetscCallAbort(PetscObjectComm((PetscObject)A),MatMult(A,x,y));
    PetscCallAbort(PetscObjectComm((PetscObject)A),VecResetArray(x));
    PetscCallAbort(PetscObjectComm((PetscObject)A),VecResetArray(y));
  }
  PetscFunctionReturnVoid();
}

#if defined(SLEPC_HAVE_PRIMME3)
static void massMatrixMatvec_PRIMME(void *xa,PRIMME_INT *ldx,void *ya,PRIMME_INT *ldy,int *blockSize,struct primme_params *primme,int *ierr)
{
  PetscInt   i;
  EPS_PRIMME *ops = (EPS_PRIMME*)primme->massMatrix;
  Vec        x = ops->x,y = ops->y;
  Mat        B = ops->B;

  PetscFunctionBegin;
  for (i=0;i<*blockSize;i++) {
    PetscCallAbort(PetscObjectComm((PetscObject)B),VecPlaceArray(x,(PetscScalar*)xa+(*ldx)*i));
    PetscCallAbort(PetscObjectComm((PetscObject)B),VecPlaceArray(y,(PetscScalar*)ya+(*ldy)*i));
    PetscCallAbort(PetscObjectComm((PetscObject)B),MatMult(B,x,y));
    PetscCallAbort(PetscObjectComm((PetscObject)B),VecResetArray(x));
    PetscCallAbort(PetscObjectComm((PetscObject)B),VecResetArray(y));
  }
  PetscFunctionReturnVoid();
}
#endif

static void applyPreconditioner_PRIMME(void *xa,PRIMME_INT *ldx,void *ya,PRIMME_INT *ldy,int *blockSize,struct primme_params *primme,int *ierr)
{
  PetscInt   i;
  EPS_PRIMME *ops = (EPS_PRIMME*)primme->matrix;
  Vec        x = ops->x,y = ops->y;

  PetscFunctionBegin;
  for (i=0;i<*blockSize;i++) {
    PetscCallAbort(PetscObjectComm((PetscObject)ops->ksp),VecPlaceArray(x,(PetscScalar*)xa+(*ldx)*i));
    PetscCallAbort(PetscObjectComm((PetscObject)ops->ksp),VecPlaceArray(y,(PetscScalar*)ya+(*ldy)*i));
    PetscCallAbort(PetscObjectComm((PetscObject)ops->ksp),KSPSolve(ops->ksp,x,y));
    PetscCallAbort(PetscObjectComm((PetscObject)ops->ksp),VecResetArray(x));
    PetscCallAbort(PetscObjectComm((PetscObject)ops->ksp),VecResetArray(y));
  }
  PetscFunctionReturnVoid();
}

static PetscErrorCode EPSSetUp_PRIMME(EPS eps)
{
  PetscMPIInt    numProcs,procID;
  EPS_PRIMME     *ops = (EPS_PRIMME*)eps->data;
  primme_params  *primme = &ops->primme;
  PetscBool      flg;

  PetscFunctionBegin;
  EPSCheckHermitianDefinite(eps);
  EPSCheckNotStructured(eps);
  if (eps->nev==0) eps->nev = 1;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)eps),&numProcs));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)eps),&procID));

  /* Check some constraints and set some default values */
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PETSC_INT_MAX;
  PetscCall(STGetMatrix(eps->st,0,&ops->A));
  if (eps->isgeneralized) {
#if defined(SLEPC_HAVE_PRIMME3)
    PetscCall(STGetMatrix(eps->st,1,&ops->B));
#else
    SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This version of PRIMME is not available for generalized problems");
#endif
  }
  EPSCheckUnsupported(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_STOPPING);
  EPSCheckIgnored(eps,EPS_FEATURE_BALANCE);
  if (!eps->which) eps->which = EPS_LARGEST_REAL;
#if !defined(SLEPC_HAVE_PRIMME2p2)
  if (eps->converged != EPSConvergedAbsolute) PetscCall(PetscInfo(eps,"Warning: using absolute convergence test\n"));
#else
  EPSCheckIgnored(eps,EPS_FEATURE_CONVERGENCE);
#endif

  /* Transfer SLEPc options to PRIMME options */
  primme_free(primme);
  primme_initialize(primme);
  primme->n                             = (PRIMME_INT)eps->n;
  primme->nLocal                        = (PRIMME_INT)eps->nloc;
  primme->numEvals                      = (int)eps->nev;
  primme->matrix                        = ops;
  primme->matrixMatvec                  = matrixMatvec_PRIMME;
#if defined(SLEPC_HAVE_PRIMME3)
  if (eps->isgeneralized) {
    primme->massMatrix                  = ops;
    primme->massMatrixMatvec            = massMatrixMatvec_PRIMME;
  }
#endif
  primme->commInfo                      = eps;
  primme->maxOuterIterations            = (PRIMME_INT)eps->max_it;
#if !defined(SLEPC_HAVE_PRIMME2p2)
  primme->eps                           = SlepcDefaultTol(eps->tol);
#endif
  primme->numProcs                      = numProcs;
  primme->procID                        = procID;
  primme->printLevel                    = 1;
  primme->correctionParams.precondition = 1;
  primme->globalSumReal                 = par_GlobalSumReal;
#if defined(SLEPC_HAVE_PRIMME3)
  primme->broadcastReal                 = par_broadcastReal;
#endif
#if defined(SLEPC_HAVE_PRIMME2p2)
  primme->convTestFun                   = convTestFun;
  primme->monitorFun                    = monitorFun;
#endif
  if (ops->bs > 0) primme->maxBlockSize = (int)ops->bs;

  switch (eps->which) {
    case EPS_LARGEST_REAL:
      primme->target = primme_largest;
      break;
    case EPS_SMALLEST_REAL:
      primme->target = primme_smallest;
      break;
    case EPS_LARGEST_MAGNITUDE:
      primme->target = primme_largest_abs;
      ops->target = 0.0;
      primme->numTargetShifts = 1;
      primme->targetShifts = &ops->target;
      break;
    case EPS_SMALLEST_MAGNITUDE:
      primme->target = primme_closest_abs;
      ops->target = 0.0;
      primme->numTargetShifts = 1;
      primme->targetShifts = &ops->target;
      break;
    case EPS_TARGET_MAGNITUDE:
    case EPS_TARGET_REAL:
      primme->target = primme_closest_abs;
      primme->numTargetShifts = 1;
      ops->target = PetscRealPart(eps->target);
      primme->targetShifts = &ops->target;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"'which' value not supported by PRIMME");
  }

  switch (eps->extraction) {
    case EPS_RITZ:
      primme->projectionParams.projection = primme_proj_RR;
      break;
    case EPS_HARMONIC:
      primme->projectionParams.projection = primme_proj_harmonic;
      break;
    case EPS_REFINED:
      primme->projectionParams.projection = primme_proj_refined;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"'extraction' value not supported by PRIMME");
  }

  /* If user sets mpd or ncv, maxBasisSize is modified */
  if (eps->mpd!=PETSC_DETERMINE) {
    primme->maxBasisSize = (int)eps->mpd;
    if (eps->ncv!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: 'ncv' is ignored by PRIMME\n"));
  } else if (eps->ncv!=PETSC_DETERMINE) primme->maxBasisSize = (int)eps->ncv;

  PetscCheck(primme_set_method(ops->method,primme)>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"PRIMME method not valid");

  eps->mpd = (PetscInt)primme->maxBasisSize;
  eps->ncv = (PetscInt)(primme->locking?eps->nev:0)+primme->maxBasisSize;
  ops->bs  = (PetscInt)primme->maxBlockSize;

  /* Set workspace */
  PetscCall(EPSAllocateSolution(eps,0));

  /* Setup the preconditioner */
  if (primme->correctionParams.precondition) {
    PetscCall(STGetKSP(eps->st,&ops->ksp));
    PetscCall(PetscObjectTypeCompare((PetscObject)ops->ksp,KSPPREONLY,&flg));
    if (!flg) PetscCall(PetscInfo(eps,"Warning: ignoring KSP, should use KSPPREONLY\n"));
    primme->preconditioner = NULL;
    primme->applyPreconditioner = applyPreconditioner_PRIMME;
  }

  /* Prepare auxiliary vectors */
  if (!ops->x) PetscCall(MatCreateVecsEmpty(ops->A,&ops->x,&ops->y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_PRIMME(EPS eps)
{
  EPS_PRIMME     *ops = (EPS_PRIMME*)eps->data;
  PetscScalar    *a;
  PetscInt       i,ld,ierrprimme;
  PetscReal      *evals,*rnorms;

  PetscFunctionBegin;
  /* Reset some parameters left from previous runs */
#if defined(SLEPC_HAVE_PRIMME2p2)
  ops->primme.aNorm    = 0.0;
#else
  /* Force PRIMME to stop by absolute error */
  ops->primme.aNorm    = 1.0;
#endif
  ops->primme.initSize = (int)eps->nini;
  ops->primme.iseed[0] = -1;
  ops->primme.iseed[1] = -1;
  ops->primme.iseed[2] = -1;
  ops->primme.iseed[3] = -1;
  PetscCall(BVGetLeadingDimension(eps->V,&ld));
  ops->primme.ldevecs  = (PRIMME_INT)ld;

  /* Call PRIMME solver */
  PetscCall(BVGetArray(eps->V,&a));
  PetscCall(PetscMalloc2(eps->ncv,&evals,eps->ncv,&rnorms));
  ierrprimme = PRIMME_DRIVER(evals,a,rnorms,&ops->primme);
  for (i=0;i<eps->ncv;i++) eps->eigr[i] = evals[i];
  for (i=0;i<eps->ncv;i++) eps->errest[i] = rnorms[i];
  PetscCall(PetscFree2(evals,rnorms));
  PetscCall(BVRestoreArray(eps->V,&a));

  eps->nconv  = ops->primme.initSize >= 0 ? (PetscInt)ops->primme.initSize : 0;
  eps->reason = eps->nconv >= eps->nev ? EPS_CONVERGED_TOL: EPS_DIVERGED_ITS;
  PetscCall(PetscIntCast(ops->primme.stats.numOuterIterations,&eps->its));

  /* Process PRIMME error code */
  switch (ierrprimme) {
    case 0: /* no error */
      break;
    case -1:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"PRIMME library failed with error code=%" PetscInt_FMT ": unexpected error",ierrprimme);
    case -2:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"PRIMME library failed with error code=%" PetscInt_FMT ": allocation error",ierrprimme);
    case -3: /* stop due to maximum number of iterations or matvecs */
      break;
    default:
      PetscCheck(ierrprimme<-39,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"PRIMME library failed with error code=%" PetscInt_FMT ": configuration error; check PRIMME's manual",ierrprimme);
      PetscCheck(ierrprimme>=-39,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"PRIMME library failed with error code=%" PetscInt_FMT ": runtime error; check PRIMME's manual",ierrprimme);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_PRIMME(EPS eps)
{
  EPS_PRIMME     *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  primme_free(&ops->primme);
  PetscCall(VecDestroy(&ops->x));
  PetscCall(VecDestroy(&ops->y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_PRIMME(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMESetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMESetMethod_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMEGetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMEGetMethod_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_PRIMME(EPS eps,PetscViewer viewer)
{
  PetscBool      isascii;
  EPS_PRIMME     *ctx = (EPS_PRIMME*)eps->data;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  block size=%" PetscInt_FMT "\n",ctx->bs));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  solver method: %s\n",EPSPRIMMEMethods[(EPSPRIMMEMethod)ctx->method]));

    /* Display PRIMME params */
    PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)eps),&rank));
    if (!rank) primme_display_params(ctx->primme);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_PRIMME(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  EPS_PRIMME      *ctx = (EPS_PRIMME*)eps->data;
  PetscInt        bs;
  EPSPRIMMEMethod meth;
  PetscBool       flg;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS PRIMME Options");

    PetscCall(PetscOptionsInt("-eps_primme_blocksize","Maximum block size","EPSPRIMMESetBlockSize",ctx->bs,&bs,&flg));
    if (flg) PetscCall(EPSPRIMMESetBlockSize(eps,bs));

    PetscCall(PetscOptionsEnum("-eps_primme_method","Method for solving the eigenproblem","EPSPRIMMESetMethod",EPSPRIMMEMethods,(PetscEnum)ctx->method,(PetscEnum*)&meth,&flg));
    if (flg) PetscCall(EPSPRIMMESetMethod(eps,meth));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPRIMMESetBlockSize_PRIMME(EPS eps,PetscInt bs)
{
  EPS_PRIMME *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  if (bs == PETSC_DEFAULT || bs == PETSC_DECIDE) ops->bs = 0;
  else {
    PetscCheck(bs>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"PRIMME: block size must be positive");
    ops->bs = bs;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPRIMMESetBlockSize - The maximum block size that PRIMME will try to use.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  bs - block size

   Options Database Key:
.  -eps_primme_blocksize - Sets the max allowed block size value

   Notes:
   If the block size is not set, the value established by primme_initialize
   is used.

   The user should set the block size based on the architecture specifics
   of the target computer, as well as any a priori knowledge of multiplicities.
   The code does NOT require bs > 1 to find multiple eigenvalues. For some
   methods, keeping bs = 1 yields the best overall performance.

   Level: advanced

.seealso: EPSPRIMMEGetBlockSize()
@*/
PetscErrorCode EPSPRIMMESetBlockSize(EPS eps,PetscInt bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,bs,2);
  PetscTryMethod(eps,"EPSPRIMMESetBlockSize_C",(EPS,PetscInt),(eps,bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPRIMMEGetBlockSize_PRIMME(EPS eps,PetscInt *bs)
{
  EPS_PRIMME *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  *bs = ops->bs;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPRIMMEGetBlockSize - Get the maximum block size the code will try to use.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  bs - returned block size

   Level: advanced

.seealso: EPSPRIMMESetBlockSize()
@*/
PetscErrorCode EPSPRIMMEGetBlockSize(EPS eps,PetscInt *bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(bs,2);
  PetscUseMethod(eps,"EPSPRIMMEGetBlockSize_C",(EPS,PetscInt*),(eps,bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPRIMMESetMethod_PRIMME(EPS eps,EPSPRIMMEMethod method)
{
  EPS_PRIMME *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  ops->method = (primme_preset_method)method;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPRIMMESetMethod - Sets the method for the PRIMME library.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  method - method that will be used by PRIMME

   Options Database Key:
.  -eps_primme_method - Sets the method for the PRIMME library

   Note:
   If not set, the method defaults to EPS_PRIMME_DEFAULT_MIN_TIME.

   Level: advanced

.seealso: EPSPRIMMEGetMethod(), EPSPRIMMEMethod
@*/
PetscErrorCode EPSPRIMMESetMethod(EPS eps,EPSPRIMMEMethod method)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,method,2);
  PetscTryMethod(eps,"EPSPRIMMESetMethod_C",(EPS,EPSPRIMMEMethod),(eps,method));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPRIMMEGetMethod_PRIMME(EPS eps,EPSPRIMMEMethod *method)
{
  EPS_PRIMME *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  *method = (EPSPRIMMEMethod)ops->method;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPRIMMEGetMethod - Gets the method for the PRIMME library.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  method - method that will be used by PRIMME

   Level: advanced

.seealso: EPSPRIMMESetMethod(), EPSPRIMMEMethod
@*/
PetscErrorCode EPSPRIMMEGetMethod(EPS eps,EPSPRIMMEMethod *method)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(method,2);
  PetscUseMethod(eps,"EPSPRIMMEGetMethod_C",(EPS,EPSPRIMMEMethod*),(eps,method));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_PRIMME(EPS eps)
{
  EPS_PRIMME     *primme;

  PetscFunctionBegin;
  PetscCall(PetscNew(&primme));
  eps->data = (void*)primme;

  primme_initialize(&primme->primme);
  primme->primme.globalSumReal = par_GlobalSumReal;
#if defined(SLEPC_HAVE_PRIMME3)
  primme->primme.broadcastReal = par_broadcastReal;
#endif
#if defined(SLEPC_HAVE_PRIMME2p2)
  primme->primme.convTestFun = convTestFun;
  primme->primme.monitorFun = monitorFun;
#endif
  primme->method = (primme_preset_method)EPS_PRIMME_DEFAULT_MIN_TIME;

  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_PRIMME;
  eps->ops->setup          = EPSSetUp_PRIMME;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->setfromoptions = EPSSetFromOptions_PRIMME;
  eps->ops->destroy        = EPSDestroy_PRIMME;
  eps->ops->reset          = EPSReset_PRIMME;
  eps->ops->view           = EPSView_PRIMME;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_GMRES;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMESetBlockSize_C",EPSPRIMMESetBlockSize_PRIMME));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMESetMethod_C",EPSPRIMMESetMethod_PRIMME));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMEGetBlockSize_C",EPSPRIMMEGetBlockSize_PRIMME));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMEGetMethod_C",EPSPRIMMEGetMethod_PRIMME));
  PetscFunctionReturn(PETSC_SUCCESS);
}
