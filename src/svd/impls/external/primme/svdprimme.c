/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the PRIMME SVD solver
*/

#include <slepc/private/svdimpl.h>    /*I "slepcsvd.h" I*/

#include <primme.h>

#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#define PRIMME_DRIVER cprimme_svds
#else
#define PRIMME_DRIVER zprimme_svds
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
#define PRIMME_DRIVER sprimme_svds
#else
#define PRIMME_DRIVER dprimme_svds
#endif
#endif

#if defined(PRIMME_VERSION_MAJOR) && PRIMME_VERSION_MAJOR*100+PRIMME_VERSION_MINOR >= 202
#define SLEPC_HAVE_PRIMME2p2
#endif

typedef struct {
  primme_svds_params        primme;   /* param struct */
  PetscInt                  bs;       /* block size */
  primme_svds_preset_method method;   /* primme method */
  SVD                       svd;      /* reference to the solver */
  Vec                       x,y;      /* auxiliary vectors */
} SVD_PRIMME;

static void multMatvec_PRIMME(void*,PRIMME_INT*,void*,PRIMME_INT*,int*,int*,struct primme_svds_params*,int*);

static void par_GlobalSumReal(void *sendBuf,void *recvBuf,int *count,primme_svds_params *primme,int *ierr)
{
  if (sendBuf == recvBuf) {
    *ierr = MPI_Allreduce(MPI_IN_PLACE,recvBuf,*count,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)primme->commInfo));
  } else {
    *ierr = MPI_Allreduce(sendBuf,recvBuf,*count,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)primme->commInfo));
  }
}

#if defined(SLEPC_HAVE_PRIMME3)
static void par_broadcastReal(void *buf,int *count,primme_svds_params *primme,int *ierr)
{
  *ierr = MPI_Bcast(buf,*count,MPIU_REAL,0/*root*/,PetscObjectComm((PetscObject)primme->commInfo));
}
#endif

#if defined(SLEPC_HAVE_PRIMME2p2)
static void convTestFun(double *sval,void *leftsvec,void *rightsvec,double *resNorm,
#if defined(SLEPC_HAVE_PRIMME3)
                        int *method,
#endif
                        int *isconv,struct primme_svds_params *primme,int *err)
{
  PetscErrorCode ierr;
  SVD            svd = (SVD)primme->commInfo;
  PetscReal      sigma = sval?*sval:0.0;
  PetscReal      r = resNorm?*resNorm:HUGE_VAL,errest;

  ierr = (*svd->converged)(svd,sigma,r,&errest,svd->convergedctx);
  if (ierr) *err = 1;
  else {
    *isconv = (errest<=svd->tol?1:0);
    *err = 0;
  }
}

static void monitorFun(void *basisSvals,int *basisSize,int *basisFlags,int *iblock,int *blockSize,void *basisNorms,int *numConverged,void *lockedSvals,int *numLocked,int *lockedFlags,void *lockedNorms,int *inner_its,void *LSRes,
#if defined(SLEPC_HAVE_PRIMME3)
                       const char *msg,double *time,
#endif
                       primme_event *event,int *stage,struct primme_svds_params *primme,int *err)
{

  PetscErrorCode ierr = PETSC_SUCCESS;
  SVD            svd = (SVD)primme->commInfo;
  PetscInt       i,k,nerrest;

  *err = 1;
  switch (*event) {
    case primme_event_outer_iteration:
      /* Update SVD */
      PetscCallVoid(PetscIntCast(primme->stats.numOuterIterations,&svd->its));
      if (numConverged) svd->nconv = *numConverged;
      k = 0;
      if (lockedSvals && numLocked) for (i=0; i<*numLocked && k<svd->ncv; i++) svd->sigma[k++] = ((PetscReal*)lockedSvals)[i];
      nerrest = k;
      if (iblock && blockSize) {
        for (i=0; i<*blockSize && k+iblock[i]<svd->ncv; i++) svd->errest[k+iblock[i]] = ((PetscReal*)basisNorms)[i];
        nerrest = k+(*blockSize>0?1+iblock[*blockSize-1]:0);
      }
      if (basisSvals && basisSize) for (i=0; i<*basisSize && k<svd->ncv; i++) svd->sigma[k++] = ((PetscReal*)basisSvals)[i];
      /* Show progress */
      ierr = SVDMonitor(svd,svd->its,numConverged?*numConverged:0,svd->sigma,svd->errest,nerrest);
      break;
#if defined(SLEPC_HAVE_PRIMME3)
    case primme_event_message:
      /* Print PRIMME information messages */
      ierr = PetscInfo(svd,"%s\n",msg);
      break;
#endif
    default:
      break;
  }
  *err = (ierr!=0)? 1: 0;
}
#endif /* SLEPC_HAVE_PRIMME2p2 */

static void multMatvec_PRIMME(void *xa,PRIMME_INT *ldx,void *ya,PRIMME_INT *ldy,int *blockSize,int *transpose,struct primme_svds_params *primme,int *ierr)
{
  PetscInt   i;
  SVD_PRIMME *ops = (SVD_PRIMME*)primme->matrix;
  Vec        x = ops->x,y = ops->y;
  SVD        svd = ops->svd;

  PetscFunctionBegin;
  for (i=0;i<*blockSize;i++) {
    if (*transpose) {
      PetscCallAbort(PetscObjectComm((PetscObject)svd),VecPlaceArray(y,(PetscScalar*)xa+(*ldx)*i));
      PetscCallAbort(PetscObjectComm((PetscObject)svd),VecPlaceArray(x,(PetscScalar*)ya+(*ldy)*i));
      PetscCallAbort(PetscObjectComm((PetscObject)svd),MatMult(svd->AT,y,x));
    } else {
      PetscCallAbort(PetscObjectComm((PetscObject)svd),VecPlaceArray(x,(PetscScalar*)xa+(*ldx)*i));
      PetscCallAbort(PetscObjectComm((PetscObject)svd),VecPlaceArray(y,(PetscScalar*)ya+(*ldy)*i));
      PetscCallAbort(PetscObjectComm((PetscObject)svd),MatMult(svd->A,x,y));
    }
    PetscCallAbort(PetscObjectComm((PetscObject)svd),VecResetArray(x));
    PetscCallAbort(PetscObjectComm((PetscObject)svd),VecResetArray(y));
  }
  PetscFunctionReturnVoid();
}

static PetscErrorCode SVDSetUp_PRIMME(SVD svd)
{
  PetscMPIInt        numProcs,procID;
  PetscInt           n,m,nloc,mloc;
  SVD_PRIMME         *ops = (SVD_PRIMME*)svd->data;
  primme_svds_params *primme = &ops->primme;

  PetscFunctionBegin;
  SVDCheckStandard(svd);
  SVDCheckDefinite(svd);
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)svd),&numProcs));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)svd),&procID));

  /* Check some constraints and set some default values */
  PetscCall(MatGetSize(svd->A,&m,&n));
  PetscCall(MatGetLocalSize(svd->A,&mloc,&nloc));
  PetscCall(SVDSetDimensions_Default(svd));
  if (svd->max_it==PETSC_DETERMINE) svd->max_it = PETSC_INT_MAX;
  svd->leftbasis = PETSC_TRUE;
  SVDCheckUnsupported(svd,SVD_FEATURE_STOPPING);
#if !defined(SLEPC_HAVE_PRIMME2p2)
  if (svd->converged != SVDConvergedAbsolute) PetscCall(PetscInfo(svd,"Warning: using absolute convergence test\n"));
#endif

  /* Transfer SLEPc options to PRIMME options */
  primme_svds_free(primme);
  primme_svds_initialize(primme);
  primme->m             = (PRIMME_INT)m;
  primme->n             = (PRIMME_INT)n;
  primme->mLocal        = (PRIMME_INT)mloc;
  primme->nLocal        = (PRIMME_INT)nloc;
  primme->numSvals      = (int)svd->nsv;
  primme->matrix        = ops;
  primme->commInfo      = svd;
  primme->maxMatvecs    = (PRIMME_INT)svd->max_it;
#if !defined(SLEPC_HAVE_PRIMME2p2)
  primme->eps           = SlepcDefaultTol(svd->tol);
#endif
  primme->numProcs      = numProcs;
  primme->procID        = procID;
  primme->printLevel    = 1;
  primme->matrixMatvec  = multMatvec_PRIMME;
  primme->globalSumReal = par_GlobalSumReal;
#if defined(SLEPC_HAVE_PRIMME3)
  primme->broadcastReal = par_broadcastReal;
#endif
#if defined(SLEPC_HAVE_PRIMME2p2)
  primme->convTestFun   = convTestFun;
  primme->monitorFun    = monitorFun;
#endif
  if (ops->bs > 0) primme->maxBlockSize = (int)ops->bs;

  switch (svd->which) {
    case SVD_LARGEST:
      primme->target = primme_svds_largest;
      break;
    case SVD_SMALLEST:
      primme->target = primme_svds_smallest;
      break;
  }

  /* If user sets mpd or ncv, maxBasisSize is modified */
  if (svd->mpd!=PETSC_DETERMINE) {
    primme->maxBasisSize = (int)svd->mpd;
    if (svd->ncv!=PETSC_DETERMINE) PetscCall(PetscInfo(svd,"Warning: 'ncv' is ignored by PRIMME\n"));
  } else if (svd->ncv!=PETSC_DETERMINE) primme->maxBasisSize = (int)svd->ncv;

  PetscCheck(primme_svds_set_method(ops->method,(primme_preset_method)EPS_PRIMME_DEFAULT_MIN_TIME,PRIMME_DEFAULT_METHOD,primme)>=0,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"PRIMME method not valid");

  svd->mpd = (PetscInt)primme->maxBasisSize;
  svd->ncv = (PetscInt)(primme->locking?svd->nsv:0)+primme->maxBasisSize;
  ops->bs  = (PetscInt)primme->maxBlockSize;

  /* Set workspace */
  PetscCall(SVDAllocateSolution(svd,0));

  /* Prepare auxiliary vectors */
  if (!ops->x) PetscCall(MatCreateVecsEmpty(svd->A,&ops->x,&ops->y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_PRIMME(SVD svd)
{
  SVD_PRIMME     *ops = (SVD_PRIMME*)svd->data;
  PetscScalar    *svecs, *a;
  PetscInt       i,ierrprimme,ld;
  PetscReal      *svals,*rnorms;

  PetscFunctionBegin;
  /* Reset some parameters left from previous runs */
  ops->primme.aNorm    = 0.0;
  ops->primme.initSize = (int)svd->nini;
  ops->primme.iseed[0] = -1;
  ops->primme.iseed[1] = -1;
  ops->primme.iseed[2] = -1;
  ops->primme.iseed[3] = -1;

  /* Allocating left and right singular vectors contiguously */
  PetscCall(PetscCalloc1(ops->primme.numSvals*(ops->primme.mLocal+ops->primme.nLocal),&svecs));

  /* Call PRIMME solver */
  PetscCall(PetscMalloc2(svd->ncv,&svals,svd->ncv,&rnorms));
  ierrprimme = PRIMME_DRIVER(svals,svecs,rnorms,&ops->primme);
  for (i=0;i<svd->ncv;i++) svd->sigma[i] = svals[i];
  for (i=0;i<svd->ncv;i++) svd->errest[i] = rnorms[i];
  PetscCall(PetscFree2(svals,rnorms));

  /* Copy left and right singular vectors into svd */
  PetscCall(BVGetLeadingDimension(svd->U,&ld));
  PetscCall(BVGetArray(svd->U,&a));
  for (i=0;i<ops->primme.initSize;i++) PetscCall(PetscArraycpy(a+i*ld,svecs+i*ops->primme.mLocal,ops->primme.mLocal));
  PetscCall(BVRestoreArray(svd->U,&a));

  PetscCall(BVGetLeadingDimension(svd->V,&ld));
  PetscCall(BVGetArray(svd->V,&a));
  for (i=0;i<ops->primme.initSize;i++) PetscCall(PetscArraycpy(a+i*ld,svecs+ops->primme.mLocal*ops->primme.initSize+i*ops->primme.nLocal,ops->primme.nLocal));
  PetscCall(BVRestoreArray(svd->V,&a));

  PetscCall(PetscFree(svecs));

  svd->nconv  = ops->primme.initSize >= 0 ? (PetscInt)ops->primme.initSize : 0;
  svd->reason = svd->nconv >= svd->nsv ? SVD_CONVERGED_TOL: SVD_DIVERGED_ITS;
  PetscCall(PetscIntCast(ops->primme.stats.numOuterIterations,&svd->its));

  /* Process PRIMME error code */
  if (ierrprimme != 0) {
    switch (ierrprimme%100) {
      case -1:
        SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"PRIMME library failed with error code=%" PetscInt_FMT ": unexpected error",ierrprimme);
      case -2:
        SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"PRIMME library failed with error code=%" PetscInt_FMT ": allocation error",ierrprimme);
      case -3: /* stop due to maximum number of iterations or matvecs */
        break;
      default:
        PetscCheck(ierrprimme<-39,PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"PRIMME library failed with error code=%" PetscInt_FMT ": configuration error; check PRIMME's manual",ierrprimme);
        PetscCheck(ierrprimme>=-39,PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"PRIMME library failed with error code=%" PetscInt_FMT ": runtime error; check PRIMME's manual",ierrprimme);
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDReset_PRIMME(SVD svd)
{
  SVD_PRIMME     *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  primme_svds_free(&ops->primme);
  PetscCall(VecDestroy(&ops->x));
  PetscCall(VecDestroy(&ops->y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDDestroy_PRIMME(SVD svd)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(svd->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMESetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMEGetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMESetMethod_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMEGetMethod_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDView_PRIMME(SVD svd,PetscViewer viewer)
{
  PetscBool      isascii;
  SVD_PRIMME     *ctx = (SVD_PRIMME*)svd->data;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  block size=%" PetscInt_FMT "\n",ctx->bs));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  solver method: %s\n",SVDPRIMMEMethods[(SVDPRIMMEMethod)ctx->method]));

    /* Display PRIMME params */
    PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)svd),&rank));
    if (!rank) primme_svds_display_params(ctx->primme);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSetFromOptions_PRIMME(SVD svd,PetscOptionItems PetscOptionsObject)
{
  SVD_PRIMME      *ctx = (SVD_PRIMME*)svd->data;
  PetscInt        bs;
  SVDPRIMMEMethod meth;
  PetscBool       flg;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"SVD PRIMME Options");

    PetscCall(PetscOptionsInt("-svd_primme_blocksize","Maximum block size","SVDPRIMMESetBlockSize",ctx->bs,&bs,&flg));
    if (flg) PetscCall(SVDPRIMMESetBlockSize(svd,bs));

    PetscCall(PetscOptionsEnum("-svd_primme_method","Method for solving the singular value problem","SVDPRIMMESetMethod",SVDPRIMMEMethods,(PetscEnum)ctx->method,(PetscEnum*)&meth,&flg));
    if (flg) PetscCall(SVDPRIMMESetMethod(svd,meth));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDPRIMMESetBlockSize_PRIMME(SVD svd,PetscInt bs)
{
  SVD_PRIMME *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  if (bs == PETSC_DEFAULT || bs == PETSC_DECIDE) ops->bs = 0;
  else {
    PetscCheck(bs>0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"PRIMME: block size must be positive");
    ops->bs = bs;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDPRIMMESetBlockSize - The maximum block size that PRIMME will try to use.

   Logically Collective

   Input Parameters:
+  svd - the singular value solver context
-  bs - block size

   Options Database Key:
.  -svd_primme_blocksize - Sets the max allowed block size value

   Notes:
   If the block size is not set, the value established by primme_svds_initialize
   is used.

   The user should set the block size based on the architecture specifics
   of the target computer, as well as any a priori knowledge of multiplicities.
   The code does NOT require bs > 1 to find multiple eigenvalues. For some
   methods, keeping bs = 1 yields the best overall performance.

   Level: advanced

.seealso: SVDPRIMMEGetBlockSize()
@*/
PetscErrorCode SVDPRIMMESetBlockSize(SVD svd,PetscInt bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,bs,2);
  PetscTryMethod(svd,"SVDPRIMMESetBlockSize_C",(SVD,PetscInt),(svd,bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDPRIMMEGetBlockSize_PRIMME(SVD svd,PetscInt *bs)
{
  SVD_PRIMME *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  *bs = ops->bs;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDPRIMMEGetBlockSize - Get the maximum block size the code will try to use.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  bs - returned block size

   Level: advanced

.seealso: SVDPRIMMESetBlockSize()
@*/
PetscErrorCode SVDPRIMMEGetBlockSize(SVD svd,PetscInt *bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(bs,2);
  PetscUseMethod(svd,"SVDPRIMMEGetBlockSize_C",(SVD,PetscInt*),(svd,bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDPRIMMESetMethod_PRIMME(SVD svd,SVDPRIMMEMethod method)
{
  SVD_PRIMME *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  ops->method = (primme_svds_preset_method)method;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDPRIMMESetMethod - Sets the method for the PRIMME SVD solver.

   Logically Collective

   Input Parameters:
+  svd - the singular value solver context
-  method - method that will be used by PRIMME

   Options Database Key:
.  -svd_primme_method - Sets the method for the PRIMME SVD solver

   Note:
   If not set, the method defaults to SVD_PRIMME_HYBRID.

   Level: advanced

.seealso: SVDPRIMMEGetMethod(), SVDPRIMMEMethod
@*/
PetscErrorCode SVDPRIMMESetMethod(SVD svd,SVDPRIMMEMethod method)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,method,2);
  PetscTryMethod(svd,"SVDPRIMMESetMethod_C",(SVD,SVDPRIMMEMethod),(svd,method));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDPRIMMEGetMethod_PRIMME(SVD svd,SVDPRIMMEMethod *method)
{
  SVD_PRIMME *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  *method = (SVDPRIMMEMethod)ops->method;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDPRIMMEGetMethod - Gets the method for the PRIMME SVD solver.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  method - method that will be used by PRIMME

   Level: advanced

.seealso: SVDPRIMMESetMethod(), SVDPRIMMEMethod
@*/
PetscErrorCode SVDPRIMMEGetMethod(SVD svd,SVDPRIMMEMethod *method)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(method,2);
  PetscUseMethod(svd,"SVDPRIMMEGetMethod_C",(SVD,SVDPRIMMEMethod*),(svd,method));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_PRIMME(SVD svd)
{
  SVD_PRIMME     *primme;

  PetscFunctionBegin;
  PetscCall(PetscNew(&primme));
  svd->data = (void*)primme;

  primme_svds_initialize(&primme->primme);
  primme->bs = 0;
  primme->method = (primme_svds_preset_method)SVD_PRIMME_HYBRID;
  primme->svd = svd;

  svd->ops->solve          = SVDSolve_PRIMME;
  svd->ops->setup          = SVDSetUp_PRIMME;
  svd->ops->setfromoptions = SVDSetFromOptions_PRIMME;
  svd->ops->destroy        = SVDDestroy_PRIMME;
  svd->ops->reset          = SVDReset_PRIMME;
  svd->ops->view           = SVDView_PRIMME;

  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMESetBlockSize_C",SVDPRIMMESetBlockSize_PRIMME));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMEGetBlockSize_C",SVDPRIMMEGetBlockSize_PRIMME));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMESetMethod_C",SVDPRIMMESetMethod_PRIMME));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMEGetMethod_C",SVDPRIMMEGetMethod_PRIMME));
  PetscFunctionReturn(PETSC_SUCCESS);
}
