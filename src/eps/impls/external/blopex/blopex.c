/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the BLOPEX package
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include "blopex.h"
#include <lobpcg.h>
#include <interpreter.h>
#include <multivector.h>
#include <temp_multivector.h>

PetscInt slepc_blopex_useconstr = -1;

typedef struct {
  lobpcg_Tolerance           tol;
  lobpcg_BLASLAPACKFunctions blap_fn;
  mv_InterfaceInterpreter    ii;
  ST                         st;
  Vec                        w;
  PetscInt                   bs;     /* block size */
} EPS_BLOPEX;

static void Precond_FnSingleVector(void *data,void *x,void *y)
{
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)data;
  MPI_Comm       comm = PetscObjectComm((PetscObject)blopex->st);
  KSP            ksp;

  PetscFunctionBegin;
  PetscCallAbort(comm,STGetKSP(blopex->st,&ksp));
  PetscCallAbort(comm,KSPSolve(ksp,(Vec)x,(Vec)y));
  PetscFunctionReturnVoid();
}

static void Precond_FnMultiVector(void *data,void *x,void *y)
{
  EPS_BLOPEX *blopex = (EPS_BLOPEX*)data;

  PetscFunctionBegin;
  blopex->ii.Eval(Precond_FnSingleVector,data,x,y);
  PetscFunctionReturnVoid();
}

static void OperatorASingleVector(void *data,void *x,void *y)
{
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)data;
  MPI_Comm       comm = PetscObjectComm((PetscObject)blopex->st);
  Mat            A,B;
  PetscScalar    sigma;
  PetscInt       nmat;

  PetscFunctionBegin;
  PetscCallAbort(comm,STGetNumMatrices(blopex->st,&nmat));
  PetscCallAbort(comm,STGetMatrix(blopex->st,0,&A));
  if (nmat>1) PetscCallAbort(comm,STGetMatrix(blopex->st,1,&B));
  PetscCallAbort(comm,MatMult(A,(Vec)x,(Vec)y));
  PetscCallAbort(comm,STGetShift(blopex->st,&sigma));
  if (sigma != 0.0) {
    if (nmat>1) PetscCallAbort(comm,MatMult(B,(Vec)x,blopex->w));
    else PetscCallAbort(comm,VecCopy((Vec)x,blopex->w));
    PetscCallAbort(comm,VecAXPY((Vec)y,-sigma,blopex->w));
  }
  PetscFunctionReturnVoid();
}

static void OperatorAMultiVector(void *data,void *x,void *y)
{
  EPS_BLOPEX *blopex = (EPS_BLOPEX*)data;

  PetscFunctionBegin;
  blopex->ii.Eval(OperatorASingleVector,data,x,y);
  PetscFunctionReturnVoid();
}

static void OperatorBSingleVector(void *data,void *x,void *y)
{
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)data;
  MPI_Comm       comm = PetscObjectComm((PetscObject)blopex->st);
  Mat            B;

  PetscFunctionBegin;
  PetscCallAbort(comm,STGetMatrix(blopex->st,1,&B));
  PetscCallAbort(comm,MatMult(B,(Vec)x,(Vec)y));
  PetscFunctionReturnVoid();
}

static void OperatorBMultiVector(void *data,void *x,void *y)
{
  EPS_BLOPEX *blopex = (EPS_BLOPEX*)data;

  PetscFunctionBegin;
  blopex->ii.Eval(OperatorBSingleVector,data,x,y);
  PetscFunctionReturnVoid();
}

static PetscErrorCode EPSSetDimensions_BLOPEX(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  EPS_BLOPEX     *ctx = (EPS_BLOPEX*)eps->data;
  PetscInt       k;

  PetscFunctionBegin;
  if (*nev==0) *nev = 1;
  k = ((*nev-1)/ctx->bs+1)*ctx->bs;
  if (*ncv!=PETSC_DETERMINE) { /* ncv set */
    PetscCheck(*ncv>=k,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv is not sufficiently large");
  } else *ncv = k;
  if (*mpd==PETSC_DETERMINE) *mpd = *ncv;
  else PetscCall(PetscInfo(eps,"Warning: given value of mpd ignored\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUp_BLOPEX(EPS eps)
{
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)eps->data;
  PetscBool      flg;
  KSP            ksp;

  PetscFunctionBegin;
  EPSCheckHermitianDefinite(eps);
  EPSCheckNotStructured(eps);
  if (!blopex->bs) blopex->bs = PetscMin(16,eps->nev?eps->nev:1);
  PetscCall(EPSSetDimensions_BLOPEX(eps,&eps->nev,&eps->ncv,&eps->mpd));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) eps->which = EPS_SMALLEST_REAL;
  PetscCheck(eps->which==EPS_SMALLEST_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver supports only smallest real eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_STOPPING);
  EPSCheckIgnored(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_EXTRACTION);

  blopex->st = eps->st;

  PetscCheck(eps->converged==EPSConvergedRelative || eps->converged==EPSConvergedAbsolute,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Convergence test not supported in this solver");
  if (eps->converged == EPSConvergedRelative) {
    blopex->tol.absolute = 0.0;
    blopex->tol.relative = SlepcDefaultTol(eps->tol);
  } else {  /* EPSConvergedAbsolute */
    blopex->tol.absolute = SlepcDefaultTol(eps->tol);
    blopex->tol.relative = 0.0;
  }

  SLEPCSetupInterpreter(&blopex->ii);

  PetscCall(STGetKSP(eps->st,&ksp));
  PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&flg));
  if (!flg) PetscCall(PetscInfo(eps,"Warning: ignoring KSP, should use KSPPREONLY\n"));

  /* allocate memory */
  if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
  PetscCall(PetscObjectTypeCompareAny((PetscObject)eps->V,&flg,BVVECS,BVCONTIGUOUS,""));
  if (!flg) {  /* blopex only works with BVVECS or BVCONTIGUOUS */
    PetscCall(BVSetType(eps->V,BVCONTIGUOUS));
  }
  PetscCall(EPSAllocateSolution(eps,0));
  if (!blopex->w) PetscCall(BVCreateVec(eps->V,&blopex->w));

#if defined(PETSC_USE_COMPLEX)
  blopex->blap_fn.zpotrf = PETSC_zpotrf_interface;
  blopex->blap_fn.zhegv = PETSC_zsygv_interface;
#else
  blopex->blap_fn.dpotrf = PETSC_dpotrf_interface;
  blopex->blap_fn.dsygv = PETSC_dsygv_interface;
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_BLOPEX(EPS eps)
{
  EPS_BLOPEX        *blopex = (EPS_BLOPEX*)eps->data;
  PetscScalar       sigma,*eigr=NULL;
  PetscReal         *errest=NULL;
  int               i,j,info,its,nconv;
  double            *residhist=NULL;
  mv_MultiVectorPtr eigenvectors,constraints;
#if defined(PETSC_USE_COMPLEX)
  komplex           *lambda=NULL,*lambdahist=NULL;
#else
  double            *lambda=NULL,*lambdahist=NULL;
#endif

  PetscFunctionBegin;
  PetscCall(STGetShift(eps->st,&sigma));
  PetscCall(PetscMalloc1(blopex->bs,&lambda));
  if (eps->numbermonitors>0) PetscCall(PetscMalloc4(blopex->bs*(eps->max_it+1),&lambdahist,eps->ncv,&eigr,blopex->bs*(eps->max_it+1),&residhist,eps->ncv,&errest));

  /* Complete the initial basis with random vectors */
  for (i=0;i<eps->nini;i++) {  /* in case the initial vectors were also set with VecSetRandom */
    PetscCall(BVSetRandomColumn(eps->V,eps->nini));
  }
  for (i=eps->nini;i<eps->ncv;i++) PetscCall(BVSetRandomColumn(eps->V,i));

  while (eps->reason == EPS_CONVERGED_ITERATING) {

    /* Create multivector of constraints from leading columns of V */
    PetscCall(PetscObjectComposedDataSetInt((PetscObject)eps->V,slepc_blopex_useconstr,1));
    PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
    constraints = mv_MultiVectorCreateFromSampleVector(&blopex->ii,eps->nds+eps->nconv,eps->V);

    /* Create multivector where eigenvectors of this run will be stored */
    PetscCall(PetscObjectComposedDataSetInt((PetscObject)eps->V,slepc_blopex_useconstr,0));
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,eps->nconv+blopex->bs));
    eigenvectors = mv_MultiVectorCreateFromSampleVector(&blopex->ii,blopex->bs,eps->V);

#if defined(PETSC_USE_COMPLEX)
    info = lobpcg_solve_complex(eigenvectors,blopex,OperatorAMultiVector,
          eps->isgeneralized?blopex:NULL,eps->isgeneralized?OperatorBMultiVector:NULL,
          blopex,Precond_FnMultiVector,constraints,
          blopex->blap_fn,blopex->tol,eps->max_it,0,&its,
          lambda,lambdahist,blopex->bs,eps->errest+eps->nconv,residhist,blopex->bs);
#else
    info = lobpcg_solve_double(eigenvectors,blopex,OperatorAMultiVector,
          eps->isgeneralized?blopex:NULL,eps->isgeneralized?OperatorBMultiVector:NULL,
          blopex,Precond_FnMultiVector,constraints,
          blopex->blap_fn,blopex->tol,eps->max_it,0,&its,
          lambda,lambdahist,blopex->bs,eps->errest+eps->nconv,residhist,blopex->bs);
#endif
    PetscCheck(info==0,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"BLOPEX failed with exit code=%d",info);
    mv_MultiVectorDestroy(constraints);
    mv_MultiVectorDestroy(eigenvectors);

    for (j=0;j<blopex->bs;j++) {
#if defined(PETSC_USE_COMPLEX)
      eps->eigr[eps->nconv+j] = PetscCMPLX(lambda[j].real,lambda[j].imag);
#else
      eps->eigr[eps->nconv+j] = lambda[j];
#endif
    }

    if (eps->numbermonitors>0) {
      for (i=0;i<its;i++) {
        nconv = 0;
        for (j=0;j<blopex->bs;j++) {
#if defined(PETSC_USE_COMPLEX)
          eigr[eps->nconv+j] = PetscCMPLX(lambdahist[j+i*blopex->bs].real,lambdahist[j+i*blopex->bs].imag);
#else
          eigr[eps->nconv+j] = lambdahist[j+i*blopex->bs];
#endif
          errest[eps->nconv+j] = residhist[j+i*blopex->bs];
          if (residhist[j+i*blopex->bs]<=eps->tol) nconv++;
        }
        PetscCall(EPSMonitor(eps,eps->its+i,eps->nconv+nconv,eigr,eps->eigi,errest,eps->nconv+blopex->bs));
      }
    }

    eps->its += its;
    if (info==-1) {
      eps->reason = EPS_DIVERGED_ITS;
      break;
    } else {
      for (i=0;i<blopex->bs;i++) {
        if (sigma != 0.0) eps->eigr[eps->nconv+i] += sigma;
      }
      eps->nconv += blopex->bs;
      if (eps->nconv>=eps->nev) eps->reason = EPS_CONVERGED_TOL;
    }
  }

  PetscCall(PetscFree(lambda));
  if (eps->numbermonitors>0) PetscCall(PetscFree4(lambdahist,eigr,residhist,errest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSBLOPEXSetBlockSize_BLOPEX(EPS eps,PetscInt bs)
{
  EPS_BLOPEX *ctx = (EPS_BLOPEX*)eps->data;

  PetscFunctionBegin;
  if (bs==PETSC_DEFAULT || bs==PETSC_DECIDE) {
    ctx->bs    = 0;
    eps->state = EPS_STATE_INITIAL;
  } else {
    PetscCheck(bs>0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Block size must be >0");
    ctx->bs = bs;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSBLOPEXSetBlockSize - Sets the block size of the BLOPEX solver.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  bs  - the block size

   Options Database Key:
.  -eps_blopex_blocksize - Sets the block size

   Level: advanced

.seealso: EPSBLOPEXGetBlockSize()
@*/
PetscErrorCode EPSBLOPEXSetBlockSize(EPS eps,PetscInt bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,bs,2);
  PetscTryMethod(eps,"EPSBLOPEXSetBlockSize_C",(EPS,PetscInt),(eps,bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSBLOPEXGetBlockSize_BLOPEX(EPS eps,PetscInt *bs)
{
  EPS_BLOPEX *ctx = (EPS_BLOPEX*)eps->data;

  PetscFunctionBegin;
  *bs = ctx->bs;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSBLOPEXGetBlockSize - Gets the block size used in the BLOPEX solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  bs - the block size

   Level: advanced

.seealso: EPSBLOPEXSetBlockSize()
@*/
PetscErrorCode EPSBLOPEXGetBlockSize(EPS eps,PetscInt *bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(bs,2);
  PetscUseMethod(eps,"EPSBLOPEXGetBlockSize_C",(EPS,PetscInt*),(eps,bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_BLOPEX(EPS eps)
{
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)eps->data;

  PetscFunctionBegin;
  PetscCall(VecDestroy(&blopex->w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_BLOPEX(EPS eps)
{
  PetscFunctionBegin;
  LOBPCG_DestroyRandomContext();
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSBLOPEXSetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSBLOPEXGetBlockSize_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_BLOPEX(EPS eps,PetscViewer viewer)
{
  EPS_BLOPEX     *ctx = (EPS_BLOPEX*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) PetscCall(PetscViewerASCIIPrintf(viewer,"  block size %" PetscInt_FMT "\n",ctx->bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_BLOPEX(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  PetscBool      flg;
  PetscInt       bs;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS BLOPEX Options");

    PetscCall(PetscOptionsInt("-eps_blopex_blocksize","Block size","EPSBLOPEXSetBlockSize",20,&bs,&flg));
    if (flg) PetscCall(EPSBLOPEXSetBlockSize(eps,bs));

  PetscOptionsHeadEnd();

  LOBPCG_SetFromOptionsRandomContext();
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_BLOPEX(EPS eps)
{
  EPS_BLOPEX     *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_BLOPEX;
  eps->ops->setup          = EPSSetUp_BLOPEX;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->setfromoptions = EPSSetFromOptions_BLOPEX;
  eps->ops->destroy        = EPSDestroy_BLOPEX;
  eps->ops->reset          = EPSReset_BLOPEX;
  eps->ops->view           = EPSView_BLOPEX;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_GMRES;

  LOBPCG_InitRandomContext(PetscObjectComm((PetscObject)eps),NULL);
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSBLOPEXSetBlockSize_C",EPSBLOPEXSetBlockSize_BLOPEX));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSBLOPEXGetBlockSize_C",EPSBLOPEXGetBlockSize_BLOPEX));
  if (slepc_blopex_useconstr < 0) PetscCall(PetscObjectComposedDataRegister(&slepc_blopex_useconstr));
  PetscFunctionReturn(PETSC_SUCCESS);
}
