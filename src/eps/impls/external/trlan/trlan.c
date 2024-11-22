/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the TRLAN package
*/

#include <slepc/private/epsimpl.h>
#include "trlan.h"

/* Nasty global variable to access EPS data from TRLan_ */
static struct {
  EPS eps;
  Vec x,y;
} globaldata;

static PetscErrorCode EPSSetUp_TRLAN(EPS eps)
{
  EPS_TRLAN      *tr = (EPS_TRLAN*)eps->data;

  PetscFunctionBegin;
  EPSCheckHermitian(eps);
  EPSCheckStandard(eps);
  EPSCheckNotStructured(eps);
  if (eps->nev==0) eps->nev = 1;
  PetscCall(PetscBLASIntCast(PetscMax(7,eps->nev+PetscMin(eps->nev,6)),&tr->maxlan));
  if (eps->ncv!=PETSC_DETERMINE) {
    PetscCheck(eps->ncv>=eps->nev,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev");
  } else eps->ncv = tr->maxlan;
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(1000,eps->n);

  if (!eps->which) eps->which = EPS_LARGEST_REAL;
  PetscCheck(eps->which==EPS_SMALLEST_REAL || eps->which==EPS_LARGEST_REAL || eps->which==EPS_TARGET_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver supports only smallest, largest or target real eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_STOPPING);
  EPSCheckIgnored(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_EXTRACTION);

  tr->restart = 0;
  if (tr->maxlan+1-eps->ncv<=0) PetscCall(PetscBLASIntCast(tr->maxlan*(tr->maxlan+10),&tr->lwork));
  else PetscCall(PetscBLASIntCast(eps->nloc*(tr->maxlan+1-eps->ncv) + tr->maxlan*(tr->maxlan+10),&tr->lwork));
  if (tr->work) PetscCall(PetscFree(tr->work));
  PetscCall(PetscMalloc1(tr->lwork,&tr->work));

  PetscCall(EPSAllocateSolution(eps,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscBLASInt MatMult_TRLAN(PetscBLASInt *n,PetscBLASInt *m,PetscReal *xin,PetscBLASInt *ldx,PetscReal *yout,PetscBLASInt *ldy)
{
  Vec            x=globaldata.x,y=globaldata.y;
  EPS            eps=globaldata.eps;
  PetscBLASInt   i;

  PetscFunctionBegin;
  for (i=0;i<*m;i++) {
    PetscCall(VecPlaceArray(x,(PetscScalar*)xin+i*(*ldx)));
    PetscCall(VecPlaceArray(y,(PetscScalar*)yout+i*(*ldy)));
    PetscCall(STApply(eps->st,x,y));
    PetscCall(BVOrthogonalizeVec(eps->V,y,NULL,NULL,NULL));
    PetscCall(VecResetArray(x));
    PetscCall(VecResetArray(y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_TRLAN(EPS eps)
{
  PetscInt       i;
  PetscBLASInt   ipar[32],n,lohi,stat,ncv;
  EPS_TRLAN      *tr = (EPS_TRLAN*)eps->data;
  PetscScalar    *pV;
  Vec            v0;
  Mat            A;
#if !defined(PETSC_HAVE_MPIUNI)
  MPI_Fint       fcomm;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(eps->ncv,&ncv));
  PetscCall(PetscBLASIntCast(eps->nloc,&n));

  PetscCheck(eps->which==EPS_LARGEST_REAL || eps->which==EPS_TARGET_REAL || eps->which==EPS_SMALLEST_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"Wrong value of eps->which");
  lohi = (eps->which==EPS_SMALLEST_REAL)? -1: 1;

  globaldata.eps = eps;
  PetscCall(STGetMatrix(eps->st,0,&A));
  PetscCall(MatCreateVecsEmpty(A,&globaldata.x,&globaldata.y));

  ipar[0]  = 0;            /* stat: error flag */
  ipar[1]  = lohi;         /* smallest (lohi<0) or largest eigenvalues (lohi>0) */
  PetscCall(PetscBLASIntCast(eps->nev,&ipar[2])); /* number of desired eigenpairs */
  ipar[3]  = 0;            /* number of eigenpairs already converged */
  ipar[4]  = tr->maxlan;   /* maximum Lanczos basis size */
  ipar[5]  = tr->restart;  /* restarting scheme */
  PetscCall(PetscBLASIntCast(eps->max_it,&ipar[6])); /* maximum number of MATVECs */
#if !defined(PETSC_HAVE_MPIUNI)
  fcomm    = MPI_Comm_c2f(PetscObjectComm((PetscObject)eps));
  ipar[7]  = fcomm;
#endif
  ipar[8]  = 0;            /* verboseness */
  ipar[9]  = 99;           /* Fortran IO unit number used to write log messages */
  ipar[10] = 1;            /* use supplied starting vector */
  ipar[11] = 0;            /* checkpointing flag */
  ipar[12] = 98;           /* Fortran IO unit number used to write checkpoint files */
  ipar[13] = 0;            /* number of flops per matvec per PE (not used) */
  tr->work[0] = eps->tol;  /* relative tolerance on residual norms */

  for (i=0;i<eps->ncv;i++) eps->eigr[i]=0.0;
  PetscCall(EPSGetStartVector(eps,0,NULL));
  PetscCall(BVSetActiveColumns(eps->V,0,0));  /* just for deflation space */
  PetscCall(BVGetColumn(eps->V,0,&v0));
  PetscCall(VecGetArray(v0,&pV));

  PetscStackCallExternalVoid("TRLan",TRLan_(MatMult_TRLAN,ipar,&n,&ncv,eps->eigr,pV,&n,tr->work,&tr->lwork));

  PetscCall(VecRestoreArray(v0,&pV));
  PetscCall(BVRestoreColumn(eps->V,0,&v0));

  stat        = ipar[0];
  eps->nconv  = ipar[3];
  eps->its    = ipar[25];
  eps->reason = EPS_CONVERGED_TOL;

  PetscCall(VecDestroy(&globaldata.x));
  PetscCall(VecDestroy(&globaldata.y));
  PetscCheck(stat==0,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Error in TRLAN (code=%" PetscBLASInt_FMT ")",stat);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_TRLAN(EPS eps)
{
  EPS_TRLAN      *tr = (EPS_TRLAN*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscFree(tr->work));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_TRLAN(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_TRLAN(EPS eps)
{
  EPS_TRLAN      *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->ops->solve          = EPSSolve_TRLAN;
  eps->ops->setup          = EPSSetUp_TRLAN;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->destroy        = EPSDestroy_TRLAN;
  eps->ops->reset          = EPSReset_TRLAN;
  eps->ops->backtransform  = EPSBackTransform_Default;
  PetscFunctionReturn(PETSC_SUCCESS);
}
