/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepceps.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define epsmonitorset_                    EPSMONITORSET
#define epsmonitorall_                    EPSMONITORALL
#define epsmonitorfirst_                  EPSMONITORFIRST
#define epsmonitorconverged_              EPSMONITORCONVERGED
#define epsmonitorconvergedcreate_        EPSMONITORCONVERGEDCREATE
#define epsmonitorconvergeddestroy_       EPSMONITORCONVERGEDDESTROY
#define epsconvergedabsolute_             EPSCONVERGEDABSOLUTE
#define epsconvergedrelative_             EPSCONVERGEDRELATIVE
#define epsconvergednorm_                 EPSCONVERGEDNORM
#define epssetconvergencetestfunction_    EPSSETCONVERGENCETESTFUNCTION
#define epsstoppingbasic_                 EPSSTOPPINGBASIC
#define epsstoppingthreshold_             EPSSTOPPINGTHRESHOLD
#define epssetstoppingtestfunction_       EPSSETSTOPPINGTESTFUNCTION
#define epsseteigenvaluecomparison_       EPSSETEIGENVALUECOMPARISON
#define epssetarbitraryselection_         EPSSETARBITRARYSELECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define epsmonitorset_                    epsmonitorset
#define epsmonitorall_                    epsmonitorall
#define epsmonitorfirst_                  epsmonitorfirst
#define epsmonitorconverged_              epsmonitorconverged
#define epsmonitorconvergedcreate_        epsmonitorconvergedcreate
#define epsmonitorconvergeddestroy_       epsmonitorconvergeddestroy
#define epsconvergedabsolute_             epsconvergedabsolute
#define epsconvergedrelative_             epsconvergedrelative
#define epsconvergednorm_                 epsconvergednorm
#define epssetconvergencetestfunction_    epssetconvergencetestfunction
#define epsstoppingbasic_                 epsstoppingbasic
#define epsstoppingthreshold_             epsstoppingthreshold
#define epssetstoppingtestfunction_       epssetstoppingtestfunction
#define epsseteigenvaluecomparison_       epsseteigenvaluecomparison
#define epssetarbitraryselection_         epssetarbitraryselection
#endif

/*
   These are not usually called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void epsmonitorall_(EPS *eps,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorAll(*eps,*it,*nconv,eigr,eigi,errest,*nest,*vf);
}

SLEPC_EXTERN void epsmonitorfirst_(EPS *eps,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorFirst(*eps,*it,*nconv,eigr,eigi,errest,*nest,*vf);
}

SLEPC_EXTERN void epsmonitorconverged_(EPS *eps,PetscInt *it,PetscInt *nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorConverged(*eps,*it,*nconv,eigr,eigi,errest,*nest,*vf);
}

SLEPC_EXTERN void epsmonitorconvergedcreate_(PetscViewer *vin,PetscViewerFormat *format,void *ctx,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(vin,v);
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = EPSMonitorConvergedCreate(v,*format,ctx,vf);
}

SLEPC_EXTERN void epsmonitorconvergeddestroy_(PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = EPSMonitorConvergedDestroy(vf);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
  PetscFortranCallbackId convergence;
  PetscFortranCallbackId convdestroy;
  PetscFortranCallbackId stopping;
  PetscFortranCallbackId stopdestroy;
  PetscFortranCallbackId comparison;
  PetscFortranCallbackId arbitrary;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(EPS eps,PetscInt i,PetscInt nc,PetscScalar *er,PetscScalar *ei,PetscReal *d,PetscInt l,void* ctx)
{
  PetscObjectUseFortranCallback(eps,_cb.monitor,(EPS*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),(&eps,&i,&nc,er,ei,d,&l,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void **ctx)
{
  EPS eps = (EPS)*ctx;
  PetscObjectUseFortranCallback(eps,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourconvergence(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscObjectUseFortranCallback(eps,_cb.convergence,(EPS*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),(&eps,&eigr,&eigi,&res,errest,_ctx,&ierr));
}

static PetscErrorCode ourconvdestroy(void **ctx)
{
  EPS eps = (EPS)*ctx;
  PetscObjectUseFortranCallback(eps,_cb.convdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourstopping(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx)
{
  PetscObjectUseFortranCallback(eps,_cb.stopping,(EPS*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,EPSConvergedReason*,void*,PetscErrorCode*),(&eps,&its,&max_it,&nconv,&nev,reason,_ctx,&ierr));
}

static PetscErrorCode ourstopdestroy(void **ctx)
{
  EPS eps = (EPS)*ctx;
  PetscObjectUseFortranCallback(eps,_cb.stopdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode oureigenvaluecomparison(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  EPS eps = (EPS)ctx;
  PetscObjectUseFortranCallback(eps,_cb.comparison,(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*,PetscErrorCode*),(&ar,&ai,&br,&bi,r,_ctx,&ierr));
}

static PetscErrorCode ourarbitraryfunc(PetscScalar er,PetscScalar ei,Vec xr,Vec xi,PetscScalar *rr,PetscScalar *ri,void *ctx)
{
  EPS eps = (EPS)ctx;
  PetscObjectUseFortranCallback(eps,_cb.arbitrary,(PetscScalar*,PetscScalar*,Vec*,Vec*,PetscScalar*,PetscScalar*,void*,PetscErrorCode*),(&er,&ei,&xr,&xi,rr,ri,_ctx,&ierr));
}

SLEPC_EXTERN void epsmonitorset_(EPS *eps,void (*monitor)(EPS*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)epsmonitorall_) {
    *ierr = EPSMonitorSet(*eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))EPSMonitorAll,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)epsmonitorconverged_) {
    *ierr = EPSMonitorSet(*eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))EPSMonitorConverged,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)EPSMonitorConvergedDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)epsmonitorfirst_) {
    *ierr = EPSMonitorSet(*eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))EPSMonitorFirst,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = EPSMonitorSet(*eps,ourmonitor,*eps,ourdestroy);
  }
}

SLEPC_EXTERN void epsconvergedabsolute_(EPS *eps,PetscScalar *eigr,PetscScalar *eigi,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSConvergedAbsolute(*eps,*eigr,*eigi,*res,errest,ctx);
}

SLEPC_EXTERN void epsconvergedrelative_(EPS *eps,PetscScalar *eigr,PetscScalar *eigi,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSConvergedRelative(*eps,*eigr,*eigi,*res,errest,ctx);
}

SLEPC_EXTERN void epsconvergednorm_(EPS *eps,PetscScalar *eigr,PetscScalar *eigi,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSConvergedNorm(*eps,*eigr,*eigi,*res,errest,ctx);
}

SLEPC_EXTERN void epssetconvergencetestfunction_(EPS *eps,void (*func)(EPS*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)epsconvergedabsolute_) {
    *ierr = EPSSetConvergenceTest(*eps,EPS_CONV_ABS);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)epsconvergedrelative_) {
    *ierr = EPSSetConvergenceTest(*eps,EPS_CONV_REL);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)epsconvergednorm_) {
    *ierr = EPSSetConvergenceTest(*eps,EPS_CONV_NORM);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convergence,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = EPSSetConvergenceTestFunction(*eps,ourconvergence,*eps,ourconvdestroy);
  }
}

SLEPC_EXTERN void epsstoppingbasic_(EPS *eps,PetscInt *its,PetscInt *max_it,PetscInt *nconv,PetscInt *nev,EPSConvergedReason *reason,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSStoppingBasic(*eps,*its,*max_it,*nconv,*nev,reason,ctx);
}

SLEPC_EXTERN void epsstoppingthreshold_(EPS *eps,PetscInt *its,PetscInt *max_it,PetscInt *nconv,PetscInt *nsv,EPSConvergedReason *reason,void *ctx,PetscErrorCode *ierr)
{
  *ierr = EPSStoppingThreshold(*eps,*its,*max_it,*nconv,*nsv,reason,ctx);
}

SLEPC_EXTERN void epssetstoppingtestfunction_(EPS *eps,void (*func)(EPS*,PetscInt,PetscInt,PetscInt,PetscInt,EPSConvergedReason*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)epsstoppingbasic_) {
    *ierr = EPSSetStoppingTest(*eps,EPS_STOP_BASIC);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)epsstoppingthreshold_) {
    *ierr = EPSSetStoppingTest(*eps,EPS_STOP_THRESHOLD);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopping,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = EPSSetStoppingTestFunction(*eps,ourstopping,*eps,ourstopdestroy);
  }
}

SLEPC_EXTERN void epsseteigenvaluecomparison_(EPS *eps,void (*func)(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*),void* ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.comparison,(PetscVoidFunction)func,ctx); if (*ierr) return;
  *ierr = EPSSetEigenvalueComparison(*eps,oureigenvaluecomparison,*eps);
}

SLEPC_EXTERN void epssetarbitraryselection_(EPS *eps,void (*func)(PetscScalar*,PetscScalar*,Vec*,Vec*,PetscScalar*,PetscScalar*,void*,PetscErrorCode*),void *ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PetscObjectSetFortranCallback((PetscObject)*eps,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.arbitrary,(PetscVoidFunction)func,ctx); if (*ierr) return;
  *ierr = EPSSetArbitrarySelection(*eps,ourarbitraryfunc,*eps);
}
