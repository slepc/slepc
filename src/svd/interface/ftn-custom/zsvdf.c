/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepcsvd.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define svdmonitorset_                    SVDMONITORSET
#define svdmonitorall_                    SVDMONITORALL
#define svdmonitorfirst_                  SVDMONITORFIRST
#define svdmonitorconditioning_           SVDMONITORCONDITIONING
#define svdmonitorconverged_              SVDMONITORCONVERGED
#define svdmonitorconvergedcreate_        SVDMONITORCONVERGEDCREATE
#define svdmonitorconvergeddestroy_       SVDMONITORCONVERGEDDESTROY
#define svdconvergedabsolute_             SVDCONVERGEDABSOLUTE
#define svdconvergedrelative_             SVDCONVERGEDRELATIVE
#define svdconvergednorm_                 SVDCONVERGEDNORM
#define svdconvergedmaxit_                SVDCONVERGEDMAXIT
#define svdsetconvergencetestfunction_    SVDSETCONVERGENCETESTFUNCTION
#define svdstoppingbasic_                 SVDSTOPPINGBASIC
#define svdstoppingthreshold_             SVDSTOPPINGTHRESHOLD
#define svdsetstoppingtestfunction_       SVDSETSTOPPINGTESTFUNCTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define svdmonitorset_                    svdmonitorset
#define svdmonitorall_                    svdmonitorall
#define svdmonitorfirst_                  svdmonitorfirst
#define svdmonitorconditioning_           svdmonitorconditioning
#define svdmonitorconverged_              svdmonitorconverged
#define svdmonitorconvergedcreate_        svdmonitorconvergedcreate
#define svdmonitorconvergeddestroy_       svdmonitorconvergeddestroy
#define svdconvergedabsolute_             svdconvergedabsolute
#define svdconvergedrelative_             svdconvergedrelative
#define svdconvergednorm_                 svdconvergednorm
#define svdconvergedmaxit_                svdconvergedmaxit
#define svdsetconvergencetestfunction_    svdsetconvergencetestfunction
#define svdstoppingbasic_                 svdstoppingbasic
#define svdstoppingthreshold_             svdstoppingthreshold
#define svdsetstoppingtestfunction_       svdsetstoppingtestfunction
#endif

/*
   These cannot be called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void svdmonitorall_(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);
SLEPC_EXTERN void svdmonitorfirst_(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);
SLEPC_EXTERN void svdmonitorconditioning_(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);
SLEPC_EXTERN void svdmonitorconverged_(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);

SLEPC_EXTERN void svdmonitorconvergedcreate_(PetscViewer *vin,PetscViewerFormat *format,void *ctx,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(vin,v);
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = SVDMonitorConvergedCreate(v,*format,ctx,vf);
}

SLEPC_EXTERN void svdmonitorconvergeddestroy_(PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorConvergedDestroy(vf);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
  PetscFortranCallbackId convergence;
  PetscFortranCallbackId convdestroy;
  PetscFortranCallbackId stopping;
  PetscFortranCallbackId stopdestroy;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(SVD svd,PetscInt i,PetscInt nc,PetscReal *sigma,PetscReal *d,PetscInt l,void *ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.monitor,(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,void*,PetscErrorCode*),(&svd,&i,&nc,sigma,d,&l,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void **ctx)
{
  SVD svd = (SVD)*ctx;
  PetscObjectUseFortranCallback(svd,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourconvergence(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.convergence,(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*),(&svd,&sigma,&res,errest,_ctx,&ierr));
}

static PetscErrorCode ourconvdestroy(void **ctx)
{
  SVD svd = (SVD)*ctx;
  PetscObjectUseFortranCallback(svd,_cb.convdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourstopping(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nsv,SVDConvergedReason *reason,void *ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.stopping,(SVD*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,SVDConvergedReason*,void*,PetscErrorCode*),(&svd,&its,&max_it,&nconv,&nsv,reason,_ctx,&ierr));
}

static PetscErrorCode ourstopdestroy(void **ctx)
{
  SVD svd = (SVD)*ctx;
  PetscObjectUseFortranCallback(svd,_cb.stopdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

SLEPC_EXTERN void svdmonitorset_(SVD *svd,SVDMonitorFn monitor,void *mctx,PetscCtxDestroyFn monitordestroy,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorall_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorAll,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorconverged_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorConverged,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)SVDMonitorConvergedDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorfirst_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorFirst,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = SVDMonitorSet(*svd,ourmonitor,*svd,ourdestroy);
  }
}

SLEPC_EXTERN void svdconvergedabsolute_(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*);
SLEPC_EXTERN void svdconvergedrelative_(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*);
SLEPC_EXTERN void svdconvergednorm_(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*);
SLEPC_EXTERN void svdconvergedmaxit_(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*);

SLEPC_EXTERN void svdsetconvergencetestfunction_(SVD *svd,SVDConvergenceTestFn func,void *ctx,PetscCtxDestroyFn destroy,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergedabsolute_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_ABS);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergedrelative_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_REL);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergednorm_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_NORM);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergedmaxit_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_MAXIT);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convergence,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = SVDSetConvergenceTestFunction(*svd,ourconvergence,*svd,ourconvdestroy);
  }
}

SLEPC_EXTERN void svdstoppingbasic_(SVD*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,SVDConvergedReason*,void*,PetscErrorCode*);
SLEPC_EXTERN void svdstoppingthreshold_(SVD*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,SVDConvergedReason*,void*,PetscErrorCode*);

SLEPC_EXTERN void svdsetstoppingtestfunction_(SVD *svd,SVDStoppingTestFn func,void *ctx,PetscCtxDestroyFn destroy,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)svdstoppingbasic_) {
    *ierr = SVDSetStoppingTest(*svd,SVD_STOP_BASIC);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)svdstoppingthreshold_) {
    *ierr = SVDSetStoppingTest(*svd,SVD_STOP_THRESHOLD);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopping,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = SVDSetStoppingTestFunction(*svd,ourstopping,*svd,ourstopdestroy);
  }
}
