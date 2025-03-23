/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepcpep.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define pepmonitorset_                    PEPMONITORSET
#define pepmonitorall_                    PEPMONITORALL
#define pepmonitorfirst_                  PEPMONITORFIRST
#define pepmonitorconverged_              PEPMONITORCONVERGED
#define pepmonitorconvergedcreate_        PEPMONITORCONVERGEDCREATE
#define pepmonitorconvergeddestroy_       PEPMONITORCONVERGEDDESTROY
#define pepconvergedabsolute_             PEPCONVERGEDABSOLUTE
#define pepconvergedrelative_             PEPCONVERGEDRELATIVE
#define pepsetconvergencetestfunction_    PEPSETCONVERGENCETESTFUNCTION
#define pepsetstoppingtestfunction_       PEPSETSTOPPINGTESTFUNCTION
#define pepseteigenvaluecomparison_       PEPSETEIGENVALUECOMPARISON
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define pepmonitorset_                    pepmonitorset
#define pepmonitorall_                    pepmonitorall
#define pepmonitorfirst_                  pepmonitorfirst
#define pepmonitorconverged_              pepmonitorconverged
#define pepmonitorconvergedcreate_        pepmonitorconvergedcreate
#define pepmonitorconvergeddestroy_       pepmonitorconvergeddestroy
#define pepconvergedabsolute_             pepconvergedabsolute
#define pepconvergedrelative_             pepconvergedrelative
#define pepsetconvergencetestfunction_    pepsetconvergencetestfunction
#define pepsetstoppingtestfunction_       pepsetstoppingtestfunction
#define pepseteigenvaluecomparison_       pepseteigenvaluecomparison
#endif

/*
   These cannot be called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void pepmonitorall_(PEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);
SLEPC_EXTERN void pepmonitorfirst_(PEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);
SLEPC_EXTERN void pepmonitorconverged_(PEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);

SLEPC_EXTERN void pepmonitorconvergedcreate_(PetscViewer *vin,PetscViewerFormat *format,void *ctx,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(vin,v);
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PEPMonitorConvergedCreate(v,*format,ctx,vf);
}

SLEPC_EXTERN void pepmonitorconvergeddestroy_(PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = PEPMonitorConvergedDestroy(vf);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
  PetscFortranCallbackId convergence;
  PetscFortranCallbackId convdestroy;
  PetscFortranCallbackId stopping;
  PetscFortranCallbackId stopdestroy;
  PetscFortranCallbackId comparison;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(PEP pep,PetscInt i,PetscInt nc,PetscScalar *er,PetscScalar *ei,PetscReal *d,PetscInt l,void* ctx)
{
  PetscObjectUseFortranCallback(pep,_cb.monitor,(PEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),(&pep,&i,&nc,er,ei,d,&l,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void **ctx)
{
  PEP pep = (PEP)*ctx;
  PetscObjectUseFortranCallback(pep,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourconvergence(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscObjectUseFortranCallback(pep,_cb.convergence,(PEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),(&pep,&eigr,&eigi,&res,errest,_ctx,&ierr));
}

static PetscErrorCode ourconvdestroy(void **ctx)
{
  PEP pep = (PEP)*ctx;
  PetscObjectUseFortranCallback(pep,_cb.convdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourstopping(PEP pep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,PEPConvergedReason *reason,void *ctx)
{
  PetscObjectUseFortranCallback(pep,_cb.stopping,(PEP*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PEPConvergedReason*,void*,PetscErrorCode*),(&pep,&its,&max_it,&nconv,&nev,reason,_ctx,&ierr));
}

static PetscErrorCode ourstopdestroy(void **ctx)
{
  PEP pep = (PEP)*ctx;
  PetscObjectUseFortranCallback(pep,_cb.stopdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode oureigenvaluecomparison(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  PEP pep = (PEP)ctx;
  PetscObjectUseFortranCallback(pep,_cb.comparison,(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*,PetscErrorCode*),(&ar,&ai,&br,&bi,r,_ctx,&ierr));
}

SLEPC_EXTERN void pepmonitorset_(PEP *pep,void (*monitor)(PEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)pepmonitorall_) {
    *ierr = PEPMonitorSet(*pep,(PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))PEPMonitorAll,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)pepmonitorconverged_) {
    *ierr = PEPMonitorSet(*pep,(PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))PEPMonitorConverged,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PEPMonitorConvergedDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)pepmonitorfirst_) {
    *ierr = PEPMonitorSet(*pep,(PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))PEPMonitorFirst,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*pep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*pep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = PEPMonitorSet(*pep,ourmonitor,*pep,ourdestroy);
  }
}

SLEPC_EXTERN void pepconvergedabsolute_(PEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*);
SLEPC_EXTERN void pepconvergedrelative_(PEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*);

SLEPC_EXTERN void pepsetconvergencetestfunction_(PEP *pep,void (*func)(PEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)pepconvergedabsolute_) {
    *ierr = PEPSetConvergenceTest(*pep,PEP_CONV_ABS);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)pepconvergedrelative_) {
    *ierr = PEPSetConvergenceTest(*pep,PEP_CONV_REL);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*pep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convergence,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*pep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = PEPSetConvergenceTestFunction(*pep,ourconvergence,*pep,ourconvdestroy);
  }
}

SLEPC_EXTERN void pepstoppingbasic_(PEP*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PEPConvergedReason*,void*,PetscErrorCode*);

SLEPC_EXTERN void pepsetstoppingtestfunction_(PEP *pep,void (*func)(PEP*,PetscInt,PetscInt,PetscInt,PetscInt,PEPConvergedReason*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)pepstoppingbasic_) {
    *ierr = PEPSetStoppingTest(*pep,PEP_STOP_BASIC);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*pep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopping,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*pep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = PEPSetStoppingTestFunction(*pep,ourstopping,*pep,ourstopdestroy);
  }
}

SLEPC_EXTERN void pepseteigenvaluecomparison_(PEP *pep,void (*func)(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*),void* ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PetscObjectSetFortranCallback((PetscObject)*pep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.comparison,(PetscVoidFunction)func,ctx); if (*ierr) return;
  *ierr = PEPSetEigenvalueComparison(*pep,oureigenvaluecomparison,*pep);
}
