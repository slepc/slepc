/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepclme.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define lmemonitordefault_                LMEMONITORDEFAULT
#define lmemonitorset_                    LMEMONITORSET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define lmemonitordefault_                lmemonitordefault
#define lmemonitorset_                    lmemonitorset
#endif

/*
   These cannot be called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void lmemonitordefault_(LME*,PetscInt*,PetscReal*,PetscViewerAndFormat*,PetscErrorCode*);

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(LME lme,PetscInt i,PetscReal d,void* ctx)
{
  PetscObjectUseFortranCallback(lme,_cb.monitor,(LME*,PetscInt*,PetscReal*,void*,PetscErrorCode*),(&lme,&i,&d,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void** ctx)
{
  LME lme = (LME)*ctx;
  PetscObjectUseFortranCallback(lme,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

SLEPC_EXTERN void lmemonitorset_(LME *lme,void (*monitor)(LME*,PetscInt*,PetscReal*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)lmemonitordefault_) {
    *ierr = LMEMonitorSet(*lme,(PetscErrorCode (*)(LME,PetscInt,PetscReal,void*))LMEMonitorDefault,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*lme,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*lme,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = LMEMonitorSet(*lme,ourmonitor,*lme,ourdestroy);
  }
}
