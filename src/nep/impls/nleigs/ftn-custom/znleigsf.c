/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepcnep.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define nepnleigssetsingularitiesfunction_ NEPNLEIGSSETSINGULARITIESFUNCTION
#define nepnleigsgetrkshifts_              NEPNLEIGSGETRKSHIFTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define nepnleigssetsingularitiesfunction_ nepnleigssetsingularitiesfunction
#define nepnleigsgetrkshifts_              nepnleigsgetrkshifts
#endif

static struct {
  PetscFortranCallbackId singularities;
} _cb;

static PetscErrorCode oursingularitiesfunc(NEP nep,PetscInt *maxnp,PetscScalar *xi,void *ctx)
{
  PetscObjectUseFortranCallback(nep,_cb.singularities,(NEP*,PetscInt*,PetscScalar*,void*,PetscErrorCode*),(&nep,maxnp,xi,_ctx,&ierr));
}

SLEPC_EXTERN void nepnleigssetsingularitiesfunction_(NEP *nep,void (*func)(NEP*,PetscInt*,PetscScalar*,void*,PetscErrorCode*),void *ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.singularities,(PetscVoidFunction)func,ctx); if (*ierr) return;
  *ierr = NEPNLEIGSSetSingularitiesFunction(*nep,oursingularitiesfunc,*nep);
}

SLEPC_EXTERN void nepnleigsgetrkshifts_(NEP *nep,PetscInt *ns,PetscScalar *pshifts,PetscErrorCode *ierr)
{
  PetscScalar *oshifts;
  PetscInt    n;

  CHKFORTRANNULLSCALAR(pshifts);
  *ierr = NEPNLEIGSGetRKShifts(*nep,&n,&oshifts); if (*ierr) return;
  if (pshifts) { *ierr = PetscArraycpy(pshifts,oshifts,n); if (*ierr) return; }
  *ns = n;
  *ierr = PetscFree(oshifts);
}
