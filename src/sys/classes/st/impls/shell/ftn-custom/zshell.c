/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepcst.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define stshellsetapply_                   STSHELLSETAPPLY
#define stshellsetapplytranspose_          STSHELLSETAPPLYTRANSPOSE
#define stshellsetapplyhermitiantranspose_ STSHELLSETAPPLYHERMITIANTRANSPOSE
#define stshellsetbacktransform_           STSHELLSETBACKTRANSFORM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define stshellsetapply_                   stshellsetapply
#define stshellsetapplytranspose_          stshellsetapplytranspose
#define stshellsetapplyhermitiantranspose_ stshellsetapplyhermitiantranspose
#define stshellsetbacktransform_           stshellsetbacktransform
#endif

static struct {
  PetscFortranCallbackId apply;
  PetscFortranCallbackId applytranspose;
  PetscFortranCallbackId applyhermtrans;
  PetscFortranCallbackId backtransform;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourshellapply(ST st,Vec x,Vec y)
{
  PetscObjectUseFortranCallback(st,_cb.apply,(ST*,Vec*,Vec*,PetscErrorCode*),(&st,&x,&y,&ierr));
}

static PetscErrorCode ourshellapplytranspose(ST st,Vec x,Vec y)
{
  PetscObjectUseFortranCallback(st,_cb.applytranspose,(ST*,Vec*,Vec*,PetscErrorCode*),(&st,&x,&y,&ierr));
}

static PetscErrorCode ourshellapplyhermitiantranspose(ST st,Vec x,Vec y)
{
  PetscObjectUseFortranCallback(st,_cb.applyhermtrans,(ST*,Vec*,Vec*,PetscErrorCode*),(&st,&x,&y,&ierr));
}

static PetscErrorCode ourshellbacktransform(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscObjectUseFortranCallback(st,_cb.backtransform,(ST*,PetscInt*,PetscScalar*,PetscScalar*,PetscErrorCode*),(&st,&n,eigr,eigi,&ierr));
}

SLEPC_EXTERN void stshellsetapply_(ST *st,void (*apply)(void*,Vec*,Vec*,PetscErrorCode*),PetscErrorCode *ierr)
{
  *ierr = PetscObjectSetFortranCallback((PetscObject)*st,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.apply,(PetscVoidFunction)apply,NULL); if (*ierr) return;
  *ierr = STShellSetApply(*st,ourshellapply);
}

SLEPC_EXTERN void stshellsetapplytranspose_(ST *st,void (*applytranspose)(void*,Vec*,Vec*,PetscErrorCode*),PetscErrorCode *ierr)
{
  *ierr = PetscObjectSetFortranCallback((PetscObject)*st,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.applytranspose,(PetscVoidFunction)applytranspose,NULL); if (*ierr) return;
  *ierr = STShellSetApplyTranspose(*st,ourshellapplytranspose);
}

SLEPC_EXTERN void stshellsetapplyhermitiantranspose_(ST *st,void (*applyhermtrans)(void*,Vec*,Vec*,PetscErrorCode*),PetscErrorCode *ierr)
{
  *ierr = PetscObjectSetFortranCallback((PetscObject)*st,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.applyhermtrans,(PetscVoidFunction)applyhermtrans,NULL); if (*ierr) return;
  *ierr = STShellSetApplyHermitianTranspose(*st,ourshellapplyhermitiantranspose);
}

SLEPC_EXTERN void stshellsetbacktransform_(ST *st,void (*backtransform)(void*,PetscScalar*,PetscScalar*,PetscErrorCode*),PetscErrorCode *ierr)
{
  *ierr = PetscObjectSetFortranCallback((PetscObject)*st,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.backtransform,(PetscVoidFunction)backtransform,NULL); if (*ierr) return;
  *ierr = STShellSetBackTransform(*st,ourshellbacktransform);
}
