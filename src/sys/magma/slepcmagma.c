/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepcimpl.h>
#include <slepcmagma.h>

static PetscBool SlepcBeganMagma = PETSC_FALSE;

static void slepc_magma_finalize(void PETSC_UNUSED *unused, magma_int_t *ierr)
{
  (void)unused;
  *ierr = magma_finalize();
  return;
}

static PetscErrorCode SlepcMagmaFinalize(void)
{
  PetscFunctionBegin;
  SlepcBeganMagma = PETSC_FALSE;
  PetscCallMAGMA(slepc_magma_finalize, NULL);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static void slepc_magma_init(void PETSC_UNUSED *unused, magma_int_t *ierr)
{
  (void)unused;
  *ierr = magma_init();
  return;
}

/*@
   SlepcMagmaInit - A helper function to cleanly initialize MAGMA from within SLEPc
   when needed.

   Logically Collective

   Level: developer

.seealso: `SlepcInitialize()`
@*/
PetscErrorCode SlepcMagmaInit(void)
{
  PetscFunctionBegin;
  if (!SlepcBeganMagma) {
    PetscCallMAGMA(slepc_magma_init, NULL);
    SlepcBeganMagma = PETSC_TRUE;
    PetscCall(PetscRegisterFinalize(SlepcMagmaFinalize));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
