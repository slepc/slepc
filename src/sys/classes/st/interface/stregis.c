/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/stimpl.h>          /*I   "slepcst.h"   I*/

SLEPC_EXTERN PetscErrorCode STCreate_Shell(ST);
SLEPC_EXTERN PetscErrorCode STCreate_Shift(ST);
SLEPC_EXTERN PetscErrorCode STCreate_Sinvert(ST);
SLEPC_EXTERN PetscErrorCode STCreate_Cayley(ST);
SLEPC_EXTERN PetscErrorCode STCreate_Precond(ST);
SLEPC_EXTERN PetscErrorCode STCreate_Filter(ST);

/*@C
   STRegisterAll - Registers all of the spectral transformations in the `ST` package.

   Not Collective

   Level: advanced

.seealso: [](ch:st), `STRegister()`
@*/
PetscErrorCode STRegisterAll(void)
{
  PetscFunctionBegin;
  if (STRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  STRegisterAllCalled = PETSC_TRUE;
  PetscCall(STRegister(STSHELL,STCreate_Shell));
  PetscCall(STRegister(STSHIFT,STCreate_Shift));
  PetscCall(STRegister(STSINVERT,STCreate_Sinvert));
  PetscCall(STRegister(STCAYLEY,STCreate_Cayley));
  PetscCall(STRegister(STPRECOND,STCreate_Precond));
  PetscCall(STRegister(STFILTER,STCreate_Filter));
  PetscFunctionReturn(PETSC_SUCCESS);
}
