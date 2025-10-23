/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/lmeimpl.h>  /*I "slepclme.h" I*/

SLEPC_EXTERN PetscErrorCode LMECreate_Krylov(LME);

/*@C
  LMERegisterAll - Registers all the matrix functions in the `LME` package.

  Not Collective

  Level: advanced

.seealso: [](ch:lme), `LMERegister()`
@*/
PetscErrorCode LMERegisterAll(void)
{
  PetscFunctionBegin;
  if (LMERegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  LMERegisterAllCalled = PETSC_TRUE;
  PetscCall(LMERegister(LMEKRYLOV,LMECreate_Krylov));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  LMEMonitorRegisterAll - Registers all the monitors in the `LME `package.

  Not Collective

  Level: advanced

.seealso: [](ch:lme), `LMEMonitorRegister()`
@*/
PetscErrorCode LMEMonitorRegisterAll(void)
{
  PetscFunctionBegin;
  if (LMEMonitorRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  LMEMonitorRegisterAllCalled = PETSC_TRUE;

  PetscCall(LMEMonitorRegister("error_estimate",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,LMEMonitorDefault,NULL,NULL));
  PetscCall(LMEMonitorRegister("error_estimate",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,LMEMonitorDefaultDrawLG,LMEMonitorDefaultDrawLGCreate,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}
