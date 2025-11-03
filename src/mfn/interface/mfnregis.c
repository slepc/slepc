/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/mfnimpl.h>  /*I "slepcmfn.h" I*/

SLEPC_EXTERN PetscErrorCode MFNCreate_Krylov(MFN);
SLEPC_EXTERN PetscErrorCode MFNCreate_Expokit(MFN);

/*@C
  MFNRegisterAll - Registers all the matrix functions in the `MFN` package.

  Not Collective

  Level: advanced

.seealso: [](ch:mfn), `MFNRegister()`
@*/
PetscErrorCode MFNRegisterAll(void)
{
  PetscFunctionBegin;
  if (MFNRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  MFNRegisterAllCalled = PETSC_TRUE;
  PetscCall(MFNRegister(MFNKRYLOV,MFNCreate_Krylov));
  PetscCall(MFNRegister(MFNEXPOKIT,MFNCreate_Expokit));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  MFNMonitorRegisterAll - Registers all the monitors in the `MFN` package.

  Not Collective

  Level: advanced

.seealso: [](ch:mfn), `MFNMonitorRegister()`
@*/
PetscErrorCode MFNMonitorRegisterAll(void)
{
  PetscFunctionBegin;
  if (MFNMonitorRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  MFNMonitorRegisterAllCalled = PETSC_TRUE;

  PetscCall(MFNMonitorRegister("error_estimate",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,MFNMonitorDefault,NULL,NULL));
  PetscCall(MFNMonitorRegister("error_estimate",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,MFNMonitorDefaultDrawLG,MFNMonitorDefaultDrawLGCreate,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}
