/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/pepimpl.h>      /*I "slepcpep.h" I*/

SLEPC_EXTERN PetscErrorCode PEPCreate_Linear(PEP);
SLEPC_EXTERN PetscErrorCode PEPCreate_QArnoldi(PEP);
SLEPC_EXTERN PetscErrorCode PEPCreate_TOAR(PEP);
SLEPC_EXTERN PetscErrorCode PEPCreate_STOAR(PEP);
SLEPC_EXTERN PetscErrorCode PEPCreate_JD(PEP);
#if defined(PETSC_USE_COMPLEX)
SLEPC_EXTERN PetscErrorCode PEPCreate_CISS(PEP);
#endif

/*@C
   PEPRegisterAll - Registers all the solvers in the `PEP` package.

   Not Collective

   Level: advanced

.seealso: [](ch:pep), `PEPRegister()`
@*/
PetscErrorCode PEPRegisterAll(void)
{
  PetscFunctionBegin;
  if (PEPRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  PEPRegisterAllCalled = PETSC_TRUE;
  PetscCall(PEPRegister(PEPLINEAR,PEPCreate_Linear));
  PetscCall(PEPRegister(PEPQARNOLDI,PEPCreate_QArnoldi));
  PetscCall(PEPRegister(PEPTOAR,PEPCreate_TOAR));
  PetscCall(PEPRegister(PEPSTOAR,PEPCreate_STOAR));
  PetscCall(PEPRegister(PEPJD,PEPCreate_JD));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PEPRegister(PEPCISS,PEPCreate_CISS));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  PEPMonitorRegisterAll - Registers all the monitors in the `PEP` package.

  Not Collective

  Level: advanced

.seealso: [](ch:pep), `PEPMonitorRegister()`
@*/
PetscErrorCode PEPMonitorRegisterAll(void)
{
  PetscFunctionBegin;
  if (PEPMonitorRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  PEPMonitorRegisterAllCalled = PETSC_TRUE;

  PetscCall(PEPMonitorRegister("first_approximation",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,PEPMonitorFirst,NULL,NULL));
  PetscCall(PEPMonitorRegister("first_approximation",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,PEPMonitorFirstDrawLG,PEPMonitorFirstDrawLGCreate,NULL));
  PetscCall(PEPMonitorRegister("all_approximations",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,PEPMonitorAll,NULL,NULL));
  PetscCall(PEPMonitorRegister("all_approximations",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,PEPMonitorAllDrawLG,PEPMonitorAllDrawLGCreate,NULL));
  PetscCall(PEPMonitorRegister("convergence_history",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,PEPMonitorConverged,PEPMonitorConvergedCreate,PEPMonitorConvergedDestroy));
  PetscCall(PEPMonitorRegister("convergence_history",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,PEPMonitorConvergedDrawLG,PEPMonitorConvergedDrawLGCreate,PEPMonitorConvergedDestroy));
  PetscFunctionReturn(PETSC_SUCCESS);
}
