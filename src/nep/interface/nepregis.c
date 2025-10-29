/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/nepimpl.h>      /*I "slepcnep.h" I*/

SLEPC_EXTERN PetscErrorCode NEPCreate_RII(NEP);
SLEPC_EXTERN PetscErrorCode NEPCreate_SLP(NEP);
SLEPC_EXTERN PetscErrorCode NEPCreate_NArnoldi(NEP);
SLEPC_EXTERN PetscErrorCode NEPCreate_Interpol(NEP);
#if defined(PETSC_USE_COMPLEX)
SLEPC_EXTERN PetscErrorCode NEPCreate_CISS(NEP);
#endif
SLEPC_EXTERN PetscErrorCode NEPCreate_NLEIGS(NEP);

/*@C
   NEPRegisterAll - Registers all the solvers in the `NEP` package.

   Not Collective

   Level: advanced

.seealso: [](ch:nep), `NEPRegister()`
@*/
PetscErrorCode NEPRegisterAll(void)
{
  PetscFunctionBegin;
  if (NEPRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  NEPRegisterAllCalled = PETSC_TRUE;
  PetscCall(NEPRegister(NEPRII,NEPCreate_RII));
  PetscCall(NEPRegister(NEPSLP,NEPCreate_SLP));
  PetscCall(NEPRegister(NEPNARNOLDI,NEPCreate_NArnoldi));
  PetscCall(NEPRegister(NEPINTERPOL,NEPCreate_Interpol));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(NEPRegister(NEPCISS,NEPCreate_CISS));
#endif
  PetscCall(NEPRegister(NEPNLEIGS,NEPCreate_NLEIGS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  NEPMonitorRegisterAll - Registers all the monitors in the `NEP` package.

  Not Collective

  Level: advanced

.seealso: [](ch:nep), `NEPMonitorRegister()`
@*/
PetscErrorCode NEPMonitorRegisterAll(void)
{
  PetscFunctionBegin;
  if (NEPMonitorRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  NEPMonitorRegisterAllCalled = PETSC_TRUE;

  PetscCall(NEPMonitorRegister("first_approximation",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,NEPMonitorFirst,NULL,NULL));
  PetscCall(NEPMonitorRegister("first_approximation",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,NEPMonitorFirstDrawLG,NEPMonitorFirstDrawLGCreate,NULL));
  PetscCall(NEPMonitorRegister("all_approximations",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,NEPMonitorAll,NULL,NULL));
  PetscCall(NEPMonitorRegister("all_approximations",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,NEPMonitorAllDrawLG,NEPMonitorAllDrawLGCreate,NULL));
  PetscCall(NEPMonitorRegister("convergence_history",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,NEPMonitorConverged,NEPMonitorConvergedCreate,NEPMonitorConvergedDestroy));
  PetscCall(NEPMonitorRegister("convergence_history",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,NEPMonitorConvergedDrawLG,NEPMonitorConvergedDrawLGCreate,NEPMonitorConvergedDestroy));
  PetscFunctionReturn(PETSC_SUCCESS);
}
