/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/lmeimpl.h>

static PetscBool LMEPackageInitialized = PETSC_FALSE;

const char *LMEProblemTypes[] = {"LYAPUNOV","SYLVESTER","GEN_LYAPUNOV","GEN_SYLVESTER","DT_LYAPUNOV","STEIN","LMEProblemType","LME_",NULL};
const char *const LMEConvergedReasons_Shifted[] = {"DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","LMEConvergedReason","LME_",NULL};
const char *const*LMEConvergedReasons = LMEConvergedReasons_Shifted + 2;

/*@C
  LMEFinalizePackage - This function destroys everything in the SLEPc interface
  to the `LME` package. It is called from `SlepcFinalize()`.

  Level: developer

.seealso: [](ch:lme), `SlepcFinalize()`, `LMEInitializePackage()`
@*/
PetscErrorCode LMEFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&LMEList));
  PetscCall(PetscFunctionListDestroy(&LMEMonitorList));
  PetscCall(PetscFunctionListDestroy(&LMEMonitorCreateList));
  PetscCall(PetscFunctionListDestroy(&LMEMonitorDestroyList));
  LMEPackageInitialized       = PETSC_FALSE;
  LMERegisterAllCalled        = PETSC_FALSE;
  LMEMonitorRegisterAllCalled = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   LMEInitializePackage - This function initializes everything in the `LME` package.
   It is called from `PetscDLLibraryRegister_slepclme()` when using dynamic libraries, and
   on the first call to `LMECreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](ch:lme), `LME`, `SlepcInitialize()`, `LMEFinalizePackage()`
@*/
PetscErrorCode LMEInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (LMEPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  LMEPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("Lin. Matrix Equation",&LME_CLASSID));
  /* Register Constructors */
  PetscCall(LMERegisterAll());
  /* Register Monitors */
  PetscCall(LMEMonitorRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("LMESetUp",LME_CLASSID,&LME_SetUp));
  PetscCall(PetscLogEventRegister("LMESolve",LME_CLASSID,&LME_Solve));
  PetscCall(PetscLogEventRegister("LMEComputeError",LME_CLASSID,&LME_ComputeError));
  /* Process Info */
  classids[0] = LME_CLASSID;
  PetscCall(PetscInfoProcessClass("lme",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("lme",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(LME_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(LMEFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the LME methods that are in the basic SLEPc libslepclme
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepclme(void)
{
  PetscFunctionBegin;
  PetscCall(LMEInitializePackage());
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */
