/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/mfnimpl.h>

static PetscBool MFNPackageInitialized = PETSC_FALSE;

const char *const MFNConvergedReasons_Shifted[] = {"DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_ITS","MFNConvergedReason","MFN_",NULL};
const char *const*MFNConvergedReasons = MFNConvergedReasons_Shifted + 2;

/*@C
  MFNFinalizePackage - This function destroys everything in the SLEPc interface
  to the `MFN` package. It is called from `SlepcFinalize()`.

  Level: developer

.seealso: [](ch:mfn), `SlepcFinalize()`, `MFNInitializePackage()`
@*/
PetscErrorCode MFNFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&MFNList));
  PetscCall(PetscFunctionListDestroy(&MFNMonitorList));
  PetscCall(PetscFunctionListDestroy(&MFNMonitorCreateList));
  PetscCall(PetscFunctionListDestroy(&MFNMonitorDestroyList));
  MFNPackageInitialized       = PETSC_FALSE;
  MFNRegisterAllCalled        = PETSC_FALSE;
  MFNMonitorRegisterAllCalled = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   MFNInitializePackage - This function initializes everything in the `MFN` package.
   It is called from `PetscDLLibraryRegister_slepcmfn()` when using dynamic libraries, and
   on the first call to `MFNCreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](ch:mfn), `MFN`, `SlepcInitialize()`, `MFNFinalizePackage()`
@*/
PetscErrorCode MFNInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (MFNPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  MFNPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("Matrix Function",&MFN_CLASSID));
  /* Register Constructors */
  PetscCall(MFNRegisterAll());
  /* Register Monitors */
  PetscCall(MFNMonitorRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("MFNSetUp",MFN_CLASSID,&MFN_SetUp));
  PetscCall(PetscLogEventRegister("MFNSolve",MFN_CLASSID,&MFN_Solve));
  /* Process Info */
  classids[0] = MFN_CLASSID;
  PetscCall(PetscInfoProcessClass("mfn",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("mfn",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(MFN_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(MFNFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the MFN methods that are in the basic SLEPc libslepcmfn
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcmfn(void)
{
  PetscFunctionBegin;
  PetscCall(MFNInitializePackage());
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */
