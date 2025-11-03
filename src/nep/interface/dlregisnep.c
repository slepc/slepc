/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/nepimpl.h>

static PetscBool NEPPackageInitialized = PETSC_FALSE;

const char *NEPErrorTypes[] = {"ABSOLUTE","RELATIVE","BACKWARD","NEPErrorType","NEP_ERROR_",NULL};
const char *NEPRefineTypes[] = {"NONE","SIMPLE","MULTIPLE","NEPRefine","NEP_REFINE_",NULL};
const char *NEPRefineSchemes[] = {"","SCHUR","MBE","EXPLICIT","NEPRefineScheme","NEP_REFINE_SCHEME_",NULL};
const char *NEPCISSExtractions[] = {"RITZ","HANKEL","CAA","NEPCISSExtraction","NEP_CISS_EXTRACTION_",NULL};
const char *const NEPConvergedReasons_Shifted[] = {"DIVERGED_SUBSPACE_EXHAUSTED","DIVERGED_LINEAR_SOLVE","","DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_USER","NEPConvergedReason","NEP_",NULL};
const char *const*NEPConvergedReasons = NEPConvergedReasons_Shifted + 5;

/*@C
  NEPFinalizePackage - This function destroys everything in the SLEPc interface
  to the `NEP` package. It is called from `SlepcFinalize()`.

  Level: developer

.seealso: [](ch:nep), `SlepcFinalize()`, `NEPInitializePackage()`
@*/
PetscErrorCode NEPFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&NEPList));
  PetscCall(PetscFunctionListDestroy(&NEPMonitorList));
  PetscCall(PetscFunctionListDestroy(&NEPMonitorCreateList));
  PetscCall(PetscFunctionListDestroy(&NEPMonitorDestroyList));
  NEPPackageInitialized       = PETSC_FALSE;
  NEPRegisterAllCalled        = PETSC_FALSE;
  NEPMonitorRegisterAllCalled = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPInitializePackage - This function initializes everything in the `NEP` package.
   It is called from `PetscDLLibraryRegister_slepcnep()` when using dynamic libraries, and
   on the first call to `NEPCreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](ch:nep), `NEP`, `SlepcInitialize()`, `NEPFinalizePackage()`
@*/
PetscErrorCode NEPInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (NEPPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  NEPPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("NEP Solver",&NEP_CLASSID));
  /* Register Constructors */
  PetscCall(NEPRegisterAll());
  /* Register Monitors */
  PetscCall(NEPMonitorRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("NEPSetUp",NEP_CLASSID,&NEP_SetUp));
  PetscCall(PetscLogEventRegister("NEPSolve",NEP_CLASSID,&NEP_Solve));
  PetscCall(PetscLogEventRegister("NEPRefine",NEP_CLASSID,&NEP_Refine));
  PetscCall(PetscLogEventRegister("NEPFunctionEval",NEP_CLASSID,&NEP_FunctionEval));
  PetscCall(PetscLogEventRegister("NEPJacobianEval",NEP_CLASSID,&NEP_JacobianEval));
  PetscCall(PetscLogEventRegister("NEPResolvent",NEP_CLASSID,&NEP_Resolvent));
  PetscCall(PetscLogEventRegister("NEPCISS_SVD",NEP_CLASSID,&NEP_CISS_SVD));
  /* Process Info */
  classids[0] = NEP_CLASSID;
  PetscCall(PetscInfoProcessClass("nep",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("nep",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(NEP_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(NEPFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the NEP methods that are in the basic SLEPc libslepcnep
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcnep(void)
{
  PetscFunctionBegin;
  PetscCall(NEPInitializePackage());
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */
