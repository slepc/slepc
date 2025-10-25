/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/pepimpl.h>

static PetscBool PEPPackageInitialized = PETSC_FALSE;

const char *PEPBasisTypes[] = {"MONOMIAL","CHEBYSHEV1","CHEBYSHEV2","LEGENDRE","LAGUERRE","HERMITE","PEPBasis","PEP_BASIS_",NULL};
const char *PEPScaleTypes[] = {"NONE","SCALAR","DIAGONAL","BOTH","PEPScale","PEP_SCALE_",NULL};
const char *PEPRefineTypes[] = {"NONE","SIMPLE","MULTIPLE","PEPRefine","PEP_REFINE_",NULL};
const char *PEPRefineSchemes[] = {"","SCHUR","MBE","EXPLICIT","PEPRefineScheme","PEP_REFINE_SCHEME_",NULL};
const char *PEPExtractTypes[] = {"","NONE","NORM","RESIDUAL","STRUCTURED","PEPExtract","PEP_EXTRACT_",NULL};
const char *PEPErrorTypes[] = {"ABSOLUTE","RELATIVE","BACKWARD","PEPErrorType","PEP_ERROR_",NULL};
const char *const PEPConvergedReasons_Shifted[] = {"","DIVERGED_SYMMETRY_LOST","DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_USER","PEPConvergedReason","PEP_",NULL};
const char *const*PEPConvergedReasons = PEPConvergedReasons_Shifted + 4;
const char *PEPJDProjectionTypes[] = {"HARMONIC","ORTHOGONAL","PEPJDProjection","PEP_JD_PROJECTION_",NULL};
const char *PEPCISSExtractions[] = {"RITZ","HANKEL","CAA","PEPCISSExtraction","PEP_CISS_EXTRACTION_",NULL};

/*@C
  PEPFinalizePackage - This function destroys everything in the SLEPc interface
  to the `PEP` package. It is called from `SlepcFinalize()`.

  Level: developer

.seealso: [](ch:pep), `SlepcFinalize()`, `PEPInitializePackage()`
@*/
PetscErrorCode PEPFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&PEPList));
  PetscCall(PetscFunctionListDestroy(&PEPMonitorList));
  PetscCall(PetscFunctionListDestroy(&PEPMonitorCreateList));
  PetscCall(PetscFunctionListDestroy(&PEPMonitorDestroyList));
  PEPPackageInitialized       = PETSC_FALSE;
  PEPRegisterAllCalled        = PETSC_FALSE;
  PEPMonitorRegisterAllCalled = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPInitializePackage - This function initializes everything in the `PEP` package.
   It is called from `PetscDLLibraryRegister_slepcpep()` when using dynamic libraries, and
   on the first call to `PEPCreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](ch:pep), `PEP`, `SlepcInitialize()`, `PEPFinalizePackage()`
@*/
PetscErrorCode PEPInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (PEPPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  PEPPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("PEP Solver",&PEP_CLASSID));
  /* Register Constructors */
  PetscCall(PEPRegisterAll());
  /* Register Monitors */
  PetscCall(PEPMonitorRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("PEPSetUp",PEP_CLASSID,&PEP_SetUp));
  PetscCall(PetscLogEventRegister("PEPSolve",PEP_CLASSID,&PEP_Solve));
  PetscCall(PetscLogEventRegister("PEPRefine",PEP_CLASSID,&PEP_Refine));
  PetscCall(PetscLogEventRegister("PEPCISS_SVD",PEP_CLASSID,&PEP_CISS_SVD));
  /* Process Info */
  classids[0] = PEP_CLASSID;
  PetscCall(PetscInfoProcessClass("pep",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("pep",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(PEP_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(PEPFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the PEP methods that are in the basic SLEPc libslepcpep
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcpep(void)
{
  PetscFunctionBegin;
  PetscCall(PEPInitializePackage());
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */
