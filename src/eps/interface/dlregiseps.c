/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/epsimpl.h>

static PetscBool EPSPackageInitialized = PETSC_FALSE;

const char *EPSBalanceTypes[] = {"NONE","ONESIDE","TWOSIDE","USER","EPSBalance","EPS_BALANCE_",NULL};
const char *EPSErrorTypes[] = {"ABSOLUTE","RELATIVE","BACKWARD","EPSErrorType","EPS_ERROR_",NULL};
const char *EPSPowerShiftTypes[] = {"CONSTANT","RAYLEIGH","WILKINSON","EPSPowerShiftType","EPS_POWER_SHIFT_",NULL};
const char *EPSKrylovSchurBSETypes[] = {"SHAO","GRUNING","PROJECTEDBSE","EPSKrylovSchurBSEType","EPS_KRYLOVSCHUR_BSE_",NULL};
const char *EPSLanczosReorthogTypes[] = {"LOCAL","FULL","SELECTIVE","PERIODIC","PARTIAL","DELAYED","EPSLanczosReorthogType","EPS_LANCZOS_REORTHOG_",NULL};
const char *EPSPRIMMEMethods[] = {"","DYNAMIC","DEFAULT_MIN_TIME","DEFAULT_MIN_MATVECS","ARNOLDI","GD","GD_PLUSK","GD_OLSEN_PLUSK","JD_OLSEN_PLUSK","RQI","JDQR","JDQMR","JDQMR_ETOL","SUBSPACE_ITERATION","LOBPCG_ORTHOBASIS","LOBPCG_ORTHOBASISW","EPSPRIMMEMethod","EPS_PRIMME_",NULL};
const char *EPSCISSQuadRules[] = {"(not set yet)","TRAPEZOIDAL","CHEBYSHEV","EPSCISSQuadRule","EPS_CISS_QUADRULE_",NULL};
const char *EPSCISSExtractions[] = {"RITZ","HANKEL","EPSCISSExtraction","EPS_CISS_EXTRACTION_",NULL};
const char *EPSEVSLDOSMethods[] = {"KPM","LANCZOS","EPSEVSLDOSMethod","EPS_EVSL_DOS_",NULL};
const char *EPSEVSLDampings[] = {"NONE","JACKSON","SIGMA","EPSEVSLDamping","EPS_EVSL_DAMPING_",NULL};
const char *const EPSConvergedReasons_Shifted[] = {"","DIVERGED_SYMMETRY_LOST","DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_USER","EPSConvergedReason","EPS_",NULL};
const char *const*EPSConvergedReasons = EPSConvergedReasons_Shifted + 4;

/*@C
  EPSFinalizePackage - This function destroys everything in the SLEPc interface
  to the `EPS` package. It is called from `SlepcFinalize()`.

  Level: developer

.seealso: [](ch:eps), `SlepcFinalize()`, `EPSInitializePackage()`
@*/
PetscErrorCode EPSFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&EPSList));
  PetscCall(PetscFunctionListDestroy(&EPSMonitorList));
  PetscCall(PetscFunctionListDestroy(&EPSMonitorCreateList));
  PetscCall(PetscFunctionListDestroy(&EPSMonitorDestroyList));
  EPSPackageInitialized       = PETSC_FALSE;
  EPSRegisterAllCalled        = PETSC_FALSE;
  EPSMonitorRegisterAllCalled = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSInitializePackage - This function initializes everything in the `EPS` package.
   It is called from `PetscDLLibraryRegister_slepceps()` when using dynamic libraries, and
   on the first call to `EPSCreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](ch:eps), `EPS`, `SlepcInitialize()`, `EPSFinalizePackage()`
@*/
PetscErrorCode EPSInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (EPSPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  EPSPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("EPS Solver",&EPS_CLASSID));
  /* Register Constructors */
  PetscCall(EPSRegisterAll());
  /* Register Monitors */
  PetscCall(EPSMonitorRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("EPSSetUp",EPS_CLASSID,&EPS_SetUp));
  PetscCall(PetscLogEventRegister("EPSSolve",EPS_CLASSID,&EPS_Solve));
  PetscCall(PetscLogEventRegister("EPSCISS_SVD",EPS_CLASSID,&EPS_CISS_SVD));
  /* Process Info */
  classids[0] = EPS_CLASSID;
  PetscCall(PetscInfoProcessClass("eps",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("eps",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(EPS_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(EPSFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the EPS methods that are in the basic SLEPc libslepceps
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepceps(void)
{
  PetscFunctionBegin;
  PetscCall(EPSInitializePackage());
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */
