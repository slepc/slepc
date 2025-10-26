/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/svdimpl.h>

static PetscBool SVDPackageInitialized = PETSC_FALSE;

const char *SVDTRLanczosGBidiags[] = {"SINGLE","UPPER","LOWER","SVDTRLanczosGBidiag","SVD_TRLANCZOS_GBIDIAG_",NULL};
const char *SVDErrorTypes[] = {"ABSOLUTE","RELATIVE","SVDErrorType","SVD_ERROR_",NULL};
const char *SVDPRIMMEMethods[] = {"","HYBRID","NORMALEQUATIONS","AUGMENTED","SVDPRIMMEMethod","SVD_PRIMME_",NULL};
const char *SVDKSVDEigenMethods[] = {"","MRRR","DC","ELPA","SVDKSVDEigenMethod","SVD_KSVD_EIGEN_",NULL};
const char *SVDKSVDPolarMethods[] = {"","QDWH","ZOLOPD","SVDKSVDPolarMethod","SVD_KSVD_POLAR_",NULL};
const char *const SVDConvergedReasons_Shifted[] = {"","DIVERGED_SYMMETRY_LOST","DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_USER","CONVERGED_MAXIT","SVDConvergedReason","SVD_",NULL};
const char *const*SVDConvergedReasons = SVDConvergedReasons_Shifted + 4;

/*@C
  SVDFinalizePackage - This function destroys everything in the SLEPc interface
  to the `SVD` package. It is called from `SlepcFinalize()`.

  Level: developer

.seealso: [](ch:svd), `SlepcFinalize()`, `SVDInitializePackage()`
@*/
PetscErrorCode SVDFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&SVDList));
  PetscCall(PetscFunctionListDestroy(&SVDMonitorList));
  PetscCall(PetscFunctionListDestroy(&SVDMonitorCreateList));
  PetscCall(PetscFunctionListDestroy(&SVDMonitorDestroyList));
  SVDPackageInitialized       = PETSC_FALSE;
  SVDRegisterAllCalled        = PETSC_FALSE;
  SVDMonitorRegisterAllCalled = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDInitializePackage - This function initializes everything in the `SVD` package.
   It is called from `PetscDLLibraryRegister_slepcsvd()` when using dynamic libraries, and
   on the first call to `SVDCreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](ch:svd), `SVD`, `SlepcInitialize()`, `SVDFinalizePackage()`
@*/
PetscErrorCode SVDInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (SVDPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  SVDPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("SVD Solver",&SVD_CLASSID));
  /* Register Constructors */
  PetscCall(SVDRegisterAll());
  /* Register Monitors */
  PetscCall(SVDMonitorRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("SVDSetUp",SVD_CLASSID,&SVD_SetUp));
  PetscCall(PetscLogEventRegister("SVDSolve",SVD_CLASSID,&SVD_Solve));
  /* Process Info */
  classids[0] = SVD_CLASSID;
  PetscCall(PetscInfoProcessClass("svd",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("svd",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(SVD_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(SVDFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the SVD methods that are in the basic SLEPc libslepcsvd
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcsvd(void)
{
  PetscFunctionBegin;
  PetscCall(SVDInitializePackage());
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */
