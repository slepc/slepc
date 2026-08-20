/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This include file contains definitions of system functions
*/

#pragma once

#include <slepcsystypes.h>

#include <slepcmath.h>   /* definition of basic math functions */
#include <slepcsc.h>     /* definition of sorting criteria */
#include <slepcvec.h>    /* utilities related to Vec, extend functionality in PETSc */
#include <slepcmat.h>    /* utilities related to Mat, extend functionality in PETSc */

/* SUBMANSEC = Sys */

/*
   Initialization of SLEPc and other system routines
*/
SLEPC_EXTERN PetscErrorCode SlepcInitialize(int*,char***,const char[],const char[]);
SLEPC_EXTERN PetscErrorCode SlepcInitializeNoPointers(int,char**,const char[],const char[]);
SLEPC_EXTERN PetscErrorCode SlepcInitializeNoArguments(void);
SLEPC_EXTERN PetscErrorCode SlepcFinalize(void);
SLEPC_EXTERN PetscErrorCode SlepcInitializeFortran(void);
SLEPC_EXTERN PetscErrorCode SlepcInitialized(PetscBool*);
SLEPC_EXTERN PetscErrorCode SlepcFinalized(PetscBool*);
SLEPC_EXTERN PetscErrorCode SlepcGetVersion(char[],size_t);
SLEPC_EXTERN PetscErrorCode SlepcGetVersionNumber(PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SlepcGetConfiguration(const char *[]);
SLEPC_EXTERN PetscErrorCode SlepcHasExternalPackage(const char[],PetscBool*);

SLEPC_EXTERN PetscErrorCode SlepcSNPrintfScalar(char[],size_t,PetscScalar,PetscBool);

SLEPC_EXTERN PetscBool SlepcInitializeCalled;
SLEPC_EXTERN PetscBool SlepcFinalizeCalled;

/*MC
   SlepcLogFlopsComplex - Log how many flops are performed in a calculation.

   Synopsis:
   #include <slepcsys.h>
   PetscReal SlepcLogFlopsComplex(PetscLogDouble n)

   Not Collective

   Input parameter:
.  n - the number of flops

   Note:
   Like `PetscLogFlops()` but will log correctly in the case of real scalar computations
   that explicitly use complex arithmetic.

   Level: developer

.seealso: `PetscLogFlops()`
M*/
#if PetscDefined(USE_COMPLEX)
#define SlepcLogFlopsComplex(a) PetscLogFlops((a))
#else
#define SlepcLogFlopsComplex(a) PetscLogFlops(4.0*(a))
#endif

/*MC
   SlepcLogGpuFlopsComplex - Log how many flops are performed in a calculation on the device.

   Synopsis:
   #include <slepcsys.h>
   PetscReal SlepcLogGpuFlopsComplex(PetscLogDouble n)

   Not Collective

   Input parameter:
.  n - the number of flops

   Note:
   Like `PetscLogGpuFlops()` but will log correctly in the case of real scalar computations
   that explicitly use complex arithmetic.

   Level: developer

.seealso: `PetscLogGpuFlops()`
M*/
#if PetscDefined(USE_COMPLEX)
#define SlepcLogGpuFlopsComplex(a) PetscLogGpuFlops((a))
#else
#define SlepcLogGpuFlopsComplex(a) PetscLogGpuFlops(4.0*(a))
#endif

/*
    Developer routines to be used with a debugger
*/
#if PetscDefined(USE_DEBUG)
SLEPC_EXTERN PetscErrorCode SlepcDebugViewMatrix(PetscInt,PetscInt,PetscScalar[],PetscScalar[],PetscInt,const char[],const char[]);
#endif
