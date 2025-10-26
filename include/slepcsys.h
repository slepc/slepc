/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This include file contains definitions of system functions. It is included
   by all other SLEPc include files.
*/

#pragma once

#include <petscsys.h>

/* SUBMANSEC = Sys */

#if defined(slepc_EXPORTS)
#define SLEPC_VISIBILITY_PUBLIC PETSC_DLLEXPORT
#else
#define SLEPC_VISIBILITY_PUBLIC PETSC_DLLIMPORT
#endif
#define SLEPC_VISIBILITY_INTERNAL PETSC_VISIBILITY_INTERNAL

/*
    Functions tagged with SLEPC_EXTERN in the header files are
  always defined as extern "C" when compiled with C++ so they may be
  used from C and are always visible in the shared libraries
*/
#if defined(__cplusplus)
#define SLEPC_EXTERN extern "C" SLEPC_VISIBILITY_PUBLIC
#define SLEPC_INTERN extern "C" SLEPC_VISIBILITY_INTERNAL
#else
#define SLEPC_EXTERN extern SLEPC_VISIBILITY_PUBLIC
#define SLEPC_INTERN extern SLEPC_VISIBILITY_INTERNAL
#endif

#if defined(PETSC_USE_SINGLE_LIBRARY)
  #define SLEPC_SINGLE_LIBRARY_VISIBILITY_INTERNAL SLEPC_VISIBILITY_INTERNAL
  #define SLEPC_SINGLE_LIBRARY_INTERN              SLEPC_INTERN
#else
  #define SLEPC_SINGLE_LIBRARY_VISIBILITY_INTERNAL SLEPC_VISIBILITY_PUBLIC
  #define SLEPC_SINGLE_LIBRARY_INTERN              SLEPC_EXTERN
#endif

/* ========================================================================== */
/*
   slepcconf.h is created by the configure script and placed in ${PETSC_ARCH}/include.
   It contains macro definitions set at configure time.
*/
#include <slepcconf.h>
/*
    slepcversion.h contains version info
*/
#include <slepcversion.h>
#define SLEPC_AUTHOR_INFO "       The SLEPc Team\n    slepc-maint@upv.es\n https://slepc.upv.es\n"

/* ========================================================================== */
/*
   The PETSc include files.
*/
#include <petscmat.h>
/*
    slepcmath.h contains definition of basic math functions
*/
#include <slepcmath.h>
/*
    slepcsc.h contains definition of sorting criterion
*/
#include <slepcsc.h>
/*
    slepcmat.h, slepcvec.h contain utilities related to Mat and Vec, extend functionality in PETSc
*/
#include <slepcmat.h>
#include <slepcvec.h>

/*
    Context for monitors of type XXXMonitorConverged
*/
typedef struct _n_SlepcConvMon* SlepcConvMon;

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
SLEPC_EXTERN PetscErrorCode SlepcHasExternalPackage(const char[],PetscBool*);

SLEPC_EXTERN PetscErrorCode SlepcSNPrintfScalar(char[],size_t,PetscScalar,PetscBool);

SLEPC_EXTERN PetscBool SlepcInitializeCalled;
SLEPC_EXTERN PetscBool SlepcFinalizeCalled;

#if defined(PETSC_USE_COMPLEX)
#define SlepcLogFlopsComplex(a) PetscLogFlops((a))
#else
#define SlepcLogFlopsComplex(a) PetscLogFlops((4.0*a))
#endif

#if defined(PETSC_USE_COMPLEX)
#define SlepcLogGpuFlopsComplex(a) PetscLogGpuFlops((a))
#else
#define SlepcLogGpuFlopsComplex(a) PetscLogGpuFlops((4.0*a))
#endif

/*
    Developer routines to be used with a debugger
*/
#if defined(PETSC_USE_DEBUG)
SLEPC_EXTERN PetscErrorCode SlepcDebugViewMatrix(PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscInt,const char*,const char*);
#endif
