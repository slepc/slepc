/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <petscmacros.h>
#include <slepcversion.h>

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

#if PetscDefined(USE_SINGLE_LIBRARY)
  #define SLEPC_SINGLE_LIBRARY_VISIBILITY_INTERNAL SLEPC_VISIBILITY_INTERNAL
  #define SLEPC_SINGLE_LIBRARY_INTERN              SLEPC_INTERN
#else
  #define SLEPC_SINGLE_LIBRARY_VISIBILITY_INTERNAL SLEPC_VISIBILITY_PUBLIC
  #define SLEPC_SINGLE_LIBRARY_INTERN              SLEPC_EXTERN
#endif

#define SLEPC_AUTHOR_INFO "       The SLEPc Team\n    slepc-maint@upv.es\n https://slepc.upv.es\n"
