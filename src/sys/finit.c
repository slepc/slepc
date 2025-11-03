/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepcimpl.h>           /*I "slepcsys.h" I*/

/*@C
   SlepcInitializeFortran - Routine that should be called from C after
   the call to `SlepcInitialize()` if one is using a C main program
   that calls Fortran routines that in turn call SLEPc routines.

   Collective on `PETSC_COMM_WORLD`

   Level: advanced

   Note:
   `SlepcInitializeFortran()` initializes some of the default SLEPc variables
   for use in Fortran if a user's main program is written in C.
   `SlepcInitializeFortran()` is NOT needed if a user's main
   program is written in Fortran; in this case, just calling
   `SlepcInitialize()` in the main (Fortran) program is sufficient.

.seealso: [](sec:fortran), `SlepcInitialize()`
@*/
PetscErrorCode SlepcInitializeFortran(void)
{
#if defined(PETSC_USE_FORTRAN_BINDINGS)
  PetscCall(PetscInitializeFortran());
#endif
  return PETSC_SUCCESS;
}

