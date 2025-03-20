/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file contains the Fortran version of SlepcInitialize()
*/

#include <slepc/private/slepcimpl.h>
#include <petsc/private/ftnimpl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define petscinitializef_             PETSCINITIALIZEF
#define petscfinalize_                PETSCFINALIZE
#define slepcinitializef_             SLEPCINITIALIZEF
#define slepcfinalize_                SLEPCFINALIZE
#define slepcinitializefortran_       SLEPCINITIALIZEFORTRAN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define petscinitializef_             petscinitializef
#define petscfinalize_                petscfinalize
#define slepcinitializef_             slepcinitializef
#define slepcfinalize_                slepcfinalize
#define slepcinitializefortran_       slepcinitializefortran
#endif

SLEPC_EXTERN void petscinitializef_(char *filename,char* help,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len,PETSC_FORTRAN_CHARLEN_T helplen);
SLEPC_EXTERN void petscfinalize_(PetscErrorCode *ierr);

/*
    SlepcInitialize - Version called from Fortran.

    Notes:
    Since this routine is called from Fortran it does not return error codes.
*/
SLEPC_EXTERN void slepcinitializef_(char *filename,char* help,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len,PETSC_FORTRAN_CHARLEN_T helplen)
{
  PetscBool flg;

  if (SlepcInitializeCalled) { *ierr = PETSC_SUCCESS; return; }

  *ierr = PetscInitialized(&flg);
  if (*ierr) { (void)(*PetscErrorPrintf)("SlepcInitialize:PetscInitialized failed");return; }
  if (!flg) {
    petscinitializef_(filename,help,ierr,len,helplen);
    if (*ierr) { (void)(*PetscErrorPrintf)("SlepcInitialize:PetscInitialize failed");return; }
    SlepcBeganPetsc = PETSC_TRUE;
  }

  *ierr = SlepcCitationsInitialize();
  if (*ierr) { (void)(*PetscErrorPrintf)("SlepcInitialize:SlepcCitationsInitialize()\n");return; }

  *ierr = SlepcInitialize_DynamicLibraries();
  if (*ierr) { (void)(*PetscErrorPrintf)("SlepcInitialize:Initializing dynamic libraries\n");return; }

  SlepcInitializeCalled = PETSC_TRUE;
  SlepcFinalizeCalled   = PETSC_FALSE;
  *ierr = PetscInfo(0,"SLEPc successfully started from Fortran\n");
  if (*ierr) { (void)(*PetscErrorPrintf)("SlepcInitialize:Calling PetscInfo()");return; }
}
