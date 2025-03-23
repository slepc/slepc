!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcsysdef
        use petscmatdef
#include <../src/sys/ftn-mod/slepcsys.h>
        end module

        module slepcsys
        use,intrinsic :: iso_c_binding
        use slepcsysdef
        use petscmat
#include <../src/sys/ftn-mod/slepcsys.h90>
#include <../ftn/sys/slepcall.h90>
        interface SlepcInitialize
          module procedure SlepcInitializeWithHelp, SlepcInitializeNoHelp, SlepcInitializeNoArguments
        end interface
      contains
#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)
!DEC$ ATTRIBUTES DLLEXPORT::SlepcInitializeWithHelp
#endif
      subroutine SlepcInitializeWithHelp(filename,help,ierr)
          character(len=*)           :: filename
          character(len=*)           :: help
          PetscErrorCode             :: ierr

          if (filename .ne. PETSC_NULL_CHARACTER) then
             call SlepcInitializeF(trim(filename),help,ierr)
             CHKERRQ(ierr)
          else
             call SlepcInitializeF(filename,help,ierr)
             CHKERRQ(ierr)
          endif
        end subroutine SlepcInitializeWithHelp

#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)
!DEC$ ATTRIBUTES DLLEXPORT::SlepcInitializeNoHelp
#endif
        subroutine SlepcInitializeNoHelp(filename,ierr)
          character(len=*)           :: filename
          PetscErrorCode             :: ierr

          if (filename .ne. PETSC_NULL_CHARACTER) then
             call SlepcInitializeF(trim(filename),PETSC_NULL_CHARACTER,ierr)
             CHKERRQ(ierr)
          else
             call SlepcInitializeF(filename,PETSC_NULL_CHARACTER,ierr)
             CHKERRQ(ierr)
          endif
        end subroutine SlepcInitializeNoHelp

#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)
!DEC$ ATTRIBUTES DLLEXPORT::SlepcInitializeNoArguments
#endif
        subroutine SlepcInitializeNoArguments(ierr)
          PetscErrorCode             :: ierr

          call SlepcInitializeF(PETSC_NULL_CHARACTER,PETSC_NULL_CHARACTER,ierr)
          CHKERRQ(ierr)
        end subroutine SlepcInitializeNoArguments

#include <../ftn/sys/slepcall.hf90>
        end module
