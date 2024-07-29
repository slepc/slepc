!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test14f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Simple example that tests solving a DSNHEP problem.
!
!  The command line options are:
!    -n <n>, where <n> = matrix size
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcds.h>
      use slepcds
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     A     problem matrix
!     ds    dense solver context

      Mat            A
      DS             ds
      PetscInt       n, i, ld, zero
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg
      PetscScalar    wr(100), wi(100)
      PetscReal      re, im
      PetscScalar, pointer :: aa(:,:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      zero = 0
      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 10
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      if (n .gt. 100) then; SETERRA(PETSC_COMM_SELF,1,'Program currently limited to n=100'); endif

      if (rank .eq. 0) then
        write(*,110) n
      endif
 110  format (/'Solve a Dense System of type NHEP, n =',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create DS object
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(DSCreate(PETSC_COMM_WORLD,ds,ierr))
      PetscCallA(DSSetType(ds,DSNHEP,ierr))
      PetscCallA(DSSetFromOptions(ds,ierr))
      ld = n
      PetscCallA(DSAllocate(ds,ld,ierr))
      PetscCallA(DSSetDimensions(ds,n,zero,zero,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Fill with Grcar matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(DSGetMat(ds,DS_MAT_A,A,ierr))
      PetscCallA(MatDenseGetArray(A,aa,ierr))
      call FillUpMatrix(n,aa)
      PetscCallA(MatDenseRestoreArray(A,aa,ierr))
      PetscCallA(DSRestoreMat(ds,DS_MAT_A,A,ierr))
      PetscCallA(DSSetState(ds,DS_STATE_INTERMEDIATE,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the problem and show eigenvalues
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(DSSolve(ds,wr,wi,ierr))
!     PetscCallA(DSSort(ds,wr,wi,PETSC_NULL_SCALAR,PETSC_NULL_SCALAR,PETSC_NULL_INTEGER,ierr))

      if (rank .eq. 0) then
        write(*,*) 'Computed eigenvalues ='
        do i=1,n
#if defined(PETSC_USE_COMPLEX)
          re = PetscRealPart(wr(i))
          im = PetscImaginaryPart(wr(i))
#else
          re = wr(i)
          im = wi(i)
#endif
          if (abs(im).lt.1.d-10) then
            write(*,120) re
          else
            write(*,130) re, im
          endif
        end do
      endif
 120  format ('  ',F8.5)
 130  format ('  ',F8.5,SP,F8.5,'i')

!     *** Clean up
      PetscCallA(DSDestroy(ds,ierr))
      PetscCallA(SlepcFinalize(ierr))
      end

! -----------------------------------------------------------------

      subroutine FillUpMatrix(n,X)
      PetscInt    n,i,j
      PetscScalar X(n,n)

      do i=2,n
        X(i,i-1) = -1.d0
      end do
      do j=0,3
        do i=1,n-j
          X(i,i+j) = 1.d0
        end do
      end do

      end

!/*TEST
!
!   test:
!      suffix: 1
!      requires: !complex
!
!TEST*/
