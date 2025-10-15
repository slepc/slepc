!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex1f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Simple example that solves an eigensystem with the EPS object.
!  The standard symmetric eigenvalue problem to be solved corresponds to the
!  Laplacian operator in 1 dimension.
!
!  The command line options are:
!    -n <n>, where <n> = number of grid points = matrix size
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepceps.h>
      program ex1f
      use slepceps
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     A      operator matrix
!     eps    eigenproblem solver context

      Mat            :: A
      EPS            :: eps
      EPSType        :: tname
      PetscInt       :: n, i, Istart, Iend, one, two, three
      PetscInt       :: nev
      PetscInt       :: row(1), col(3)
      PetscMPIInt    :: rank
      PetscErrorCode :: ierr
      PetscBool      :: flg, terse
      PetscScalar    :: val(3)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      one = 1
      two = 2
      three = 3
      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,"ex1f test"//c_new_line,ierr))
      if (ierr/=0) then
        print*,'SlepcInitialize failed'
        stop
      end if
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 30
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))

      if (rank==0) then
        write(*,100) n
      end if
 100  format (/'1-D Laplacian Eigenproblem, n =',I4,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the operator matrix that defines the eigensystem, Ax=kx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))

      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      if (Istart==0) then
        row(1) = 0
        col(1) = 0
        col(2) = 1
        val(1) =  2.0
        val(2) = -1.0
        PetscCallA(MatSetValues(A,one,row,two,col,val,INSERT_VALUES,ierr))
        Istart = Istart+1
      end if
      if (Iend==n) then
        row(1) = n-1
        col(1) = n-2
        col(2) = n-1
        val(1) = -1.0
        val(2) =  2.0
        PetscCallA(MatSetValues(A,one,row,two,col,val,INSERT_VALUES,ierr))
        Iend = Iend-1
      end if
      val(1) = -1.0
      val(2) =  2.0
      val(3) = -1.0
      do i=Istart,Iend-1
        row(1) = i
        col(1) = i-1
        col(2) = i
        col(3) = i+1
        PetscCallA(MatSetValues(A,one,row,three,col,val,INSERT_VALUES,ierr))
      end do

      PetscCallA(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      PetscCallA(EPSCreate(PETSC_COMM_WORLD,eps,ierr))

!     ** Set operators. In this case, it is a standard eigenvalue problem
      PetscCallA(EPSSetOperators(eps,A,PETSC_NULL_MAT,ierr))
      PetscCallA(EPSSetProblemType(eps,EPS_HEP,ierr))

!     ** Set solver parameters at runtime
      PetscCallA(EPSSetFromOptions(eps,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSSolve(eps,ierr))

!     ** Optional: Get some information from the solver and display it
      PetscCallA(EPSGetType(eps,tname,ierr))
      if (rank==0) then
        write(*,120) tname
      end if
 120  format (' Solution method: ',A)
      PetscCallA(EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
      if (rank==0) then
        write(*,130) nev
      end if
 130  format (' Number of requested eigenvalues:',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** show detailed info unless -terse option is given by user
      PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-terse',terse,ierr))
      if (terse) then
        PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      else
        PetscCallA(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr))
        PetscCallA(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr))
      end if
      PetscCallA(EPSDestroy(eps,ierr))
      PetscCallA(MatDestroy(A,ierr))

      PetscCallA(SlepcFinalize(ierr))
      end program ex1f

!/*TEST
!
!   test:
!      args: -eps_nev 4 -terse
!      filter: sed -e "s/3.83791/3.83792/"
!
!TEST*/
