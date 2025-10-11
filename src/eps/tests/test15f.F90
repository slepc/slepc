!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test15f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Tests custom monitors from Fortran.
!
!  The command line options are:
!    -n <n>, where <n> = number of grid points = matrix size
!    -my_eps_monitor, activates the custom monitor
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepceps.h>
      use slepceps
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     A     operator matrix
!     eps   eigenproblem solver context

      Mat            A
      EPS            eps
      EPSType        tname
      PetscInt       n, i, Istart, Iend
      PetscInt       nev
      PetscInt       col(3)
      PetscInt       i1,i2,i3
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg
      PetscScalar    value(3)

!  Note: Any user-defined Fortran routines (such as MyEPSMonitor)
!  MUST be declared as external.

      external MyEPSMonitor

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 30
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))

      if (rank .eq. 0) then
        write(*,100) n
      end if
 100  format (/'1-D Laplacian Eigenproblem, n =',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the operator matrix that defines the eigensystem, Ax=kx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))

      i1 = 1
      i2 = 2
      i3 = 3
      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      if (Istart .eq. 0) then
        i = 0
        col(1) = 0
        col(2) = 1
        value(1) =  2.0
        value(2) = -1.0
        PetscCallA(MatSetValues(A,i1,[i],i2,col,value,INSERT_VALUES,ierr))
        Istart = Istart+1
      end if
      if (Iend .eq. n) then
        i = n-1
        col(1) = n-2
        col(2) = n-1
        value(1) = -1.0
        value(2) =  2.0
        PetscCallA(MatSetValues(A,i1,[i],i2,col,value,INSERT_VALUES,ierr))
        Iend = Iend-1
      end if
      value(1) = -1.0
      value(2) =  2.0
      value(3) = -1.0
      do i=Istart,Iend-1
        col(1) = i-1
        col(2) = i
        col(3) = i+1
        PetscCallA(MatSetValues(A,i1,[i],i3,col,value,INSERT_VALUES,ierr))
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

!     ** Set user-defined monitor
      PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-my_eps_monitor',flg,ierr))
      if (flg) then
        PetscCallA(EPSMonitorSet(eps,MyEPSMonitor,0,PETSC_NULL_FUNCTION,ierr))
      end if

!     ** Set solver parameters at runtime
      PetscCallA(EPSSetFromOptions(eps,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSSolve(eps,ierr))

!     ** Optional: Get some information from the solver and display it
      PetscCallA(EPSGetType(eps,tname,ierr))
      if (rank .eq. 0) then
        write(*,120) tname
      end if
 120  format (' Solution method: ',A)
      PetscCallA(EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
      if (rank .eq. 0) then
        write(*,130) nev
      end if
 130  format (' Number of requested eigenvalues:',I2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      PetscCallA(EPSDestroy(eps,ierr))
      PetscCallA(MatDestroy(A,ierr))

      PetscCallA(SlepcFinalize(ierr))
      end

! --------------------------------------------------------------
!
!  MyEPSMonitor - This is a user-defined routine for monitoring
!  the EPS iterative solvers.
!
!  Input Parameters:
!    eps   - eigensolver context
!    its   - iteration number
!    nconv - number of converged eigenpairs
!    eigr  - real part of the eigenvalues
!    eigi  - imaginary part of the eigenvalues
!    errest- relative error estimates for each eigenpair
!    nest  - number of error estimates
!    dummy - optional user-defined monitor context (unused here)
!
      subroutine MyEPSMonitor(eps,its,nconv,eigr,eigi,errest,nest,dummy,ierr)
#include <slepc/finclude/slepceps.h>
      use slepceps
      implicit none

      EPS            eps
      PetscErrorCode ierr
      PetscInt       its,nconv,nest,dummy
      PetscScalar    eigr(*),eigi(*)
      PetscReal      re,errest(*)
      PetscMPIInt    rank

      PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      if (its .gt. 0 .and. rank .eq. 0) then
        re = PetscRealPart(eigr(nconv+1))
        write(6,140) its,nconv,re,errest(nconv+1)
      end if

 140  format(i3,' EPS nconv=',i2,' first unconverged value (error) ',f7.4,' (',g10.3,')')
      ierr = 0
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -my_eps_monitor
!      requires: double
!
!TEST*/
