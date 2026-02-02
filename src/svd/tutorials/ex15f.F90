!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex15f [-help] [-n <n>] [-mu <mu>] [all SLEPc options]
!
!  Description: Singular value decomposition of the Lauchli matrix.
!
!  The command line options are:
!    -n <n>, where <n> = matrix dimension.
!    -mu <mu>, where <mu> = subdiagonal value.
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepcsvd.h>
program ex15f
  use slepcsvd
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Mat            :: A    ! operator matrix
  SVD            :: svd  ! singular value solver context
  SVDType        :: tname
  PetscReal      :: tol, error, sigma, mu
  PetscInt       :: n, i, j, Istart, Iend
  PetscInt       :: nsv, maxit, its, nconv
  PetscMPIInt    :: rank
  PetscErrorCode :: ierr
  PetscBool      :: flg
  PetscScalar    :: alpha
  PetscScalar, parameter :: one = 1.0

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
  n = 100
  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-n', n, flg, ierr))
  mu = PETSC_SQRT_MACHINE_EPSILON
  PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-mu', mu, flg, ierr))

  if (rank == 0) then
    write (*, '(/a,i3,a,e12.4,a)') 'Lauchli SVD, n =', n, ', mu=', mu, ' (Fortran)'
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Build the Lauchli matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(MatCreate(PETSC_COMM_WORLD, A, ierr))
  PetscCallA(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n + 1, n, ierr))
  PetscCallA(MatSetFromOptions(A, ierr))

  PetscCallA(MatGetOwnershipRange(A, Istart, Iend, ierr))
  do i = Istart, Iend - 1
    if (i == 0) then
      do j = 0, n - 1
        PetscCallA(MatSetValue(A, i, j, one, INSERT_VALUES, ierr))
      end do
    else
      alpha = mu
      PetscCallA(MatSetValue(A, i, i - 1, alpha, INSERT_VALUES, ierr))
    end if
  end do

  PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create the singular value solver and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ** Create singular value solver context
  PetscCallA(SVDCreate(PETSC_COMM_WORLD, svd, ierr))

! ** Set operators and problem type
  PetscCallA(SVDSetOperators(svd, A, PETSC_NULL_MAT, ierr))
  PetscCallA(SVDSetProblemType(svd, SVD_STANDARD, ierr))

! ** Use thick-restart Lanczos as default solver
  PetscCallA(SVDSetType(svd, SVDTRLANCZOS, ierr))

! ** Set solver parameters at runtime
  PetscCallA(SVDSetFromOptions(svd, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Solve the singular value system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SVDSolve(svd, ierr))
  PetscCallA(SVDGetIterationNumber(svd, its, ierr))
  if (rank == 0) then
    write (*, '(/a,i4)') ' Number of iterations of the method:', its
  end if

! ** Optional: Get some information from the solver and display it
  PetscCallA(SVDGetType(svd, tname, ierr))
  if (rank == 0) then
    write (*, '(a,a)') ' Solution method: ', tname
  end if
  PetscCallA(SVDGetDimensions(svd, nsv, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Number of requested singular values:', nsv
  end if
  PetscCallA(SVDGetTolerances(svd, tol, maxit, ierr))
  if (rank == 0) then
    write (*, '(a,1pe11.4,a,i4)') ' Stopping condition: tol=', tol, ', maxit=', maxit
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ** Get number of converged singular triplets
  PetscCallA(SVDGetConverged(svd, nconv, ierr))
  if (rank == 0) then
    write (*, '(a,i2/)') ' Number of converged approximate singular triplets:', nconv
  end if

! ** Display singular values and relative errors
  if (nconv > 0) then
    if (rank == 0) then
      write (*, *) '       sigma          relative error'
      write (*, *) ' ----------------- ------------------'
    end if
    do i = 0, nconv - 1
!     ** Get i-th singular value
      PetscCallA(SVDGetSingularTriplet(svd, i, sigma, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr))

!     ** Compute the relative error for each singular triplet
      PetscCallA(SVDComputeError(svd, i, SVD_ERROR_RELATIVE, error, ierr))
      if (rank == 0) then
        write (*, '(1p,a,e12.4,a,e12.4)') '   ', sigma, '       ', error
      end if

    end do
    if (rank == 0) then
      write (*, *)
    end if
  end if

! ** Free work space
  PetscCallA(SVDDestroy(svd, ierr))
  PetscCallA(MatDestroy(A, ierr))

  PetscCallA(SlepcFinalize(ierr))
end program ex15f

!/*TEST
!
!   test:
!      suffix: 1
!      filter: sed -e "s/[0-9]\.[0-9]*E[+-]\([0-9]*\)/removed/g"
!
!TEST*/
