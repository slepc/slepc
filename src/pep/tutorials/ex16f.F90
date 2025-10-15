!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex16f [-help] [-n <n>] [-m <m>] [SLEPc opts]
!
!  Description: Simple example that solves a quadratic eigensystem with PEP.
!  This is the Fortran90 equivalent to ex16.c
!
!  The command line options are:
!    -n <n>, where <n> = number of grid subdivisions in x dimension
!    -m <m>, where <m> = number of grid subdivisions in y dimension
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepcpep.h>
program ex16f
  use slepcpep
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Mat            :: M, C, K, A(3)   ! problem matrices
  PEP            :: pep             ! polynomial eigenproblem solver context
  PEPType        :: tname
  PetscInt       :: N, nx, ny, i, j, Istart, Iend, II
  PetscInt       :: nev, ithree
  PetscMPIInt    :: rank
  PetscErrorCode :: ierr
  PetscBool      :: flg, terse
  PetscScalar    :: mone, two, four, val

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
  nx = 10
  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-n', nx, flg, ierr))
  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-m', ny, flg, ierr))
  if (.not. flg) then
    ny = nx
  end if
  N = nx*ny
  if (rank == 0) then
    write (*, 100) N, nx, ny
  end if
100 format(/'Quadratic Eigenproblem, N=', I6, ' (', I4, 'x', I4, ' grid)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ** K is the 2-D Laplacian
  PetscCallA(MatCreate(PETSC_COMM_WORLD, K, ierr))
  PetscCallA(MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, N, N, ierr))
  PetscCallA(MatSetFromOptions(K, ierr))
  PetscCallA(MatGetOwnershipRange(K, Istart, Iend, ierr))
  mone = -1.0
  four = 4.0
  do II = Istart, Iend - 1
    i = II/nx
    j = II - i*nx
    if (i > 0) then
      PetscCallA(MatSetValue(K, II, II - nx, mone, INSERT_VALUES, ierr))
    end if
    if (i < ny - 1) then
      PetscCallA(MatSetValue(K, II, II + nx, mone, INSERT_VALUES, ierr))
    end if
    if (j > 0) then
      PetscCallA(MatSetValue(K, II, II - 1, mone, INSERT_VALUES, ierr))
    end if
    if (j < nx - 1) then
      PetscCallA(MatSetValue(K, II, II + 1, mone, INSERT_VALUES, ierr))
    end if
    PetscCallA(MatSetValue(K, II, II, four, INSERT_VALUES, ierr))
  end do
  PetscCallA(MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY, ierr))

! ** C is the 1-D Laplacian on horizontal lines
  PetscCallA(MatCreate(PETSC_COMM_WORLD, C, ierr))
  PetscCallA(MatSetSizes(C, PETSC_DECIDE, PETSC_DECIDE, N, N, ierr))
  PetscCallA(MatSetFromOptions(C, ierr))
  PetscCallA(MatGetOwnershipRange(C, Istart, Iend, ierr))
  two = 2.0
  do II = Istart, Iend - 1
    i = II/nx
    j = II - i*nx
    if (j > 0) then
      PetscCallA(MatSetValue(C, II, II - 1, mone, INSERT_VALUES, ierr))
    end if
    if (j < nx - 1) then
      PetscCallA(MatSetValue(C, II, II + 1, mone, INSERT_VALUES, ierr))
    end if
    PetscCallA(MatSetValue(C, II, II, two, INSERT_VALUES, ierr))
  end do
  PetscCallA(MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY, ierr))

! ** M is a diagonal matrix
  PetscCallA(MatCreate(PETSC_COMM_WORLD, M, ierr))
  PetscCallA(MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, N, N, ierr))
  PetscCallA(MatSetFromOptions(M, ierr))
  PetscCallA(MatGetOwnershipRange(M, Istart, Iend, ierr))
  do II = Istart, Iend - 1
    val = II + 1
    PetscCallA(MatSetValue(M, II, II, val, INSERT_VALUES, ierr))
  end do
  PetscCallA(MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ** Create eigensolver context
  PetscCallA(PEPCreate(PETSC_COMM_WORLD, pep, ierr))

! ** Set matrices and problem type
  A(1) = K
  A(2) = C
  A(3) = M
  ithree = 3
  PetscCallA(PEPSetOperators(pep, ithree, A, ierr))
  PetscCallA(PEPSetProblemType(pep, PEP_GENERAL, ierr))

! ** Set solver parameters at runtime
  PetscCallA(PEPSetFromOptions(pep, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(PEPSolve(pep, ierr))

! ** Optional: Get some information from the solver and display it
  PetscCallA(PEPGetType(pep, tname, ierr))
  if (rank == 0) then
    write (*, 120) tname
  end if
120 format(' Solution method: ', A)
  PetscCallA(PEPGetDimensions(pep, nev, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr))
  if (rank == 0) then
    write (*, 130) nev
  end if
130 format(' Number of requested eigenvalues:', I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ** show detailed info unless -terse option is given by user
  PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-terse', terse, ierr))
  if (terse) then
    PetscCallA(PEPErrorView(pep, PEP_ERROR_BACKWARD, PETSC_NULL_VIEWER, ierr))
  else
    PetscCallA(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL, ierr))
    PetscCallA(PEPConvergedReasonView(pep, PETSC_VIEWER_STDOUT_WORLD, ierr))
    PetscCallA(PEPErrorView(pep, PEP_ERROR_BACKWARD, PETSC_VIEWER_STDOUT_WORLD, ierr))
    PetscCallA(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD, ierr))
  end if
  PetscCallA(PEPDestroy(pep, ierr))
  PetscCallA(MatDestroy(K, ierr))
  PetscCallA(MatDestroy(C, ierr))
  PetscCallA(MatDestroy(M, ierr))
  PetscCallA(SlepcFinalize(ierr))
end program ex16f

!/*TEST
!
!   test:
!      suffix: 1
!      args: -pep_nev 4 -pep_ncv 19 -terse
!      requires: !complex
!
!TEST*/
