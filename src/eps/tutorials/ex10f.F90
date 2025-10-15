!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex10f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Illustrates the use of shell spectral transformations.
!  The problem to be solved is the same as ex1.c and corresponds to the
!  Laplacian operator in 1 dimension
!
!  The command line options are:
!    nm <n>, where <n> is the number of grid subdivisions = matrix dimension
! ----------------------------------------------------------------------

! Module contains data needed by shell ST
!
#include <slepc/finclude/slepceps.h>
module ex10fmodule
  use slepceps
  implicit none

  KSP myksp

contains
  ! -------------------------------------------------------------------
  ! STApply_User - This routine demonstrates the use of a user-provided spectral
  ! transformation. The transformation implemented in this code is just OP=A^-1.
  !
  ! Input Parameters:
  !   st - spectral transformation context
  !   x - input vector
  !
  ! Output Parameter:
  !   y - output vector
  !
  subroutine STApply_User(st, x, y, ierr)
    use slepceps
    implicit none

    ST             :: st
    Vec            :: x, y
    PetscErrorCode :: ierr

    PetscCall(KSPSolve(myksp, x, y, ierr))
  end subroutine

  ! -------------------------------------------------------------------
  ! STApplyTranspose_User - This is not required unless using a two-sided eigensolver
  !
  ! Input Parameters:
  !   st - spectral transformation context
  !   x - input vector
  !
  ! Output Parameter:
  !   y - output vector
  !
  subroutine STApplyTranspose_User(st, x, y, ierr)
    use slepceps
    implicit none

    ST             :: st
    Vec            :: x, y
    PetscErrorCode :: ierr

    PetscCall(KSPSolveTranspose(myksp, x, y, ierr))
  end subroutine

#if defined(PETSC_USE_COMPLEX)
  ! -------------------------------------------------------------------
  ! STApplyHermitianTranspose_User - This is not required unless using a two-sided eigensolver
  ! in complex scalars
  !
  ! Input Parameters:
  !   st - spectral transformation context
  !   x - input vector
  !
  ! Output Parameter:
  !   y - output vector
  !
  subroutine STApplyHermitianTranspose_User(st, x, y, ierr)
    use slepceps
    implicit none

    ST             :: st
    Vec            :: x, y, w
    PetscErrorCode :: ierr

    PetscCall(VecDuplicate(x, w, ierr))
    PetscCall(VecCopy(x, w, ierr))
    PetscCall(VecConjugate(w, ierr))
    PetscCall(KSPSolveTranspose(myksp, w, y, ierr))
    PetscCall(VecConjugate(y, ierr))
    PetscCall(VecDestroy(w, ierr))
  end subroutine
#endif

  ! -------------------------------------------------------------------
  ! STBackTransform_User - This routine demonstrates the use of a user-provided spectral
  ! transformation
  !
  ! Input Parameters:
  !   st - spectral transformation context
  !   n  - number of eigenvalues to transform
  !
  ! Output Parameters:
  !   eigr - real part of eigenvalues
  !   eigi - imaginary part of eigenvalues
  !
  subroutine STBackTransform_User(st, n, eigr, eigi, ierr)
    use slepceps
    implicit none

    ST             :: st
    PetscInt       :: n, j
    PetscScalar    :: eigr(*), eigi(*)
    PetscErrorCode :: ierr

    do j = 1, n
      eigr(j) = 1.0/eigr(j)
    end do
    ierr = 0
  end subroutine

end module ex10fmodule

! ----------------------------------------------------------------------

program ex10f
  use slepceps
  use ex10fmodule
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  Mat            :: A    ! operator matrix
  EPS            :: eps  ! eigenproblem solver context
  ST             :: st
  EPSType        :: tname
  PetscInt       :: n, i, Istart, Iend, one, two, three
  PetscInt       :: nev, row(1), col(3)
  PetscScalar    :: val(3)
  PetscBool      :: flg, isShell, terse
  PetscMPIInt    :: rank
  PetscErrorCode :: ierr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  one = 1
  two = 2
  three = 3
  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
  n = 30
  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-n', n, flg, ierr))

  if (rank == 0) then
    write (*, '(/A,I6/)') '1-D Laplacian Eigenproblem (shell-enabled), n=', n
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute the operator matrix that defines the eigensystem, Ax=kx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(MatCreate(PETSC_COMM_WORLD, A, ierr))
  PetscCallA(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr))
  PetscCallA(MatSetFromOptions(A, ierr))

  PetscCallA(MatGetOwnershipRange(A, Istart, Iend, ierr))
  if (Istart == 0) then
    row(1) = 0
    col(1) = 0
    col(2) = 1
    val(1) = 2.0
    val(2) = -1.0
    PetscCallA(MatSetValues(A, one, row, two, col, val, INSERT_VALUES, ierr))
    Istart = Istart + 1
  end if
  if (Iend == n) then
    row(1) = n - 1
    col(1) = n - 2
    col(2) = n - 1
    val(1) = -1.0
    val(2) = 2.0
    PetscCallA(MatSetValues(A, one, row, two, col, val, INSERT_VALUES, ierr))
    Iend = Iend - 1
  end if
  val(1) = -1.0
  val(2) = 2.0
  val(3) = -1.0
  do i = Istart, Iend - 1
    row(1) = i
    col(1) = i - 1
    col(2) = i
    col(3) = i + 1
    PetscCallA(MatSetValues(A, one, row, three, col, val, INSERT_VALUES, ierr))
  end do

  PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ** Create eigensolver context
  PetscCallA(EPSCreate(PETSC_COMM_WORLD, eps, ierr))

! ** Set operators. In this case, it is a standard eigenvalue problem
  PetscCallA(EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr))
  PetscCallA(EPSSetProblemType(eps, EPS_NHEP, ierr))

! ** Set solver parameters at runtime
  PetscCallA(EPSSetFromOptions(eps, ierr))

! ** Initialize shell spectral transformation if selected by user
  PetscCallA(EPSGetST(eps, st, ierr))
  PetscCallA(PetscObjectTypeCompare(st, STSHELL, isShell, ierr))

  if (isShell) then
!   ** Change sorting criterion since this ST example computes values
!   ** closest to 0
    PetscCallA(EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL, ierr))

!   ** In Fortran, instead of a context for the user-defined spectral transform
!   ** we use a module containing any application-specific data, initialized here
    PetscCallA(KSPCreate(PETSC_COMM_WORLD, myksp, ierr))
    PetscCallA(KSPAppendOptionsPrefix(myksp, "st_", ierr))

!   ** (Required) Set the user-defined routine for applying the operator
    PetscCallA(STShellSetApply(st, STApply_User, ierr))

!   ** (Optional) Set the user-defined routine for applying the transposed operator
    PetscCallA(STShellSetApplyTranspose(st, STApplyTranspose_User, ierr))

#if defined(PETSC_USE_COMPLEX)
!   ** (Optional) Set the user-defined routine for applying the conjugate-transposed operator
    PetscCallA(STShellSetApplyHermitianTranspose(st, STApplyHermitianTranspose_User, ierr))
#endif

!   ** (Optional) Set the user-defined routine for back-transformation
    PetscCallA(STShellSetBackTransform(st, STBackTransform_User, ierr))

!   ** (Optional) Set a name for the transformation, used for STView()
    PetscCallA(PetscObjectSetName(st, 'MyTransformation', ierr))

!   ** (Optional) Do any setup required for the new transformation
    PetscCallA(KSPSetOperators(myksp, A, A, ierr))
    PetscCallA(KSPSetFromOptions(myksp, ierr))
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(EPSSolve(eps, ierr))

! ** Optional: Get some information from the solver and display it
  PetscCallA(EPSGetType(eps, tname, ierr))
  if (rank == 0) then
    write (*, '(A,A,/)') ' Solution method: ', tname
  end if
  PetscCallA(EPSGetDimensions(eps, nev, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr))
  if (rank == 0) then
    write (*, '(A,I2)') ' Number of requested eigenvalues:', nev
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ** show detailed info unless -terse option is given by user
  PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-terse', terse, ierr))
  if (terse) then
    PetscCallA(EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_NULL_VIEWER, ierr))
  else
    PetscCallA(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL, ierr))
    PetscCallA(EPSConvergedReasonView(eps, PETSC_VIEWER_STDOUT_WORLD, ierr))
    PetscCallA(EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD, ierr))
    PetscCallA(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD, ierr))
  end if
  if (isShell) then
    PetscCallA(KSPDestroy(myksp, ierr))
  end if
  PetscCallA(EPSDestroy(eps, ierr))
  PetscCallA(MatDestroy(A, ierr))
  PetscCallA(SlepcFinalize(ierr))
end program ex10f

!/*TEST
!
!   testset:
!      args: -eps_nev 5 -eps_non_hermitian -terse
!      output_file: output/ex10_1.out
!      requires: !single
!      test:
!         suffix: 1_sinvert
!         args: -st_type sinvert
!      test:
!         suffix: 1_shell
!         args: -st_type shell -eps_two_sided {{0 1}}
!
!TEST*/
