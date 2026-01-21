!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex6f [-help] [-m <m>] [all SLEPc options]
!
!  Description: Eigensystem from the Ising model for ferromagnetic materials.
!  Information about the model can be found at the following
!  site https://math.nist.gov/MatrixMarket/data/NEP
!
!  The command line options are:
!    -m <m>, where <m> is the number of 2x2 blocks, i.e. matrix size N=2*m
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepceps.h>

module ex6fmodule
  use slepceps
  implicit none

contains
  ! -------------------------------------------------------------------
  ! The actual routine for the matrix-vector product
  ! See https://math.nist.gov/MatrixMarket/data/NEP/mvmisg/mvmisg.html
  !
  ! Computes Y(:,1:M) = op(A)*X(:,1:M)
  ! where op(A) is A or A' and A is the Ising matrix.
  !
  !   trans   (input) integer
  !           If trans = 0, compute y(:,1:M) = A*x(:,1:M)
  !           If trans = 1, compute y(:,1:M) = A'*x(:,1:M)
  !
  !   n       (input) integer
  !           The order of the matrix, it has to be an even number.
  !
  !   m       (input) integeR
  !           The number of columns of x to multiply.
  !
  !   x       (input) double precision array, dimension (ldx, m)
  !           x contains the matrix (vectors) x.
  !
  !   ldx     (input) integer
  !           The leading dimension of array x, ldx >= max(1, n)
  !
  !   y       (output) double precision array, dimension (ldx, m)
  !           contains the product of the matrix op(A) with x.
  !
  !   ldy     (input) integer
  !           The leading dimension of array y, ldy >= max(1, n)
  !
  subroutine mvmisg(trans, n, m, x, ldx, y, ldy)
    use petscsys
    implicit none

    PetscInt    :: ldy, ldx, m, n, trans
    PetscScalar :: y(ldy, *), x(ldx, *)
    PetscInt    :: i, k
    PetscReal   :: alpha, beta, cosa, cosb, sina, sinb
    PetscScalar :: temp, temp1

    alpha = PETSC_PI/4
    beta = PETSC_PI/4
    cosa = cos(alpha)
    sina = sin(alpha)
    cosb = cos(beta)
    sinb = sin(beta)

    if (trans == 0) then

!     ** Compute y(:,1:m) = A*x(:,1:m)

      do k = 1, m
        y(1, k) = cosb*x(1, k) - sinb*x(n, k)
        do i = 2, n - 1, 2
          y(i, k) = cosb*x(i, k) + sinb*x(i + 1, k)
          y(i + 1, k) = -sinb*x(i, k) + cosb*x(i + 1, k)
        end do
        y(n, k) = sinb*x(1, k) + cosb*x(n, k)
        do i = 1, n, 2
          temp = cosa*y(i, k) + sina*y(i + 1, k)
          y(i + 1, k) = -sina*y(i, k) + cosa*y(i + 1, k)
          y(i, k) = temp
        end do
      end do

    else if (trans == 1) then

!     ** Compute y(:1:m) = A'*x(:,1:m)

      do k = 1, m
        do i = 1, n, 2
          y(i, k) = cosa*x(i, k) - sina*x(i + 1, k)
          y(i + 1, k) = sina*x(i, k) + cosa*x(i + 1, k)
        end do
        temp = cosb*y(1, k) + sinb*y(n, k)
        do i = 2, n - 1, 2
          temp1 = cosb*y(i, k) - sinb*y(i + 1, k)
          y(i + 1, k) = sinb*y(i, k) + cosb*y(i + 1, k)
          y(i, k) = temp1
        end do
        y(n, k) = -sinb*y(1, k) + cosb*y(n, k)
        y(1, k) = temp
      end do

    end if
  end subroutine

  ! -------------------------------------------------------------------
  ! MatMult_Ising - user provided matrix-vector multiply
  !
  ! Input Parameters:
  !   A - matrix
  !   x - input vector
  !
  ! Output Parameter:
  !   y - output vector
  !
  subroutine MatMult_Ising(A, x, y, ierr)
    use petscmat
    implicit none

    Mat                  :: A
    Vec                  :: x, y
    PetscInt             :: trans, N
    PetscScalar, pointer :: xx(:), yy(:)
    PetscErrorCode, intent(out) :: ierr

    PetscCall(MatGetSize(A, N, PETSC_NULL_INTEGER, ierr))
    PetscCall(VecGetArrayRead(x, xx, ierr))
    PetscCall(VecGetArray(y, yy, ierr))

    trans = 0
    call mvmisg(trans, N, 1_PETSC_INT_KIND, xx, N, yy, N)

    PetscCall(VecRestoreArrayRead(x, xx, ierr))
    PetscCall(VecRestoreArray(y, yy, ierr))
  end subroutine

  ! -------------------------------------------------------------------
  ! MatMultTranspose_Ising - user provided transpose matrix-vector multiply
  !
  ! Input Parameters:
  !   A - matrix
  !   x - input vector
  !
  ! Output Parameter:
  !   y - output vector
  !
  subroutine MatMultTranspose_Ising(A, x, y, ierr)
    use petscmat
    implicit none

    Mat                  :: A
    Vec                  :: x, y
    PetscInt             :: trans, N
    PetscScalar, pointer :: xx(:), yy(:)
    PetscErrorCode, intent(out) :: ierr

    PetscCall(MatGetSize(A, N, PETSC_NULL_INTEGER, ierr))
    PetscCall(VecGetArrayRead(x, xx, ierr))
    PetscCall(VecGetArray(y, yy, ierr))

    trans = 1
    call mvmisg(trans, N, 1_PETSC_INT_KIND, xx, N, yy, N)

    PetscCall(VecRestoreArrayRead(x, xx, ierr))
    PetscCall(VecRestoreArray(y, yy, ierr))
  end subroutine

end module ex6fmodule

program ex6f
  use slepceps
  use ex6fmodule
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Mat            :: A    ! operator matrix
  EPS            :: eps  ! aeigenproblem solver context
  EPSType        :: tname
  PetscReal      :: tol
  PetscInt       :: N, m
  PetscInt       :: nev, maxit, its
  PetscMPIInt    :: sz, rank
  PetscErrorCode :: ierr
  PetscBool      :: flg, terse

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER, ierr))
#if defined(PETSC_USE_COMPLEX)
  SETERRA(PETSC_COMM_SELF, PETSC_ERR_SUP, 'This example requires real numbers')
#endif
  PetscCallMPIA(MPI_Comm_size(PETSC_COMM_WORLD, sz, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
  PetscCheckA(sz == 1, PETSC_COMM_SELF, PETSC_ERR_WRONG_MPI_SIZE, 'This is a uniprocessor example only!')
  m = 30
  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-m', m, flg, ierr))
  N = 2*m

  if (rank == 0) then
    write (*, '(/a,i6,a/)') 'Ising Model Eigenproblem, m=', m, ', (N=2*m)'
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Register the matrix-vector subroutine for the operator that defines
!  the eigensystem, Ax=kx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(MatCreateShell(PETSC_COMM_WORLD, N, N, N, N, PETSC_NULL_INTEGER, A, ierr))
  PetscCallA(MatShellSetOperation(A, MATOP_MULT, MatMult_Ising, ierr))
  PetscCallA(MatShellSetOperation(A, MATOP_MULT_TRANSPOSE, MatMultTranspose_Ising, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create the eigensolver and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ** Create eigensolver context
  PetscCallA(EPSCreate(PETSC_COMM_WORLD, eps, ierr))

! ** Set operators. In this case, it is a standard eigenvalue problem
  PetscCallA(EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr))
  PetscCallA(EPSSetProblemType(eps, EPS_NHEP, ierr))

! ** Set solver parameters at runtime
  PetscCallA(EPSSetFromOptions(eps, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(EPSSolve(eps, ierr))
  PetscCallA(EPSGetIterationNumber(eps, its, ierr))
  if (rank == 0) then
    write (*, '(a,i4)') ' Number of iterations of the method: ', its
  end if

! ** Optional: Get some information from the solver and display it
  PetscCallA(EPSGetType(eps, tname, ierr))
  if (rank == 0) then
    write (*, '(a,a)') ' Solution method: ', tname
  end if
  PetscCallA(EPSGetDimensions(eps, nev, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Number of requested eigenvalues:', nev
  end if
  PetscCallA(EPSGetTolerances(eps, tol, maxit, ierr))
  if (rank == 0) then
    write (*, '(a,1pe11.4,a,i6)') ' Stopping condition: tol=', tol, ', maxit=', maxit
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
  PetscCallA(EPSDestroy(eps, ierr))
  PetscCallA(MatDestroy(A, ierr))
  PetscCallA(SlepcFinalize(ierr))
end program ex6f

!/*TEST
!
!   testset:
!      args: -eps_max_it 1000 -eps_ncv 12 -eps_tol 1e-5 -eps_nev 4 -eps_largest_imaginary -terse
!      requires: !complex
!      output_file: output/ex6f_1.out
!      filter: grep -v iterations | sed -e 's/-0.00000/0.00000/g'
!      test:
!         suffix: 1
!      test:
!         suffix: 1_ts
!         args: -eps_two_sided
!         requires: !single
!
!TEST*/
