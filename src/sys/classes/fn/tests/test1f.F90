!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test1f [-help]
!
!  Description: Test rational function in Fortran.
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepcfn.h>
program test1f
  use slepcfn
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FN             :: fn
  PetscInt       :: i, n, na, nb
  PetscMPIInt    :: rank
  PetscErrorCode :: ierr
  PetscScalar    :: x, y, yp, p(10), q(10), five
  PetscScalar    :: pp(10), qq(10), tau, eta

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
  PetscCallA(FNCreate(PETSC_COMM_WORLD, fn, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Polynomial p(x)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  na = 5
  p(1) = -3.1
  p(2) = 1.1
  p(3) = 1.0
  p(4) = -2.0
  p(5) = 3.5
  PetscCallA(FNSetType(fn, FNRATIONAL, ierr))
  PetscCallA(FNRationalSetNumerator(fn, na, p, ierr))
  PetscCallA(FNView(fn, PETSC_NULL_VIEWER, ierr))
  x = 2.2
  PetscCallA(FNEvaluateFunction(fn, x, y, ierr))
  PetscCallA(FNEvaluateDerivative(fn, x, yp, ierr))
  call PrintInfo(x, y, yp)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Inverse of polynomial 1/q(x)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  na = 0
  nb = 3
  q(1) = -3.1
  q(2) = 1.1
  q(3) = 1.0
  PetscCallA(FNSetType(fn, FNRATIONAL, ierr))
  PetscCallA(FNRationalSetNumerator(fn, na, PETSC_NULL_SCALAR_ARRAY, ierr))
  PetscCallA(FNRationalSetDenominator(fn, nb, q, ierr))
  PetscCallA(FNView(fn, PETSC_NULL_VIEWER, ierr))
  x = 2.2
  PetscCallA(FNEvaluateFunction(fn, x, y, ierr))
  PetscCallA(FNEvaluateDerivative(fn, x, yp, ierr))
  call PrintInfo(x, y, yp)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Rational p(x)/q(x)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  na = 2
  nb = 3
  p(1) = 1.1
  p(2) = 1.1
  q(1) = 1.0
  q(2) = -2.0
  q(3) = 3.5
  PetscCallA(FNSetType(fn, FNRATIONAL, ierr))
  PetscCallA(FNRationalSetNumerator(fn, na, p, ierr))
  PetscCallA(FNRationalSetDenominator(fn, nb, q, ierr))
  tau = 1.2
  eta = 0.5
  PetscCallA(FNSetScale(fn, tau, eta, ierr))
  PetscCallA(FNView(fn, PETSC_NULL_VIEWER, ierr))
  x = 2.2
  PetscCallA(FNEvaluateFunction(fn, x, y, ierr))
  PetscCallA(FNEvaluateDerivative(fn, x, yp, ierr))
  call PrintInfo(x, y, yp)

  PetscCallA(FNRationalGetNumerator(fn, n, pp, ierr))
  if (rank == 0) then
    write (*, 100) 'Numerator', (PetscRealPart(pp(i)), i=1, n)
  end if
  PetscCallA(FNRationalGetDenominator(fn, n, qq, ierr))
  if (rank == 0) then
    write (*, 100) 'Denominator', (PetscRealPart(qq(i)), i=1, n)
  end if
100 format(A15, 10F6.1)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Constant
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  na = 1
  nb = 0
  five = 5.0
  PetscCallA(FNSetType(fn, FNRATIONAL, ierr))
  PetscCallA(FNRationalSetNumerator(fn, na, [five], ierr))
  PetscCallA(FNRationalSetDenominator(fn, nb, PETSC_NULL_SCALAR_ARRAY, ierr))
  PetscCallA(FNView(fn, PETSC_NULL_VIEWER, ierr))
  x = 2.2
  PetscCallA(FNEvaluateFunction(fn, x, y, ierr))
  PetscCallA(FNEvaluateDerivative(fn, x, yp, ierr))
  call PrintInfo(x, y, yp)

! *** Clean up
  PetscCallA(FNDestroy(fn, ierr))
  PetscCallA(SlepcFinalize(ierr))

contains

  subroutine PrintInfo(x, y, yp)
    use slepcfn
    implicit none
    PetscScalar    :: x, y, yp
    PetscReal      :: re, im
    PetscMPIInt    :: rank
    PetscErrorCode :: ierr

    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
    if (rank == 0) then
      re = PetscRealPart(y)
      im = PetscImaginaryPart(y)
      if (abs(im) < 1.d-10) then
        write (*, 110) 'f', PetscRealPart(x), re
      else
        write (*, 120) 'f', PetscRealPart(x), re, im
      end if
      re = PetscRealPart(yp)
      im = PetscImaginaryPart(yp)
      if (abs(im) < 1.d-10) then
        write (*, 110) 'f''', PetscRealPart(x), re
      else
        write (*, 120) 'f''', PetscRealPart(x), re, im
      end if
    end if
110 format(A2, '(', F4.1, ') = ', F10.5)
120 format(A2, '(', F4.1, ') = ', F10.5, SP, F9.5, 'i')

  end subroutine

end program test1f

!/*TEST
!
!   test:
!      suffix: 1
!      nsize: 1
!      requires: !single
!
!TEST*/
