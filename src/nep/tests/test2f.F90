!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Simple example to test the NEP Fortran interface.
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepcnep.h>
program test2f
  use slepcnep
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Mat                  :: A(3), B
  FN                   :: f(3), g
  NEP                  :: nep
  DS                   :: ds
  RG                   :: rg
  PetscReal            :: tol
  PetscScalar          :: coeffs(2), tget, val
  PetscInt             :: n, i, its, Istart, Iend
  PetscInt             :: nev, ncv, mpd, nterm
  PetscInt             :: nc, np
  NEPWhich             :: which
  NEPConvergedReason   :: reason
  NEPType              :: tname
  NEPRefine            :: refine
  NEPRefineScheme      :: rscheme
  NEPConv              :: conv
  NEPStop              :: stp
  NEPProblemType       :: ptype
  MatStructure         :: mstr
  PetscMPIInt          :: rank
  PetscErrorCode       :: ierr
  PetscViewerAndFormat :: vf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
  n = 20
  if (rank == 0) then
    write (*, '(/a,i3,a)') 'Diagonal Nonlinear Eigenproblem, n =', n, ' (Fortran)'
  end if

! Matrices
  PetscCallA(MatCreate(PETSC_COMM_WORLD, A(1), ierr))
  PetscCallA(MatSetSizes(A(1), PETSC_DECIDE, PETSC_DECIDE, n, n, ierr))
  PetscCallA(MatSetFromOptions(A(1), ierr))
  PetscCallA(MatGetOwnershipRange(A(1), Istart, Iend, ierr))
  do i = Istart, Iend - 1
    val = i + 1
    PetscCallA(MatSetValue(A(1), i, i, val, INSERT_VALUES, ierr))
  end do
  PetscCallA(MatAssemblyBegin(A(1), MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A(1), MAT_FINAL_ASSEMBLY, ierr))

  PetscCallA(MatCreate(PETSC_COMM_WORLD, A(2), ierr))
  PetscCallA(MatSetSizes(A(2), PETSC_DECIDE, PETSC_DECIDE, n, n, ierr))
  PetscCallA(MatSetFromOptions(A(2), ierr))
  PetscCallA(MatGetOwnershipRange(A(2), Istart, Iend, ierr))
  do i = Istart, Iend - 1
    val = 1
    PetscCallA(MatSetValue(A(2), i, i, val, INSERT_VALUES, ierr))
  end do
  PetscCallA(MatAssemblyBegin(A(2), MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A(2), MAT_FINAL_ASSEMBLY, ierr))

  PetscCallA(MatCreate(PETSC_COMM_WORLD, A(3), ierr))
  PetscCallA(MatSetSizes(A(3), PETSC_DECIDE, PETSC_DECIDE, n, n, ierr))
  PetscCallA(MatSetFromOptions(A(3), ierr))
  PetscCallA(MatGetOwnershipRange(A(3), Istart, Iend, ierr))
  do i = Istart, Iend - 1
    val = real(n)/real(i + 1)
    PetscCallA(MatSetValue(A(3), i, i, val, INSERT_VALUES, ierr))
  end do
  PetscCallA(MatAssemblyBegin(A(3), MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A(3), MAT_FINAL_ASSEMBLY, ierr))

! Functions: f0=-lambda, f1=1.0, f2=sqrt(lambda)
  PetscCallA(FNCreate(PETSC_COMM_WORLD, f(1), ierr))
  PetscCallA(FNSetType(f(1), FNRATIONAL, ierr))
  nc = 2
  coeffs(1) = -1.0
  coeffs(2) = 0.0
  PetscCallA(FNRationalSetNumerator(f(1), nc, coeffs, ierr))

  PetscCallA(FNCreate(PETSC_COMM_WORLD, f(2), ierr))
  PetscCallA(FNSetType(f(2), FNRATIONAL, ierr))
  nc = 1
  coeffs(1) = 1.0
  PetscCallA(FNRationalSetNumerator(f(2), nc, coeffs, ierr))

  PetscCallA(FNCreate(PETSC_COMM_WORLD, f(3), ierr))
  PetscCallA(FNSetType(f(3), FNSQRT, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create eigensolver and test interface functions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(NEPCreate(PETSC_COMM_WORLD, nep, ierr))
  nterm = 3
  mstr = SAME_NONZERO_PATTERN
  PetscCallA(NEPSetSplitOperator(nep, nterm, A, f, mstr, ierr))
  PetscCallA(NEPGetSplitOperatorInfo(nep, nterm, mstr, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a)') ' Nonlinear function with ', nterm, ' terms'
  end if
  i = 0
  PetscCallA(NEPGetSplitOperatorTerm(nep, i, B, g, ierr))
  PetscCallA(MatView(B, PETSC_NULL_VIEWER, ierr))
  PetscCallA(FNView(g, PETSC_NULL_VIEWER, ierr))

  PetscCallA(NEPSetType(nep, NEPRII, ierr))
  PetscCallA(NEPGetType(nep, tname, ierr))
  if (rank == 0) then
    write (*, '(a,a)') ' Type set to ', tname
  end if

  PetscCallA(NEPGetProblemType(nep, ptype, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Problem type before changing = ', ptype
  end if
  PetscCallA(NEPSetProblemType(nep, NEP_RATIONAL, ierr))
  PetscCallA(NEPGetProblemType(nep, ptype, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' ... changed to ', ptype
  end if

  np = 1
  tol = 1e-9
  its = 2
  refine = NEP_REFINE_SIMPLE
  rscheme = NEP_REFINE_SCHEME_EXPLICIT
  PetscCallA(NEPSetRefine(nep, refine, np, tol, its, rscheme, ierr))
  PetscCallA(NEPGetRefine(nep, refine, np, tol, its, rscheme, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,f12.9,a,i2,a,i2)') ' Refinement: ', refine, ', tol=', tol, ', its=', its, ', scheme=', rscheme
  end if

  tget = 1.1
  PetscCallA(NEPSetTarget(nep, tget, ierr))
  PetscCallA(NEPGetTarget(nep, tget, ierr))
  PetscCallA(NEPSetWhichEigenpairs(nep, NEP_TARGET_MAGNITUDE, ierr))
  PetscCallA(NEPGetWhichEigenpairs(nep, which, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,f4.1)') ' Which = ', which, ', target = ', PetscRealPart(tget)
  end if

  nev = 1
  ncv = 12
  PetscCallA(NEPSetDimensions(nep, nev, ncv, PETSC_DETERMINE_INTEGER, ierr))
  PetscCallA(NEPGetDimensions(nep, nev, ncv, mpd, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,i2,a,i2)') ' Dimensions: nev=', nev, ', ncv=', ncv, ', mpd=', mpd
  end if

  tol = 1.0e-6
  its = 200
  PetscCallA(NEPSetTolerances(nep, tol, its, ierr))
  PetscCallA(NEPGetTolerances(nep, tol, its, ierr))
  if (rank == 0) then
    write (*, '(a,f9.6,a,i4)') ' Tolerance =', tol, ', max_its =', its
  end if

  PetscCallA(NEPSetConvergenceTest(nep, NEP_CONV_ABS, ierr))
  PetscCallA(NEPGetConvergenceTest(nep, conv, ierr))
  PetscCallA(NEPSetStoppingTest(nep, NEP_STOP_BASIC, ierr))
  PetscCallA(NEPGetStoppingTest(nep, stp, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,i2)') ' Convergence test =', conv, ', stopping test =', stp
  end if

  PetscCallA(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr))
  PetscCallA(NEPMonitorSet(nep, NEPMONITORFIRST, vf, PetscViewerAndFormatDestroy, ierr))
  PetscCallA(NEPMonitorConvergedCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, PETSC_NULL_VEC, vf, ierr))
  PetscCallA(NEPMonitorSet(nep, NEPMONITORCONVERGED, vf, PetscViewerAndFormatDestroy, ierr))
  PetscCallA(NEPMonitorCancel(nep, ierr))

  PetscCallA(NEPGetDS(nep, ds, ierr))
  PetscCallA(DSView(ds, PETSC_NULL_VIEWER, ierr))
  PetscCallA(NEPSetFromOptions(nep, ierr))

  PetscCallA(NEPGetRG(nep, rg, ierr))
  PetscCallA(RGView(rg, PETSC_NULL_VIEWER, ierr))

  PetscCallA(NEPSolve(nep, ierr))
  PetscCallA(NEPGetConvergedReason(nep, reason, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Finished - converged reason =', reason
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  PetscCallA(NEPErrorView(nep, NEP_ERROR_RELATIVE, PETSC_NULL_VIEWER, ierr))
  PetscCallA(NEPDestroy(nep, ierr))
  PetscCallA(MatDestroy(A(1), ierr))
  PetscCallA(MatDestroy(A(2), ierr))
  PetscCallA(MatDestroy(A(3), ierr))
  PetscCallA(FNDestroy(f(1), ierr))
  PetscCallA(FNDestroy(f(2), ierr))
  PetscCallA(FNDestroy(f(3), ierr))

  PetscCallA(SlepcFinalize(ierr))
end program test2f

!/*TEST
!
!   test:
!      suffix: 1
!      requires: !single
!
!TEST*/
