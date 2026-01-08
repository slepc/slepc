!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Simple example to test the EPS Fortran interface.
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepceps.h>
program test14f
  use slepceps
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Mat                  :: A, B
  EPS                  :: eps
  ST                   :: st
  KSP                  :: ksp
  DS                   :: ds
  PetscReal            :: cut, tol, tolabs
  PetscScalar          :: tget, val
  PetscInt             :: n, i, its, Istart, Iend
  PetscInt             :: nev, ncv, mpd
  PetscBool            :: flg
  EPSConvergedReason   :: reason
  EPSType              :: tname
  EPSExtraction        :: extr
  EPSBalance           :: bal
  EPSWhich             :: which
  EPSConv              :: conv
  EPSStop              :: stp
  EPSProblemType       :: ptype
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
    write (*, '(/a,i3,a)') 'Diagonal Eigenproblem, n =', n, ' (Fortran)'
  end if

  PetscCallA(MatCreate(PETSC_COMM_WORLD, A, ierr))
  PetscCallA(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr))
  PetscCallA(MatSetFromOptions(A, ierr))
  PetscCallA(MatGetOwnershipRange(A, Istart, Iend, ierr))
  do i = Istart, Iend - 1
    val = i + 1
    PetscCallA(MatSetValue(A, i, i, val, INSERT_VALUES, ierr))
  end do
  PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create eigensolver and test interface functions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(EPSCreate(PETSC_COMM_WORLD, eps, ierr))
  PetscCallA(EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr))
  PetscCallA(EPSGetOperators(eps, B, PETSC_NULL_MAT, ierr))
  PetscCallA(MatView(B, PETSC_NULL_VIEWER, ierr))

  PetscCallA(EPSSetType(eps, EPSKRYLOVSCHUR, ierr))
  PetscCallA(EPSGetType(eps, tname, ierr))
  if (rank == 0) then
    write (*, '(a,a)') ' Type set to ', tname
  end if

  PetscCallA(EPSGetProblemType(eps, ptype, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Problem type before changing = ', ptype
  end if
  PetscCallA(EPSSetProblemType(eps, EPS_HEP, ierr))
  PetscCallA(EPSGetProblemType(eps, ptype, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' ... changed to ', ptype
  end if
  PetscCallA(EPSIsGeneralized(eps, flg, ierr))
  if (flg .and. rank == 0) then
    write (*, *) 'generalized'
  end if
  PetscCallA(EPSIsHermitian(eps, flg, ierr))
  if (flg .and. rank == 0) then
    write (*, *) 'hermitian'
  end if
  PetscCallA(EPSIsPositive(eps, flg, ierr))
  if (flg .and. rank == 0) then
    write (*, *) 'positive'
  end if

  PetscCallA(EPSGetExtraction(eps, extr, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Extraction before changing = ', extr
  end if
  PetscCallA(EPSSetExtraction(eps, EPS_HARMONIC, ierr))
  PetscCallA(EPSGetExtraction(eps, extr, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' ... changed to ', extr
  end if

  its = 8
  cut = 2.0e-6
  bal = EPS_BALANCE_ONESIDE
  PetscCallA(EPSSetBalance(eps, bal, its, cut, ierr))
  PetscCallA(EPSGetBalance(eps, bal, its, cut, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,i2,a,f9.6)') ' Balance: ', bal, ', its=', its, ', cutoff=', cut
  end if

  tget = 4.8
  PetscCallA(EPSSetTarget(eps, tget, ierr))
  PetscCallA(EPSGetTarget(eps, tget, ierr))
  PetscCallA(EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE, ierr))
  PetscCallA(EPSGetWhichEigenpairs(eps, which, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,f4.1)') ' Which = ', which, ', target = ', PetscRealPart(tget)
  end if

  nev = 4
  PetscCallA(EPSSetDimensions(eps, nev, PETSC_DETERMINE_INTEGER, PETSC_DETERMINE_INTEGER, ierr))
  PetscCallA(EPSGetDimensions(eps, nev, ncv, mpd, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,i2,a,i2)') ' Dimensions: nev=', nev, ', ncv=', ncv, ', mpd=', mpd
  end if

  tol = 2.2e-4
  its = 200
  PetscCallA(EPSSetTolerances(eps, tol, its, ierr))
  PetscCallA(EPSGetTolerances(eps, tol, its, ierr))
  if (rank == 0) then
    write (*, '(a,f8.5,a,i4)') ' Tolerance =', tol, ', max_its =', its
  end if

  PetscCallA(EPSSetConvergenceTest(eps, EPS_CONV_ABS, ierr))
  PetscCallA(EPSGetConvergenceTest(eps, conv, ierr))
  PetscCallA(EPSSetStoppingTest(eps, EPS_STOP_BASIC, ierr))
  PetscCallA(EPSGetStoppingTest(eps, stp, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,i2)') ' Convergence test =', conv, ', stopping test =', stp
  end if

  PetscCallA(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr))
  PetscCallA(EPSMonitorSet(eps, EPSMONITORFIRST, vf, PetscViewerAndFormatDestroy, ierr))
  PetscCallA(EPSMonitorConvergedCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, PETSC_NULL, vf, ierr))
  PetscCallA(EPSMonitorSet(eps, EPSMONITORCONVERGED, vf, PetscViewerAndFormatDestroy, ierr))
  PetscCallA(EPSMonitorCancel(eps, ierr))

  PetscCallA(EPSGetST(eps, st, ierr))
  PetscCallA(STGetKSP(st, ksp, ierr))
  tol = 1.e-8
  tolabs = 1.e-35
  PetscCallA(KSPSetTolerances(ksp, tol, tolabs, PETSC_CURRENT_REAL, PETSC_CURRENT_INTEGER, ierr))
  PetscCallA(STView(st, PETSC_NULL_VIEWER, ierr))
  PetscCallA(EPSGetDS(eps, ds, ierr))
  PetscCallA(DSView(ds, PETSC_NULL_VIEWER, ierr))

  PetscCallA(EPSSetFromOptions(eps, ierr))
  PetscCallA(EPSSolve(eps, ierr))
  PetscCallA(EPSGetConvergedReason(eps, reason, ierr))
  PetscCallA(EPSGetIterationNumber(eps, its, ierr))
  if (rank == 0) then
    write (*, '(a,i2,a,i4)') ' Finished - converged reason =', reason, ', its = ', its
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  PetscCallA(EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_NULL_VIEWER, ierr))
  PetscCallA(EPSDestroy(eps, ierr))
  PetscCallA(MatDestroy(A, ierr))

  PetscCallA(SlepcFinalize(ierr))
end program test14f

!/*TEST
!
!   test:
!      suffix: 1
!      args: -eps_ncv 14
!      filter: sed -e "s/00001/00000/" | sed -e "s/4.99999/5.00000/" | sed -e "s/5.99999/6.00000/"
!
!TEST*/
