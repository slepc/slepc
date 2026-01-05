!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test4f [-help] [-n <n>] [-m <m>] [all SLEPc options]
!
!  Description: Singular value decomposition of a bidiagonal matrix.
!
!               |  1  2                     |
!               |     1  2                  |
!               |        1  2               |
!           A = |          .  .             |
!               |             .  .          |
!               |                1  2       |
!               |                   1  2    |
!
!  The command line options are:
!    -m <m>, where <m> = matrix rows.
!    -n <n>, where <n> = matrix columns (defaults to m+2).
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepcsvd.h>
program test4f
  use slepcsvd
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  Mat                  :: A, B
  SVD                  :: svd
  SVDConv              :: conv
  SVDStop              :: stp
  SVDWhich             :: which
  SVDConvergedReason   :: reason
  PetscInt             :: m, n, i, Istart
  PetscInt             :: col(2), its, Iend
  PetscScalar          :: val(2)
  SVDProblemType       :: ptype
  PetscMPIInt          :: rank
  PetscErrorCode       :: ierr
  PetscBool            :: flg, tmode
  PetscViewerAndFormat :: vf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
  m = 20
  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-m', m, flg, ierr))
  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-n', n, flg, ierr))
  if (.not. flg) n = m + 2

  if (rank == 0) then
    write (*, '(/a,i3,a,i3,a)') 'Bidiagonal matrix, m =', m, ', n=', n, ' (Fortran)'
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Build the Lauchli matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(MatCreate(PETSC_COMM_WORLD, A, ierr))
  PetscCallA(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m, n, ierr))
  PetscCallA(MatSetFromOptions(A, ierr))

  PetscCallA(MatGetOwnershipRange(A, Istart, Iend, ierr))
  val(1) = 1.0
  val(2) = 2.0
  do i = Istart, Iend - 1
    col(1) = i
    col(2) = i + 1
    if (i < n) then
      PetscCallA(MatSetValue(A, i, col(1), val(1), INSERT_VALUES, ierr))
    end if
    if (i < n - 1) then
      PetscCallA(MatSetValue(A, i, col(2), val(2), INSERT_VALUES, ierr))
    end if
  end do

  PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute singular values
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SVDCreate(PETSC_COMM_WORLD, svd, ierr))
  PetscCallA(SVDSetOperators(svd, A, PETSC_NULL_MAT, ierr))

! ** test some interface functions
  PetscCallA(SVDGetOperators(svd, B, PETSC_NULL_MAT, ierr))
  PetscCallA(MatView(B, PETSC_VIEWER_STDOUT_WORLD, ierr))
  PetscCallA(SVDSetConvergenceTest(svd, SVD_CONV_ABS, ierr))
  PetscCallA(SVDSetStoppingTest(svd, SVD_STOP_BASIC, ierr))

! ** query properties and print them
  PetscCallA(SVDGetProblemType(svd, ptype, ierr))
  if (rank == 0) then
    write (*, '(/a,i2)') ' Problem type = ', ptype
  end if
  PetscCallA(SVDIsGeneralized(svd, flg, ierr))
  if (flg .and. rank == 0) then
    write (*, *) 'generalized'
  end if
  PetscCallA(SVDGetImplicitTranspose(svd, tmode, ierr))
  if (rank == 0) then
    if (tmode) then
      write (*, *) ' Transpose mode is implicit'
    else
      write (*, *) ' Transpose mode is explicit'
    end if
  end if
  PetscCallA(SVDGetConvergenceTest(svd, conv, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Convergence test is', conv
  end if
  PetscCallA(SVDGetStoppingTest(svd, stp, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Stopping test is', stp
  end if
  PetscCallA(SVDGetWhichSingularTriplets(svd, which, ierr))
  if (rank == 0) then
    if (which == SVD_LARGEST) then
      write (*, *) ' Which = largest'
    else
      write (*, *) ' Which = smallest'
    end if
  end if

  PetscCallA(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr))
  PetscCallA(SVDMonitorSet(svd, SVDMONITORFIRST, vf, PetscViewerAndFormatDestroy, ierr))
  PetscCallA(SVDMonitorConvergedCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, PETSC_NULL, vf, ierr))
  PetscCallA(SVDMonitorSet(svd, SVDMONITORCONVERGED, vf, PetscViewerAndFormatDestroy, ierr))
  PetscCallA(SVDMonitorCancel(svd, ierr))

! ** call the solver
  PetscCallA(SVDSetFromOptions(svd, ierr))
  PetscCallA(SVDSolve(svd, ierr))
  PetscCallA(SVDGetConvergedReason(svd, reason, ierr))
  if (rank == 0) then
    write (*, '(a,i2)') ' Converged reason:', reason
  end if
  PetscCallA(SVDGetIterationNumber(svd, its, ierr))
! if (rank==0) then
!   write(*,'(a,i4)') ' Number of iterations of the method:', its
! end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SVDErrorView(svd, SVD_ERROR_RELATIVE, PETSC_NULL_VIEWER, ierr))
  PetscCallA(SVDDestroy(svd, ierr))
  PetscCallA(MatDestroy(A, ierr))

  PetscCallA(SlepcFinalize(ierr))
end program test4f

!/*TEST
!
!   test:
!      suffix: 1
!      args: -svd_type {{lanczos trlanczos cross cyclic randomized}}
!      filter: sed -e 's/2.99255/2.99254/'
!
!TEST*/
