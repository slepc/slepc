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
      program main
#include <slepc/finclude/slepcsvd.h>
      use slepcsvd
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      Mat                A, B
      SVD                svd
      SVDConv            conv
      SVDStop            stp
      SVDWhich           which
      SVDConvergedReason reason
      PetscInt           m, n, i, Istart
      PetscInt           col(2), its, Iend
      PetscScalar        val(2)
      SVDProblemType     ptype
      PetscMPIInt        rank
      PetscErrorCode     ierr
      PetscBool          flg, tmode
      PetscViewerAndFormat vf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      m = 20
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-m',m,flg,ierr))
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      if (.not. flg) n = m+2

      if (rank .eq. 0) then
        write(*,100) m, n
      endif
 100  format (/'Bidiagonal matrix, m =',I3,', n=',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Build the Lauchli matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))

      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      val(1) = 1.0
      val(2) = 2.0
      do i=Istart,Iend-1
        col(1) = i
        col(2) = i+1
        if (i .le. n-1) then
          PetscCallA(MatSetValue(A,i,col(1),val(1),INSERT_VALUES,ierr))
        end if
        if (i .lt. n-1) then
          PetscCallA(MatSetValue(A,i,col(2),val(2),INSERT_VALUES,ierr))
        end if
      enddo

      PetscCallA(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute singular values
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SVDCreate(PETSC_COMM_WORLD,svd,ierr))
      PetscCallA(SVDSetOperators(svd,A,PETSC_NULL_MAT,ierr))

!     ** test some interface functions
      PetscCallA(SVDGetOperators(svd,B,PETSC_NULL_MAT,ierr))
      PetscCallA(MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr))
      PetscCallA(SVDSetConvergenceTest(svd,SVD_CONV_ABS,ierr))
      PetscCallA(SVDSetStoppingTest(svd,SVD_STOP_BASIC,ierr))

!     ** query properties and print them
      PetscCallA(SVDGetProblemType(svd,ptype,ierr))
      if (rank .eq. 0) then
        write(*,105) ptype
      endif
 105  format (/' Problem type = ',I2)
      PetscCallA(SVDIsGeneralized(svd,flg,ierr))
      if (flg .and. rank .eq. 0) then
        write(*,*) 'generalized'
      endif
      PetscCallA(SVDGetImplicitTranspose(svd,tmode,ierr))
      if (rank .eq. 0) then
        if (tmode) then
          write(*,110) 'implicit'
        else
          write(*,110) 'explicit'
        endif
      endif
 110  format (' Transpose mode is',A9)
      PetscCallA(SVDGetConvergenceTest(svd,conv,ierr))
      if (rank .eq. 0) then
        write(*,120) conv
      endif
 120  format (' Convergence test is',I2)
      PetscCallA(SVDGetStoppingTest(svd,stp,ierr))
      if (rank .eq. 0) then
        write(*,130) stp
      endif
 130  format (' Stopping test is',I2)
      PetscCallA(SVDGetWhichSingularTriplets(svd,which,ierr))
      if (rank .eq. 0) then
        if (which .eq. SVD_LARGEST) then
          write(*,140) 'largest'
        else
          write(*,140) 'smallest'
        endif
      endif
 140  format (' Which =',A9)

      PetscCallA(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr))
      PetscCallA(SVDMonitorSet(svd,SVDMONITORFIRST,vf,PetscViewerAndFormatDestroy,ierr))
      PetscCallA(SVDMonitorConvergedCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,PETSC_NULL_VEC,vf,ierr))
      PetscCallA(SVDMonitorSet(svd,SVDMONITORCONVERGED,vf,SVDMonitorConvergedDestroy,ierr))
      PetscCallA(SVDMonitorCancel(svd,ierr))

!     ** call the solver
      PetscCallA(SVDSetFromOptions(svd,ierr))
      PetscCallA(SVDSolve(svd,ierr))
      PetscCallA(SVDGetConvergedReason(svd,reason,ierr))
      if (rank .eq. 0) then
        write(*,150) reason
      endif
 150  format (' Converged reason:',I2)
      PetscCallA(SVDGetIterationNumber(svd,its,ierr))
!     if (rank .eq. 0) then
!       write(*,160) its
!     endif
!160  format (' Number of iterations of the method:',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SVDErrorView(svd,SVD_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      PetscCallA(SVDDestroy(svd,ierr))
      PetscCallA(MatDestroy(A,ierr))

      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -svd_type {{lanczos trlanczos cross cyclic randomized}}
!      filter: sed -e 's/2.99255/2.99254/'
!
!TEST*/
