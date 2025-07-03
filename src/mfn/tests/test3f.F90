!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test3f [-help] [-n <n>] [all SLEPc options]
!
!  Description: square root of the 2-D Laplacian.
!
!  The command line options are:
!    -n <n>, where <n> = matrix rows and columns
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcmfn.h>
      use slepcmfn
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      Mat                A, B
      MFN                mfn
      FN                 f
      MFNConvergedReason reason
      Vec                v, y
      PetscInt           Nt, n, i, j, II
      PetscInt           Istart, maxit, ncv
      PetscInt           col, its, Iend
      PetscScalar        val
      PetscReal          tol, norm
      PetscMPIInt        rank
      PetscErrorCode     ierr
      PetscBool          flg

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 4
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      Nt = n*n

      if (rank .eq. 0) then
        write(*,100) n
      endif
 100  format (/'nSquare root of Laplacian, n=',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the discrete 2-D Laplacian
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Nt,Nt,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))

      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      do II=Istart,Iend-1
        i = II/n
        j = II-i*n
        val = -1.0
        if (i .gt. 0) then
          col = II-n
          PetscCallA(MatSetValue(A,II,col,val,INSERT_VALUES,ierr))
        end if
        if (i .lt. n-1) then
          col = II+n
          PetscCallA(MatSetValue(A,II,col,val,INSERT_VALUES,ierr))
        end if
        if (j .gt. 0) then
          col = II-1
          PetscCallA(MatSetValue(A,II,col,val,INSERT_VALUES,ierr))
        end if
        if (j .lt. n-1) then
          col = II+1
          PetscCallA(MatSetValue(A,II,col,val,INSERT_VALUES,ierr))
        end if
        val = 4.0
        PetscCallA(MatSetValue(A,II,II,val,INSERT_VALUES,ierr))
      enddo

      PetscCallA(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))

      PetscCallA(MatCreateVecs(A,PETSC_NULL_VEC,v,ierr))
      i = 0
      val = 1.0
      PetscCallA(VecSetValue(v,i,val,INSERT_VALUES,ierr))
      PetscCallA(VecAssemblyBegin(v,ierr))
      PetscCallA(VecAssemblyEnd(v,ierr))
      PetscCallA(VecDuplicate(v,y,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute singular values
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MFNCreate(PETSC_COMM_WORLD,mfn,ierr))
      PetscCallA(MFNSetOperator(mfn,A,ierr))
      PetscCallA(MFNGetFN(mfn,f,ierr))
      PetscCallA(FNSetType(f,FNSQRT,ierr))

!     ** test some interface functions
      PetscCallA(MFNGetOperator(mfn,B,ierr))
      PetscCallA(MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr))
      PetscCallA(MFNSetOptionsPrefix(mfn,'myprefix_',ierr))
      tol = 1e-4
      maxit = 500
      PetscCallA(MFNSetTolerances(mfn,tol,maxit,ierr))
      ncv = 6
      PetscCallA(MFNSetDimensions(mfn,ncv,ierr))
      PetscCallA(MFNSetErrorIfNotConverged(mfn,PETSC_TRUE,ierr))
      PetscCallA(MFNSetFromOptions(mfn,ierr))

!     ** query properties and print them
      PetscCallA(MFNGetTolerances(mfn,tol,maxit,ierr))
      if (rank .eq. 0) then
        write(*,110) tol,maxit
      endif
 110  format (/' Tolerance: ',F7.4,', maxit: ',I4)
      PetscCallA(MFNGetDimensions(mfn,ncv,ierr))
      if (rank .eq. 0) then
        write(*,120) ncv
      endif
 120  format (' Subspace dimension: ',I3)
      PetscCallA(MFNGetErrorIfNotConverged(mfn,flg,ierr))
      if (rank .eq. 0 .and. flg) then
        write(*,*) 'Erroring out if convergence fails'
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Call the solver
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MFNSolve(mfn,v,y,ierr))
      PetscCallA(MFNGetConvergedReason(mfn,reason,ierr))
      if (rank .eq. 0) then
        write(*,130) reason
      endif
 130  format (' Converged reason:',I2)
      PetscCallA(MFNGetIterationNumber(mfn,its,ierr))
!     if (rank .eq. 0) then
!       write(*,140) its
!     endif
!140  format (' Number of iterations of the method:',I4)

      PetscCallA(VecNorm(y,NORM_2,norm,ierr))
      if (rank .eq. 0) then
        write(*,150) norm
      endif
 150  format (' sqrt(A)*v has norm ',F7.4)

      PetscCallA(MFNDestroy(mfn,ierr))
      PetscCallA(MatDestroy(A,ierr))
      PetscCallA(VecDestroy(v,ierr))
      PetscCallA(VecDestroy(y,ierr))

      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -log_exclude mfn
!
!TEST*/
