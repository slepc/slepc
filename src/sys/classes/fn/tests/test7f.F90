!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test7f [-help] [-n <n>] [-verbose] [-inplace]
!
!  Description: Simple example that tests the matrix square root.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcfn.h>
      use slepcfn
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Mat            A,S,R
      FN             fn
      PetscInt       n
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg,verbose,inplace
      PetscReal      re,im,nrm
      PetscScalar    tau,eta,alpha,x,y,yp
      PetscScalar, pointer :: aa(:,:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 10
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-verbose',verbose,ierr))
      PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-inplace',inplace,ierr))

      if (rank .eq. 0) then
        write(*,100) n
      endif
 100  format (/'Matrix square root, n =',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create FN object and matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(FNCreate(PETSC_COMM_WORLD,fn,ierr))
      PetscCallA(FNSetType(fn,FNSQRT,ierr))
      tau = 0.15
      eta = 1.0
      PetscCallA(FNSetScale(fn,tau,eta,ierr))
      PetscCallA(FNSetFromOptions(fn,ierr))
      PetscCallA(FNGetScale(fn,tau,eta,ierr))
      PetscCallA(FNView(fn,PETSC_NULL_VIEWER,ierr))

      PetscCallA(MatCreateSeqDense(PETSC_COMM_SELF,n,n,PETSC_NULL_SCALAR_ARRAY,A,ierr))
      PetscCallA(PetscObjectSetName(A,'A',ierr))
      PetscCallA(MatDenseGetArray(A,aa,ierr))
      call FillUpMatrix(n,aa)
      PetscCallA(MatDenseRestoreArray(A,aa,ierr))
      PetscCallA(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Scalar evaluation
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      x = 2.2
      PetscCallA(FNEvaluateFunction(fn,x,y,ierr))
      PetscCallA(FNEvaluateDerivative(fn,x,yp,ierr))

      if (rank .eq. 0) then
        re = PetscRealPart(y)
        im = PetscImaginaryPart(y)
        if (abs(im).lt.1.d-10) then
          write(*,110) 'f', PetscRealPart(x), re
        else
          write(*,120) 'f', PetscRealPart(x), re, im
        endif
        re = PetscRealPart(yp)
        im = PetscImaginaryPart(yp)
        if (abs(im).lt.1.d-10) then
          write(*,110) 'f''', PetscRealPart(x), re
        else
          write(*,120) 'f''', PetscRealPart(x), re, im
        endif
      endif
 110  format (A2,'(',F4.1,') = ',F8.5)
 120  format (A2,'(',F4.1,') = ',F8.5,SP,F8.5,'i')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute matrix square root
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreateSeqDense(PETSC_COMM_SELF,n,n,PETSC_NULL_SCALAR_ARRAY,S,ierr))
      PetscCallA(PetscObjectSetName(S,'S',ierr))
      if (inplace) then
        PetscCallA(MatCopy(A,S,SAME_NONZERO_PATTERN,ierr))
        PetscCallA(MatSetOption(S,MAT_HERMITIAN,PETSC_TRUE,ierr))
        PetscCallA(FNEvaluateFunctionMat(fn,S,PETSC_NULL_MAT,ierr))
      else
        PetscCallA(FNEvaluateFunctionMat(fn,A,S,ierr))
      endif
      if (verbose) then
        if (rank .eq. 0) write (*,*) 'Matrix A - - - - - - - -'
        PetscCallA(MatView(A,PETSC_NULL_VIEWER,ierr))
        if (rank .eq. 0) write (*,*) 'Computed sqrtm(A) - - - - - - - -'
        PetscCallA(MatView(S,PETSC_NULL_VIEWER,ierr))
      endif

!     *** check error ||S*S-A||_F
      PetscCallA(MatMatMult(S,S,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,R,ierr))
      if (eta .ne. 1.0) then
        alpha = 1.0/(eta*eta)
        PetscCallA(MatScale(R,alpha,ierr))
      endif
      alpha = -tau
      PetscCallA(MatAXPY(R,alpha,A,SAME_NONZERO_PATTERN,ierr))
      PetscCallA(MatNorm(R,NORM_FROBENIUS,nrm,ierr))
      if (nrm<100*PETSC_MACHINE_EPSILON) then
        write (*,*) '||S*S-A||_F < 100*eps'
      else
        write (*,130) nrm
      endif
 130  format ('||S*S-A||_F = ',F8.5)

!     *** Clean up
      PetscCallA(MatDestroy(S,ierr))
      PetscCallA(MatDestroy(R,ierr))
      PetscCallA(MatDestroy(A,ierr))
      PetscCallA(FNDestroy(fn,ierr))
      PetscCallA(SlepcFinalize(ierr))
      end

! -----------------------------------------------------------------

      subroutine FillUpMatrix(n,X)
      PetscInt    n,i,j
      PetscScalar X(n,n)

      do i=1,n
        X(i,i) = 2.5
      end do
      do j=1,2
        do i=1,n-j
          X(i,i+j) = 1.d0
          X(i+j,i) = 1.d0
        end do
      end do

      end

!/*TEST
!
!   test:
!      suffix: 1
!      nsize: 1
!      args: -fn_scale .13,2 -n 19 -fn_method {{0 1 2 3}shared output}
!      filter: grep -v "computing matrix functions"
!      output_file: output/test7f_1.out
!
!TEST*/
