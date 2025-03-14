!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex22f [-n <n>] [-tau <tau>] [SLEPc opts]
!
!  Description: Delay differential equation. Fortran90 equivalent of ex22.c
!
!  The command line options are:
!    -n <n>, where <n> = number of grid subdivisions
!    -tau <tau>, where <tau> = delay parameter
!
! ----------------------------------------------------------------------
!  Solve parabolic partial differential equation with time delay tau
!
!           u_t = u_xx + aa*u(t) + bb*u(t-tau)
!           u(0,t) = u(pi,t) = 0
!
!  with aa = 20 and bb(x) = -4.1+x*(1-exp(x-pi)).
!
!  Discretization leads to a DDE of dimension n
!
!           -u' = A*u(t) + B*u(t-tau)
!
!  which results in the nonlinear eigenproblem
!
!           (-lambda*I + A + exp(-tau*lambda)*B)*u = 0
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcnep.h>
      use slepcnep
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     nep       nonlinear eigensolver context
!     Id,A,B    problem matrices
!     f1,f2,f3  functions to define the nonlinear operator

      Mat            Id, A, B, mats(3)
      FN             f1, f2, f3, funs(3)
      NEP            nep
      NEPType        tname
      PetscScalar    one, bb, coeffs(2), scal
      PetscReal      tau, h, aa, xi, tol
      PetscInt       n, i, k, nev, Istart, Iend
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg, terse

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 128
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      tau = 0.001
      PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-tau',tau,flg,ierr))
      if (rank .eq. 0) then
        write(*,100) n, tau
      endif
 100  format (/'Delay Eigenproblem, n =',I4,', tau =',F6.3)

      one = 1.0
      aa = 20.0
      h = PETSC_PI/real(n+1)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create problem matrices
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Id is the identity matrix
      PetscCallA(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,one,Id,ierr))
      PetscCallA(MatSetOption(Id,MAT_HERMITIAN,PETSC_TRUE,ierr))

!     ** A = 1/h^2*tridiag(1,-2,1) + aa*I
      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))
      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      coeffs(1) = 1.0/(h*h)
      coeffs(2) = -2.0/(h*h)+aa
      do i=Istart,Iend-1
        if (i .gt. 0) then
          PetscCallA(MatSetValue(A,i,i-1,coeffs(1),INSERT_VALUES,ierr))
        endif
        if (i .lt. n-1) then
          PetscCallA(MatSetValue(A,i,i+1,coeffs(1),INSERT_VALUES,ierr))
        endif
        PetscCallA(MatSetValue(A,i,i,coeffs(2),INSERT_VALUES,ierr))
      end do
      PetscCallA(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE,ierr))

!     ** B = diag(bb(xi))
      PetscCallA(MatCreate(PETSC_COMM_WORLD,B,ierr))
      PetscCallA(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(B,ierr))
      PetscCallA(MatGetOwnershipRange(B,Istart,Iend,ierr))
      do i=Istart,Iend-1
        xi = (i+1)*h
        bb  = -4.1+xi*(1.0-exp(xi-PETSC_PI))
        PetscCallA(MatSetValue(B,i,i,bb,INSERT_VALUES,ierr))
      end do
      PetscCallA(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatSetOption(B,MAT_HERMITIAN,PETSC_TRUE,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create problem functions, f1=-lambda, f2=1.0, f3=exp(-tau*lambda)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(FNCreate(PETSC_COMM_WORLD,f1,ierr))
      PetscCallA(FNSetType(f1,FNRATIONAL,ierr))
      k = 2
      coeffs(1) = -1.0
      coeffs(2) = 0.0
      PetscCallA(FNRationalSetNumerator(f1,k,coeffs,ierr))

      PetscCallA(FNCreate(PETSC_COMM_WORLD,f2,ierr))
      PetscCallA(FNSetType(f2,FNRATIONAL,ierr))
      k = 1
      coeffs(1) = 1.0
      PetscCallA(FNRationalSetNumerator(f2,k,coeffs,ierr))

      PetscCallA(FNCreate(PETSC_COMM_WORLD,f3,ierr))
      PetscCallA(FNSetType(f3,FNEXP,ierr))
      scal = -tau
      PetscCallA(FNSetScale(f3,scal,one,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      PetscCallA(NEPCreate(PETSC_COMM_WORLD,nep,ierr))

!     ** Set the split operator. Note that A is passed first so that
!     ** SUBSET_NONZERO_PATTERN can be used
      k = 3
      mats(1) = A
      mats(2) = Id
      mats(3) = B
      funs(1) = f2
      funs(2) = f1
      funs(3) = f3
      PetscCallA(NEPSetSplitOperator(nep,k,mats,funs,SUBSET_NONZERO_PATTERN,ierr))
      PetscCallA(NEPSetProblemType(nep,NEP_GENERAL,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Customize nonlinear solver; set runtime options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      tol = 1e-9
      PetscCallA(NEPSetTolerances(nep,tol,PETSC_CURRENT_INTEGER,ierr))
      k = 1
      PetscCallA(NEPSetDimensions(nep,k,PETSC_DETERMINE_INTEGER,PETSC_DETERMINE_INTEGER,ierr))
      k = 0
      PetscCallA(NEPRIISetLagPreconditioner(nep,k,ierr))

!     ** Set solver parameters at runtime
      PetscCallA(NEPSetFromOptions(nep,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(NEPSolve(nep,ierr))

!     ** Optional: Get some information from the solver and display it
      PetscCallA(NEPGetType(nep,tname,ierr))
      if (rank .eq. 0) then
        write(*,120) tname
      endif
 120  format (' Solution method: ',A)
      PetscCallA(NEPGetDimensions(nep,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
      if (rank .eq. 0) then
        write(*,130) nev
      endif
 130  format (' Number of requested eigenvalues:',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** show detailed info unless -terse option is given by user
      PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-terse',terse,ierr))
      if (terse) then
        PetscCallA(NEPErrorView(nep,NEP_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      else
        PetscCallA(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr))
        PetscCallA(NEPConvergedReasonView(nep,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(NEPErrorView(nep,NEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr))
      endif
      PetscCallA(NEPDestroy(nep,ierr))
      PetscCallA(MatDestroy(Id,ierr))
      PetscCallA(MatDestroy(A,ierr))
      PetscCallA(MatDestroy(B,ierr))
      PetscCallA(FNDestroy(f1,ierr))
      PetscCallA(FNDestroy(f2,ierr))
      PetscCallA(FNDestroy(f3,ierr))
      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -terse
!      requires: !single
!      filter: sed -e "s/[+-]0\.0*i//g"
!
!TEST*/
