!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test7f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Simple example that solves an eigensystem with the EPS object.
!  Same problem as ex1f but with simplified output.
!
!  The command line options are:
!    -n <n>, where <n> = number of grid points = matrix size
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepceps.h>
      use slepceps
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     A     operator matrix
!     eps   eigenproblem solver context

      Mat            A
      EPS            eps
      EPSType        tname
      PetscInt       n, i, Istart, Iend
      PetscInt       nev, nini
      PetscInt       col(3)
      PetscInt       i0, i1, i2, i3
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg
      PetscScalar    value(3), one
      Vec            v(1)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 30
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))

      if (rank .eq. 0) then
        write(*,100) n
      endif
 100  format (/'1-D Laplacian Eigenproblem, n =',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the operator matrix that defines the eigensystem, Ax=kx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))

      i0 = 0
      i1 = 1
      i2 = 2
      i3 = 3
      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      if (Istart .eq. 0) then
        i = 0
        col(1) = 0
        col(2) = 1
        value(1) =  2.0
        value(2) = -1.0
        PetscCallA(MatSetValues(A,i1,[i],i2,col,value,INSERT_VALUES,ierr))
        Istart = Istart+1
      endif
      if (Iend .eq. n) then
        i = n-1
        col(1) = n-2
        col(2) = n-1
        value(1) = -1.0
        value(2) =  2.0
        PetscCallA(MatSetValues(A,i1,[i],i2,col,value,INSERT_VALUES,ierr))
        Iend = Iend-1
      endif
      value(1) = -1.0
      value(2) =  2.0
      value(3) = -1.0
      do i=Istart,Iend-1
        col(1) = i-1
        col(2) = i
        col(3) = i+1
        PetscCallA(MatSetValues(A,i1,[i],i3,col,value,INSERT_VALUES,ierr))
      enddo

      PetscCallA(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))

      PetscCallA(MatCreateVecs(A,v(1),PETSC_NULL_VEC,ierr))
      one = 1.0
      if (Istart .eq. 0) then
        PetscCallA(VecSetValue(v(1),i0,one,INSERT_VALUES,ierr))
      endif
      PetscCallA(VecAssemblyBegin(v(1),ierr))
      PetscCallA(VecAssemblyEnd(v(1),ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      PetscCallA(EPSCreate(PETSC_COMM_WORLD,eps,ierr))

!     ** Set operators. In this case, it is a standard eigenvalue problem
      PetscCallA(EPSSetOperators(eps,A,PETSC_NULL_MAT,ierr))
      PetscCallA(EPSSetProblemType(eps,EPS_HEP,ierr))

!     ** Set solver parameters at runtime
      PetscCallA(EPSSetFromOptions(eps,ierr))

!     ** Set initial vectors
      nini = 1
      PetscCallA(EPSSetInitialSpace(eps,nini,v,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSSolve(eps,ierr))

!     ** Optional: Get some information from the solver and display it
      PetscCallA(EPSGetType(eps,tname,ierr))
      if (rank .eq. 0) then
        write(*,120) tname
      endif
 120  format (' Solution method: ',A)
      PetscCallA(EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
      if (rank .eq. 0) then
        write(*,130) nev
      endif
 130  format (' Number of requested eigenvalues:',I2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      PetscCallA(EPSDestroy(eps,ierr))
      PetscCallA(MatDestroy(A,ierr))
      PetscCallA(VecDestroy(v(1),ierr))

      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -eps_nev 4 -eps_ncv 19
!      filter: sed -e "s/83791/83792/"
!
!TEST*/
