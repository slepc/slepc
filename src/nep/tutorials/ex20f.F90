!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex20f [-n <n>] [SLEPc opts]
!
!  Description: Simple 1-D nonlinear eigenproblem. Fortran90 equivalent of ex20.c
!
!  The command line options are:
!    -n <n>, where <n> = number of grid subdivisions
!
! ----------------------------------------------------------------------
!  Solve 1-D PDE
!           -u'' = lambda*u
!  on [0,1] subject to
!           u(0)=0, u'(1)=u(1)*lambda*kappa/(kappa-lambda)
! ----------------------------------------------------------------------
!

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     User-defined application context
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      module UserModule
#include <slepc/finclude/slepcnep.h>
      use slepcnep
      type User
        PetscScalar kappa
        PetscReal   h
      end type User
      end module UserModule

      program main
#include <slepc/finclude/slepcnep.h>
      use UserModule
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     nep       nonlinear eigensolver context
!     x         eigenvector
!     lambda    eigenvalue
!     F,J       Function and Jacobian matrices
!     ctx       user-defined context

      NEP            nep
      Vec            x, v(1)
      PetscScalar    lambda
      Mat            F, J
      type(User)     ctx
      NEPType        tname
      PetscInt       n, i, k, nev, its, maxit, nconv, three, one
      PetscReal      tol, norm
      PetscScalar    alpha
      PetscMPIInt    rank
      PetscBool      flg
      PetscErrorCode ierr
!  Note: Any user-defined Fortran routines (such as FormJacobian)
!  MUST be declared as external.
      external       FormFunction, FormJacobian

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 128
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      if (rank .eq. 0) then
        write(*,'(/A,I4)') 'Nonlinear Eigenproblem, n =',n
      endif

      ctx%h = 1.0/real(n)
      ctx%kappa = 1.0

      three = 3
      one = 1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create matrix data structure to hold the Function and the Jacobian
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreate(PETSC_COMM_WORLD,F,ierr))
      PetscCallA(MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(F,ierr))
      PetscCallA(MatSeqAIJSetPreallocation(F,three,PETSC_NULL_INTEGER_ARRAY,ierr))
      PetscCallA(MatMPIAIJSetPreallocation(F,three,PETSC_NULL_INTEGER_ARRAY,one,PETSC_NULL_INTEGER_ARRAY,ierr))

      PetscCallA(MatCreate(PETSC_COMM_WORLD,J,ierr))
      PetscCallA(MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(J,ierr))
      PetscCallA(MatSeqAIJSetPreallocation(J,three,PETSC_NULL_INTEGER_ARRAY,ierr))
      PetscCallA(MatMPIAIJSetPreallocation(J,three,PETSC_NULL_INTEGER_ARRAY,one,PETSC_NULL_INTEGER_ARRAY,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      PetscCallA(NEPCreate(PETSC_COMM_WORLD,nep,ierr))

!     ** Set routines for evaluation of Function and Jacobian
      PetscCallA(NEPSetFunction(nep,F,F,FormFunction,ctx,ierr))
      PetscCallA(NEPSetJacobian(nep,J,FormJacobian,ctx,ierr))

!     ** Customize nonlinear solver
      tol = 1e-9
      PetscCallA(NEPSetTolerances(nep,tol,PETSC_CURRENT_INTEGER,ierr))
      k = 1
      PetscCallA(NEPSetDimensions(nep,k,PETSC_DETERMINE_INTEGER,PETSC_DETERMINE_INTEGER,ierr))

!     ** Set solver parameters at runtime
      PetscCallA(NEPSetFromOptions(nep,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Evaluate initial guess
      PetscCallA(MatCreateVecs(F,x,PETSC_NULL_VEC,ierr))
      PetscCallA(VecDuplicate(x,v(1),ierr))
      alpha = 1.0
      PetscCallA(VecSet(v(1),alpha,ierr))
      k = 1
      PetscCallA(NEPSetInitialSpace(nep,k,v,ierr))
      PetscCallA(VecDestroy(v(1),ierr))

!     ** Call the solver
      PetscCallA(NEPSolve(nep,ierr))
      PetscCallA(NEPGetIterationNumber(nep,its,ierr))
      if (rank .eq. 0) then
        write(*,'(A,I3)') ' Number of NEP iterations =',its
      endif

!     ** Optional: Get some information from the solver and display it
      PetscCallA(NEPGetType(nep,tname,ierr))
      if (rank .eq. 0) then
        write(*,'(A,A10)') ' Solution method: ',tname
      endif
      PetscCallA(NEPGetDimensions(nep,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
      if (rank .eq. 0) then
        write(*,'(A,I4)') ' Number of requested eigenvalues:',nev
      endif
      PetscCallA(NEPGetTolerances(nep,tol,maxit,ierr))
      if (rank .eq. 0) then
        write(*,'(A,F12.9,A,I5)') ' Stopping condition: tol=',tol,', maxit=',maxit
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(NEPGetConverged(nep,nconv,ierr))
      if (rank .eq. 0) then
        write(*,'(A,I2/)') ' Number of converged approximate eigenpairs:',nconv
      endif

!     ** Display eigenvalues and relative errors
      if (nconv .gt. 0) then
        if (rank .eq. 0) then
          write(*,*) '        k              ||T(k)x||'
          write(*,*) '----------------- ------------------'
        endif
        do i=0,nconv-1
!         ** Get converged eigenpairs: (in this example they are always real)
          PetscCallA(NEPGetEigenpair(nep,i,lambda,PETSC_NULL_SCALAR,x,PETSC_NULL_VEC,ierr))

!         ** Compute residual norm and error
          PetscCallA(NEPComputeError(nep,i,NEP_ERROR_RELATIVE,norm,ierr))
          if (rank .eq. 0) then
            write(*,'(1P,E15.4,E18.4)') PetscRealPart(lambda), norm
          endif
        enddo
        if (rank .eq. 0) then
          write(*,*)
        endif
      endif

      PetscCallA(NEPDestroy(nep,ierr))
      PetscCallA(MatDestroy(F,ierr))
      PetscCallA(MatDestroy(J,ierr))
      PetscCallA(VecDestroy(x,ierr))
      PetscCallA(SlepcFinalize(ierr))
      end

! ---------------  Evaluate Function matrix  T(lambda)  ----------------

      subroutine FormFunction(nep,lambda,fun,B,ctx,ierr)
      use UserModule
      implicit none
      NEP            nep
      PetscScalar    lambda, A(3), c, d
      Mat            fun,B
      type(User)     ctx
      PetscReal      h
      PetscInt       i, n, j(3), Istart, Iend, one, two, three
      PetscErrorCode ierr

!     ** Compute Function entries and insert into matrix
      PetscCall(MatGetSize(fun,n,PETSC_NULL_INTEGER,ierr))
      PetscCall(MatGetOwnershipRange(fun,Istart,Iend,ierr))
      h = ctx%h
      c = ctx%kappa/(lambda-ctx%kappa)
      d = n
      one = 1
      two = 2
      three = 3

!     ** Boundary points
      if (Istart .eq. 0) then
        i = 0
        j(1) = 0
        j(2) = 1
        A(1) = 2.0*(d-lambda*h/3.0)
        A(2) = -d-lambda*h/6.0
        PetscCall(MatSetValues(fun,one,[i],two,j,A,INSERT_VALUES,ierr))
        Istart = Istart + 1
      endif

      if (Iend .eq. n) then
        i = n-1
        j(1) = n-2
        j(2) = n-1
        A(1) = -d-lambda*h/6.0
        A(2) = d-lambda*h/3.0+lambda*c
        PetscCall(MatSetValues(fun,one,[i],two,j,A,INSERT_VALUES,ierr))
        Iend = Iend - 1
      endif

!     ** Interior grid points
      do i=Istart,Iend-1
        j(1) = i-1
        j(2) = i
        j(3) = i+1
        A(1) = -d-lambda*h/6.0
        A(2) = 2.0*(d-lambda*h/3.0)
        A(3) = -d-lambda*h/6.0
        PetscCall(MatSetValues(fun,one,[i],three,j,A,INSERT_VALUES,ierr))
      enddo

!     ** Assemble matrix
      PetscCall(MatAssemblyBegin(fun,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatAssemblyEnd(fun,MAT_FINAL_ASSEMBLY,ierr))

      end

! ---------------  Evaluate Jacobian matrix  T'(lambda)  ---------------

      subroutine FormJacobian(nep,lambda,jac,ctx,ierr)
      use UserModule
      implicit none
      NEP            nep
      PetscScalar    lambda, A(3), c
      Mat            jac
      type(User)     ctx
      PetscReal      h
      PetscInt       i, n, j(3), Istart, Iend, one, two, three
      PetscErrorCode ierr

!     ** Compute Jacobian entries and insert into matrix
      PetscCall(MatGetSize(jac,n,PETSC_NULL_INTEGER,ierr))
      PetscCall(MatGetOwnershipRange(jac,Istart,Iend,ierr))
      h = ctx%h
      c = ctx%kappa/(lambda-ctx%kappa)
      one = 1
      two = 2
      three = 3

!     ** Boundary points
      if (Istart .eq. 0) then
        i = 0
        j(1) = 0
        j(2) = 1
        A(1) = -2.0*h/3.0
        A(2) = -h/6.0
        PetscCall(MatSetValues(jac,one,[i],two,j,A,INSERT_VALUES,ierr))
        Istart = Istart + 1
      endif

      if (Iend .eq. n) then
        i = n-1
        j(1) = n-2
        j(2) = n-1
        A(1) = -h/6.0
        A(2) = -h/3.0-c*c
        PetscCall(MatSetValues(jac,one,[i],two,j,A,INSERT_VALUES,ierr))
        Iend = Iend - 1
      endif

!     ** Interior grid points
      do i=Istart,Iend-1
        j(1) = i-1
        j(2) = i
        j(3) = i+1
        A(1) = -h/6.0
        A(2) = -2.0*h/3.0
        A(3) = -h/6.0
        PetscCall(MatSetValues(jac,one,[i],three,j,A,INSERT_VALUES,ierr))
      enddo

!     ** Assemble matrix
      PetscCall(MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr))

      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -nep_target 4
!      filter: sed -e "s/[0-9]\.[0-9]*E-[0-9]*/removed/g" -e "s/ Number of NEP iterations = [ 0-9]*/ Number of NEP iterations = /"
!      requires: !single
!
!TEST*/
