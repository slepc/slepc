!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex54f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Illustrates use of shell matrices in callback interface from Fortran.
!  Similar to ex21.c. This one solves a simple diagonal linear eigenproblem as a NEP.
!
!  The command line options are:
!    -n <n>, where <n> = matrix dimension

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    Modules needed to pass and get the context to/from the Mat
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   module shell_ctx
#include <petsc/finclude/petscmat.h>
      use petscmat
      implicit none
      type :: MatCtx
         PetscScalar :: lambda
      end type MatCtx
   end module shell_ctx

   module shell_ctx_interfaces
      use shell_ctx
      implicit none

      interface MatCreateShell
        subroutine MatCreateShell(comm,mloc,nloc,m,n,ctx,mat,ierr)
          use shell_ctx
          MPI_Comm       :: comm
          PetscInt       :: mloc,nloc,m,n
          type(MatCtx)   :: ctx
          Mat            :: mat
          PetscErrorCode :: ierr
        end subroutine MatCreateShell
      end interface MatCreateShell

      interface MatShellSetContext
        subroutine MatShellSetContext(mat,ctx,ierr)
          use shell_ctx
          Mat            :: mat
          type(MatCtx)   :: ctx
          PetscErrorCode :: ierr
        end subroutine MatShellSetContext
      end interface MatShellSetContext

      interface MatShellGetContext
        subroutine MatShellGetContext(mat,ctx,ierr)
          use shell_ctx
          Mat                   :: mat
          type(MatCtx), pointer :: ctx
          PetscErrorCode        :: ierr
        end subroutine matShellGetContext
      end interface MatShellGetContext

   end module shell_ctx_interfaces

!=================================================================================================

   program main
#include <slepc/finclude/slepcnep.h>
      use slepcnep
      use shell_ctx_interfaces
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      NEP            :: nep
      Mat            :: A,B
      PetscInt       :: n=400,nev=3,nconv
      PetscErrorCode :: ierr
      PetscScalar    :: sigma
      PetscBool      :: flg,terse
      PetscMPIInt    :: rank
      type(MatCtx)   :: ctxA,ctxB

      external MyNEPFunction,MyNEPJacobian,MatMult_A,MatDuplicate_A,MatMult_B

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      if (rank .eq. 0) then
         write(*,'(/A,I4)') 'Nonlinear eigenproblem with shell matrices, n =',n
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create matrix-free operators for A and B corresponding to T and T'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ctxA%lambda = 0.0
      PetscCallA(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,ctxA,A,ierr))
      PetscCallA(MatShellSetOperation(A,MATOP_MULT,MatMult_A,ierr))
      PetscCallA(MatShellSetOperation(A,MATOP_DUPLICATE,MatDuplicate_A,ierr))

      ctxB%lambda = 0.0   ! unused
      PetscCallA(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,ctxB,B,ierr))
      PetscCallA(MatShellSetOperation(B,MATOP_MULT,MatMult_B,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(NEPCreate(PETSC_COMM_WORLD,nep,ierr))
      PetscCallA(NEPSetFunction(nep,A,A,MyNEPFunction,PETSC_NULL_INTEGER,ierr))
      PetscCallA(NEPSetJacobian(nep,B,MyNEPJacobian,PETSC_NULL_INTEGER,ierr))
      PetscCallA(NEPSetDimensions(nep,nev,PETSC_DETERMINE_INTEGER,PETSC_DETERMINE_INTEGER,ierr))
      sigma = 1.05
      PetscCallA(NEPSetTarget(nep,sigma,ierr))
      PetscCallA(NEPSetFromOptions(nep,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem, display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(NEPSolve(nep, ierr))

      PetscCallA(NEPGetConverged(nep,nconv,ierr))
      if (rank .eq. 0) then
         write(*,'(A,I2/)') ' Number of converged approximate eigenpairs:',nconv
      endif
!
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
      PetscCallA(MatDestroy(A,ierr))
      PetscCallA(MatDestroy(B,ierr))
      PetscCallA(SlepcFinalize(ierr))

   end program main

! --------------------------------------------------------------
!
!  MyNEPFunction - Computes Function matrix  T(lambda)
!
   subroutine MyNEPFunction(nep,lambda,T,P,ctx,ierr)
      use slepcnep
      use shell_ctx_interfaces
      implicit none

      NEP                   :: nep
      PetscScalar           :: lambda
      Mat                   :: T,P
      PetscInt              :: ctx
      PetscErrorCode        :: ierr
      type(MatCtx), pointer :: ctxT

      PetscCall(MatShellGetContext(T,ctxT,ierr))
      ctxT%lambda = lambda

   end subroutine MyNEPFunction

! --------------------------------------------------------------
!
!  MyNEPJacobian - Computes Jacobian matrix  T'(lambda)
!
   subroutine MyNEPJacobian(nep,lambda,T,ctx,ierr)
      use slepcnep
      use shell_ctx_interfaces
      implicit none

      NEP                   :: nep
      PetscScalar           :: lambda
      Mat                   :: T
      PetscInt              :: ctx
      PetscErrorCode        :: ierr
      type(MatCtx), pointer :: ctxT

      PetscCall(MatShellGetContext(T,ctxT,ierr))
      ctxT%lambda = lambda

   end subroutine MyNEPJacobian

! --------------------------------------------------------------
!
!  MatMult_A - Shell matrix operation, multiples y=A*x
!  Here A=(D-lambda*I) where D is a diagonal matrix
!
   subroutine MatMult_A(A,x,y,ierr)
      use shell_ctx_interfaces
      implicit none

      Mat                   :: A
      Vec                   :: x,y
      PetscErrorCode        :: ierr
      PetscInt              :: i,istart,iend
      PetscScalar           :: val
      type(MatCtx), pointer :: ctxA
!
      PetscCall(VecGetOwnershipRange(x,istart,iend,ierr))
      do i=istart,iend-1
         val = i+1
         val = 1.0/val
         PetscCall(VecSetValue(y,i,val,INSERT_VALUES,ierr))
      end do
      PetscCall(VecAssemblyBegin(y,ierr))
      PetscCall(VecAssemblyEnd(y,ierr))
      PetscCall(VecPointwiseMult(y,y,x,ierr))
      PetscCall(MatShellGetContext(A,ctxA,ierr))
      PetscCall(VecAXPY(y,-ctxA%lambda,x,ierr))

   end subroutine MatMult_A

! --------------------------------------------------------------
!
!  MatDuplicate_A - Shell matrix operation, duplicates A
!
   subroutine MatDuplicate_A(A,opt,M,ierr)
      use shell_ctx_interfaces
      implicit none

      Mat                   :: A,M
      MatDuplicateOption    :: opt
      PetscErrorCode        :: ierr
      PetscInt              :: ml,nl
      type(MatCtx), pointer :: ctxM,ctxA

      external MatMult_A,MatDestroy_A

      PetscCall(MatGetLocalSize(A,ml,nl,ierr));
      PetscCall(MatShellGetContext(A,ctxA,ierr))
      allocate(ctxM)
      ctxM%lambda = ctxA%lambda
      PetscCall(MatCreateShell(PETSC_COMM_WORLD,ml,nl,PETSC_DETERMINE,PETSC_DETERMINE,ctxM,M,ierr))
      PetscCall(MatShellSetOperation(M,MATOP_MULT,MatMult_A,ierr))
      PetscCall(MatShellSetOperation(M,MATOP_DESTROY,MatDestroy_A,ierr))

   end subroutine MatDuplicate_A

! --------------------------------------------------------------
!
!  MatDestroy_A - Shell matrix operation, destroys A
!
   subroutine MatDestroy_A(A,ierr)
      use shell_ctx_interfaces
      implicit none

      Mat                   :: A
      PetscErrorCode        :: ierr
      type(MatCtx), pointer :: ctxA

      PetscCall(MatShellGetContext(A,ctxA,ierr))
      deallocate(ctxA)

   end subroutine MatDestroy_A

! --------------------------------------------------------------
!
!  MatMult_B - Shell matrix operation, multiples y=B*x
!  Here B=-I
!
   subroutine MatMult_B(B,x,y,ierr)
      use petscmat
      implicit none

      Mat            :: B
      Vec            :: x
      Vec            :: y
      PetscErrorCode :: ierr
      PetscScalar    :: mone

      PetscCall(VecCopy(x,y,ierr))
      mone = -1.0
      PetscCall(VecScale(y,mone,ierr))

   end subroutine MatMult_B

!/*TEST
!
!   testset:
!      args: -terse
!      output_file: output/ex54f_1.out
!      filter: grep -v approximate | sed -e "s/[+-]0\.0*i//g"
!      test:
!         suffix: 1_slp
!         args: -nep_type slp -nep_slp_ksp_type gmres -nep_slp_pc_type none
!         requires: double
!      test:
!         suffix: 1_nleigs
!         args: -nep_type nleigs -rg_interval_endpoints 0.2,1.1 -nep_nleigs_ksp_type gmres -nep_nleigs_pc_type none
!         requires: !complex
!      test:
!         suffix: 1_nleigs_complex
!         args: -nep_type nleigs -rg_interval_endpoints 0.2,1.1,-.1,.1 -nep_nleigs_ksp_type gmres -nep_nleigs_pc_type none
!         requires: complex
!
!TEST*/
