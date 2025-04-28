!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex27f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Simple NLEIGS example. Fortran90 equivalent of ex27.c
!
!  The command line options are:
!    -n <n>, where <n> = matrix dimension
!
! ----------------------------------------------------------------------
!   Solve T(lambda)x=0 using NLEIGS solver
!      with T(lambda) = -D+sqrt(lambda)*I
!      where D is the Laplacian operator in 1 dimension
!      and with the interpolation interval [.01,16]
! ----------------------------------------------------------------------
!
PROGRAM main
#include <slepc/finclude/slepcnep.h>
  USE slepcnep
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NEP                :: nep
  Mat                :: A(2),F,J
  NEPType            :: ntype
  PetscInt           :: n=100,nev,Istart,Iend,i,col,one,two,three
  PetscErrorCode     :: ierr
  PetscBool          :: terse,flg,split=PETSC_TRUE
  PetscReal          :: ia,ib,ic,id
  RG                 :: rg
  FN                 :: fn(2)
  PetscScalar        :: coeffs,sigma,done
  CHARACTER(LEN=128) :: string

  ! NOTE: Any user-defined Fortran routines (such as ComputeSingularities)
  !       MUST be declared as external.
  external ComputeSingularities, FormFunction, FormJacobian

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-n",n,flg,ierr))
  PetscCallA(PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-split",split,flg,ierr))
  if (split) then
     write(string,*) 'Square root eigenproblem, n=',n,' (split-form)\n'
  else
     write(string,*) 'Square root eigenproblem, n=',n,'\n'
  end if
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,trim(string),ierr))
  done  = 1.0
  one   = 1
  two   = 2
  three = 3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create nonlinear eigensolver context and set options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(NEPCreate(PETSC_COMM_WORLD,nep,ierr))
  PetscCallA(NEPSetType(nep,NEPNLEIGS,ierr))
  PetscCallA(NEPNLEIGSSetSingularitiesFunction(nep,ComputeSingularities,0,ierr))
  PetscCallA(NEPGetRG(nep,rg,ierr))
  PetscCallA(RGSetType(rg,RGINTERVAL,ierr))
  ia = 0.01
  ib = 16.0
#if defined(PETSC_USE_COMPLEX)
  ic = -0.001
  id = 0.001
#else
  ic = 0.0
  id = 0.0
#endif
  PetscCallA(RGIntervalSetEndpoints(rg,ia,ib,ic,id,ierr))
  sigma = 1.1
  PetscCallA(NEPSetTarget(nep,sigma,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Define the nonlinear problem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (split) then
     ! ** Create matrices for the split form
     PetscCallA(MatCreate(PETSC_COMM_WORLD,A(1),ierr))
     PetscCallA(MatSetSizes(A(1),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
     PetscCallA(MatSetFromOptions(A(1),ierr))
     PetscCallA(MatGetOwnershipRange(A(1),Istart,Iend,ierr))
     coeffs = -2.0
     do i=Istart,Iend-1
        if (i.gt.0) then
           col = i-1
           PetscCallA(MatSetValue(A(1),i,col,done,INSERT_VALUES,ierr))
        end if
        if (i.lt.n-1) then
           col = i+1
           PetscCallA(MatSetValue(A(1),i,col,done,INSERT_VALUES,ierr))
        end if
        PetscCallA(MatSetValue(A(1),i,i,coeffs,INSERT_VALUES,ierr))
     end do
     PetscCallA(MatAssemblyBegin(A(1),MAT_FINAL_ASSEMBLY,ierr))
     PetscCallA(MatAssemblyEnd(A(1),MAT_FINAL_ASSEMBLY,ierr))

     PetscCallA(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,done,A(2),ierr))

     ! ** Define functions for the split form
     PetscCallA(FNCreate(PETSC_COMM_WORLD,fn(1),ierr))
     PetscCallA(FNSetType(fn(1),FNRATIONAL,ierr))
     PetscCallA(FNRationalSetNumerator(fn(1),one,[done],ierr))
     PetscCallA(FNCreate(PETSC_COMM_WORLD,fn(2),ierr))
     PetscCallA(FNSetType(fn(2),FNSQRT,ierr))
     PetscCallA(NEPSetSplitOperator(nep,two,A,fn,SUBSET_NONZERO_PATTERN,ierr))
  else
    ! ** Callback form: create matrix and set Function evaluation routine
    PetscCallA(MatCreate(PETSC_COMM_WORLD,F,ierr))
    PetscCallA(MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
    PetscCallA(MatSetFromOptions(F,ierr))
    PetscCallA(MatSeqAIJSetPreallocation(F,three,PETSC_NULL_INTEGER_ARRAY,ierr))
    PetscCallA(MatMPIAIJSetPreallocation(F,three,PETSC_NULL_INTEGER_ARRAY,one,PETSC_NULL_INTEGER_ARRAY,ierr))
    PetscCallA(NEPSetFunction(nep,F,F,FormFunction,PETSC_NULL_INTEGER,ierr))

    PetscCallA(MatCreate(PETSC_COMM_WORLD,J,ierr))
    PetscCallA(MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
    PetscCallA(MatSetFromOptions(J,ierr))
    PetscCallA(MatSeqAIJSetPreallocation(J,one,PETSC_NULL_INTEGER_ARRAY,ierr))
    PetscCallA(MatMPIAIJSetPreallocation(J,one,PETSC_NULL_INTEGER_ARRAY,one,PETSC_NULL_INTEGER_ARRAY,ierr))
    PetscCallA(NEPSetJacobian(nep,J,FormJacobian,PETSC_NULL_INTEGER,ierr))
  end if

  PetscCallA(NEPSetFromOptions(nep,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(NEPSolve(nep,ierr))
  PetscCallA(NEPGetType(nep,ntype,ierr))
  write(string,*) 'Solution method: ',ntype,'\n'
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,trim(string),ierr))
  PetscCallA(NEPGetDimensions(nep,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
  write(string,*) 'Number of requested eigenvalues:',nev,'\n'
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,trim(string),ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! ** show detailed info unless -terse option is given by user
  PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-terse',terse,ierr))
  if (terse) then
    PetscCallA(NEPErrorView(nep,NEP_ERROR_BACKWARD,PETSC_NULL_VIEWER,ierr))
  else
    PetscCallA(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr))
    PetscCallA(NEPConvergedReasonView(nep,PETSC_VIEWER_STDOUT_WORLD,ierr))
    PetscCallA(NEPErrorView(nep,NEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD,ierr))
    PetscCallA(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr))
  end if

  if (split) then
    PetscCallA(MatDestroy(A(1),ierr))
    PetscCallA(MatDestroy(A(2),ierr))
    PetscCallA(FNDestroy(fn(1),ierr))
    PetscCallA(FNDestroy(fn(2),ierr))
  else
    PetscCallA(MatDestroy(F,ierr))
    PetscCallA(MatDestroy(J,ierr))
  end if
  PetscCallA(NEPDestroy(nep,ierr))
  PetscCallA(SlepcFinalize(ierr))

END PROGRAM main

! --------------------------------------------------------------
!
!   FormFunction - Computes Function matrix  T(lambda)
!
SUBROUTINE FormFunction(nep,lambda,fun,B,ctx,ierr)
#include <slepc/finclude/slepcnep.h>
  use slepcnep
  implicit none

  NEP            :: nep
  PetscScalar    :: lambda,val(0:2),t
  Mat            :: fun,B
  PetscInt       :: ctx,i,n,col(0:2),Istart,Iend,Istart0,Iend0,one,two,three
  PetscErrorCode :: ierr
  PetscBool      :: FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE

  one   = 1
  two   = 2
  three = 3

  ! ** Compute Function entries and insert into matrix
  t = sqrt(lambda)
  PetscCall(MatGetSize(fun,n,PETSC_NULL_INTEGER,ierr))
  PetscCall(MatGetOwnershipRange(fun,Istart,Iend,ierr))
  if (Istart.eq.0) FirstBlock=PETSC_TRUE;
  if (Iend.eq.n) LastBlock=PETSC_TRUE;
  val(0)=1.0; val(1)=t-2.0; val(2)=1.0;

  Istart0 = Istart
  if (FirstBlock) Istart0 = Istart+1
  Iend0 = Iend
  if (LastBlock) Iend0 = Iend-1

  do i=Istart0,Iend0-1
     col(0) = i-1
     col(1) = i
     col(2) = i+1
     PetscCall(MatSetValues(fun,one,[i],three,col,val,INSERT_VALUES,ierr))
  end do

  if (LastBlock) then
     i = n-1
     col(0) = n-2
     col(1) = n-1
     val(0) = 1.0
     val(1) = t-2.0
     PetscCall(MatSetValues(fun,one,[i],two,col,val,INSERT_VALUES,ierr))
  end if

  if (FirstBlock) then
     i = 0
     col(0) = 0
     col(1) = 1
     val(0) = t-2.0
     val(1) = 1.0
     PetscCall(MatSetValues(fun,one,[i],two,col,val,INSERT_VALUES,ierr))
  end if

  ! ** Assemble matrix
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr))
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr))
  PetscCall(MatAssemblyBegin(fun,MAT_FINAL_ASSEMBLY,ierr))
  PetscCall(MatAssemblyEnd(fun,MAT_FINAL_ASSEMBLY,ierr))

END SUBROUTINE FormFunction

! --------------------------------------------------------------
!
!   FormJacobian - Computes Jacobian matrix  T'(lambda)
!
SUBROUTINE FormJacobian(nep,lambda,jac,ctx,ierr)
#include <slepc/finclude/slepcnep.h>
  USE slepcnep
  implicit none

  NEP            :: nep
  PetscScalar    :: lambda,t
  Mat            :: jac
  PetscInt       :: ctx
  PetscErrorCode :: ierr
  Vec            :: d

  PetscCall(MatCreateVecs(jac,d,PETSC_NULL_VEC,ierr))
  t = 0.5/sqrt(lambda)
  PetscCall(VecSet(d,t,ierr))
  PetscCall(MatDiagonalSet(jac,d,INSERT_VALUES,ierr))
  PetscCall(VecDestroy(d,ierr))

END SUBROUTINE FormJacobian

! --------------------------------------------------------------
!
!  ComputeSingularities - This is a user-defined routine to compute maxnp
!  points (at most) in the complex plane where the function T(.) is not analytic.
!
!  In this case, we discretize the singularity region (-inf,0)~(-1e+6,-1e-5)
!
!  Input Parameters:
!    nep   - nonlinear eigensolver context
!    maxnp - on input number of requested points in the discretization (can be set)
!    xi    - computed values of the discretization
!    dummy - optional user-defined monitor context (unused here)
!
SUBROUTINE ComputeSingularities(nep,maxnp,xi,dummy,ierr)
#include <slepc/finclude/slepcnep.h>
  use slepcnep
  implicit none

  NEP            :: nep
  PetscInt       :: maxnp, dummy
  PetscScalar    :: xi(0:maxnp-1)
  PetscErrorCode :: ierr
  PetscReal      :: h
  PetscInt       :: i

  h = 11.0/real(maxnp-1)
  xi(0) = -1e-5
  xi(maxnp-1) = -1e+6
  do i=1,maxnp-2
     xi(i) = -10**(-5+h*i)
  end do
  ierr = 0

END SUBROUTINE ComputeSingularities

!/*TEST
!
!   test:
!      suffix: 1
!      args: -nep_nev 3 -nep_nleigs_interpolation_degree 90 -terse
!      requires: !single
!      filter: sed -e "s/[+-]0\.0*i//g"
!
!   test:
!      suffix: 2
!      args: -split 0 -nep_nev 3 -nep_nleigs_interpolation_degree 90 -terse
!      requires: !single
!      filter: sed -e "s/[+-]0\.0*i//g"
!
!TEST*/
