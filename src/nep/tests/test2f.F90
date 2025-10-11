!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Simple example to test the NEP Fortran interface.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcnep.h>
      use slepcnep
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Mat                A(3),B
      FN                 f(3),g
      NEP                nep
      DS                 ds
      RG                 rg
      PetscReal          tol
      PetscScalar        coeffs(2),tget,val
      PetscInt           n,i,its,Istart,Iend
      PetscInt           nev,ncv,mpd,nterm
      PetscInt           nc,np
      NEPWhich           which
      NEPConvergedReason reason
      NEPType            tname
      NEPRefine          refine
      NEPRefineScheme    rscheme
      NEPConv            conv
      NEPStop            stp
      NEPProblemType     ptype
      MatStructure       mstr
      PetscMPIInt        rank
      PetscErrorCode     ierr
      PetscViewerAndFormat vf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 20
      if (rank==0) then
        write(*,100) n
      end if
 100  format (/'Diagonal Nonlinear Eigenproblem, n =',I3,' (Fortran)')

!     Matrices
      PetscCallA(MatCreate(PETSC_COMM_WORLD,A(1),ierr))
      PetscCallA(MatSetSizes(A(1),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(A(1),ierr))
      PetscCallA(MatGetOwnershipRange(A(1),Istart,Iend,ierr))
      do i=Istart,Iend-1
        val = i+1
        PetscCallA(MatSetValue(A(1),i,i,val,INSERT_VALUES,ierr))
      end do
      PetscCallA(MatAssemblyBegin(A(1),MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A(1),MAT_FINAL_ASSEMBLY,ierr))

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A(2),ierr))
      PetscCallA(MatSetSizes(A(2),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(A(2),ierr))
      PetscCallA(MatGetOwnershipRange(A(2),Istart,Iend,ierr))
      do i=Istart,Iend-1
        val = 1
        PetscCallA(MatSetValue(A(2),i,i,val,INSERT_VALUES,ierr))
      end do
      PetscCallA(MatAssemblyBegin(A(2),MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A(2),MAT_FINAL_ASSEMBLY,ierr))

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A(3),ierr))
      PetscCallA(MatSetSizes(A(3),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(A(3),ierr))
      PetscCallA(MatGetOwnershipRange(A(3),Istart,Iend,ierr))
      do i=Istart,Iend-1
        val = real(n)/real(i+1)
        PetscCallA(MatSetValue(A(3),i,i,val,INSERT_VALUES,ierr))
      end do
      PetscCallA(MatAssemblyBegin(A(3),MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A(3),MAT_FINAL_ASSEMBLY,ierr))

!     Functions: f0=-lambda, f1=1.0, f2=sqrt(lambda)
      PetscCallA(FNCreate(PETSC_COMM_WORLD,f(1),ierr))
      PetscCallA(FNSetType(f(1),FNRATIONAL,ierr))
      nc = 2
      coeffs(1) = -1.0
      coeffs(2) = 0.0
      PetscCallA(FNRationalSetNumerator(f(1),nc,coeffs,ierr))

      PetscCallA(FNCreate(PETSC_COMM_WORLD,f(2),ierr))
      PetscCallA(FNSetType(f(2),FNRATIONAL,ierr))
      nc = 1
      coeffs(1) = 1.0
      PetscCallA(FNRationalSetNumerator(f(2),nc,coeffs,ierr))

      PetscCallA(FNCreate(PETSC_COMM_WORLD,f(3),ierr))
      PetscCallA(FNSetType(f(3),FNSQRT,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create eigensolver and test interface functions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(NEPCreate(PETSC_COMM_WORLD,nep,ierr))
      nterm = 3
      mstr = SAME_NONZERO_PATTERN
      PetscCallA(NEPSetSplitOperator(nep,nterm,A,f,mstr,ierr))
      PetscCallA(NEPGetSplitOperatorInfo(nep,nterm,mstr,ierr))
      if (rank==0) then
        write(*,110) nterm
      end if
 110  format (' Nonlinear function with ',I2,' terms')
      i = 0
      PetscCallA(NEPGetSplitOperatorTerm(nep,i,B,g,ierr))
      PetscCallA(MatView(B,PETSC_NULL_VIEWER,ierr))
      PetscCallA(FNView(g,PETSC_NULL_VIEWER,ierr))

      PetscCallA(NEPSetType(nep,NEPRII,ierr))
      PetscCallA(NEPGetType(nep,tname,ierr))
      if (rank==0) then
        write(*,120) tname
      end if
 120  format (' Type set to ',A)

      PetscCallA(NEPGetProblemType(nep,ptype,ierr))
      if (rank==0) then
        write(*,130) ptype
      end if
 130  format (' Problem type before changing = ',I2)
      PetscCallA(NEPSetProblemType(nep,NEP_RATIONAL,ierr))
      PetscCallA(NEPGetProblemType(nep,ptype,ierr))
      if (rank==0) then
        write(*,140) ptype
      end if
 140  format (' ... changed to ',I2)

      np = 1
      tol = 1e-9
      its = 2
      refine = NEP_REFINE_SIMPLE
      rscheme = NEP_REFINE_SCHEME_EXPLICIT
      PetscCallA(NEPSetRefine(nep,refine,np,tol,its,rscheme,ierr))
      PetscCallA(NEPGetRefine(nep,refine,np,tol,its,rscheme,ierr))
      if (rank==0) then
        write(*,190) refine,tol,its,rscheme
      end if
 190  format (' Refinement: ',I2,', tol=',F12.9,', its=',I2,', scheme=',I2)

      tget = 1.1
      PetscCallA(NEPSetTarget(nep,tget,ierr))
      PetscCallA(NEPGetTarget(nep,tget,ierr))
      PetscCallA(NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE,ierr))
      PetscCallA(NEPGetWhichEigenpairs(nep,which,ierr))
      if (rank==0) then
        write(*,200) which,PetscRealPart(tget)
      end if
 200  format (' Which = ',I2,', target = ',F4.1)

      nev = 1
      ncv = 12
      PetscCallA(NEPSetDimensions(nep,nev,ncv,PETSC_DETERMINE_INTEGER,ierr))
      PetscCallA(NEPGetDimensions(nep,nev,ncv,mpd,ierr))
      if (rank==0) then
        write(*,210) nev,ncv,mpd
      end if
 210  format (' Dimensions: nev=',I2,', ncv=',I2,', mpd=',I2)

      tol = 1.0e-6
      its = 200
      PetscCallA(NEPSetTolerances(nep,tol,its,ierr))
      PetscCallA(NEPGetTolerances(nep,tol,its,ierr))
      if (rank==0) then
        write(*,220) tol,its
      end if
 220  format (' Tolerance =',F9.6,', max_its =',I4)

      PetscCallA(NEPSetConvergenceTest(nep,NEP_CONV_ABS,ierr))
      PetscCallA(NEPGetConvergenceTest(nep,conv,ierr))
      PetscCallA(NEPSetStoppingTest(nep,NEP_STOP_BASIC,ierr))
      PetscCallA(NEPGetStoppingTest(nep,stp,ierr))
      if (rank==0) then
        write(*,230) conv,stp
      end if
 230  format (' Convergence test =',I2,', stopping test =',I2)

      PetscCallA(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr))
      PetscCallA(NEPMonitorSet(nep,NEPMONITORFIRST,vf,PetscViewerAndFormatDestroy,ierr))
      PetscCallA(NEPMonitorConvergedCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,PETSC_NULL_VEC,vf,ierr))
      PetscCallA(NEPMonitorSet(nep,NEPMONITORCONVERGED,vf,NEPMonitorConvergedDestroy,ierr))
      PetscCallA(NEPMonitorCancel(nep,ierr))

      PetscCallA(NEPGetDS(nep,ds,ierr))
      PetscCallA(DSView(ds,PETSC_NULL_VIEWER,ierr))
      PetscCallA(NEPSetFromOptions(nep,ierr))

      PetscCallA(NEPGetRG(nep,rg,ierr))
      PetscCallA(RGView(rg,PETSC_NULL_VIEWER,ierr))

      PetscCallA(NEPSolve(nep,ierr))
      PetscCallA(NEPGetConvergedReason(nep,reason,ierr))
      if (rank==0) then
        write(*,240) reason
      end if
 240  format (' Finished - converged reason =',I2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PetscCallA(NEPErrorView(nep,NEP_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      PetscCallA(NEPDestroy(nep,ierr))
      PetscCallA(MatDestroy(A(1),ierr))
      PetscCallA(MatDestroy(A(2),ierr))
      PetscCallA(MatDestroy(A(3),ierr))
      PetscCallA(FNDestroy(f(1),ierr))
      PetscCallA(FNDestroy(f(2),ierr))
      PetscCallA(FNDestroy(f(3),ierr))

      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      requires: !single
!
!TEST*/
