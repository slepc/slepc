!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Simple example to test the PEP Fortran interface.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcpep.h>
      use slepcpep
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Mat                A(3),B
      PEP                pep
      ST                 st
      KSP                ksp
      DS                 ds
      PetscReal          tol,tolabs,alpha,lambda
      PetscScalar        tget,val
      PetscInt           n,i,its,Istart,Iend
      PetscInt           nev,ncv,mpd,nmat,np
      PEPWhich           which
      PEPConvergedReason reason
      PEPType            tname
      PEPExtract         extr
      PEPBasis           basis
      PEPScale           scal
      PEPRefine          refine
      PEPRefineScheme    rscheme
      PEPConv            conv
      PEPStop            stp
      PEPProblemType     ptype
      PetscMPIInt        rank
      PetscErrorCode     ierr
      PetscViewerAndFormat vf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 20
      if (rank .eq. 0) then
        write(*,100) n
      end if
 100  format (/'Diagonal Quadratic Eigenproblem, n =',I3,' (Fortran)')

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

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create eigensolver and test interface functions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(PEPCreate(PETSC_COMM_WORLD,pep,ierr))
      nmat = 3
      PetscCallA(PEPSetOperators(pep,nmat,A,ierr))
      PetscCallA(PEPGetNumMatrices(pep,nmat,ierr))
      if (rank .eq. 0) then
        write(*,110) nmat-1
      end if
 110  format (' Polynomial of degree ',I2)
      i = 0
      PetscCallA(PEPGetOperators(pep,i,B,ierr))
      PetscCallA(MatView(B,PETSC_NULL_VIEWER,ierr))

      PetscCallA(PEPSetType(pep,PEPTOAR,ierr))
      PetscCallA(PEPGetType(pep,tname,ierr))
      if (rank .eq. 0) then
        write(*,120) tname
      end if
 120  format (' Type set to ',A)

      PetscCallA(PEPGetProblemType(pep,ptype,ierr))
      if (rank .eq. 0) then
        write(*,130) ptype
      end if
 130  format (' Problem type before changing = ',I2)
      PetscCallA(PEPSetProblemType(pep,PEP_HERMITIAN,ierr))
      PetscCallA(PEPGetProblemType(pep,ptype,ierr))
      if (rank .eq. 0) then
        write(*,140) ptype
      end if
 140  format (' ... changed to ',I2)

      PetscCallA(PEPGetExtract(pep,extr,ierr))
      if (rank .eq. 0) then
        write(*,150) extr
      end if
 150  format (' Extraction before changing = ',I2)
      PetscCallA(PEPSetExtract(pep,PEP_EXTRACT_STRUCTURED,ierr))
      PetscCallA(PEPGetExtract(pep,extr,ierr))
      if (rank .eq. 0) then
        write(*,160) extr
      end if
 160  format (' ... changed to ',I2)

      alpha = .1
      its = 5
      lambda = 1.
      scal = PEP_SCALE_SCALAR
      PetscCallA(PEPSetScale(pep,scal,alpha,PETSC_NULL_VEC,PETSC_NULL_VEC,its,lambda,ierr))
      PetscCallA(PEPGetScale(pep,scal,alpha,PETSC_NULL_VEC,PETSC_NULL_VEC,its,lambda,ierr))
      if (rank .eq. 0) then
        write(*,170) scal,alpha,its
      end if
 170  format (' Scaling: ',I2,', alpha=',F7.4,', its=',I2)

      basis = PEP_BASIS_CHEBYSHEV1
      PetscCallA(PEPSetBasis(pep,basis,ierr))
      PetscCallA(PEPGetBasis(pep,basis,ierr))
      if (rank .eq. 0) then
        write(*,180) basis
      end if
 180  format (' Polynomial basis: ',I2)

      np = 1
      tol = 1e-9
      its = 2
      refine = PEP_REFINE_SIMPLE
      rscheme = PEP_REFINE_SCHEME_SCHUR
      PetscCallA(PEPSetRefine(pep,refine,np,tol,its,rscheme,ierr))
      PetscCallA(PEPGetRefine(pep,refine,np,tol,its,rscheme,ierr))
      if (rank .eq. 0) then
        write(*,190) refine,tol,its,rscheme
      end if
 190  format (' Refinement: ',I2,', tol=',F7.4,', its=',I2,', schem=',I2)

      tget = 4.8
      PetscCallA(PEPSetTarget(pep,tget,ierr))
      PetscCallA(PEPGetTarget(pep,tget,ierr))
      PetscCallA(PEPSetWhichEigenpairs(pep,PEP_TARGET_MAGNITUDE,ierr))
      PetscCallA(PEPGetWhichEigenpairs(pep,which,ierr))
      if (rank .eq. 0) then
        write(*,200) which,PetscRealPart(tget)
      end if
 200  format (' Which = ',I2,', target = ',F4.1)

      nev = 4
      PetscCallA(PEPSetDimensions(pep,nev,PETSC_DETERMINE_INTEGER,PETSC_DETERMINE_INTEGER,ierr))
      PetscCallA(PEPGetDimensions(pep,nev,ncv,mpd,ierr))
      if (rank .eq. 0) then
        write(*,210) nev,ncv,mpd
      end if
 210  format (' Dimensions: nev=',I2,', ncv=',I2,', mpd=',I2)

      tol = 2.2e-4
      its = 200
      PetscCallA(PEPSetTolerances(pep,tol,its,ierr))
      PetscCallA(PEPGetTolerances(pep,tol,its,ierr))
      if (rank .eq. 0) then
        write(*,220) tol,its
      end if
 220  format (' Tolerance =',F8.5,', max_its =',I4)

      PetscCallA(PEPSetConvergenceTest(pep,PEP_CONV_ABS,ierr))
      PetscCallA(PEPGetConvergenceTest(pep,conv,ierr))
      PetscCallA(PEPSetStoppingTest(pep,PEP_STOP_BASIC,ierr))
      PetscCallA(PEPGetStoppingTest(pep,stp,ierr))
      if (rank .eq. 0) then
        write(*,230) conv,stp
      end if
 230  format (' Convergence test =',I2,', stopping test =',I2)

      PetscCallA(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr))
      PetscCallA(PEPMonitorSet(pep,PEPMONITORFIRST,vf,PetscViewerAndFormatDestroy,ierr))
      PetscCallA(PEPMonitorConvergedCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,PETSC_NULL_VEC,vf,ierr))
      PetscCallA(PEPMonitorSet(pep,PEPMONITORCONVERGED,vf,PEPMonitorConvergedDestroy,ierr))
      PetscCallA(PEPMonitorCancel(pep,ierr))

      PetscCallA(PEPGetST(pep,st,ierr))
      PetscCallA(STGetKSP(st,ksp,ierr))
      tol = 1.e-8
      tolabs = 1.e-35
      PetscCallA(KSPSetTolerances(ksp,tol,tolabs,PETSC_CURRENT_REAL,PETSC_CURRENT_INTEGER,ierr))
      PetscCallA(STView(st,PETSC_NULL_VIEWER,ierr))
      PetscCallA(PEPGetDS(pep,ds,ierr))
      PetscCallA(DSView(ds,PETSC_NULL_VIEWER,ierr))

      PetscCallA(PEPSetFromOptions(pep,ierr))
      PetscCallA(PEPSolve(pep,ierr))
      PetscCallA(PEPGetConvergedReason(pep,reason,ierr))
      if (rank .eq. 0) then
        write(*,240) reason
      end if
 240  format (' Finished - converged reason =',I2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PetscCallA(PEPErrorView(pep,PEP_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      PetscCallA(PEPDestroy(pep,ierr))
      PetscCallA(MatDestroy(A(1),ierr))
      PetscCallA(MatDestroy(A(2),ierr))
      PetscCallA(MatDestroy(A(3),ierr))

      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -pep_tol 1e-6 -pep_ncv 22
!      filter: sed -e "s/[+-]0\.0*i//g"
!
!TEST*/
