!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Simple example to test the EPS Fortran interface.
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
      Mat                A,B
      EPS                eps
      ST                 st
      KSP                ksp
      DS                 ds
      PetscReal          cut,tol,tolabs
      PetscScalar        tget,value
      PetscInt           n,i,its,Istart,Iend
      PetscInt           nev,ncv,mpd
      PetscBool          flg
      EPSConvergedReason reason
      EPSType            tname
      EPSExtraction      extr
      EPSBalance         bal
      EPSWhich           which
      EPSConv            conv
      EPSStop            stp
      EPSProblemType     ptype
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
 100  format (/'Diagonal Eigenproblem, n =',I3,' (Fortran)')

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))
      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      do i=Istart,Iend-1
        value = i+1
        PetscCallA(MatSetValue(A,i,i,value,INSERT_VALUES,ierr))
      end do
      PetscCallA(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create eigensolver and test interface functions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSCreate(PETSC_COMM_WORLD,eps,ierr))
      PetscCallA(EPSSetOperators(eps,A,PETSC_NULL_MAT,ierr))
      PetscCallA(EPSGetOperators(eps,B,PETSC_NULL_MAT,ierr))
      PetscCallA(MatView(B,PETSC_NULL_VIEWER,ierr))

      PetscCallA(EPSSetType(eps,EPSKRYLOVSCHUR,ierr))
      PetscCallA(EPSGetType(eps,tname,ierr))
      if (rank .eq. 0) then
        write(*,110) tname
      end if
 110  format (' Type set to ',A)

      PetscCallA(EPSGetProblemType(eps,ptype,ierr))
      if (rank .eq. 0) then
        write(*,120) ptype
      end if
 120  format (' Problem type before changing = ',I2)
      PetscCallA(EPSSetProblemType(eps,EPS_HEP,ierr))
      PetscCallA(EPSGetProblemType(eps,ptype,ierr))
      if (rank .eq. 0) then
        write(*,130) ptype
      end if
 130  format (' ... changed to ',I2)
      PetscCallA(EPSIsGeneralized(eps,flg,ierr))
      if (flg .and. rank .eq. 0) then
        write(*,*) 'generalized'
      end if
      PetscCallA(EPSIsHermitian(eps,flg,ierr))
      if (flg .and. rank .eq. 0) then
        write(*,*) 'hermitian'
      end if
      PetscCallA(EPSIsPositive(eps,flg,ierr))
      if (flg .and. rank .eq. 0) then
        write(*,*) 'positive'
      end if

      PetscCallA(EPSGetExtraction(eps,extr,ierr))
      if (rank .eq. 0) then
        write(*,140) extr
      end if
 140  format (' Extraction before changing = ',I2)
      PetscCallA(EPSSetExtraction(eps,EPS_HARMONIC,ierr))
      PetscCallA(EPSGetExtraction(eps,extr,ierr))
      if (rank .eq. 0) then
        write(*,150) extr
      end if
 150  format (' ... changed to ',I2)

      its = 8
      cut = 2.0e-6
      bal = EPS_BALANCE_ONESIDE
      PetscCallA(EPSSetBalance(eps,bal,its,cut,ierr))
      PetscCallA(EPSGetBalance(eps,bal,its,cut,ierr))
      if (rank .eq. 0) then
        write(*,160) bal,its,cut
      end if
 160  format (' Balance: ',I2,', its=',I2,', cutoff=',F9.6)

      tget = 4.8
      PetscCallA(EPSSetTarget(eps,tget,ierr))
      PetscCallA(EPSGetTarget(eps,tget,ierr))
      PetscCallA(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE,ierr))
      PetscCallA(EPSGetWhichEigenpairs(eps,which,ierr))
      if (rank .eq. 0) then
        write(*,170) which,PetscRealPart(tget)
      end if
 170  format (' Which = ',I2,', target = ',F4.1)

      nev = 4
      PetscCallA(EPSSetDimensions(eps,nev,PETSC_DETERMINE_INTEGER,PETSC_DETERMINE_INTEGER,ierr))
      PetscCallA(EPSGetDimensions(eps,nev,ncv,mpd,ierr))
      if (rank .eq. 0) then
        write(*,180) nev,ncv,mpd
      end if
 180  format (' Dimensions: nev=',I2,', ncv=',I2,', mpd=',I2)

      tol = 2.2e-4
      its = 200
      PetscCallA(EPSSetTolerances(eps,tol,its,ierr))
      PetscCallA(EPSGetTolerances(eps,tol,its,ierr))
      if (rank .eq. 0) then
        write(*,190) tol,its
      end if
 190  format (' Tolerance =',F8.5,', max_its =',I4)

      PetscCallA(EPSSetConvergenceTest(eps,EPS_CONV_ABS,ierr))
      PetscCallA(EPSGetConvergenceTest(eps,conv,ierr))
      PetscCallA(EPSSetStoppingTest(eps,EPS_STOP_BASIC,ierr))
      PetscCallA(EPSGetStoppingTest(eps,stp,ierr))
      if (rank .eq. 0) then
        write(*,200) conv,stp
      end if
 200  format (' Convergence test =',I2,', stopping test =',I2)

      PetscCallA(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr))
      PetscCallA(EPSMonitorSet(eps,EPSMONITORFIRST,vf,PetscViewerAndFormatDestroy,ierr))
      PetscCallA(EPSMonitorConvergedCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,PETSC_NULL_VEC,vf,ierr))
      PetscCallA(EPSMonitorSet(eps,EPSMONITORCONVERGED,vf,EPSMonitorConvergedDestroy,ierr))
      PetscCallA(EPSMonitorCancel(eps,ierr))

      PetscCallA(EPSGetST(eps,st,ierr))
      PetscCallA(STGetKSP(st,ksp,ierr))
      tol = 1.e-8
      tolabs = 1.e-35
      PetscCallA(KSPSetTolerances(ksp,tol,tolabs,PETSC_CURRENT_REAL,PETSC_CURRENT_INTEGER,ierr))
      PetscCallA(STView(st,PETSC_NULL_VIEWER,ierr))
      PetscCallA(EPSGetDS(eps,ds,ierr))
      PetscCallA(DSView(ds,PETSC_NULL_VIEWER,ierr))

      PetscCallA(EPSSetFromOptions(eps,ierr))
      PetscCallA(EPSSolve(eps,ierr))
      PetscCallA(EPSGetConvergedReason(eps,reason,ierr))
      PetscCallA(EPSGetIterationNumber(eps,its,ierr))
      if (rank .eq. 0) then
        write(*,210) reason,its
      end if
 210  format (' Finished - converged reason =',I2,', its=',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      PetscCallA(EPSDestroy(eps,ierr))
      PetscCallA(MatDestroy(A,ierr))

      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -eps_ncv 14
!      filter: sed -e "s/00001/00000/" | sed -e "s/4.99999/5.00000/" | sed -e "s/5.99999/6.00000/"
!
!TEST*/
