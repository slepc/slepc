!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Test Fortran interface of spectrum-slicing Krylov-Schur.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepceps.h>
      use slepceps
      implicit none

#define MAXSUB 16
#define MAXSHI 16

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Mat            A,B,As,Bs,Au
      EPS            eps
      ST             st
      KSP            ksp
      PC             pc
      Vec            v
      PetscScalar    value
      PetscInt       n,m,i,j,k,Istart,Iend
      PetscInt       nev,ncv,mpd,nval
      PetscInt       row,col,nloc,nlocs,mlocs
      PetscInt       II,npart,inertias(MAXSHI)
      PetscBool      flg,lock
      PetscMPIInt    nprc,rank
      PetscReal      int0,int1,keep,subint(MAXSUB)
      PetscReal      shifts(MAXSHI)
      PetscScalar    eval,one,mone,zero
      PetscErrorCode ierr
      MPI_Comm       comm

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      if (ierr/=0) then
        print*,'SlepcInitialize failed'
        stop
      end if
      PetscCallMPIA(MPI_Comm_size(PETSC_COMM_WORLD,nprc,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      n = 35
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      m = n*n
      if (rank==0) then
        write(*,100) n
      end if
 100  format (/'Spectrum-slicing test, n =',I3,' (Fortran)'/)

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))
      PetscCallA(MatCreate(PETSC_COMM_WORLD,B,ierr))
      PetscCallA(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,m,m,ierr))
      PetscCallA(MatSetFromOptions(B,ierr))
      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      do II=Istart,Iend-1
        i = II/n
        j = II-i*n
        value = -1.0
        row = II
        if (i>0) then
          col = II-n
          PetscCallA(MatSetValue(A,row,col,value,INSERT_VALUES,ierr))
        end if
        if (i<n-1) then
          col = II+n
          PetscCallA(MatSetValue(A,row,col,value,INSERT_VALUES,ierr))
        end if
        if (j>0) then
          col = II-1
          PetscCallA(MatSetValue(A,row,col,value,INSERT_VALUES,ierr))
        end if
        if (j<n-1) then
          col = II+1
          PetscCallA(MatSetValue(A,row,col,value,INSERT_VALUES,ierr))
        end if
        col = II
        value = 4.0
        PetscCallA(MatSetValue(A,row,col,value,INSERT_VALUES,ierr))
        value = 2.0
        PetscCallA(MatSetValue(B,row,col,value,INSERT_VALUES,ierr))
      end do
      if (Istart==0) then
        row = 0
        col = 0
        value = 6.0
        PetscCallA(MatSetValue(B,row,col,value,INSERT_VALUES,ierr))
        row = 0
        col = 1
        value = -1.0
        PetscCallA(MatSetValue(B,row,col,value,INSERT_VALUES,ierr))
        row = 1
        col = 0
        value = -1.0
        PetscCallA(MatSetValue(B,row,col,value,INSERT_VALUES,ierr))
        row = 1
        col = 1
        value = 1.0
        PetscCallA(MatSetValue(B,row,col,value,INSERT_VALUES,ierr))
      end if
      PetscCallA(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSCreate(PETSC_COMM_WORLD,eps,ierr))
      PetscCallA(EPSSetOperators(eps,A,B,ierr))
      PetscCallA(EPSSetProblemType(eps,EPS_GHEP,ierr))
      PetscCallA(EPSSetType(eps,EPSKRYLOVSCHUR,ierr))

!     Set interval and other settings for spectrum slicing

      PetscCallA(EPSSetWhichEigenpairs(eps,EPS_ALL,ierr))
      int0 = 1.1
      int1 = 1.3
      PetscCallA(EPSSetInterval(eps,int0,int1,ierr))
      PetscCallA(EPSGetST(eps,st,ierr))
      PetscCallA(STSetType(st,STSINVERT,ierr))
      if (nprc>0) then
        npart = nprc
        PetscCallA(EPSKrylovSchurSetPartitions(eps,npart,ierr))
      end if
      PetscCallA(EPSKrylovSchurGetKSP(eps,ksp,ierr))
      PetscCallA(KSPGetPC(ksp,pc,ierr))
      PetscCallA(KSPSetType(ksp,KSPPREONLY,ierr))
      PetscCallA(PCSetType(pc,PCCHOLESKY,ierr))

!     Test interface functions of Krylov-Schur solver

      PetscCallA(EPSKrylovSchurGetRestart(eps,keep,ierr))
      if (rank==0) then
        write(*,110) keep
      end if
 110  format (' Restart parameter before changing = ',f7.4)
      keep = 0.4
      PetscCallA(EPSKrylovSchurSetRestart(eps,keep,ierr))
      PetscCallA(EPSKrylovSchurGetRestart(eps,keep,ierr))
      if (rank==0) then
        write(*,120) keep
      end if
 120  format (' ... changed to ',f7.4)

      PetscCallA(EPSKrylovSchurGetLocking(eps,lock,ierr))
      if (rank==0) then
        write(*,130) lock
      end if
 130  format (' Locking flag before changing = ',L4)
      PetscCallA(EPSKrylovSchurSetLocking(eps,PETSC_FALSE,ierr))
      PetscCallA(EPSKrylovSchurGetLocking(eps,lock,ierr))
      if (rank==0) then
        write(*,140) lock
      end if
 140  format (' ... changed to ',L4)

      PetscCallA(EPSKrylovSchurGetDimensions(eps,nev,ncv,mpd,ierr))
      if (rank==0) then
        write(*,150) nev,ncv,mpd
      end if
 150  format (' Sub-solve dimensions before changing: nev=',I2,', ncv=',I2,', mpd=',I2)
      nev = 30
      ncv = 60
      mpd = 60
      PetscCallA(EPSKrylovSchurSetDimensions(eps,nev,ncv,mpd,ierr))
      PetscCallA(EPSKrylovSchurGetDimensions(eps,nev,ncv,mpd,ierr))
      if (rank==0) then
        write(*,160) nev,ncv,mpd
      end if
 160  format (' ... changed to: nev=',I2,', ncv=',I2,', mpd=',I2)

      if (nprc>0) then
        PetscCallA(EPSKrylovSchurGetPartitions(eps,npart,ierr))
        if (rank==0) then
          write(*,170) npart
        end if
 170    format (' Using ',I2,' partitions')
        PetscCheckA(npart<=MAXSUB,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,'Too many subintervals')
        subint(1) = int0
        subint(npart+1) = int1
        do i=2,npart
          subint(i) = int0+(i-1)*(int1-int0)/npart
        end do
        PetscCallA(EPSKrylovSchurSetSubintervals(eps,subint,ierr))
        PetscCallA(EPSKrylovSchurGetSubintervals(eps,subint,ierr))
        if (rank==0) then
          write(*,*) 'Using sub-interval separations ='
          do i=2,npart
            write(*,180) subint(i)
          end do
        end if
 180    format (f7.4)
      end if

      PetscCallA(EPSSetFromOptions(eps,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute all eigenvalues in interval and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSSetUp(eps,ierr))
      PetscCallA(EPSKrylovSchurGetInertias(eps,k,PETSC_NULL_REAL_ARRAY,PETSC_NULL_INTEGER_ARRAY,ierr))
      PetscCheckA(k<=MAXSHI,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,'Too many shifts')
      PetscCallA(EPSKrylovSchurGetInertias(eps,k,shifts,inertias,ierr))
      if (rank==0) then
        write(*,*) 'Inertias after EPSSetUp:'
        do i=1,k
          write(*,185) shifts(i),inertias(i)
        end do
      end if
 185  format (' .. ',f4.1,' (',I3,')')

      PetscCallA(EPSSolve(eps,ierr))
      PetscCallA(EPSGetDimensions(eps,nev,ncv,mpd,ierr))
      PetscCallA(EPSGetInterval(eps,int0,int1,ierr))
      if (rank==0) then
        write(*,190) nev,int0,int1
      end if
 190  format (' Found ',I2,' eigenvalues in interval [',f7.4,',',f7.4,']')

      if (nprc>0) then
        PetscCallA(EPSKrylovSchurGetSubcommInfo(eps,k,nval,v,ierr))
        if (rank==0) then
          write(*,200) rank,k,nval
          do i=0,nval-1
            PetscCallA(EPSKrylovSchurGetSubcommPairs(eps,i,eval,v,ierr))
            write(*,210) PetscRealPart(eval)
          end do
        end if
 200    format (' Process ',I2,' has worked in sub-interval ',I2,', containing ',I2,' eigenvalues')
 210    format (f7.4)
        PetscCallA(VecDestroy(v,ierr))

        PetscCallA(EPSKrylovSchurGetSubcommMats(eps,As,Bs,ierr))
        PetscCallA(MatGetLocalSize(A,nloc,PETSC_NULL_INTEGER,ierr))
        PetscCallA(MatGetLocalSize(As,nlocs,mlocs,ierr))
        if (rank==0) then
          write(*,220) rank,nloc,nlocs
        end if
 220    format (' Process ',I2,' owns ',I5,', rows of the global',' matrices, and ',I5,' rows in the subcommunicator')

!       modify A on subcommunicators
        PetscCallA(PetscObjectGetComm(As,comm,ierr))
        PetscCallA(MatCreate(comm,Au,ierr))
        PetscCallA(MatSetSizes(Au,nlocs,mlocs,m,m,ierr))
        PetscCallA(MatSetFromOptions(Au,ierr))
        PetscCallA(MatGetOwnershipRange(Au,Istart,Iend,ierr))
        do II=Istart,Iend-1
          value = 0.5
          PetscCallA(MatSetValue(Au,II,II,value,INSERT_VALUES,ierr))
        end do
        PetscCallA(MatAssemblyBegin(Au,MAT_FINAL_ASSEMBLY,ierr))
        PetscCallA(MatAssemblyEnd(Au,MAT_FINAL_ASSEMBLY,ierr))
        one = 1.0
        mone = -1.0
        zero = 0.0
        PetscCallA(EPSKrylovSchurUpdateSubcommMats(eps,one,mone,Au,zero,zero,PETSC_NULL_MAT,DIFFERENT_NONZERO_PATTERN,PETSC_TRUE,ierr))
        PetscCallA(MatDestroy(Au,ierr))
      end if

      PetscCallA(EPSDestroy(eps,ierr))
      PetscCallA(MatDestroy(A,ierr))
      PetscCallA(MatDestroy(B,ierr))

      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      nsize: 2
!
!TEST*/
