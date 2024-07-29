!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test1f [-help]
!
!  Description: Simple example that tests BV interface functions.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcbv.h>
      use slepcbv
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define KMAX 35

      Vec            t,v
      Mat            Q,M
      BV             X,Y;
      PetscMPIInt    rank
      PetscInt       i,j,n,k,l,izero,ione
      PetscScalar    z(KMAX),val
      PetscScalar, pointer :: qq(:,:)
      PetscScalar    one,mone,two,zero
      PetscReal      nrm
      PetscBool      flg
      PetscErrorCode ierr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      n = 10
      k = 5
      l = 3
      one = 1.0
      mone = -1.0
      two = 2.0
      zero = 0.0
      izero = 0
      ione = 1
      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-k',k,flg,ierr))
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-l',l,flg,ierr))
      if (k .gt. KMAX) then; SETERRA(PETSC_COMM_SELF,1,'Program currently limited to k=35'); endif
      if (rank .eq. 0) then
        write(*,110) k,n
      endif
 110  format (/'Test BV with',I3,' columns of length',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initialize data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create template vector
      PetscCallA(VecCreate(PETSC_COMM_WORLD,t,ierr))
      PetscCallA(VecSetSizes(t,PETSC_DECIDE,n,ierr))
      PetscCallA(VecSetFromOptions(t,ierr))

!     ** Create BV object X
      PetscCallA(BVCreate(PETSC_COMM_WORLD,X,ierr))
      PetscCallA(BVSetSizesFromVec(X,t,k,ierr))
      PetscCallA(BVSetFromOptions(X,ierr))

!     ** Fill X entries
      do j=0,k-1
        PetscCallA(BVGetColumn(X,j,v,ierr))
        PetscCallA(VecSet(v,zero,ierr))
        do i=0,3
          if (i+j<n) then
            val = 3*i+j-2
            PetscCallA(VecSetValue(v,i+j,val,INSERT_VALUES,ierr))
          end if
        end do
        PetscCallA(VecAssemblyBegin(v,ierr))
        PetscCallA(VecAssemblyEnd(v,ierr))
        PetscCallA(BVRestoreColumn(X,j,v,ierr))
      end do

!     ** Create BV object Y
      PetscCallA(BVCreate(PETSC_COMM_WORLD,Y,ierr))
      PetscCallA(BVSetSizesFromVec(Y,t,l,ierr))
      PetscCallA(BVSetFromOptions(Y,ierr))

!     ** Fill Y entries
      do j=0,l-1
        PetscCallA(BVGetColumn(Y,j,v,ierr))
        val = real(j+1)/4.0
        PetscCallA(VecSet(v,val,ierr))
        PetscCallA(BVRestoreColumn(Y,j,v,ierr))
      end do

!     ** Create Mat
      PetscCallA(MatCreateSeqDense(PETSC_COMM_SELF,k,l,PETSC_NULL_SCALAR_ARRAY,Q,ierr))
      PetscCallA(MatDenseGetArray(Q,qq,ierr))
      do i=1,k
        do j=1,l
          if (i<j) then
            qq(i,j) = 2.0
          else
            qq(i,j) = -0.5
          end if
        end do
      end do
      PetscCallA(MatDenseRestoreArray(Q,qq,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Test several operations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Test BVMult
      PetscCallA(BVMult(Y,two,one,X,Q,ierr))

!     ** Test BVMultVec
      PetscCallA(BVGetColumn(Y,izero,v,ierr))
      z(1) = 2.0
      do i=2,k
        z(i) = -0.5*z(i-1)
      end do
      PetscCallA(BVMultVec(X,mone,one,v,z,ierr))
      PetscCallA(BVRestoreColumn(Y,izero,v,ierr))

!     ** Test BVDot
      PetscCallA(MatCreateSeqDense(PETSC_COMM_SELF,l,k,PETSC_NULL_SCALAR_ARRAY,M,ierr))
      PetscCallA(BVDot(X,Y,M,ierr))

!     ** Test BVDotVec
      PetscCallA(BVGetColumn(Y,izero,v,ierr))
      PetscCallA(BVDotVec(X,v,z,ierr))
      PetscCallA(BVRestoreColumn(Y,izero,v,ierr))

!     ** Test BVMultInPlace and BVScale
      PetscCallA(BVMultInPlace(X,Q,ione,l,ierr))
      PetscCallA(BVScale(X,two,ierr))

!     ** Test BVNorm
      PetscCallA(BVNormColumn(X,izero,NORM_2,nrm,ierr))
      if (rank .eq. 0) then
        write(*,120) nrm
      endif
 120  format ('2-Norm of X[0] = ',f8.4)
      PetscCallA(BVNorm(X,NORM_FROBENIUS,nrm,ierr))
      if (rank .eq. 0) then
        write(*,130) nrm
      endif
 130  format ('Frobenius Norm of X = ',f8.4)

!     *** Clean up
      PetscCallA(BVDestroy(X,ierr))
      PetscCallA(BVDestroy(Y,ierr))
      PetscCallA(VecDestroy(t,ierr))
      PetscCallA(MatDestroy(Q,ierr))
      PetscCallA(MatDestroy(M,ierr))
      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      nsize: 1
!      args: -bv_type {{vecs contiguous svec mat}separate output}
!      output_file: output/test1f_1.out
!
!   test:
!      suffix: 2
!      nsize: 2
!      args: -bv_type {{vecs contiguous svec mat}separate output}
!      output_file: output/test1f_1.out
!
!TEST*/
