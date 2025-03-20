!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex23f [-help] [-t <t>] [-m <m>] [SLEPc opts]
!
!  Description: Computes exp(t*A)*v for a matrix from a Markov model.
!  This is the Fortran90 equivalent to ex23.c
!
!  The command line options are:
!    -t <t>, where <t> = time parameter (multiplies the matrix)
!    -m <m>, where <m> = number of grid subdivisions in each dimension
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcmfn.h>
      use slepcmfn
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     A      problem matrix
!     mfn    matrix function solver context

      Mat            A
      MFN            mfn
      FN             f
      PetscReal      tol, norm, cst, pd, pu
      PetscScalar    t, z
      Vec            v, y
      PetscInt       N, m, ncv, maxit, its, ii, jj
      PetscInt       i, j, jmax, ix, Istart, Iend
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg
      MFNConvergedReason reason

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      m = 15
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-m',m,flg,ierr))
      t = 2.0
      PetscCallA(PetscOptionsGetScalar(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-t',t,flg,ierr))
      N = m*(m+1)/2
      if (rank .eq. 0) then
        write(*,100) N, m
      endif
 100  format (/'Markov y=exp(t*A)*e_1, N=',I6,' (m=',I4,')')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the transition probability matrix, A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
      PetscCallA(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr))
      PetscCallA(MatSetFromOptions(A,ierr))
      PetscCallA(MatGetOwnershipRange(A,Istart,Iend,ierr))
      ix = 0
      cst = 0.5/real(m-1)
      do i=1,m
        jmax = m-i+1
        do j=1,jmax
          ix = ix + 1
          ii = ix - 1
          if (ix-1.ge.Istart .and. ix.le.Iend) then
            if (j.ne.jmax) then
              pd = cst*(i+j-1)
              !** north
              if (i.eq.1) then
                z = 2.0*pd
              else
                z = pd
              end if
              PetscCallA(MatSetValue(A,ii,ix,z,INSERT_VALUES,ierr))
              !** east
              if (j.eq.1) then
                z = 2.0*pd
              else
                z = pd
              end if
              jj = ix+jmax-1
              PetscCallA(MatSetValue(A,ii,jj,z,INSERT_VALUES,ierr))
            end if

            !** south
            pu = 0.5 - cst*(i+j-3)
            z = pu
            if (j.gt.1) then
              jj = ix-2
              PetscCallA(MatSetValue(A,ii,jj,z,INSERT_VALUES,ierr))
            end if
            !** west
            if (i.gt.1) then
              jj = ix-jmax-2
              PetscCallA(MatSetValue(A,ii,jj,z,INSERT_VALUES,ierr))
            end if
          end if
        end do
      end do
      PetscCallA(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCallA(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))

!     ** Set v = e_1
      PetscCallA(MatCreateVecs(A,y,v,ierr))
      ii = 0
      z = 1.0
      PetscCallA(VecSetValue(v,ii,z,INSERT_VALUES,ierr))
      PetscCallA(VecAssemblyBegin(v,ierr))
      PetscCallA(VecAssemblyEnd(v,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create matrix function solver context
      PetscCallA(MFNCreate(PETSC_COMM_WORLD,mfn,ierr))

!     ** Set operator matrix, the function to compute, and other options
      PetscCallA(MFNSetOperator(mfn,A,ierr))
      PetscCallA(MFNGetFN(mfn,f,ierr))
      PetscCallA(FNSetType(f,FNEXP,ierr))
      z = 1.0
      PetscCallA(FNSetScale(f,t,z,ierr))
      tol = 0.0000001
      PetscCallA(MFNSetTolerances(mfn,tol,PETSC_CURRENT_INTEGER,ierr))

!     ** Set solver parameters at runtime
      PetscCallA(MFNSetFromOptions(mfn,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the problem, y=exp(t*A)*v
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MFNSolve(mfn,v,y,ierr))
      PetscCallA(MFNGetConvergedReason(mfn,reason,ierr))
      if (reason%v.lt.0) then; SETERRA(PETSC_COMM_WORLD,1,'Solver did not converge'); endif
      PetscCallA(VecNorm(y,NORM_2,norm,ierr))

!     ** Optional: Get some information from the solver and display it
      PetscCallA(MFNGetIterationNumber(mfn,its,ierr))
      if (rank .eq. 0) then
        write(*,120) its
      endif
 120  format (' Number of iterations of the method: ',I4)
      PetscCallA(MFNGetDimensions(mfn,ncv,ierr))
      if (rank .eq. 0) then
        write(*,130) ncv
      endif
 130  format (' Subspace dimension:',I4)
      PetscCallA(MFNGetTolerances(mfn,tol,maxit,ierr))
      if (rank .eq. 0) then
        write(*,140) tol,maxit
      endif
 140  format (' Stopping condition: tol=',f10.7,' maxit=',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (rank .eq. 0) then
        write(*,150) PetscRealPart(t),norm
      endif
 150  format (' Computed vector at time t=',f4.1,' has norm ',f8.5)

      PetscCallA(MFNDestroy(mfn,ierr))
      PetscCallA(MatDestroy(A,ierr))
      PetscCallA(VecDestroy(v,ierr))
      PetscCallA(VecDestroy(y,ierr))
      PetscCallA(SlepcFinalize(ierr))
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -mfn_ncv 6
!
!TEST*/
