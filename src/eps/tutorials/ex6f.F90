!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex6f [-help] [-m <m>] [all SLEPc options]
!
!  Description: Eigensystem from the Ising model for ferromagnetic materials.
!  Information about the model can be found at the following
!  site https://math.nist.gov/MatrixMarket/data/NEP
!
!  The command line options are:
!    -m <m>, where <m> is the number of 2x2 blocks, i.e. matrix size N=2*m
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
      PetscReal      tol
      PetscInt       N, m
      PetscInt       nev, maxit, its
      PetscMPIInt    sz, rank
      PetscErrorCode ierr
      PetscBool      flg, terse

!     This is the routine to use for matrix-free approach
!
      external MatMult_Ising, MatMultTranspose_Ising

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))
#if defined(PETSC_USE_COMPLEX)
      SETERRA(PETSC_COMM_SELF,PETSC_ERR_SUP,'This example requires real numbers')
#endif
      PetscCallMPIA(MPI_Comm_size(PETSC_COMM_WORLD,sz,ierr))
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      if (sz .ne. 1) then; SETERRA(PETSC_COMM_SELF,PETSC_ERR_WRONG_MPI_SIZE,'This is a uniprocessor example only!'); endif
      m = 30
      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-m',m,flg,ierr))
      N = 2*m

      if (rank .eq. 0) then
        write(*,'(/A,I6,A/)') 'Ising Model Eigenproblem, m=',m,', (N=2*m)'
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Register the matrix-vector subroutine for the operator that defines
!     the eigensystem, Ax=kx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,PETSC_NULL_INTEGER,A,ierr))
      PetscCallA(MatShellSetOperation(A,MATOP_MULT,MatMult_Ising,ierr))
      PetscCallA(MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,MatMultTranspose_Ising,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      PetscCallA(EPSCreate(PETSC_COMM_WORLD,eps,ierr))

!     ** Set operators. In this case, it is a standard eigenvalue problem
      PetscCallA(EPSSetOperators(eps,A,PETSC_NULL_MAT,ierr))
      PetscCallA(EPSSetProblemType(eps,EPS_NHEP,ierr))

!     ** Set solver parameters at runtime
      PetscCallA(EPSSetFromOptions(eps,ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      PetscCallA(EPSSolve(eps,ierr))
      PetscCallA(EPSGetIterationNumber(eps,its,ierr))
      if (rank .eq. 0) then
        write(*,'(A,I4)') ' Number of iterations of the method: ',its
      endif

!     ** Optional: Get some information from the solver and display it
      PetscCallA(EPSGetType(eps,tname,ierr))
      if (rank .eq. 0) then
        write(*,'(A,A)') ' Solution method: ', tname
      endif
      PetscCallA(EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
      if (rank .eq. 0) then
        write(*,'(A,I2)') ' Number of requested eigenvalues:',nev
      endif
      PetscCallA(EPSGetTolerances(eps,tol,maxit,ierr))
      if (rank .eq. 0) then
        write(*,'(A,1PE11.4,A,I6)') ' Stopping condition: tol=',tol,', maxit=', maxit
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** show detailed info unless -terse option is given by user
      PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-terse',terse,ierr))
      if (terse) then
        PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr))
      else
        PetscCallA(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr))
        PetscCallA(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr))
      endif
      PetscCallA(EPSDestroy(eps,ierr))
      PetscCallA(MatDestroy(A,ierr))
      PetscCallA(SlepcFinalize(ierr))
      end

! -------------------------------------------------------------------
!
!   MatMult_Ising - user provided matrix-vector multiply
!
!   Input Parameters:
!   A - matrix
!   x - input vector
!
!   Output Parameter:
!   y - output vector
!
      subroutine MatMult_Ising(A,x,y,ierr)
#include <petsc/finclude/petscmat.h>
      use petscmat
      implicit none

      Mat            A
      Vec            x,y
      PetscInt       trans,one,N
      PetscScalar, pointer :: xx(:),yy(:)
      PetscErrorCode ierr

!     The actual routine for the matrix-vector product
      external mvmisg

      PetscCall(MatGetSize(A,N,PETSC_NULL_INTEGER,ierr))
      PetscCall(VecGetArrayRead(x,xx,ierr))
      PetscCall(VecGetArray(y,yy,ierr))

      trans = 0
      one = 1
      call mvmisg(trans,N,one,xx,N,yy,N)

      PetscCall(VecRestoreArrayRead(x,xx,ierr))
      PetscCall(VecRestoreArray(y,yy,ierr))

      end

! -------------------------------------------------------------------
!
!   MatMultTranspose_Ising - user provided transpose matrix-vector multiply
!
!   Input Parameters:
!   A - matrix
!   x - input vector
!
!   Output Parameter:
!   y - output vector
!
      subroutine MatMultTranspose_Ising(A,x,y,ierr)
#include <petsc/finclude/petscmat.h>
      use petscmat
      implicit none

      Mat            A
      Vec            x,y
      PetscInt       trans,one,N
      PetscScalar, pointer :: xx(:),yy(:)
      PetscErrorCode ierr

!     The actual routine for the transpose matrix-vector product
      external mvmisg

      PetscCall(MatGetSize(A,N,PETSC_NULL_INTEGER,ierr))
      PetscCall(VecGetArrayRead(x,xx,ierr))
      PetscCall(VecGetArray(y,yy,ierr))

      trans = 1
      one = 1
      call mvmisg(trans,N,one,xx,N,yy,N)

      PetscCall(VecRestoreArrayRead(x,xx,ierr))
      PetscCall(VecRestoreArray(y,yy,ierr))

      end

! -------------------------------------------------------------------
!     The actual routine for the matrix-vector product
!     See https://math.nist.gov/MatrixMarket/data/NEP/mvmisg/mvmisg.html

      SUBROUTINE MVMISG(TRANS, N, M, X, LDX, Y, LDY)
#include <petsc/finclude/petscsys.h>
      use petscsys
!     ..
!     .. Scalar Arguments ..
      PetscInt     LDY, LDX, M, N, TRANS
!     ..
!     .. Array Arguments ..
      PetscScalar  Y(LDY, *), X(LDX, *)
!     ..
!
!  Purpose
!  =======
!
!  Compute
!
!               Y(:,1:M) = op(A)*X(:,1:M)
!
!  where op(A) is A or A' (the transpose of A). The A is the Ising
!  matrix.
!
!  Arguments
!  =========
!
!  TRANS   (input) INTEGER
!          If TRANS = 0, compute Y(:,1:M) = A*X(:,1:M)
!          If TRANS = 1, compute Y(:,1:M) = A'*X(:,1:M)
!
!  N       (input) INTEGER
!          The order of the matrix A. N has to be an even number.
!
!  M       (input) INTEGER
!          The number of columns of X to multiply.
!
!  X       (input) DOUBLE PRECISION array, dimension (LDX, M)
!          X contains the matrix (vectors) X.
!
!  LDX     (input) INTEGER
!          The leading dimension of array X, LDX >= max(1, N)
!
!  Y       (output) DOUBLE PRECISION array, dimension (LDX, M)
!          contains the product of the matrix op(A) with X.
!
!  LDY     (input) INTEGER
!          The leading dimension of array Y, LDY >= max(1, N)
!
!  ===================================================================
!
!

!     .. Local Variables ..
      PetscInt    I, K
      PetscReal   ALPHA, BETA
      PetscReal   COSA, COSB, SINA, SINB
      PetscScalar TEMP, TEMP1
!
!     .. Intrinsic functions ..
      INTRINSIC   COS, SIN
!
      ALPHA = PETSC_PI/4
      BETA = PETSC_PI/4
      COSA = COS(ALPHA)
      SINA = SIN(ALPHA)
      COSB = COS(BETA)
      SINB = SIN(BETA)
!
      IF (TRANS.EQ.0) THEN
!
!     Compute Y(:,1:M) = A*X(:,1:M)

         DO 30 K = 1, M
!
            Y(1, K) = COSB*X(1, K) - SINB*X(N, K)
            DO 10 I = 2, N-1, 2
               Y(I, K)   =  COSB*X(I, K) + SINB*X(I+1, K)
               Y(I+1, K) = -SINB*X(I, K) + COSB*X(I+1, K)
   10       CONTINUE
            Y(N, K) = SINB*X(1, K) + COSB*X(N, K)
!
            DO 20 I = 1, N, 2
               TEMP      =  COSA*Y(I, K) + SINA*Y(I+1, K)
               Y(I+1, K) = -SINA*Y(I, K) + COSA*Y(I+1, K)
               Y(I, K)   = TEMP
   20       CONTINUE
!
   30    CONTINUE
!
      ELSE IF (TRANS.EQ.1) THEN
!
!        Compute Y(:1:M) = A'*X(:,1:M)
!
         DO 60 K = 1, M
!
            DO 40 I = 1, N, 2
               Y(I, K)   =  COSA*X(I, K) - SINA*X(I+1, K)
               Y(I+1, K) =  SINA*X(I, K) + COSA*X(I+1, K)
   40       CONTINUE
            TEMP  = COSB*Y(1,K) + SINB*Y(N,K)
            DO 50 I = 2, N-1, 2
               TEMP1     =  COSB*Y(I, K) - SINB*Y(I+1, K)
               Y(I+1, K) =  SINB*Y(I, K) + COSB*Y(I+1, K)
               Y(I, K)   =  TEMP1
   50       CONTINUE
            Y(N, K) = -SINB*Y(1, K) + COSB*Y(N, K)
            Y(1, K) = TEMP
!
   60    CONTINUE
!
      END IF
!
      RETURN
!
!     END OF MVMISG
      END

!/*TEST
!
!   testset:
!      args: -eps_max_it 1000 -eps_ncv 12 -eps_tol 1e-5 -eps_nev 4 -eps_largest_imaginary -terse
!      requires: !complex
!      output_file: output/ex6f_1.out
!      filter: grep -v iterations | sed -e 's/-0.00000/0.00000/g'
!      test:
!         suffix: 1
!      test:
!         suffix: 1_ts
!         args: -eps_two_sided
!         requires: !single
!
!TEST*/
