!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex10f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Illustrates the use of shell spectral transformations.
!  The problem to be solved is the same as ex1.c and corresponds to the
!  Laplacian operator in 1 dimension
!
!  The command line options are:
!    nm <n>, where <n> is the number of grid subdivisions = matrix dimension
!
!  Note: this example illustrates old error checking with CHKERRA instead
!  of PetscCallA()
! ----------------------------------------------------------------------
!
!     Module contains data needed by shell ST
!
      module mymoduleex10f
#include <slepc/finclude/slepceps.h>
      use slepceps
      implicit none

      KSP myksp
      end module

! ----------------------------------------------------------------------

      program main
#include <slepc/finclude/slepceps.h>
      use slepceps
      use mymoduleex10f
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
      ST             st
      EPSType        tname
      PetscInt       n, i, Istart, Iend, one, two, three
      PetscInt       nev, row(1), col(3)
      PetscScalar    val(3)
      PetscBool      flg, isShell, terse
      PetscMPIInt    rank
      PetscErrorCode ierr

!     Note: Any user-defined Fortran routines MUST be declared as external.
      external STApply_User, STApplyTranspose_User, STApplyHermitianTranspose_User, STBackTransform_User

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      one = 1
      two = 2
      three = 3
      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRMPIA(ierr)
      n = 30
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRA(ierr)

      if (rank .eq. 0) then
        write(*,'(/A,I6/)') '1-D Laplacian Eigenproblem (shell-enabled), n=',n
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the operator matrix that defines the eigensystem, Ax=kx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call MatCreate(PETSC_COMM_WORLD,A,ierr);CHKERRA(ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
      call MatSetFromOptions(A,ierr);CHKERRA(ierr)

      call MatGetOwnershipRange(A,Istart,Iend,ierr);CHKERRA(ierr)
      if (Istart .eq. 0) then
        row(1) = 0
        col(1) = 0
        col(2) = 1
        val(1) =  2.0
        val(2) = -1.0
        call MatSetValues(A,one,row,two,col,val,INSERT_VALUES,ierr);CHKERRA(ierr)
        Istart = Istart+1
      endif
      if (Iend .eq. n) then
        row(1) = n-1
        col(1) = n-2
        col(2) = n-1
        val(1) = -1.0
        val(2) =  2.0
        call MatSetValues(A,one,row,two,col,val,INSERT_VALUES,ierr);CHKERRA(ierr)
        Iend = Iend-1
      endif
      val(1) = -1.0
      val(2) =  2.0
      val(3) = -1.0
      do i=Istart,Iend-1
        row(1) = i
        col(1) = i-1
        col(2) = i
        col(3) = i+1
        call MatSetValues(A,one,row,three,col,val,INSERT_VALUES,ierr);CHKERRA(ierr)
      enddo

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      call EPSCreate(PETSC_COMM_WORLD,eps,ierr);CHKERRA(ierr)

!     ** Set operators. In this case, it is a standard eigenvalue problem
      call EPSSetOperators(eps,A,PETSC_NULL_MAT,ierr);CHKERRA(ierr)
      call EPSSetProblemType(eps,EPS_NHEP,ierr);CHKERRA(ierr)

!     ** Set solver parameters at runtime
      call EPSSetFromOptions(eps,ierr);CHKERRA(ierr)

!     ** Initialize shell spectral transformation if selected by user
      call EPSGetST(eps,st,ierr);CHKERRA(ierr)
      call PetscObjectTypeCompare(st,STSHELL,isShell,ierr);CHKERRA(ierr)

      if (isShell) then
!       ** Change sorting criterion since this ST example computes values
!       ** closest to 0
        call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL,ierr);CHKERRA(ierr)

!       ** In Fortran, instead of a context for the user-defined spectral transform
!       ** we use a module containing any application-specific data, initialized here
        call KSPCreate(PETSC_COMM_WORLD,myksp,ierr);CHKERRA(ierr)
        call KSPAppendOptionsPrefix(myksp,"st_",ierr);CHKERRA(ierr)

!       ** (Required) Set the user-defined routine for applying the operator
        call STShellSetApply(st,STApply_User,ierr);CHKERRA(ierr)

!       ** (Optional) Set the user-defined routine for applying the transposed operator
        call STShellSetApplyTranspose(st,STApplyTranspose_User,ierr);CHKERRA(ierr)

#if defined(PETSC_USE_COMPLEX)
!       ** (Optional) Set the user-defined routine for applying the conjugate-transposed operator
        call STShellSetApplyHermitianTranspose(st,STApplyHermitianTranspose_User,ierr);CHKERRA(ierr)
#endif

!       ** (Optional) Set the user-defined routine for back-transformation
        call STShellSetBackTransform(st,STBackTransform_User,ierr);CHKERRA(ierr)

!       ** (Optional) Set a name for the transformation, used for STView()
        call PetscObjectSetName(st,'MyTransformation',ierr);CHKERRA(ierr)

!       ** (Optional) Do any setup required for the new transformation
        call KSPSetOperators(myksp,A,A,ierr);CHKERRA(ierr)
        call KSPSetFromOptions(myksp,ierr);CHKERRA(ierr)
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call EPSSolve(eps,ierr);CHKERRA(ierr)

!     ** Optional: Get some information from the solver and display it
      call EPSGetType(eps,tname,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,'(A,A,/)') ' Solution method: ', tname
      endif
      call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,'(A,I2)') ' Number of requested eigenvalues:',nev
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** show detailed info unless -terse option is given by user
      call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-terse',terse,ierr);CHKERRA(ierr)
      if (terse) then
        call EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr);CHKERRA(ierr)
      else
        call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr);CHKERRA(ierr)
        call EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
        call EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
        call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
      endif
      if (isShell) then
        call KSPDestroy(myksp,ierr);CHKERRA(ierr)
      endif
      call EPSDestroy(eps,ierr);CHKERRA(ierr)
      call MatDestroy(A,ierr);CHKERRA(ierr)
      call SlepcFinalize(ierr)
      end

! -------------------------------------------------------------------
!
!   STApply_User - This routine demonstrates the use of a user-provided spectral
!   transformation. The transformation implemented in this code is just OP=A^-1.
!
!   Input Parameters:
!   st - spectral transformation context
!   x - input vector
!
!   Output Parameter:
!   y - output vector
!
      subroutine STApply_User(st,x,y,ierr)
#include <slepc/finclude/slepceps.h>
      use slepceps
      use mymoduleex10f
      implicit none

      ST             st
      Vec            x,y
      PetscErrorCode ierr

      call KSPSolve(myksp,x,y,ierr);CHKERRQ(ierr)

      end

! -------------------------------------------------------------------
!
!   STApplyTranspose_User - This is not required unless using a two-sided eigensolver
!
!   Input Parameters:
!   st - spectral transformation context
!   x - input vector
!
!   Output Parameter:
!   y - output vector
!
      subroutine STApplyTranspose_User(st,x,y,ierr)
#include <slepc/finclude/slepceps.h>
      use slepceps
      use mymoduleex10f
      implicit none

      ST             st
      Vec            x,y
      PetscErrorCode ierr

      call KSPSolveTranspose(myksp,x,y,ierr);CHKERRQ(ierr)

      end

#if defined(PETSC_USE_COMPLEX)
! -------------------------------------------------------------------
!
!   STApplyHermitianTranspose_User - This is not required unless using a two-sided eigensolver
!   in complex scalars
!
!   Input Parameters:
!   st - spectral transformation context
!   x - input vector
!
!   Output Parameter:
!   y - output vector
!
      subroutine STApplyHermitianTranspose_User(st,x,y,ierr)
#include <slepc/finclude/slepceps.h>
      use slepceps
      use mymoduleex10f
      implicit none

      ST             st
      Vec            x,y,w
      PetscErrorCode ierr

      call VecDuplicate(x,w,ierr);CHKERRQ(ierr)
      call VecCopy(x,w,ierr);CHKERRQ(ierr)
      call VecConjugate(w,ierr);CHKERRQ(ierr)
      call KSPSolveTranspose(myksp,w,y,ierr);CHKERRQ(ierr)
      call VecConjugate(y,ierr);CHKERRQ(ierr)
      call VecDestroy(w,ierr);CHKERRQ(ierr)

      end
#endif

! -------------------------------------------------------------------
!
!   STBackTransform_User - This routine demonstrates the use of a user-provided spectral
!   transformation
!
!   Input Parameters:
!   st - spectral transformation context
!   n  - number of eigenvalues to transform
!
!   Output Parameters:
!   eigr - real part of eigenvalues
!   eigi - imaginary part of eigenvalues
!
      subroutine STBackTransform_User(st,n,eigr,eigi,ierr)
#include <slepc/finclude/slepceps.h>
      use slepceps
      use mymoduleex10f
      implicit none

      ST             st
      PetscInt       n, j
      PetscScalar    eigr(*), eigi(*)
      PetscErrorCode ierr

      do j=1,n
        eigr(j) = 1.0 / eigr(j)
      enddo
      ierr = 0

      end

!/*TEST
!
!   testset:
!      args: -eps_nev 5 -eps_non_hermitian -terse
!      output_file: output/ex10_1.out
!      requires: !single
!      test:
!         suffix: 1_sinvert
!         args: -st_type sinvert
!      test:
!         suffix: 1_shell
!         args: -st_type shell -eps_two_sided {{0 1}}
!
!TEST*/
