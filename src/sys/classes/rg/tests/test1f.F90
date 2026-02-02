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
!  Description: Simple example that tests RG interface functions.
!
! ----------------------------------------------------------------------
!
#include <slepc/finclude/slepcrg.h>
program test1f
  use slepcrg
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  RG             :: rg
  PetscInt       :: i, n
  PetscInt       :: inside(1)
  PetscMPIInt    :: rank
  PetscErrorCode :: ierr
  PetscReal      :: re, im
  PetscScalar    :: ar, ai, cr(10), ci(10)
  PetscScalar    :: vr(7), vi(7)
  PetscScalar    :: center
  PetscReal      :: radius, vscale, a, b, c, d

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
  PetscCallA(RGCreate(PETSC_COMM_WORLD, rg, ierr))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Ellipse
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(RGSetType(rg, RGELLIPSE, ierr))
  center = 1.1
  radius = 2
  vscale = 0.1
  PetscCallA(RGEllipseSetParameters(rg, center, radius, vscale, ierr))
  PetscCallA(RGSetFromOptions(rg, ierr))
  PetscCallA(RGView(rg, PETSC_NULL_VIEWER, ierr))
  re = 0.1
  im = 0.3
#if defined(PETSC_USE_COMPLEX)
  ar = re + im*PETSC_i
  ai = 0.0
#else
  ar = re
  ai = im
#endif
  PetscCallA(RGCheckInside(rg, 1_PETSC_INT_KIND, [ar], [ai], inside, ierr))
  if (rank == 0) then
    if (inside(1) >= 0) then
      write (*, '(a,f4.1,a,f4.1,a)') 'Point (', re, ',', im, ') is inside the region'
    else
      write (*, '(a,f4.1,a,f4.1,a)') 'Point (', re, ',', im, ') is outside the region'
    end if
  end if

  PetscCallA(RGComputeBoundingBox(rg, a, b, c, d, ierr))
  if (rank == 0) then
    write (*, '(a,f4.1,a,f4.1,a,f4.1,a,f4.1,a)') 'Bounding box: [', a, ',', b, ']x[', c, ',', d, ']'
  end if

  if (rank == 0) then
    write (*, *) 'Contour points:'
  end if
  n = 10
  PetscCallA(RGComputeContour(rg, n, cr, ci, ierr))
  do i = 1, n
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(cr(i))
    im = PetscImaginaryPart(cr(i))
#else
    re = cr(i)
    im = ci(i)
#endif
    if (rank == 0) then
      write (*, '(a,f7.4,a,f7.4,a)') '(', re, ',', im, ')'
    end if
  end do

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Interval
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscCallA(RGSetType(rg, RGINTERVAL, ierr))
  a = -1
  b = 1
  c = -0.1
  d = 0.1
  PetscCallA(RGIntervalSetEndpoints(rg, a, b, c, d, ierr))
  PetscCallA(RGSetFromOptions(rg, ierr))
  PetscCallA(RGView(rg, PETSC_NULL_VIEWER, ierr))
  re = 0.2
  im = 0
#if defined(PETSC_USE_COMPLEX)
  ar = re + im*PETSC_i
  ai = 0.0
#else
  ar = re
  ai = im
#endif
  PetscCallA(RGCheckInside(rg, 1_PETSC_INT_KIND, [ar], [ai], inside, ierr))
  if (rank == 0) then
    if (inside(1) >= 0) then
      write (*, '(a,f4.1,a,f4.1,a)') 'Point (', re, ',', im, ') is inside the region'
    else
      write (*, '(a,f4.1,a,f4.1,a)') 'Point (', re, ',', im, ') is outside the region'
    end if
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Polygon
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if defined(PETSC_USE_COMPLEX)
  vr(1) = (0.0, 2.0)
  vr(2) = (1.0, 4.0)
  vr(3) = (2.0, 5.0)
  vr(4) = (4.0, 3.0)
  vr(5) = (5.0, 4.0)
  vr(6) = (6.0, 1.0)
  vr(7) = (2.0, 0.0)
#else
  vr(1) = 0.0
  vi(1) = 1.0
  vr(2) = 0.0
  vi(2) = -1.0
  vr(3) = 0.6
  vi(3) = -0.8
  vr(4) = 1.0
  vi(4) = -1.0
  vr(5) = 2.0
  vi(5) = 0.0
  vr(6) = 1.0
  vi(6) = 1.0
  vr(7) = 0.6
  vi(7) = 0.8
#endif
  PetscCallA(RGSetType(rg, RGPOLYGON, ierr))
  n = 7
  PetscCallA(RGPolygonSetVertices(rg, n, vr, vi, ierr))
  PetscCallA(RGSetFromOptions(rg, ierr))
  PetscCallA(RGView(rg, PETSC_NULL_VIEWER, ierr))
  re = 5
  im = 0.9
#if defined(PETSC_USE_COMPLEX)
  ar = re + im*PETSC_i
  ai = 0.0
#else
  ar = re
  ai = im
#endif
  PetscCallA(RGCheckInside(rg, 1_PETSC_INT_KIND, [ar], [ai], inside, ierr))
  if (rank == 0) then
    if (inside(1) >= 0) then
      write (*, '(a,f4.1,a,f4.1,a)') 'Point (', re, ',', im, ') is inside the region'
    else
      write (*, '(a,f4.1,a,f4.1,a)') 'Point (', re, ',', im, ') is outside the region'
    end if
  end if

! *** Clean up
  PetscCallA(RGDestroy(rg, ierr))
  PetscCallA(SlepcFinalize(ierr))
end program test1f

!/*TEST
!
!   test:
!      suffix: 1
!      requires: !complex
!
!   test:
!      suffix: 1_complex
!      requires: complex
!
!TEST*/
