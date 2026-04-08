/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepcrg.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define rgpolygongetvertices0_ RGPOLYGONGETVERTICES0
#define rgpolygongetvertices1_ RGPOLYGONGETVERTICES1
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define rgpolygongetvertices0_ rgpolygongetvertices0
#define rgpolygongetvertices1_ rgpolygongetvertices1
#endif

SLEPC_EXTERN void rgpolygongetvertices_(RG *rg,PetscInt *n,PetscScalar *vr,PetscScalar *vi,PetscErrorCode *ierr)
{
  PetscScalar *ovr,*ovi;
  PetscInt    n_;

  CHKFORTRANNULLINTEGER(n);
  CHKFORTRANNULLSCALAR(vr);
  CHKFORTRANNULLSCALAR(vi);
  *ierr = RGPolygonGetVertices(*rg,&n_,&ovr,&ovi); if (*ierr) return;
  if (vr && ovr) { *ierr = PetscArraycpy(vr,ovr,n_); if (*ierr) return; }
  if (vi && ovi) { *ierr = PetscArraycpy(vi,ovi,n_); if (*ierr) return; }
  if (n) *n = n_;
  *ierr = PetscFree(ovr); if (*ierr) return;
  *ierr = PetscFree(ovi);
}

SLEPC_EXTERN void rgpolygongetvertices0_(RG *rg,PetscInt *n,PetscScalar *vr,PetscScalar *vi,PetscErrorCode *ierr)
{
  rgpolygongetvertices_(rg,n,vr,vi,ierr);
}

SLEPC_EXTERN void rgpolygongetvertices1_(RG *rg,PetscInt *n,PetscScalar *vr,PetscScalar *vi,PetscErrorCode *ierr)
{
  rgpolygongetvertices_(rg,n,vr,vi,ierr);
}
