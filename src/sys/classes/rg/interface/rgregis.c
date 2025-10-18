/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/

SLEPC_EXTERN PetscErrorCode RGCreate_Interval(RG);
SLEPC_EXTERN PetscErrorCode RGCreate_Ellipse(RG);
SLEPC_EXTERN PetscErrorCode RGCreate_Ring(RG);
SLEPC_EXTERN PetscErrorCode RGCreate_Polygon(RG);

/*@C
   RGRegisterAll - Registers all of the regions in the `RG` package.

   Not Collective

   Level: advanced

.seealso: [](sec:rg), `RGRegister()`
@*/
PetscErrorCode RGRegisterAll(void)
{
  PetscFunctionBegin;
  if (RGRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  RGRegisterAllCalled = PETSC_TRUE;
  PetscCall(RGRegister(RGINTERVAL,RGCreate_Interval));
  PetscCall(RGRegister(RGELLIPSE,RGCreate_Ellipse));
  PetscCall(RGRegister(RGRING,RGCreate_Ring));
  PetscCall(RGRegister(RGPOLYGON,RGCreate_Polygon));
  PetscFunctionReturn(PETSC_SUCCESS);
}
