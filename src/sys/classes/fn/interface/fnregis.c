/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/

SLEPC_EXTERN PetscErrorCode FNCreate_Combine(FN);
SLEPC_EXTERN PetscErrorCode FNCreate_Rational(FN);
SLEPC_EXTERN PetscErrorCode FNCreate_Exp(FN);
SLEPC_EXTERN PetscErrorCode FNCreate_Log(FN);
SLEPC_EXTERN PetscErrorCode FNCreate_Phi(FN);
SLEPC_EXTERN PetscErrorCode FNCreate_Sqrt(FN);
SLEPC_EXTERN PetscErrorCode FNCreate_Invsqrt(FN);

/*@C
   FNRegisterAll - Registers all of the math functions in the `FN` package.

   Not Collective

   Level: advanced

.seealso: [](sec:fn), `FNRegister()`
@*/
PetscErrorCode FNRegisterAll(void)
{
  PetscFunctionBegin;
  if (FNRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  FNRegisterAllCalled = PETSC_TRUE;
  PetscCall(FNRegister(FNCOMBINE,FNCreate_Combine));
  PetscCall(FNRegister(FNRATIONAL,FNCreate_Rational));
  PetscCall(FNRegister(FNEXP,FNCreate_Exp));
  PetscCall(FNRegister(FNLOG,FNCreate_Log));
  PetscCall(FNRegister(FNPHI,FNCreate_Phi));
  PetscCall(FNRegister(FNSQRT,FNCreate_Sqrt));
  PetscCall(FNRegister(FNINVSQRT,FNCreate_Invsqrt));
  PetscFunctionReturn(PETSC_SUCCESS);
}
