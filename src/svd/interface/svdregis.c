/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/svdimpl.h>       /*I "slepcsvd.h" I*/

SLEPC_EXTERN PetscErrorCode SVDCreate_Cross(SVD);
SLEPC_EXTERN PetscErrorCode SVDCreate_Cyclic(SVD);
SLEPC_EXTERN PetscErrorCode SVDCreate_LAPACK(SVD);
SLEPC_EXTERN PetscErrorCode SVDCreate_Lanczos(SVD);
SLEPC_EXTERN PetscErrorCode SVDCreate_TRLanczos(SVD);
SLEPC_EXTERN PetscErrorCode SVDCreate_Randomized(SVD);
#if defined(SLEPC_HAVE_SCALAPACK)
SLEPC_EXTERN PetscErrorCode SVDCreate_ScaLAPACK(SVD);
#endif
#if defined(SLEPC_HAVE_KSVD)
SLEPC_EXTERN PetscErrorCode SVDCreate_KSVD(SVD);
#endif
#if defined(SLEPC_HAVE_ELEMENTAL)
SLEPC_EXTERN PetscErrorCode SVDCreate_Elemental(SVD);
#endif
#if defined(SLEPC_HAVE_PRIMME)
SLEPC_EXTERN PetscErrorCode SVDCreate_PRIMME(SVD);
#endif

/*@C
   SVDRegisterAll - Registers all the singular value solvers in the `SVD` package.

   Not Collective

   Level: advanced

.seealso: [](ch:svd), `SVDRegister()`
@*/
PetscErrorCode SVDRegisterAll(void)
{
  PetscFunctionBegin;
  if (SVDRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  SVDRegisterAllCalled = PETSC_TRUE;
  PetscCall(SVDRegister(SVDCROSS,SVDCreate_Cross));
  PetscCall(SVDRegister(SVDCYCLIC,SVDCreate_Cyclic));
  PetscCall(SVDRegister(SVDLAPACK,SVDCreate_LAPACK));
  PetscCall(SVDRegister(SVDLANCZOS,SVDCreate_Lanczos));
  PetscCall(SVDRegister(SVDTRLANCZOS,SVDCreate_TRLanczos));
  PetscCall(SVDRegister(SVDRANDOMIZED,SVDCreate_Randomized));
#if defined(SLEPC_HAVE_SCALAPACK)
  PetscCall(SVDRegister(SVDSCALAPACK,SVDCreate_ScaLAPACK));
#endif
#if defined(SLEPC_HAVE_KSVD)
  PetscCall(SVDRegister(SVDKSVD,SVDCreate_KSVD));
#endif
#if defined(SLEPC_HAVE_ELEMENTAL)
  PetscCall(SVDRegister(SVDELEMENTAL,SVDCreate_Elemental));
#endif
#if defined(SLEPC_HAVE_PRIMME)
  PetscCall(SVDRegister(SVDPRIMME,SVDCreate_PRIMME));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  SVDMonitorRegisterAll - Registers all the monitors in the `SVD` package.

  Not Collective

  Level: advanced

.seealso: [](ch:svd), `SVDMonitorRegister()`
@*/
PetscErrorCode SVDMonitorRegisterAll(void)
{
  PetscFunctionBegin;
  if (SVDMonitorRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  SVDMonitorRegisterAllCalled = PETSC_TRUE;

  PetscCall(SVDMonitorRegister("first_approximation",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,SVDMonitorFirst,NULL,NULL));
  PetscCall(SVDMonitorRegister("first_approximation",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,SVDMonitorFirstDrawLG,SVDMonitorFirstDrawLGCreate,NULL));
  PetscCall(SVDMonitorRegister("all_approximations",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,SVDMonitorAll,NULL,NULL));
  PetscCall(SVDMonitorRegister("all_approximations",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,SVDMonitorAllDrawLG,SVDMonitorAllDrawLGCreate,NULL));
  PetscCall(SVDMonitorRegister("convergence_history",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,SVDMonitorConverged,SVDMonitorConvergedCreate,SVDMonitorConvergedDestroy));
  PetscCall(SVDMonitorRegister("convergence_history",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,SVDMonitorConvergedDrawLG,SVDMonitorConvergedDrawLGCreate,SVDMonitorConvergedDestroy));
  PetscCall(SVDMonitorRegister("conditioning",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,SVDMonitorConditioning,NULL,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}
