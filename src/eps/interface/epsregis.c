/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/epsimpl.h>  /*I "slepceps.h" I*/

SLEPC_EXTERN PetscErrorCode EPSCreate_Power(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_Subspace(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_Arnoldi(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_Lanczos(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_KrylovSchur(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_GD(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_JD(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_RQCG(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_LOBPCG(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_CISS(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_LyapII(EPS);
SLEPC_EXTERN PetscErrorCode EPSCreate_LAPACK(EPS);
#if defined(SLEPC_HAVE_ARPACK)
SLEPC_EXTERN PetscErrorCode EPSCreate_ARPACK(EPS);
#endif
#if defined(SLEPC_HAVE_TRLAN)
SLEPC_EXTERN PetscErrorCode EPSCreate_TRLAN(EPS);
#endif
#if defined(SLEPC_HAVE_BLOPEX)
SLEPC_EXTERN PetscErrorCode EPSCreate_BLOPEX(EPS);
#endif
#if defined(SLEPC_HAVE_PRIMME)
SLEPC_EXTERN PetscErrorCode EPSCreate_PRIMME(EPS);
#endif
#if defined(SLEPC_HAVE_FEAST)
SLEPC_EXTERN PetscErrorCode EPSCreate_FEAST(EPS);
#endif
#if defined(SLEPC_HAVE_SCALAPACK)
SLEPC_EXTERN PetscErrorCode EPSCreate_ScaLAPACK(EPS);
#endif
#if defined(SLEPC_HAVE_ELPA)
SLEPC_EXTERN PetscErrorCode EPSCreate_ELPA(EPS);
#endif
#if defined(SLEPC_HAVE_ELEMENTAL)
SLEPC_EXTERN PetscErrorCode EPSCreate_Elemental(EPS);
#endif
#if defined(SLEPC_HAVE_EVSL)
SLEPC_EXTERN PetscErrorCode EPSCreate_EVSL(EPS);
#endif
#if defined(SLEPC_HAVE_CHASE)
SLEPC_EXTERN PetscErrorCode EPSCreate_ChASE(EPS);
#endif

/*@C
  EPSRegisterAll - Registers all the eigenvalue solvers in the EPS package.

  Not Collective

  Level: advanced

.seealso: EPSRegister()
@*/
PetscErrorCode EPSRegisterAll(void)
{
  PetscFunctionBegin;
  if (EPSRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  EPSRegisterAllCalled = PETSC_TRUE;
  PetscCall(EPSRegister(EPSPOWER,EPSCreate_Power));
  PetscCall(EPSRegister(EPSSUBSPACE,EPSCreate_Subspace));
  PetscCall(EPSRegister(EPSARNOLDI,EPSCreate_Arnoldi));
  PetscCall(EPSRegister(EPSLANCZOS,EPSCreate_Lanczos));
  PetscCall(EPSRegister(EPSKRYLOVSCHUR,EPSCreate_KrylovSchur));
  PetscCall(EPSRegister(EPSGD,EPSCreate_GD));
  PetscCall(EPSRegister(EPSJD,EPSCreate_JD));
  PetscCall(EPSRegister(EPSRQCG,EPSCreate_RQCG));
  PetscCall(EPSRegister(EPSLOBPCG,EPSCreate_LOBPCG));
  PetscCall(EPSRegister(EPSCISS,EPSCreate_CISS));
  PetscCall(EPSRegister(EPSLYAPII,EPSCreate_LyapII));
  PetscCall(EPSRegister(EPSLAPACK,EPSCreate_LAPACK));
#if defined(SLEPC_HAVE_ARPACK)
  PetscCall(EPSRegister(EPSARPACK,EPSCreate_ARPACK));
#endif
#if defined(SLEPC_HAVE_TRLAN)
  PetscCall(EPSRegister(EPSTRLAN,EPSCreate_TRLAN));
#endif
#if defined(SLEPC_HAVE_BLOPEX)
  PetscCall(EPSRegister(EPSBLOPEX,EPSCreate_BLOPEX));
#endif
#if defined(SLEPC_HAVE_PRIMME)
  PetscCall(EPSRegister(EPSPRIMME,EPSCreate_PRIMME));
#endif
#if defined(SLEPC_HAVE_FEAST)
  PetscCall(EPSRegister(EPSFEAST,EPSCreate_FEAST));
#endif
#if defined(SLEPC_HAVE_SCALAPACK)
  PetscCall(EPSRegister(EPSSCALAPACK,EPSCreate_ScaLAPACK));
#endif
#if defined(SLEPC_HAVE_ELPA)
  PetscCall(EPSRegister(EPSELPA,EPSCreate_ELPA));
#endif
#if defined(SLEPC_HAVE_ELEMENTAL)
  PetscCall(EPSRegister(EPSELEMENTAL,EPSCreate_Elemental));
#endif
#if defined(SLEPC_HAVE_EVSL)
  PetscCall(EPSRegister(EPSEVSL,EPSCreate_EVSL));
#endif
#if defined(SLEPC_HAVE_CHASE)
  PetscCall(EPSRegister(EPSCHASE,EPSCreate_ChASE));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  EPSMonitorRegisterAll - Registers all the monitors in the EPS package.

  Not Collective

  Level: advanced

.seealso: EPSMonitorRegister()
@*/
PetscErrorCode EPSMonitorRegisterAll(void)
{
  PetscFunctionBegin;
  if (EPSMonitorRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  EPSMonitorRegisterAllCalled = PETSC_TRUE;

  PetscCall(EPSMonitorRegister("first_approximation",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,EPSMonitorFirst,NULL,NULL));
  PetscCall(EPSMonitorRegister("first_approximation",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,EPSMonitorFirstDrawLG,EPSMonitorFirstDrawLGCreate,NULL));
  PetscCall(EPSMonitorRegister("all_approximations",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,EPSMonitorAll,NULL,NULL));
  PetscCall(EPSMonitorRegister("all_approximations",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,EPSMonitorAllDrawLG,EPSMonitorAllDrawLGCreate,NULL));
  PetscCall(EPSMonitorRegister("convergence_history",PETSCVIEWERASCII,PETSC_VIEWER_DEFAULT,EPSMonitorConverged,EPSMonitorConvergedCreate,EPSMonitorConvergedDestroy));
  PetscCall(EPSMonitorRegister("convergence_history",PETSCVIEWERDRAW,PETSC_VIEWER_DRAW_LG,EPSMonitorConvergedDrawLG,EPSMonitorConvergedDrawLGCreate,EPSMonitorConvergedDestroy));
  PetscFunctionReturn(PETSC_SUCCESS);
}
