/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Simple default routines for common NEP operations
*/

#include <slepc/private/nepimpl.h>     /*I "slepcnep.h" I*/

/*@
   NEPSetWorkVecs - Sets a number of work vectors into a `NEP` object

   Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  nw  - number of work vectors to allocate

   Developer Note:
   This is `SLEPC_EXTERN` because it may be required by user plugin `NEP`
   implementations.

   Level: developer

.seealso: [](ch:nep), `NEPSetUp()`
@*/
PetscErrorCode NEPSetWorkVecs(NEP nep,PetscInt nw)
{
  Vec            t;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,nw,2);
  PetscCheck(nw>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"nw must be > 0: nw = %" PetscInt_FMT,nw);
  if (nep->nwork < nw) {
    PetscCall(VecDestroyVecs(nep->nwork,&nep->work));
    nep->nwork = nw;
    PetscCall(BVGetColumn(nep->V,0,&t));
    PetscCall(VecDuplicateVecs(t,nw,&nep->work));
    PetscCall(BVRestoreColumn(nep->V,0,&t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  NEPGetDefaultShift - Return the value of sigma to start the nonlinear iteration.
 */
PetscErrorCode NEPGetDefaultShift(NEP nep,PetscScalar *sigma)
{
  PetscFunctionBegin;
  PetscAssertPointer(sigma,2);
  switch (nep->which) {
    case NEP_LARGEST_MAGNITUDE:
    case NEP_LARGEST_IMAGINARY:
    case NEP_ALL:
    case NEP_WHICH_USER:
      *sigma = 1.0;   /* arbitrary value */
      break;
    case NEP_SMALLEST_MAGNITUDE:
    case NEP_SMALLEST_IMAGINARY:
      *sigma = 0.0;
      break;
    case NEP_LARGEST_REAL:
      *sigma = PETSC_MAX_REAL;
      break;
    case NEP_SMALLEST_REAL:
      *sigma = PETSC_MIN_REAL;
      break;
    case NEP_TARGET_MAGNITUDE:
    case NEP_TARGET_REAL:
    case NEP_TARGET_IMAGINARY:
      *sigma = nep->target;
      break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  NEPConvergedRelative - Checks convergence relative to the eigenvalue.
*/
PetscErrorCode NEPConvergedRelative(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal w;

  PetscFunctionBegin;
  w = SlepcAbsEigenvalue(eigr,eigi);
  *errest = (w!=0.0)? res/w: PETSC_MAX_REAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  NEPConvergedAbsolute - Checks convergence absolutely.
*/
PetscErrorCode NEPConvergedAbsolute(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = res;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  NEPConvergedNorm - Checks convergence relative to the matrix norms.
*/
PetscErrorCode NEPConvergedNorm(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscScalar    s;
  PetscReal      w=0.0;
  PetscInt       j;
  PetscBool      flg;

  PetscFunctionBegin;
  if (nep->fui!=NEP_USER_INTERFACE_SPLIT) {
    PetscCall(NEPComputeFunction(nep,eigr,nep->function,nep->function));
    PetscCall(MatHasOperation(nep->function,MATOP_NORM,&flg));
    PetscCheck(flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
    PetscCall(MatNorm(nep->function,NORM_INFINITY,&w));
  } else {
    /* initialization of matrix norms */
    if (!nep->nrma[0]) {
      for (j=0;j<nep->nt;j++) {
        PetscCall(MatHasOperation(nep->A[j],MATOP_NORM,&flg));
        PetscCheck(flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"The convergence test related to the matrix norms requires a matrix norm operation");
        PetscCall(MatNorm(nep->A[j],NORM_INFINITY,&nep->nrma[j]));
      }
    }
    for (j=0;j<nep->nt;j++) {
      PetscCall(FNEvaluateFunction(nep->f[j],eigr,&s));
      w = w + nep->nrma[j]*PetscAbsScalar(s);
    }
  }
  *errest = res/w;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPStoppingBasic - Default routine to determine whether the outer eigensolver
   iteration must be stopped.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged eigenpairs
.  nev    - number of requested eigenpairs
-  ctx    - context (not used here)

   Output Parameter:
.  reason - result of the stopping test

   Notes:
   `NEPStoppingBasic()` will stop if all requested eigenvalues are converged, or if
   the maximum number of iterations has been reached.

   This is the default stopping test.
   Use `NEPSetStoppingTest()` to provide your own test instead of using this one.

   Level: advanced

.seealso: [](ch:nep), `NEPSetStoppingTest()`, `NEPConvergedReason`, `NEPGetConvergedReason()`
@*/
PetscErrorCode NEPStoppingBasic(NEP nep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,NEPConvergedReason *reason,void *ctx)
{
  PetscFunctionBegin;
  *reason = NEP_CONVERGED_ITERATING;
  if (nconv >= nev) {
    PetscCall(PetscInfo(nep,"Nonlinear eigensolver finished successfully: %" PetscInt_FMT " eigenpairs converged at iteration %" PetscInt_FMT "\n",nconv,its));
    *reason = NEP_CONVERGED_TOL;
  } else if (its >= max_it) {
    *reason = NEP_DIVERGED_ITS;
    PetscCall(PetscInfo(nep,"Nonlinear eigensolver iteration reached maximum number of iterations (%" PetscInt_FMT ")\n",its));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPComputeVectors_Schur(NEP nep)
{
  Mat            Z;

  PetscFunctionBegin;
  PetscCall(DSVectors(nep->ds,DS_MAT_X,NULL,NULL));
  PetscCall(DSGetMat(nep->ds,DS_MAT_X,&Z));
  PetscCall(BVMultInPlace(nep->V,Z,0,nep->nconv));
  PetscCall(DSRestoreMat(nep->ds,DS_MAT_X,&Z));
  PetscCall(BVNormalize(nep->V,nep->eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}
