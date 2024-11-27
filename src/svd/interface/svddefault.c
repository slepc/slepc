/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Simple default routines for common SVD operations
*/

#include <slepc/private/svdimpl.h>      /*I "slepcsvd.h" I*/

/*
  SVDConvergedAbsolute - Checks convergence absolutely.
*/
PetscErrorCode SVDConvergedAbsolute(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = res;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  SVDConvergedRelative - Checks convergence relative to the singular value.
*/
PetscErrorCode SVDConvergedRelative(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = (sigma!=0.0)? res/sigma: PETSC_MAX_REAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  SVDConvergedNorm - Checks convergence relative to the matrix norms.
*/
PetscErrorCode SVDConvergedNorm(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = res/PetscMax(svd->nrma,svd->nrmb);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  SVDConvergedMaxIt - Always returns Inf to force reaching the maximum number of iterations.
*/
PetscErrorCode SVDConvergedMaxIt(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = PETSC_MAX_REAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDStoppingBasic - Default routine to determine whether the outer singular value
   solver iteration must be stopped.

   Collective

   Input Parameters:
+  svd    - singular value solver context obtained from SVDCreate()
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged singular triplets
.  nsv    - number of requested singular triplets
-  ctx    - context (not used here)

   Output Parameter:
.  reason - result of the stopping test

   Notes:
   A positive value of reason indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, reason must be set to SVD_CONVERGED_ITERATING
   (zero).

   SVDStoppingBasic() will stop if all requested singular values are converged, or if
   the maximum number of iterations has been reached.

   Use SVDSetStoppingTest() to provide your own test instead of using this one.

   Level: advanced

.seealso: SVDSetStoppingTest(), SVDStoppingThreshold(), SVDConvergedReason, SVDGetConvergedReason()
@*/
PetscErrorCode SVDStoppingBasic(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nsv,SVDConvergedReason *reason,void *ctx)
{
  PetscFunctionBegin;
  *reason = SVD_CONVERGED_ITERATING;
  if (nconv >= nsv) {
    PetscCall(PetscInfo(svd,"Singular value solver finished successfully: %" PetscInt_FMT " singular triplets converged at iteration %" PetscInt_FMT "\n",nconv,its));
    *reason = SVD_CONVERGED_TOL;
  } else if (its >= max_it) {
    if (svd->conv == SVD_CONV_MAXIT) *reason = SVD_CONVERGED_MAXIT;
    else {
      *reason = SVD_DIVERGED_ITS;
      PetscCall(PetscInfo(svd,"Singular value solver iteration reached maximum number of iterations (%" PetscInt_FMT ")\n",its));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDStoppingThreshold - Routine to determine whether the outer singular value
   solver iteration must be stopped, according to some threshold for the computed values.

   Collective

   Input Parameters:
+  svd    - singular value solver context obtained from SVDCreate()
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged singular triplets (ignored here)
.  nsv    - number of requested singular triplets (ignored here)
-  ctx    - context containing additional data (SVDStoppingCtx)

   Output Parameter:
.  reason - result of the stopping test

   Notes:
   A positive value of reason indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, reason must be set to SVD_CONVERGED_ITERATING
   (zero).

   SVDStoppingThreshold() will stop when one of the computed singular values is not
   above/below the threshold given at SVDSetThreshold(). If a number of wanted singular
   values has been specified via SVDSetDimensions() then it is also taken into account,
   and the solver will stop when one of the two conditions (threshold or number of
   converged values) is met.

   Use SVDSetStoppingTest() to provide your own test instead of using this one.

   Level: advanced

.seealso: SVDSetStoppingTest(), SVDStoppingBasic(), SVDSetThreshold(), SVDSetDimensions(), SVDConvergedReason, SVDGetConvergedReason()
@*/
PetscErrorCode SVDStoppingThreshold(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nsv,SVDConvergedReason *reason,void *ctx)
{
  PetscReal thres,firstsv,lastsv;
  PetscBool rel;
  SVDWhich  which;

  PetscFunctionBegin;
  *reason = SVD_CONVERGED_ITERATING;
  firstsv = ((SVDStoppingCtx)ctx)->firstsv;
  lastsv  = ((SVDStoppingCtx)ctx)->lastsv;
  thres   = ((SVDStoppingCtx)ctx)->thres;
  rel     = ((SVDStoppingCtx)ctx)->threlative;
  which   = ((SVDStoppingCtx)ctx)->which;
  if (nconv && ((which==SVD_LARGEST && ((rel && lastsv<thres*firstsv) || (!rel && lastsv<thres))) || (which==SVD_SMALLEST && lastsv>thres))) {
    if (which==SVD_SMALLEST) PetscCall(PetscInfo(svd,"Singular value solver finished successfully: singular value %g is above the threshold %g\n",(double)lastsv,(double)thres));
    else if (!rel) PetscCall(PetscInfo(svd,"Singular value solver finished successfully: singular value %g is below the threshold %g\n",(double)lastsv,(double)thres));
    else PetscCall(PetscInfo(svd,"Singular value solver finished successfully: the ratio %g/%g is below the threshold %g\n",(double)lastsv,(double)firstsv,(double)thres));
    *reason = SVD_CONVERGED_TOL;
  } else if (nsv && nconv >= nsv) {
    PetscCall(PetscInfo(svd,"Singular value solver finished successfully: %" PetscInt_FMT " singular triplets converged at iteration %" PetscInt_FMT "\n",nconv,its));
    *reason = SVD_CONVERGED_TOL;
  } else if (its >= max_it) {
    if (svd->conv == SVD_CONV_MAXIT) *reason = SVD_CONVERGED_MAXIT;
    else {
      *reason = SVD_DIVERGED_ITS;
      PetscCall(PetscInfo(svd,"Singular value solver iteration reached maximum number of iterations (%" PetscInt_FMT ")\n",its));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetWorkVecs - Sets a number of work vectors into an SVD object.

   Collective

   Input Parameters:
+  svd    - singular value solver context
.  nleft  - number of work vectors of dimension equal to left singular vector
-  nright - number of work vectors of dimension equal to right singular vector

   Developer Notes:
   This is SLEPC_EXTERN because it may be required by user plugin SVD
   implementations.

   Level: developer

.seealso: SVDSetUp()
@*/
PetscErrorCode SVDSetWorkVecs(SVD svd,PetscInt nleft,PetscInt nright)
{
  Vec            t;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,nleft,2);
  PetscValidLogicalCollectiveInt(svd,nright,3);
  PetscCheck(nleft>=0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"nleft must be >= 0: nleft = %" PetscInt_FMT,nleft);
  PetscCheck(nright>=0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"nright must be >= 0: nright = %" PetscInt_FMT,nright);
  PetscCheck(nleft>0 || nright>0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"nleft and nright cannot be both zero");
  if (svd->nworkl < nleft) {
    PetscCall(VecDestroyVecs(svd->nworkl,&svd->workl));
    svd->nworkl = nleft;
    if (svd->isgeneralized) PetscCall(SVDCreateLeftTemplate(svd,&t));
    else PetscCall(MatCreateVecsEmpty(svd->OP,NULL,&t));
    PetscCall(VecDuplicateVecs(t,nleft,&svd->workl));
    PetscCall(VecDestroy(&t));
  }
  if (svd->nworkr < nright) {
    PetscCall(VecDestroyVecs(svd->nworkr,&svd->workr));
    svd->nworkr = nright;
    PetscCall(MatCreateVecsEmpty(svd->OP,&t,NULL));
    PetscCall(VecDuplicateVecs(t,nright,&svd->workr));
    PetscCall(VecDestroy(&t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
