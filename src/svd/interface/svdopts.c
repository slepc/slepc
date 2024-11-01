/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SVD routines for setting solver options
*/

#include <slepc/private/svdimpl.h>      /*I "slepcsvd.h" I*/
#include <petscdraw.h>

/*@
   SVDSetImplicitTranspose - Indicates how to handle the transpose of the matrix
   associated with the singular value problem.

   Logically Collective

   Input Parameters:
+  svd  - the singular value solver context
-  impl - how to handle the transpose (implicitly or not)

   Options Database Key:
.  -svd_implicittranspose - Activate the implicit transpose mode.

   Notes:
   By default, the transpose of the matrix is explicitly built (if the matrix
   has defined the MatTranspose operation).

   If this flag is set to true, the solver does not build the transpose, but
   handles it implicitly via MatMultTranspose() (or MatMultHermitianTranspose()
   in the complex case) operations. This is likely to be more inefficient
   than the default behaviour, both in sequential and in parallel, but
   requires less storage.

   Level: advanced

.seealso: SVDGetImplicitTranspose(), SVDSolve(), SVDSetOperators()
@*/
PetscErrorCode SVDSetImplicitTranspose(SVD svd,PetscBool impl)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,impl,2);
  if (svd->impltrans!=impl) {
    svd->impltrans = impl;
    svd->state     = SVD_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetImplicitTranspose - Gets the mode used to handle the transpose
   of the matrix associated with the singular value problem.

   Not Collective

   Input Parameter:
.  svd  - the singular value solver context

   Output Parameter:
.  impl - how to handle the transpose (implicitly or not)

   Level: advanced

.seealso: SVDSetImplicitTranspose(), SVDSolve(), SVDSetOperators()
@*/
PetscErrorCode SVDGetImplicitTranspose(SVD svd,PetscBool *impl)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(impl,2);
  *impl = svd->impltrans;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetTolerances - Sets the tolerance and maximum
   iteration count used by the default SVD convergence testers.

   Logically Collective

   Input Parameters:
+  svd - the singular value solver context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -svd_tol <tol> - Sets the convergence tolerance
-  -svd_max_it <maxits> - Sets the maximum number of iterations allowed

   Note:
   Use PETSC_CURRENT to retain the current value of any of the parameters.
   Use PETSC_DETERMINE for either argument to assign a default value computed
   internally (may be different in each solver).
   For maxits use PETSC_UMLIMITED to indicate there is no upper bound on this value.

   Level: intermediate

.seealso: SVDGetTolerances()
@*/
PetscErrorCode SVDSetTolerances(SVD svd,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveReal(svd,tol,2);
  PetscValidLogicalCollectiveInt(svd,maxits,3);
  if (tol == (PetscReal)PETSC_DETERMINE) {
    svd->tol   = PETSC_DETERMINE;
    svd->state = SVD_STATE_INITIAL;
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    svd->tol = tol;
  }
  if (maxits == PETSC_DETERMINE) {
    svd->max_it = PETSC_DETERMINE;
    svd->state  = SVD_STATE_INITIAL;
  } else if (maxits == PETSC_UNLIMITED) {
    svd->max_it = PETSC_INT_MAX;
  } else if (maxits != PETSC_CURRENT) {
    PetscCheck(maxits>0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    svd->max_it = maxits;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetTolerances - Gets the tolerance and maximum
   iteration count used by the default SVD convergence tests.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: SVDSetTolerances()
@*/
PetscErrorCode SVDGetTolerances(SVD svd,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (tol)    *tol    = svd->tol;
  if (maxits) *maxits = svd->max_it;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetThreshold - Sets the threshold used in the threshold stopping test.

   Logically Collective

   Input Parameters:
+  svd   - the singular value solver context
.  thres - the threshold value
-  rel   - whether the threshold is relative or not

   Options Database Keys:
+  -svd_threshold_absolute <thres> - Sets an absolute threshold
-  -svd_threshold_relative <thres> - Sets a relative threshold

   Notes:
   This function internally calls SVDSetStoppingTest() to set a special stopping
   test based on the threshold, where singular values are computed in sequence
   until one of the computed singular values is below the threshold.

   If the solver is configured to compute smallest singular values, then the
   threshold must be interpreted in the opposite direction, i.e., the computation
   will stop when one of the computed singular values is above the threshold.

   In the case of largest singular values, the threshold can be made relative
   with respect to the largest singular value (i.e., the matrix norm).

   The test against the threshold is done for converged singular values, which
   implies that the final number of converged singular values will be at least
   one more than the actual number of values below/above the threshold.

   Since the number of computed singular values is not known a priori, the solver
   will need to reallocate the basis of vectors internally, to have enough room
   to accommodate all the singular vectors. Hence, this option must be used with
   caution to avoid out-of-memory problems. The recommendation is to set the value
   of ncv to be larger than the estimated number of singular values, to minimize
   the number of reallocations.

   This functionality is most useful when computing largest singular values. A
   typical use case is to compute a low rank approximation of a matrix. Suppose
   we know that singular values decay abruptly around a certain index k, which
   is unknown. Then using a small relative threshold such as 0.2 will guarantee that
   the computed singular vectors capture the numerical rank k. However, if the matrix
   does not have low rank, i.e., singular values decay progressively, then a
   value of 0.2 will imply a very high cost, both computationally and in memory.

   If a number of wanted singular values has been set with SVDSetDimensions()
   it is also taken into account and the solver will stop when one of the two
   conditions (threshold or number of converged values) is met.

   Use SVDSetStoppingTest() to return to the usual computation of a fixed number
   of singular values.

   Level: advanced

.seealso: SVDGetThreshold(), SVDSetStoppingTest(), SVDSetDimensions(), SVDSetWhichSingularTriplets()
@*/
PetscErrorCode SVDSetThreshold(SVD svd,PetscReal thres,PetscBool rel)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveReal(svd,thres,2);
  PetscValidLogicalCollectiveBool(svd,rel,3);
  PetscCheck(thres>0.0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of the threshold. Must be > 0");
  if (svd->thres != thres || svd->threlative != rel) {
    svd->thres = thres;
    svd->threlative = rel;
    svd->state = SVD_STATE_INITIAL;
    PetscCall(SVDSetStoppingTest(svd,SVD_STOP_THRESHOLD));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetThreshold - Gets the threshold used by the threshold stopping test.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameters:
+  thres - the threshold
-  rel   - whether the threshold is relative or not

   Level: advanced

.seealso: SVDSetThreshold()
@*/
PetscErrorCode SVDGetThreshold(SVD svd,PetscReal *thres,PetscBool *rel)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (thres) *thres = svd->thres;
  if (rel)   *rel   = svd->threlative;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetDimensions - Sets the number of singular values to compute
   and the dimension of the subspace.

   Logically Collective

   Input Parameters:
+  svd - the singular value solver context
.  nsv - number of singular values to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Keys:
+  -svd_nsv <nsv> - Sets the number of singular values
.  -svd_ncv <ncv> - Sets the dimension of the subspace
-  -svd_mpd <mpd> - Sets the maximum projected dimension

   Notes:
   Use PETSC_DETERMINE for ncv and mpd to assign a reasonably good value, which is
   dependent on the solution method and the number of singular values required. For
   any of the arguments, use PETSC_CURRENT to preserve the current value.

   The parameters ncv and mpd are intimately related, so that the user is advised
   to set one of them at most. Normal usage is that
   (a) in cases where nsv is small, the user sets ncv (a reasonable default is 2*nsv); and
   (b) in cases where nsv is large, the user sets mpd.

   The value of ncv should always be between nsv and (nsv+mpd), typically
   ncv=nsv+mpd. If nsv is not too large, mpd=nsv is a reasonable choice, otherwise
   a smaller value should be used.

   Level: intermediate

.seealso: SVDGetDimensions()
@*/
PetscErrorCode SVDSetDimensions(SVD svd,PetscInt nsv,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,nsv,2);
  PetscValidLogicalCollectiveInt(svd,ncv,3);
  PetscValidLogicalCollectiveInt(svd,mpd,4);
  if (nsv != PETSC_CURRENT) {
    PetscCheck(nsv>0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nsv. Must be > 0");
    svd->nsv = nsv;
  }
  if (ncv == PETSC_DETERMINE) {
    svd->ncv = PETSC_DETERMINE;
  } else if (ncv != PETSC_CURRENT) {
    PetscCheck(ncv>0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    svd->ncv = ncv;
  }
  if (mpd == PETSC_DETERMINE) {
    svd->mpd = PETSC_DETERMINE;
  } else if (mpd != PETSC_CURRENT) {
    PetscCheck(mpd>0,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    svd->mpd = mpd;
  }
  svd->state = SVD_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetDimensions - Gets the number of singular values to compute
   and the dimension of the subspace.

   Not Collective

   Input Parameter:
.  svd - the singular value context

   Output Parameters:
+  nsv - number of singular values to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: SVDSetDimensions()
@*/
PetscErrorCode SVDGetDimensions(SVD svd,PetscInt *nsv,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (nsv) *nsv = svd->nsv? svd->nsv: 1;
  if (ncv) *ncv = svd->ncv;
  if (mpd) *mpd = svd->mpd;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
    SVDSetWhichSingularTriplets - Specifies which singular triplets are
    to be sought.

    Logically Collective

    Input Parameter:
.   svd - singular value solver context obtained from SVDCreate()

    Output Parameter:
.   which - which singular triplets are to be sought

    Options Database Keys:
+   -svd_largest  - Sets largest singular values
-   -svd_smallest - Sets smallest singular values

    Notes:
    The parameter 'which' can have one of these values

+     SVD_LARGEST  - largest singular values
-     SVD_SMALLEST - smallest singular values

    Level: intermediate

.seealso: SVDGetWhichSingularTriplets(), SVDWhich
@*/
PetscErrorCode SVDSetWhichSingularTriplets(SVD svd,SVDWhich which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,which,2);
  switch (which) {
    case SVD_LARGEST:
    case SVD_SMALLEST:
      if (svd->which != which) {
        svd->state = SVD_STATE_INITIAL;
        svd->which = which;
      }
      break;
  default:
    SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'which' parameter");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
    SVDGetWhichSingularTriplets - Returns which singular triplets are
    to be sought.

    Not Collective

    Input Parameter:
.   svd - singular value solver context obtained from SVDCreate()

    Output Parameter:
.   which - which singular triplets are to be sought

    Notes:
    See SVDSetWhichSingularTriplets() for possible values of which

    Level: intermediate

.seealso: SVDSetWhichSingularTriplets(), SVDWhich
@*/
PetscErrorCode SVDGetWhichSingularTriplets(SVD svd,SVDWhich *which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(which,2);
  *which = svd->which;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDSetConvergenceTestFunction - Sets a function to compute the error estimate
   used in the convergence test.

   Logically Collective

   Input Parameters:
+  svd     - singular value solver context obtained from SVDCreate()
.  conv    - the convergence test function, see SVDConvergenceTestFn for the calling sequence
.  ctx     - context for private data for the convergence routine (may be NULL)
-  destroy - a routine for destroying the context (may be NULL), see PetscCtxDestroyFn for the calling sequence

   Note:
   If the error estimate returned by the convergence test function is less than
   the tolerance, then the singular value is accepted as converged.

   Level: advanced

.seealso: SVDSetConvergenceTest(), SVDSetTolerances()
@*/
PetscErrorCode SVDSetConvergenceTestFunction(SVD svd,SVDConvergenceTestFn *conv,void* ctx,PetscCtxDestroyFn *destroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (svd->convergeddestroy) PetscCall((*svd->convergeddestroy)(&svd->convergedctx));
  svd->convergeduser    = conv;
  svd->convergeddestroy = destroy;
  svd->convergedctx     = ctx;
  if (conv == SVDConvergedAbsolute) svd->conv = SVD_CONV_ABS;
  else if (conv == SVDConvergedRelative) svd->conv = SVD_CONV_REL;
  else if (conv == SVDConvergedNorm) svd->conv = SVD_CONV_NORM;
  else if (conv == SVDConvergedMaxIt) svd->conv = SVD_CONV_MAXIT;
  else {
    svd->conv      = SVD_CONV_USER;
    svd->converged = svd->convergeduser;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetConvergenceTest - Specifies how to compute the error estimate
   used in the convergence test.

   Logically Collective

   Input Parameters:
+  svd  - singular value solver context obtained from SVDCreate()
-  conv - the type of convergence test

   Options Database Keys:
+  -svd_conv_abs   - Sets the absolute convergence test
.  -svd_conv_rel   - Sets the convergence test relative to the singular value
.  -svd_conv_norm  - Sets the convergence test relative to the matrix norm
.  -svd_conv_maxit - Forces the maximum number of iterations as set by -svd_max_it
-  -svd_conv_user  - Selects the user-defined convergence test

   Notes:
   The parameter 'conv' can have one of these values
+     SVD_CONV_ABS   - absolute error ||r||
.     SVD_CONV_REL   - error relative to the singular value sigma, ||r||/sigma
.     SVD_CONV_NORM  - error relative to the matrix norms, ||r||/||Z||, with Z=A or Z=[A;B]
.     SVD_CONV_MAXIT - no convergence until maximum number of iterations has been reached
-     SVD_CONV_USER  - function set by SVDSetConvergenceTestFunction()

   The default in standard SVD is SVD_CONV_REL, while in GSVD the default is SVD_CONV_NORM.

   Level: intermediate

.seealso: SVDGetConvergenceTest(), SVDSetConvergenceTestFunction(), SVDSetStoppingTest(), SVDConv
@*/
PetscErrorCode SVDSetConvergenceTest(SVD svd,SVDConv conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,conv,2);
  switch (conv) {
    case SVD_CONV_ABS:   svd->converged = SVDConvergedAbsolute; break;
    case SVD_CONV_REL:   svd->converged = SVDConvergedRelative; break;
    case SVD_CONV_NORM:  svd->converged = SVDConvergedNorm; break;
    case SVD_CONV_MAXIT: svd->converged = SVDConvergedMaxIt; break;
    case SVD_CONV_USER:
      PetscCheck(svd->convergeduser,PetscObjectComm((PetscObject)svd),PETSC_ERR_ORDER,"Must call SVDSetConvergenceTestFunction() first");
      svd->converged = svd->convergeduser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'conv' value");
  }
  svd->conv = conv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetConvergenceTest - Gets the method used to compute the error estimate
   used in the convergence test.

   Not Collective

   Input Parameters:
.  svd   - singular value solver context obtained from SVDCreate()

   Output Parameters:
.  conv  - the type of convergence test

   Level: intermediate

.seealso: SVDSetConvergenceTest(), SVDConv
@*/
PetscErrorCode SVDGetConvergenceTest(SVD svd,SVDConv *conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(conv,2);
  *conv = svd->conv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDSetStoppingTestFunction - Sets a function to decide when to stop the outer
   iteration of the singular value solver.

   Logically Collective

   Input Parameters:
+  svd     - singular value solver context obtained from SVDCreate()
.  stop    - the stopping test function, see SVDStoppingTestFn for the calling sequence
.  ctx     - context for private data for the stopping routine (may be NULL)
-  destroy - a routine for destroying the context (may be NULL), see PetscCtxDestroyFn for the calling sequence

   Note:
   Normal usage is to first call the default routine SVDStoppingBasic() and then
   set reason to SVD_CONVERGED_USER if some user-defined conditions have been
   met. To let the singular value solver continue iterating, the result must be
   left as SVD_CONVERGED_ITERATING.

   Level: advanced

.seealso: SVDSetStoppingTest(), SVDStoppingBasic()
@*/
PetscErrorCode SVDSetStoppingTestFunction(SVD svd,SVDStoppingTestFn *stop,void* ctx,PetscCtxDestroyFn *destroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (svd->stoppingdestroy) PetscCall((*svd->stoppingdestroy)(&svd->stoppingctx));
  svd->stoppinguser    = stop;
  svd->stoppingdestroy = destroy;
  svd->stoppingctx     = ctx;
  if (stop == SVDStoppingBasic) PetscCall(SVDSetStoppingTest(svd,SVD_STOP_BASIC));
  else if (stop == SVDStoppingThreshold) PetscCall(SVDSetStoppingTest(svd,SVD_STOP_THRESHOLD));
  else {
    svd->stop     = SVD_STOP_USER;
    svd->stopping = svd->stoppinguser;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetStoppingTest - Specifies how to decide the termination of the outer
   loop of the singular value solver.

   Logically Collective

   Input Parameters:
+  svd  - singular value solver context obtained from SVDCreate()
-  stop - the type of stopping test

   Options Database Keys:
+  -svd_stop_basic     - Sets the default stopping test
.  -svd_stop_threshold - Sets the threshold stopping test
-  -svd_stop_user      - Selects the user-defined stopping test

   Note:
   The parameter 'stop' can have one of these values
+     SVD_STOP_BASIC     - default stopping test
.     SVD_STOP_THRESHOLD - threshold stopping test
-     SVD_STOP_USER      - function set by SVDSetStoppingTestFunction()

   Level: advanced

.seealso: SVDGetStoppingTest(), SVDSetStoppingTestFunction(), SVDSetConvergenceTest(), SVDStop
@*/
PetscErrorCode SVDSetStoppingTest(SVD svd,SVDStop stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,stop,2);
  switch (stop) {
    case SVD_STOP_BASIC: svd->stopping = SVDStoppingBasic; break;
    case SVD_STOP_THRESHOLD: svd->stopping = SVDStoppingThreshold; break;
    case SVD_STOP_USER:
      PetscCheck(svd->stoppinguser,PetscObjectComm((PetscObject)svd),PETSC_ERR_ORDER,"Must call SVDSetStoppingTestFunction() first");
      svd->stopping = svd->stoppinguser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'stop' value");
  }
  svd->stop = stop;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetStoppingTest - Gets the method used to decide the termination of the outer
   loop of the singular value solver.

   Not Collective

   Input Parameters:
.  svd   - singular value solver context obtained from SVDCreate()

   Output Parameters:
.  stop  - the type of stopping test

   Level: advanced

.seealso: SVDSetStoppingTest(), SVDStop
@*/
PetscErrorCode SVDGetStoppingTest(SVD svd,SVDStop *stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(stop,2);
  *stop = svd->stop;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective

   Input Parameters:
+  svd      - the singular value solver context
.  opt      - the command line option for this monitor
.  name     - the monitor type one is seeking
.  ctx      - an optional user context for the monitor, or NULL
-  trackall - whether this monitor tracks all singular values or not

   Level: developer

.seealso: SVDMonitorSet(), SVDSetTrackAll()
@*/
PetscErrorCode SVDMonitorSetFromOptions(SVD svd,const char opt[],const char name[],void *ctx,PetscBool trackall)
{
  PetscErrorCode       (*mfunc)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);
  PetscErrorCode       (*cfunc)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
  PetscErrorCode       (*dfunc)(PetscViewerAndFormat**);
  PetscViewerAndFormat *vf;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerType      vtype;
  char                 key[PETSC_MAX_PATH_LEN];
  PetscBool            flg;

  PetscFunctionBegin;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,opt,&viewer,&format,&flg));
  if (!flg) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetType(viewer,&vtype));
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListFind(SVDMonitorList,key,&mfunc));
  PetscCheck(mfunc,PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Specified viewer and format not supported");
  PetscCall(PetscFunctionListFind(SVDMonitorCreateList,key,&cfunc));
  PetscCall(PetscFunctionListFind(SVDMonitorDestroyList,key,&dfunc));
  if (!cfunc) cfunc = PetscViewerAndFormatCreate_Internal;
  if (!dfunc) dfunc = PetscViewerAndFormatDestroy;

  PetscCall((*cfunc)(viewer,format,ctx,&vf));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(SVDMonitorSet(svd,mfunc,vf,(PetscCtxDestroyFn*)dfunc));
  if (trackall) PetscCall(SVDSetTrackAll(svd,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetFromOptions - Sets SVD options from the options database.
   This routine must be called before SVDSetUp() if the user is to be
   allowed to set the solver type.

   Collective

   Input Parameters:
.  svd - the singular value solver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner

.seealso: SVDSetOptionsPrefix()
@*/
PetscErrorCode SVDSetFromOptions(SVD svd)
{
  char           type[256];
  PetscBool      set,flg,val,flg1,flg2,flg3;
  PetscInt       i,j,k;
  PetscReal      r;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscCall(SVDRegisterAll());
  PetscObjectOptionsBegin((PetscObject)svd);
    PetscCall(PetscOptionsFList("-svd_type","SVD solver method","SVDSetType",SVDList,(char*)(((PetscObject)svd)->type_name?((PetscObject)svd)->type_name:SVDCROSS),type,sizeof(type),&flg));
    if (flg) PetscCall(SVDSetType(svd,type));
    else if (!((PetscObject)svd)->type_name) PetscCall(SVDSetType(svd,SVDCROSS));

    PetscCall(PetscOptionsBoolGroupBegin("-svd_standard","Singular value decomposition (SVD)","SVDSetProblemType",&flg));
    if (flg) PetscCall(SVDSetProblemType(svd,SVD_STANDARD));
    PetscCall(PetscOptionsBoolGroup("-svd_generalized","Generalized singular value decomposition (GSVD)","SVDSetProblemType",&flg));
    if (flg) PetscCall(SVDSetProblemType(svd,SVD_GENERALIZED));
    PetscCall(PetscOptionsBoolGroupEnd("-svd_hyperbolic","Hyperbolic singular value decomposition (HSVD)","SVDSetProblemType",&flg));
    if (flg) PetscCall(SVDSetProblemType(svd,SVD_HYPERBOLIC));

    PetscCall(PetscOptionsBool("-svd_implicittranspose","Handle matrix transpose implicitly","SVDSetImplicitTranspose",svd->impltrans,&val,&flg));
    if (flg) PetscCall(SVDSetImplicitTranspose(svd,val));

    i = svd->max_it;
    PetscCall(PetscOptionsInt("-svd_max_it","Maximum number of iterations","SVDSetTolerances",svd->max_it,&i,&flg1));
    r = svd->tol;
    PetscCall(PetscOptionsReal("-svd_tol","Tolerance","SVDSetTolerances",SlepcDefaultTol(svd->tol),&r,&flg2));
    if (flg1 || flg2) PetscCall(SVDSetTolerances(svd,r,i));

    r = svd->thres;
    PetscCall(PetscOptionsReal("-svd_threshold_absolute","Absolute threshold","SVDSetThreshold",r,&r,&flg));
    if (flg) PetscCall(SVDSetThreshold(svd,r,PETSC_FALSE));
    PetscCall(PetscOptionsReal("-svd_threshold_relative","Relative threshold","SVDSetThreshold",r,&r,&flg));
    if (flg) PetscCall(SVDSetThreshold(svd,r,PETSC_TRUE));

    PetscCall(PetscOptionsBoolGroupBegin("-svd_conv_abs","Absolute error convergence test","SVDSetConvergenceTest",&flg));
    if (flg) PetscCall(SVDSetConvergenceTest(svd,SVD_CONV_ABS));
    PetscCall(PetscOptionsBoolGroup("-svd_conv_rel","Relative error convergence test","SVDSetConvergenceTest",&flg));
    if (flg) PetscCall(SVDSetConvergenceTest(svd,SVD_CONV_REL));
    PetscCall(PetscOptionsBoolGroup("-svd_conv_norm","Convergence test relative to the matrix norms","SVDSetConvergenceTest",&flg));
    if (flg) PetscCall(SVDSetConvergenceTest(svd,SVD_CONV_NORM));
    PetscCall(PetscOptionsBoolGroup("-svd_conv_maxit","Maximum iterations convergence test","SVDSetConvergenceTest",&flg));
    if (flg) PetscCall(SVDSetConvergenceTest(svd,SVD_CONV_MAXIT));
    PetscCall(PetscOptionsBoolGroupEnd("-svd_conv_user","User-defined convergence test","SVDSetConvergenceTest",&flg));
    if (flg) PetscCall(SVDSetConvergenceTest(svd,SVD_CONV_USER));

    PetscCall(PetscOptionsBoolGroupBegin("-svd_stop_basic","Stop iteration if all singular values converged or max_it reached","SVDSetStoppingTest",&flg));
    if (flg) PetscCall(SVDSetStoppingTest(svd,SVD_STOP_BASIC));
    PetscCall(PetscOptionsBoolGroup("-svd_stop_threshold","Stop iteration if a converged singular value is below/above the threshold","SVDSetStoppingTest",&flg));
    if (flg) PetscCall(SVDSetStoppingTest(svd,SVD_STOP_THRESHOLD));
    PetscCall(PetscOptionsBoolGroupEnd("-svd_stop_user","User-defined stopping test","SVDSetStoppingTest",&flg));
    if (flg) PetscCall(SVDSetStoppingTest(svd,SVD_STOP_USER));

    i = svd->nsv;
    PetscCall(PetscOptionsInt("-svd_nsv","Number of singular values to compute","SVDSetDimensions",svd->nsv,&i,&flg1));
    if (!flg1) i = PETSC_CURRENT;
    j = svd->ncv;
    PetscCall(PetscOptionsInt("-svd_ncv","Number of basis vectors","SVDSetDimensions",svd->ncv,&j,&flg2));
    k = svd->mpd;
    PetscCall(PetscOptionsInt("-svd_mpd","Maximum dimension of projected problem","SVDSetDimensions",svd->mpd,&k,&flg3));
    if (flg1 || flg2 || flg3) PetscCall(SVDSetDimensions(svd,i,j,k));

    PetscCall(PetscOptionsBoolGroupBegin("-svd_largest","Compute largest singular values","SVDSetWhichSingularTriplets",&flg));
    if (flg) PetscCall(SVDSetWhichSingularTriplets(svd,SVD_LARGEST));
    PetscCall(PetscOptionsBoolGroupEnd("-svd_smallest","Compute smallest singular values","SVDSetWhichSingularTriplets",&flg));
    if (flg) PetscCall(SVDSetWhichSingularTriplets(svd,SVD_SMALLEST));

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to SVDSetFromOptions()
    */
    PetscCall(PetscOptionsBool("-svd_monitor_cancel","Remove any hardwired monitor routines","SVDMonitorCancel",PETSC_FALSE,&flg,&set));
    if (set && flg) PetscCall(SVDMonitorCancel(svd));
    PetscCall(SVDMonitorSetFromOptions(svd,"-svd_monitor","first_approximation",NULL,PETSC_FALSE));
    PetscCall(SVDMonitorSetFromOptions(svd,"-svd_monitor_all","all_approximations",NULL,PETSC_TRUE));
    PetscCall(SVDMonitorSetFromOptions(svd,"-svd_monitor_conv","convergence_history",NULL,PETSC_FALSE));
    PetscCall(SVDMonitorSetFromOptions(svd,"-svd_monitor_conditioning","conditioning",NULL,PETSC_FALSE));

    /* -----------------------------------------------------------------------*/
    PetscCall(PetscOptionsName("-svd_view","Print detailed information on solver used","SVDView",&set));
    PetscCall(PetscOptionsName("-svd_view_vectors","View computed singular vectors","SVDVectorsView",&set));
    PetscCall(PetscOptionsName("-svd_view_values","View computed singular values","SVDValuesView",&set));
    PetscCall(PetscOptionsName("-svd_converged_reason","Print reason for convergence, and number of iterations","SVDConvergedReasonView",&set));
    PetscCall(PetscOptionsName("-svd_error_absolute","Print absolute errors of each singular triplet","SVDErrorView",&set));
    PetscCall(PetscOptionsName("-svd_error_relative","Print relative errors of each singular triplet","SVDErrorView",&set));
    PetscCall(PetscOptionsName("-svd_error_norm","Print errors relative to the matrix norms of each singular triplet","SVDErrorView",&set));

    PetscTryTypeMethod(svd,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)svd,PetscOptionsObject));
  PetscOptionsEnd();

  if (!svd->V) PetscCall(SVDGetBV(svd,&svd->V,NULL));
  PetscCall(BVSetFromOptions(svd->V));
  if (!svd->U) PetscCall(SVDGetBV(svd,NULL,&svd->U));
  PetscCall(BVSetFromOptions(svd->U));
  if (!svd->ds) PetscCall(SVDGetDS(svd,&svd->ds));
  PetscCall(SVDSetDSType(svd));
  PetscCall(DSSetFromOptions(svd->ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetProblemType - Specifies the type of the singular value problem.

   Logically Collective

   Input Parameters:
+  svd  - the singular value solver context
-  type - a known type of singular value problem

   Options Database Keys:
+  -svd_standard    - standard singular value decomposition (SVD)
.  -svd_generalized - generalized singular value problem (GSVD)
-  -svd_hyperbolic  - hyperbolic singular value problem (HSVD)

   Notes:
   The GSVD requires that two matrices have been passed via SVDSetOperators().
   The HSVD requires that a signature matrix has been passed via SVDSetSignature().

   Level: intermediate

.seealso: SVDSetOperators(), SVDSetSignature(), SVDSetType(), SVDGetProblemType(), SVDProblemType
@*/
PetscErrorCode SVDSetProblemType(SVD svd,SVDProblemType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,type,2);
  if (type == svd->problem_type) PetscFunctionReturn(PETSC_SUCCESS);
  switch (type) {
    case SVD_STANDARD:
      svd->isgeneralized = PETSC_FALSE;
      svd->ishyperbolic  = PETSC_FALSE;
      break;
    case SVD_GENERALIZED:
      svd->isgeneralized = PETSC_TRUE;
      svd->ishyperbolic  = PETSC_FALSE;
      break;
    case SVD_HYPERBOLIC:
      svd->isgeneralized = PETSC_FALSE;
      svd->ishyperbolic  = PETSC_TRUE;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_WRONG,"Unknown singular value problem type");
  }
  svd->problem_type = type;
  svd->state = SVD_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetProblemType - Gets the problem type from the SVD object.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  type - the problem type

   Level: intermediate

.seealso: SVDSetProblemType(), SVDProblemType
@*/
PetscErrorCode SVDGetProblemType(SVD svd,SVDProblemType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = svd->problem_type;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDIsGeneralized - Ask if the SVD object corresponds to a generalized
   singular value problem.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  is - the answer

   Level: intermediate

.seealso: SVDIsHyperbolic()
@*/
PetscErrorCode SVDIsGeneralized(SVD svd,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(is,2);
  *is = svd->isgeneralized;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDIsHyperbolic - Ask if the SVD object corresponds to a hyperbolic
   singular value problem.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  is - the answer

   Level: intermediate

.seealso: SVDIsGeneralized()
@*/
PetscErrorCode SVDIsHyperbolic(SVD svd,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(is,2);
  *is = svd->ishyperbolic;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetTrackAll - Specifies if the solver must compute the residual norm of all
   approximate singular value or not.

   Logically Collective

   Input Parameters:
+  svd      - the singular value solver context
-  trackall - whether to compute all residuals or not

   Notes:
   If the user sets trackall=PETSC_TRUE then the solver computes (or estimates)
   the residual norm for each singular value approximation. Computing the residual is
   usually an expensive operation and solvers commonly compute only the residual
   associated to the first unconverged singular value.

   The option '-svd_monitor_all' automatically activates this option.

   Level: developer

.seealso: SVDGetTrackAll()
@*/
PetscErrorCode SVDSetTrackAll(SVD svd,PetscBool trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,trackall,2);
  svd->trackall = trackall;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetTrackAll - Returns the flag indicating whether all residual norms must
   be computed or not.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  trackall - the returned flag

   Level: developer

.seealso: SVDSetTrackAll()
@*/
PetscErrorCode SVDGetTrackAll(SVD svd,PetscBool *trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(trackall,2);
  *trackall = svd->trackall;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetOptionsPrefix - Sets the prefix used for searching for all
   SVD options in the database.

   Logically Collective

   Input Parameters:
+  svd - the singular value solver context
-  prefix - the prefix string to prepend to all SVD option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different SVD contexts, one could call
.vb
      SVDSetOptionsPrefix(svd1,"svd1_")
      SVDSetOptionsPrefix(svd2,"svd2_")
.ve

   Level: advanced

.seealso: SVDAppendOptionsPrefix(), SVDGetOptionsPrefix()
@*/
PetscErrorCode SVDSetOptionsPrefix(SVD svd,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!svd->V) PetscCall(SVDGetBV(svd,&svd->V,&svd->U));
  PetscCall(BVSetOptionsPrefix(svd->V,prefix));
  PetscCall(BVSetOptionsPrefix(svd->U,prefix));
  if (!svd->ds) PetscCall(SVDGetDS(svd,&svd->ds));
  PetscCall(DSSetOptionsPrefix(svd->ds,prefix));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)svd,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDAppendOptionsPrefix - Appends to the prefix used for searching for all
   SVD options in the database.

   Logically Collective

   Input Parameters:
+  svd - the singular value solver context
-  prefix - the prefix string to prepend to all SVD option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: SVDSetOptionsPrefix(), SVDGetOptionsPrefix()
@*/
PetscErrorCode SVDAppendOptionsPrefix(SVD svd,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!svd->V) PetscCall(SVDGetBV(svd,&svd->V,&svd->U));
  PetscCall(BVAppendOptionsPrefix(svd->V,prefix));
  PetscCall(BVAppendOptionsPrefix(svd->U,prefix));
  if (!svd->ds) PetscCall(SVDGetDS(svd,&svd->ds));
  PetscCall(DSAppendOptionsPrefix(svd->ds,prefix));
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)svd,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetOptionsPrefix - Gets the prefix used for searching for all
   SVD options in the database.

   Not Collective

   Input Parameters:
.  svd - the singular value solver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: SVDSetOptionsPrefix(), SVDAppendOptionsPrefix()
@*/
PetscErrorCode SVDGetOptionsPrefix(SVD svd,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)svd,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}
