/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to options that can be set via the command-line
   or procedurally
*/

#include <slepc/private/nepimpl.h>       /*I "slepcnep.h" I*/
#include <petscdraw.h>

/*@C
   NEPMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective

   Input Parameters:
+  nep      - the nonlinear eigensolver context
.  opt      - the command line option for this monitor
.  name     - the monitor type one is seeking
.  ctx      - an optional user context for the monitor, or NULL
-  trackall - whether this monitor tracks all eigenvalues or not

   Level: developer

.seealso: NEPMonitorSet(), NEPSetTrackAll()
@*/
PetscErrorCode NEPMonitorSetFromOptions(NEP nep,const char opt[],const char name[],void *ctx,PetscBool trackall)
{
  PetscErrorCode       (*mfunc)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
  PetscErrorCode       (*cfunc)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
  PetscErrorCode       (*dfunc)(PetscViewerAndFormat**);
  PetscViewerAndFormat *vf;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerType      vtype;
  char                 key[PETSC_MAX_PATH_LEN];
  PetscBool            flg;

  PetscFunctionBegin;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,opt,&viewer,&format,&flg));
  if (!flg) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetType(viewer,&vtype));
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListFind(NEPMonitorList,key,&mfunc));
  PetscCheck(mfunc,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Specified viewer and format not supported");
  PetscCall(PetscFunctionListFind(NEPMonitorCreateList,key,&cfunc));
  PetscCall(PetscFunctionListFind(NEPMonitorDestroyList,key,&dfunc));
  if (!cfunc) cfunc = PetscViewerAndFormatCreate_Internal;
  if (!dfunc) dfunc = PetscViewerAndFormatDestroy;

  PetscCall((*cfunc)(viewer,format,ctx,&vf));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(NEPMonitorSet(nep,mfunc,vf,(PetscCtxDestroyFn*)dfunc));
  if (trackall) PetscCall(NEPSetTrackAll(nep,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetFromOptions - Sets NEP options from the options database.
   This routine must be called before NEPSetUp() if the user is to be
   allowed to set the solver type.

   Collective

   Input Parameters:
.  nep - the nonlinear eigensolver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner

.seealso: NEPSetOptionsPrefix()
@*/
PetscErrorCode NEPSetFromOptions(NEP nep)
{
  char            type[256];
  PetscBool       set,flg,flg1,flg2,flg3,flg4,flg5,bval;
  PetscReal       r;
  PetscScalar     s;
  PetscInt        i,j,k;
  NEPRefine       refine;
  NEPRefineScheme scheme;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscCall(NEPRegisterAll());
  PetscObjectOptionsBegin((PetscObject)nep);
    PetscCall(PetscOptionsFList("-nep_type","Nonlinear eigensolver method","NEPSetType",NEPList,(char*)(((PetscObject)nep)->type_name?((PetscObject)nep)->type_name:NEPRII),type,sizeof(type),&flg));
    if (flg) PetscCall(NEPSetType(nep,type));
    else if (!((PetscObject)nep)->type_name) PetscCall(NEPSetType(nep,NEPRII));

    PetscCall(PetscOptionsBoolGroupBegin("-nep_general","General nonlinear eigenvalue problem","NEPSetProblemType",&flg));
    if (flg) PetscCall(NEPSetProblemType(nep,NEP_GENERAL));
    PetscCall(PetscOptionsBoolGroupEnd("-nep_rational","Rational eigenvalue problem","NEPSetProblemType",&flg));
    if (flg) PetscCall(NEPSetProblemType(nep,NEP_RATIONAL));

    refine = nep->refine;
    PetscCall(PetscOptionsEnum("-nep_refine","Iterative refinement method","NEPSetRefine",NEPRefineTypes,(PetscEnum)refine,(PetscEnum*)&refine,&flg1));
    i = nep->npart;
    PetscCall(PetscOptionsInt("-nep_refine_partitions","Number of partitions of the communicator for iterative refinement","NEPSetRefine",nep->npart,&i,&flg2));
    r = nep->rtol;
    PetscCall(PetscOptionsReal("-nep_refine_tol","Tolerance for iterative refinement","NEPSetRefine",nep->rtol==(PetscReal)PETSC_DETERMINE?SLEPC_DEFAULT_TOL/1000:nep->rtol,&r,&flg3));
    j = nep->rits;
    PetscCall(PetscOptionsInt("-nep_refine_its","Maximum number of iterations for iterative refinement","NEPSetRefine",nep->rits,&j,&flg4));
    scheme = nep->scheme;
    PetscCall(PetscOptionsEnum("-nep_refine_scheme","Scheme used for linear systems within iterative refinement","NEPSetRefine",NEPRefineSchemes,(PetscEnum)scheme,(PetscEnum*)&scheme,&flg5));
    if (flg1 || flg2 || flg3 || flg4 || flg5) PetscCall(NEPSetRefine(nep,refine,i,r,j,scheme));

    i = nep->max_it;
    PetscCall(PetscOptionsInt("-nep_max_it","Maximum number of iterations","NEPSetTolerances",nep->max_it,&i,&flg1));
    r = nep->tol;
    PetscCall(PetscOptionsReal("-nep_tol","Tolerance","NEPSetTolerances",SlepcDefaultTol(nep->tol),&r,&flg2));
    if (flg1 || flg2) PetscCall(NEPSetTolerances(nep,r,i));

    PetscCall(PetscOptionsBoolGroupBegin("-nep_conv_rel","Relative error convergence test","NEPSetConvergenceTest",&flg));
    if (flg) PetscCall(NEPSetConvergenceTest(nep,NEP_CONV_REL));
    PetscCall(PetscOptionsBoolGroup("-nep_conv_norm","Convergence test relative to the matrix norms","NEPSetConvergenceTest",&flg));
    if (flg) PetscCall(NEPSetConvergenceTest(nep,NEP_CONV_NORM));
    PetscCall(PetscOptionsBoolGroup("-nep_conv_abs","Absolute error convergence test","NEPSetConvergenceTest",&flg));
    if (flg) PetscCall(NEPSetConvergenceTest(nep,NEP_CONV_ABS));
    PetscCall(PetscOptionsBoolGroupEnd("-nep_conv_user","User-defined convergence test","NEPSetConvergenceTest",&flg));
    if (flg) PetscCall(NEPSetConvergenceTest(nep,NEP_CONV_USER));

    PetscCall(PetscOptionsBoolGroupBegin("-nep_stop_basic","Stop iteration if all eigenvalues converged or max_it reached","NEPSetStoppingTest",&flg));
    if (flg) PetscCall(NEPSetStoppingTest(nep,NEP_STOP_BASIC));
    PetscCall(PetscOptionsBoolGroupEnd("-nep_stop_user","User-defined stopping test","NEPSetStoppingTest",&flg));
    if (flg) PetscCall(NEPSetStoppingTest(nep,NEP_STOP_USER));

    i = nep->nev;
    PetscCall(PetscOptionsInt("-nep_nev","Number of eigenvalues to compute","NEPSetDimensions",nep->nev,&i,&flg1));
    j = nep->ncv;
    PetscCall(PetscOptionsInt("-nep_ncv","Number of basis vectors","NEPSetDimensions",nep->ncv,&j,&flg2));
    k = nep->mpd;
    PetscCall(PetscOptionsInt("-nep_mpd","Maximum dimension of projected problem","NEPSetDimensions",nep->mpd,&k,&flg3));
    if (flg1 || flg2 || flg3) PetscCall(NEPSetDimensions(nep,i,j,k));

    PetscCall(PetscOptionsBoolGroupBegin("-nep_largest_magnitude","Compute largest eigenvalues in magnitude","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_LARGEST_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-nep_smallest_magnitude","Compute smallest eigenvalues in magnitude","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_SMALLEST_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-nep_largest_real","Compute eigenvalues with largest real parts","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_LARGEST_REAL));
    PetscCall(PetscOptionsBoolGroup("-nep_smallest_real","Compute eigenvalues with smallest real parts","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_SMALLEST_REAL));
    PetscCall(PetscOptionsBoolGroup("-nep_largest_imaginary","Compute eigenvalues with largest imaginary parts","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_LARGEST_IMAGINARY));
    PetscCall(PetscOptionsBoolGroup("-nep_smallest_imaginary","Compute eigenvalues with smallest imaginary parts","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_SMALLEST_IMAGINARY));
    PetscCall(PetscOptionsBoolGroup("-nep_target_magnitude","Compute eigenvalues closest to target","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-nep_target_real","Compute eigenvalues with real parts closest to target","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_TARGET_REAL));
    PetscCall(PetscOptionsBoolGroup("-nep_target_imaginary","Compute eigenvalues with imaginary parts closest to target","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_TARGET_IMAGINARY));
    PetscCall(PetscOptionsBoolGroupEnd("-nep_all","Compute all eigenvalues in a region","NEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(NEPSetWhichEigenpairs(nep,NEP_ALL));

    PetscCall(PetscOptionsScalar("-nep_target","Value of the target","NEPSetTarget",nep->target,&s,&flg));
    if (flg) {
      if (nep->which!=NEP_TARGET_REAL && nep->which!=NEP_TARGET_IMAGINARY) PetscCall(NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE));
      PetscCall(NEPSetTarget(nep,s));
    }

    PetscCall(PetscOptionsBool("-nep_two_sided","Use two-sided variant (to compute left eigenvectors)","NEPSetTwoSided",nep->twosided,&bval,&flg));
    if (flg) PetscCall(NEPSetTwoSided(nep,bval));

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to NEPSetFromOptions()
    */
    PetscCall(PetscOptionsBool("-nep_monitor_cancel","Remove any hardwired monitor routines","NEPMonitorCancel",PETSC_FALSE,&flg,&set));
    if (set && flg) PetscCall(NEPMonitorCancel(nep));
    PetscCall(NEPMonitorSetFromOptions(nep,"-nep_monitor","first_approximation",NULL,PETSC_FALSE));
    PetscCall(NEPMonitorSetFromOptions(nep,"-nep_monitor_all","all_approximations",NULL,PETSC_TRUE));
    PetscCall(NEPMonitorSetFromOptions(nep,"-nep_monitor_conv","convergence_history",NULL,PETSC_FALSE));

    /* -----------------------------------------------------------------------*/
    PetscCall(PetscOptionsName("-nep_view","Print detailed information on solver used","NEPView",&set));
    PetscCall(PetscOptionsName("-nep_view_vectors","View computed eigenvectors","NEPVectorsView",&set));
    PetscCall(PetscOptionsName("-nep_view_values","View computed eigenvalues","NEPValuesView",&set));
    PetscCall(PetscOptionsName("-nep_converged_reason","Print reason for convergence, and number of iterations","NEPConvergedReasonView",&set));
    PetscCall(PetscOptionsName("-nep_error_absolute","Print absolute errors of each eigenpair","NEPErrorView",&set));
    PetscCall(PetscOptionsName("-nep_error_relative","Print relative errors of each eigenpair","NEPErrorView",&set));

    PetscTryTypeMethod(nep,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)nep,PetscOptionsObject));
  PetscOptionsEnd();

  if (!nep->V) PetscCall(NEPGetBV(nep,&nep->V));
  PetscCall(BVSetFromOptions(nep->V));
  if (!nep->rg) PetscCall(NEPGetRG(nep,&nep->rg));
  PetscCall(RGSetFromOptions(nep->rg));
  if (nep->useds) {
    if (!nep->ds) PetscCall(NEPGetDS(nep,&nep->ds));
    PetscCall(NEPSetDSType(nep));
    PetscCall(DSSetFromOptions(nep->ds));
  }
  if (!nep->refineksp) PetscCall(NEPRefineGetKSP(nep,&nep->refineksp));
  PetscCall(KSPSetFromOptions(nep->refineksp));
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) for (i=0;i<nep->nt;i++) PetscCall(FNSetFromOptions(nep->f[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetTolerances - Gets the tolerance and maximum iteration count used
   by the NEP convergence tests.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: NEPSetTolerances()
@*/
PetscErrorCode NEPGetTolerances(NEP nep,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (tol)    *tol    = nep->tol;
  if (maxits) *maxits = nep->max_it;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetTolerances - Sets the tolerance and maximum iteration count used
   by the NEP convergence tests.

   Logically Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  tol    - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -nep_tol <tol> - Sets the convergence tolerance
-  -nep_max_it <maxits> - Sets the maximum number of iterations allowed

   Notes:
   Use PETSC_CURRENT to retain the current value of any of the parameters.
   Use PETSC_DETERMINE for either argument to assign a default value computed
   internally (may be different in each solver).
   For maxits use PETSC_UMLIMITED to indicate there is no upper bound on this value.

   Level: intermediate

.seealso: NEPGetTolerances()
@*/
PetscErrorCode NEPSetTolerances(NEP nep,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,tol,2);
  PetscValidLogicalCollectiveInt(nep,maxits,3);
  if (tol == (PetscReal)PETSC_DETERMINE) {
    nep->tol   = PETSC_DETERMINE;
    nep->state = NEP_STATE_INITIAL;
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    nep->tol = tol;
  }
  if (maxits == PETSC_DETERMINE) {
    nep->max_it = PETSC_DETERMINE;
    nep->state  = NEP_STATE_INITIAL;
  } else if (maxits == PETSC_UNLIMITED) {
    nep->max_it = PETSC_INT_MAX;
  } else if (maxits != PETSC_CURRENT) {
    PetscCheck(maxits>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    nep->max_it = maxits;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetDimensions - Gets the number of eigenvalues to compute
   and the dimension of the subspace.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: NEPSetDimensions()
@*/
PetscErrorCode NEPGetDimensions(NEP nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (nev) *nev = nep->nev;
  if (ncv) *ncv = nep->ncv;
  if (mpd) *mpd = nep->mpd;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetDimensions - Sets the number of eigenvalues to compute
   and the dimension of the subspace.

   Logically Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Keys:
+  -nep_nev <nev> - Sets the number of eigenvalues
.  -nep_ncv <ncv> - Sets the dimension of the subspace
-  -nep_mpd <mpd> - Sets the maximum projected dimension

   Notes:
   Use PETSC_DETERMINE for ncv and mpd to assign a reasonably good value, which is
   dependent on the solution method. For any of the arguments, use PETSC_CURRENT
   to preserve the current value.

   The parameters ncv and mpd are intimately related, so that the user is advised
   to set one of them at most. Normal usage is that
   (a) in cases where nev is small, the user sets ncv (a reasonable default is 2*nev); and
   (b) in cases where nev is large, the user sets mpd.

   The value of ncv should always be between nev and (nev+mpd), typically
   ncv=nev+mpd. If nev is not too large, mpd=nev is a reasonable choice, otherwise
   a smaller value should be used.

   Level: intermediate

.seealso: NEPGetDimensions()
@*/
PetscErrorCode NEPSetDimensions(NEP nep,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,nev,2);
  PetscValidLogicalCollectiveInt(nep,ncv,3);
  PetscValidLogicalCollectiveInt(nep,mpd,4);
  if (nev != PETSC_CURRENT) {
    PetscCheck(nev>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
    nep->nev = nev;
  }
  if (ncv == PETSC_DETERMINE) {
    nep->ncv = PETSC_DETERMINE;
  } else if (ncv != PETSC_CURRENT) {
    PetscCheck(ncv>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    nep->ncv = ncv;
  }
  if (mpd == PETSC_DETERMINE) {
    nep->mpd = PETSC_DETERMINE;
  } else if (mpd != PETSC_CURRENT) {
    PetscCheck(mpd>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    nep->mpd = mpd;
  }
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
    NEPSetWhichEigenpairs - Specifies which portion of the spectrum is
    to be sought.

    Logically Collective

    Input Parameters:
+   nep   - eigensolver context obtained from NEPCreate()
-   which - the portion of the spectrum to be sought

    Options Database Keys:
+   -nep_largest_magnitude - Sets largest eigenvalues in magnitude
.   -nep_smallest_magnitude - Sets smallest eigenvalues in magnitude
.   -nep_largest_real - Sets largest real parts
.   -nep_smallest_real - Sets smallest real parts
.   -nep_largest_imaginary - Sets largest imaginary parts
.   -nep_smallest_imaginary - Sets smallest imaginary parts
.   -nep_target_magnitude - Sets eigenvalues closest to target
.   -nep_target_real - Sets real parts closest to target
.   -nep_target_imaginary - Sets imaginary parts closest to target
-   -nep_all - Sets all eigenvalues in a region

    Notes:
    The parameter 'which' can have one of these values

+     NEP_LARGEST_MAGNITUDE - largest eigenvalues in magnitude (default)
.     NEP_SMALLEST_MAGNITUDE - smallest eigenvalues in magnitude
.     NEP_LARGEST_REAL - largest real parts
.     NEP_SMALLEST_REAL - smallest real parts
.     NEP_LARGEST_IMAGINARY - largest imaginary parts
.     NEP_SMALLEST_IMAGINARY - smallest imaginary parts
.     NEP_TARGET_MAGNITUDE - eigenvalues closest to the target (in magnitude)
.     NEP_TARGET_REAL - eigenvalues with real part closest to target
.     NEP_TARGET_IMAGINARY - eigenvalues with imaginary part closest to target
.     NEP_ALL - all eigenvalues contained in a given region
-     NEP_WHICH_USER - user defined ordering set with NEPSetEigenvalueComparison()

    Not all eigensolvers implemented in NEP account for all the possible values
    stated above. If SLEPc is compiled for real numbers NEP_LARGEST_IMAGINARY
    and NEP_SMALLEST_IMAGINARY use the absolute value of the imaginary part
    for eigenvalue selection.

    The target is a scalar value provided with NEPSetTarget().

    NEP_ALL is intended for use in the context of the CISS solver for
    computing all eigenvalues in a region.

    Level: intermediate

.seealso: NEPGetWhichEigenpairs(), NEPSetTarget(), NEPSetEigenvalueComparison(), NEPWhich
@*/
PetscErrorCode NEPSetWhichEigenpairs(NEP nep,NEPWhich which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(nep,which,2);
  switch (which) {
    case NEP_LARGEST_MAGNITUDE:
    case NEP_SMALLEST_MAGNITUDE:
    case NEP_LARGEST_REAL:
    case NEP_SMALLEST_REAL:
    case NEP_LARGEST_IMAGINARY:
    case NEP_SMALLEST_IMAGINARY:
    case NEP_TARGET_MAGNITUDE:
    case NEP_TARGET_REAL:
#if defined(PETSC_USE_COMPLEX)
    case NEP_TARGET_IMAGINARY:
#endif
    case NEP_ALL:
    case NEP_WHICH_USER:
      if (nep->which != which) {
        nep->state = NEP_STATE_INITIAL;
        nep->which = which;
      }
      break;
#if !defined(PETSC_USE_COMPLEX)
    case NEP_TARGET_IMAGINARY:
      SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NEP_TARGET_IMAGINARY can be used only with complex scalars");
#endif
    default:
      SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'which' value");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
    NEPGetWhichEigenpairs - Returns which portion of the spectrum is to be
    sought.

    Not Collective

    Input Parameter:
.   nep - eigensolver context obtained from NEPCreate()

    Output Parameter:
.   which - the portion of the spectrum to be sought

    Notes:
    See NEPSetWhichEigenpairs() for possible values of 'which'.

    Level: intermediate

.seealso: NEPSetWhichEigenpairs(), NEPWhich
@*/
PetscErrorCode NEPGetWhichEigenpairs(NEP nep,NEPWhich *which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(which,2);
  *which = nep->which;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPSetEigenvalueComparison - Specifies the eigenvalue comparison function
   when NEPSetWhichEigenpairs() is set to NEP_WHICH_USER.

   Logically Collective

   Input Parameters:
+  nep  - eigensolver context obtained from NEPCreate()
.  comp - a pointer to the comparison function
-  ctx  - a context pointer (the last parameter to the comparison function)

   Calling sequence of comp:
$  PetscErrorCode comp(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx)
+   ar     - real part of the 1st eigenvalue
.   ai     - imaginary part of the 1st eigenvalue
.   br     - real part of the 2nd eigenvalue
.   bi     - imaginary part of the 2nd eigenvalue
.   res    - result of comparison
-   ctx    - optional context, as set by NEPSetEigenvalueComparison()

   Note:
   The returning parameter 'res' can be
+  negative - if the 1st eigenvalue is preferred to the 2st one
.  zero     - if both eigenvalues are equally preferred
-  positive - if the 2st eigenvalue is preferred to the 1st one

   Level: advanced

.seealso: NEPSetWhichEigenpairs(), NEPWhich
@*/
PetscErrorCode NEPSetEigenvalueComparison(NEP nep,PetscErrorCode (*comp)(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx),void* ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  nep->sc->comparison    = comp;
  nep->sc->comparisonctx = ctx;
  nep->which             = NEP_WHICH_USER;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetProblemType - Specifies the type of the nonlinear eigenvalue problem.

   Logically Collective

   Input Parameters:
+  nep  - the nonlinear eigensolver context
-  type - a known type of nonlinear eigenvalue problem

   Options Database Keys:
+  -nep_general - general problem with no particular structure
-  -nep_rational - a rational eigenvalue problem defined in split form with all f_i rational

   Notes:
   Allowed values for the problem type are general (NEP_GENERAL), and rational
   (NEP_RATIONAL).

   This function is used to provide a hint to the NEP solver to exploit certain
   properties of the nonlinear eigenproblem. This hint may be used or not,
   depending on the solver. By default, no particular structure is assumed.

   Level: intermediate

.seealso: NEPSetType(), NEPGetProblemType(), NEPProblemType
@*/
PetscErrorCode NEPSetProblemType(NEP nep,NEPProblemType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(nep,type,2);
  PetscCheck(type==NEP_GENERAL || type==NEP_RATIONAL,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"Unknown eigenvalue problem type");
  if (type != nep->problem_type) {
    nep->problem_type = type;
    nep->state = NEP_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetProblemType - Gets the problem type from the NEP object.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  type - the problem type

   Level: intermediate

.seealso: NEPSetProblemType(), NEPProblemType
@*/
PetscErrorCode NEPGetProblemType(NEP nep,NEPProblemType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = nep->problem_type;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetTwoSided - Sets the solver to use a two-sided variant so that left
   eigenvectors are also computed.

   Logically Collective

   Input Parameters:
+  nep      - the eigensolver context
-  twosided - whether the two-sided variant is to be used or not

   Options Database Keys:
.  -nep_two_sided <boolean> - Sets/resets the twosided flag

   Notes:
   If the user sets twosided=PETSC_TRUE then the solver uses a variant of
   the algorithm that computes both right and left eigenvectors. This is
   usually much more costly. This option is not available in all solvers.

   When using two-sided solvers, the problem matrices must have both the
   MatMult and MatMultTranspose operations defined.

   Level: advanced

.seealso: NEPGetTwoSided(), NEPGetLeftEigenvector()
@*/
PetscErrorCode NEPSetTwoSided(NEP nep,PetscBool twosided)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,twosided,2);
  if (twosided!=nep->twosided) {
    nep->twosided = twosided;
    nep->state    = NEP_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetTwoSided - Returns the flag indicating whether a two-sided variant
   of the algorithm is being used or not.

   Not Collective

   Input Parameter:
.  nep - the eigensolver context

   Output Parameter:
.  twosided - the returned flag

   Level: advanced

.seealso: NEPSetTwoSided()
@*/
PetscErrorCode NEPGetTwoSided(NEP nep,PetscBool *twosided)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(twosided,2);
  *twosided = nep->twosided;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPSetConvergenceTestFunction - Sets a function to compute the error estimate
   used in the convergence test.

   Logically Collective

   Input Parameters:
+  nep     - nonlinear eigensolver context obtained from NEPCreate()
.  conv    - convergence test function, see NEPConvergenceTestFn for the calling sequence
.  ctx     - context for private data for the convergence routine (may be NULL)
-  destroy - a routine for destroying the context (may be NULL), see PetscCtxDestroyFn for the calling sequence

   Note:
   If the error estimate returned by the convergence test function is less than
   the tolerance, then the eigenvalue is accepted as converged.

   Level: advanced

.seealso: NEPSetConvergenceTest(), NEPSetTolerances()
@*/
PetscErrorCode NEPSetConvergenceTestFunction(NEP nep,NEPConvergenceTestFn *conv,void* ctx,PetscCtxDestroyFn *destroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (nep->convergeddestroy) PetscCall((*nep->convergeddestroy)(&nep->convergedctx));
  nep->convergeduser    = conv;
  nep->convergeddestroy = destroy;
  nep->convergedctx     = ctx;
  if (conv == NEPConvergedRelative) nep->conv = NEP_CONV_REL;
  else if (conv == NEPConvergedNorm) nep->conv = NEP_CONV_NORM;
  else if (conv == NEPConvergedAbsolute) nep->conv = NEP_CONV_ABS;
  else {
    nep->conv      = NEP_CONV_USER;
    nep->converged = nep->convergeduser;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetConvergenceTest - Specifies how to compute the error estimate
   used in the convergence test.

   Logically Collective

   Input Parameters:
+  nep  - nonlinear eigensolver context obtained from NEPCreate()
-  conv - the type of convergence test

   Options Database Keys:
+  -nep_conv_abs  - Sets the absolute convergence test
.  -nep_conv_rel  - Sets the convergence test relative to the eigenvalue
-  -nep_conv_user - Selects the user-defined convergence test

   Note:
   The parameter 'conv' can have one of these values
+     NEP_CONV_ABS  - absolute error ||r||
.     NEP_CONV_REL  - error relative to the eigenvalue l, ||r||/|l|
.     NEP_CONV_NORM - error relative matrix norms, ||r||/sum_i(|f_i(l)|*||A_i||)
-     NEP_CONV_USER - function set by NEPSetConvergenceTestFunction()

   Level: intermediate

.seealso: NEPGetConvergenceTest(), NEPSetConvergenceTestFunction(), NEPSetStoppingTest(), NEPConv
@*/
PetscErrorCode NEPSetConvergenceTest(NEP nep,NEPConv conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(nep,conv,2);
  switch (conv) {
    case NEP_CONV_ABS:  nep->converged = NEPConvergedAbsolute; break;
    case NEP_CONV_REL:  nep->converged = NEPConvergedRelative; break;
    case NEP_CONV_NORM: nep->converged = NEPConvergedNorm; break;
    case NEP_CONV_USER:
      PetscCheck(nep->convergeduser,PetscObjectComm((PetscObject)nep),PETSC_ERR_ORDER,"Must call NEPSetConvergenceTestFunction() first");
      nep->converged = nep->convergeduser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'conv' value");
  }
  nep->conv = conv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetConvergenceTest - Gets the method used to compute the error estimate
   used in the convergence test.

   Not Collective

   Input Parameters:
.  nep   - nonlinear eigensolver context obtained from NEPCreate()

   Output Parameters:
.  conv  - the type of convergence test

   Level: intermediate

.seealso: NEPSetConvergenceTest(), NEPConv
@*/
PetscErrorCode NEPGetConvergenceTest(NEP nep,NEPConv *conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(conv,2);
  *conv = nep->conv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPSetStoppingTestFunction - Sets a function to decide when to stop the outer
   iteration of the eigensolver.

   Logically Collective

   Input Parameters:
+  nep     - nonlinear eigensolver context obtained from NEPCreate()
.  stop    - the stopping test function, see NEPStoppingTestFn for the calling sequence
.  ctx     - context for private data for the stopping routine (may be NULL)
-  destroy - a routine for destroying the context (may be NULL), see PetscCtxDestroyFn for the calling sequence

   Note:
   Normal usage is to first call the default routine NEPStoppingBasic() and then
   set reason to NEP_CONVERGED_USER if some user-defined conditions have been
   met. To let the eigensolver continue iterating, the result must be left as
   NEP_CONVERGED_ITERATING.

   Level: advanced

.seealso: NEPSetStoppingTest(), NEPStoppingBasic()
@*/
PetscErrorCode NEPSetStoppingTestFunction(NEP nep,NEPStoppingTestFn *stop,void* ctx,PetscCtxDestroyFn *destroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (nep->stoppingdestroy) PetscCall((*nep->stoppingdestroy)(&nep->stoppingctx));
  nep->stoppinguser    = stop;
  nep->stoppingdestroy = destroy;
  nep->stoppingctx     = ctx;
  if (stop == NEPStoppingBasic) nep->stop = NEP_STOP_BASIC;
  else {
    nep->stop     = NEP_STOP_USER;
    nep->stopping = nep->stoppinguser;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetStoppingTest - Specifies how to decide the termination of the outer
   loop of the eigensolver.

   Logically Collective

   Input Parameters:
+  nep  - nonlinear eigensolver context obtained from NEPCreate()
-  stop - the type of stopping test

   Options Database Keys:
+  -nep_stop_basic - Sets the default stopping test
-  -nep_stop_user  - Selects the user-defined stopping test

   Note:
   The parameter 'stop' can have one of these values
+     NEP_STOP_BASIC - default stopping test
-     NEP_STOP_USER  - function set by NEPSetStoppingTestFunction()

   Level: advanced

.seealso: NEPGetStoppingTest(), NEPSetStoppingTestFunction(), NEPSetConvergenceTest(), NEPStop
@*/
PetscErrorCode NEPSetStoppingTest(NEP nep,NEPStop stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(nep,stop,2);
  switch (stop) {
    case NEP_STOP_BASIC: nep->stopping = NEPStoppingBasic; break;
    case NEP_STOP_USER:
      PetscCheck(nep->stoppinguser,PetscObjectComm((PetscObject)nep),PETSC_ERR_ORDER,"Must call NEPSetStoppingTestFunction() first");
      nep->stopping = nep->stoppinguser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'stop' value");
  }
  nep->stop = stop;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetStoppingTest - Gets the method used to decide the termination of the outer
   loop of the eigensolver.

   Not Collective

   Input Parameters:
.  nep   - nonlinear eigensolver context obtained from NEPCreate()

   Output Parameters:
.  stop  - the type of stopping test

   Level: advanced

.seealso: NEPSetStoppingTest(), NEPStop
@*/
PetscErrorCode NEPGetStoppingTest(NEP nep,NEPStop *stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(stop,2);
  *stop = nep->stop;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetTrackAll - Specifies if the solver must compute the residual of all
   approximate eigenpairs or not.

   Logically Collective

   Input Parameters:
+  nep      - the eigensolver context
-  trackall - whether compute all residuals or not

   Notes:
   If the user sets trackall=PETSC_TRUE then the solver explicitly computes
   the residual for each eigenpair approximation. Computing the residual is
   usually an expensive operation and solvers commonly compute the associated
   residual to the first unconverged eigenpair.

   The option '-nep_monitor_all' automatically activates this option.

   Level: developer

.seealso: NEPGetTrackAll()
@*/
PetscErrorCode NEPSetTrackAll(NEP nep,PetscBool trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,trackall,2);
  nep->trackall = trackall;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetTrackAll - Returns the flag indicating whether all residual norms must
   be computed or not.

   Not Collective

   Input Parameter:
.  nep - the eigensolver context

   Output Parameter:
.  trackall - the returned flag

   Level: developer

.seealso: NEPSetTrackAll()
@*/
PetscErrorCode NEPGetTrackAll(NEP nep,PetscBool *trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(trackall,2);
  *trackall = nep->trackall;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetRefine - Specifies the refinement type (and options) to be used
   after the solve.

   Logically Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  refine - refinement type
.  npart  - number of partitions of the communicator
.  tol    - the convergence tolerance
.  its    - maximum number of refinement iterations
-  scheme - which scheme to be used for solving the involved linear systems

   Options Database Keys:
+  -nep_refine <type> - refinement type, one of <none,simple,multiple>
.  -nep_refine_partitions <n> - the number of partitions
.  -nep_refine_tol <tol> - the tolerance
.  -nep_refine_its <its> - number of iterations
-  -nep_refine_scheme - to set the scheme for the linear solves

   Notes:
   By default, iterative refinement is disabled, since it may be very
   costly. There are two possible refinement strategies, simple and multiple.
   The simple approach performs iterative refinement on each of the
   converged eigenpairs individually, whereas the multiple strategy works
   with the invariant pair as a whole, refining all eigenpairs simultaneously.
   The latter may be required for the case of multiple eigenvalues.

   In some cases, especially when using direct solvers within the
   iterative refinement method, it may be helpful for improved scalability
   to split the communicator in several partitions. The npart parameter
   indicates how many partitions to use (defaults to 1).

   The tol and its parameters specify the stopping criterion. In the simple
   method, refinement continues until the residual of each eigenpair is
   below the tolerance (tol defaults to the NEP tol, but may be set to a
   different value). In contrast, the multiple method simply performs its
   refinement iterations (just one by default).

   The scheme argument is used to change the way in which linear systems are
   solved. Possible choices are explicit, mixed block elimination (MBE),
   and Schur complement.

   Use PETSC_CURRENT to retain the current value of npart, tol or its. Use
   PETSC_DETERMINE to assign a default value.

   Level: intermediate

.seealso: NEPGetRefine()
@*/
PetscErrorCode NEPSetRefine(NEP nep,NEPRefine refine,PetscInt npart,PetscReal tol,PetscInt its,NEPRefineScheme scheme)
{
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(nep,refine,2);
  PetscValidLogicalCollectiveInt(nep,npart,3);
  PetscValidLogicalCollectiveReal(nep,tol,4);
  PetscValidLogicalCollectiveInt(nep,its,5);
  PetscValidLogicalCollectiveEnum(nep,scheme,6);
  nep->refine = refine;
  if (refine) {  /* process parameters only if not REFINE_NONE */
    if (npart!=nep->npart) {
      PetscCall(PetscSubcommDestroy(&nep->refinesubc));
      PetscCall(KSPDestroy(&nep->refineksp));
    }
    if (npart == PETSC_DETERMINE) {
      nep->npart = 1;
    } else if (npart != PETSC_CURRENT) {
      PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)nep),&size));
      PetscCheck(npart>0 && npart<=size,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of npart");
      nep->npart = npart;
    }
    if (tol == (PetscReal)PETSC_DETERMINE) {
      nep->rtol = PETSC_DETERMINE;
    } else if (tol != (PetscReal)PETSC_CURRENT) {
      PetscCheck(tol>0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
      nep->rtol = tol;
    }
    if (its==PETSC_DETERMINE) {
      nep->rits = PETSC_DETERMINE;
    } else if (its != PETSC_CURRENT) {
      PetscCheck(its>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of its. Must be >= 0");
      nep->rits = its;
    }
    nep->scheme = scheme;
  }
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetRefine - Gets the refinement strategy used by the NEP object, and the
   associated parameters.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  refine - refinement type
.  npart  - number of partitions of the communicator
.  tol    - the convergence tolerance
.  its    - maximum number of refinement iterations
-  scheme - the scheme used for solving linear systems

   Level: intermediate

   Note:
   The user can specify NULL for any parameter that is not needed.

.seealso: NEPSetRefine()
@*/
PetscErrorCode NEPGetRefine(NEP nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (refine) *refine = nep->refine;
  if (npart)  *npart  = nep->npart;
  if (tol)    *tol    = nep->rtol;
  if (its)    *its    = nep->rits;
  if (scheme) *scheme = nep->scheme;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetOptionsPrefix - Sets the prefix used for searching for all
   NEP options in the database.

   Logically Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  prefix - the prefix string to prepend to all NEP option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different NEP contexts, one could call
.vb
      NEPSetOptionsPrefix(nep1,"neig1_")
      NEPSetOptionsPrefix(nep2,"neig2_")
.ve

   Level: advanced

.seealso: NEPAppendOptionsPrefix(), NEPGetOptionsPrefix()
@*/
PetscErrorCode NEPSetOptionsPrefix(NEP nep,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!nep->V) PetscCall(NEPGetBV(nep,&nep->V));
  PetscCall(BVSetOptionsPrefix(nep->V,prefix));
  if (!nep->ds) PetscCall(NEPGetDS(nep,&nep->ds));
  PetscCall(DSSetOptionsPrefix(nep->ds,prefix));
  if (!nep->rg) PetscCall(NEPGetRG(nep,&nep->rg));
  PetscCall(RGSetOptionsPrefix(nep->rg,prefix));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)nep,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPAppendOptionsPrefix - Appends to the prefix used for searching for all
   NEP options in the database.

   Logically Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  prefix - the prefix string to prepend to all NEP option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: NEPSetOptionsPrefix(), NEPGetOptionsPrefix()
@*/
PetscErrorCode NEPAppendOptionsPrefix(NEP nep,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!nep->V) PetscCall(NEPGetBV(nep,&nep->V));
  PetscCall(BVAppendOptionsPrefix(nep->V,prefix));
  if (!nep->ds) PetscCall(NEPGetDS(nep,&nep->ds));
  PetscCall(DSAppendOptionsPrefix(nep->ds,prefix));
  if (!nep->rg) PetscCall(NEPGetRG(nep,&nep->rg));
  PetscCall(RGAppendOptionsPrefix(nep->rg,prefix));
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)nep,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetOptionsPrefix - Gets the prefix used for searching for all
   NEP options in the database.

   Not Collective

   Input Parameters:
.  nep - the nonlinear eigensolver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: NEPSetOptionsPrefix(), NEPAppendOptionsPrefix()
@*/
PetscErrorCode NEPGetOptionsPrefix(NEP nep,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)nep,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}
