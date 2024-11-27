/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   EPS routines related to options that can be set via the command-line
   or procedurally.
*/

#include <slepc/private/epsimpl.h>   /*I "slepceps.h" I*/
#include <petscdraw.h>

/*@C
   EPSMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective

   Input Parameters:
+  eps      - the eigensolver context
.  opt      - the command line option for this monitor
.  name     - the monitor type one is seeking
.  ctx      - an optional user context for the monitor, or NULL
-  trackall - whether this monitor tracks all eigenvalues or not

   Level: developer

.seealso: EPSMonitorSet(), EPSSetTrackAll()
@*/
PetscErrorCode EPSMonitorSetFromOptions(EPS eps,const char opt[],const char name[],void *ctx,PetscBool trackall)
{
  PetscErrorCode       (*mfunc)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
  PetscErrorCode       (*cfunc)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
  PetscErrorCode       (*dfunc)(PetscViewerAndFormat**);
  PetscViewerAndFormat *vf;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerType      vtype;
  char                 key[PETSC_MAX_PATH_LEN];
  PetscBool            flg;

  PetscFunctionBegin;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,opt,&viewer,&format,&flg));
  if (!flg) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetType(viewer,&vtype));
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListFind(EPSMonitorList,key,&mfunc));
  PetscCheck(mfunc,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Specified viewer and format not supported");
  PetscCall(PetscFunctionListFind(EPSMonitorCreateList,key,&cfunc));
  PetscCall(PetscFunctionListFind(EPSMonitorDestroyList,key,&dfunc));
  if (!cfunc) cfunc = PetscViewerAndFormatCreate_Internal;
  if (!dfunc) dfunc = PetscViewerAndFormatDestroy;

  PetscCall((*cfunc)(viewer,format,ctx,&vf));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(EPSMonitorSet(eps,mfunc,vf,(PetscCtxDestroyFn*)dfunc));
  if (trackall) PetscCall(EPSSetTrackAll(eps,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetFromOptions - Sets EPS options from the options database.
   This routine must be called before EPSSetUp() if the user is to be
   allowed to set the solver type.

   Collective

   Input Parameters:
.  eps - the eigensolver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner

.seealso: EPSSetOptionsPrefix()
@*/
PetscErrorCode EPSSetFromOptions(EPS eps)
{
  char           type[256];
  PetscBool      set,flg,flg1,flg2,flg3,bval;
  PetscReal      r,array[2]={0,0};
  PetscScalar    s;
  PetscInt       i,j,k;
  EPSBalance     bal;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscCall(EPSRegisterAll());
  PetscObjectOptionsBegin((PetscObject)eps);
    PetscCall(PetscOptionsFList("-eps_type","Eigensolver method","EPSSetType",EPSList,(char*)(((PetscObject)eps)->type_name?((PetscObject)eps)->type_name:EPSKRYLOVSCHUR),type,sizeof(type),&flg));
    if (flg) PetscCall(EPSSetType(eps,type));
    else if (!((PetscObject)eps)->type_name) PetscCall(EPSSetType(eps,EPSKRYLOVSCHUR));

    PetscCall(PetscOptionsBoolGroupBegin("-eps_hermitian","Hermitian eigenvalue problem","EPSSetProblemType",&flg));
    if (flg) PetscCall(EPSSetProblemType(eps,EPS_HEP));
    PetscCall(PetscOptionsBoolGroup("-eps_gen_hermitian","Generalized Hermitian eigenvalue problem","EPSSetProblemType",&flg));
    if (flg) PetscCall(EPSSetProblemType(eps,EPS_GHEP));
    PetscCall(PetscOptionsBoolGroup("-eps_non_hermitian","Non-Hermitian eigenvalue problem","EPSSetProblemType",&flg));
    if (flg) PetscCall(EPSSetProblemType(eps,EPS_NHEP));
    PetscCall(PetscOptionsBoolGroup("-eps_gen_non_hermitian","Generalized non-Hermitian eigenvalue problem","EPSSetProblemType",&flg));
    if (flg) PetscCall(EPSSetProblemType(eps,EPS_GNHEP));
    PetscCall(PetscOptionsBoolGroup("-eps_pos_gen_non_hermitian","Generalized non-Hermitian eigenvalue problem with positive semi-definite B","EPSSetProblemType",&flg));
    if (flg) PetscCall(EPSSetProblemType(eps,EPS_PGNHEP));
    PetscCall(PetscOptionsBoolGroup("-eps_gen_indefinite","Generalized Hermitian-indefinite eigenvalue problem","EPSSetProblemType",&flg));
    if (flg) PetscCall(EPSSetProblemType(eps,EPS_GHIEP));
    PetscCall(PetscOptionsBoolGroupEnd("-eps_bse","Structured Bethe-Salpeter eigenvalue problem","EPSSetProblemType",&flg));
    if (flg) PetscCall(EPSSetProblemType(eps,EPS_BSE));

    PetscCall(PetscOptionsBoolGroupBegin("-eps_ritz","Rayleigh-Ritz extraction","EPSSetExtraction",&flg));
    if (flg) PetscCall(EPSSetExtraction(eps,EPS_RITZ));
    PetscCall(PetscOptionsBoolGroup("-eps_harmonic","Harmonic Ritz extraction","EPSSetExtraction",&flg));
    if (flg) PetscCall(EPSSetExtraction(eps,EPS_HARMONIC));
    PetscCall(PetscOptionsBoolGroup("-eps_harmonic_relative","Relative harmonic Ritz extraction","EPSSetExtraction",&flg));
    if (flg) PetscCall(EPSSetExtraction(eps,EPS_HARMONIC_RELATIVE));
    PetscCall(PetscOptionsBoolGroup("-eps_harmonic_right","Right harmonic Ritz extraction","EPSSetExtraction",&flg));
    if (flg) PetscCall(EPSSetExtraction(eps,EPS_HARMONIC_RIGHT));
    PetscCall(PetscOptionsBoolGroup("-eps_harmonic_largest","Largest harmonic Ritz extraction","EPSSetExtraction",&flg));
    if (flg) PetscCall(EPSSetExtraction(eps,EPS_HARMONIC_LARGEST));
    PetscCall(PetscOptionsBoolGroup("-eps_refined","Refined Ritz extraction","EPSSetExtraction",&flg));
    if (flg) PetscCall(EPSSetExtraction(eps,EPS_REFINED));
    PetscCall(PetscOptionsBoolGroupEnd("-eps_refined_harmonic","Refined harmonic Ritz extraction","EPSSetExtraction",&flg));
    if (flg) PetscCall(EPSSetExtraction(eps,EPS_REFINED_HARMONIC));

    bal = eps->balance;
    PetscCall(PetscOptionsEnum("-eps_balance","Balancing method","EPSSetBalance",EPSBalanceTypes,(PetscEnum)bal,(PetscEnum*)&bal,&flg1));
    j = eps->balance_its;
    PetscCall(PetscOptionsInt("-eps_balance_its","Number of iterations in balancing","EPSSetBalance",eps->balance_its,&j,&flg2));
    r = eps->balance_cutoff;
    PetscCall(PetscOptionsReal("-eps_balance_cutoff","Cutoff value in balancing","EPSSetBalance",eps->balance_cutoff,&r,&flg3));
    if (flg1 || flg2 || flg3) PetscCall(EPSSetBalance(eps,bal,j,r));

    i = eps->max_it;
    PetscCall(PetscOptionsInt("-eps_max_it","Maximum number of iterations","EPSSetTolerances",eps->max_it,&i,&flg1));
    r = eps->tol;
    PetscCall(PetscOptionsReal("-eps_tol","Tolerance","EPSSetTolerances",SlepcDefaultTol(eps->tol),&r,&flg2));
    if (flg1 || flg2) PetscCall(EPSSetTolerances(eps,r,i));

    r = eps->thres;
    PetscCall(PetscOptionsReal("-eps_threshold_absolute","Absolute threshold","EPSSetThreshold",r,&r,&flg));
    if (flg) PetscCall(EPSSetThreshold(eps,r,PETSC_FALSE));
    PetscCall(PetscOptionsReal("-eps_threshold_relative","Relative threshold","EPSSetThreshold",r,&r,&flg));
    if (flg) PetscCall(EPSSetThreshold(eps,r,PETSC_TRUE));

    PetscCall(PetscOptionsBoolGroupBegin("-eps_conv_rel","Relative error convergence test","EPSSetConvergenceTest",&flg));
    if (flg) PetscCall(EPSSetConvergenceTest(eps,EPS_CONV_REL));
    PetscCall(PetscOptionsBoolGroup("-eps_conv_norm","Convergence test relative to the eigenvalue and the matrix norms","EPSSetConvergenceTest",&flg));
    if (flg) PetscCall(EPSSetConvergenceTest(eps,EPS_CONV_NORM));
    PetscCall(PetscOptionsBoolGroup("-eps_conv_abs","Absolute error convergence test","EPSSetConvergenceTest",&flg));
    if (flg) PetscCall(EPSSetConvergenceTest(eps,EPS_CONV_ABS));
    PetscCall(PetscOptionsBoolGroupEnd("-eps_conv_user","User-defined convergence test","EPSSetConvergenceTest",&flg));
    if (flg) PetscCall(EPSSetConvergenceTest(eps,EPS_CONV_USER));

    PetscCall(PetscOptionsBoolGroupBegin("-eps_stop_basic","Stop iteration if all eigenvalues converged or max_it reached","EPSSetStoppingTest",&flg));
    if (flg) PetscCall(EPSSetStoppingTest(eps,EPS_STOP_BASIC));
    PetscCall(PetscOptionsBoolGroup("-eps_stop_threshold","Stop iteration if a converged eigenvalue is below/above the threshold","EPSSetStoppingTest",&flg));
    if (flg) PetscCall(EPSSetStoppingTest(eps,EPS_STOP_THRESHOLD));
    PetscCall(PetscOptionsBoolGroupEnd("-eps_stop_user","User-defined stopping test","EPSSetStoppingTest",&flg));
    if (flg) PetscCall(EPSSetStoppingTest(eps,EPS_STOP_USER));

    i = eps->nev;
    PetscCall(PetscOptionsInt("-eps_nev","Number of eigenvalues to compute","EPSSetDimensions",eps->nev,&i,&flg1));
    if (!flg1) i = PETSC_CURRENT;
    j = eps->ncv;
    PetscCall(PetscOptionsInt("-eps_ncv","Number of basis vectors","EPSSetDimensions",eps->ncv,&j,&flg2));
    k = eps->mpd;
    PetscCall(PetscOptionsInt("-eps_mpd","Maximum dimension of projected problem","EPSSetDimensions",eps->mpd,&k,&flg3));
    if (flg1 || flg2 || flg3) PetscCall(EPSSetDimensions(eps,i,j,k));

    PetscCall(PetscOptionsBoolGroupBegin("-eps_largest_magnitude","Compute largest eigenvalues in magnitude","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-eps_smallest_magnitude","Compute smallest eigenvalues in magnitude","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-eps_largest_real","Compute eigenvalues with largest real parts","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
    PetscCall(PetscOptionsBoolGroup("-eps_smallest_real","Compute eigenvalues with smallest real parts","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));
    PetscCall(PetscOptionsBoolGroup("-eps_largest_imaginary","Compute eigenvalues with largest imaginary parts","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_IMAGINARY));
    PetscCall(PetscOptionsBoolGroup("-eps_smallest_imaginary","Compute eigenvalues with smallest imaginary parts","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_IMAGINARY));
    PetscCall(PetscOptionsBoolGroup("-eps_target_magnitude","Compute eigenvalues closest to target","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-eps_target_real","Compute eigenvalues with real parts closest to target","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_REAL));
    PetscCall(PetscOptionsBoolGroup("-eps_target_imaginary","Compute eigenvalues with imaginary parts closest to target","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_IMAGINARY));
    PetscCall(PetscOptionsBoolGroupEnd("-eps_all","Compute all eigenvalues in an interval or a region","EPSSetWhichEigenpairs",&flg));
    if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_ALL));

    PetscCall(PetscOptionsScalar("-eps_target","Value of the target","EPSSetTarget",eps->target,&s,&flg));
    if (flg) {
      if (eps->which!=EPS_TARGET_REAL && eps->which!=EPS_TARGET_IMAGINARY) PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));
      PetscCall(EPSSetTarget(eps,s));
    }

    k = 2;
    PetscCall(PetscOptionsRealArray("-eps_interval","Computational interval (two real values separated with a comma without spaces)","EPSSetInterval",array,&k,&flg));
    if (flg) {
      PetscCheck(k>1,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_SIZ,"Must pass two values in -eps_interval (comma-separated without spaces)");
      PetscCall(EPSSetWhichEigenpairs(eps,EPS_ALL));
      PetscCall(EPSSetInterval(eps,array[0],array[1]));
    }

    PetscCall(PetscOptionsBool("-eps_true_residual","Compute true residuals explicitly","EPSSetTrueResidual",eps->trueres,&eps->trueres,NULL));
    PetscCall(PetscOptionsBool("-eps_purify","Postprocess eigenvectors for purification","EPSSetPurify",eps->purify,&bval,&flg));
    if (flg) PetscCall(EPSSetPurify(eps,bval));
    PetscCall(PetscOptionsBool("-eps_two_sided","Use two-sided variant (to compute left eigenvectors)","EPSSetTwoSided",eps->twosided,&bval,&flg));
    if (flg) PetscCall(EPSSetTwoSided(eps,bval));

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to EPSSetFromOptions()
    */
    PetscCall(PetscOptionsBool("-eps_monitor_cancel","Remove any hardwired monitor routines","EPSMonitorCancel",PETSC_FALSE,&flg,&set));
    if (set && flg) PetscCall(EPSMonitorCancel(eps));
    PetscCall(EPSMonitorSetFromOptions(eps,"-eps_monitor","first_approximation",NULL,PETSC_FALSE));
    PetscCall(EPSMonitorSetFromOptions(eps,"-eps_monitor_all","all_approximations",NULL,PETSC_TRUE));
    PetscCall(EPSMonitorSetFromOptions(eps,"-eps_monitor_conv","convergence_history",NULL,PETSC_FALSE));

    /* -----------------------------------------------------------------------*/
    PetscCall(PetscOptionsName("-eps_view","Print detailed information on solver used","EPSView",&set));
    PetscCall(PetscOptionsName("-eps_view_vectors","View computed eigenvectors","EPSVectorsView",&set));
    PetscCall(PetscOptionsName("-eps_view_values","View computed eigenvalues","EPSValuesView",&set));
    PetscCall(PetscOptionsName("-eps_converged_reason","Print reason for convergence, and number of iterations","EPSConvergedReasonView",&set));
    PetscCall(PetscOptionsName("-eps_error_absolute","Print absolute errors of each eigenpair","EPSErrorView",&set));
    PetscCall(PetscOptionsName("-eps_error_relative","Print relative errors of each eigenpair","EPSErrorView",&set));
    PetscCall(PetscOptionsName("-eps_error_backward","Print backward errors of each eigenpair","EPSErrorView",&set));

    PetscTryTypeMethod(eps,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)eps,PetscOptionsObject));
  PetscOptionsEnd();

  if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
  PetscCall(BVSetFromOptions(eps->V));
  if (!eps->rg) PetscCall(EPSGetRG(eps,&eps->rg));
  PetscCall(RGSetFromOptions(eps->rg));
  if (eps->useds) {
    if (!eps->ds) PetscCall(EPSGetDS(eps,&eps->ds));
    PetscCall(EPSSetDSType(eps));
    PetscCall(DSSetFromOptions(eps->ds));
  }
  if (!eps->st) PetscCall(EPSGetST(eps,&eps->st));
  PetscCall(EPSSetDefaultST(eps));
  PetscCall(STSetFromOptions(eps->st));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetTolerances - Gets the tolerance and maximum iteration count used
   by the EPS convergence tests.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: EPSSetTolerances()
@*/
PetscErrorCode EPSGetTolerances(EPS eps,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (tol)    *tol    = eps->tol;
  if (maxits) *maxits = eps->max_it;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetTolerances - Sets the tolerance and maximum iteration count used
   by the EPS convergence tests.

   Logically Collective

   Input Parameters:
+  eps - the eigensolver context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -eps_tol <tol> - Sets the convergence tolerance
-  -eps_max_it <maxits> - Sets the maximum number of iterations allowed

   Notes:
   Use PETSC_CURRENT to retain the current value of any of the parameters.
   Use PETSC_DETERMINE for either argument to assign a default value computed
   internally (may be different in each solver).
   For maxits use PETSC_UMLIMITED to indicate there is no upper bound on this value.

   Level: intermediate

.seealso: EPSGetTolerances()
@*/
PetscErrorCode EPSSetTolerances(EPS eps,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,tol,2);
  PetscValidLogicalCollectiveInt(eps,maxits,3);
  if (tol == (PetscReal)PETSC_DETERMINE) {
    eps->tol   = PETSC_DETERMINE;
    eps->state = EPS_STATE_INITIAL;
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    eps->tol = tol;
  }
  if (maxits == PETSC_DETERMINE) {
    eps->max_it = PETSC_DETERMINE;
    eps->state  = EPS_STATE_INITIAL;
  } else if (maxits == PETSC_UNLIMITED) {
    eps->max_it = PETSC_INT_MAX;
  } else if (maxits != PETSC_CURRENT) {
    PetscCheck(maxits>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    eps->max_it = maxits;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetDimensions - Gets the number of eigenvalues to compute
   and the dimension of the subspace.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Level: intermediate

.seealso: EPSSetDimensions()
@*/
PetscErrorCode EPSGetDimensions(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (nev) *nev = eps->nev? eps->nev: 1;
  if (ncv) *ncv = eps->ncv;
  if (mpd) *mpd = eps->mpd;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetDimensions - Sets the number of eigenvalues to compute
   and the dimension of the subspace.

   Logically Collective

   Input Parameters:
+  eps - the eigensolver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Keys:
+  -eps_nev <nev> - Sets the number of eigenvalues
.  -eps_ncv <ncv> - Sets the dimension of the subspace
-  -eps_mpd <mpd> - Sets the maximum projected dimension

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

   When computing all eigenvalues in an interval, see EPSSetInterval(), these
   parameters lose relevance, and tuning must be done with
   EPSKrylovSchurSetDimensions().

   Level: intermediate

.seealso: EPSGetDimensions(), EPSSetInterval(), EPSKrylovSchurSetDimensions()
@*/
PetscErrorCode EPSSetDimensions(EPS eps,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nev,2);
  PetscValidLogicalCollectiveInt(eps,ncv,3);
  PetscValidLogicalCollectiveInt(eps,mpd,4);
  if (nev != PETSC_CURRENT) {
    PetscCheck(nev>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
    eps->nev = nev;
  }
  if (ncv == PETSC_DETERMINE) {
    eps->ncv = PETSC_DETERMINE;
  } else if (ncv != PETSC_CURRENT) {
    PetscCheck(ncv>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    eps->ncv = ncv;
  }
  if (mpd == PETSC_DETERMINE) {
    eps->mpd = PETSC_DETERMINE;
  } else if (mpd != PETSC_CURRENT) {
    PetscCheck(mpd>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    eps->mpd = mpd;
  }
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetWhichEigenpairs - Specifies which portion of the spectrum is
   to be sought.

   Logically Collective

   Input Parameters:
+  eps   - eigensolver context obtained from EPSCreate()
-  which - the portion of the spectrum to be sought

   Options Database Keys:
+   -eps_largest_magnitude - Sets largest eigenvalues in magnitude
.   -eps_smallest_magnitude - Sets smallest eigenvalues in magnitude
.   -eps_largest_real - Sets largest real parts
.   -eps_smallest_real - Sets smallest real parts
.   -eps_largest_imaginary - Sets largest imaginary parts
.   -eps_smallest_imaginary - Sets smallest imaginary parts
.   -eps_target_magnitude - Sets eigenvalues closest to target
.   -eps_target_real - Sets real parts closest to target
.   -eps_target_imaginary - Sets imaginary parts closest to target
-   -eps_all - Sets all eigenvalues in an interval or region

   Notes:
   The parameter 'which' can have one of these values

+     EPS_LARGEST_MAGNITUDE - largest eigenvalues in magnitude (default)
.     EPS_SMALLEST_MAGNITUDE - smallest eigenvalues in magnitude
.     EPS_LARGEST_REAL - largest real parts
.     EPS_SMALLEST_REAL - smallest real parts
.     EPS_LARGEST_IMAGINARY - largest imaginary parts
.     EPS_SMALLEST_IMAGINARY - smallest imaginary parts
.     EPS_TARGET_MAGNITUDE - eigenvalues closest to the target (in magnitude)
.     EPS_TARGET_REAL - eigenvalues with real part closest to target
.     EPS_TARGET_IMAGINARY - eigenvalues with imaginary part closest to target
.     EPS_ALL - all eigenvalues contained in a given interval or region
-     EPS_WHICH_USER - user defined ordering set with EPSSetEigenvalueComparison()

   Not all eigensolvers implemented in EPS account for all the possible values
   stated above. Also, some values make sense only for certain types of
   problems. If SLEPc is compiled for real numbers EPS_LARGEST_IMAGINARY
   and EPS_SMALLEST_IMAGINARY use the absolute value of the imaginary part
   for eigenvalue selection.

   The target is a scalar value provided with EPSSetTarget().

   The criterion EPS_TARGET_IMAGINARY is available only in case PETSc and
   SLEPc have been built with complex scalars.

   EPS_ALL is intended for use in combination with an interval (see
   EPSSetInterval()), when all eigenvalues within the interval are requested,
   or in the context of the CISS solver for computing all eigenvalues in a region.
   In those cases, the number of eigenvalues is unknown, so the nev parameter
   has a different sense, see EPSSetDimensions().

   Level: intermediate

.seealso: EPSGetWhichEigenpairs(), EPSSetTarget(), EPSSetInterval(),
          EPSSetDimensions(), EPSSetEigenvalueComparison(), EPSWhich
@*/
PetscErrorCode EPSSetWhichEigenpairs(EPS eps,EPSWhich which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,which,2);
  switch (which) {
    case EPS_LARGEST_MAGNITUDE:
    case EPS_SMALLEST_MAGNITUDE:
    case EPS_LARGEST_REAL:
    case EPS_SMALLEST_REAL:
    case EPS_LARGEST_IMAGINARY:
    case EPS_SMALLEST_IMAGINARY:
    case EPS_TARGET_MAGNITUDE:
    case EPS_TARGET_REAL:
#if defined(PETSC_USE_COMPLEX)
    case EPS_TARGET_IMAGINARY:
#endif
    case EPS_ALL:
    case EPS_WHICH_USER:
      if (eps->which != which) {
        eps->state = EPS_STATE_INITIAL;
        eps->which = which;
      }
      break;
#if !defined(PETSC_USE_COMPLEX)
    case EPS_TARGET_IMAGINARY:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"EPS_TARGET_IMAGINARY can be used only with complex scalars");
#endif
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'which' value");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetWhichEigenpairs - Returns which portion of the spectrum is to be
   sought.

   Not Collective

   Input Parameter:
.  eps - eigensolver context obtained from EPSCreate()

   Output Parameter:
.  which - the portion of the spectrum to be sought

   Notes:
   See EPSSetWhichEigenpairs() for possible values of 'which'.

   Level: intermediate

.seealso: EPSSetWhichEigenpairs(), EPSWhich
@*/
PetscErrorCode EPSGetWhichEigenpairs(EPS eps,EPSWhich *which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(which,2);
  *which = eps->which;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetThreshold - Sets the threshold used in the threshold stopping test.

   Logically Collective

   Input Parameters:
+  eps   - the eigenvalue solver context
.  thres - the threshold value
-  rel   - whether the threshold is relative or not

   Options Database Keys:
+  -eps_threshold_absolute <thres> - Sets an absolute threshold
-  -eps_threshold_relative <thres> - Sets a relative threshold

   Notes:
   This function internally calls EPSSetStoppingTest() to set a special stopping
   test based on the threshold, where eigenvalues are computed in sequence until
   one of the computed eigenvalues is below the threshold (in magnitude). This is
   the interpretation in case of searching for largest eigenvalues in magnitude,
   see EPSSetWhichEigenpairs().

   If the solver is configured to compute smallest magnitude eigenvalues, then the
   threshold must be interpreted in the opposite direction, i.e., the computation
   will stop when one of the computed values is above the threshold (in magnitude).

   The threshold can also be used when computing largest/smallest real eigenvalues
   (i.e, rightmost or leftmost), in which case the threshold is allowed to be
   negative. The solver will stop when one of the computed eigenvalues is above
   or below the threshold (considering the real part of the eigenvalue). This mode
   is allowed only in problem types whose eigenvalues are always real (e.g., HEP).

   In the case of largest magnitude eigenvalues, the threshold can be made relative
   with respect to the dominant eigenvalue. Otherwise, the argument rel should be
   PETSC_FALSE.

   An additional use case is with target magnitude selection of eigenvalues (e.g.,
   with shift-and-invert), but this must be used with caution to avoid unexpected
   behaviour. With an absolute threshold, the solver will assume that leftmost
   eigenvalues are being computed (e.g., with target=0 for a problem with real
   positive eigenvalues). In case of a relative threshold, a value of threshold<1
   implies that the wanted eigenvalues are the largest ones, and otherwise the
   solver assumes that smallest eigenvalues are being computed.

   The test against the threshold is done for converged eigenvalues, which
   implies that the final number of converged eigenvalues will be at least
   one more than the actual number of values below/above the threshold.

   Since the number of computed eigenvalues is not known a priori, the solver
   will need to reallocate the basis of vectors internally, to have enough room
   to accommodate all the eigenvectors. Hence, this option must be used with
   caution to avoid out-of-memory problems. The recommendation is to set the value
   of ncv to be larger than the estimated number of eigenvalues, to minimize the
   number of reallocations.

   If a number of wanted eigenvalues has been set with EPSSetDimensions()
   it is also taken into account and the solver will stop when one of the two
   conditions (threshold or number of converged values) is met.

   Use EPSSetStoppingTest() to return to the usual computation of a fixed number
   of eigenvalues.

   Level: advanced

.seealso: EPSGetThreshold(), EPSSetStoppingTest(), EPSSetDimensions(), EPSSetWhichEigenpairs(), EPSSetProblemType()
@*/
PetscErrorCode EPSSetThreshold(EPS eps,PetscReal thres,PetscBool rel)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,thres,2);
  PetscValidLogicalCollectiveBool(eps,rel,3);
  if (eps->thres != thres || eps->threlative != rel) {
    eps->thres = thres;
    eps->threlative = rel;
    eps->state = EPS_STATE_INITIAL;
    PetscCall(EPSSetStoppingTest(eps,EPS_STOP_THRESHOLD));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetThreshold - Gets the threshold used by the threshold stopping test.

   Not Collective

   Input Parameter:
.  eps - the eigenvalue solver context

   Output Parameters:
+  thres - the threshold
-  rel   - whether the threshold is relative or not

   Level: advanced

.seealso: EPSSetThreshold()
@*/
PetscErrorCode EPSGetThreshold(EPS eps,PetscReal *thres,PetscBool *rel)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (thres) *thres = eps->thres;
  if (rel)   *rel   = eps->threlative;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSSetEigenvalueComparison - Specifies the eigenvalue comparison function
   when EPSSetWhichEigenpairs() is set to EPS_WHICH_USER.

   Logically Collective

   Input Parameters:
+  eps  - eigensolver context obtained from EPSCreate()
.  func - the comparison function, see EPSEigenvalueComparisonFn for the calling sequence
-  ctx  - a context pointer (the last parameter to the comparison function)

   Note:
   The returning parameter 'res' can be
+  negative - if the 1st eigenvalue is preferred to the 2st one
.  zero     - if both eigenvalues are equally preferred
-  positive - if the 2st eigenvalue is preferred to the 1st one

   Level: advanced

.seealso: EPSSetWhichEigenpairs(), EPSWhich
@*/
PetscErrorCode EPSSetEigenvalueComparison(EPS eps,SlepcEigenvalueComparisonFn *func,void* ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  eps->sc->comparison    = func;
  eps->sc->comparisonctx = ctx;
  eps->which             = EPS_WHICH_USER;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSSetArbitrarySelection - Specifies a function intended to look for
   eigenvalues according to an arbitrary selection criterion. This criterion
   can be based on a computation involving the current eigenvector approximation.

   Logically Collective

   Input Parameters:
+  eps  - eigensolver context obtained from EPSCreate()
.  func - the arbitrary selection function, see SlepcArbitrarySelectionFn for a calling sequence
-  ctx  - a context pointer (the last parameter to the arbitrary selection function)

   Notes:
   This provides a mechanism to select eigenpairs by evaluating a user-defined
   function. When a function has been provided, the default selection based on
   sorting the eigenvalues is replaced by the sorting of the results of this
   function (with the same sorting criterion given in EPSSetWhichEigenpairs()).

   For instance, suppose you want to compute those eigenvectors that maximize
   a certain computable expression. Then implement the computation using
   the arguments xr and xi, and return the result in rr. Then set the standard
   sorting by magnitude so that the eigenpair with largest value of rr is
   selected.

   This evaluation function is collective, that is, all processes call it and
   it can use collective operations; furthermore, the computed result must
   be the same in all processes.

   The result of func is expressed as a complex number so that it is possible to
   use the standard eigenvalue sorting functions, but normally only rr is used.
   Set ri to zero unless it is meaningful in your application.

   Level: advanced

.seealso: EPSSetWhichEigenpairs()
@*/
PetscErrorCode EPSSetArbitrarySelection(EPS eps,SlepcArbitrarySelectionFn *func,void* ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  eps->arbitrary    = func;
  eps->arbitraryctx = ctx;
  eps->state        = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSSetConvergenceTestFunction - Sets a function to compute the error estimate
   used in the convergence test.

   Logically Collective

   Input Parameters:
+  eps     - eigensolver context obtained from EPSCreate()
.  func    - convergence test function, see EPSConvergenceTestFn for the calling sequence
.  ctx     - context for private data for the convergence routine (may be NULL)
-  destroy - a routine for destroying the context (may be NULL), see PetscCtxDestroyFn for the calling sequence

   Note:
   If the error estimate returned by the convergence test function is less than
   the tolerance, then the eigenvalue is accepted as converged.

   Level: advanced

.seealso: EPSSetConvergenceTest(), EPSSetTolerances()
@*/
PetscErrorCode EPSSetConvergenceTestFunction(EPS eps,EPSConvergenceTestFn *func,void* ctx,PetscCtxDestroyFn *destroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (eps->convergeddestroy) PetscCall((*eps->convergeddestroy)(&eps->convergedctx));
  eps->convergeduser    = func;
  eps->convergeddestroy = destroy;
  eps->convergedctx     = ctx;
  if (func == EPSConvergedRelative) eps->conv = EPS_CONV_REL;
  else if (func == EPSConvergedNorm) eps->conv = EPS_CONV_NORM;
  else if (func == EPSConvergedAbsolute) eps->conv = EPS_CONV_ABS;
  else {
    eps->conv      = EPS_CONV_USER;
    eps->converged = eps->convergeduser;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetConvergenceTest - Specifies how to compute the error estimate
   used in the convergence test.

   Logically Collective

   Input Parameters:
+  eps  - eigensolver context obtained from EPSCreate()
-  conv - the type of convergence test

   Options Database Keys:
+  -eps_conv_abs  - Sets the absolute convergence test
.  -eps_conv_rel  - Sets the convergence test relative to the eigenvalue
.  -eps_conv_norm - Sets the convergence test relative to the matrix norms
-  -eps_conv_user - Selects the user-defined convergence test

   Note:
   The parameter 'conv' can have one of these values
+     EPS_CONV_ABS  - absolute error ||r||
.     EPS_CONV_REL  - error relative to the eigenvalue l, ||r||/|l|
.     EPS_CONV_NORM - error relative to the matrix norms, ||r||/(||A||+|l|*||B||)
-     EPS_CONV_USER - function set by EPSSetConvergenceTestFunction()

   Level: intermediate

.seealso: EPSGetConvergenceTest(), EPSSetConvergenceTestFunction(), EPSSetStoppingTest(), EPSConv
@*/
PetscErrorCode EPSSetConvergenceTest(EPS eps,EPSConv conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,conv,2);
  switch (conv) {
    case EPS_CONV_ABS:  eps->converged = EPSConvergedAbsolute; break;
    case EPS_CONV_REL:  eps->converged = EPSConvergedRelative; break;
    case EPS_CONV_NORM: eps->converged = EPSConvergedNorm; break;
    case EPS_CONV_USER:
      PetscCheck(eps->convergeduser,PetscObjectComm((PetscObject)eps),PETSC_ERR_ORDER,"Must call EPSSetConvergenceTestFunction() first");
      eps->converged = eps->convergeduser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'conv' value");
  }
  eps->conv = conv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetConvergenceTest - Gets the method used to compute the error estimate
   used in the convergence test.

   Not Collective

   Input Parameters:
.  eps   - eigensolver context obtained from EPSCreate()

   Output Parameters:
.  conv  - the type of convergence test

   Level: intermediate

.seealso: EPSSetConvergenceTest(), EPSConv
@*/
PetscErrorCode EPSGetConvergenceTest(EPS eps,EPSConv *conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(conv,2);
  *conv = eps->conv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSSetStoppingTestFunction - Sets a function to decide when to stop the outer
   iteration of the eigensolver.

   Logically Collective

   Input Parameters:
+  eps     - eigensolver context obtained from EPSCreate()
.  func    - stopping test function, see EPSStoppingTestFn for the calling sequence
.  ctx     - context for private data for the stopping routine (may be NULL)
-  destroy - a routine for destroying the context (may be NULL), see PetscCtxDestroyFn for the calling sequence

   Note:
   Normal usage is to first call the default routine EPSStoppingBasic() and then
   set reason to EPS_CONVERGED_USER if some user-defined conditions have been
   met. To let the eigensolver continue iterating, the result must be left as
   EPS_CONVERGED_ITERATING.

   Level: advanced

.seealso: EPSSetStoppingTest(), EPSStoppingBasic()
@*/
PetscErrorCode EPSSetStoppingTestFunction(EPS eps,EPSStoppingTestFn *func,void* ctx,PetscCtxDestroyFn *destroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (eps->stoppingdestroy) PetscCall((*eps->stoppingdestroy)(&eps->stoppingctx));
  eps->stoppinguser    = func;
  eps->stoppingdestroy = destroy;
  eps->stoppingctx     = ctx;
  if (func == EPSStoppingBasic) PetscCall(EPSSetStoppingTest(eps,EPS_STOP_BASIC));
  else if (func == EPSStoppingThreshold) PetscCall(EPSSetStoppingTest(eps,EPS_STOP_THRESHOLD));
  else {
    eps->stop     = EPS_STOP_USER;
    eps->stopping = eps->stoppinguser;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetStoppingTest - Specifies how to decide the termination of the outer
   loop of the eigensolver.

   Logically Collective

   Input Parameters:
+  eps  - eigensolver context obtained from EPSCreate()
-  stop - the type of stopping test

   Options Database Keys:
+  -eps_stop_basic     - Sets the default stopping test
.  -eps_stop_threshold - Sets the threshold stopping test
-  -eps_stop_user      - Selects the user-defined stopping test

   Note:
   The parameter 'stop' can have one of these values
+     EPS_STOP_BASIC     - default stopping test
.     EPS_STOP_THRESHOLD - threshold stopping test)
-     EPS_STOP_USER      - function set by EPSSetStoppingTestFunction()

   Level: advanced

.seealso: EPSGetStoppingTest(), EPSSetStoppingTestFunction(), EPSSetConvergenceTest(), EPSStop
@*/
PetscErrorCode EPSSetStoppingTest(EPS eps,EPSStop stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,stop,2);
  switch (stop) {
    case EPS_STOP_BASIC: eps->stopping = EPSStoppingBasic; break;
    case EPS_STOP_THRESHOLD: eps->stopping = EPSStoppingThreshold; break;
    case EPS_STOP_USER:
      PetscCheck(eps->stoppinguser,PetscObjectComm((PetscObject)eps),PETSC_ERR_ORDER,"Must call EPSSetStoppingTestFunction() first");
      eps->stopping = eps->stoppinguser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'stop' value");
  }
  eps->stop = stop;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetStoppingTest - Gets the method used to decide the termination of the outer
   loop of the eigensolver.

   Not Collective

   Input Parameters:
.  eps   - eigensolver context obtained from EPSCreate()

   Output Parameters:
.  stop  - the type of stopping test

   Level: advanced

.seealso: EPSSetStoppingTest(), EPSStop
@*/
PetscErrorCode EPSGetStoppingTest(EPS eps,EPSStop *stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(stop,2);
  *stop = eps->stop;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetProblemType - Specifies the type of the eigenvalue problem.

   Logically Collective

   Input Parameters:
+  eps      - the eigensolver context
-  type     - a known type of eigenvalue problem

   Options Database Keys:
+  -eps_hermitian - Hermitian eigenvalue problem
.  -eps_gen_hermitian - generalized Hermitian eigenvalue problem
.  -eps_non_hermitian - non-Hermitian eigenvalue problem
.  -eps_gen_non_hermitian - generalized non-Hermitian eigenvalue problem
.  -eps_pos_gen_non_hermitian - generalized non-Hermitian eigenvalue problem
   with positive semi-definite B
.  -eps_gen_indefinite - generalized Hermitian-indefinite eigenvalue problem
-  -eps_bse - structured Bethe-Salpeter eigenvalue problem

   Notes:
   This function must be used to instruct SLEPc to exploit symmetry or other
   kind of structure. If no
   problem type is specified, by default a non-Hermitian problem is assumed
   (either standard or generalized). If the user knows that the problem is
   Hermitian (i.e. A=A^H) or generalized Hermitian (i.e. A=A^H, B=B^H, and
   B positive definite) then it is recommended to set the problem type so
   that eigensolver can exploit these properties.

   If the user does not call this function, the solver will use a reasonable
   guess.

   For structured problem types such as EPS_BSE, the matrices passed in via
   EPSSetOperators() must have been created with the corresponding helper
   function, i.e., MatCreateBSE().

   Level: intermediate

.seealso: EPSSetOperators(), EPSSetType(), EPSGetProblemType(), EPSProblemType
@*/
PetscErrorCode EPSSetProblemType(EPS eps,EPSProblemType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,type,2);
  if (type == eps->problem_type) PetscFunctionReturn(PETSC_SUCCESS);
  switch (type) {
    case EPS_HEP:
      eps->isgeneralized = PETSC_FALSE;
      eps->ishermitian = PETSC_TRUE;
      eps->ispositive = PETSC_FALSE;
      eps->isstructured = PETSC_FALSE;
      break;
    case EPS_NHEP:
      eps->isgeneralized = PETSC_FALSE;
      eps->ishermitian = PETSC_FALSE;
      eps->ispositive = PETSC_FALSE;
      eps->isstructured = PETSC_FALSE;
      break;
    case EPS_GHEP:
      eps->isgeneralized = PETSC_TRUE;
      eps->ishermitian = PETSC_TRUE;
      eps->ispositive = PETSC_TRUE;
      eps->isstructured = PETSC_FALSE;
      break;
    case EPS_GNHEP:
      eps->isgeneralized = PETSC_TRUE;
      eps->ishermitian = PETSC_FALSE;
      eps->ispositive = PETSC_FALSE;
      eps->isstructured = PETSC_FALSE;
      break;
    case EPS_PGNHEP:
      eps->isgeneralized = PETSC_TRUE;
      eps->ishermitian = PETSC_FALSE;
      eps->ispositive = PETSC_TRUE;
      eps->isstructured = PETSC_FALSE;
      break;
    case EPS_GHIEP:
      eps->isgeneralized = PETSC_TRUE;
      eps->ishermitian = PETSC_TRUE;
      eps->ispositive = PETSC_FALSE;
      eps->isstructured = PETSC_FALSE;
      break;
    case EPS_BSE:
      eps->isgeneralized = PETSC_FALSE;
      eps->ishermitian = PETSC_FALSE;
      eps->ispositive = PETSC_FALSE;
      eps->isstructured = PETSC_TRUE;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Unknown eigenvalue problem type");
  }
  eps->problem_type = type;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetProblemType - Gets the problem type from the EPS object.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  type - the problem type

   Level: intermediate

.seealso: EPSSetProblemType(), EPSProblemType
@*/
PetscErrorCode EPSGetProblemType(EPS eps,EPSProblemType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = eps->problem_type;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetExtraction - Specifies the type of extraction technique to be employed
   by the eigensolver.

   Logically Collective

   Input Parameters:
+  eps  - the eigensolver context
-  extr - a known type of extraction

   Options Database Keys:
+  -eps_ritz - Rayleigh-Ritz extraction
.  -eps_harmonic - harmonic Ritz extraction
.  -eps_harmonic_relative - harmonic Ritz extraction relative to the eigenvalue
.  -eps_harmonic_right - harmonic Ritz extraction for rightmost eigenvalues
.  -eps_harmonic_largest - harmonic Ritz extraction for largest magnitude
   (without target)
.  -eps_refined - refined Ritz extraction
-  -eps_refined_harmonic - refined harmonic Ritz extraction

   Notes:
   Not all eigensolvers support all types of extraction. See the SLEPc
   Users Manual for details.

   By default, a standard Rayleigh-Ritz extraction is used. Other extractions
   may be useful when computing interior eigenvalues.

   Harmonic-type extractions are used in combination with a 'target'.

   Level: advanced

.seealso: EPSSetTarget(), EPSGetExtraction(), EPSExtraction
@*/
PetscErrorCode EPSSetExtraction(EPS eps,EPSExtraction extr)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,extr,2);
  if (eps->extraction != extr) {
    eps->state      = EPS_STATE_INITIAL;
    eps->extraction = extr;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetExtraction - Gets the extraction type used by the EPS object.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  extr - name of extraction type

   Level: advanced

.seealso: EPSSetExtraction(), EPSExtraction
@*/
PetscErrorCode EPSGetExtraction(EPS eps,EPSExtraction *extr)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(extr,2);
  *extr = eps->extraction;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetBalance - Specifies the balancing technique to be employed by the
   eigensolver, and some parameters associated to it.

   Logically Collective

   Input Parameters:
+  eps    - the eigensolver context
.  bal    - the balancing method, one of EPS_BALANCE_NONE, EPS_BALANCE_ONESIDE,
            EPS_BALANCE_TWOSIDE, or EPS_BALANCE_USER
.  its    - number of iterations of the balancing algorithm
-  cutoff - cutoff value

   Options Database Keys:
+  -eps_balance <method> - the balancing method, where <method> is one of
                           'none', 'oneside', 'twoside', or 'user'
.  -eps_balance_its <its> - number of iterations
-  -eps_balance_cutoff <cutoff> - cutoff value

   Notes:
   When balancing is enabled, the solver works implicitly with matrix DAD^-1,
   where D is an appropriate diagonal matrix. This improves the accuracy of
   the computed results in some cases. See the SLEPc Users Manual for details.

   Balancing makes sense only for non-Hermitian problems when the required
   precision is high (i.e. a small tolerance such as 1e-15).

   By default, balancing is disabled. The two-sided method is much more
   effective than the one-sided counterpart, but it requires the system
   matrices to have the MatMultTranspose operation defined.

   The parameter 'its' is the number of iterations performed by the method. The
   cutoff value is used only in the two-side variant. Use PETSC_DETERMINE to assign
   a reasonably good value, or PETSC_CURRENT to leave the value unchanged.

   User-defined balancing is allowed provided that the corresponding matrix
   is set via STSetBalanceMatrix.

   Level: intermediate

.seealso: EPSGetBalance(), EPSBalance, STSetBalanceMatrix()
@*/
PetscErrorCode EPSSetBalance(EPS eps,EPSBalance bal,PetscInt its,PetscReal cutoff)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,bal,2);
  PetscValidLogicalCollectiveInt(eps,its,3);
  PetscValidLogicalCollectiveReal(eps,cutoff,4);
  switch (bal) {
    case EPS_BALANCE_NONE:
    case EPS_BALANCE_ONESIDE:
    case EPS_BALANCE_TWOSIDE:
    case EPS_BALANCE_USER:
      if (eps->balance != bal) {
        eps->state = EPS_STATE_INITIAL;
        eps->balance = bal;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid value of argument 'bal'");
  }
  if (its==PETSC_DETERMINE) eps->balance_its = 5;
  else if (its!=PETSC_CURRENT) {
    PetscCheck(its>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of its. Must be > 0");
    eps->balance_its = its;
  }
  if (cutoff==(PetscReal)PETSC_DETERMINE) eps->balance_cutoff = 1e-8;
  else if (cutoff!=(PetscReal)PETSC_CURRENT) {
    PetscCheck(cutoff>0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of cutoff. Must be > 0");
    eps->balance_cutoff = cutoff;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetBalance - Gets the balancing type used by the EPS object, and the
   associated parameters.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  bal    - the balancing method
.  its    - number of iterations of the balancing algorithm
-  cutoff - cutoff value

   Level: intermediate

   Note:
   The user can specify NULL for any parameter that is not needed.

.seealso: EPSSetBalance(), EPSBalance
@*/
PetscErrorCode EPSGetBalance(EPS eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (bal)    *bal = eps->balance;
  if (its)    *its = eps->balance_its;
  if (cutoff) *cutoff = eps->balance_cutoff;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetTwoSided - Sets the solver to use a two-sided variant so that left
   eigenvectors are also computed.

   Logically Collective

   Input Parameters:
+  eps      - the eigensolver context
-  twosided - whether the two-sided variant is to be used or not

   Options Database Keys:
.  -eps_two_sided <boolean> - Sets/resets the twosided flag

   Notes:
   If the user sets twosided=PETSC_TRUE then the solver uses a variant of
   the algorithm that computes both right and left eigenvectors. This is
   usually much more costly. This option is not available in all solvers.

   When using two-sided solvers, the problem matrices must have both the
   MatMult and MatMultTranspose operations defined.

   Level: advanced

.seealso: EPSGetTwoSided(), EPSGetLeftEigenvector()
@*/
PetscErrorCode EPSSetTwoSided(EPS eps,PetscBool twosided)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,twosided,2);
  if (twosided!=eps->twosided) {
    eps->twosided = twosided;
    eps->state    = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetTwoSided - Returns the flag indicating whether a two-sided variant
   of the algorithm is being used or not.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  twosided - the returned flag

   Level: advanced

.seealso: EPSSetTwoSided()
@*/
PetscErrorCode EPSGetTwoSided(EPS eps,PetscBool *twosided)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(twosided,2);
  *twosided = eps->twosided;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetTrueResidual - Specifies if the solver must compute the true residual
   explicitly or not.

   Logically Collective

   Input Parameters:
+  eps     - the eigensolver context
-  trueres - whether true residuals are required or not

   Options Database Keys:
.  -eps_true_residual <boolean> - Sets/resets the boolean flag 'trueres'

   Notes:
   If the user sets trueres=PETSC_TRUE then the solver explicitly computes
   the true residual for each eigenpair approximation, and uses it for
   convergence testing. Computing the residual is usually an expensive
   operation. Some solvers (e.g., Krylov solvers) can avoid this computation
   by using a cheap estimate of the residual norm, but this may sometimes
   give inaccurate results (especially if a spectral transform is being
   used). On the contrary, preconditioned eigensolvers (e.g., Davidson solvers)
   do rely on computing the true residual, so this option is irrelevant for them.

   Level: advanced

.seealso: EPSGetTrueResidual()
@*/
PetscErrorCode EPSSetTrueResidual(EPS eps,PetscBool trueres)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,trueres,2);
  eps->trueres = trueres;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetTrueResidual - Returns the flag indicating whether true
   residuals must be computed explicitly or not.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  trueres - the returned flag

   Level: advanced

.seealso: EPSSetTrueResidual()
@*/
PetscErrorCode EPSGetTrueResidual(EPS eps,PetscBool *trueres)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(trueres,2);
  *trueres = eps->trueres;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetTrackAll - Specifies if the solver must compute the residual norm of all
   approximate eigenpairs or not.

   Logically Collective

   Input Parameters:
+  eps      - the eigensolver context
-  trackall - whether to compute all residuals or not

   Notes:
   If the user sets trackall=PETSC_TRUE then the solver computes (or estimates)
   the residual norm for each eigenpair approximation. Computing the residual is
   usually an expensive operation and solvers commonly compute only the residual
   associated to the first unconverged eigenpair.

   The option '-eps_monitor_all' automatically activates this option.

   Level: developer

.seealso: EPSGetTrackAll()
@*/
PetscErrorCode EPSSetTrackAll(EPS eps,PetscBool trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,trackall,2);
  eps->trackall = trackall;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetTrackAll - Returns the flag indicating whether all residual norms must
   be computed or not.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  trackall - the returned flag

   Level: developer

.seealso: EPSSetTrackAll()
@*/
PetscErrorCode EPSGetTrackAll(EPS eps,PetscBool *trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(trackall,2);
  *trackall = eps->trackall;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetPurify - Deactivate eigenvector purification (which is activated by default).

   Logically Collective

   Input Parameters:
+  eps    - the eigensolver context
-  purify - whether purification is required or not

   Options Database Keys:
.  -eps_purify <boolean> - Sets/resets the boolean flag 'purify'

   Notes:
   By default, eigenvectors of generalized symmetric eigenproblems are purified
   in order to purge directions in the nullspace of matrix B. If the user knows
   that B is non-singular, then purification can be safely deactivated and some
   computational cost is avoided (this is particularly important in interval computations).

   Level: intermediate

.seealso: EPSGetPurify(), EPSSetInterval()
@*/
PetscErrorCode EPSSetPurify(EPS eps,PetscBool purify)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,purify,2);
  if (purify!=eps->purify) {
    eps->purify = purify;
    eps->state  = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetPurify - Returns the flag indicating whether purification is activated
   or not.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  purify - the returned flag

   Level: intermediate

.seealso: EPSSetPurify()
@*/
PetscErrorCode EPSGetPurify(EPS eps,PetscBool *purify)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(purify,2);
  *purify = eps->purify;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetOptionsPrefix - Sets the prefix used for searching for all
   EPS options in the database.

   Logically Collective

   Input Parameters:
+  eps - the eigensolver context
-  prefix - the prefix string to prepend to all EPS option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different EPS contexts, one could call
.vb
      EPSSetOptionsPrefix(eps1,"eig1_")
      EPSSetOptionsPrefix(eps2,"eig2_")
.ve

   Level: advanced

.seealso: EPSAppendOptionsPrefix(), EPSGetOptionsPrefix()
@*/
PetscErrorCode EPSSetOptionsPrefix(EPS eps,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!eps->st) PetscCall(EPSGetST(eps,&eps->st));
  PetscCall(STSetOptionsPrefix(eps->st,prefix));
  if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
  PetscCall(BVSetOptionsPrefix(eps->V,prefix));
  if (!eps->ds) PetscCall(EPSGetDS(eps,&eps->ds));
  PetscCall(DSSetOptionsPrefix(eps->ds,prefix));
  if (!eps->rg) PetscCall(EPSGetRG(eps,&eps->rg));
  PetscCall(RGSetOptionsPrefix(eps->rg,prefix));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)eps,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSAppendOptionsPrefix - Appends to the prefix used for searching for all
   EPS options in the database.

   Logically Collective

   Input Parameters:
+  eps - the eigensolver context
-  prefix - the prefix string to prepend to all EPS option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: EPSSetOptionsPrefix(), EPSGetOptionsPrefix()
@*/
PetscErrorCode EPSAppendOptionsPrefix(EPS eps,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!eps->st) PetscCall(EPSGetST(eps,&eps->st));
  PetscCall(STAppendOptionsPrefix(eps->st,prefix));
  if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
  PetscCall(BVAppendOptionsPrefix(eps->V,prefix));
  if (!eps->ds) PetscCall(EPSGetDS(eps,&eps->ds));
  PetscCall(DSAppendOptionsPrefix(eps->ds,prefix));
  if (!eps->rg) PetscCall(EPSGetRG(eps,&eps->rg));
  PetscCall(RGAppendOptionsPrefix(eps->rg,prefix));
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)eps,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetOptionsPrefix - Gets the prefix used for searching for all
   EPS options in the database.

   Not Collective

   Input Parameters:
.  eps - the eigensolver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: EPSSetOptionsPrefix(), EPSAppendOptionsPrefix()
@*/
PetscErrorCode EPSGetOptionsPrefix(EPS eps,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)eps,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}
