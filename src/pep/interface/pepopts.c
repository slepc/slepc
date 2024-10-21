/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   PEP routines related to options that can be set via the command-line
   or procedurally
*/

#include <slepc/private/pepimpl.h>       /*I "slepcpep.h" I*/
#include <petscdraw.h>

/*@C
   PEPMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective

   Input Parameters:
+  pep      - the polynomial eigensolver context
.  opt      - the command line option for this monitor
.  name     - the monitor type one is seeking
.  ctx      - an optional user context for the monitor, or NULL
-  trackall - whether this monitor tracks all eigenvalues or not

   Level: developer

.seealso: PEPMonitorSet(), PEPSetTrackAll()
@*/
PetscErrorCode PEPMonitorSetFromOptions(PEP pep,const char opt[],const char name[],void *ctx,PetscBool trackall)
{
  PetscErrorCode       (*mfunc)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
  PetscErrorCode       (*cfunc)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
  PetscErrorCode       (*dfunc)(PetscViewerAndFormat**);
  PetscViewerAndFormat *vf;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerType      vtype;
  char                 key[PETSC_MAX_PATH_LEN];
  PetscBool            flg;

  PetscFunctionBegin;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,opt,&viewer,&format,&flg));
  if (!flg) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetType(viewer,&vtype));
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListFind(PEPMonitorList,key,&mfunc));
  PetscCheck(mfunc,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Specified viewer and format not supported");
  PetscCall(PetscFunctionListFind(PEPMonitorCreateList,key,&cfunc));
  PetscCall(PetscFunctionListFind(PEPMonitorDestroyList,key,&dfunc));
  if (!cfunc) cfunc = PetscViewerAndFormatCreate_Internal;
  if (!dfunc) dfunc = PetscViewerAndFormatDestroy;

  PetscCall((*cfunc)(viewer,format,ctx,&vf));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(PEPMonitorSet(pep,mfunc,vf,(PetscCtxDestroyFn*)dfunc));
  if (trackall) PetscCall(PEPSetTrackAll(pep,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetFromOptions - Sets PEP options from the options database.
   This routine must be called before PEPSetUp() if the user is to be
   allowed to set the solver type.

   Collective

   Input Parameters:
.  pep - the polynomial eigensolver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner

.seealso: PEPSetOptionsPrefix()
@*/
PetscErrorCode PEPSetFromOptions(PEP pep)
{
  char            type[256];
  PetscBool       set,flg,flg1,flg2,flg3,flg4,flg5;
  PetscReal       r,t,array[2]={0,0};
  PetscScalar     s;
  PetscInt        i,j,k;
  PEPScale        scale;
  PEPRefine       refine;
  PEPRefineScheme scheme;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscCall(PEPRegisterAll());
  PetscObjectOptionsBegin((PetscObject)pep);
    PetscCall(PetscOptionsFList("-pep_type","Polynomial eigensolver method","PEPSetType",PEPList,(char*)(((PetscObject)pep)->type_name?((PetscObject)pep)->type_name:PEPTOAR),type,sizeof(type),&flg));
    if (flg) PetscCall(PEPSetType(pep,type));
    else if (!((PetscObject)pep)->type_name) PetscCall(PEPSetType(pep,PEPTOAR));

    PetscCall(PetscOptionsBoolGroupBegin("-pep_general","General polynomial eigenvalue problem","PEPSetProblemType",&flg));
    if (flg) PetscCall(PEPSetProblemType(pep,PEP_GENERAL));
    PetscCall(PetscOptionsBoolGroup("-pep_hermitian","Hermitian polynomial eigenvalue problem","PEPSetProblemType",&flg));
    if (flg) PetscCall(PEPSetProblemType(pep,PEP_HERMITIAN));
    PetscCall(PetscOptionsBoolGroup("-pep_hyperbolic","Hyperbolic polynomial eigenvalue problem","PEPSetProblemType",&flg));
    if (flg) PetscCall(PEPSetProblemType(pep,PEP_HYPERBOLIC));
    PetscCall(PetscOptionsBoolGroupEnd("-pep_gyroscopic","Gyroscopic polynomial eigenvalue problem","PEPSetProblemType",&flg));
    if (flg) PetscCall(PEPSetProblemType(pep,PEP_GYROSCOPIC));

    scale = pep->scale;
    PetscCall(PetscOptionsEnum("-pep_scale","Scaling strategy","PEPSetScale",PEPScaleTypes,(PetscEnum)scale,(PetscEnum*)&scale,&flg1));
    r = pep->sfactor;
    PetscCall(PetscOptionsReal("-pep_scale_factor","Scale factor","PEPSetScale",pep->sfactor,&r,&flg2));
    if (!flg2 && r==1.0) r = PETSC_DETERMINE;
    j = pep->sits;
    PetscCall(PetscOptionsInt("-pep_scale_its","Number of iterations in diagonal scaling","PEPSetScale",pep->sits,&j,&flg3));
    t = pep->slambda;
    PetscCall(PetscOptionsReal("-pep_scale_lambda","Estimate of eigenvalue (modulus) for diagonal scaling","PEPSetScale",pep->slambda,&t,&flg4));
    if (flg1 || flg2 || flg3 || flg4) PetscCall(PEPSetScale(pep,scale,r,NULL,NULL,j,t));

    PetscCall(PetscOptionsEnum("-pep_extract","Extraction method","PEPSetExtract",PEPExtractTypes,(PetscEnum)pep->extract,(PetscEnum*)&pep->extract,NULL));

    refine = pep->refine;
    PetscCall(PetscOptionsEnum("-pep_refine","Iterative refinement method","PEPSetRefine",PEPRefineTypes,(PetscEnum)refine,(PetscEnum*)&refine,&flg1));
    i = pep->npart;
    PetscCall(PetscOptionsInt("-pep_refine_partitions","Number of partitions of the communicator for iterative refinement","PEPSetRefine",pep->npart,&i,&flg2));
    r = pep->rtol;
    PetscCall(PetscOptionsReal("-pep_refine_tol","Tolerance for iterative refinement","PEPSetRefine",pep->rtol==(PetscReal)PETSC_DETERMINE?SLEPC_DEFAULT_TOL/1000:pep->rtol,&r,&flg3));
    j = pep->rits;
    PetscCall(PetscOptionsInt("-pep_refine_its","Maximum number of iterations for iterative refinement","PEPSetRefine",pep->rits,&j,&flg4));
    scheme = pep->scheme;
    PetscCall(PetscOptionsEnum("-pep_refine_scheme","Scheme used for linear systems within iterative refinement","PEPSetRefine",PEPRefineSchemes,(PetscEnum)scheme,(PetscEnum*)&scheme,&flg5));
    if (flg1 || flg2 || flg3 || flg4 || flg5) PetscCall(PEPSetRefine(pep,refine,i,r,j,scheme));

    i = pep->max_it;
    PetscCall(PetscOptionsInt("-pep_max_it","Maximum number of iterations","PEPSetTolerances",pep->max_it,&i,&flg1));
    r = pep->tol;
    PetscCall(PetscOptionsReal("-pep_tol","Tolerance","PEPSetTolerances",SlepcDefaultTol(pep->tol),&r,&flg2));
    if (flg1 || flg2) PetscCall(PEPSetTolerances(pep,r,i));

    PetscCall(PetscOptionsBoolGroupBegin("-pep_conv_rel","Relative error convergence test","PEPSetConvergenceTest",&flg));
    if (flg) PetscCall(PEPSetConvergenceTest(pep,PEP_CONV_REL));
    PetscCall(PetscOptionsBoolGroup("-pep_conv_norm","Convergence test relative to the matrix norms","PEPSetConvergenceTest",&flg));
    if (flg) PetscCall(PEPSetConvergenceTest(pep,PEP_CONV_NORM));
    PetscCall(PetscOptionsBoolGroup("-pep_conv_abs","Absolute error convergence test","PEPSetConvergenceTest",&flg));
    if (flg) PetscCall(PEPSetConvergenceTest(pep,PEP_CONV_ABS));
    PetscCall(PetscOptionsBoolGroupEnd("-pep_conv_user","User-defined convergence test","PEPSetConvergenceTest",&flg));
    if (flg) PetscCall(PEPSetConvergenceTest(pep,PEP_CONV_USER));

    PetscCall(PetscOptionsBoolGroupBegin("-pep_stop_basic","Stop iteration if all eigenvalues converged or max_it reached","PEPSetStoppingTest",&flg));
    if (flg) PetscCall(PEPSetStoppingTest(pep,PEP_STOP_BASIC));
    PetscCall(PetscOptionsBoolGroupEnd("-pep_stop_user","User-defined stopping test","PEPSetStoppingTest",&flg));
    if (flg) PetscCall(PEPSetStoppingTest(pep,PEP_STOP_USER));

    i = pep->nev;
    PetscCall(PetscOptionsInt("-pep_nev","Number of eigenvalues to compute","PEPSetDimensions",pep->nev,&i,&flg1));
    j = pep->ncv;
    PetscCall(PetscOptionsInt("-pep_ncv","Number of basis vectors","PEPSetDimensions",pep->ncv,&j,&flg2));
    k = pep->mpd;
    PetscCall(PetscOptionsInt("-pep_mpd","Maximum dimension of projected problem","PEPSetDimensions",pep->mpd,&k,&flg3));
    if (flg1 || flg2 || flg3) PetscCall(PEPSetDimensions(pep,i,j,k));

    PetscCall(PetscOptionsEnum("-pep_basis","Polynomial basis","PEPSetBasis",PEPBasisTypes,(PetscEnum)pep->basis,(PetscEnum*)&pep->basis,NULL));

    PetscCall(PetscOptionsBoolGroupBegin("-pep_largest_magnitude","Compute largest eigenvalues in magnitude","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_LARGEST_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-pep_smallest_magnitude","Compute smallest eigenvalues in magnitude","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_SMALLEST_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-pep_largest_real","Compute eigenvalues with largest real parts","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_LARGEST_REAL));
    PetscCall(PetscOptionsBoolGroup("-pep_smallest_real","Compute eigenvalues with smallest real parts","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_SMALLEST_REAL));
    PetscCall(PetscOptionsBoolGroup("-pep_largest_imaginary","Compute eigenvalues with largest imaginary parts","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_LARGEST_IMAGINARY));
    PetscCall(PetscOptionsBoolGroup("-pep_smallest_imaginary","Compute eigenvalues with smallest imaginary parts","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_SMALLEST_IMAGINARY));
    PetscCall(PetscOptionsBoolGroup("-pep_target_magnitude","Compute eigenvalues closest to target","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_TARGET_MAGNITUDE));
    PetscCall(PetscOptionsBoolGroup("-pep_target_real","Compute eigenvalues with real parts closest to target","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_TARGET_REAL));
    PetscCall(PetscOptionsBoolGroup("-pep_target_imaginary","Compute eigenvalues with imaginary parts closest to target","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_TARGET_IMAGINARY));
    PetscCall(PetscOptionsBoolGroupEnd("-pep_all","Compute all eigenvalues in an interval or a region","PEPSetWhichEigenpairs",&flg));
    if (flg) PetscCall(PEPSetWhichEigenpairs(pep,PEP_ALL));

    PetscCall(PetscOptionsScalar("-pep_target","Value of the target","PEPSetTarget",pep->target,&s,&flg));
    if (flg) {
      if (pep->which!=PEP_TARGET_REAL && pep->which!=PEP_TARGET_IMAGINARY) PetscCall(PEPSetWhichEigenpairs(pep,PEP_TARGET_MAGNITUDE));
      PetscCall(PEPSetTarget(pep,s));
    }

    k = 2;
    PetscCall(PetscOptionsRealArray("-pep_interval","Computational interval (two real values separated with a comma without spaces)","PEPSetInterval",array,&k,&flg));
    if (flg) {
      PetscCheck(k>1,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_SIZ,"Must pass two values in -pep_interval (comma-separated without spaces)");
      PetscCall(PEPSetWhichEigenpairs(pep,PEP_ALL));
      PetscCall(PEPSetInterval(pep,array[0],array[1]));
    }

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to PEPSetFromOptions()
    */
    PetscCall(PetscOptionsBool("-pep_monitor_cancel","Remove any hardwired monitor routines","PEPMonitorCancel",PETSC_FALSE,&flg,&set));
    if (set && flg) PetscCall(PEPMonitorCancel(pep));
    PetscCall(PEPMonitorSetFromOptions(pep,"-pep_monitor","first_approximation",NULL,PETSC_FALSE));
    PetscCall(PEPMonitorSetFromOptions(pep,"-pep_monitor_all","all_approximations",NULL,PETSC_TRUE));
    PetscCall(PEPMonitorSetFromOptions(pep,"-pep_monitor_conv","convergence_history",NULL,PETSC_FALSE));

    /* -----------------------------------------------------------------------*/
    PetscCall(PetscOptionsName("-pep_view","Print detailed information on solver used","PEPView",&set));
    PetscCall(PetscOptionsName("-pep_view_vectors","View computed eigenvectors","PEPVectorsView",&set));
    PetscCall(PetscOptionsName("-pep_view_values","View computed eigenvalues","PEPValuesView",&set));
    PetscCall(PetscOptionsName("-pep_converged_reason","Print reason for convergence, and number of iterations","PEPConvergedReasonView",&set));
    PetscCall(PetscOptionsName("-pep_error_absolute","Print absolute errors of each eigenpair","PEPErrorView",&set));
    PetscCall(PetscOptionsName("-pep_error_relative","Print relative errors of each eigenpair","PEPErrorView",&set));
    PetscCall(PetscOptionsName("-pep_error_backward","Print backward errors of each eigenpair","PEPErrorView",&set));

    PetscTryTypeMethod(pep,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)pep,PetscOptionsObject));
  PetscOptionsEnd();

  if (!pep->V) PetscCall(PEPGetBV(pep,&pep->V));
  PetscCall(BVSetFromOptions(pep->V));
  if (!pep->rg) PetscCall(PEPGetRG(pep,&pep->rg));
  PetscCall(RGSetFromOptions(pep->rg));
  if (!pep->ds) PetscCall(PEPGetDS(pep,&pep->ds));
  PetscCall(PEPSetDSType(pep));
  PetscCall(DSSetFromOptions(pep->ds));
  if (!pep->st) PetscCall(PEPGetST(pep,&pep->st));
  PetscCall(PEPSetDefaultST(pep));
  PetscCall(STSetFromOptions(pep->st));
  if (!pep->refineksp) PetscCall(PEPRefineGetKSP(pep,&pep->refineksp));
  PetscCall(KSPSetFromOptions(pep->refineksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetTolerances - Gets the tolerance and maximum iteration count used
   by the PEP convergence tests.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: PEPSetTolerances()
@*/
PetscErrorCode PEPGetTolerances(PEP pep,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (tol)    *tol    = pep->tol;
  if (maxits) *maxits = pep->max_it;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetTolerances - Sets the tolerance and maximum iteration count used
   by the PEP convergence tests.

   Logically Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -pep_tol <tol> - Sets the convergence tolerance
-  -pep_max_it <maxits> - Sets the maximum number of iterations allowed

   Notes:
   Use PETSC_CURRENT to retain the current value of any of the parameters.
   Use PETSC_DETERMINE for either argument to assign a default value computed
   internally (may be different in each solver).
   For maxits use PETSC_UMLIMITED to indicate there is no upper bound on this value.

   Level: intermediate

.seealso: PEPGetTolerances()
@*/
PetscErrorCode PEPSetTolerances(PEP pep,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,tol,2);
  PetscValidLogicalCollectiveInt(pep,maxits,3);
  if (tol == (PetscReal)PETSC_DETERMINE) {
    pep->tol   = PETSC_DETERMINE;
    pep->state = PEP_STATE_INITIAL;
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    pep->tol = tol;
  }
  if (maxits == PETSC_DETERMINE) {
    pep->max_it = PETSC_DETERMINE;
    pep->state  = PEP_STATE_INITIAL;
  } else if (maxits == PETSC_UNLIMITED) {
    pep->max_it = PETSC_INT_MAX;
  } else if (maxits != PETSC_CURRENT) {
    PetscCheck(maxits>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    pep->max_it = maxits;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetDimensions - Gets the number of eigenvalues to compute
   and the dimension of the subspace.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: PEPSetDimensions()
@*/
PetscErrorCode PEPGetDimensions(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (nev) *nev = pep->nev;
  if (ncv) *ncv = pep->ncv;
  if (mpd) *mpd = pep->mpd;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetDimensions - Sets the number of eigenvalues to compute
   and the dimension of the subspace.

   Logically Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Keys:
+  -pep_nev <nev> - Sets the number of eigenvalues
.  -pep_ncv <ncv> - Sets the dimension of the subspace
-  -pep_mpd <mpd> - Sets the maximum projected dimension

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

   When computing all eigenvalues in an interval, see PEPSetInterval(), these
   parameters lose relevance, and tuning must be done with PEPSTOARSetDimensions().

   Level: intermediate

.seealso: PEPGetDimensions(), PEPSetInterval(), PEPSTOARSetDimensions()
@*/
PetscErrorCode PEPSetDimensions(PEP pep,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,nev,2);
  PetscValidLogicalCollectiveInt(pep,ncv,3);
  PetscValidLogicalCollectiveInt(pep,mpd,4);
  if (nev != PETSC_CURRENT) {
    PetscCheck(nev>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
    pep->nev = nev;
  }
  if (ncv == PETSC_DETERMINE) {
    pep->ncv = PETSC_DETERMINE;
  } else if (ncv != PETSC_CURRENT) {
    PetscCheck(ncv>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    pep->ncv = ncv;
  }
  if (mpd == PETSC_DETERMINE) {
    pep->mpd = PETSC_DETERMINE;
  } else if (mpd != PETSC_CURRENT) {
    PetscCheck(mpd>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    pep->mpd = mpd;
  }
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetWhichEigenpairs - Specifies which portion of the spectrum is
   to be sought.

   Logically Collective

   Input Parameters:
+  pep   - eigensolver context obtained from PEPCreate()
-  which - the portion of the spectrum to be sought

   Options Database Keys:
+   -pep_largest_magnitude - Sets largest eigenvalues in magnitude
.   -pep_smallest_magnitude - Sets smallest eigenvalues in magnitude
.   -pep_largest_real - Sets largest real parts
.   -pep_smallest_real - Sets smallest real parts
.   -pep_largest_imaginary - Sets largest imaginary parts
.   -pep_smallest_imaginary - Sets smallest imaginary parts
.   -pep_target_magnitude - Sets eigenvalues closest to target
.   -pep_target_real - Sets real parts closest to target
.   -pep_target_imaginary - Sets imaginary parts closest to target
-   -pep_all - Sets all eigenvalues in an interval or region

   Notes:
   The parameter 'which' can have one of these values

+     PEP_LARGEST_MAGNITUDE - largest eigenvalues in magnitude (default)
.     PEP_SMALLEST_MAGNITUDE - smallest eigenvalues in magnitude
.     PEP_LARGEST_REAL - largest real parts
.     PEP_SMALLEST_REAL - smallest real parts
.     PEP_LARGEST_IMAGINARY - largest imaginary parts
.     PEP_SMALLEST_IMAGINARY - smallest imaginary parts
.     PEP_TARGET_MAGNITUDE - eigenvalues closest to the target (in magnitude)
.     PEP_TARGET_REAL - eigenvalues with real part closest to target
.     PEP_TARGET_IMAGINARY - eigenvalues with imaginary part closest to target
.     PEP_ALL - all eigenvalues contained in a given interval or region
-     PEP_WHICH_USER - user defined ordering set with PEPSetEigenvalueComparison()

   Not all eigensolvers implemented in PEP account for all the possible values
   stated above. If SLEPc is compiled for real numbers PEP_LARGEST_IMAGINARY
   and PEP_SMALLEST_IMAGINARY use the absolute value of the imaginary part
   for eigenvalue selection.

   The target is a scalar value provided with PEPSetTarget().

   The criterion PEP_TARGET_IMAGINARY is available only in case PETSc and
   SLEPc have been built with complex scalars.

   PEP_ALL is intended for use in combination with an interval (see
   PEPSetInterval()), when all eigenvalues within the interval are requested,
   and also for computing all eigenvalues in a region with the CISS solver.
   In both cases, the number of eigenvalues is unknown, so the nev parameter
   has a different sense, see PEPSetDimensions().

   Level: intermediate

.seealso: PEPGetWhichEigenpairs(), PEPSetTarget(), PEPSetInterval(),
          PEPSetDimensions(), PEPSetEigenvalueComparison(), PEPWhich
@*/
PetscErrorCode PEPSetWhichEigenpairs(PEP pep,PEPWhich which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,which,2);
  switch (which) {
    case PEP_LARGEST_MAGNITUDE:
    case PEP_SMALLEST_MAGNITUDE:
    case PEP_LARGEST_REAL:
    case PEP_SMALLEST_REAL:
    case PEP_LARGEST_IMAGINARY:
    case PEP_SMALLEST_IMAGINARY:
    case PEP_TARGET_MAGNITUDE:
    case PEP_TARGET_REAL:
#if defined(PETSC_USE_COMPLEX)
    case PEP_TARGET_IMAGINARY:
#endif
    case PEP_ALL:
    case PEP_WHICH_USER:
      if (pep->which != which) {
        pep->state = PEP_STATE_INITIAL;
        pep->which = which;
      }
      break;
#if !defined(PETSC_USE_COMPLEX)
    case PEP_TARGET_IMAGINARY:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"PEP_TARGET_IMAGINARY can be used only with complex scalars");
#endif
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'which' value");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
    PEPGetWhichEigenpairs - Returns which portion of the spectrum is to be
    sought.

    Not Collective

    Input Parameter:
.   pep - eigensolver context obtained from PEPCreate()

    Output Parameter:
.   which - the portion of the spectrum to be sought

    Notes:
    See PEPSetWhichEigenpairs() for possible values of 'which'.

    Level: intermediate

.seealso: PEPSetWhichEigenpairs(), PEPWhich
@*/
PetscErrorCode PEPGetWhichEigenpairs(PEP pep,PEPWhich *which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(which,2);
  *which = pep->which;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPSetEigenvalueComparison - Specifies the eigenvalue comparison function
   when PEPSetWhichEigenpairs() is set to PEP_WHICH_USER.

   Logically Collective

   Input Parameters:
+  pep  - eigensolver context obtained from PEPCreate()
.  comp - a pointer to the comparison function
-  ctx  - a context pointer (the last parameter to the comparison function)

   Calling sequence of comp:
$  PetscErrorCode comp(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx)
+   ar     - real part of the 1st eigenvalue
.   ai     - imaginary part of the 1st eigenvalue
.   br     - real part of the 2nd eigenvalue
.   bi     - imaginary part of the 2nd eigenvalue
.   res    - result of comparison
-   ctx    - optional context, as set by PEPSetEigenvalueComparison()

   Note:
   The returning parameter 'res' can be
+  negative - if the 1st eigenvalue is preferred to the 2st one
.  zero     - if both eigenvalues are equally preferred
-  positive - if the 2st eigenvalue is preferred to the 1st one

   Level: advanced

.seealso: PEPSetWhichEigenpairs(), PEPWhich
@*/
PetscErrorCode PEPSetEigenvalueComparison(PEP pep,PetscErrorCode (*comp)(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx),void* ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  pep->sc->comparison    = comp;
  pep->sc->comparisonctx = ctx;
  pep->which             = PEP_WHICH_USER;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetProblemType - Specifies the type of the polynomial eigenvalue problem.

   Logically Collective

   Input Parameters:
+  pep  - the polynomial eigensolver context
-  type - a known type of polynomial eigenvalue problem

   Options Database Keys:
+  -pep_general - general problem with no particular structure
.  -pep_hermitian - problem whose coefficient matrices are Hermitian
.  -pep_hyperbolic - Hermitian problem that satisfies the definition of hyperbolic
-  -pep_gyroscopic - problem with Hamiltonian structure

   Notes:
   Allowed values for the problem type are general (PEP_GENERAL), Hermitian
   (PEP_HERMITIAN), hyperbolic (PEP_HYPERBOLIC), and gyroscopic (PEP_GYROSCOPIC).

   This function is used to instruct SLEPc to exploit certain structure in
   the polynomial eigenproblem. By default, no particular structure is assumed.

   If the problem matrices are Hermitian (symmetric in the real case) or
   Hermitian/skew-Hermitian then the solver can exploit this fact to perform
   less operations or provide better stability. Hyperbolic problems are a
   particular case of Hermitian problems, some solvers may treat them simply as
   Hermitian.

   Level: intermediate

.seealso: PEPSetOperators(), PEPSetType(), PEPGetProblemType(), PEPProblemType
@*/
PetscErrorCode PEPSetProblemType(PEP pep,PEPProblemType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,type,2);
  PetscCheck(type==PEP_GENERAL || type==PEP_HERMITIAN || type==PEP_HYPERBOLIC || type==PEP_GYROSCOPIC,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"Unknown eigenvalue problem type");
  if (type != pep->problem_type) {
    pep->problem_type = type;
    pep->state = PEP_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetProblemType - Gets the problem type from the PEP object.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameter:
.  type - the problem type

   Level: intermediate

.seealso: PEPSetProblemType(), PEPProblemType
@*/
PetscErrorCode PEPGetProblemType(PEP pep,PEPProblemType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = pep->problem_type;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetBasis - Specifies the type of polynomial basis used to describe the
   polynomial eigenvalue problem.

   Logically Collective

   Input Parameters:
+  pep   - the polynomial eigensolver context
-  basis - the type of polynomial basis

   Options Database Key:
.  -pep_basis <basis> - Select the basis type

   Notes:
   By default, the coefficient matrices passed via PEPSetOperators() are
   expressed in the monomial basis, i.e.
   P(lambda) = A_0 + lambda*A_1 + lambda^2*A_2 + ... + lambda^d*A_d.
   Other polynomial bases may have better numerical behaviour, but the user
   must then pass the coefficient matrices accordingly.

   Level: intermediate

.seealso: PEPSetOperators(), PEPGetBasis(), PEPBasis
@*/
PetscErrorCode PEPSetBasis(PEP pep,PEPBasis basis)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,basis,2);
  pep->basis = basis;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetBasis - Gets the type of polynomial basis from the PEP object.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameter:
.  basis - the polynomial basis

   Level: intermediate

.seealso: PEPSetBasis(), PEPBasis
@*/
PetscErrorCode PEPGetBasis(PEP pep,PEPBasis *basis)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(basis,2);
  *basis = pep->basis;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetTrackAll - Specifies if the solver must compute the residual of all
   approximate eigenpairs or not.

   Logically Collective

   Input Parameters:
+  pep      - the eigensolver context
-  trackall - whether compute all residuals or not

   Notes:
   If the user sets trackall=PETSC_TRUE then the solver explicitly computes
   the residual for each eigenpair approximation. Computing the residual is
   usually an expensive operation and solvers commonly compute the associated
   residual to the first unconverged eigenpair.

   The option '-pep_monitor_all' automatically activates this option.

   Level: developer

.seealso: PEPGetTrackAll()
@*/
PetscErrorCode PEPSetTrackAll(PEP pep,PetscBool trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,trackall,2);
  pep->trackall = trackall;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetTrackAll - Returns the flag indicating whether all residual norms must
   be computed or not.

   Not Collective

   Input Parameter:
.  pep - the eigensolver context

   Output Parameter:
.  trackall - the returned flag

   Level: developer

.seealso: PEPSetTrackAll()
@*/
PetscErrorCode PEPGetTrackAll(PEP pep,PetscBool *trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(trackall,2);
  *trackall = pep->trackall;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPSetConvergenceTestFunction - Sets a function to compute the error estimate
   used in the convergence test.

   Logically Collective

   Input Parameters:
+  pep     - eigensolver context obtained from PEPCreate()
.  conv    - convergence test function, see PEPConvergenceTestFn for the calling sequence
.  ctx     - context for private data for the convergence routine (may be NULL)
-  destroy - a routine for destroying the context (may be NULL), see PetscCtxDestroyFn for the calling sequence

   Note:
   If the error estimate returned by the convergence test function is less than
   the tolerance, then the eigenvalue is accepted as converged.

   Level: advanced

.seealso: PEPSetConvergenceTest(), PEPSetTolerances()
@*/
PetscErrorCode PEPSetConvergenceTestFunction(PEP pep,PEPConvergenceTestFn *conv,void* ctx,PetscCtxDestroyFn *destroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (pep->convergeddestroy) PetscCall((*pep->convergeddestroy)(&pep->convergedctx));
  pep->convergeduser    = conv;
  pep->convergeddestroy = destroy;
  pep->convergedctx     = ctx;
  if (conv == PEPConvergedRelative) pep->conv = PEP_CONV_REL;
  else if (conv == PEPConvergedNorm) pep->conv = PEP_CONV_NORM;
  else if (conv == PEPConvergedAbsolute) pep->conv = PEP_CONV_ABS;
  else {
    pep->conv      = PEP_CONV_USER;
    pep->converged = pep->convergeduser;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetConvergenceTest - Specifies how to compute the error estimate
   used in the convergence test.

   Logically Collective

   Input Parameters:
+  pep  - eigensolver context obtained from PEPCreate()
-  conv - the type of convergence test

   Options Database Keys:
+  -pep_conv_abs    - Sets the absolute convergence test
.  -pep_conv_rel    - Sets the convergence test relative to the eigenvalue
.  -pep_conv_norm   - Sets the convergence test relative to the matrix norms
-  -pep_conv_user   - Selects the user-defined convergence test

   Note:
   The parameter 'conv' can have one of these values
+     PEP_CONV_ABS    - absolute error ||r||
.     PEP_CONV_REL    - error relative to the eigenvalue l, ||r||/|l|
.     PEP_CONV_NORM   - error relative matrix norms, ||r||/sum_i(l^i*||A_i||)
-     PEP_CONV_USER   - function set by PEPSetConvergenceTestFunction()

   Level: intermediate

.seealso: PEPGetConvergenceTest(), PEPSetConvergenceTestFunction(), PEPSetStoppingTest(), PEPConv
@*/
PetscErrorCode PEPSetConvergenceTest(PEP pep,PEPConv conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,conv,2);
  switch (conv) {
    case PEP_CONV_ABS:  pep->converged = PEPConvergedAbsolute; break;
    case PEP_CONV_REL:  pep->converged = PEPConvergedRelative; break;
    case PEP_CONV_NORM: pep->converged = PEPConvergedNorm; break;
    case PEP_CONV_USER:
      PetscCheck(pep->convergeduser,PetscObjectComm((PetscObject)pep),PETSC_ERR_ORDER,"Must call PEPSetConvergenceTestFunction() first");
      pep->converged = pep->convergeduser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'conv' value");
  }
  pep->conv = conv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetConvergenceTest - Gets the method used to compute the error estimate
   used in the convergence test.

   Not Collective

   Input Parameters:
.  pep   - eigensolver context obtained from PEPCreate()

   Output Parameters:
.  conv  - the type of convergence test

   Level: intermediate

.seealso: PEPSetConvergenceTest(), PEPConv
@*/
PetscErrorCode PEPGetConvergenceTest(PEP pep,PEPConv *conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(conv,2);
  *conv = pep->conv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPSetStoppingTestFunction - Sets a function to decide when to stop the outer
   iteration of the eigensolver.

   Logically Collective

   Input Parameters:
+  pep     - eigensolver context obtained from PEPCreate()
.  stop    - stopping test function, see PEPStoppingTestFn for the calling sequence
.  ctx     - context for private data for the stopping routine (may be NULL)
-  destroy - a routine for destroying the context (may be NULL), see PetscCtxDestroyFn for the calling sequence

   Note:
   Normal usage is to first call the default routine PEPStoppingBasic() and then
   set reason to PEP_CONVERGED_USER if some user-defined conditions have been
   met. To let the eigensolver continue iterating, the result must be left as
   PEP_CONVERGED_ITERATING.

   Level: advanced

.seealso: PEPSetStoppingTest(), PEPStoppingBasic()
@*/
PetscErrorCode PEPSetStoppingTestFunction(PEP pep,PEPStoppingTestFn *stop,void* ctx,PetscCtxDestroyFn *destroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (pep->stoppingdestroy) PetscCall((*pep->stoppingdestroy)(&pep->stoppingctx));
  pep->stoppinguser    = stop;
  pep->stoppingdestroy = destroy;
  pep->stoppingctx     = ctx;
  if (stop == PEPStoppingBasic) pep->stop = PEP_STOP_BASIC;
  else {
    pep->stop     = PEP_STOP_USER;
    pep->stopping = pep->stoppinguser;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetStoppingTest - Specifies how to decide the termination of the outer
   loop of the eigensolver.

   Logically Collective

   Input Parameters:
+  pep  - eigensolver context obtained from PEPCreate()
-  stop - the type of stopping test

   Options Database Keys:
+  -pep_stop_basic - Sets the default stopping test
-  -pep_stop_user  - Selects the user-defined stopping test

   Note:
   The parameter 'stop' can have one of these values
+     PEP_STOP_BASIC - default stopping test
-     PEP_STOP_USER  - function set by PEPSetStoppingTestFunction()

   Level: advanced

.seealso: PEPGetStoppingTest(), PEPSetStoppingTestFunction(), PEPSetConvergenceTest(), PEPStop
@*/
PetscErrorCode PEPSetStoppingTest(PEP pep,PEPStop stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,stop,2);
  switch (stop) {
    case PEP_STOP_BASIC: pep->stopping = PEPStoppingBasic; break;
    case PEP_STOP_USER:
      PetscCheck(pep->stoppinguser,PetscObjectComm((PetscObject)pep),PETSC_ERR_ORDER,"Must call PEPSetStoppingTestFunction() first");
      pep->stopping = pep->stoppinguser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'stop' value");
  }
  pep->stop = stop;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetStoppingTest - Gets the method used to decide the termination of the outer
   loop of the eigensolver.

   Not Collective

   Input Parameters:
.  pep   - eigensolver context obtained from PEPCreate()

   Output Parameters:
.  stop  - the type of stopping test

   Level: advanced

.seealso: PEPSetStoppingTest(), PEPStop
@*/
PetscErrorCode PEPGetStoppingTest(PEP pep,PEPStop *stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(stop,2);
  *stop = pep->stop;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetScale - Specifies the scaling strategy to be used.

   Collective

   Input Parameters:
+  pep    - the eigensolver context
.  scale  - scaling strategy
.  alpha  - the scaling factor used in the scalar strategy
.  Dl     - the left diagonal matrix of the diagonal scaling algorithm
.  Dr     - the right diagonal matrix of the diagonal scaling algorithm
.  its    - number of iterations of the diagonal scaling algorithm
-  lambda - approximation to wanted eigenvalues (modulus)

   Options Database Keys:
+  -pep_scale <type> - scaling type, one of <none,scalar,diagonal,both>
.  -pep_scale_factor <alpha> - the scaling factor
.  -pep_scale_its <its> - number of iterations
-  -pep_scale_lambda <lambda> - approximation to eigenvalues

   Notes:
   There are two non-exclusive scaling strategies, scalar and diagonal.

   In the scalar strategy, scaling is applied to the eigenvalue, that is,
   mu = lambda/alpha is the new eigenvalue and all matrices are scaled
   accordingly. After solving the scaled problem, the original lambda is
   recovered. Parameter 'alpha' must be positive. Use PETSC_DETERMINE to let
   the solver compute a reasonable scaling factor, and PETSC_CURRENT to
   retain a previously set value.

   In the diagonal strategy, the solver works implicitly with matrix Dl*A*Dr,
   where Dl and Dr are appropriate diagonal matrices. This improves the accuracy
   of the computed results in some cases. The user may provide the Dr and Dl
   matrices represented as Vec objects storing diagonal elements. If not
   provided, these matrices are computed internally. This option requires
   that the polynomial coefficient matrices are of MATAIJ type.
   The parameter 'its' is the number of iterations performed by the method.
   Parameter 'lambda' must be positive. Use PETSC_DETERMINE or set lambda = 1.0
   if no information about eigenvalues is available. PETSC_CURRENT can also
   be used to leave its and lambda unchanged.

   Level: intermediate

.seealso: PEPGetScale()
@*/
PetscErrorCode PEPSetScale(PEP pep,PEPScale scale,PetscReal alpha,Vec Dl,Vec Dr,PetscInt its,PetscReal lambda)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,scale,2);
  pep->scale = scale;
  if (scale==PEP_SCALE_SCALAR || scale==PEP_SCALE_BOTH) {
    PetscValidLogicalCollectiveReal(pep,alpha,3);
    if (alpha == (PetscReal)PETSC_DETERMINE) {
      pep->sfactor = 0.0;
      pep->sfactor_set = PETSC_FALSE;
    } else if (alpha != (PetscReal)PETSC_CURRENT) {
      PetscCheck(alpha>0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of alpha. Must be > 0");
      pep->sfactor = alpha;
      pep->sfactor_set = PETSC_TRUE;
    }
  }
  if (scale==PEP_SCALE_DIAGONAL || scale==PEP_SCALE_BOTH) {
    if (Dl) {
      PetscValidHeaderSpecific(Dl,VEC_CLASSID,4);
      PetscCheckSameComm(pep,1,Dl,4);
      PetscCall(PetscObjectReference((PetscObject)Dl));
      PetscCall(VecDestroy(&pep->Dl));
      pep->Dl = Dl;
    }
    if (Dr) {
      PetscValidHeaderSpecific(Dr,VEC_CLASSID,5);
      PetscCheckSameComm(pep,1,Dr,5);
      PetscCall(PetscObjectReference((PetscObject)Dr));
      PetscCall(VecDestroy(&pep->Dr));
      pep->Dr = Dr;
    }
    PetscValidLogicalCollectiveInt(pep,its,6);
    PetscValidLogicalCollectiveReal(pep,lambda,7);
    if (its==PETSC_DETERMINE) pep->sits = 5;
    else if (its!=PETSC_CURRENT) {
      PetscCheck(its>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of its. Must be > 0");
      pep->sits = its;
    }
    if (lambda == (PetscReal)PETSC_DETERMINE) pep->slambda = 1.0;
    else if (lambda != (PetscReal)PETSC_CURRENT) {
      PetscCheck(lambda>0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of lambda. Must be > 0");
      pep->slambda = lambda;
    }
  }
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetScale - Gets the scaling strategy used by the PEP object, and the
   associated parameters.

   Not Collective

   Input Parameter:
.  pep - the eigensolver context

   Output Parameters:
+  scale  - scaling strategy
.  alpha  - the scaling factor used in the scalar strategy
.  Dl     - the left diagonal matrix of the diagonal scaling algorithm
.  Dr     - the right diagonal matrix of the diagonal scaling algorithm
.  its    - number of iterations of the diagonal scaling algorithm
-  lambda - approximation to wanted eigenvalues (modulus)

   Level: intermediate

   Note:
   The user can specify NULL for any parameter that is not needed.

   If Dl or Dr were not set by the user, then the ones computed internally are
   returned (or a null pointer if called before PEPSetUp).

.seealso: PEPSetScale(), PEPSetUp()
@*/
PetscErrorCode PEPGetScale(PEP pep,PEPScale *scale,PetscReal *alpha,Vec *Dl,Vec *Dr,PetscInt *its,PetscReal *lambda)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (scale)  *scale  = pep->scale;
  if (alpha)  *alpha  = pep->sfactor;
  if (Dl)     *Dl     = pep->Dl;
  if (Dr)     *Dr     = pep->Dr;
  if (its)    *its    = pep->sits;
  if (lambda) *lambda = pep->slambda;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetExtract - Specifies the extraction strategy to be used.

   Logically Collective

   Input Parameters:
+  pep     - the eigensolver context
-  extract - extraction strategy

   Options Database Keys:
.  -pep_extract <type> - extraction type, one of <none,norm,residual,structured>

   Level: intermediate

.seealso: PEPGetExtract()
@*/
PetscErrorCode PEPSetExtract(PEP pep,PEPExtract extract)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,extract,2);
  pep->extract = extract;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetExtract - Gets the extraction strategy used by the PEP object.

   Not Collective

   Input Parameter:
.  pep - the eigensolver context

   Output Parameter:
.  extract - extraction strategy

   Level: intermediate

.seealso: PEPSetExtract(), PEPExtract
@*/
PetscErrorCode PEPGetExtract(PEP pep,PEPExtract *extract)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(extract,2);
  *extract = pep->extract;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetRefine - Specifies the refinement type (and options) to be used
   after the solve.

   Logically Collective

   Input Parameters:
+  pep    - the polynomial eigensolver context
.  refine - refinement type
.  npart  - number of partitions of the communicator
.  tol    - the convergence tolerance
.  its    - maximum number of refinement iterations
-  scheme - which scheme to be used for solving the involved linear systems

   Options Database Keys:
+  -pep_refine <type> - refinement type, one of <none,simple,multiple>
.  -pep_refine_partitions <n> - the number of partitions
.  -pep_refine_tol <tol> - the tolerance
.  -pep_refine_its <its> - number of iterations
-  -pep_refine_scheme - to set the scheme for the linear solves

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
   below the tolerance (tol defaults to the PEP tol, but may be set to a
   different value). In contrast, the multiple method simply performs its
   refinement iterations (just one by default).

   The scheme argument is used to change the way in which linear systems are
   solved. Possible choices are explicit, mixed block elimination (MBE),
   and Schur complement.

   Use PETSC_CURRENT to retain the current value of npart, tol or its. Use
   PETSC_DETERMINE to assign a default value.

   Level: intermediate

.seealso: PEPGetRefine()
@*/
PetscErrorCode PEPSetRefine(PEP pep,PEPRefine refine,PetscInt npart,PetscReal tol,PetscInt its,PEPRefineScheme scheme)
{
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,refine,2);
  PetscValidLogicalCollectiveInt(pep,npart,3);
  PetscValidLogicalCollectiveReal(pep,tol,4);
  PetscValidLogicalCollectiveInt(pep,its,5);
  PetscValidLogicalCollectiveEnum(pep,scheme,6);
  pep->refine = refine;
  if (refine) {  /* process parameters only if not REFINE_NONE */
    if (npart!=pep->npart) {
      PetscCall(PetscSubcommDestroy(&pep->refinesubc));
      PetscCall(KSPDestroy(&pep->refineksp));
    }
    if (npart == PETSC_DETERMINE) {
      pep->npart = 1;
    } else if (npart != PETSC_CURRENT) {
      PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pep),&size));
      PetscCheck(npart>0 && npart<=size,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of npart");
      pep->npart = npart;
    }
    if (tol == (PetscReal)PETSC_DETERMINE) {
      pep->rtol = PETSC_DETERMINE;
    } else if (tol != (PetscReal)PETSC_CURRENT) {
      PetscCheck(tol>0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
      pep->rtol = tol;
    }
    if (its==PETSC_DETERMINE) {
      pep->rits = PETSC_DETERMINE;
    } else if (its != PETSC_CURRENT) {
      PetscCheck(its>=0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of its. Must be >= 0");
      pep->rits = its;
    }
    pep->scheme = scheme;
  }
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetRefine - Gets the refinement strategy used by the PEP object, and the
   associated parameters.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  refine - refinement type
.  npart  - number of partitions of the communicator
.  tol    - the convergence tolerance
.  its    - maximum number of refinement iterations
-  scheme - the scheme used for solving linear systems

   Level: intermediate

   Note:
   The user can specify NULL for any parameter that is not needed.

.seealso: PEPSetRefine()
@*/
PetscErrorCode PEPGetRefine(PEP pep,PEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,PEPRefineScheme *scheme)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (refine) *refine = pep->refine;
  if (npart)  *npart  = pep->npart;
  if (tol)    *tol    = pep->rtol;
  if (its)    *its    = pep->rits;
  if (scheme) *scheme = pep->scheme;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetOptionsPrefix - Sets the prefix used for searching for all
   PEP options in the database.

   Logically Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
-  prefix - the prefix string to prepend to all PEP option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different PEP contexts, one could call
.vb
      PEPSetOptionsPrefix(pep1,"qeig1_")
      PEPSetOptionsPrefix(pep2,"qeig2_")
.ve

   Level: advanced

.seealso: PEPAppendOptionsPrefix(), PEPGetOptionsPrefix()
@*/
PetscErrorCode PEPSetOptionsPrefix(PEP pep,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!pep->st) PetscCall(PEPGetST(pep,&pep->st));
  PetscCall(STSetOptionsPrefix(pep->st,prefix));
  if (!pep->V) PetscCall(PEPGetBV(pep,&pep->V));
  PetscCall(BVSetOptionsPrefix(pep->V,prefix));
  if (!pep->ds) PetscCall(PEPGetDS(pep,&pep->ds));
  PetscCall(DSSetOptionsPrefix(pep->ds,prefix));
  if (!pep->rg) PetscCall(PEPGetRG(pep,&pep->rg));
  PetscCall(RGSetOptionsPrefix(pep->rg,prefix));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)pep,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPAppendOptionsPrefix - Appends to the prefix used for searching for all
   PEP options in the database.

   Logically Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
-  prefix - the prefix string to prepend to all PEP option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: PEPSetOptionsPrefix(), PEPGetOptionsPrefix()
@*/
PetscErrorCode PEPAppendOptionsPrefix(PEP pep,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!pep->st) PetscCall(PEPGetST(pep,&pep->st));
  PetscCall(STAppendOptionsPrefix(pep->st,prefix));
  if (!pep->V) PetscCall(PEPGetBV(pep,&pep->V));
  PetscCall(BVAppendOptionsPrefix(pep->V,prefix));
  if (!pep->ds) PetscCall(PEPGetDS(pep,&pep->ds));
  PetscCall(DSAppendOptionsPrefix(pep->ds,prefix));
  if (!pep->rg) PetscCall(PEPGetRG(pep,&pep->rg));
  PetscCall(RGAppendOptionsPrefix(pep->rg,prefix));
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)pep,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetOptionsPrefix - Gets the prefix used for searching for all
   PEP options in the database.

   Not Collective

   Input Parameters:
.  pep - the polynomial eigensolver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: PEPSetOptionsPrefix(), PEPAppendOptionsPrefix()
@*/
PetscErrorCode PEPGetOptionsPrefix(PEP pep,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)pep,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}
