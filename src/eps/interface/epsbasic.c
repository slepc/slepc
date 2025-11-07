/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic EPS routines
*/

#include <slepc/private/epsimpl.h>      /*I "slepceps.h" I*/

/* Logging support */
PetscClassId      EPS_CLASSID = 0;
PetscLogEvent     EPS_SetUp = 0,EPS_Solve = 0,EPS_CISS_SVD = 0;

/* List of registered EPS routines */
PetscFunctionList EPSList = NULL;
PetscBool         EPSRegisterAllCalled = PETSC_FALSE;

/* List of registered EPS monitors */
PetscFunctionList EPSMonitorList              = NULL;
PetscFunctionList EPSMonitorCreateList        = NULL;
PetscFunctionList EPSMonitorDestroyList       = NULL;
PetscBool         EPSMonitorRegisterAllCalled = PETSC_FALSE;

/*@
   EPSCreate - Creates the `EPS` context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  outeps - location to put the `EPS` context

   Note:
   The default `EPS` type is `EPSKRYLOVSCHUR`.

   Level: beginner

.seealso: [](ch:eps), `EPSSetUp()`, `EPSSolve()`, `EPSDestroy()`, `EPS`
@*/
PetscErrorCode EPSCreate(MPI_Comm comm,EPS *outeps)
{
  EPS            eps;

  PetscFunctionBegin;
  PetscAssertPointer(outeps,2);
  PetscCall(EPSInitializePackage());
  PetscCall(SlepcHeaderCreate(eps,EPS_CLASSID,"EPS","Eigenvalue Problem Solver","EPS",comm,EPSDestroy,EPSView));

  eps->max_it          = PETSC_DETERMINE;
  eps->nev             = 0;
  eps->ncv             = PETSC_DETERMINE;
  eps->mpd             = PETSC_DETERMINE;
  eps->nini            = 0;
  eps->nds             = 0;
  eps->target          = 0.0;
  eps->tol             = PETSC_DETERMINE;
  eps->thres           = PETSC_MIN_REAL;
  eps->threlative      = PETSC_FALSE;
  eps->conv            = EPS_CONV_REL;
  eps->stop            = EPS_STOP_BASIC;
  eps->which           = (EPSWhich)0;
  eps->inta            = 0.0;
  eps->intb            = 0.0;
  eps->problem_type    = (EPSProblemType)0;
  eps->extraction      = EPS_RITZ;
  eps->balance         = EPS_BALANCE_NONE;
  eps->balance_its     = 5;
  eps->balance_cutoff  = 1e-8;
  eps->trueres         = PETSC_FALSE;
  eps->trackall        = PETSC_FALSE;
  eps->purify          = PETSC_TRUE;
  eps->twosided        = PETSC_FALSE;

  eps->converged       = EPSConvergedRelative;
  eps->convergeduser   = NULL;
  eps->convergeddestroy= NULL;
  eps->stopping        = EPSStoppingBasic;
  eps->stoppinguser    = NULL;
  eps->stoppingdestroy = NULL;
  eps->arbitrary       = NULL;
  eps->convergedctx    = NULL;
  eps->stoppingctx     = NULL;
  eps->arbitraryctx    = NULL;
  eps->numbermonitors  = 0;

  eps->st              = NULL;
  eps->ds              = NULL;
  eps->V               = NULL;
  eps->W               = NULL;
  eps->rg              = NULL;
  eps->D               = NULL;
  eps->IS              = NULL;
  eps->ISL             = NULL;
  eps->defl            = NULL;
  eps->eigr            = NULL;
  eps->eigi            = NULL;
  eps->errest          = NULL;
  eps->rr              = NULL;
  eps->ri              = NULL;
  eps->perm            = NULL;
  eps->nwork           = 0;
  eps->work            = NULL;
  eps->data            = NULL;

  eps->state           = EPS_STATE_INITIAL;
  eps->categ           = EPS_CATEGORY_KRYLOV;
  eps->nconv           = 0;
  eps->its             = 0;
  eps->nloc            = 0;
  eps->nrma            = 0.0;
  eps->nrmb            = 0.0;
  eps->useds           = PETSC_FALSE;
  eps->isgeneralized   = PETSC_FALSE;
  eps->ispositive      = PETSC_FALSE;
  eps->ishermitian     = PETSC_FALSE;
  eps->isstructured    = PETSC_FALSE;
  eps->reason          = EPS_CONVERGED_ITERATING;

  PetscCall(PetscNew(&eps->sc));
  *outeps = eps;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetType - Selects the particular solver to be used in the `EPS` object.

   Logically Collective

   Input Parameters:
+  eps  - the linear eigensolver context
-  type - a known method

   Options Database Key:
.  -eps_type \<type\> - sets the method; use `-help` for a list of available methods

   Notes:
   See `EPSType` for available methods. The default is `EPSKRYLOVSCHUR`.

   Normally, it is best to use the `EPSSetFromOptions()` command and
   then set the `EPS` type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The `EPSSetType()` routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: [](ch:eps), `STSetType()`, `EPSType`
@*/
PetscErrorCode EPSSetType(EPS eps,EPSType type)
{
  PetscErrorCode (*r)(EPS);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)eps,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscFunctionListFind(EPSList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown EPS type given: %s",type);

  PetscTryTypeMethod(eps,destroy);
  PetscCall(PetscMemzero(eps->ops,sizeof(struct _EPSOps)));

  eps->state = EPS_STATE_INITIAL;
  PetscCall(PetscObjectChangeTypeName((PetscObject)eps,type));
  PetscCall((*r)(eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetType - Gets the `EPS` type as a string from the `EPS` object.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  type - name of `EPS` method

   Level: intermediate

.seealso: [](ch:eps), `EPSSetType()`
@*/
PetscErrorCode EPSGetType(EPS eps,EPSType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)eps)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSRegister - Adds a method to the eigenproblem solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Note:
   `EPSRegister()` may be called multiple times to add several user-defined solvers.

   Example Usage:
.vb
   EPSRegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
.vb
   EPSSetType(eps,"my_solver")
.ve
   or at runtime via the option `-eps_type my_solver`.

   Level: advanced

.seealso: [](ch:eps), `EPSRegisterAll()`
@*/
PetscErrorCode EPSRegister(const char *name,PetscErrorCode (*function)(EPS))
{
  PetscFunctionBegin;
  PetscCall(EPSInitializePackage());
  PetscCall(PetscFunctionListAdd(&EPSList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorRegister - Registers an `EPS` monitor routine that may be accessed with
   `EPSMonitorSetFromOptions()`.

   Not Collective

   Input Parameters:
+  name    - name of a new monitor routine
.  vtype   - a `PetscViewerType` for the output
.  format  - a `PetscViewerFormat` for the output
.  monitor - monitor routine, see `EPSMonitorRegisterFn`
.  create  - creation routine, or `NULL`
-  destroy - destruction routine, or `NULL`

   Notes:
   `EPSMonitorRegister()` may be called multiple times to add several user-defined monitors.

   The calling sequence for the given function matches the calling sequence of `EPSMonitorFn`
   functions passed to `EPSMonitorSet()` with the additional requirement that its final argument
   be a `PetscViewerAndFormat`.

   Example Usage:
.vb
   EPSMonitorRegister("my_monitor",PETSCVIEWERASCII,PETSC_VIEWER_ASCII_INFO_DETAIL,MyMonitor,NULL,NULL);
.ve

   Then, your monitor can be chosen with the procedural interface via
.vb
   EPSMonitorSetFromOptions(eps,"-eps_monitor_my_monitor","my_monitor",NULL)
.ve
   or at runtime via the option `-eps_monitor_my_monitor`.

   Level: advanced

.seealso: [](ch:eps), `EPSMonitorSet()`, `EPSMonitorRegisterAll()`, `EPSMonitorSetFromOptions()`
@*/
PetscErrorCode EPSMonitorRegister(const char name[],PetscViewerType vtype,PetscViewerFormat format,EPSMonitorRegisterFn *monitor,EPSMonitorRegisterCreateFn *create,EPSMonitorRegisterDestroyFn *destroy)
{
  char           key[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  PetscCall(EPSInitializePackage());
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListAdd(&EPSMonitorList,key,monitor));
  if (create)  PetscCall(PetscFunctionListAdd(&EPSMonitorCreateList,key,create));
  if (destroy) PetscCall(PetscFunctionListAdd(&EPSMonitorDestroyList,key,destroy));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSReset - Resets the `EPS` context to the initial state (prior to setup)
   and destroys any allocated `Vec`s and `Mat`s.

   Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Note:
   This can be used when a problem of different matrix size wants to be solved.
   All options that have previously been set are preserved, so in a next use
   the solver configuration is the same, but new sizes for matrices and vectors
   are allowed.

   Level: advanced

.seealso: [](ch:eps), `EPSDestroy()`
@*/
PetscErrorCode EPSReset(EPS eps)
{
  PetscFunctionBegin;
  if (eps) PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!eps) PetscFunctionReturn(PETSC_SUCCESS);
  PetscTryTypeMethod(eps,reset);
  if (eps->st) PetscCall(STReset(eps->st));
  PetscCall(VecDestroy(&eps->D));
  PetscCall(BVDestroy(&eps->V));
  PetscCall(BVDestroy(&eps->W));
  PetscCall(VecDestroyVecs(eps->nwork,&eps->work));
  eps->nwork = 0;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSDestroy - Destroys the `EPS` context.

   Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Level: beginner

.seealso: [](ch:eps), `EPSCreate()`, `EPSSetUp()`, `EPSSolve()`
@*/
PetscErrorCode EPSDestroy(EPS *eps)
{
  PetscFunctionBegin;
  if (!*eps) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*eps,EPS_CLASSID,1);
  if (--((PetscObject)*eps)->refct > 0) { *eps = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscCall(EPSReset(*eps));
  PetscTryTypeMethod(*eps,destroy);
  if ((*eps)->eigr) PetscCall(PetscFree4((*eps)->eigr,(*eps)->eigi,(*eps)->errest,(*eps)->perm));
  if ((*eps)->rr) PetscCall(PetscFree2((*eps)->rr,(*eps)->ri));
  PetscCall(STDestroy(&(*eps)->st));
  PetscCall(RGDestroy(&(*eps)->rg));
  PetscCall(DSDestroy(&(*eps)->ds));
  PetscCall(PetscFree((*eps)->sc));
  /* just in case the initial vectors have not been used */
  PetscCall(SlepcBasisDestroy_Private(&(*eps)->nds,&(*eps)->defl));
  PetscCall(SlepcBasisDestroy_Private(&(*eps)->nini,&(*eps)->IS));
  PetscCall(SlepcBasisDestroy_Private(&(*eps)->ninil,&(*eps)->ISL));
  if ((*eps)->convergeddestroy) PetscCall((*(*eps)->convergeddestroy)(&(*eps)->convergedctx));
  if ((*eps)->stoppingdestroy) PetscCall((*(*eps)->stoppingdestroy)(&(*eps)->stoppingctx));
  PetscCall(EPSMonitorCancel(*eps));
  PetscCall(PetscHeaderDestroy(eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetTarget - Sets the value of the target.

   Logically Collective

   Input Parameters:
+  eps    - the linear eigensolver context
-  target - the value of the target

   Options Database Key:
.  -eps_target \<target\> - the value of the target

   Notes:
   The target is a scalar value used to determine the portion of the spectrum
   of interest. It is used in combination with `EPSSetWhichEigenpairs()`.

   When PETSc is built with real scalars, it is not possible to specify a
   complex target.

   In the case of complex scalars, a complex value can be provided in the
   command line with `[+/-][realnumber][+/-]realnumberi` with no spaces, e.g.
   `-eps_target 1.0+2.0i`.

   Level: intermediate

.seealso: [](ch:eps), `EPSGetTarget()`, `EPSSetWhichEigenpairs()`
@*/
PetscErrorCode EPSSetTarget(EPS eps,PetscScalar target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveScalar(eps,target,2);
  eps->target = target;
  if (!eps->st) PetscCall(EPSGetST(eps,&eps->st));
  PetscCall(STSetDefaultShift(eps->st,target));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetTarget - Gets the value of the target.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  target - the value of the target

   Note:
   If the target was not set by the user, then zero is returned.

   Level: intermediate

.seealso: [](ch:eps), `EPSSetTarget()`
@*/
PetscErrorCode EPSGetTarget(EPS eps,PetscScalar* target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(target,2);
  *target = eps->target;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetInterval - Defines the computational interval for spectrum slicing.

   Logically Collective

   Input Parameters:
+  eps  - the linear eigensolver context
.  inta - left end of the interval
-  intb - right end of the interval

   Options Database Key:
.  -eps_interval <a,b> - set $[a,b]$ as the interval of interest

   Notes:
   Spectrum slicing is a technique employed for computing all eigenvalues of
   symmetric eigenproblems in a given interval, see section [](#sec:slice).
   This function provides the interval to be considered. It must be used in
   combination with `EPS_ALL`, see `EPSSetWhichEigenpairs()`.

   A computational interval is also needed when using polynomial filters,
   see `STFILTER` and section [](#sec:filter).

   In the command-line option, two values must be provided. For an open interval,
   one can give an infinite, e.g., `-eps_interval 1.0,inf` or `-eps_interval -inf,1.0`.
   An open interval in the programmatic interface can be specified with
   `PETSC_MAX_REAL` and -`PETSC_MAX_REAL`.

   Level: intermediate

.seealso: [](ch:eps), [](#sec:slice), [](#sec:filter), `EPSGetInterval()`, `EPSSetWhichEigenpairs()`, `STFILTER`
@*/
PetscErrorCode EPSSetInterval(EPS eps,PetscReal inta,PetscReal intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,inta,2);
  PetscValidLogicalCollectiveReal(eps,intb,3);
  PetscCheck(inta<intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be inta<intb");
  if (eps->inta != inta || eps->intb != intb) {
    eps->inta = inta;
    eps->intb = intb;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetInterval - Gets the computational interval for spectrum slicing.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameters:
+  inta - left end of the interval
-  intb - right end of the interval

   Level: intermediate

   Note:
   If the interval was not set by the user, then zeros are returned.

.seealso: [](ch:eps), `EPSSetInterval()`
@*/
PetscErrorCode EPSGetInterval(EPS eps,PetscReal* inta,PetscReal* intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (inta) *inta = eps->inta;
  if (intb) *intb = eps->intb;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetST - Associates a spectral transformation object to the eigensolver.

   Collective

   Input Parameters:
+  eps - the linear eigensolver context
-  st   - the spectral transformation object

   Note:
   Use `EPSGetST()` to retrieve the spectral transformation context at a later time
   (for example, to free it at the end of the computations).

   Level: advanced

.seealso: [](ch:eps), `EPSGetST()`
@*/
PetscErrorCode EPSSetST(EPS eps,ST st)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(st,ST_CLASSID,2);
  PetscCheckSameComm(eps,1,st,2);
  PetscCall(PetscObjectReference((PetscObject)st));
  PetscCall(STDestroy(&eps->st));
  eps->st = st;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetST - Obtain the spectral transformation (`ST`) object associated
   to the eigensolver object.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  st - spectral transformation context

   Level: intermediate

.seealso: [](ch:eps), `EPSSetST()`
@*/
PetscErrorCode EPSGetST(EPS eps,ST *st)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(st,2);
  if (!eps->st) {
    PetscCall(STCreate(PetscObjectComm((PetscObject)eps),&eps->st));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)eps->st,(PetscObject)eps,0));
    PetscCall(PetscObjectSetOptions((PetscObject)eps->st,((PetscObject)eps)->options));
  }
  *st = eps->st;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetBV - Associates a basis vectors object to the eigensolver.

   Collective

   Input Parameters:
+  eps - the linear eigensolver context
-  V   - the basis vectors object

   Level: advanced

.seealso: [](ch:eps), `EPSGetBV()`
@*/
PetscErrorCode EPSSetBV(EPS eps,BV V)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(V,BV_CLASSID,2);
  PetscCheckSameComm(eps,1,V,2);
  PetscCall(PetscObjectReference((PetscObject)V));
  PetscCall(BVDestroy(&eps->V));
  eps->V = V;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetBV - Obtain the basis vectors object associated to the eigensolver object.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  V - basis vectors context

   Level: advanced

.seealso: [](ch:eps), `EPSSetBV()`
@*/
PetscErrorCode EPSGetBV(EPS eps,BV *V)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(V,2);
  if (!eps->V) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)eps),&eps->V));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)eps->V,(PetscObject)eps,0));
    PetscCall(PetscObjectSetOptions((PetscObject)eps->V,((PetscObject)eps)->options));
  }
  *V = eps->V;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetRG - Associates a region object to the eigensolver.

   Collective

   Input Parameters:
+  eps - the linear eigensolver context
-  rg  - the region object

   Note:
   Use `EPSGetRG()` to retrieve the region context at a later time (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: [](ch:eps), `EPSGetRG()`
@*/
PetscErrorCode EPSSetRG(EPS eps,RG rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (rg) {
    PetscValidHeaderSpecific(rg,RG_CLASSID,2);
    PetscCheckSameComm(eps,1,rg,2);
  }
  PetscCall(PetscObjectReference((PetscObject)rg));
  PetscCall(RGDestroy(&eps->rg));
  eps->rg = rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetRG - Obtain the region object associated to the eigensolver.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  rg - region context

   Level: advanced

.seealso: [](ch:eps), `EPSSetRG()`
@*/
PetscErrorCode EPSGetRG(EPS eps,RG *rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(rg,2);
  if (!eps->rg) {
    PetscCall(RGCreate(PetscObjectComm((PetscObject)eps),&eps->rg));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)eps->rg,(PetscObject)eps,0));
    PetscCall(PetscObjectSetOptions((PetscObject)eps->rg,((PetscObject)eps)->options));
  }
  *rg = eps->rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetDS - Associates a direct solver object to the eigensolver.

   Collective

   Input Parameters:
+  eps - the linear eigensolver context
-  ds  - the direct solver object

   Note:
   Use `EPSGetDS()` to retrieve the direct solver context at a later time (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: [](ch:eps), `EPSGetDS()`
@*/
PetscErrorCode EPSSetDS(EPS eps,DS ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(ds,DS_CLASSID,2);
  PetscCheckSameComm(eps,1,ds,2);
  PetscCall(PetscObjectReference((PetscObject)ds));
  PetscCall(DSDestroy(&eps->ds));
  eps->ds = ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetDS - Obtain the direct solver object associated to the eigensolver object.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  ds - direct solver context

   Level: advanced

.seealso: [](ch:eps), `EPSSetDS()`
@*/
PetscErrorCode EPSGetDS(EPS eps,DS *ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(ds,2);
  if (!eps->ds) {
    PetscCall(DSCreate(PetscObjectComm((PetscObject)eps),&eps->ds));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)eps->ds,(PetscObject)eps,0));
    PetscCall(PetscObjectSetOptions((PetscObject)eps->ds,((PetscObject)eps)->options));
  }
  *ds = eps->ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSIsGeneralized - Ask if the `EPS` object corresponds to a generalized
   eigenvalue problem.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  is - `PETSC_TRUE` if the problem is generalized

   Level: intermediate

.seealso: [](ch:eps), `EPSIsHermitian()`, `EPSIsPositive()`, `EPSIsStructured()`
@*/
PetscErrorCode EPSIsGeneralized(EPS eps,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(is,2);
  *is = eps->isgeneralized;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSIsHermitian - Ask if the `EPS` object corresponds to a Hermitian
   eigenvalue problem.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  is - `PETSC_TRUE` if the problem is Hermitian

   Level: intermediate

.seealso: [](ch:eps), `EPSIsGeneralized()`, `EPSIsPositive()`, `EPSIsStructured()`
@*/
PetscErrorCode EPSIsHermitian(EPS eps,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(is,2);
  *is = eps->ishermitian;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSIsPositive - Ask if the `EPS` object corresponds to an eigenvalue
   problem type that requires a positive (semi-) definite matrix $B$.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  is - `PETSC_TRUE` if the problem is positive (semi-) definite

   Level: intermediate

.seealso: [](ch:eps), `EPSIsGeneralized()`, `EPSIsHermitian()`, `EPSIsStructured()`
@*/
PetscErrorCode EPSIsPositive(EPS eps,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(is,2);
  *is = eps->ispositive;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSIsStructured - Ask if the `EPS` object corresponds to a structured
   eigenvalue problem.

   Not Collective

   Input Parameter:
.  eps - the linear eigensolver context

   Output Parameter:
.  is - `PETSC_TRUE` if the problem is structured

   Note:
   The result will be true if the problem type has been set to some
   structured type such as `EPS_BSE`. This is independent of whether the
   input matrix has been built with a certain structure with a helper function.

   Level: intermediate

.seealso: [](ch:eps), `EPSIsGeneralized()`, `EPSIsHermitian()`, `EPSIsPositive()`, `EPSSetProblemType()`
@*/
PetscErrorCode EPSIsStructured(EPS eps,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(is,2);
  *is = eps->isstructured;
  PetscFunctionReturn(PETSC_SUCCESS);
}
