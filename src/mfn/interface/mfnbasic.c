/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic MFN routines
*/

#include <slepc/private/mfnimpl.h>      /*I "slepcmfn.h" I*/

/* Logging support */
PetscClassId      MFN_CLASSID = 0;
PetscLogEvent     MFN_SetUp = 0,MFN_Solve = 0;

/* List of registered MFN routines */
PetscFunctionList MFNList = NULL;
PetscBool         MFNRegisterAllCalled = PETSC_FALSE;

/* List of registered MFN monitors */
PetscFunctionList MFNMonitorList              = NULL;
PetscFunctionList MFNMonitorCreateList        = NULL;
PetscFunctionList MFNMonitorDestroyList       = NULL;
PetscBool         MFNMonitorRegisterAllCalled = PETSC_FALSE;

/*@
   MFNView - Prints the `MFN` data structure.

   Collective

   Input Parameters:
+  mfn - the matrix function solver context
-  viewer - optional visualization context

   Options Database Key:
.  -mfn_view - calls `MFNView()` at end of `MFNSolve()`

   Notes:
   The available visualization contexts include
+     `PETSC_VIEWER_STDOUT_SELF` - standard output (default)
-     `PETSC_VIEWER_STDOUT_WORLD` - synchronized standard output where only the
         first process opens the file; all other processes send their data to the
         first one to print

   The user can open an alternative visualization context with `PetscViewerASCIIOpen()`
   to output to a specified file.

   Level: beginner

.seealso: [](ch:mfn), `MFNCreate()`, `MFNViewFromOptions()`
@*/
PetscErrorCode MFNView(MFN mfn,PetscViewer viewer)
{
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)mfn),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(mfn,1,viewer,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)mfn,viewer));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscTryTypeMethod(mfn,view,viewer);
    PetscCall(PetscViewerASCIIPopTab(viewer));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %" PetscInt_FMT "\n",mfn->ncv));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %" PetscInt_FMT "\n",mfn->max_it));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)mfn->tol));
  } else PetscTryTypeMethod(mfn,view,viewer);
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO));
  if (!mfn->V) PetscCall(MFNGetFN(mfn,&mfn->fn));
  PetscCall(FNView(mfn->fn,viewer));
  if (!mfn->V) PetscCall(MFNGetBV(mfn,&mfn->V));
  PetscCall(BVView(mfn->V,viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNViewFromOptions - View (print) an `MFN` object based on values in the options database.

   Collective

   Input Parameters:
+  mfn  - the matrix function solver context
.  obj  - optional object that provides the options prefix used to query the options database
-  name - command line option

   Level: intermediate

.seealso: [](ch:mfn), `MFNView()`, `MFNCreate()`, `PetscObjectViewFromOptions()`
@*/
PetscErrorCode MFNViewFromOptions(MFN mfn,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)mfn,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*@
   MFNConvergedReasonView - Displays the reason an `MFN` solve converged or diverged.

   Collective

   Input Parameters:
+  mfn - the matrix function solver context
-  viewer - the viewer to display the reason

   Options Database Key:
.  -mfn_converged_reason - print reason for convergence/divergence, and number of iterations

   Notes:
   Use `MFNConvergedReasonViewFromOptions()` to display the reason based on values
   in the options database.

   To change the format of the output call `PetscViewerPushFormat()` before this
   call. Use `PETSC_VIEWER_DEFAULT` for the default, or `PETSC_VIEWER_FAILED` to only
   display a reason if it fails. The latter can be set in the command line with
   `-mfn_converged_reason ::failed`.

   Level: intermediate

.seealso: [](ch:mfn), `MFNSetTolerances()`, `MFNGetIterationNumber()`, `MFNConvergedReasonViewFromOptions()`
@*/
PetscErrorCode MFNConvergedReasonView(MFN mfn,PetscViewer viewer)
{
  PetscBool         isAscii;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (!viewer) viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)mfn));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii));
  if (isAscii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)mfn)->tablevel));
    if (mfn->reason > 0 && format != PETSC_VIEWER_FAILED) PetscCall(PetscViewerASCIIPrintf(viewer,"%s Matrix function solve converged due to %s; iterations %" PetscInt_FMT "\n",((PetscObject)mfn)->prefix?((PetscObject)mfn)->prefix:"",MFNConvergedReasons[mfn->reason],mfn->its));
    else if (mfn->reason <= 0) PetscCall(PetscViewerASCIIPrintf(viewer,"%s Matrix function solve did not converge due to %s; iterations %" PetscInt_FMT "\n",((PetscObject)mfn)->prefix?((PetscObject)mfn)->prefix:"",MFNConvergedReasons[mfn->reason],mfn->its));
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)mfn)->tablevel));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNConvergedReasonViewFromOptions - Processes command line options to determine if/how
   the `MFN` converged reason is to be viewed.

   Collective

   Input Parameter:
.  mfn - the matrix function solver context

   Level: intermediate

.seealso: [](ch:mfn), `MFNConvergedReasonView()`
@*/
PetscErrorCode MFNConvergedReasonViewFromOptions(MFN mfn)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)mfn),((PetscObject)mfn)->options,((PetscObject)mfn)->prefix,"-mfn_converged_reason",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(MFNConvergedReasonView(mfn,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNCreate - Creates the `MFN` context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  outmfn - location to put the `MFN` context

   Note:
   The default `MFN` type is `MFNKRYLOV`.

   Level: beginner

.seealso: [](ch:mfn), `MFNSetUp()`, `MFNSolve()`, `MFNDestroy()`, `MFN`
@*/
PetscErrorCode MFNCreate(MPI_Comm comm,MFN *outmfn)
{
  MFN            mfn;

  PetscFunctionBegin;
  PetscAssertPointer(outmfn,2);
  PetscCall(MFNInitializePackage());
  PetscCall(SlepcHeaderCreate(mfn,MFN_CLASSID,"MFN","Matrix Function","MFN",comm,MFNDestroy,MFNView));

  mfn->A               = NULL;
  mfn->fn              = NULL;
  mfn->max_it          = PETSC_DETERMINE;
  mfn->ncv             = PETSC_DETERMINE;
  mfn->tol             = PETSC_DETERMINE;
  mfn->errorifnotconverged = PETSC_FALSE;

  mfn->numbermonitors  = 0;

  mfn->V               = NULL;
  mfn->nwork           = 0;
  mfn->work            = NULL;
  mfn->data            = NULL;

  mfn->its             = 0;
  mfn->nv              = 0;
  mfn->errest          = 0;
  mfn->setupcalled     = 0;
  mfn->reason          = MFN_CONVERGED_ITERATING;

  *outmfn = mfn;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetType - Selects the particular solver to be used in the `MFN` object.

   Logically Collective

   Input Parameters:
+  mfn  - the matrix function solver context
-  type - a known method

   Options Database Key:
.  -mfn_type \<type\> - sets the method; use `-help` for a list of available methods

   Notes:
   See `MFNType` for available methods. The default is `MFNKRYLOV`.

   Normally, it is best to use the `MFNSetFromOptions()` command and
   then set the `MFN` type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The `MFNSetType()` routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: [](ch:mfn), `MFNType`
@*/
PetscErrorCode MFNSetType(MFN mfn,MFNType type)
{
  PetscErrorCode (*r)(MFN);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)mfn,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscFunctionListFind(MFNList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown MFN type given: %s",type);

  PetscTryTypeMethod(mfn,destroy);
  PetscCall(PetscMemzero(mfn->ops,sizeof(struct _MFNOps)));

  mfn->setupcalled = 0;
  PetscCall(PetscObjectChangeTypeName((PetscObject)mfn,type));
  PetscCall((*r)(mfn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNGetType - Gets the `MFN` type as a string from the `MFN` object.

   Not Collective

   Input Parameter:
.  mfn - the matrix function solver context

   Output Parameter:
.  type - name of `MFN` method

   Level: intermediate

.seealso: [](ch:mfn), `MFNSetType()`
@*/
PetscErrorCode MFNGetType(MFN mfn,MFNType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)mfn)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   MFNRegister - Adds a method to the matrix function solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Note:
   `MFNRegister()` may be called multiple times to add several user-defined solvers.

   Example Usage:
.vb
   MFNRegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
.vb
   MFNSetType(mfn,"my_solver")
.ve
   or at runtime via the option `-mfn_type my_solver`.

   Level: advanced

.seealso: [](ch:mfn), `MFNRegisterAll()`
@*/
PetscErrorCode MFNRegister(const char *name,PetscErrorCode (*function)(MFN))
{
  PetscFunctionBegin;
  PetscCall(MFNInitializePackage());
  PetscCall(PetscFunctionListAdd(&MFNList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   MFNMonitorRegister - Registers an `MFN` monitor routine that may be accessed with
   `MFNMonitorSetFromOptions()`.

   Not Collective

   Input Parameters:
+  name    - name of a new monitor routine
.  vtype   - a `PetscViewerType` for the output
.  format  - a `PetscViewerFormat` for the output
.  monitor - monitor routine, see `MFNMonitorRegisterFn`
.  create  - creation routine, or `NULL`
-  destroy - destruction routine, or `NULL`

   Notes:
   `MFNMonitorRegister()` may be called multiple times to add several user-defined monitors.

   The calling sequence for the given function matches the calling sequence of `MFNMonitorFn`
   functions passed to `MFNMonitorSet()` with the additional requirement that its final argument
   be a `PetscViewerAndFormat`.

   Example Usage:
.vb
   MFNMonitorRegister("my_monitor",PETSCVIEWERASCII,PETSC_VIEWER_ASCII_INFO_DETAIL,MyMonitor,NULL,NULL);
.ve

   Then, your monitor can be chosen with the procedural interface via
.vb
   MFNMonitorSetFromOptions(mfn,"-mfn_monitor_my_monitor","my_monitor",NULL);
.ve
   or at runtime via the option `-mfn_monitor_my_monitor`.

   Level: advanced

.seealso: [](ch:mfn), `MFNMonitorSet()`, `MFNMonitorRegisterAll()`, `MFNMonitorSetFromOptions()`
@*/
PetscErrorCode MFNMonitorRegister(const char name[],PetscViewerType vtype,PetscViewerFormat format,MFNMonitorRegisterFn *monitor,MFNMonitorRegisterCreateFn *create,MFNMonitorRegisterDestroyFn *destroy)
{
  char           key[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  PetscCall(MFNInitializePackage());
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListAdd(&MFNMonitorList,key,monitor));
  if (create)  PetscCall(PetscFunctionListAdd(&MFNMonitorCreateList,key,create));
  if (destroy) PetscCall(PetscFunctionListAdd(&MFNMonitorDestroyList,key,destroy));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNReset - Resets the `MFN` context to the initial state (prior to setup)
   and destroys any allocated `Vec`s and `Mat`s.

   Collective

   Input Parameter:
.  mfn - the matrix function solver context

   Level: advanced

.seealso: [](ch:mfn), `MFNDestroy()`
@*/
PetscErrorCode MFNReset(MFN mfn)
{
  PetscFunctionBegin;
  if (mfn) PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  if (!mfn) PetscFunctionReturn(PETSC_SUCCESS);
  PetscTryTypeMethod(mfn,reset);
  PetscCall(MatDestroy(&mfn->A));
  PetscCall(BVDestroy(&mfn->V));
  PetscCall(VecDestroyVecs(mfn->nwork,&mfn->work));
  mfn->nwork = 0;
  mfn->setupcalled = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNDestroy - Destroys the `MFN` context.

   Collective

   Input Parameter:
.  mfn - the matrix function solver context

   Level: beginner

.seealso: [](ch:mfn), `MFNCreate()`, `MFNSetUp()`, `MFNSolve()`
@*/
PetscErrorCode MFNDestroy(MFN *mfn)
{
  PetscFunctionBegin;
  if (!*mfn) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*mfn,MFN_CLASSID,1);
  if (--((PetscObject)*mfn)->refct > 0) { *mfn = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscCall(MFNReset(*mfn));
  PetscTryTypeMethod(*mfn,destroy);
  PetscCall(FNDestroy(&(*mfn)->fn));
  PetscCall(MatDestroy(&(*mfn)->AT));
  PetscCall(MFNMonitorCancel(*mfn));
  PetscCall(PetscHeaderDestroy(mfn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetBV - Associates a basis vectors object to the matrix function solver.

   Collective

   Input Parameters:
+  mfn - the matrix function solver context
-  bv  - the basis vectors object

   Note:
   Use `MFNGetBV()` to retrieve the basis vectors context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: [](ch:mfn), `MFNGetBV()`
@*/
PetscErrorCode MFNSetBV(MFN mfn,BV bv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidHeaderSpecific(bv,BV_CLASSID,2);
  PetscCheckSameComm(mfn,1,bv,2);
  PetscCall(PetscObjectReference((PetscObject)bv));
  PetscCall(BVDestroy(&mfn->V));
  mfn->V = bv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNGetBV - Obtain the basis vectors object associated to the matrix
   function solver.

   Not Collective

   Input Parameter:
.  mfn - the matrix function solver context

   Output Parameter:
.  bv - basis vectors context

   Level: advanced

.seealso: [](ch:mfn), `MFNSetBV()`
@*/
PetscErrorCode MFNGetBV(MFN mfn,BV *bv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(bv,2);
  if (!mfn->V) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)mfn),&mfn->V));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)mfn->V,(PetscObject)mfn,0));
    PetscCall(PetscObjectSetOptions((PetscObject)mfn->V,((PetscObject)mfn)->options));
  }
  *bv = mfn->V;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetFN - Specifies the function to be computed.

   Collective

   Input Parameters:
+  mfn - the matrix function solver context
-  fn  - the math function object

   Notes:
   At a later time, use `MFNGetFN()` to retrieve the math function context
   (for example, to free it at the end of the computations).

   This function is not called in normal usage. Instead, it is easier to
   extract the internal `FN` object with `MFNGetFN()` and modify it.

   Level: beginner

.seealso: [](ch:mfn), `MFNGetFN()`
@*/
PetscErrorCode MFNSetFN(MFN mfn,FN fn)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidHeaderSpecific(fn,FN_CLASSID,2);
  PetscCheckSameComm(mfn,1,fn,2);
  PetscCall(PetscObjectReference((PetscObject)fn));
  PetscCall(FNDestroy(&mfn->fn));
  mfn->fn = fn;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNGetFN - Obtain the math function object associated to the `MFN` object.

   Not Collective

   Input Parameter:
.  mfn - the matrix function solver context

   Output Parameter:
.  fn - math function context

   Note:
   This is the usual way to specify the function that needs to be applied
   to a given vector in `MFNSolve()`.

   Level: beginner

.seealso: [](ch:mfn), `MFNSetFN()`. `MFNSolve()`
@*/
PetscErrorCode MFNGetFN(MFN mfn,FN *fn)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(fn,2);
  if (!mfn->fn) {
    PetscCall(FNCreate(PetscObjectComm((PetscObject)mfn),&mfn->fn));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)mfn->fn,(PetscObject)mfn,0));
    PetscCall(PetscObjectSetOptions((PetscObject)mfn->fn,((PetscObject)mfn)->options));
  }
  *fn = mfn->fn;
  PetscFunctionReturn(PETSC_SUCCESS);
}
