/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   MFN routines related to options that can be set via the command-line
   or procedurally
*/

#include <slepc/private/mfnimpl.h>   /*I "slepcmfn.h" I*/
#include <petscdraw.h>

/*@C
   MFNMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective

   Input Parameters:
+  mfn  - the matrix function context
.  opt  - the command line option for this monitor
.  name - the monitor type one is seeking
-  ctx  - an optional user context for the monitor, or NULL

   Level: developer

.seealso: MFNMonitorSet()
@*/
PetscErrorCode MFNMonitorSetFromOptions(MFN mfn,const char opt[],const char name[],void *ctx)
{
  PetscErrorCode       (*mfunc)(MFN,PetscInt,PetscReal,void*);
  PetscErrorCode       (*cfunc)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
  PetscErrorCode       (*dfunc)(PetscViewerAndFormat**);
  PetscViewerAndFormat *vf;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerType      vtype;
  char                 key[PETSC_MAX_PATH_LEN];
  PetscBool            flg;

  PetscFunctionBegin;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)mfn),((PetscObject)mfn)->options,((PetscObject)mfn)->prefix,opt,&viewer,&format,&flg));
  if (!flg) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetType(viewer,&vtype));
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListFind(MFNMonitorList,key,&mfunc));
  PetscCheck(mfunc,PetscObjectComm((PetscObject)mfn),PETSC_ERR_SUP,"Specified viewer and format not supported");
  PetscCall(PetscFunctionListFind(MFNMonitorCreateList,key,&cfunc));
  PetscCall(PetscFunctionListFind(MFNMonitorDestroyList,key,&dfunc));
  if (!cfunc) cfunc = PetscViewerAndFormatCreate_Internal;
  if (!dfunc) dfunc = PetscViewerAndFormatDestroy;

  PetscCall((*cfunc)(viewer,format,ctx,&vf));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(MFNMonitorSet(mfn,mfunc,vf,(PetscCtxDestroyFn*)dfunc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetFromOptions - Sets MFN options from the options database.
   This routine must be called before MFNSetUp() if the user is to be
   allowed to set the solver type.

   Collective

   Input Parameters:
.  mfn - the matrix function context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner

.seealso: MFNSetOptionsPrefix()
@*/
PetscErrorCode MFNSetFromOptions(MFN mfn)
{
  char           type[256];
  PetscBool      set,flg,flg1,flg2;
  PetscReal      r;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscCall(MFNRegisterAll());
  PetscObjectOptionsBegin((PetscObject)mfn);
    PetscCall(PetscOptionsFList("-mfn_type","Matrix Function method","MFNSetType",MFNList,(char*)(((PetscObject)mfn)->type_name?((PetscObject)mfn)->type_name:MFNKRYLOV),type,sizeof(type),&flg));
    if (flg) PetscCall(MFNSetType(mfn,type));
    else if (!((PetscObject)mfn)->type_name) PetscCall(MFNSetType(mfn,MFNKRYLOV));

    i = mfn->max_it;
    PetscCall(PetscOptionsInt("-mfn_max_it","Maximum number of iterations","MFNSetTolerances",mfn->max_it,&i,&flg1));
    if (!flg1) i = PETSC_DETERMINE;
    r = mfn->tol;
    PetscCall(PetscOptionsReal("-mfn_tol","Tolerance","MFNSetTolerances",SlepcDefaultTol(mfn->tol),&r,&flg2));
    if (flg1 || flg2) PetscCall(MFNSetTolerances(mfn,r,i));

    PetscCall(PetscOptionsInt("-mfn_ncv","Number of basis vectors","MFNSetDimensions",mfn->ncv,&i,&flg));
    if (flg) PetscCall(MFNSetDimensions(mfn,i));

    PetscCall(PetscOptionsBool("-mfn_error_if_not_converged","Generate error if solver does not converge","MFNSetErrorIfNotConverged",mfn->errorifnotconverged,&mfn->errorifnotconverged,NULL));

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to MFNSetFromOptions()
    */
    PetscCall(PetscOptionsBool("-mfn_monitor_cancel","Remove any hardwired monitor routines","MFNMonitorCancel",PETSC_FALSE,&flg,&set));
    if (set && flg) PetscCall(MFNMonitorCancel(mfn));
    PetscCall(MFNMonitorSetFromOptions(mfn,"-mfn_monitor","error_estimate",NULL));

    /* -----------------------------------------------------------------------*/
    PetscCall(PetscOptionsName("-mfn_view","Print detailed information on solver used","MFNView",&set));

    PetscTryTypeMethod(mfn,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)mfn,PetscOptionsObject));
  PetscOptionsEnd();

  if (!mfn->V) PetscCall(MFNGetBV(mfn,&mfn->V));
  PetscCall(BVSetFromOptions(mfn->V));
  if (!mfn->fn) PetscCall(MFNGetFN(mfn,&mfn->fn));
  PetscCall(FNSetFromOptions(mfn->fn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNGetTolerances - Gets the tolerance and maximum iteration count used
   by the MFN convergence tests.

   Not Collective

   Input Parameter:
.  mfn - the matrix function context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: MFNSetTolerances()
@*/
PetscErrorCode MFNGetTolerances(MFN mfn,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  if (tol)    *tol    = mfn->tol;
  if (maxits) *maxits = mfn->max_it;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetTolerances - Sets the tolerance and maximum iteration count used
   by the MFN convergence tests.

   Logically Collective

   Input Parameters:
+  mfn - the matrix function context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -mfn_tol <tol> - Sets the convergence tolerance
-  -mfn_max_it <maxits> - Sets the maximum number of iterations allowed

   Notes:
   Use PETSC_CURRENT to retain the current value of any of the parameters.
   Use PETSC_DETERMINE for either argument to assign a default value computed
   internally (may be different in each solver).
   For maxits use PETSC_UMLIMITED to indicate there is no upper bound on this value.

   Level: intermediate

.seealso: MFNGetTolerances()
@*/
PetscErrorCode MFNSetTolerances(MFN mfn,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidLogicalCollectiveReal(mfn,tol,2);
  PetscValidLogicalCollectiveInt(mfn,maxits,3);
  if (tol == (PetscReal)PETSC_DETERMINE) {
    mfn->tol = PETSC_DETERMINE;
    mfn->setupcalled = 0;
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    mfn->tol = tol;
  }
  if (maxits == PETSC_DETERMINE) {
    mfn->max_it = PETSC_DETERMINE;
    mfn->setupcalled = 0;
  } else if (maxits == PETSC_UNLIMITED) {
    mfn->max_it = PETSC_INT_MAX;
  } else if (maxits != PETSC_CURRENT) {
    PetscCheck(maxits>0,PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    mfn->max_it = maxits;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNGetDimensions - Gets the dimension of the subspace used by the solver.

   Not Collective

   Input Parameter:
.  mfn - the matrix function context

   Output Parameter:
.  ncv - the maximum dimension of the subspace to be used by the solver

   Level: intermediate

.seealso: MFNSetDimensions()
@*/
PetscErrorCode MFNGetDimensions(MFN mfn,PetscInt *ncv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(ncv,2);
  *ncv = mfn->ncv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetDimensions - Sets the dimension of the subspace to be used by the solver.

   Logically Collective

   Input Parameters:
+  mfn - the matrix function context
-  ncv - the maximum dimension of the subspace to be used by the solver

   Options Database Keys:
.  -mfn_ncv <ncv> - Sets the dimension of the subspace

   Notes:
   Use PETSC_DETERMINE for ncv to assign a reasonably good value, which is
   dependent on the solution method.

   Level: intermediate

.seealso: MFNGetDimensions()
@*/
PetscErrorCode MFNSetDimensions(MFN mfn,PetscInt ncv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidLogicalCollectiveInt(mfn,ncv,2);
  if (ncv == PETSC_DECIDE || ncv == PETSC_DEFAULT) {
    mfn->ncv = PETSC_DETERMINE;
  } else {
    PetscCheck(ncv>0,PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    mfn->ncv = ncv;
  }
  mfn->setupcalled = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetErrorIfNotConverged - Causes MFNSolve() to generate an error if the
   solver has not converged.

   Logically Collective

   Input Parameters:
+  mfn - the matrix function context
-  flg - PETSC_TRUE indicates you want the error generated

   Options Database Keys:
.  -mfn_error_if_not_converged - this takes an optional truth value (0/1/no/yes/true/false)

   Level: intermediate

   Note:
   Normally SLEPc continues if the solver fails to converge, you can call
   MFNGetConvergedReason() after a MFNSolve() to determine if it has converged.

.seealso: MFNGetErrorIfNotConverged()
@*/
PetscErrorCode MFNSetErrorIfNotConverged(MFN mfn,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidLogicalCollectiveBool(mfn,flg,2);
  mfn->errorifnotconverged = flg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNGetErrorIfNotConverged - Return a flag indicating whether MFNSolve() will
   generate an error if the solver does not converge.

   Not Collective

   Input Parameter:
.  mfn - the matrix function context

   Output Parameter:
.  flag - PETSC_TRUE if it will generate an error, else PETSC_FALSE

   Level: intermediate

.seealso: MFNSetErrorIfNotConverged()
@*/
PetscErrorCode MFNGetErrorIfNotConverged(MFN mfn,PetscBool *flag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(flag,2);
  *flag = mfn->errorifnotconverged;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetOptionsPrefix - Sets the prefix used for searching for all
   MFN options in the database.

   Logically Collective

   Input Parameters:
+  mfn - the matrix function context
-  prefix - the prefix string to prepend to all MFN option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different MFN contexts, one could call
.vb
      MFNSetOptionsPrefix(mfn1,"fun1_")
      MFNSetOptionsPrefix(mfn2,"fun2_")
.ve

   Level: advanced

.seealso: MFNAppendOptionsPrefix(), MFNGetOptionsPrefix()
@*/
PetscErrorCode MFNSetOptionsPrefix(MFN mfn,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  if (!mfn->V) PetscCall(MFNGetBV(mfn,&mfn->V));
  PetscCall(BVSetOptionsPrefix(mfn->V,prefix));
  if (!mfn->fn) PetscCall(MFNGetFN(mfn,&mfn->fn));
  PetscCall(FNSetOptionsPrefix(mfn->fn,prefix));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)mfn,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNAppendOptionsPrefix - Appends to the prefix used for searching for all
   MFN options in the database.

   Logically Collective

   Input Parameters:
+  mfn - the matrix function context
-  prefix - the prefix string to prepend to all MFN option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: MFNSetOptionsPrefix(), MFNGetOptionsPrefix()
@*/
PetscErrorCode MFNAppendOptionsPrefix(MFN mfn,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  if (!mfn->V) PetscCall(MFNGetBV(mfn,&mfn->V));
  PetscCall(BVAppendOptionsPrefix(mfn->V,prefix));
  if (!mfn->fn) PetscCall(MFNGetFN(mfn,&mfn->fn));
  PetscCall(FNAppendOptionsPrefix(mfn->fn,prefix));
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)mfn,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNGetOptionsPrefix - Gets the prefix used for searching for all
   MFN options in the database.

   Not Collective

   Input Parameters:
.  mfn - the matrix function context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: MFNSetOptionsPrefix(), MFNAppendOptionsPrefix()
@*/
PetscErrorCode MFNGetOptionsPrefix(MFN mfn,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)mfn,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}
