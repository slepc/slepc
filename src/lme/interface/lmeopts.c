/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   LME routines related to options that can be set via the command-line
   or procedurally
*/

#include <slepc/private/lmeimpl.h>   /*I "slepclme.h" I*/
#include <petscdraw.h>

/*@C
   LMEMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective

   Input Parameters:
+  lme  - the linear matrix equation solver context
.  opt  - the command line option for this monitor
.  name - the monitor type one is seeking
-  ctx  - an optional user context for the monitor, or `NULL`

   Level: developer

.seealso: [](ch:lme), `LMEMonitorSet()`
@*/
PetscErrorCode LMEMonitorSetFromOptions(LME lme,const char opt[],const char name[],void *ctx)
{
  PetscErrorCode       (*mfunc)(LME,PetscInt,PetscReal,void*);
  PetscErrorCode       (*cfunc)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);
  PetscErrorCode       (*dfunc)(PetscViewerAndFormat**);
  PetscViewerAndFormat *vf;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerType      vtype;
  char                 key[PETSC_MAX_PATH_LEN];
  PetscBool            flg;

  PetscFunctionBegin;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)lme),((PetscObject)lme)->options,((PetscObject)lme)->prefix,opt,&viewer,&format,&flg));
  if (!flg) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetType(viewer,&vtype));
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListFind(LMEMonitorList,key,&mfunc));
  PetscCheck(mfunc,PetscObjectComm((PetscObject)lme),PETSC_ERR_SUP,"Specified viewer and format not supported");
  PetscCall(PetscFunctionListFind(LMEMonitorCreateList,key,&cfunc));
  PetscCall(PetscFunctionListFind(LMEMonitorDestroyList,key,&dfunc));
  if (!cfunc) cfunc = PetscViewerAndFormatCreate_Internal;
  if (!dfunc) dfunc = PetscViewerAndFormatDestroy;

  PetscCall((*cfunc)(viewer,format,ctx,&vf));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(LMEMonitorSet(lme,mfunc,vf,(PetscErrorCode(*)(void **))dfunc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetFromOptions - Sets `LME` options from the options database.
   This routine must be called before `LMESetUp()` if the user is to be
   allowed to configure the solver.

   Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Note:
   To see all options, run your program with the `-help` option.

   Level: beginner

.seealso: [](ch:lme), `LMESetOptionsPrefix()`
@*/
PetscErrorCode LMESetFromOptions(LME lme)
{
  char           type[256];
  PetscBool      set,flg,flg1,flg2;
  PetscReal      r;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscCall(LMERegisterAll());
  PetscObjectOptionsBegin((PetscObject)lme);
    PetscCall(PetscOptionsFList("-lme_type","Linear matrix equation","LMESetType",LMEList,(char*)(((PetscObject)lme)->type_name?((PetscObject)lme)->type_name:LMEKRYLOV),type,sizeof(type),&flg));
    if (flg) PetscCall(LMESetType(lme,type));
    else if (!((PetscObject)lme)->type_name) PetscCall(LMESetType(lme,LMEKRYLOV));

    PetscCall(PetscOptionsBoolGroupBegin("-lme_lyapunov","Continuous-time Lyapunov equation","LMESetProblemType",&flg));
    if (flg) PetscCall(LMESetProblemType(lme,LME_LYAPUNOV));
    PetscCall(PetscOptionsBoolGroup("-lme_sylvester","Continuous-time Sylvester equation","LMESetProblemType",&flg));
    if (flg) PetscCall(LMESetProblemType(lme,LME_SYLVESTER));
    PetscCall(PetscOptionsBoolGroup("-lme_gen_lyapunov","Generalized Lyapunov equation","LMESetProblemType",&flg));
    if (flg) PetscCall(LMESetProblemType(lme,LME_GEN_LYAPUNOV));
    PetscCall(PetscOptionsBoolGroup("-lme_gen_sylvester","Generalized Sylvester equation","LMESetProblemType",&flg));
    if (flg) PetscCall(LMESetProblemType(lme,LME_GEN_SYLVESTER));
    PetscCall(PetscOptionsBoolGroup("-lme_dt_lyapunov","Discrete-time Lyapunov equation","LMESetProblemType",&flg));
    if (flg) PetscCall(LMESetProblemType(lme,LME_DT_LYAPUNOV));
    PetscCall(PetscOptionsBoolGroupEnd("-lme_stein","Stein equation","LMESetProblemType",&flg));
    if (flg) PetscCall(LMESetProblemType(lme,LME_STEIN));

    i = lme->max_it;
    PetscCall(PetscOptionsInt("-lme_max_it","Maximum number of iterations","LMESetTolerances",lme->max_it,&i,&flg1));
    if (!flg1) i = PETSC_DETERMINE;
    r = lme->tol;
    PetscCall(PetscOptionsReal("-lme_tol","Tolerance","LMESetTolerances",SlepcDefaultTol(lme->tol),&r,&flg2));
    if (flg1 || flg2) PetscCall(LMESetTolerances(lme,r,i));

    PetscCall(PetscOptionsInt("-lme_ncv","Number of basis vectors","LMESetDimensions",lme->ncv,&i,&flg));
    if (flg) PetscCall(LMESetDimensions(lme,i));

    PetscCall(PetscOptionsBool("-lme_error_if_not_converged","Generate error if solver does not converge","LMESetErrorIfNotConverged",lme->errorifnotconverged,&lme->errorifnotconverged,NULL));

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to LMESetFromOptions()
    */
    PetscCall(PetscOptionsBool("-lme_monitor_cancel","Remove any hardwired monitor routines","LMEMonitorCancel",PETSC_FALSE,&flg,&set));
    if (set && flg) PetscCall(LMEMonitorCancel(lme));
    PetscCall(LMEMonitorSetFromOptions(lme,"-lme_monitor","error_estimate",NULL));

    /* -----------------------------------------------------------------------*/
    PetscCall(PetscOptionsName("-lme_view","Print detailed information on solver used","LMEView",&set));

    PetscTryTypeMethod(lme,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)lme,PetscOptionsObject));
  PetscOptionsEnd();

  if (!lme->V) PetscCall(LMEGetBV(lme,&lme->V));
  PetscCall(BVSetFromOptions(lme->V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetProblemType - Specifies the type of matrix equation to be solved.

   Logically Collective

   Input Parameters:
+  lme  - the linear matrix equation solver context
-  type - a known type of matrix equation

   Options Database Keys:
+  -lme_lyapunov      - continuous-time Lyapunov equation
.  -lme_sylvester     - Sylvester equation
.  -lme_gen_lyapunov  - generalized Lyapunov equation
.  -lme_gen_sylvester - generalized Sylvester equation
.  -lme_dt_lyapunov   - discrete-time Lyapunov equation
-  -lme_stein         - Stein equation

   Notes:
   The coefficient matrices $A$, $B$, $D$, $E$ must be provided via `LMESetCoefficients()`,
   but some of them are optional depending on the matrix equation.

   Problem Type             | Equation         |`LMEProblemType`   | $A$ | $B$ | $D$ | $E$
   -------------------------|------------------|-------------------|-----|-----|-----|-----
   Continuous-Time Lyapunov | $AX+XA^*=-C$     |`LME_LYAPUNOV`     | yes |$A^*$|  -  |  -
   Sylvester                | $AX+XB=C$        |`LME_SYLVESTER`    | yes | yes |  -  |  -
   Generalized Lyapunov     | $AXD^*+DXA^*=-C$ |`LME_GEN_LYAPUNOV` | yes |$A^*$| yes |$D^*$
   Generalized Sylvester    | $AXE+DXB=C$      |`LME_GEN_SYLVESTER`| yes | yes | yes | yes
   Discrete-Time Lyapunov   | $AXA^*-X=-C$     |`LME_DT_LYAPUNOV`  | yes |  -  |  -  |$A^*$
   Stein                    | $AXE-X=-C$       |`LME_STEIN`        | yes |  -  |  -  | yes

   In the above table, the notation $A^*$ means that this matrix need
   not be passed, but the user may choose to pass an explicit transpose
   of matrix $A$ (for improved efficiency).

   Also note that some of the equation types impose restrictions on the
   properties of the coefficient matrices and possibly on the right-hand
   side $C$.

   Level: beginner

.seealso: [](ch:lme), `LMESetCoefficients()`, `LMESetType()`, `LMEGetProblemType()`, `LMEProblemType`
@*/
PetscErrorCode LMESetProblemType(LME lme,LMEProblemType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveEnum(lme,type,2);
  if (type == lme->problem_type) PetscFunctionReturn(PETSC_SUCCESS);
  switch (type) {
    case LME_LYAPUNOV:
    case LME_SYLVESTER:
    case LME_GEN_LYAPUNOV:
    case LME_GEN_SYLVESTER:
    case LME_DT_LYAPUNOV:
    case LME_STEIN:
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONG,"Unknown matrix equation type");
  }
  lme->problem_type = type;
  lme->setupcalled  = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetProblemType - Gets the matrix equation type from the `LME` object.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  type - the problem type

   Level: beginner

.seealso: [](ch:lme), `LMESetProblemType()`, `LMEProblemType`
@*/
PetscErrorCode LMEGetProblemType(LME lme,LMEProblemType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = lme->problem_type;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetTolerances - Gets the tolerance and maximum iteration count used
   by the `LME` convergence tests.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Note:
   The user can specify `NULL` for any parameter that is not needed.

   Level: intermediate

.seealso: [](ch:lme), `LMESetTolerances()`
@*/
PetscErrorCode LMEGetTolerances(LME lme,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (tol)    *tol    = lme->tol;
  if (maxits) *maxits = lme->max_it;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetTolerances - Sets the tolerance and maximum iteration count used
   by the `LME` convergence tests.

   Logically Collective

   Input Parameters:
+  lme - the linear matrix equation solver context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -lme_tol \<tol\>       - sets the convergence tolerance
-  -lme_max_it \<maxits\> - sets the maximum number of iterations allowed

   Notes:
   Use `PETSC_CURRENT` to retain the current value of any of the parameters.
   Use `PETSC_DETERMINE` for either argument to assign a default value computed
   internally (may be different in each solver).
   For `maxits` use `PETSC_UNLIMITED` to indicate there is no upper bound on this value.

   Level: intermediate

.seealso: [](ch:lme), `LMEGetTolerances()`
@*/
PetscErrorCode LMESetTolerances(LME lme,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveReal(lme,tol,2);
  PetscValidLogicalCollectiveInt(lme,maxits,3);
  if (tol == (PetscReal)PETSC_DETERMINE) {
    lme->tol = PETSC_DETERMINE;
    lme->setupcalled = 0;
  } else if (tol != (PetscReal)PETSC_CURRENT) {
    PetscCheck(tol>0.0,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    lme->tol = tol;
  }
  if (maxits == PETSC_DETERMINE) {
    lme->max_it = PETSC_DETERMINE;
    lme->setupcalled = 0;
  } else if (maxits == PETSC_UNLIMITED) {
    lme->max_it = PETSC_INT_MAX;
  } else if (maxits != PETSC_CURRENT) {
    PetscCheck(maxits>0,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    lme->max_it = maxits;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetDimensions - Gets the dimension of the subspace used by the solver.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  ncv - the maximum dimension of the subspace to be used by the solver

   Level: intermediate

.seealso: [](ch:lme), `LMESetDimensions()`
@*/
PetscErrorCode LMEGetDimensions(LME lme,PetscInt *ncv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(ncv,2);
  *ncv = lme->ncv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetDimensions - Sets the dimension of the subspace to be used by the solver.

   Logically Collective

   Input Parameters:
+  lme - the linear matrix equation solver context
-  ncv - the maximum dimension of the subspace to be used by the solver

   Options Database Key:
.  -lme_ncv \<ncv\> - sets the dimension of the subspace

   Notes:
   Use `PETSC_DETERMINE` for `ncv` to assign a reasonably good value, which is
   dependent on the solution method.

   Level: intermediate

.seealso: [](ch:lme), `LMEGetDimensions()`
@*/
PetscErrorCode LMESetDimensions(LME lme,PetscInt ncv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveInt(lme,ncv,2);
  if (ncv == PETSC_DECIDE || ncv == PETSC_DEFAULT) {
    lme->ncv = PETSC_DETERMINE;
  } else {
    PetscCheck(ncv>0,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    lme->ncv = ncv;
  }
  lme->setupcalled = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetErrorIfNotConverged - Causes `LMESolve()` to generate an error if the
   solver has not converged.

   Logically Collective

   Input Parameters:
+  lme - the linear matrix equation solver context
-  flg - `PETSC_TRUE` indicates you want the error generated

   Options Database Key:
.  -lme_error_if_not_converged - generate an error and stop the program

   Note:
   Normally SLEPc continues if the solver fails to converge, you can call
   `LMEGetConvergedReason()` after a `LMESolve()` to determine if it has converged.

   Level: intermediate

.seealso: [](ch:lme), `LMEGetErrorIfNotConverged()`
@*/
PetscErrorCode LMESetErrorIfNotConverged(LME lme,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveBool(lme,flg,2);
  lme->errorifnotconverged = flg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetErrorIfNotConverged - Return a flag indicating whether `LMESolve()` will
   generate an error if the solver does not converge.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  flag - `PETSC_TRUE` if it will generate an error, else `PETSC_FALSE`

   Level: intermediate

.seealso: [](ch:lme), `LMESetErrorIfNotConverged()`
@*/
PetscErrorCode LMEGetErrorIfNotConverged(LME lme,PetscBool *flag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(flag,2);
  *flag = lme->errorifnotconverged;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetOptionsPrefix - Sets the prefix used for searching for all
   `LME` options in the database.

   Logically Collective

   Input Parameters:
+  lme    - the linear matrix equation solver context
-  prefix - the prefix string to prepend to all `LME` option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different `LME` contexts, one could call
.vb
   LMESetOptionsPrefix(lme1,"fun1_")
   LMESetOptionsPrefix(lme2,"fun2_")
.ve

   Level: advanced

.seealso: [](ch:lme), `LMEAppendOptionsPrefix()`, `LMEGetOptionsPrefix()`
@*/
PetscErrorCode LMESetOptionsPrefix(LME lme,const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (!lme->V) PetscCall(LMEGetBV(lme,&lme->V));
  PetscCall(BVSetOptionsPrefix(lme->V,prefix));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)lme,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEAppendOptionsPrefix - Appends to the prefix used for searching for all
   `LME` options in the database.

   Logically Collective

   Input Parameters:
+  lme    - the linear matrix equation solver context
-  prefix - the prefix string to prepend to all `LME` option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: [](ch:lme), `LMESetOptionsPrefix()`, `LMEGetOptionsPrefix()`
@*/
PetscErrorCode LMEAppendOptionsPrefix(LME lme,const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (!lme->V) PetscCall(LMEGetBV(lme,&lme->V));
  PetscCall(BVAppendOptionsPrefix(lme->V,prefix));
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)lme,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetOptionsPrefix - Gets the prefix used for searching for all
   `LME` options in the database.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  prefix - pointer to the prefix string used is returned

   Level: advanced

.seealso: [](ch:lme), `LMESetOptionsPrefix()`, `LMEAppendOptionsPrefix()`
@*/
PetscErrorCode LMEGetOptionsPrefix(LME lme,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)lme,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}
