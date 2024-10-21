/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic NEP routines
*/

#include <slepc/private/nepimpl.h>      /*I "slepcnep.h" I*/

/* Logging support */
PetscClassId      NEP_CLASSID = 0;
PetscLogEvent     NEP_SetUp = 0,NEP_Solve = 0,NEP_Refine = 0,NEP_FunctionEval = 0,NEP_JacobianEval = 0,NEP_Resolvent = 0,NEP_CISS_SVD = 0;

/* List of registered NEP routines */
PetscFunctionList NEPList = NULL;
PetscBool         NEPRegisterAllCalled = PETSC_FALSE;

/* List of registered NEP monitors */
PetscFunctionList NEPMonitorList              = NULL;
PetscFunctionList NEPMonitorCreateList        = NULL;
PetscFunctionList NEPMonitorDestroyList       = NULL;
PetscBool         NEPMonitorRegisterAllCalled = PETSC_FALSE;

/*@
   NEPCreate - Creates the default NEP context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  outnep - location to put the NEP context

   Level: beginner

.seealso: NEPSetUp(), NEPSolve(), NEPDestroy(), NEP
@*/
PetscErrorCode NEPCreate(MPI_Comm comm,NEP *outnep)
{
  NEP            nep;

  PetscFunctionBegin;
  PetscAssertPointer(outnep,2);
  PetscCall(NEPInitializePackage());
  PetscCall(SlepcHeaderCreate(nep,NEP_CLASSID,"NEP","Nonlinear Eigenvalue Problem","NEP",comm,NEPDestroy,NEPView));

  nep->max_it          = PETSC_DETERMINE;
  nep->nev             = 1;
  nep->ncv             = PETSC_DETERMINE;
  nep->mpd             = PETSC_DETERMINE;
  nep->nini            = 0;
  nep->target          = 0.0;
  nep->tol             = PETSC_DETERMINE;
  nep->conv            = NEP_CONV_REL;
  nep->stop            = NEP_STOP_BASIC;
  nep->which           = (NEPWhich)0;
  nep->problem_type    = (NEPProblemType)0;
  nep->refine          = NEP_REFINE_NONE;
  nep->npart           = 1;
  nep->rtol            = PETSC_DETERMINE;
  nep->rits            = PETSC_DETERMINE;
  nep->scheme          = (NEPRefineScheme)0;
  nep->trackall        = PETSC_FALSE;
  nep->twosided        = PETSC_FALSE;

  nep->computefunction = NULL;
  nep->computejacobian = NULL;
  nep->functionctx     = NULL;
  nep->jacobianctx     = NULL;
  nep->converged       = NEPConvergedRelative;
  nep->convergeduser   = NULL;
  nep->convergeddestroy= NULL;
  nep->stopping        = NEPStoppingBasic;
  nep->stoppinguser    = NULL;
  nep->stoppingdestroy = NULL;
  nep->convergedctx    = NULL;
  nep->stoppingctx     = NULL;
  nep->numbermonitors  = 0;

  nep->ds              = NULL;
  nep->V               = NULL;
  nep->W               = NULL;
  nep->rg              = NULL;
  nep->function        = NULL;
  nep->function_pre    = NULL;
  nep->jacobian        = NULL;
  nep->A               = NULL;
  nep->f               = NULL;
  nep->nt              = 0;
  nep->mstr            = UNKNOWN_NONZERO_PATTERN;
  nep->P               = NULL;
  nep->mstrp           = UNKNOWN_NONZERO_PATTERN;
  nep->IS              = NULL;
  nep->eigr            = NULL;
  nep->eigi            = NULL;
  nep->errest          = NULL;
  nep->perm            = NULL;
  nep->nwork           = 0;
  nep->work            = NULL;
  nep->data            = NULL;

  nep->state           = NEP_STATE_INITIAL;
  nep->nconv           = 0;
  nep->its             = 0;
  nep->n               = 0;
  nep->nloc            = 0;
  nep->nrma            = NULL;
  nep->fui             = (NEPUserInterface)0;
  nep->useds           = PETSC_FALSE;
  nep->resolvent       = NULL;
  nep->reason          = NEP_CONVERGED_ITERATING;

  PetscCall(PetscNew(&nep->sc));
  *outnep = nep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetType - Selects the particular solver to be used in the NEP object.

   Logically Collective

   Input Parameters:
+  nep      - the nonlinear eigensolver context
-  type     - a known method

   Options Database Key:
.  -nep_type <method> - Sets the method; use -help for a list
    of available methods

   Notes:
   See "slepc/include/slepcnep.h" for available methods.

   Normally, it is best to use the NEPSetFromOptions() command and
   then set the NEP type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The NEPSetType() routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: NEPType
@*/
PetscErrorCode NEPSetType(NEP nep,NEPType type)
{
  PetscErrorCode (*r)(NEP);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)nep,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscFunctionListFind(NEPList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown NEP type given: %s",type);

  PetscTryTypeMethod(nep,destroy);
  PetscCall(PetscMemzero(nep->ops,sizeof(struct _NEPOps)));

  nep->state = NEP_STATE_INITIAL;
  PetscCall(PetscObjectChangeTypeName((PetscObject)nep,type));
  PetscCall((*r)(nep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetType - Gets the NEP type as a string from the NEP object.

   Not Collective

   Input Parameter:
.  nep - the eigensolver context

   Output Parameter:
.  type - name of NEP method

   Level: intermediate

.seealso: NEPSetType()
@*/
PetscErrorCode NEPGetType(NEP nep,NEPType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)nep)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPRegister - Adds a method to the nonlinear eigenproblem solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Notes:
   NEPRegister() may be called multiple times to add several user-defined solvers.

   Example Usage:
.vb
    NEPRegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     NEPSetType(nep,"my_solver")
   or at runtime via the option
$     -nep_type my_solver

   Level: advanced

.seealso: NEPRegisterAll()
@*/
PetscErrorCode NEPRegister(const char *name,PetscErrorCode (*function)(NEP))
{
  PetscFunctionBegin;
  PetscCall(NEPInitializePackage());
  PetscCall(PetscFunctionListAdd(&NEPList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorRegister - Adds NEP monitor routine.

   Not Collective

   Input Parameters:
+  name    - name of a new monitor routine
.  vtype   - a PetscViewerType for the output
.  format  - a PetscViewerFormat for the output
.  monitor - monitor routine
.  create  - creation routine, or NULL
-  destroy - destruction routine, or NULL

   Notes:
   NEPMonitorRegister() may be called multiple times to add several user-defined monitors.

   Example Usage:
.vb
   NEPMonitorRegister("my_monitor",PETSCVIEWERASCII,PETSC_VIEWER_ASCII_INFO_DETAIL,MyMonitor,NULL,NULL);
.ve

   Then, your monitor can be chosen with the procedural interface via
$      NEPMonitorSetFromOptions(nep,"-nep_monitor_my_monitor","my_monitor",NULL)
   or at runtime via the option
$      -nep_monitor_my_monitor

   Level: advanced

.seealso: NEPMonitorRegisterAll()
@*/
PetscErrorCode NEPMonitorRegister(const char name[],PetscViewerType vtype,PetscViewerFormat format,PetscErrorCode (*monitor)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscErrorCode (*create)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**),PetscErrorCode (*destroy)(PetscViewerAndFormat**))
{
  char           key[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  PetscCall(NEPInitializePackage());
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListAdd(&NEPMonitorList,key,monitor));
  if (create)  PetscCall(PetscFunctionListAdd(&NEPMonitorCreateList,key,create));
  if (destroy) PetscCall(PetscFunctionListAdd(&NEPMonitorDestroyList,key,destroy));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   NEPReset_Problem - Destroys the problem matrices.
*/
PetscErrorCode NEPReset_Problem(NEP nep)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscCall(MatDestroy(&nep->function));
  PetscCall(MatDestroy(&nep->function_pre));
  PetscCall(MatDestroy(&nep->jacobian));
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    PetscCall(MatDestroyMatrices(nep->nt,&nep->A));
    for (i=0;i<nep->nt;i++) PetscCall(FNDestroy(&nep->f[i]));
    PetscCall(PetscFree(nep->f));
    PetscCall(PetscFree(nep->nrma));
    if (nep->P) PetscCall(MatDestroyMatrices(nep->nt,&nep->P));
    nep->nt = 0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*@
   NEPReset - Resets the NEP context to the initial state (prior to setup)
   and destroys any allocated Vecs and Mats.

   Collective

   Input Parameter:
.  nep - eigensolver context obtained from NEPCreate()

   Level: advanced

.seealso: NEPDestroy()
@*/
PetscErrorCode NEPReset(NEP nep)
{
  PetscFunctionBegin;
  if (nep) PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!nep) PetscFunctionReturn(PETSC_SUCCESS);
  PetscTryTypeMethod(nep,reset);
  if (nep->refineksp) PetscCall(KSPReset(nep->refineksp));
  PetscCall(NEPReset_Problem(nep));
  PetscCall(BVDestroy(&nep->V));
  PetscCall(BVDestroy(&nep->W));
  PetscCall(VecDestroyVecs(nep->nwork,&nep->work));
  PetscCall(MatDestroy(&nep->resolvent));
  nep->nwork = 0;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPDestroy - Destroys the NEP context.

   Collective

   Input Parameter:
.  nep - eigensolver context obtained from NEPCreate()

   Level: beginner

.seealso: NEPCreate(), NEPSetUp(), NEPSolve()
@*/
PetscErrorCode NEPDestroy(NEP *nep)
{
  PetscFunctionBegin;
  if (!*nep) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*nep,NEP_CLASSID,1);
  if (--((PetscObject)*nep)->refct > 0) { *nep = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscCall(NEPReset(*nep));
  PetscTryTypeMethod(*nep,destroy);
  if ((*nep)->eigr) PetscCall(PetscFree4((*nep)->eigr,(*nep)->eigi,(*nep)->errest,(*nep)->perm));
  PetscCall(RGDestroy(&(*nep)->rg));
  PetscCall(DSDestroy(&(*nep)->ds));
  PetscCall(KSPDestroy(&(*nep)->refineksp));
  PetscCall(PetscSubcommDestroy(&(*nep)->refinesubc));
  PetscCall(PetscFree((*nep)->sc));
  /* just in case the initial vectors have not been used */
  PetscCall(SlepcBasisDestroy_Private(&(*nep)->nini,&(*nep)->IS));
  if ((*nep)->convergeddestroy) PetscCall((*(*nep)->convergeddestroy)(&(*nep)->convergedctx));
  if ((*nep)->stoppingdestroy) PetscCall((*(*nep)->stoppingdestroy)(&(*nep)->stoppingctx));
  PetscCall(NEPMonitorCancel(*nep));
  PetscCall(PetscHeaderDestroy(nep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetBV - Associates a basis vectors object to the nonlinear eigensolver.

   Collective

   Input Parameters:
+  nep - eigensolver context obtained from NEPCreate()
-  bv  - the basis vectors object

   Note:
   Use NEPGetBV() to retrieve the basis vectors context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: NEPGetBV()
@*/
PetscErrorCode NEPSetBV(NEP nep,BV bv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(bv,BV_CLASSID,2);
  PetscCheckSameComm(nep,1,bv,2);
  PetscCall(PetscObjectReference((PetscObject)bv));
  PetscCall(BVDestroy(&nep->V));
  nep->V = bv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetBV - Obtain the basis vectors object associated to the nonlinear
   eigensolver object.

   Not Collective

   Input Parameters:
.  nep - eigensolver context obtained from NEPCreate()

   Output Parameter:
.  bv - basis vectors context

   Level: advanced

.seealso: NEPSetBV()
@*/
PetscErrorCode NEPGetBV(NEP nep,BV *bv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(bv,2);
  if (!nep->V) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)nep),&nep->V));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)nep->V,(PetscObject)nep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)nep->V,((PetscObject)nep)->options));
  }
  *bv = nep->V;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetRG - Associates a region object to the nonlinear eigensolver.

   Collective

   Input Parameters:
+  nep - eigensolver context obtained from NEPCreate()
-  rg  - the region object

   Note:
   Use NEPGetRG() to retrieve the region context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: NEPGetRG()
@*/
PetscErrorCode NEPSetRG(NEP nep,RG rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (rg) {
    PetscValidHeaderSpecific(rg,RG_CLASSID,2);
    PetscCheckSameComm(nep,1,rg,2);
  }
  PetscCall(PetscObjectReference((PetscObject)rg));
  PetscCall(RGDestroy(&nep->rg));
  nep->rg = rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetRG - Obtain the region object associated to the
   nonlinear eigensolver object.

   Not Collective

   Input Parameters:
.  nep - eigensolver context obtained from NEPCreate()

   Output Parameter:
.  rg - region context

   Level: advanced

.seealso: NEPSetRG()
@*/
PetscErrorCode NEPGetRG(NEP nep,RG *rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(rg,2);
  if (!nep->rg) {
    PetscCall(RGCreate(PetscObjectComm((PetscObject)nep),&nep->rg));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)nep->rg,(PetscObject)nep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)nep->rg,((PetscObject)nep)->options));
  }
  *rg = nep->rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetDS - Associates a direct solver object to the nonlinear eigensolver.

   Collective

   Input Parameters:
+  nep - eigensolver context obtained from NEPCreate()
-  ds  - the direct solver object

   Note:
   Use NEPGetDS() to retrieve the direct solver context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: NEPGetDS()
@*/
PetscErrorCode NEPSetDS(NEP nep,DS ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(ds,DS_CLASSID,2);
  PetscCheckSameComm(nep,1,ds,2);
  PetscCall(PetscObjectReference((PetscObject)ds));
  PetscCall(DSDestroy(&nep->ds));
  nep->ds = ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetDS - Obtain the direct solver object associated to the
   nonlinear eigensolver object.

   Not Collective

   Input Parameters:
.  nep - eigensolver context obtained from NEPCreate()

   Output Parameter:
.  ds - direct solver context

   Level: advanced

.seealso: NEPSetDS()
@*/
PetscErrorCode NEPGetDS(NEP nep,DS *ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(ds,2);
  if (!nep->ds) {
    PetscCall(DSCreate(PetscObjectComm((PetscObject)nep),&nep->ds));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)nep->ds,(PetscObject)nep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)nep->ds,((PetscObject)nep)->options));
  }
  *ds = nep->ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPRefineGetKSP - Obtain the ksp object used by the eigensolver
   object in the refinement phase.

   Collective

   Input Parameters:
.  nep - eigensolver context obtained from NEPCreate()

   Output Parameter:
.  ksp - ksp context

   Level: advanced

.seealso: NEPSetRefine()
@*/
PetscErrorCode NEPRefineGetKSP(NEP nep,KSP *ksp)
{
  MPI_Comm       comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(ksp,2);
  if (!nep->refineksp) {
    if (nep->npart>1) {
      /* Split in subcomunicators */
      PetscCall(PetscSubcommCreate(PetscObjectComm((PetscObject)nep),&nep->refinesubc));
      PetscCall(PetscSubcommSetNumber(nep->refinesubc,nep->npart));
      PetscCall(PetscSubcommSetType(nep->refinesubc,PETSC_SUBCOMM_CONTIGUOUS));
      PetscCall(PetscSubcommGetChild(nep->refinesubc,&comm));
    } else PetscCall(PetscObjectGetComm((PetscObject)nep,&comm));
    PetscCall(KSPCreate(comm,&nep->refineksp));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)nep->refineksp,(PetscObject)nep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)nep->refineksp,((PetscObject)nep)->options));
    PetscCall(KSPSetOptionsPrefix(*ksp,((PetscObject)nep)->prefix));
    PetscCall(KSPAppendOptionsPrefix(*ksp,"nep_refine_"));
    PetscCall(KSPSetTolerances(nep->refineksp,SlepcDefaultTol(nep->rtol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
  }
  *ksp = nep->refineksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetTarget - Sets the value of the target.

   Logically Collective

   Input Parameters:
+  nep    - eigensolver context
-  target - the value of the target

   Options Database Key:
.  -nep_target <scalar> - the value of the target

   Notes:
   The target is a scalar value used to determine the portion of the spectrum
   of interest. It is used in combination with NEPSetWhichEigenpairs().

   In the case of complex scalars, a complex value can be provided in the
   command line with [+/-][realnumber][+/-]realnumberi with no spaces, e.g.
   -nep_target 1.0+2.0i

   Level: intermediate

.seealso: NEPGetTarget(), NEPSetWhichEigenpairs()
@*/
PetscErrorCode NEPSetTarget(NEP nep,PetscScalar target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,target,2);
  nep->target = target;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetTarget - Gets the value of the target.

   Not Collective

   Input Parameter:
.  nep - eigensolver context

   Output Parameter:
.  target - the value of the target

   Note:
   If the target was not set by the user, then zero is returned.

   Level: intermediate

.seealso: NEPSetTarget()
@*/
PetscErrorCode NEPGetTarget(NEP nep,PetscScalar* target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(target,2);
  *target = nep->target;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPSetFunction - Sets the function to compute the nonlinear Function T(lambda)
   as well as the location to store the matrix.

   Collective

   Input Parameters:
+  nep - the NEP context
.  A   - Function matrix
.  B   - preconditioner matrix (usually same as A)
.  fun - Function evaluation routine (if NULL then NEP retains any
         previously set value), see NEPFunctionFn for the calling sequence
-  ctx - [optional] user-defined context for private data for the Function
         evaluation routine (may be NULL) (if NULL then NEP retains any
         previously set value)

   Level: beginner

.seealso: NEPGetFunction(), NEPSetJacobian()
@*/
PetscErrorCode NEPSetFunction(NEP nep,Mat A,Mat B,NEPFunctionFn *fun,void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  if (B) PetscValidHeaderSpecific(B,MAT_CLASSID,3);
  if (A) PetscCheckSameComm(nep,1,A,2);
  if (B) PetscCheckSameComm(nep,1,B,3);

  if (nep->state) PetscCall(NEPReset(nep));
  else if (nep->fui && nep->fui!=NEP_USER_INTERFACE_CALLBACK) PetscCall(NEPReset_Problem(nep));

  if (fun) nep->computefunction = fun;
  if (ctx) nep->functionctx     = ctx;
  if (A) {
    PetscCall(PetscObjectReference((PetscObject)A));
    PetscCall(MatDestroy(&nep->function));
    nep->function = A;
  }
  if (B) {
    PetscCall(PetscObjectReference((PetscObject)B));
    PetscCall(MatDestroy(&nep->function_pre));
    nep->function_pre = B;
  }
  nep->fui   = NEP_USER_INTERFACE_CALLBACK;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPGetFunction - Returns the Function matrix and optionally the user
   provided context for evaluating the Function.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  A   - location to stash Function matrix (or NULL)
.  B   - location to stash preconditioner matrix (or NULL)
.  fun - location to put Function function (or NULL)
-  ctx - location to stash Function context (or NULL)

   Level: advanced

.seealso: NEPSetFunction()
@*/
PetscErrorCode NEPGetFunction(NEP nep,Mat *A,Mat *B,NEPFunctionFn **fun,void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckCallback(nep,1);
  if (A)   *A   = nep->function;
  if (B)   *B   = nep->function_pre;
  if (fun) *fun = nep->computefunction;
  if (ctx) *ctx = nep->functionctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPSetJacobian - Sets the function to compute the Jacobian T'(lambda) as well
   as the location to store the matrix.

   Collective

   Input Parameters:
+  nep - the NEP context
.  A   - Jacobian matrix
.  jac - Jacobian evaluation routine (if NULL then NEP retains any
         previously set value), see NEPJacobianFn for the calling sequence
-  ctx - [optional] user-defined context for private data for the Jacobian
         evaluation routine (may be NULL) (if NULL then NEP retains any
         previously set value)

   Level: beginner

.seealso: NEPSetFunction(), NEPGetJacobian()
@*/
PetscErrorCode NEPSetJacobian(NEP nep,Mat A,NEPJacobianFn *jac,void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  if (A) PetscCheckSameComm(nep,1,A,2);

  if (nep->state) PetscCall(NEPReset(nep));
  else if (nep->fui && nep->fui!=NEP_USER_INTERFACE_CALLBACK) PetscCall(NEPReset_Problem(nep));

  if (jac) nep->computejacobian = jac;
  if (ctx) nep->jacobianctx     = ctx;
  if (A) {
    PetscCall(PetscObjectReference((PetscObject)A));
    PetscCall(MatDestroy(&nep->jacobian));
    nep->jacobian = A;
  }
  nep->fui   = NEP_USER_INTERFACE_CALLBACK;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPGetJacobian - Returns the Jacobian matrix and optionally the user
   provided routine and context for evaluating the Jacobian.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  A   - location to stash Jacobian matrix (or NULL)
.  jac - location to put Jacobian function (or NULL)
-  ctx - location to stash Jacobian context (or NULL)

   Level: advanced

.seealso: NEPSetJacobian()
@*/
PetscErrorCode NEPGetJacobian(NEP nep,Mat *A,NEPJacobianFn **jac,void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckCallback(nep,1);
  if (A)   *A   = nep->jacobian;
  if (jac) *jac = nep->computejacobian;
  if (ctx) *ctx = nep->jacobianctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetSplitOperator - Sets the operator of the nonlinear eigenvalue problem
   in split form.

   Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
.  nt  - number of terms in the split form
.  A   - array of matrices
.  f   - array of functions
-  str - structure flag for matrices

   Notes:
   The nonlinear operator is written as T(lambda) = sum_i A_i*f_i(lambda),
   for i=1,...,n. The derivative T'(lambda) can be obtained using the
   derivatives of f_i.

   The structure flag provides information about A_i's nonzero pattern
   (see MatStructure enum). If all matrices have the same pattern, then
   use SAME_NONZERO_PATTERN. If the patterns are different but contained
   in the pattern of the first one, then use SUBSET_NONZERO_PATTERN. If
   patterns are known to be different, use DIFFERENT_NONZERO_PATTERN.
   If set to UNKNOWN_NONZERO_PATTERN, the patterns will be compared to
   determine if they are equal.

   This function must be called before NEPSetUp(). If it is called again
   after NEPSetUp() then the NEP object is reset.

   Level: beginner

.seealso: NEPGetSplitOperatorTerm(), NEPGetSplitOperatorInfo(), NEPSetSplitPreconditioner()
@*/
PetscErrorCode NEPSetSplitOperator(NEP nep,PetscInt nt,Mat A[],FN f[],MatStructure str)
{
  PetscInt       i,n=0,m,m0=0,mloc,nloc,mloc0=0;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,nt,2);
  PetscCheck(nt>0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Must have one or more terms, you have %" PetscInt_FMT,nt);
  PetscAssertPointer(A,3);
  PetscAssertPointer(f,4);
  PetscValidLogicalCollectiveEnum(nep,str,5);

  for (i=0;i<nt;i++) {
    PetscValidHeaderSpecific(A[i],MAT_CLASSID,3);
    PetscCheckSameComm(nep,1,A[i],3);
    PetscValidHeaderSpecific(f[i],FN_CLASSID,4);
    PetscCheckSameComm(nep,1,f[i],4);
    PetscCall(MatGetSize(A[i],&m,&n));
    PetscCall(MatGetLocalSize(A[i],&mloc,&nloc));
    PetscCheck(m==n,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"A[%" PetscInt_FMT "] is a non-square matrix (%" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",i,m,n);
    PetscCheck(mloc==nloc,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"A[%" PetscInt_FMT "] does not have equal row and column local sizes (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,mloc,nloc);
    if (!i) { m0 = m; mloc0 = mloc; }
    PetscCheck(m==m0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_INCOMP,"Dimensions of A[%" PetscInt_FMT "] do not match with previous matrices (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,m,m0);
    PetscCheck(mloc==mloc0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_INCOMP,"Local dimensions of A[%" PetscInt_FMT "] do not match with previous matrices (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,mloc,mloc0);
    PetscCall(PetscObjectReference((PetscObject)A[i]));
    PetscCall(PetscObjectReference((PetscObject)f[i]));
  }

  if (nep->state && (n!=nep->n || nloc!=nep->nloc)) PetscCall(NEPReset(nep));
  else PetscCall(NEPReset_Problem(nep));

  /* allocate space and copy matrices and functions */
  PetscCall(PetscMalloc1(nt,&nep->A));
  for (i=0;i<nt;i++) nep->A[i] = A[i];
  PetscCall(PetscMalloc1(nt,&nep->f));
  for (i=0;i<nt;i++) nep->f[i] = f[i];
  PetscCall(PetscCalloc1(nt,&nep->nrma));
  nep->nt    = nt;
  nep->mstr  = str;
  nep->fui   = NEP_USER_INTERFACE_SPLIT;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetSplitOperatorTerm - Gets the matrices and functions associated with
   the nonlinear operator in split form.

   Not Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  k   - the index of the requested term (starting in 0)

   Output Parameters:
+  A - the matrix of the requested term
-  f - the function of the requested term

   Level: intermediate

.seealso: NEPSetSplitOperator(), NEPGetSplitOperatorInfo()
@*/
PetscErrorCode NEPGetSplitOperatorTerm(NEP nep,PetscInt k,Mat *A,FN *f)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,k,2);
  NEPCheckSplit(nep,1);
  PetscCheck(k>=0 && k<nep->nt,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,nep->nt-1);
  if (A) *A = nep->A[k];
  if (f) *f = nep->f[k];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetSplitOperatorInfo - Returns the number of terms of the split form of
   the nonlinear operator, as well as the structure flag for matrices.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  n   - the number of terms passed in NEPSetSplitOperator()
-  str - the matrix structure flag passed in NEPSetSplitOperator()

   Level: intermediate

.seealso: NEPSetSplitOperator(), NEPGetSplitOperatorTerm()
@*/
PetscErrorCode NEPGetSplitOperatorInfo(NEP nep,PetscInt *n,MatStructure *str)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckSplit(nep,1);
  if (n)   *n = nep->nt;
  if (str) *str = nep->mstr;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetSplitPreconditioner - Sets an operator in split form from which
   to build the preconditioner to be used when solving the nonlinear
   eigenvalue problem in split form.

   Collective

   Input Parameters:
+  nep  - the nonlinear eigensolver context
.  ntp  - number of terms in the split preconditioner
.  P    - array of matrices
-  strp - structure flag for matrices

   Notes:
   The matrix for the preconditioner is expressed as P(lambda) =
   sum_i P_i*f_i(lambda), for i=1,...,n, where the f_i functions
   are the same as in NEPSetSplitOperator(). It is not necessary to call
   this function. If it is not invoked, then the preconditioner is
   built from T(lambda), i.e., both matrices and functions passed in
   NEPSetSplitOperator().

   The structure flag provides information about P_i's nonzero pattern
   in the same way as in NEPSetSplitOperator().

   If the functions defining the preconditioner operator were different
   from the ones given in NEPSetSplitOperator(), then the split form
   cannot be used. Use the callback interface instead.

   Use ntp=0 to reset a previously set split preconditioner.

   Level: advanced

.seealso: NEPGetSplitPreconditionerTerm(), NEPGetSplitPreconditionerInfo(), NEPSetSplitOperator()
@*/
PetscErrorCode NEPSetSplitPreconditioner(NEP nep,PetscInt ntp,Mat P[],MatStructure strp)
{
  PetscInt       i,n=0,m,m0=0,mloc,nloc,mloc0=0;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,ntp,2);
  PetscCheck(ntp>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Negative value of ntp = %" PetscInt_FMT,ntp);
  PetscCheck(nep->fui==NEP_USER_INTERFACE_SPLIT,PetscObjectComm((PetscObject)nep),PETSC_ERR_ORDER,"Must call NEPSetSplitOperator first");
  PetscCheck(ntp==0 || nep->nt==ntp,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The number of terms must be the same as in NEPSetSplitOperator()");
  if (ntp) PetscAssertPointer(P,3);
  PetscValidLogicalCollectiveEnum(nep,strp,4);

  for (i=0;i<ntp;i++) {
    PetscValidHeaderSpecific(P[i],MAT_CLASSID,3);
    PetscCheckSameComm(nep,1,P[i],3);
    PetscCall(MatGetSize(P[i],&m,&n));
    PetscCall(MatGetLocalSize(P[i],&mloc,&nloc));
    PetscCheck(m==n,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"P[%" PetscInt_FMT "] is a non-square matrix (%" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",i,m,n);
    PetscCheck(mloc==nloc,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"P[%" PetscInt_FMT "] does not have equal row and column local sizes (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,mloc,nloc);
    if (!i) { m0 = m; mloc0 = mloc; }
    PetscCheck(m==m0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_INCOMP,"Dimensions of P[%" PetscInt_FMT "] do not match with previous matrices (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,m,m0);
    PetscCheck(mloc==mloc0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_INCOMP,"Local dimensions of P[%" PetscInt_FMT "] do not match with previous matrices (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,mloc,mloc0);
    PetscCall(PetscObjectReference((PetscObject)P[i]));
  }

  PetscCheck(!nep->state,PetscObjectComm((PetscObject)nep),PETSC_ERR_ORDER,"To call this function after NEPSetUp(), you must call NEPSetSplitOperator() again");
  if (nep->P) PetscCall(MatDestroyMatrices(nep->nt,&nep->P));

  /* allocate space and copy matrices */
  if (ntp) {
    PetscCall(PetscMalloc1(ntp,&nep->P));
    for (i=0;i<ntp;i++) nep->P[i] = P[i];
  }
  nep->mstrp = strp;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetSplitPreconditionerTerm - Gets the matrices associated with
   the split preconditioner.

   Not Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  k   - the index of the requested term (starting in 0)

   Output Parameter:
.  P  - the matrix of the requested term

   Level: advanced

.seealso: NEPSetSplitPreconditioner(), NEPGetSplitPreconditionerInfo()
@*/
PetscErrorCode NEPGetSplitPreconditionerTerm(NEP nep,PetscInt k,Mat *P)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,k,2);
  PetscAssertPointer(P,3);
  NEPCheckSplit(nep,1);
  PetscCheck(k>=0 && k<nep->nt,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,nep->nt-1);
  PetscCheck(nep->P,PetscObjectComm((PetscObject)nep),PETSC_ERR_ORDER,"You have not called NEPSetSplitPreconditioner()");
  *P = nep->P[k];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetSplitPreconditionerInfo - Returns the number of terms of the split
   preconditioner, as well as the structure flag for matrices.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  n    - the number of terms passed in NEPSetSplitPreconditioner()
-  strp - the matrix structure flag passed in NEPSetSplitPreconditioner()

   Level: advanced

.seealso: NEPSetSplitPreconditioner(), NEPGetSplitPreconditionerTerm()
@*/
PetscErrorCode NEPGetSplitPreconditionerInfo(NEP nep,PetscInt *n,MatStructure *strp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckSplit(nep,1);
  if (n)    *n    = nep->P? nep->nt: 0;
  if (strp) *strp = nep->mstrp;
  PetscFunctionReturn(PETSC_SUCCESS);
}
