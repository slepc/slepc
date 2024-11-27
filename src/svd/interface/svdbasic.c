/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic SVD routines
*/

#include <slepc/private/svdimpl.h>      /*I "slepcsvd.h" I*/

/* Logging support */
PetscClassId      SVD_CLASSID = 0;
PetscLogEvent     SVD_SetUp = 0,SVD_Solve = 0;

/* List of registered SVD routines */
PetscFunctionList SVDList = NULL;
PetscBool         SVDRegisterAllCalled = PETSC_FALSE;

/* List of registered SVD monitors */
PetscFunctionList SVDMonitorList              = NULL;
PetscFunctionList SVDMonitorCreateList        = NULL;
PetscFunctionList SVDMonitorDestroyList       = NULL;
PetscBool         SVDMonitorRegisterAllCalled = PETSC_FALSE;

/*@
   SVDCreate - Creates the default SVD context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  outsvd - location to put the SVD context

   Note:
   The default SVD type is SVDCROSS

   Level: beginner

.seealso: SVDSetUp(), SVDSolve(), SVDDestroy(), SVD
@*/
PetscErrorCode SVDCreate(MPI_Comm comm,SVD *outsvd)
{
  SVD            svd;

  PetscFunctionBegin;
  PetscAssertPointer(outsvd,2);
  PetscCall(SVDInitializePackage());
  PetscCall(SlepcHeaderCreate(svd,SVD_CLASSID,"SVD","Singular Value Decomposition","SVD",comm,SVDDestroy,SVDView));

  svd->OP               = NULL;
  svd->OPb              = NULL;
  svd->omega            = NULL;
  svd->max_it           = PETSC_DETERMINE;
  svd->nsv              = 0;
  svd->ncv              = PETSC_DETERMINE;
  svd->mpd              = PETSC_DETERMINE;
  svd->nini             = 0;
  svd->ninil            = 0;
  svd->tol              = PETSC_DETERMINE;
  svd->thres            = 0.0;
  svd->threlative       = PETSC_FALSE;
  svd->conv             = (SVDConv)-1;
  svd->stop             = SVD_STOP_BASIC;
  svd->which            = SVD_LARGEST;
  svd->problem_type     = (SVDProblemType)0;
  svd->impltrans        = PETSC_FALSE;
  svd->trackall         = PETSC_FALSE;

  svd->converged        = NULL;
  svd->convergeduser    = NULL;
  svd->convergeddestroy = NULL;
  svd->stopping         = SVDStoppingBasic;
  svd->stoppinguser     = NULL;
  svd->stoppingdestroy  = NULL;
  svd->convergedctx     = NULL;
  svd->stoppingctx      = NULL;
  svd->numbermonitors   = 0;

  svd->ds               = NULL;
  svd->U                = NULL;
  svd->V                = NULL;
  svd->A                = NULL;
  svd->B                = NULL;
  svd->AT               = NULL;
  svd->BT               = NULL;
  svd->IS               = NULL;
  svd->ISL              = NULL;
  svd->sigma            = NULL;
  svd->errest           = NULL;
  svd->sign             = NULL;
  svd->perm             = NULL;
  svd->nworkl           = 0;
  svd->nworkr           = 0;
  svd->workl            = NULL;
  svd->workr            = NULL;
  svd->data             = NULL;

  svd->state            = SVD_STATE_INITIAL;
  svd->nconv            = 0;
  svd->its              = 0;
  svd->leftbasis        = PETSC_FALSE;
  svd->swapped          = PETSC_FALSE;
  svd->expltrans        = PETSC_FALSE;
  svd->nrma             = 0.0;
  svd->nrmb             = 0.0;
  svd->isgeneralized    = PETSC_FALSE;
  svd->reason           = SVD_CONVERGED_ITERATING;

  PetscCall(PetscNew(&svd->sc));
  *outsvd = svd;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDReset - Resets the SVD context to the initial state (prior to setup)
   and destroys any allocated Vecs and Mats.

   Collective

   Input Parameter:
.  svd - singular value solver context obtained from SVDCreate()

   Level: advanced

.seealso: SVDDestroy()
@*/
PetscErrorCode SVDReset(SVD svd)
{
  PetscFunctionBegin;
  if (svd) PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!svd) PetscFunctionReturn(PETSC_SUCCESS);
  PetscTryTypeMethod(svd,reset);
  PetscCall(MatDestroy(&svd->OP));
  PetscCall(MatDestroy(&svd->OPb));
  PetscCall(VecDestroy(&svd->omega));
  PetscCall(MatDestroy(&svd->A));
  PetscCall(MatDestroy(&svd->B));
  PetscCall(MatDestroy(&svd->AT));
  PetscCall(MatDestroy(&svd->BT));
  PetscCall(BVDestroy(&svd->U));
  PetscCall(BVDestroy(&svd->V));
  PetscCall(VecDestroyVecs(svd->nworkl,&svd->workl));
  svd->nworkl = 0;
  PetscCall(VecDestroyVecs(svd->nworkr,&svd->workr));
  svd->nworkr = 0;
  svd->swapped = PETSC_FALSE;
  svd->state = SVD_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDDestroy - Destroys the SVD context.

   Collective

   Input Parameter:
.  svd - singular value solver context obtained from SVDCreate()

   Level: beginner

.seealso: SVDCreate(), SVDSetUp(), SVDSolve()
@*/
PetscErrorCode SVDDestroy(SVD *svd)
{
  PetscFunctionBegin;
  if (!*svd) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*svd,SVD_CLASSID,1);
  if (--((PetscObject)*svd)->refct > 0) { *svd = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscCall(SVDReset(*svd));
  PetscTryTypeMethod(*svd,destroy);
  if ((*svd)->sigma) PetscCall(PetscFree3((*svd)->sigma,(*svd)->perm,(*svd)->errest));
  if ((*svd)->sign) PetscCall(PetscFree((*svd)->sign));
  PetscCall(DSDestroy(&(*svd)->ds));
  PetscCall(PetscFree((*svd)->sc));
  /* just in case the initial vectors have not been used */
  PetscCall(SlepcBasisDestroy_Private(&(*svd)->nini,&(*svd)->IS));
  PetscCall(SlepcBasisDestroy_Private(&(*svd)->ninil,&(*svd)->ISL));
  if ((*svd)->convergeddestroy) PetscCall((*(*svd)->convergeddestroy)(&(*svd)->convergedctx));
  if ((*svd)->stoppingdestroy) PetscCall((*(*svd)->stoppingdestroy)(&(*svd)->stoppingctx));
  PetscCall(SVDMonitorCancel(*svd));
  PetscCall(PetscHeaderDestroy(svd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetType - Selects the particular solver to be used in the SVD object.

   Logically Collective

   Input Parameters:
+  svd      - the singular value solver context
-  type     - a known method

   Options Database Key:
.  -svd_type <method> - Sets the method; use -help for a list
    of available methods

   Notes:
   See "slepc/include/slepcsvd.h" for available methods. The default
   is SVDCROSS.

   Normally, it is best to use the SVDSetFromOptions() command and
   then set the SVD type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The SVDSetType() routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: SVDType
@*/
PetscErrorCode SVDSetType(SVD svd,SVDType type)
{
  PetscErrorCode (*r)(SVD);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)svd,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscFunctionListFind(SVDList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown SVD type given: %s",type);

  PetscTryTypeMethod(svd,destroy);
  PetscCall(PetscMemzero(svd->ops,sizeof(struct _SVDOps)));

  svd->state = SVD_STATE_INITIAL;
  PetscCall(PetscObjectChangeTypeName((PetscObject)svd,type));
  PetscCall((*r)(svd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetType - Gets the SVD type as a string from the SVD object.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  type - name of SVD method

   Level: intermediate

.seealso: SVDSetType()
@*/
PetscErrorCode SVDGetType(SVD svd,SVDType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)svd)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDRegister - Adds a method to the singular value solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Notes:
   SVDRegister() may be called multiple times to add several user-defined solvers.

   Example Usage:
.vb
    SVDRegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     SVDSetType(svd,"my_solver")
   or at runtime via the option
$     -svd_type my_solver

   Level: advanced

.seealso: SVDRegisterAll()
@*/
PetscErrorCode SVDRegister(const char *name,PetscErrorCode (*function)(SVD))
{
  PetscFunctionBegin;
  PetscCall(SVDInitializePackage());
  PetscCall(PetscFunctionListAdd(&SVDList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorRegister - Adds SVD monitor routine.

   Not Collective

   Input Parameters:
+  name    - name of a new monitor routine
.  vtype   - a PetscViewerType for the output
.  format  - a PetscViewerFormat for the output
.  monitor - monitor routine
.  create  - creation routine, or NULL
-  destroy - destruction routine, or NULL

   Notes:
   SVDMonitorRegister() may be called multiple times to add several user-defined monitors.

   Example Usage:
.vb
   SVDMonitorRegister("my_monitor",PETSCVIEWERASCII,PETSC_VIEWER_ASCII_INFO_DETAIL,MyMonitor,NULL,NULL);
.ve

   Then, your monitor can be chosen with the procedural interface via
$      SVDMonitorSetFromOptions(svd,"-svd_monitor_my_monitor","my_monitor",NULL)
   or at runtime via the option
$      -svd_monitor_my_monitor

   Level: advanced

.seealso: SVDMonitorRegisterAll()
@*/
PetscErrorCode SVDMonitorRegister(const char name[],PetscViewerType vtype,PetscViewerFormat format,PetscErrorCode (*monitor)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscErrorCode (*create)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**),PetscErrorCode (*destroy)(PetscViewerAndFormat**))
{
  char           key[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  PetscCall(SVDInitializePackage());
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListAdd(&SVDMonitorList,key,monitor));
  if (create)  PetscCall(PetscFunctionListAdd(&SVDMonitorCreateList,key,create));
  if (destroy) PetscCall(PetscFunctionListAdd(&SVDMonitorDestroyList,key,destroy));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetBV - Associates basis vectors objects to the singular value solver.

   Collective

   Input Parameters:
+  svd - singular value solver context obtained from SVDCreate()
.  V   - the basis vectors object for right singular vectors
-  U   - the basis vectors object for left singular vectors

   Note:
   Use SVDGetBV() to retrieve the basis vectors contexts (for example,
   to free them at the end of the computations).

   Level: advanced

.seealso: SVDGetBV()
@*/
PetscErrorCode SVDSetBV(SVD svd,BV V,BV U)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (V) {
    PetscValidHeaderSpecific(V,BV_CLASSID,2);
    PetscCheckSameComm(svd,1,V,2);
    PetscCall(PetscObjectReference((PetscObject)V));
    PetscCall(BVDestroy(&svd->V));
    svd->V = V;
  }
  if (U) {
    PetscValidHeaderSpecific(U,BV_CLASSID,3);
    PetscCheckSameComm(svd,1,U,3);
    PetscCall(PetscObjectReference((PetscObject)U));
    PetscCall(BVDestroy(&svd->U));
    svd->U = U;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetBV - Obtain the basis vectors objects associated to the singular
   value solver object.

   Not Collective

   Input Parameter:
.  svd - singular value solver context obtained from SVDCreate()

   Output Parameters:
+  V - basis vectors context for right singular vectors
-  U - basis vectors context for left singular vectors

   Level: advanced

.seealso: SVDSetBV()
@*/
PetscErrorCode SVDGetBV(SVD svd,BV *V,BV *U)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (V) {
    if (!svd->V) {
      PetscCall(BVCreate(PetscObjectComm((PetscObject)svd),&svd->V));
      PetscCall(PetscObjectIncrementTabLevel((PetscObject)svd->V,(PetscObject)svd,0));
      PetscCall(PetscObjectSetOptions((PetscObject)svd->V,((PetscObject)svd)->options));
    }
    *V = svd->V;
  }
  if (U) {
    if (!svd->U) {
      PetscCall(BVCreate(PetscObjectComm((PetscObject)svd),&svd->U));
      PetscCall(PetscObjectIncrementTabLevel((PetscObject)svd->U,(PetscObject)svd,0));
      PetscCall(PetscObjectSetOptions((PetscObject)svd->U,((PetscObject)svd)->options));
    }
    *U = svd->U;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDSetDS - Associates a direct solver object to the singular value solver.

   Collective

   Input Parameters:
+  svd - singular value solver context obtained from SVDCreate()
-  ds  - the direct solver object

   Note:
   Use SVDGetDS() to retrieve the direct solver context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: SVDGetDS()
@*/
PetscErrorCode SVDSetDS(SVD svd,DS ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(ds,DS_CLASSID,2);
  PetscCheckSameComm(svd,1,ds,2);
  PetscCall(PetscObjectReference((PetscObject)ds));
  PetscCall(DSDestroy(&svd->ds));
  svd->ds = ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDGetDS - Obtain the direct solver object associated to the singular value
   solver object.

   Not Collective

   Input Parameters:
.  svd - singular value solver context obtained from SVDCreate()

   Output Parameter:
.  ds - direct solver context

   Level: advanced

.seealso: SVDSetDS()
@*/
PetscErrorCode SVDGetDS(SVD svd,DS *ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(ds,2);
  if (!svd->ds) {
    PetscCall(DSCreate(PetscObjectComm((PetscObject)svd),&svd->ds));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)svd->ds,(PetscObject)svd,0));
    PetscCall(PetscObjectSetOptions((PetscObject)svd->ds,((PetscObject)svd)->options));
  }
  *ds = svd->ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}
