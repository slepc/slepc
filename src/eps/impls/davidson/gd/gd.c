/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "gd"

   Method: Generalized Davidson

   Algorithm:

       Generalized Davidson with various subspace extraction and
       restart techniques.

   References:

       [1] E.R. Davidson, "Super-matrix methods", Comput. Phys. Commun.
           53(2):49-60, 1989.

       [2] E. Romero and J.E. Roman, "A parallel implementation of
           Davidson methods for large-scale eigenvalue problems in
           SLEPc", ACM Trans. Math. Software 40(2), Article 13, 2014.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <../src/eps/impls/davidson/davidson.h>

static PetscErrorCode EPSSetFromOptions_GD(EPS eps,PetscOptionItems PetscOptionsObject)
{
  PetscBool      flg,flg2,op,orth;
  PetscInt       opi,opi0;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS Generalized Davidson (GD) Options");

    PetscCall(EPSGDGetKrylovStart(eps,&op));
    PetscCall(PetscOptionsBool("-eps_gd_krylov_start","Start the search subspace with a Krylov basis","EPSGDSetKrylovStart",op,&op,&flg));
    if (flg) PetscCall(EPSGDSetKrylovStart(eps,op));

    PetscCall(EPSGDGetBOrth(eps,&orth));
    PetscCall(PetscOptionsBool("-eps_gd_borth","Use B-orthogonalization in the search subspace","EPSGDSetBOrth",op,&op,&flg));
    if (flg) PetscCall(EPSGDSetBOrth(eps,op));

    PetscCall(EPSGDGetBlockSize(eps,&opi));
    PetscCall(PetscOptionsInt("-eps_gd_blocksize","Number of vectors to add to the search subspace","EPSGDSetBlockSize",opi,&opi,&flg));
    if (flg) PetscCall(EPSGDSetBlockSize(eps,opi));

    PetscCall(EPSGDGetRestart(eps,&opi,&opi0));
    PetscCall(PetscOptionsInt("-eps_gd_minv","Size of the search subspace after restarting","EPSGDSetRestart",opi,&opi,&flg));
    PetscCall(PetscOptionsInt("-eps_gd_plusk","Number of eigenvectors saved from the previous iteration when restarting","EPSGDSetRestart",opi0,&opi0,&flg2));
    if (flg || flg2) PetscCall(EPSGDSetRestart(eps,opi,opi0));

    PetscCall(EPSGDGetInitialSize(eps,&opi));
    PetscCall(PetscOptionsInt("-eps_gd_initial_size","Initial size of the search subspace","EPSGDSetInitialSize",opi,&opi,&flg));
    if (flg) PetscCall(EPSGDSetInitialSize(eps,opi));

    PetscCall(PetscOptionsBool("-eps_gd_double_expansion","Use the doble-expansion variant of GD","EPSGDSetDoubleExpansion",PETSC_FALSE,&op,&flg));
    if (flg) PetscCall(EPSGDSetDoubleExpansion(eps,op));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUp_GD(EPS eps)
{
  PetscBool      t;
  KSP            ksp;

  PetscFunctionBegin;
  /* Setup common for all davidson solvers */
  PetscCall(EPSSetUp_XD(eps));

  /* Check some constraints */
  PetscCall(STGetKSP(eps->st,&ksp));
  PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&t));
  PetscCheck(t,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"EPSGD only works with KSPPREONLY");
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_GD(EPS eps,PetscViewer viewer)
{
  PetscBool      isascii,opb;
  PetscInt       opi,opi0;
  PetscBool      borth;
  EPS_DAVIDSON   *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (data->doubleexp) PetscCall(PetscViewerASCIIPrintf(viewer,"  using double expansion variant (GD2)\n"));
    PetscCall(EPSXDGetBOrth_XD(eps,&borth));
    if (borth) PetscCall(PetscViewerASCIIPrintf(viewer,"  search subspace is B-orthogonalized\n"));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  search subspace is orthogonalized\n"));
    PetscCall(EPSXDGetBlockSize_XD(eps,&opi));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  block size=%" PetscInt_FMT "\n",opi));
    PetscCall(EPSXDGetKrylovStart_XD(eps,&opb));
    if (!opb) PetscCall(PetscViewerASCIIPrintf(viewer,"  type of the initial subspace: non-Krylov\n"));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  type of the initial subspace: Krylov\n"));
    PetscCall(EPSXDGetRestart_XD(eps,&opi,&opi0));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  size of the subspace after restarting: %" PetscInt_FMT "\n",opi));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of vectors after restarting from the previous iteration: %" PetscInt_FMT "\n",opi0));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_GD(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetKrylovStart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetKrylovStart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetBOrth_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetBOrth_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetInitialSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetInitialSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetDoubleExpansion_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetDoubleExpansion_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDSetKrylovStart - Activates or deactivates starting the searching
   subspace with a Krylov basis.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  krylovstart - boolean flag

   Options Database Key:
.  -eps_gd_krylov_start - Activates starting the searching subspace with a
    Krylov basis

   Level: advanced

.seealso: EPSGDGetKrylovStart()
@*/
PetscErrorCode EPSGDSetKrylovStart(EPS eps,PetscBool krylovstart)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,krylovstart,2);
  PetscTryMethod(eps,"EPSGDSetKrylovStart_C",(EPS,PetscBool),(eps,krylovstart));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDGetKrylovStart - Returns a flag indicating if the search subspace is started with a
   Krylov basis.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  krylovstart - boolean flag indicating if the search subspace is started
   with a Krylov basis

   Level: advanced

.seealso: EPSGDSetKrylovStart()
@*/
PetscErrorCode EPSGDGetKrylovStart(EPS eps,PetscBool *krylovstart)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(krylovstart,2);
  PetscUseMethod(eps,"EPSGDGetKrylovStart_C",(EPS,PetscBool*),(eps,krylovstart));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDSetBlockSize - Sets the number of vectors to be added to the searching space
   in every iteration.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  blocksize - number of vectors added to the search space in every iteration

   Options Database Key:
.  -eps_gd_blocksize - number of vectors added to the search space in every iteration

   Level: advanced

.seealso: EPSGDSetKrylovStart()
@*/
PetscErrorCode EPSGDSetBlockSize(EPS eps,PetscInt blocksize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,blocksize,2);
  PetscTryMethod(eps,"EPSGDSetBlockSize_C",(EPS,PetscInt),(eps,blocksize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDGetBlockSize - Returns the number of vectors to be added to the searching space
   in every iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  blocksize - number of vectors added to the search space in every iteration

   Level: advanced

.seealso: EPSGDSetBlockSize()
@*/
PetscErrorCode EPSGDGetBlockSize(EPS eps,PetscInt *blocksize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(blocksize,2);
  PetscUseMethod(eps,"EPSGDGetBlockSize_C",(EPS,PetscInt*),(eps,blocksize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDSetRestart - Sets the number of vectors of the searching space after
   restarting and the number of vectors saved from the previous iteration.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  minv - number of vectors of the searching subspace after restarting
-  plusk - number of vectors saved from the previous iteration

   Options Database Keys:
+  -eps_gd_minv - number of vectors of the searching subspace after restarting
-  -eps_gd_plusk - number of vectors saved from the previous iteration

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: EPSGDGetRestart()
@*/
PetscErrorCode EPSGDSetRestart(EPS eps,PetscInt minv,PetscInt plusk)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,minv,2);
  PetscValidLogicalCollectiveInt(eps,plusk,3);
  PetscTryMethod(eps,"EPSGDSetRestart_C",(EPS,PetscInt,PetscInt),(eps,minv,plusk));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDGetRestart - Gets the number of vectors of the searching space after
   restarting and the number of vectors saved from the previous iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  minv - number of vectors of the searching subspace after restarting
-  plusk - number of vectors saved from the previous iteration

   Level: advanced

.seealso: EPSGDSetRestart()
@*/
PetscErrorCode EPSGDGetRestart(EPS eps,PetscInt *minv,PetscInt *plusk)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSGDGetRestart_C",(EPS,PetscInt*,PetscInt*),(eps,minv,plusk));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDSetInitialSize - Sets the initial size of the searching space.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  initialsize - number of vectors of the initial searching subspace

   Options Database Key:
.  -eps_gd_initial_size - number of vectors of the initial searching subspace

   Notes:
   If EPSGDGetKrylovStart() is PETSC_FALSE and the user provides vectors with
   EPSSetInitialSpace(), up to initialsize vectors will be used; and if the
   provided vectors are not enough, the solver completes the subspace with
   random vectors. In the case of EPSGDGetKrylovStart() being PETSC_TRUE, the solver
   gets the first vector provided by the user or, if not available, a random vector,
   and expands the Krylov basis up to initialsize vectors.

   Level: advanced

.seealso: EPSGDGetInitialSize(), EPSGDGetKrylovStart()
@*/
PetscErrorCode EPSGDSetInitialSize(EPS eps,PetscInt initialsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,initialsize,2);
  PetscTryMethod(eps,"EPSGDSetInitialSize_C",(EPS,PetscInt),(eps,initialsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDGetInitialSize - Returns the initial size of the searching space.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  initialsize - number of vectors of the initial searching subspace

   Notes:
   If EPSGDGetKrylovStart() is PETSC_FALSE and the user provides vectors with
   EPSSetInitialSpace(), up to initialsize vectors will be used; and if the
   provided vectors are not enough, the solver completes the subspace with
   random vectors. In the case of EPSGDGetKrylovStart() being PETSC_TRUE, the solver
   gets the first vector provided by the user or, if not available, a random vector,
   and expands the Krylov basis up to initialsize vectors.

   Level: advanced

.seealso: EPSGDSetInitialSize(), EPSGDGetKrylovStart()
@*/
PetscErrorCode EPSGDGetInitialSize(EPS eps,PetscInt *initialsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(initialsize,2);
  PetscUseMethod(eps,"EPSGDGetInitialSize_C",(EPS,PetscInt*),(eps,initialsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDSetBOrth - Selects the orthogonalization that will be used in the search
   subspace in case of generalized Hermitian problems.

   Logically Collective

   Input Parameters:
+  eps   - the eigenproblem solver context
-  borth - whether to B-orthogonalize the search subspace

   Options Database Key:
.  -eps_gd_borth - Set the orthogonalization used in the search subspace

   Level: advanced

.seealso: EPSGDGetBOrth()
@*/
PetscErrorCode EPSGDSetBOrth(EPS eps,PetscBool borth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,borth,2);
  PetscTryMethod(eps,"EPSGDSetBOrth_C",(EPS,PetscBool),(eps,borth));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDGetBOrth - Returns the orthogonalization used in the search
   subspace in case of generalized Hermitian problems.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  borth - whether to B-orthogonalize the search subspace

   Level: advanced

.seealso: EPSGDSetBOrth()
@*/
PetscErrorCode EPSGDGetBOrth(EPS eps,PetscBool *borth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(borth,2);
  PetscUseMethod(eps,"EPSGDGetBOrth_C",(EPS,PetscBool*),(eps,borth));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSGDSetDoubleExpansion_GD(EPS eps,PetscBool doubleexp)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  data->doubleexp = doubleexp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDSetDoubleExpansion - Activate the double expansion variant of GD.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  doubleexp - the boolean flag

   Options Database Keys:
.  -eps_gd_double_expansion - activate the double-expansion variant of GD

   Notes:
   In the double expansion variant the search subspace is expanded with K*[A*x B*x]
   instead of the classic K*r, where K is the preconditioner, x the selected
   approximate eigenvector and r its associated residual vector.

   Level: advanced

.seealso: EPSGDGetDoubleExpansion()
@*/
PetscErrorCode EPSGDSetDoubleExpansion(EPS eps,PetscBool doubleexp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,doubleexp,2);
  PetscTryMethod(eps,"EPSGDSetDoubleExpansion_C",(EPS,PetscBool),(eps,doubleexp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSGDGetDoubleExpansion_GD(EPS eps,PetscBool *doubleexp)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *doubleexp = data->doubleexp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGDGetDoubleExpansion - Gets a flag indicating whether the double
   expansion variant has been activated or not.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  doubleexp - the flag

   Level: advanced

.seealso: EPSGDSetDoubleExpansion()
@*/
PetscErrorCode EPSGDGetDoubleExpansion(EPS eps,PetscBool *doubleexp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(doubleexp,2);
  PetscUseMethod(eps,"EPSGDGetDoubleExpansion_C",(EPS,PetscBool*),(eps,doubleexp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_GD(EPS eps)
{
  EPS_DAVIDSON    *data;

  PetscFunctionBegin;
  PetscCall(PetscNew(&data));
  eps->data = (void*)data;

  data->blocksize   = 1;
  data->initialsize = 0;
  data->minv        = 0;
  data->plusk       = PETSC_DETERMINE;
  data->ipB         = PETSC_TRUE;
  data->fix         = 0.0;
  data->krylovstart = PETSC_FALSE;
  data->dynamic     = PETSC_FALSE;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_XD;
  eps->ops->setup          = EPSSetUp_GD;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->setfromoptions = EPSSetFromOptions_GD;
  eps->ops->destroy        = EPSDestroy_GD;
  eps->ops->reset          = EPSReset_XD;
  eps->ops->view           = EPSView_GD;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_XD;
  eps->ops->setdefaultst   = EPSSetDefaultST_Precond;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetKrylovStart_C",EPSXDSetKrylovStart_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetKrylovStart_C",EPSXDGetKrylovStart_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetBOrth_C",EPSXDSetBOrth_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetBOrth_C",EPSXDGetBOrth_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetBlockSize_C",EPSXDSetBlockSize_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetBlockSize_C",EPSXDGetBlockSize_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetRestart_C",EPSXDSetRestart_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetRestart_C",EPSXDGetRestart_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetInitialSize_C",EPSXDSetInitialSize_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetInitialSize_C",EPSXDGetInitialSize_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetDoubleExpansion_C",EPSGDSetDoubleExpansion_GD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetDoubleExpansion_C",EPSGDGetDoubleExpansion_GD));
  PetscFunctionReturn(PETSC_SUCCESS);
}
