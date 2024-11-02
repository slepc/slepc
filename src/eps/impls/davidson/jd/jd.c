/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "jd"

   Method: Jacobi-Davidson

   Algorithm:

       Jacobi-Davidson with various subspace extraction and
       restart techniques.

   References:

       [1] G.L.G. Sleijpen and H.A. van der Vorst, "A Jacobi-Davidson
           iteration method for linear eigenvalue problems", SIAM J.
           Matrix Anal. Appl. 17(2):401-425, 1996.

       [2] E. Romero and J.E. Roman, "A parallel implementation of
           Davidson methods for large-scale eigenvalue problems in
           SLEPc", ACM Trans. Math. Software 40(2), Article 13, 2014.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <../src/eps/impls/davidson/davidson.h>

static PetscErrorCode EPSSetFromOptions_JD(EPS eps,PetscOptionItems PetscOptionsObject)
{
  PetscBool      flg,flg2,op,orth;
  PetscInt       opi,opi0;
  PetscReal      opf;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS Jacobi-Davidson (JD) Options");

    PetscCall(EPSJDGetKrylovStart(eps,&op));
    PetscCall(PetscOptionsBool("-eps_jd_krylov_start","Start the search subspace with a Krylov basis","EPSJDSetKrylovStart",op,&op,&flg));
    if (flg) PetscCall(EPSJDSetKrylovStart(eps,op));

    PetscCall(EPSJDGetBOrth(eps,&orth));
    PetscCall(PetscOptionsBool("-eps_jd_borth","Use B-orthogonalization in the search subspace","EPSJDSetBOrth",op,&op,&flg));
    if (flg) PetscCall(EPSJDSetBOrth(eps,op));

    PetscCall(EPSJDGetBlockSize(eps,&opi));
    PetscCall(PetscOptionsInt("-eps_jd_blocksize","Number of vectors to add to the search subspace","EPSJDSetBlockSize",opi,&opi,&flg));
    if (flg) PetscCall(EPSJDSetBlockSize(eps,opi));

    PetscCall(EPSJDGetRestart(eps,&opi,&opi0));
    PetscCall(PetscOptionsInt("-eps_jd_minv","Size of the search subspace after restarting","EPSJDSetRestart",opi,&opi,&flg));
    PetscCall(PetscOptionsInt("-eps_jd_plusk","Number of eigenvectors saved from the previous iteration when restarting","EPSJDSetRestart",opi0,&opi0,&flg2));
    if (flg || flg2) PetscCall(EPSJDSetRestart(eps,opi,opi0));

    PetscCall(EPSJDGetInitialSize(eps,&opi));
    PetscCall(PetscOptionsInt("-eps_jd_initial_size","Initial size of the search subspace","EPSJDSetInitialSize",opi,&opi,&flg));
    if (flg) PetscCall(EPSJDSetInitialSize(eps,opi));

    PetscCall(EPSJDGetFix(eps,&opf));
    PetscCall(PetscOptionsReal("-eps_jd_fix","Tolerance for changing the target in the correction equation","EPSJDSetFix",opf,&opf,&flg));
    if (flg) PetscCall(EPSJDSetFix(eps,opf));

    PetscCall(EPSJDGetConstCorrectionTol(eps,&op));
    PetscCall(PetscOptionsBool("-eps_jd_const_correction_tol","Disable the dynamic stopping criterion when solving the correction equation","EPSJDSetConstCorrectionTol",op,&op,&flg));
    if (flg) PetscCall(EPSJDSetConstCorrectionTol(eps,op));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetDefaultST_JD(EPS eps)
{
  KSP            ksp;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) {
    PetscCall(STSetType(eps->st,STPRECOND));
    PetscCall(STPrecondSetKSPHasMat(eps->st,PETSC_TRUE));
  }
  PetscCall(STGetKSP(eps->st,&ksp));
  if (!((PetscObject)ksp)->type_name) {
    PetscCall(KSPSetType(ksp,KSPBCGSL));
    PetscCall(KSPSetTolerances(ksp,1e-4,PETSC_CURRENT,PETSC_CURRENT,90));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUp_JD(EPS eps)
{
  PetscBool      t;
  KSP            ksp;

  PetscFunctionBegin;
  /* Setup common for all davidson solvers */
  PetscCall(EPSSetUp_XD(eps));

  /* Check some constraints */
  PetscCall(STGetKSP(eps->st,&ksp));
  PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&t));
  PetscCheck(!t,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"EPSJD does not work with KSPPREONLY");
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_JD(EPS eps,PetscViewer viewer)
{
  PetscBool      isascii,opb;
  PetscReal      opf;
  PetscInt       opi,opi0;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(EPSXDGetBOrth_XD(eps,&opb));
    if (opb) PetscCall(PetscViewerASCIIPrintf(viewer,"  search subspace is B-orthogonalized\n"));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  search subspace is orthogonalized\n"));
    PetscCall(EPSXDGetBlockSize_XD(eps,&opi));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  block size=%" PetscInt_FMT "\n",opi));
    PetscCall(EPSXDGetKrylovStart_XD(eps,&opb));
    if (!opb) PetscCall(PetscViewerASCIIPrintf(viewer,"  type of the initial subspace: non-Krylov\n"));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  type of the initial subspace: Krylov\n"));
    PetscCall(EPSXDGetRestart_XD(eps,&opi,&opi0));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  size of the subspace after restarting: %" PetscInt_FMT "\n",opi));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of vectors after restarting from the previous iteration: %" PetscInt_FMT "\n",opi0));

    PetscCall(EPSJDGetFix_JD(eps,&opf));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  threshold for changing the target in the correction equation (fix): %g\n",(double)opf));

    PetscCall(EPSJDGetConstCorrectionTol_JD(eps,&opb));
    if (!opb) PetscCall(PetscViewerASCIIPrintf(viewer,"  using dynamic tolerance for the correction equation\n"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_JD(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetKrylovStart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetKrylovStart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetInitialSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetInitialSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetFix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetFix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetConstCorrectionTol_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetConstCorrectionTol_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetBOrth_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetBOrth_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDSetKrylovStart - Activates or deactivates starting the searching
   subspace with a Krylov basis.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  krylovstart - boolean flag

   Options Database Key:
.  -eps_jd_krylov_start - Activates starting the searching subspace with a
    Krylov basis

   Level: advanced

.seealso: EPSJDGetKrylovStart()
@*/
PetscErrorCode EPSJDSetKrylovStart(EPS eps,PetscBool krylovstart)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,krylovstart,2);
  PetscTryMethod(eps,"EPSJDSetKrylovStart_C",(EPS,PetscBool),(eps,krylovstart));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDGetKrylovStart - Returns a flag indicating if the searching subspace is started with a
   Krylov basis.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  krylovstart - boolean flag indicating if the searching subspace is started
   with a Krylov basis

   Level: advanced

.seealso: EPSJDSetKrylovStart()
@*/
PetscErrorCode EPSJDGetKrylovStart(EPS eps,PetscBool *krylovstart)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(krylovstart,2);
  PetscUseMethod(eps,"EPSJDGetKrylovStart_C",(EPS,PetscBool*),(eps,krylovstart));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDSetBlockSize - Sets the number of vectors to be added to the searching space
   in every iteration.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  blocksize - number of vectors added to the search space in every iteration

   Options Database Key:
.  -eps_jd_blocksize - number of vectors added to the searching space every iteration

   Level: advanced

.seealso: EPSJDSetKrylovStart()
@*/
PetscErrorCode EPSJDSetBlockSize(EPS eps,PetscInt blocksize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,blocksize,2);
  PetscTryMethod(eps,"EPSJDSetBlockSize_C",(EPS,PetscInt),(eps,blocksize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDGetBlockSize - Returns the number of vectors to be added to the searching space
   in every iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  blocksize - number of vectors added to the search space in every iteration

   Level: advanced

.seealso: EPSJDSetBlockSize()
@*/
PetscErrorCode EPSJDGetBlockSize(EPS eps,PetscInt *blocksize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(blocksize,2);
  PetscUseMethod(eps,"EPSJDGetBlockSize_C",(EPS,PetscInt*),(eps,blocksize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDSetRestart - Sets the number of vectors of the searching space after
   restarting and the number of vectors saved from the previous iteration.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  minv - number of vectors of the searching subspace after restarting
-  plusk - number of vectors saved from the previous iteration

   Options Database Keys:
+  -eps_jd_minv - number of vectors of the searching subspace after restarting
-  -eps_jd_plusk - number of vectors saved from the previous iteration

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: advanced

.seealso: EPSJDGetRestart()
@*/
PetscErrorCode EPSJDSetRestart(EPS eps,PetscInt minv,PetscInt plusk)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,minv,2);
  PetscValidLogicalCollectiveInt(eps,plusk,3);
  PetscTryMethod(eps,"EPSJDSetRestart_C",(EPS,PetscInt,PetscInt),(eps,minv,plusk));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDGetRestart - Gets the number of vectors of the searching space after
   restarting and the number of vectors saved from the previous iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  minv - number of vectors of the searching subspace after restarting
-  plusk - number of vectors saved from the previous iteration

   Level: advanced

.seealso: EPSJDSetRestart()
@*/
PetscErrorCode EPSJDGetRestart(EPS eps,PetscInt *minv,PetscInt *plusk)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSJDGetRestart_C",(EPS,PetscInt*,PetscInt*),(eps,minv,plusk));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDSetInitialSize - Sets the initial size of the searching space.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  initialsize - number of vectors of the initial searching subspace

   Options Database Key:
.  -eps_jd_initial_size - number of vectors of the initial searching subspace

   Notes:
   If EPSJDGetKrylovStart() is PETSC_FALSE and the user provides vectors with
   EPSSetInitialSpace(), up to initialsize vectors will be used; and if the
   provided vectors are not enough, the solver completes the subspace with
   random vectors. In the case of EPSJDGetKrylovStart() being PETSC_TRUE, the solver
   gets the first vector provided by the user or, if not available, a random vector,
   and expands the Krylov basis up to initialsize vectors.

   Level: advanced

.seealso: EPSJDGetInitialSize(), EPSJDGetKrylovStart()
@*/
PetscErrorCode EPSJDSetInitialSize(EPS eps,PetscInt initialsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,initialsize,2);
  PetscTryMethod(eps,"EPSJDSetInitialSize_C",(EPS,PetscInt),(eps,initialsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDGetInitialSize - Returns the initial size of the searching space.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  initialsize - number of vectors of the initial searching subspace

   Notes:
   If EPSJDGetKrylovStart() is PETSC_FALSE and the user provides vectors with
   EPSSetInitialSpace(), up to initialsize vectors will be used; and if the
   provided vectors are not enough, the solver completes the subspace with
   random vectors. In the case of EPSJDGetKrylovStart() being PETSC_TRUE, the solver
   gets the first vector provided by the user or, if not available, a random vector,
   and expands the Krylov basis up to initialsize vectors.

   Level: advanced

.seealso: EPSJDSetInitialSize(), EPSJDGetKrylovStart()
@*/
PetscErrorCode EPSJDGetInitialSize(EPS eps,PetscInt *initialsize)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(initialsize,2);
  PetscUseMethod(eps,"EPSJDGetInitialSize_C",(EPS,PetscInt*),(eps,initialsize));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSJDSetFix_JD(EPS eps,PetscReal fix)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (fix == (PetscReal)PETSC_DEFAULT || fix == (PetscReal)PETSC_DECIDE) fix = 0.01;
  PetscCheck(fix>0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid fix value, must be >0");
  data->fix = fix;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDSetFix - Sets the threshold for changing the target in the correction
   equation.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  fix - threshold for changing the target

   Options Database Key:
.  -eps_jd_fix - the fix value

   Note:
   The target in the correction equation is fixed at the first iterations.
   When the norm of the residual vector is lower than the fix value,
   the target is set to the corresponding eigenvalue.

   Level: advanced

.seealso: EPSJDGetFix()
@*/
PetscErrorCode EPSJDSetFix(EPS eps,PetscReal fix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,fix,2);
  PetscTryMethod(eps,"EPSJDSetFix_C",(EPS,PetscReal),(eps,fix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSJDGetFix_JD(EPS eps,PetscReal *fix)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *fix = data->fix;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDGetFix - Returns the threshold for changing the target in the correction
   equation.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  fix - threshold for changing the target

   Note:
   The target in the correction equation is fixed at the first iterations.
   When the norm of the residual vector is lower than the fix value,
   the target is set to the corresponding eigenvalue.

   Level: advanced

.seealso: EPSJDSetFix()
@*/
PetscErrorCode EPSJDGetFix(EPS eps,PetscReal *fix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(fix,2);
  PetscUseMethod(eps,"EPSJDGetFix_C",(EPS,PetscReal*),(eps,fix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSJDSetConstCorrectionTol_JD(EPS eps,PetscBool constant)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  data->dynamic = PetscNot(constant);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDSetConstCorrectionTol - If true, deactivates the dynamic stopping criterion
   (also called Newton) that sets the KSP relative tolerance
   to 0.5**i, where i is the number of EPS iterations from the last converged value.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  constant - if false, the KSP relative tolerance is set to 0.5**i.

   Options Database Key:
.  -eps_jd_const_correction_tol - Deactivates the dynamic stopping criterion.

   Level: advanced

.seealso: EPSJDGetConstCorrectionTol()
@*/
PetscErrorCode EPSJDSetConstCorrectionTol(EPS eps,PetscBool constant)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,constant,2);
  PetscTryMethod(eps,"EPSJDSetConstCorrectionTol_C",(EPS,PetscBool),(eps,constant));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSJDGetConstCorrectionTol_JD(EPS eps,PetscBool *constant)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *constant = PetscNot(data->dynamic);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDGetConstCorrectionTol - Returns a flag indicating if the dynamic stopping is being used for
   solving the correction equation.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  constant - boolean flag indicating if the dynamic stopping criterion is not being used.

   Notes:
   If the flag is false the KSP relative tolerance is set to 0.5**i, where i is the number
   of EPS iterations from the last converged value.

   Level: advanced

.seealso: EPSJDSetConstCorrectionTol()
@*/
PetscErrorCode EPSJDGetConstCorrectionTol(EPS eps,PetscBool *constant)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(constant,2);
  PetscUseMethod(eps,"EPSJDGetConstCorrectionTol_C",(EPS,PetscBool*),(eps,constant));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDSetBOrth - Selects the orthogonalization that will be used in the search
   subspace in case of generalized Hermitian problems.

   Logically Collective

   Input Parameters:
+  eps   - the eigenproblem solver context
-  borth - whether to B-orthogonalize the search subspace

   Options Database Key:
.  -eps_jd_borth - Set the orthogonalization used in the search subspace

   Level: advanced

.seealso: EPSJDGetBOrth()
@*/
PetscErrorCode EPSJDSetBOrth(EPS eps,PetscBool borth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,borth,2);
  PetscTryMethod(eps,"EPSJDSetBOrth_C",(EPS,PetscBool),(eps,borth));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSJDGetBOrth - Returns the orthogonalization used in the search
   subspace in case of generalized Hermitian problems.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  borth - whether to B-orthogonalize the search subspace

   Level: advanced

.seealso: EPSJDSetBOrth()
@*/
PetscErrorCode EPSJDGetBOrth(EPS eps,PetscBool *borth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(borth,2);
  PetscUseMethod(eps,"EPSJDGetBOrth_C",(EPS,PetscBool*),(eps,borth));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_JD(EPS eps)
{
  EPS_DAVIDSON   *data;

  PetscFunctionBegin;
  PetscCall(PetscNew(&data));
  eps->data = (void*)data;

  data->blocksize   = 1;
  data->initialsize = 0;
  data->minv        = 0;
  data->plusk       = PETSC_DETERMINE;
  data->ipB         = PETSC_TRUE;
  data->fix         = 0.01;
  data->krylovstart = PETSC_FALSE;
  data->dynamic     = PETSC_FALSE;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_XD;
  eps->ops->setup          = EPSSetUp_JD;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->setfromoptions = EPSSetFromOptions_JD;
  eps->ops->destroy        = EPSDestroy_JD;
  eps->ops->reset          = EPSReset_XD;
  eps->ops->view           = EPSView_JD;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_XD;
  eps->ops->setdefaultst   = EPSSetDefaultST_JD;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetKrylovStart_C",EPSXDSetKrylovStart_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetKrylovStart_C",EPSXDGetKrylovStart_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetBlockSize_C",EPSXDSetBlockSize_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetBlockSize_C",EPSXDGetBlockSize_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetRestart_C",EPSXDSetRestart_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetRestart_C",EPSXDGetRestart_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetInitialSize_C",EPSXDSetInitialSize_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetInitialSize_C",EPSXDGetInitialSize_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetFix_C",EPSJDSetFix_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetFix_C",EPSJDGetFix_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetConstCorrectionTol_C",EPSJDSetConstCorrectionTol_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetConstCorrectionTol_C",EPSJDGetConstCorrectionTol_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetBOrth_C",EPSXDSetBOrth_XD));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetBOrth_C",EPSXDGetBOrth_XD));
  PetscFunctionReturn(PETSC_SUCCESS);
}
