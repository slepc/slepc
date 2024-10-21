/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic PEP routines
*/

#include <slepc/private/pepimpl.h>      /*I "slepcpep.h" I*/

/* Logging support */
PetscClassId      PEP_CLASSID = 0;
PetscLogEvent     PEP_SetUp = 0,PEP_Solve = 0,PEP_Refine = 0,PEP_CISS_SVD = 0;

/* List of registered PEP routines */
PetscFunctionList PEPList = NULL;
PetscBool         PEPRegisterAllCalled = PETSC_FALSE;

/* List of registered PEP monitors */
PetscFunctionList PEPMonitorList              = NULL;
PetscFunctionList PEPMonitorCreateList        = NULL;
PetscFunctionList PEPMonitorDestroyList       = NULL;
PetscBool         PEPMonitorRegisterAllCalled = PETSC_FALSE;

/*@
   PEPCreate - Creates the default PEP context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  outpep - location to put the PEP context

   Note:
   The default PEP type is PEPTOAR

   Level: beginner

.seealso: PEPSetUp(), PEPSolve(), PEPDestroy(), PEP
@*/
PetscErrorCode PEPCreate(MPI_Comm comm,PEP *outpep)
{
  PEP            pep;

  PetscFunctionBegin;
  PetscAssertPointer(outpep,2);
  PetscCall(PEPInitializePackage());
  PetscCall(SlepcHeaderCreate(pep,PEP_CLASSID,"PEP","Polynomial Eigenvalue Problem","PEP",comm,PEPDestroy,PEPView));

  pep->max_it          = PETSC_DETERMINE;
  pep->nev             = 1;
  pep->ncv             = PETSC_DETERMINE;
  pep->mpd             = PETSC_DETERMINE;
  pep->nini            = 0;
  pep->target          = 0.0;
  pep->tol             = PETSC_DETERMINE;
  pep->conv            = PEP_CONV_REL;
  pep->stop            = PEP_STOP_BASIC;
  pep->which           = (PEPWhich)0;
  pep->basis           = PEP_BASIS_MONOMIAL;
  pep->problem_type    = (PEPProblemType)0;
  pep->scale           = PEP_SCALE_NONE;
  pep->sfactor         = 1.0;
  pep->dsfactor        = 1.0;
  pep->sits            = 5;
  pep->slambda         = 1.0;
  pep->refine          = PEP_REFINE_NONE;
  pep->npart           = 1;
  pep->rtol            = PETSC_DETERMINE;
  pep->rits            = PETSC_DETERMINE;
  pep->scheme          = (PEPRefineScheme)0;
  pep->extract         = (PEPExtract)0;
  pep->trackall        = PETSC_FALSE;

  pep->converged       = PEPConvergedRelative;
  pep->convergeduser   = NULL;
  pep->convergeddestroy= NULL;
  pep->stopping        = PEPStoppingBasic;
  pep->stoppinguser    = NULL;
  pep->stoppingdestroy = NULL;
  pep->convergedctx    = NULL;
  pep->stoppingctx     = NULL;
  pep->numbermonitors  = 0;

  pep->st              = NULL;
  pep->ds              = NULL;
  pep->V               = NULL;
  pep->rg              = NULL;
  pep->A               = NULL;
  pep->nmat            = 0;
  pep->Dl              = NULL;
  pep->Dr              = NULL;
  pep->IS              = NULL;
  pep->eigr            = NULL;
  pep->eigi            = NULL;
  pep->errest          = NULL;
  pep->perm            = NULL;
  pep->pbc             = NULL;
  pep->solvematcoeffs  = NULL;
  pep->nwork           = 0;
  pep->work            = NULL;
  pep->refineksp       = NULL;
  pep->refinesubc      = NULL;
  pep->data            = NULL;

  pep->state           = PEP_STATE_INITIAL;
  pep->nconv           = 0;
  pep->its             = 0;
  pep->n               = 0;
  pep->nloc            = 0;
  pep->nrma            = NULL;
  pep->sfactor_set     = PETSC_FALSE;
  pep->lineariz        = PETSC_FALSE;
  pep->reason          = PEP_CONVERGED_ITERATING;

  PetscCall(PetscNew(&pep->sc));
  *outpep = pep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetType - Selects the particular solver to be used in the PEP object.

   Logically Collective

   Input Parameters:
+  pep      - the polynomial eigensolver context
-  type     - a known method

   Options Database Key:
.  -pep_type <method> - Sets the method; use -help for a list
    of available methods

   Notes:
   See "slepc/include/slepcpep.h" for available methods. The default
   is PEPTOAR.

   Normally, it is best to use the PEPSetFromOptions() command and
   then set the PEP type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The PEPSetType() routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: PEPType
@*/
PetscErrorCode PEPSetType(PEP pep,PEPType type)
{
  PetscErrorCode (*r)(PEP);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)pep,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscFunctionListFind(PEPList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown PEP type given: %s",type);

  PetscTryTypeMethod(pep,destroy);
  PetscCall(PetscMemzero(pep->ops,sizeof(struct _PEPOps)));

  pep->state = PEP_STATE_INITIAL;
  PetscCall(PetscObjectChangeTypeName((PetscObject)pep,type));
  PetscCall((*r)(pep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetType - Gets the PEP type as a string from the PEP object.

   Not Collective

   Input Parameter:
.  pep - the eigensolver context

   Output Parameter:
.  type - name of PEP method

   Level: intermediate

.seealso: PEPSetType()
@*/
PetscErrorCode PEPGetType(PEP pep,PEPType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)pep)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPRegister - Adds a method to the polynomial eigenproblem solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Notes:
   PEPRegister() may be called multiple times to add several user-defined solvers.

   Example Usage:
.vb
    PEPRegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     PEPSetType(pep,"my_solver")
   or at runtime via the option
$     -pep_type my_solver

   Level: advanced

.seealso: PEPRegisterAll()
@*/
PetscErrorCode PEPRegister(const char *name,PetscErrorCode (*function)(PEP))
{
  PetscFunctionBegin;
  PetscCall(PEPInitializePackage());
  PetscCall(PetscFunctionListAdd(&PEPList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorRegister - Adds PEP monitor routine.

   Not Collective

   Input Parameters:
+  name    - name of a new monitor routine
.  vtype   - a PetscViewerType for the output
.  format  - a PetscViewerFormat for the output
.  monitor - monitor routine
.  create  - creation routine, or NULL
-  destroy - destruction routine, or NULL

   Notes:
   PEPMonitorRegister() may be called multiple times to add several user-defined monitors.

   Example Usage:
.vb
   PEPMonitorRegister("my_monitor",PETSCVIEWERASCII,PETSC_VIEWER_ASCII_INFO_DETAIL,MyMonitor,NULL,NULL);
.ve

   Then, your monitor can be chosen with the procedural interface via
$      PEPMonitorSetFromOptions(pep,"-pep_monitor_my_monitor","my_monitor",NULL)
   or at runtime via the option
$      -pep_monitor_my_monitor

   Level: advanced

.seealso: PEPMonitorRegisterAll()
@*/
PetscErrorCode PEPMonitorRegister(const char name[],PetscViewerType vtype,PetscViewerFormat format,PetscErrorCode (*monitor)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscErrorCode (*create)(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**),PetscErrorCode (*destroy)(PetscViewerAndFormat**))
{
  char           key[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  PetscCall(PEPInitializePackage());
  PetscCall(SlepcMonitorMakeKey_Internal(name,vtype,format,key));
  PetscCall(PetscFunctionListAdd(&PEPMonitorList,key,monitor));
  if (create)  PetscCall(PetscFunctionListAdd(&PEPMonitorCreateList,key,create));
  if (destroy) PetscCall(PetscFunctionListAdd(&PEPMonitorDestroyList,key,destroy));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPReset - Resets the PEP context to the initial state (prior to setup)
   and destroys any allocated Vecs and Mats.

   Collective

   Input Parameter:
.  pep - eigensolver context obtained from PEPCreate()

   Level: advanced

.seealso: PEPDestroy()
@*/
PetscErrorCode PEPReset(PEP pep)
{
  PetscFunctionBegin;
  if (pep) PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!pep) PetscFunctionReturn(PETSC_SUCCESS);
  PetscTryTypeMethod(pep,reset);
  if (pep->st) PetscCall(STReset(pep->st));
  if (pep->refineksp) PetscCall(KSPReset(pep->refineksp));
  if (pep->nmat) {
    PetscCall(MatDestroyMatrices(pep->nmat,&pep->A));
    PetscCall(PetscFree2(pep->pbc,pep->nrma));
    PetscCall(PetscFree(pep->solvematcoeffs));
    pep->nmat = 0;
  }
  PetscCall(VecDestroy(&pep->Dl));
  PetscCall(VecDestroy(&pep->Dr));
  PetscCall(BVDestroy(&pep->V));
  PetscCall(VecDestroyVecs(pep->nwork,&pep->work));
  pep->nwork = 0;
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPDestroy - Destroys the PEP context.

   Collective

   Input Parameter:
.  pep - eigensolver context obtained from PEPCreate()

   Level: beginner

.seealso: PEPCreate(), PEPSetUp(), PEPSolve()
@*/
PetscErrorCode PEPDestroy(PEP *pep)
{
  PetscFunctionBegin;
  if (!*pep) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*pep,PEP_CLASSID,1);
  if (--((PetscObject)*pep)->refct > 0) { *pep = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscCall(PEPReset(*pep));
  PetscTryTypeMethod(*pep,destroy);
  if ((*pep)->eigr) PetscCall(PetscFree4((*pep)->eigr,(*pep)->eigi,(*pep)->errest,(*pep)->perm));
  PetscCall(STDestroy(&(*pep)->st));
  PetscCall(RGDestroy(&(*pep)->rg));
  PetscCall(DSDestroy(&(*pep)->ds));
  PetscCall(KSPDestroy(&(*pep)->refineksp));
  PetscCall(PetscSubcommDestroy(&(*pep)->refinesubc));
  PetscCall(PetscFree((*pep)->sc));
  /* just in case the initial vectors have not been used */
  PetscCall(SlepcBasisDestroy_Private(&(*pep)->nini,&(*pep)->IS));
  if ((*pep)->convergeddestroy) PetscCall((*(*pep)->convergeddestroy)(&(*pep)->convergedctx));
  if ((*pep)->stoppingdestroy) PetscCall((*(*pep)->stoppingdestroy)(&(*pep)->stoppingctx));
  PetscCall(PEPMonitorCancel(*pep));
  PetscCall(PetscHeaderDestroy(pep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetBV - Associates a basis vectors object to the polynomial eigensolver.

   Collective

   Input Parameters:
+  pep - eigensolver context obtained from PEPCreate()
-  bv  - the basis vectors object

   Note:
   Use PEPGetBV() to retrieve the basis vectors context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: PEPGetBV()
@*/
PetscErrorCode PEPSetBV(PEP pep,BV bv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(bv,BV_CLASSID,2);
  PetscCheckSameComm(pep,1,bv,2);
  PetscCall(PetscObjectReference((PetscObject)bv));
  PetscCall(BVDestroy(&pep->V));
  pep->V = bv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetBV - Obtain the basis vectors object associated to the polynomial
   eigensolver object.

   Not Collective

   Input Parameters:
.  pep - eigensolver context obtained from PEPCreate()

   Output Parameter:
.  bv - basis vectors context

   Level: advanced

.seealso: PEPSetBV()
@*/
PetscErrorCode PEPGetBV(PEP pep,BV *bv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(bv,2);
  if (!pep->V) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)pep),&pep->V));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)pep->V,(PetscObject)pep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)pep->V,((PetscObject)pep)->options));
  }
  *bv = pep->V;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetRG - Associates a region object to the polynomial eigensolver.

   Collective

   Input Parameters:
+  pep - eigensolver context obtained from PEPCreate()
-  rg  - the region object

   Note:
   Use PEPGetRG() to retrieve the region context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: PEPGetRG()
@*/
PetscErrorCode PEPSetRG(PEP pep,RG rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (rg) {
    PetscValidHeaderSpecific(rg,RG_CLASSID,2);
    PetscCheckSameComm(pep,1,rg,2);
  }
  PetscCall(PetscObjectReference((PetscObject)rg));
  PetscCall(RGDestroy(&pep->rg));
  pep->rg = rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetRG - Obtain the region object associated to the
   polynomial eigensolver object.

   Not Collective

   Input Parameters:
.  pep - eigensolver context obtained from PEPCreate()

   Output Parameter:
.  rg - region context

   Level: advanced

.seealso: PEPSetRG()
@*/
PetscErrorCode PEPGetRG(PEP pep,RG *rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(rg,2);
  if (!pep->rg) {
    PetscCall(RGCreate(PetscObjectComm((PetscObject)pep),&pep->rg));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)pep->rg,(PetscObject)pep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)pep->rg,((PetscObject)pep)->options));
  }
  *rg = pep->rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetDS - Associates a direct solver object to the polynomial eigensolver.

   Collective

   Input Parameters:
+  pep - eigensolver context obtained from PEPCreate()
-  ds  - the direct solver object

   Note:
   Use PEPGetDS() to retrieve the direct solver context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: PEPGetDS()
@*/
PetscErrorCode PEPSetDS(PEP pep,DS ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(ds,DS_CLASSID,2);
  PetscCheckSameComm(pep,1,ds,2);
  PetscCall(PetscObjectReference((PetscObject)ds));
  PetscCall(DSDestroy(&pep->ds));
  pep->ds = ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetDS - Obtain the direct solver object associated to the
   polynomial eigensolver object.

   Not Collective

   Input Parameters:
.  pep - eigensolver context obtained from PEPCreate()

   Output Parameter:
.  ds - direct solver context

   Level: advanced

.seealso: PEPSetDS()
@*/
PetscErrorCode PEPGetDS(PEP pep,DS *ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(ds,2);
  if (!pep->ds) {
    PetscCall(DSCreate(PetscObjectComm((PetscObject)pep),&pep->ds));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)pep->ds,(PetscObject)pep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)pep->ds,((PetscObject)pep)->options));
  }
  *ds = pep->ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetST - Associates a spectral transformation object to the eigensolver.

   Collective

   Input Parameters:
+  pep - eigensolver context obtained from PEPCreate()
-  st   - the spectral transformation object

   Note:
   Use PEPGetST() to retrieve the spectral transformation context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: PEPGetST()
@*/
PetscErrorCode PEPSetST(PEP pep,ST st)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(st,ST_CLASSID,2);
  PetscCheckSameComm(pep,1,st,2);
  PetscCall(PetscObjectReference((PetscObject)st));
  PetscCall(STDestroy(&pep->st));
  pep->st = st;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetST - Obtain the spectral transformation (ST) object associated
   to the eigensolver object.

   Not Collective

   Input Parameters:
.  pep - eigensolver context obtained from PEPCreate()

   Output Parameter:
.  st - spectral transformation context

   Level: intermediate

.seealso: PEPSetST()
@*/
PetscErrorCode PEPGetST(PEP pep,ST *st)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(st,2);
  if (!pep->st) {
    PetscCall(STCreate(PetscObjectComm((PetscObject)pep),&pep->st));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)pep->st,(PetscObject)pep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)pep->st,((PetscObject)pep)->options));
  }
  *st = pep->st;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPRefineGetKSP - Obtain the ksp object used by the eigensolver
   object in the refinement phase.

   Collective

   Input Parameters:
.  pep - eigensolver context obtained from PEPCreate()

   Output Parameter:
.  ksp - ksp context

   Level: advanced

.seealso: PEPSetRefine()
@*/
PetscErrorCode PEPRefineGetKSP(PEP pep,KSP *ksp)
{
  MPI_Comm       comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(ksp,2);
  if (!pep->refineksp) {
    if (pep->npart>1) {
      /* Split in subcomunicators */
      PetscCall(PetscSubcommCreate(PetscObjectComm((PetscObject)pep),&pep->refinesubc));
      PetscCall(PetscSubcommSetNumber(pep->refinesubc,pep->npart));
      PetscCall(PetscSubcommSetType(pep->refinesubc,PETSC_SUBCOMM_CONTIGUOUS));
      PetscCall(PetscSubcommGetChild(pep->refinesubc,&comm));
    } else PetscCall(PetscObjectGetComm((PetscObject)pep,&comm));
    PetscCall(KSPCreate(comm,&pep->refineksp));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)pep->refineksp,(PetscObject)pep,0));
    PetscCall(PetscObjectSetOptions((PetscObject)pep->refineksp,((PetscObject)pep)->options));
    PetscCall(KSPSetOptionsPrefix(*ksp,((PetscObject)pep)->prefix));
    PetscCall(KSPAppendOptionsPrefix(*ksp,"pep_refine_"));
    PetscCall(KSPSetTolerances(pep->refineksp,SlepcDefaultTol(pep->rtol),PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
  }
  *ksp = pep->refineksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetTarget - Sets the value of the target.

   Logically Collective

   Input Parameters:
+  pep    - eigensolver context
-  target - the value of the target

   Options Database Key:
.  -pep_target <scalar> - the value of the target

   Notes:
   The target is a scalar value used to determine the portion of the spectrum
   of interest. It is used in combination with PEPSetWhichEigenpairs().

   In the case of complex scalars, a complex value can be provided in the
   command line with [+/-][realnumber][+/-]realnumberi with no spaces, e.g.
   -pep_target 1.0+2.0i

   Level: intermediate

.seealso: PEPGetTarget(), PEPSetWhichEigenpairs()
@*/
PetscErrorCode PEPSetTarget(PEP pep,PetscScalar target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(pep,target,2);
  pep->target = target;
  if (!pep->st) PetscCall(PEPGetST(pep,&pep->st));
  PetscCall(STSetDefaultShift(pep->st,target));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetTarget - Gets the value of the target.

   Not Collective

   Input Parameter:
.  pep - eigensolver context

   Output Parameter:
.  target - the value of the target

   Note:
   If the target was not set by the user, then zero is returned.

   Level: intermediate

.seealso: PEPSetTarget()
@*/
PetscErrorCode PEPGetTarget(PEP pep,PetscScalar* target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(target,2);
  *target = pep->target;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetInterval - Defines the computational interval for spectrum slicing.

   Logically Collective

   Input Parameters:
+  pep  - eigensolver context
.  inta - left end of the interval
-  intb - right end of the interval

   Options Database Key:
.  -pep_interval <a,b> - set [a,b] as the interval of interest

   Notes:
   Spectrum slicing is a technique employed for computing all eigenvalues of
   symmetric eigenproblems in a given interval. This function provides the
   interval to be considered. It must be used in combination with PEP_ALL, see
   PEPSetWhichEigenpairs().

   In the command-line option, two values must be provided. For an open interval,
   one can give an infinite, e.g., -pep_interval 1.0,inf or -pep_interval -inf,1.0.
   An open interval in the programmatic interface can be specified with
   PETSC_MAX_REAL and -PETSC_MAX_REAL.

   Level: intermediate

.seealso: PEPGetInterval(), PEPSetWhichEigenpairs()
@*/
PetscErrorCode PEPSetInterval(PEP pep,PetscReal inta,PetscReal intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,inta,2);
  PetscValidLogicalCollectiveReal(pep,intb,3);
  PetscCheck(inta<intb,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be inta<intb");
  if (pep->inta != inta || pep->intb != intb) {
    pep->inta = inta;
    pep->intb = intb;
    pep->state = PEP_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetInterval - Gets the computational interval for spectrum slicing.

   Not Collective

   Input Parameter:
.  pep - eigensolver context

   Output Parameters:
+  inta - left end of the interval
-  intb - right end of the interval

   Level: intermediate

   Note:
   If the interval was not set by the user, then zeros are returned.

.seealso: PEPSetInterval()
@*/
PetscErrorCode PEPGetInterval(PEP pep,PetscReal* inta,PetscReal* intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (inta) *inta = pep->inta;
  if (intb) *intb = pep->intb;
  PetscFunctionReturn(PETSC_SUCCESS);
}
