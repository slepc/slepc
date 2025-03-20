/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to problem setup
*/

#include <slepc/private/nepimpl.h>       /*I "slepcnep.h" I*/

/*@
   NEPSetDSType - Sets the type of the internal DS object based on the current
   settings of the nonlinear eigensolver.

   Collective

   Input Parameter:
.  nep - nonlinear eigensolver context

   Note:
   This function need not be called explicitly, since it will be called at
   both NEPSetFromOptions() and NEPSetUp().

   Level: developer

.seealso: NEPSetFromOptions(), NEPSetUp()
@*/
PetscErrorCode NEPSetDSType(NEP nep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscTryTypeMethod(nep,setdstype);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetUp - Sets up all the internal data structures necessary for the
   execution of the NEP solver.

   Collective

   Input Parameter:
.  nep   - solver context

   Notes:
   This function need not be called explicitly in most cases, since NEPSolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: NEPCreate(), NEPSolve(), NEPDestroy()
@*/
PetscErrorCode NEPSetUp(NEP nep)
{
  PetscInt       k;
  SlepcSC        sc;
  Mat            T;
  PetscBool      flg;
  KSP            ksp;
  PC             pc;
  PetscMPIInt    size;
  MatSolverType  stype;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckProblem(nep,1);
  if (nep->state) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscLogEventBegin(NEP_SetUp,nep,0,0,0));

  /* reset the convergence flag from the previous solves */
  nep->reason = NEP_CONVERGED_ITERATING;

  /* set default solver type (NEPSetFromOptions was not called) */
  if (!((PetscObject)nep)->type_name) PetscCall(NEPSetType(nep,NEPRII));
  if (nep->useds && !nep->ds) PetscCall(NEPGetDS(nep,&nep->ds));
  if (nep->useds) PetscCall(NEPSetDSType(nep));
  if (!nep->rg) PetscCall(NEPGetRG(nep,&nep->rg));
  if (!((PetscObject)nep->rg)->type_name) PetscCall(RGSetType(nep->rg,RGINTERVAL));

  /* set problem dimensions */
  switch (nep->fui) {
  case NEP_USER_INTERFACE_CALLBACK:
    PetscCall(NEPGetFunction(nep,&T,NULL,NULL,NULL));
    PetscCall(MatGetSize(T,&nep->n,NULL));
    PetscCall(MatGetLocalSize(T,&nep->nloc,NULL));
    break;
  case NEP_USER_INTERFACE_SPLIT:
    PetscCall(MatDuplicate(nep->A[0],MAT_DO_NOT_COPY_VALUES,&nep->function));
    if (nep->P) PetscCall(MatDuplicate(nep->P[0],MAT_DO_NOT_COPY_VALUES,&nep->function_pre));
    PetscCall(MatDuplicate(nep->A[0],MAT_DO_NOT_COPY_VALUES,&nep->jacobian));
    PetscCall(MatGetSize(nep->A[0],&nep->n,NULL));
    PetscCall(MatGetLocalSize(nep->A[0],&nep->nloc,NULL));
    break;
  }

  /* set default problem type */
  if (!nep->problem_type) PetscCall(NEPSetProblemType(nep,NEP_GENERAL));

  /* check consistency of refinement options */
  if (nep->refine) {
    PetscCheck(nep->fui==NEP_USER_INTERFACE_SPLIT,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Iterative refinement only implemented in split form");
    if (!nep->scheme) {  /* set default scheme */
      PetscCall(NEPRefineGetKSP(nep,&ksp));
      PetscCall(KSPGetPC(ksp,&pc));
      PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&flg));
      if (flg) PetscCall(PetscObjectTypeCompareAny((PetscObject)pc,&flg,PCLU,PCCHOLESKY,""));
      nep->scheme = flg? NEP_REFINE_SCHEME_MBE: NEP_REFINE_SCHEME_SCHUR;
    }
    if (nep->scheme==NEP_REFINE_SCHEME_MBE) {
      PetscCall(NEPRefineGetKSP(nep,&ksp));
      PetscCall(KSPGetPC(ksp,&pc));
      PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&flg));
      if (flg) PetscCall(PetscObjectTypeCompareAny((PetscObject)pc,&flg,PCLU,PCCHOLESKY,""));
      PetscCheck(flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The MBE scheme for refinement requires a direct solver in KSP");
      PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pc),&size));
      if (size>1) {   /* currently selected PC is a factorization */
        PetscCall(PCFactorGetMatSolverType(pc,&stype));
        PetscCall(PetscStrcmp(stype,MATSOLVERPETSC,&flg));
        PetscCheck(!flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"For Newton refinement, you chose to solve linear systems with a factorization, but in parallel runs you need to select an external package");
      }
    }
    if (nep->scheme==NEP_REFINE_SCHEME_SCHUR) {
      PetscCheck(nep->npart==1,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The Schur scheme for refinement does not support subcommunicators");
    }
  }
  /* call specific solver setup */
  PetscUseTypeMethod(nep,setup);

  /* set tolerance if not yet set */
  if (nep->tol==(PetscReal)PETSC_DETERMINE) nep->tol = SLEPC_DEFAULT_TOL;
  if (nep->refine) {
    if (nep->rtol==(PetscReal)PETSC_DETERMINE) nep->rtol = PetscMax(nep->tol/1000,PETSC_MACHINE_EPSILON);
    if (nep->rits==(PetscReal)PETSC_DETERMINE) nep->rits = (nep->refine==NEP_REFINE_SIMPLE)? 10: 1;
  }

  /* fill sorting criterion context */
  switch (nep->which) {
    case NEP_LARGEST_MAGNITUDE:
      nep->sc->comparison    = SlepcCompareLargestMagnitude;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_SMALLEST_MAGNITUDE:
      nep->sc->comparison    = SlepcCompareSmallestMagnitude;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_LARGEST_REAL:
      nep->sc->comparison    = SlepcCompareLargestReal;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_SMALLEST_REAL:
      nep->sc->comparison    = SlepcCompareSmallestReal;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_LARGEST_IMAGINARY:
      nep->sc->comparison    = SlepcCompareLargestImaginary;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_SMALLEST_IMAGINARY:
      nep->sc->comparison    = SlepcCompareSmallestImaginary;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_TARGET_MAGNITUDE:
      nep->sc->comparison    = SlepcCompareTargetMagnitude;
      nep->sc->comparisonctx = &nep->target;
      break;
    case NEP_TARGET_REAL:
      nep->sc->comparison    = SlepcCompareTargetReal;
      nep->sc->comparisonctx = &nep->target;
      break;
    case NEP_TARGET_IMAGINARY:
      nep->sc->comparison    = SlepcCompareTargetImaginary;
      nep->sc->comparisonctx = &nep->target;
      break;
    case NEP_ALL:
      nep->sc->comparison    = SlepcCompareSmallestReal;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_WHICH_USER:
      break;
  }

  nep->sc->map    = NULL;
  nep->sc->mapobj = NULL;

  /* fill sorting criterion for DS */
  if (nep->useds) {
    PetscCall(DSGetSlepcSC(nep->ds,&sc));
    sc->comparison    = nep->sc->comparison;
    sc->comparisonctx = nep->sc->comparisonctx;
    PetscCall(PetscObjectTypeCompare((PetscObject)nep,NEPNLEIGS,&flg));
    if (!flg) {
      sc->map    = NULL;
      sc->mapobj = NULL;
    }
  }
  PetscCheck(nep->nev<=nep->ncv,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"nev bigger than ncv");

  /* process initial vectors */
  if (nep->nini<0) {
    k = -nep->nini;
    PetscCheck(k<=nep->ncv,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"The number of initial vectors is larger than ncv");
    PetscCall(BVInsertVecs(nep->V,0,&k,nep->IS,PETSC_TRUE));
    PetscCall(SlepcBasisDestroy_Private(&nep->nini,&nep->IS));
    nep->nini = k;
  }
  PetscCall(PetscLogEventEnd(NEP_SetUp,nep,0,0,0));
  nep->state = NEP_STATE_SETUP;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSetInitialSpace - Specify a basis of vectors that constitute the initial
   space, that is, the subspace from which the solver starts to iterate.

   Collective

   Input Parameters:
+  nep   - the nonlinear eigensolver context
.  n     - number of vectors
-  is    - set of basis vectors of the initial space

   Notes:
   Some solvers start to iterate on a single vector (initial vector). In that case,
   the other vectors are ignored.

   These vectors do not persist from one NEPSolve() call to the other, so the
   initial space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Common usage of this function is when the user can provide a rough approximation
   of the wanted eigenspace. Then, convergence may be faster.

   Level: intermediate

.seealso: NEPSetUp()
@*/
PetscErrorCode NEPSetInitialSpace(NEP nep,PetscInt n,Vec is[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,n,2);
  PetscCheck(n>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative");
  if (n>0) {
    PetscAssertPointer(is,3);
    PetscValidHeaderSpecific(*is,VEC_CLASSID,3);
  }
  PetscCall(SlepcBasisReference_Private(n,is,&nep->nini,&nep->IS));
  if (n>0) nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  NEPSetDimensions_Default - Set reasonable values for ncv, mpd if not set
  by the user. This is called at setup.
 */
PetscErrorCode NEPSetDimensions_Default(NEP nep,PetscInt nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  if (*ncv!=PETSC_DETERMINE) { /* ncv set */
    PetscCheck(*ncv>=nev,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev");
  } else if (*mpd!=PETSC_DETERMINE) { /* mpd set */
    *ncv = PetscMin(nep->n,nev+(*mpd));
  } else { /* neither set: defaults depend on nev being small or large */
    if (nev<500) *ncv = PetscMin(nep->n,PetscMax(2*nev,nev+15));
    else {
      *mpd = 500;
      *ncv = PetscMin(nep->n,nev+(*mpd));
    }
  }
  if (*mpd==PETSC_DETERMINE) *mpd = *ncv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPAllocateSolution - Allocate memory storage for common variables such
   as eigenvalues and eigenvectors.

   Collective

   Input Parameters:
+  nep   - eigensolver context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developer Notes:
   This is SLEPC_EXTERN because it may be required by user plugin NEP
   implementations.

   Level: developer

.seealso: PEPSetUp()
@*/
PetscErrorCode NEPAllocateSolution(NEP nep,PetscInt extra)
{
  PetscInt       oldsize,requested;
  PetscRandom    rand;
  Mat            T;
  Vec            t;

  PetscFunctionBegin;
  requested = nep->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  PetscCall(BVGetSizes(nep->V,NULL,NULL,&oldsize));

  /* allocate space for eigenvalues and friends */
  if (requested != oldsize || !nep->eigr) {
    PetscCall(PetscFree4(nep->eigr,nep->eigi,nep->errest,nep->perm));
    PetscCall(PetscMalloc4(requested,&nep->eigr,requested,&nep->eigi,requested,&nep->errest,requested,&nep->perm));
  }

  /* allocate V */
  if (!nep->V) PetscCall(NEPGetBV(nep,&nep->V));
  if (!oldsize) {
    if (!((PetscObject)nep->V)->type_name) PetscCall(BVSetType(nep->V,BVMAT));
    if (nep->fui==NEP_USER_INTERFACE_SPLIT) T = nep->A[0];
    else PetscCall(NEPGetFunction(nep,&T,NULL,NULL,NULL));
    PetscCall(MatCreateVecsEmpty(T,&t,NULL));
    PetscCall(BVSetSizesFromVec(nep->V,t,requested));
    PetscCall(VecDestroy(&t));
  } else PetscCall(BVResize(nep->V,requested,PETSC_FALSE));

  /* allocate W */
  if (nep->twosided) {
    PetscCall(BVGetRandomContext(nep->V,&rand));  /* make sure the random context is available when duplicating */
    PetscCall(BVDestroy(&nep->W));
    PetscCall(BVDuplicate(nep->V,&nep->W));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
