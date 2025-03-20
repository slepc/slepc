/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   PEP routines related to problem setup
*/

#include <slepc/private/pepimpl.h>       /*I "slepcpep.h" I*/

/*
   Let the solver choose the ST type that should be used by default,
   otherwise set it to SHIFT.
   This is called at PEPSetFromOptions (before STSetFromOptions)
   and also at PEPSetUp (in case PEPSetFromOptions was not called).
*/
PetscErrorCode PEPSetDefaultST(PEP pep)
{
  PetscFunctionBegin;
  PetscTryTypeMethod(pep,setdefaultst);
  if (!((PetscObject)pep->st)->type_name) PetscCall(STSetType(pep->st,STSHIFT));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   This is used in Q-Arnoldi and STOAR to set the transform flag by
   default, otherwise the user has to explicitly run with -st_transform
*/
PetscErrorCode PEPSetDefaultST_Transform(PEP pep)
{
  PetscFunctionBegin;
  PetscCall(STSetTransform(pep->st,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetDSType - Sets the type of the internal DS object based on the current
   settings of the polynomial eigensolver.

   Collective

   Input Parameter:
.  pep - polynomial eigensolver context

   Note:
   This function need not be called explicitly, since it will be called at
   both PEPSetFromOptions() and PEPSetUp().

   Level: developer

.seealso: PEPSetFromOptions(), PEPSetUp()
@*/
PetscErrorCode PEPSetDSType(PEP pep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscTryTypeMethod(pep,setdstype);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetUp - Sets up all the internal data structures necessary for the
   execution of the PEP solver.

   Collective

   Input Parameter:
.  pep   - solver context

   Notes:
   This function need not be called explicitly in most cases, since PEPSolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: PEPCreate(), PEPSolve(), PEPDestroy()
@*/
PetscErrorCode PEPSetUp(PEP pep)
{
  SlepcSC        sc;
  PetscBool      istrivial,flg;
  PetscInt       k;
  KSP            ksp;
  PC             pc;
  PetscMPIInt    size;
  MatSolverType  stype;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (pep->state) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscLogEventBegin(PEP_SetUp,pep,0,0,0));

  /* reset the convergence flag from the previous solves */
  pep->reason = PEP_CONVERGED_ITERATING;

  /* set default solver type (PEPSetFromOptions was not called) */
  if (!((PetscObject)pep)->type_name) PetscCall(PEPSetType(pep,PEPTOAR));
  if (!pep->st) PetscCall(PEPGetST(pep,&pep->st));
  PetscCall(PEPSetDefaultST(pep));
  if (!pep->ds) PetscCall(PEPGetDS(pep,&pep->ds));
  PetscCall(PEPSetDSType(pep));
  if (!pep->rg) PetscCall(PEPGetRG(pep,&pep->rg));
  if (!((PetscObject)pep->rg)->type_name) PetscCall(RGSetType(pep->rg,RGINTERVAL));

  /* check matrices, transfer them to ST */
  PetscCheck(pep->A,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"PEPSetOperators must be called first");
  PetscCall(STSetMatrices(pep->st,pep->nmat,pep->A));

  /* set problem dimensions */
  PetscCall(MatGetSize(pep->A[0],&pep->n,NULL));
  PetscCall(MatGetLocalSize(pep->A[0],&pep->nloc,NULL));

  /* set default problem type */
  if (!pep->problem_type) PetscCall(PEPSetProblemType(pep,PEP_GENERAL));
  if (pep->nev > (pep->nmat-1)*pep->n) pep->nev = (pep->nmat-1)*pep->n;
  if (pep->ncv > (pep->nmat-1)*pep->n) pep->ncv = (pep->nmat-1)*pep->n;

  /* check consistency of refinement options */
  if (pep->refine) {
    if (!pep->scheme) {  /* set default scheme */
      PetscCall(PEPRefineGetKSP(pep,&ksp));
      PetscCall(KSPGetPC(ksp,&pc));
      PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&flg));
      if (flg) PetscCall(PetscObjectTypeCompareAny((PetscObject)pc,&flg,PCLU,PCCHOLESKY,""));
      pep->scheme = flg? PEP_REFINE_SCHEME_MBE: PEP_REFINE_SCHEME_SCHUR;
    }
    if (pep->scheme==PEP_REFINE_SCHEME_MBE) {
      PetscCall(PEPRefineGetKSP(pep,&ksp));
      PetscCall(KSPGetPC(ksp,&pc));
      PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&flg));
      if (flg) PetscCall(PetscObjectTypeCompareAny((PetscObject)pc,&flg,PCLU,PCCHOLESKY,""));
      PetscCheck(flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"The MBE scheme for refinement requires a direct solver in KSP");
      PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pc),&size));
      if (size>1) {   /* currently selected PC is a factorization */
        PetscCall(PCFactorGetMatSolverType(pc,&stype));
        PetscCall(PetscStrcmp(stype,MATSOLVERPETSC,&flg));
        PetscCheck(!flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"For Newton refinement, you chose to solve linear systems with a factorization, but in parallel runs you need to select an external package");
      }
    }
    if (pep->scheme==PEP_REFINE_SCHEME_SCHUR) {
      PetscCheck(pep->npart==1,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"The Schur scheme for refinement does not support subcommunicators");
    }
  }
  /* call specific solver setup */
  PetscUseTypeMethod(pep,setup);

  /* set tolerance if not yet set */
  if (pep->tol==(PetscReal)PETSC_DETERMINE) pep->tol = SLEPC_DEFAULT_TOL;
  if (pep->refine) {
    if (pep->rtol==(PetscReal)PETSC_DETERMINE) pep->rtol = PetscMax(pep->tol/1000,PETSC_MACHINE_EPSILON);
    if (pep->rits==PETSC_DETERMINE) pep->rits = (pep->refine==PEP_REFINE_SIMPLE)? 10: 1;
  }

  /* set default extraction */
  if (!pep->extract) {
    pep->extract = (pep->basis==PEP_BASIS_MONOMIAL)? PEP_EXTRACT_NORM: PEP_EXTRACT_NONE;
  }

  /* fill sorting criterion context */
  switch (pep->which) {
    case PEP_LARGEST_MAGNITUDE:
      pep->sc->comparison    = SlepcCompareLargestMagnitude;
      pep->sc->comparisonctx = NULL;
      break;
    case PEP_SMALLEST_MAGNITUDE:
      pep->sc->comparison    = SlepcCompareSmallestMagnitude;
      pep->sc->comparisonctx = NULL;
      break;
    case PEP_LARGEST_REAL:
      pep->sc->comparison    = SlepcCompareLargestReal;
      pep->sc->comparisonctx = NULL;
      break;
    case PEP_SMALLEST_REAL:
      pep->sc->comparison    = SlepcCompareSmallestReal;
      pep->sc->comparisonctx = NULL;
      break;
    case PEP_LARGEST_IMAGINARY:
      pep->sc->comparison    = SlepcCompareLargestImaginary;
      pep->sc->comparisonctx = NULL;
      break;
    case PEP_SMALLEST_IMAGINARY:
      pep->sc->comparison    = SlepcCompareSmallestImaginary;
      pep->sc->comparisonctx = NULL;
      break;
    case PEP_TARGET_MAGNITUDE:
      pep->sc->comparison    = SlepcCompareTargetMagnitude;
      pep->sc->comparisonctx = &pep->target;
      break;
    case PEP_TARGET_REAL:
      pep->sc->comparison    = SlepcCompareTargetReal;
      pep->sc->comparisonctx = &pep->target;
      break;
    case PEP_TARGET_IMAGINARY:
      pep->sc->comparison    = SlepcCompareTargetImaginary;
      pep->sc->comparisonctx = &pep->target;
      break;
    case PEP_ALL:
      pep->sc->comparison    = SlepcCompareSmallestReal;
      pep->sc->comparisonctx = NULL;
      break;
    case PEP_WHICH_USER:
      break;
  }
  pep->sc->map    = NULL;
  pep->sc->mapobj = NULL;

  /* fill sorting criterion for DS */
  if (pep->which!=PEP_ALL) {
    PetscCall(DSGetSlepcSC(pep->ds,&sc));
    PetscCall(RGIsTrivial(pep->rg,&istrivial));
    sc->rg            = istrivial? NULL: pep->rg;
    sc->comparison    = pep->sc->comparison;
    sc->comparisonctx = pep->sc->comparisonctx;
    sc->map           = SlepcMap_ST;
    sc->mapobj        = (PetscObject)pep->st;
  }
  /* setup ST */
  PetscCall(STSetUp(pep->st));

  /* compute matrix coefficients */
  PetscCall(STGetTransform(pep->st,&flg));
  if (!flg) {
    if (pep->which!=PEP_ALL && pep->solvematcoeffs) PetscCall(STMatSetUp(pep->st,1.0,pep->solvematcoeffs));
  } else {
    PetscCheck(pep->basis==PEP_BASIS_MONOMIAL,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Cannot use ST-transform with non-monomial basis in PEP");
  }

  /* compute scale factor if no set by user */
  PetscCall(PEPComputeScaleFactor(pep));

  /* build balancing matrix if required */
  if (pep->scale==PEP_SCALE_DIAGONAL || pep->scale==PEP_SCALE_BOTH) {
    if (!pep->Dl) PetscCall(BVCreateVec(pep->V,&pep->Dl));
    if (!pep->Dr) PetscCall(BVCreateVec(pep->V,&pep->Dr));
    PetscCall(PEPBuildDiagonalScaling(pep));
  }

  /* process initial vectors */
  if (pep->nini<0) {
    k = -pep->nini;
    PetscCheck(k<=pep->ncv,PetscObjectComm((PetscObject)pep),PETSC_ERR_USER_INPUT,"The number of initial vectors is larger than ncv");
    PetscCall(BVInsertVecs(pep->V,0,&k,pep->IS,PETSC_TRUE));
    PetscCall(SlepcBasisDestroy_Private(&pep->nini,&pep->IS));
    pep->nini = k;
  }
  PetscCall(PetscLogEventEnd(PEP_SetUp,pep,0,0,0));
  pep->state = PEP_STATE_SETUP;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetOperators - Sets the coefficient matrices associated with the polynomial
   eigenvalue problem.

   Collective

   Input Parameters:
+  pep  - the eigenproblem solver context
.  nmat - number of matrices in array A
-  A    - the array of matrices associated with the eigenproblem

   Notes:
   The polynomial eigenproblem is defined as P(l)*x=0, where l is
   the eigenvalue, x is the eigenvector, and P(l) is defined as
   P(l) = A_0 + l*A_1 + ... + l^d*A_d, with d=nmat-1 (the degree of P).
   For non-monomial bases, this expression is different.

   Level: beginner

.seealso: PEPSolve(), PEPGetOperators(), PEPGetNumMatrices(), PEPSetBasis()
@*/
PetscErrorCode PEPSetOperators(PEP pep,PetscInt nmat,Mat A[])
{
  PetscInt       i,n=0,m,m0=0,mloc,nloc,mloc0=0;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,nmat,2);
  PetscCheck(nmat>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Non-positive value of nmat: %" PetscInt_FMT,nmat);
  PetscCheck(nmat>2,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Cannot solve linear eigenproblems with PEP; use EPS instead");
  PetscAssertPointer(A,3);

  for (i=0;i<nmat;i++) {
    PetscValidHeaderSpecific(A[i],MAT_CLASSID,3);
    PetscCheckSameComm(pep,1,A[i],3);
    PetscCall(MatGetSize(A[i],&m,&n));
    PetscCall(MatGetLocalSize(A[i],&mloc,&nloc));
    PetscCheck(m==n,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"A[%" PetscInt_FMT "] is a non-square matrix (%" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",i,m,n);
    PetscCheck(mloc==nloc,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"A[%" PetscInt_FMT "] does not have equal row and column local sizes (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,mloc,nloc);
    if (!i) { m0 = m; mloc0 = mloc; }
    PetscCheck(m==m0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_INCOMP,"Dimensions of A[%" PetscInt_FMT "] do not match with previous matrices (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,m,m0);
    PetscCheck(mloc==mloc0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_INCOMP,"Local dimensions of A[%" PetscInt_FMT "] do not match with previous matrices (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,mloc,mloc0);
    PetscCall(PetscObjectReference((PetscObject)A[i]));
  }

  if (pep->state && (n!=pep->n || nloc!=pep->nloc)) PetscCall(PEPReset(pep));
  else if (pep->nmat) {
    PetscCall(MatDestroyMatrices(pep->nmat,&pep->A));
    PetscCall(PetscFree2(pep->pbc,pep->nrma));
    PetscCall(PetscFree(pep->solvematcoeffs));
  }

  PetscCall(PetscMalloc1(nmat,&pep->A));
  PetscCall(PetscCalloc2(3*nmat,&pep->pbc,nmat,&pep->nrma));
  for (i=0;i<nmat;i++) {
    pep->A[i]   = A[i];
    pep->pbc[i] = 1.0;  /* default to monomial basis */
  }
  pep->nmat = nmat;
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetOperators - Gets the matrices associated with the polynomial eigensystem.

   Not Collective

   Input Parameters:
+  pep - the PEP context
-  k   - the index of the requested matrix (starting in 0)

   Output Parameter:
.  A - the requested matrix

   Level: intermediate

.seealso: PEPSolve(), PEPSetOperators(), PEPGetNumMatrices()
@*/
PetscErrorCode PEPGetOperators(PEP pep,PetscInt k,Mat *A)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(A,3);
  PetscCheck(k>=0 && k<pep->nmat,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,pep->nmat-1);
  *A = pep->A[k];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetNumMatrices - Returns the number of matrices stored in the PEP.

   Not Collective

   Input Parameter:
.  pep - the PEP context

   Output Parameters:
.  nmat - the number of matrices passed in PEPSetOperators()

   Level: intermediate

.seealso: PEPSetOperators()
@*/
PetscErrorCode PEPGetNumMatrices(PEP pep,PetscInt *nmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(nmat,2);
  *nmat = pep->nmat;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSetInitialSpace - Specify a basis of vectors that constitute the initial
   space, that is, the subspace from which the solver starts to iterate.

   Collective

   Input Parameters:
+  pep   - the polynomial eigensolver context
.  n     - number of vectors
-  is    - set of basis vectors of the initial space

   Notes:
   Some solvers start to iterate on a single vector (initial vector). In that case,
   the other vectors are ignored.

   These vectors do not persist from one PEPSolve() call to the other, so the
   initial space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Common usage of this function is when the user can provide a rough approximation
   of the wanted eigenspace. Then, convergence may be faster.

   Level: intermediate

.seealso: PEPSetUp()
@*/
PetscErrorCode PEPSetInitialSpace(PEP pep,PetscInt n,Vec is[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,n,2);
  PetscCheck(n>=0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative");
  if (n>0) {
    PetscAssertPointer(is,3);
    PetscValidHeaderSpecific(*is,VEC_CLASSID,3);
  }
  PetscCall(SlepcBasisReference_Private(n,is,&pep->nini,&pep->IS));
  if (n>0) pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  PEPSetDimensions_Default - Set reasonable values for ncv, mpd if not set
  by the user. This is called at setup.
 */
PetscErrorCode PEPSetDimensions_Default(PEP pep,PetscInt nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscBool      krylov;
  PetscInt       dim;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompareAny((PetscObject)pep,&krylov,PEPTOAR,PEPSTOAR,PEPQARNOLDI,""));
  dim = (pep->nmat-1)*pep->n;
  if (*ncv!=PETSC_DETERMINE) { /* ncv set */
    if (krylov) {
      PetscCheck(*ncv>nev || (*ncv==nev && *ncv==dim),PetscObjectComm((PetscObject)pep),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev+1");
    } else {
      PetscCheck(*ncv>=nev,PetscObjectComm((PetscObject)pep),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev");
    }
  } else if (*mpd!=PETSC_DETERMINE) { /* mpd set */
    *ncv = PetscMin(dim,nev+(*mpd));
  } else { /* neither set: defaults depend on nev being small or large */
    if (nev<500) *ncv = PetscMin(dim,PetscMax(2*nev,nev+15));
    else {
      *mpd = 500;
      *ncv = PetscMin(dim,nev+(*mpd));
    }
  }
  if (*mpd==PETSC_DETERMINE) *mpd = *ncv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPAllocateSolution - Allocate memory storage for common variables such
   as eigenvalues and eigenvectors.

   Collective

   Input Parameters:
+  pep   - eigensolver context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developer Notes:
   This is SLEPC_EXTERN because it may be required by user plugin PEP
   implementations.

   Level: developer

.seealso: PEPSetUp()
@*/
PetscErrorCode PEPAllocateSolution(PEP pep,PetscInt extra)
{
  PetscInt       oldsize,requested,requestedbv;
  Vec            t;

  PetscFunctionBegin;
  requested = (pep->lineariz? pep->ncv: pep->ncv*(pep->nmat-1)) + extra;
  requestedbv = pep->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  PetscCall(BVGetSizes(pep->V,NULL,NULL,&oldsize));

  /* allocate space for eigenvalues and friends */
  if (requested != oldsize || !pep->eigr) {
    PetscCall(PetscFree4(pep->eigr,pep->eigi,pep->errest,pep->perm));
    PetscCall(PetscMalloc4(requested,&pep->eigr,requested,&pep->eigi,requested,&pep->errest,requested,&pep->perm));
  }

  /* allocate V */
  if (!pep->V) PetscCall(PEPGetBV(pep,&pep->V));
  if (!oldsize) {
    if (!((PetscObject)pep->V)->type_name) PetscCall(BVSetType(pep->V,BVMAT));
    PetscCall(STMatCreateVecsEmpty(pep->st,&t,NULL));
    PetscCall(BVSetSizesFromVec(pep->V,t,requestedbv));
    PetscCall(VecDestroy(&t));
  } else PetscCall(BVResize(pep->V,requestedbv,PETSC_FALSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}
