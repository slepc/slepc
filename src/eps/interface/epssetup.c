/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   EPS routines related to problem setup
*/

#include <slepc/private/epsimpl.h>       /*I "slepceps.h" I*/

/*
   Let the solver choose the ST type that should be used by default,
   otherwise set it to SHIFT.
   This is called at EPSSetFromOptions (before STSetFromOptions)
   and also at EPSSetUp (in case EPSSetFromOptions was not called).
*/
PetscErrorCode EPSSetDefaultST(EPS eps)
{
  PetscFunctionBegin;
  PetscTryTypeMethod(eps,setdefaultst);
  if (!((PetscObject)eps->st)->type_name) PetscCall(STSetType(eps->st,STSHIFT));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   This is done by preconditioned eigensolvers that use the PC only.
   It sets STPRECOND with KSPPREONLY.
*/
PetscErrorCode EPSSetDefaultST_Precond(EPS eps)
{
  KSP            ksp;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) PetscCall(STSetType(eps->st,STPRECOND));
  PetscCall(STGetKSP(eps->st,&ksp));
  if (!((PetscObject)ksp)->type_name) PetscCall(KSPSetType(ksp,KSPPREONLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   This is done by preconditioned eigensolvers that can also use the KSP.
   It sets STPRECOND with the default KSP (GMRES) and maxit=5.
*/
PetscErrorCode EPSSetDefaultST_GMRES(EPS eps)
{
  KSP            ksp;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) PetscCall(STSetType(eps->st,STPRECOND));
  PetscCall(STPrecondSetKSPHasMat(eps->st,PETSC_TRUE));
  PetscCall(STGetKSP(eps->st,&ksp));
  if (!((PetscObject)ksp)->type_name) {
    PetscCall(KSPSetType(ksp,KSPGMRES));
    PetscCall(KSPSetTolerances(ksp,PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT,5));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(SLEPC_HAVE_SCALAPACK) || defined(SLEPC_HAVE_ELPA) || defined(SLEPC_HAVE_ELEMENTAL) || defined(SLEPC_HAVE_EVSL)
/*
   This is for direct eigensolvers that work with A and B directly, so
   no need to factorize B.
*/
PetscErrorCode EPSSetDefaultST_NoFactor(EPS eps)
{
  KSP            ksp;
  PC             pc;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) PetscCall(STSetType(eps->st,STSHIFT));
  PetscCall(STGetKSP(eps->st,&ksp));
  if (!((PetscObject)ksp)->type_name) PetscCall(KSPSetType(ksp,KSPPREONLY));
  PetscCall(KSPGetPC(ksp,&pc));
  if (!((PetscObject)pc)->type_name) PetscCall(PCSetType(pc,PCNONE));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

/*
   Check that the ST selected by the user is compatible with the EPS solver and options
*/
static PetscErrorCode EPSCheckCompatibleST(EPS eps)
{
  PetscBool      precond,shift,sinvert,cayley,lyapii;
#if defined(PETSC_USE_COMPLEX)
  PetscScalar    sigma;
#endif

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STPRECOND,&precond));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&shift));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&sinvert));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STCAYLEY,&cayley));

  /* preconditioned eigensolvers */
  PetscCheck(eps->categ!=EPS_CATEGORY_PRECOND || precond,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver requires ST=PRECOND");
  PetscCheck(eps->categ==EPS_CATEGORY_PRECOND || !precond,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"STPRECOND is intended for preconditioned eigensolvers only");

  /* harmonic extraction */
  PetscCheck(precond || shift || !eps->extraction || eps->extraction==EPS_RITZ,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Cannot use a spectral transformation combined with harmonic extraction");

  /* real shifts in Hermitian problems */
#if defined(PETSC_USE_COMPLEX)
  PetscCall(STGetShift(eps->st,&sigma));
  PetscCheck(!eps->ishermitian || PetscImaginaryPart(sigma)==0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Hermitian problems are not compatible with complex shifts");
#endif

  /* Cayley with PGNHEP */
  PetscCheck(!cayley || eps->problem_type!=EPS_PGNHEP,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Cayley spectral transformation is not compatible with PGNHEP");

  /* make sure that the user does not specify smallest magnitude with shift-and-invert */
  if ((cayley || sinvert) && (eps->categ==EPS_CATEGORY_KRYLOV || eps->categ==EPS_CATEGORY_OTHER)) {
    PetscCall(PetscObjectTypeCompare((PetscObject)eps,EPSLYAPII,&lyapii));
    PetscCheck(lyapii || eps->which==EPS_TARGET_MAGNITUDE || eps->which==EPS_TARGET_REAL || eps->which==EPS_TARGET_IMAGINARY || eps->which==EPS_ALL || eps->which==EPS_WHICH_USER,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"Shift-and-invert requires a target 'which' (see EPSSetWhichEigenpairs), for instance -st_type sinvert -eps_target 0 -eps_target_magnitude");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   MatEstimateSpectralRange_EPS: estimate the spectral range [left,right] of a
   symmetric/Hermitian matrix A using an auxiliary EPS object
*/
PetscErrorCode MatEstimateSpectralRange_EPS(Mat A,PetscReal *left,PetscReal *right)
{
  PetscInt       nconv;
  PetscScalar    eig0;
  PetscReal      tol=1e-3,errest=tol;
  EPS            eps;

  PetscFunctionBegin;
  *left = 0.0; *right = 0.0;
  PetscCall(EPSCreate(PetscObjectComm((PetscObject)A),&eps));
  PetscCall(EPSSetOptionsPrefix(eps,"eps_filter_"));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetTolerances(eps,tol,50));
  PetscCall(EPSSetConvergenceTest(eps,EPS_CONV_ABS));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));
  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetConverged(eps,&nconv));
  if (nconv>0) {
    PetscCall(EPSGetEigenvalue(eps,0,&eig0,NULL));
    PetscCall(EPSGetErrorEstimate(eps,0,&errest));
  } else eig0 = eps->eigr[0];
  *left = PetscRealPart(eig0)-errest;
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetConverged(eps,&nconv));
  if (nconv>0) {
    PetscCall(EPSGetEigenvalue(eps,0,&eig0,NULL));
    PetscCall(EPSGetErrorEstimate(eps,0,&errest));
  } else eig0 = eps->eigr[0];
  *right = PetscRealPart(eig0)+errest;
  PetscCall(EPSDestroy(&eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSSetUpSort_Basic: configure the EPS sorting criterion according to 'which'
*/
PetscErrorCode EPSSetUpSort_Basic(EPS eps)
{
  PetscFunctionBegin;
  switch (eps->which) {
    case EPS_LARGEST_MAGNITUDE:
      eps->sc->comparison    = SlepcCompareLargestMagnitude;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_SMALLEST_MAGNITUDE:
      eps->sc->comparison    = SlepcCompareSmallestMagnitude;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_LARGEST_REAL:
      eps->sc->comparison    = SlepcCompareLargestReal;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_SMALLEST_REAL:
      eps->sc->comparison    = SlepcCompareSmallestReal;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_LARGEST_IMAGINARY:
      eps->sc->comparison    = SlepcCompareLargestImaginary;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_SMALLEST_IMAGINARY:
      eps->sc->comparison    = SlepcCompareSmallestImaginary;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_TARGET_MAGNITUDE:
      eps->sc->comparison    = SlepcCompareTargetMagnitude;
      eps->sc->comparisonctx = &eps->target;
      break;
    case EPS_TARGET_REAL:
      eps->sc->comparison    = SlepcCompareTargetReal;
      eps->sc->comparisonctx = &eps->target;
      break;
    case EPS_TARGET_IMAGINARY:
#if defined(PETSC_USE_COMPLEX)
      eps->sc->comparison    = SlepcCompareTargetImaginary;
      eps->sc->comparisonctx = &eps->target;
#endif
      break;
    case EPS_ALL:
      eps->sc->comparison    = SlepcCompareSmallestReal;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_WHICH_USER:
      break;
  }
  eps->sc->map    = NULL;
  eps->sc->mapobj = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSSetUpSort_Default: configure both EPS and DS sorting criterion
*/
PetscErrorCode EPSSetUpSort_Default(EPS eps)
{
  SlepcSC        sc;
  PetscBool      istrivial;

  PetscFunctionBegin;
  /* fill sorting criterion context */
  PetscCall(EPSSetUpSort_Basic(eps));
  /* fill sorting criterion for DS */
  PetscCall(DSGetSlepcSC(eps->ds,&sc));
  PetscCall(RGIsTrivial(eps->rg,&istrivial));
  sc->rg            = istrivial? NULL: eps->rg;
  sc->comparison    = eps->sc->comparison;
  sc->comparisonctx = eps->sc->comparisonctx;
  sc->map           = SlepcMap_ST;
  sc->mapobj        = (PetscObject)eps->st;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetDSType - Sets the type of the internal DS object based on the current
   settings of the eigenvalue solver.

   Collective

   Input Parameter:
.  eps - eigenproblem solver context

   Note:
   This function need not be called explicitly, since it will be called at
   both EPSSetFromOptions() and EPSSetUp().

   Level: developer

.seealso: EPSSetFromOptions(), EPSSetUp()
@*/
PetscErrorCode EPSSetDSType(EPS eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscTryTypeMethod(eps,setdstype);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetUp - Sets up all the internal data structures necessary for the
   execution of the eigensolver. Then calls STSetUp() for any set-up
   operations associated to the ST object.

   Collective

   Input Parameter:
.  eps   - eigenproblem solver context

   Notes:
   This function need not be called explicitly in most cases, since EPSSolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: EPSCreate(), EPSSolve(), EPSDestroy(), STSetUp(), EPSSetInitialSpace()
@*/
PetscErrorCode EPSSetUp(EPS eps)
{
  Mat            A,B;
  PetscInt       k,nmat;
  PetscBool      flg;
  EPSStoppingCtx ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (eps->state) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscLogEventBegin(EPS_SetUp,eps,0,0,0));

  /* reset the convergence flag from the previous solves */
  eps->reason = EPS_CONVERGED_ITERATING;

  /* Set default solver type (EPSSetFromOptions was not called) */
  if (!((PetscObject)eps)->type_name) PetscCall(EPSSetType(eps,EPSKRYLOVSCHUR));
  if (!eps->st) PetscCall(EPSGetST(eps,&eps->st));
  PetscCall(EPSSetDefaultST(eps));

  PetscCall(STSetTransform(eps->st,PETSC_TRUE));
  PetscCall(STSetStructured(eps->st,PETSC_FALSE));
  if (eps->useds && !eps->ds) PetscCall(EPSGetDS(eps,&eps->ds));
  if (eps->useds) PetscCall(EPSSetDSType(eps));
  if (eps->twosided) {
    PetscCheck(!eps->ishermitian || (eps->isgeneralized && !eps->ispositive),PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Two-sided methods are not intended for %s problems",SLEPC_STRING_HERMITIAN);
  }
  if (!eps->rg) PetscCall(EPSGetRG(eps,&eps->rg));
  if (!((PetscObject)eps->rg)->type_name) PetscCall(RGSetType(eps->rg,RGINTERVAL));

  /* Set problem dimensions */
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCheck(nmat,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"EPSSetOperators must be called first");
  PetscCall(STMatGetSize(eps->st,&eps->n,NULL));
  PetscCall(STMatGetLocalSize(eps->st,&eps->nloc,NULL));

  /* Set default problem type */
  if (!eps->problem_type) {
    if (nmat==1) PetscCall(EPSSetProblemType(eps,EPS_NHEP));
    else PetscCall(EPSSetProblemType(eps,EPS_GNHEP));
  } else if (nmat==1 && eps->isgeneralized) {
    PetscCall(PetscInfo(eps,"Eigenproblem set as generalized but no matrix B was provided; reverting to a standard eigenproblem\n"));
    eps->isgeneralized = PETSC_FALSE;
    eps->problem_type = eps->ishermitian? EPS_HEP: EPS_NHEP;
  } else PetscCheck(nmat==1 || eps->isgeneralized,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_INCOMP,"Inconsistent EPS state: the problem type does not match the number of matrices");

  if (eps->isstructured) {
    /* make sure the user has set the appropriate matrix */
    PetscCall(STGetMatrix(eps->st,0,&A));
    if (eps->problem_type==EPS_BSE) PetscCall(SlepcCheckMatStruct(A,SLEPC_MAT_STRUCT_BSE,NULL));
  }

  /* safeguard for small problems */
  if (eps->isstructured) {
    if (2*eps->nev > eps->n) eps->nev = eps->n/2;
    if (2*eps->ncv > eps->n) eps->ncv = eps->n/2;
  } else {
    if (eps->nev > eps->n) eps->nev = eps->n;
    if (eps->ncv > eps->n) eps->ncv = eps->n;
  }

  /* check some combinations of eps->which */
  PetscCheck(!eps->ishermitian || (eps->isgeneralized && !eps->ispositive) || (eps->which!=EPS_LARGEST_IMAGINARY && eps->which!=EPS_SMALLEST_IMAGINARY && eps->which!=EPS_TARGET_IMAGINARY),PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Sorting the eigenvalues along the imaginary axis is not allowed when all eigenvalues are real");

  /* initialization of matrix norms */
  if (eps->conv==EPS_CONV_NORM) {
    if (!eps->nrma) {
      PetscCall(STGetMatrix(eps->st,0,&A));
      PetscCall(MatNorm(A,NORM_INFINITY,&eps->nrma));
    }
    if (nmat>1 && !eps->nrmb) {
      PetscCall(STGetMatrix(eps->st,1,&B));
      PetscCall(MatNorm(B,NORM_INFINITY,&eps->nrmb));
    }
  }

  /* call specific solver setup */
  PetscUseTypeMethod(eps,setup);

  /* threshold stopping test */
  if (eps->stop==EPS_STOP_THRESHOLD) {
    PetscCheck(eps->thres!=PETSC_MIN_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"Must give a threshold value with EPSSetThreshold()");
    PetscCheck(eps->which==EPS_LARGEST_MAGNITUDE || eps->which==EPS_SMALLEST_MAGNITUDE || eps->which==EPS_LARGEST_REAL || eps->which==EPS_SMALLEST_REAL || eps->which==EPS_TARGET_MAGNITUDE,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Threshold stopping test can only be used with largest/smallest/target magnitude or largest/smallest real selection of eigenvalues");
    if (eps->which==EPS_LARGEST_REAL || eps->which==EPS_SMALLEST_REAL) PetscCheck(eps->problem_type==EPS_HEP || eps->problem_type==EPS_GHEP || eps->problem_type==EPS_BSE,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Threshold stopping test with largest/smallest real can only be used in problems that have all eigenvaues real");
    else PetscCheck(eps->thres>0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"In case of largest/smallest/target magnitude the threshold value must be positive");
    PetscCheck(eps->which==EPS_LARGEST_MAGNITUDE || eps->which==EPS_TARGET_MAGNITUDE || !eps->threlative,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"Can only use a relative threshold with largest/target magnitude selection of eigenvalues");
    PetscCall(PetscNew(&ctx));
    ctx->thres      = eps->thres;
    ctx->threlative = eps->threlative;
    ctx->which      = eps->which;
    PetscCall(EPSSetStoppingTestFunction(eps,EPSStoppingThreshold,ctx,PetscCtxDestroyDefault));
  }

  /* if purification is set, check that it really makes sense */
  if (eps->purify) {
    if (eps->categ==EPS_CATEGORY_PRECOND || eps->categ==EPS_CATEGORY_CONTOUR) eps->purify = PETSC_FALSE;
    else {
      if (!eps->isgeneralized) eps->purify = PETSC_FALSE;
      else if (!eps->ishermitian && !eps->ispositive) eps->purify = PETSC_FALSE;
      else {
        PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STCAYLEY,&flg));
        if (flg) eps->purify = PETSC_FALSE;
      }
    }
  }

  /* set tolerance if not yet set */
  if (eps->tol==(PetscReal)PETSC_DETERMINE) eps->tol = SLEPC_DEFAULT_TOL;

  /* set up sorting criterion */
  PetscTryTypeMethod(eps,setupsort);

  /* Build balancing matrix if required */
  if (eps->balance!=EPS_BALANCE_USER) {
    PetscCall(STSetBalanceMatrix(eps->st,NULL));
    if (!eps->ishermitian && (eps->balance==EPS_BALANCE_ONESIDE || eps->balance==EPS_BALANCE_TWOSIDE)) {
      if (!eps->D) PetscCall(BVCreateVec(eps->V,&eps->D));
      PetscCall(EPSBuildBalance_Krylov(eps));
      PetscCall(STSetBalanceMatrix(eps->st,eps->D));
    }
  }

  /* Setup ST */
  PetscCall(STSetUp(eps->st));
  PetscCall(EPSCheckCompatibleST(eps));

  /* process deflation and initial vectors */
  if (eps->nds<0) {
    PetscCheck(!eps->isstructured,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Deflation space is not supported in structured eigenproblems");
    k = -eps->nds;
    PetscCall(BVInsertConstraints(eps->V,&k,eps->defl));
    PetscCall(SlepcBasisDestroy_Private(&eps->nds,&eps->defl));
    eps->nds = k;
    PetscCall(STCheckNullSpace(eps->st,eps->V));
  }
  if (eps->nini<0) {
    k = -eps->nini;
    PetscCheck(k<=eps->ncv,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The number of initial vectors is larger than ncv");
    PetscCall(BVInsertVecs(eps->V,0,&k,eps->IS,PETSC_TRUE));
    PetscCall(SlepcBasisDestroy_Private(&eps->nini,&eps->IS));
    eps->nini = k;
  }
  if (eps->twosided && eps->ninil<0) {
    k = -eps->ninil;
    PetscCheck(k<=eps->ncv,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The number of left initial vectors is larger than ncv");
    PetscCall(BVInsertVecs(eps->W,0,&k,eps->ISL,PETSC_TRUE));
    PetscCall(SlepcBasisDestroy_Private(&eps->ninil,&eps->ISL));
    eps->ninil = k;
  }

  PetscCall(PetscLogEventEnd(EPS_SetUp,eps,0,0,0));
  eps->state = EPS_STATE_SETUP;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetOperators - Sets the matrices associated with the eigenvalue problem.

   Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  A  - the matrix associated with the eigensystem
-  B  - the second matrix in the case of generalized eigenproblems

   Notes:
   To specify a standard eigenproblem, use NULL for parameter B.

   It must be called before EPSSetUp(). If it is called again after EPSSetUp() and
   the matrix sizes have changed then the EPS object is reset.

   For structured eigenproblem types such as EPS_BSE (see EPSSetProblemType()), the
   provided matrices must have been created with the corresponding helper function,
   i.e., MatCreateBSE().

   Level: beginner

.seealso: EPSSolve(), EPSSetUp(), EPSReset(), EPSGetST(), STGetMatrix(), EPSSetProblemType()
@*/
PetscErrorCode EPSSetOperators(EPS eps,Mat A,Mat B)
{
  PetscInt       m,n,m0,mloc,nloc,mloc0,nmat;
  Mat            mat[2];

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  if (B) PetscValidHeaderSpecific(B,MAT_CLASSID,3);
  PetscCheckSameComm(eps,1,A,2);
  if (B) PetscCheckSameComm(eps,1,B,3);

  /* Check matrix sizes */
  PetscCall(MatGetSize(A,&m,&n));
  PetscCall(MatGetLocalSize(A,&mloc,&nloc));
  PetscCheck(m==n,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"A is a non-square matrix (%" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",m,n);
  PetscCheck(mloc==nloc,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"A does not have equal row and column sizes (%" PetscInt_FMT ", %" PetscInt_FMT ")",mloc,nloc);
  if (B) {
    PetscCall(MatGetSize(B,&m0,&n));
    PetscCall(MatGetLocalSize(B,&mloc0,&nloc));
    PetscCheck(m0==n,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"B is a non-square matrix (%" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",m0,n);
    PetscCheck(mloc0==nloc,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"B does not have equal row and column local sizes (%" PetscInt_FMT ", %" PetscInt_FMT ")",mloc0,nloc);
    PetscCheck(m==m0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_INCOMP,"Dimensions of A and B do not match (%" PetscInt_FMT ", %" PetscInt_FMT ")",m,m0);
    PetscCheck(mloc==mloc0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_INCOMP,"Local dimensions of A and B do not match (%" PetscInt_FMT ", %" PetscInt_FMT ")",mloc,mloc0);
  }
  if (eps->state && (n!=eps->n || nloc!=eps->nloc)) PetscCall(EPSReset(eps));
  eps->nrma = 0.0;
  eps->nrmb = 0.0;
  if (!eps->st) PetscCall(EPSGetST(eps,&eps->st));
  mat[0] = A;
  if (B) {
    mat[1] = B;
    nmat = 2;
  } else nmat = 1;
  PetscCall(STSetMatrices(eps->st,nmat,mat));
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetOperators - Gets the matrices associated with the eigensystem.

   Collective

   Input Parameter:
.  eps - the EPS context

   Output Parameters:
+  A  - the matrix associated with the eigensystem
-  B  - the second matrix in the case of generalized eigenproblems

   Note:
   Does not increase the reference count of the matrices, so you should not destroy them.

   Level: intermediate

.seealso: EPSSolve(), EPSGetST(), STGetMatrix(), STSetMatrices()
@*/
PetscErrorCode EPSGetOperators(EPS eps,Mat *A,Mat *B)
{
  ST             st;
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STGetNumMatrices(st,&k));
  if (A) {
    if (k<1) *A = NULL;
    else PetscCall(STGetMatrix(st,0,A));
  }
  if (B) {
    if (k<2) *B = NULL;
    else PetscCall(STGetMatrix(st,1,B));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetDeflationSpace - Specify a basis of vectors that constitute the deflation
   space.

   Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  n   - number of vectors
-  v   - set of basis vectors of the deflation space

   Notes:
   When a deflation space is given, the eigensolver seeks the eigensolution
   in the restriction of the problem to the orthogonal complement of this
   space. This can be used for instance in the case that an invariant
   subspace is known beforehand (such as the nullspace of the matrix).

   These vectors do not persist from one EPSSolve() call to the other, so the
   deflation space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Level: intermediate

.seealso: EPSSetInitialSpace()
@*/
PetscErrorCode EPSSetDeflationSpace(EPS eps,PetscInt n,Vec v[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,n,2);
  PetscCheck(n>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative");
  if (n>0) {
    PetscAssertPointer(v,3);
    PetscValidHeaderSpecific(*v,VEC_CLASSID,3);
  }
  PetscCall(SlepcBasisReference_Private(n,v,&eps->nds,&eps->defl));
  if (n>0) eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetInitialSpace - Specify a basis of vectors that constitute the initial
   space, that is, the subspace from which the solver starts to iterate.

   Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  n   - number of vectors
-  is  - set of basis vectors of the initial space

   Notes:
   Some solvers start to iterate on a single vector (initial vector). In that case,
   the other vectors are ignored.

   These vectors do not persist from one EPSSolve() call to the other, so the
   initial space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Common usage of this function is when the user can provide a rough approximation
   of the wanted eigenspace. Then, convergence may be faster.

   Level: intermediate

.seealso: EPSSetLeftInitialSpace(), EPSSetDeflationSpace()
@*/
PetscErrorCode EPSSetInitialSpace(EPS eps,PetscInt n,Vec is[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,n,2);
  PetscCheck(n>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative");
  if (n>0) {
    PetscAssertPointer(is,3);
    PetscValidHeaderSpecific(*is,VEC_CLASSID,3);
  }
  PetscCall(SlepcBasisReference_Private(n,is,&eps->nini,&eps->IS));
  if (n>0) eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetLeftInitialSpace - Specify a basis of vectors that constitute the left
   initial space, used by two-sided solvers to start the left subspace.

   Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  n   - number of vectors
-  isl - set of basis vectors of the left initial space

   Notes:
   Left initial vectors are used to initiate the left search space in two-sided
   eigensolvers. Users should pass here an approximation of the left eigenspace,
   if available.

   The same comments in EPSSetInitialSpace() are applicable here.

   Level: intermediate

.seealso: EPSSetInitialSpace(), EPSSetTwoSided()
@*/
PetscErrorCode EPSSetLeftInitialSpace(EPS eps,PetscInt n,Vec isl[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,n,2);
  PetscCheck(n>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative");
  if (n>0) {
    PetscAssertPointer(isl,3);
    PetscValidHeaderSpecific(*isl,VEC_CLASSID,3);
  }
  PetscCall(SlepcBasisReference_Private(n,isl,&eps->ninil,&eps->ISL));
  if (n>0) eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSSetDimensions_Default - Set reasonable values for ncv, mpd if not set
  by the user. This is called at setup.
 */
PetscErrorCode EPSSetDimensions_Default(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscBool      krylov;
  PetscInt       n = eps->isstructured? eps->n/2: eps->n;

  PetscFunctionBegin;
  if (*nev==0 && eps->stop!=EPS_STOP_THRESHOLD) *nev = 1;
  if (*ncv!=PETSC_DETERMINE) { /* ncv set */
    PetscCall(PetscObjectTypeCompareAny((PetscObject)eps,&krylov,EPSKRYLOVSCHUR,EPSARNOLDI,EPSLANCZOS,""));
    if (krylov) {
      PetscCheck(*ncv>=*nev+1 || (*ncv==*nev && *ncv==n),PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev+1");
    } else {
      PetscCheck(*ncv>=*nev,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev");
    }
  } else if (*mpd!=PETSC_DETERMINE) { /* mpd set */
    *ncv = PetscMin(n,*nev+(*mpd));
  } else { /* neither set: defaults depend on nev being small or large */
    if (*nev<500) *ncv = PetscMin(n,PetscMax(2*(*nev),*nev+15));
    else {
      *mpd = 500;
      *ncv = PetscMin(n,*nev+(*mpd));
    }
  }
  if (*mpd==PETSC_DETERMINE) *mpd = *ncv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSAllocateSolution - Allocate memory storage for common variables such
   as eigenvalues and eigenvectors.

   Collective

   Input Parameters:
+  eps   - eigensolver context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developer Notes:
   This is SLEPC_EXTERN because it may be required by user plugin EPS
   implementations.

   Level: developer

.seealso: EPSSetUp()
@*/
PetscErrorCode EPSAllocateSolution(EPS eps,PetscInt extra)
{
  PetscInt       oldsize,requested;
  PetscRandom    rand;
  Vec            t;

  PetscFunctionBegin;
  requested = eps->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  PetscCall(BVGetSizes(eps->V,NULL,NULL,&oldsize));

  /* allocate space for eigenvalues and friends */
  if (requested != oldsize || !eps->eigr) {
    PetscCall(PetscFree4(eps->eigr,eps->eigi,eps->errest,eps->perm));
    PetscCall(PetscMalloc4(requested,&eps->eigr,requested,&eps->eigi,requested,&eps->errest,requested,&eps->perm));
  }

  /* workspace for the case of arbitrary selection */
  if (eps->arbitrary) {
    if (eps->rr) PetscCall(PetscFree2(eps->rr,eps->ri));
    PetscCall(PetscMalloc2(requested,&eps->rr,requested,&eps->ri));
  }

  /* allocate V */
  if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
  if (!oldsize) {
    if (!((PetscObject)eps->V)->type_name) PetscCall(BVSetType(eps->V,BVMAT));
    PetscCall(STMatCreateVecsEmpty(eps->st,&t,NULL));
    PetscCall(BVSetSizesFromVec(eps->V,t,requested));
    PetscCall(VecDestroy(&t));
  } else PetscCall(BVResize(eps->V,requested,PETSC_FALSE));

  /* allocate W */
  if (eps->twosided) {
    PetscCall(BVGetRandomContext(eps->V,&rand));  /* make sure the random context is available when duplicating */
    PetscCall(BVDestroy(&eps->W));
    PetscCall(BVDuplicate(eps->V,&eps->W));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSReallocateSolution - Reallocate memory storage for common variables such
   as the eigenvalues and the basis vectors.

   Collective

   Input Parameters:
+  eps     - eigensolver context
-  newsize - new size

   Developer Notes:
   This is SLEPC_EXTERN because it may be required by user plugin EPS
   implementations.

   This is called during the iteration in case the threshold stopping test has
   been selected.

   Level: developer

.seealso: EPSAllocateSolution(), EPSSetThreshold()
@*/
PetscErrorCode EPSReallocateSolution(EPS eps,PetscInt newsize)
{
  PetscInt    oldsize,*nperm;
  PetscReal   *nerrest;
  PetscScalar *neigr,*neigi;

  PetscFunctionBegin;
  PetscCall(BVGetSizes(eps->V,NULL,NULL,&oldsize));
  if (oldsize>=newsize) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscInfo(eps,"Reallocating basis vectors to %" PetscInt_FMT " columns\n",newsize));

  /* reallocate eigenvalues */
  PetscCall(PetscMalloc4(newsize,&neigr,newsize,&neigi,newsize,&nerrest,newsize,&nperm));
  PetscCall(PetscArraycpy(neigr,eps->eigr,oldsize));
  PetscCall(PetscArraycpy(neigi,eps->eigi,oldsize));
  PetscCall(PetscArraycpy(nerrest,eps->errest,oldsize));
  PetscCall(PetscArraycpy(nperm,eps->perm,oldsize));
  PetscCall(PetscFree4(eps->eigr,eps->eigi,eps->errest,eps->perm));
  eps->eigr   = neigr;
  eps->eigi   = neigi;
  eps->errest = nerrest;
  eps->perm   = nperm;
  /* reallocate V,W */
  PetscCall(BVResize(eps->V,newsize,PETSC_TRUE));
  if (eps->twosided) PetscCall(BVResize(eps->W,newsize,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}
