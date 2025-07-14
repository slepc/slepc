/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   The ST interface routines, callable by users
*/

#include <slepc/private/stimpl.h>            /*I "slepcst.h" I*/

PetscClassId     ST_CLASSID = 0;
PetscLogEvent    ST_SetUp = 0,ST_ComputeOperator = 0,ST_Apply = 0,ST_ApplyTranspose = 0,ST_ApplyHermitianTranspose = 0,ST_MatSetUp = 0,ST_MatMult = 0,ST_MatMultTranspose = 0,ST_MatSolve = 0,ST_MatSolveTranspose = 0;
static PetscBool STPackageInitialized = PETSC_FALSE;

const char *STMatModes[] = {"COPY","INPLACE","SHELL","STMatMode","ST_MATMODE_",NULL};

/*@C
   STFinalizePackage - This function destroys everything in the Slepc interface
   to the ST package. It is called from SlepcFinalize().

   Level: developer

.seealso: SlepcFinalize()
@*/
PetscErrorCode STFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&STList));
  STPackageInitialized = PETSC_FALSE;
  STRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   STInitializePackage - This function initializes everything in the ST package.
   It is called from PetscDLLibraryRegister() when using dynamic libraries, and
   on the first call to STCreate() when using static libraries.

   Level: developer

.seealso: SlepcInitialize()
@*/
PetscErrorCode STInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (STPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  STPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("Spectral Transform",&ST_CLASSID));
  /* Register Constructors */
  PetscCall(STRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("STSetUp",ST_CLASSID,&ST_SetUp));
  PetscCall(PetscLogEventRegister("STComputeOperatr",ST_CLASSID,&ST_ComputeOperator));
  PetscCall(PetscLogEventRegister("STApply",ST_CLASSID,&ST_Apply));
  PetscCall(PetscLogEventRegister("STApplyTranspose",ST_CLASSID,&ST_ApplyTranspose));
  PetscCall(PetscLogEventRegister("STApplyHermTrans",ST_CLASSID,&ST_ApplyHermitianTranspose));
  PetscCall(PetscLogEventRegister("STMatSetUp",ST_CLASSID,&ST_MatSetUp));
  PetscCall(PetscLogEventRegister("STMatMult",ST_CLASSID,&ST_MatMult));
  PetscCall(PetscLogEventRegister("STMatMultTranspose",ST_CLASSID,&ST_MatMultTranspose));
  PetscCall(PetscLogEventRegister("STMatSolve",ST_CLASSID,&ST_MatSolve));
  PetscCall(PetscLogEventRegister("STMatSolveTranspose",ST_CLASSID,&ST_MatSolveTranspose));
  /* Process Info */
  classids[0] = ST_CLASSID;
  PetscCall(PetscInfoProcessClass("st",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("st",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(ST_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(STFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STReset - Resets the ST context to the initial state (prior to setup)
   and destroys any allocated Vecs and Mats.

   Collective

   Input Parameter:
.  st - the spectral transformation context

   Level: advanced

.seealso: STDestroy()
@*/
PetscErrorCode STReset(ST st)
{
  PetscFunctionBegin;
  if (st) PetscValidHeaderSpecific(st,ST_CLASSID,1);
  if (!st) PetscFunctionReturn(PETSC_SUCCESS);
  STCheckNotSeized(st,1);
  PetscTryTypeMethod(st,reset);
  if (st->ksp) PetscCall(KSPReset(st->ksp));
  PetscCall(MatDestroyMatrices(PetscMax(2,st->nmat),&st->T));
  PetscCall(MatDestroyMatrices(PetscMax(2,st->nmat),&st->A));
  st->nmat = 0;
  PetscCall(PetscFree(st->Astate));
  PetscCall(MatDestroy(&st->Op));
  PetscCall(MatDestroy(&st->P));
  PetscCall(MatDestroy(&st->Pmat));
  PetscCall(MatDestroyMatrices(st->nsplit,&st->Psplit));
  st->nsplit = 0;
  PetscCall(VecDestroyVecs(st->nwork,&st->work));
  st->nwork = 0;
  PetscCall(VecDestroy(&st->wb));
  PetscCall(VecDestroy(&st->wht));
  PetscCall(VecDestroy(&st->D));
  st->state   = ST_STATE_INITIAL;
  st->opready = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STDestroy - Destroys ST context that was created with STCreate().

   Collective

   Input Parameter:
.  st - the spectral transformation context

   Level: beginner

.seealso: STCreate(), STSetUp()
@*/
PetscErrorCode STDestroy(ST *st)
{
  PetscFunctionBegin;
  if (!*st) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*st,ST_CLASSID,1);
  if (--((PetscObject)*st)->refct > 0) { *st = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscCall(STReset(*st));
  PetscTryTypeMethod(*st,destroy);
  PetscCall(KSPDestroy(&(*st)->ksp));
  PetscCall(PetscHeaderDestroy(st));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STCreate - Creates a spectral transformation context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  newst - location to put the spectral transformation context

   Level: beginner

.seealso: STSetUp(), STApply(), STDestroy(), ST
@*/
PetscErrorCode STCreate(MPI_Comm comm,ST *newst)
{
  ST             st;

  PetscFunctionBegin;
  PetscAssertPointer(newst,2);
  PetscCall(STInitializePackage());
  PetscCall(SlepcHeaderCreate(st,ST_CLASSID,"ST","Spectral Transformation","ST",comm,STDestroy,STView));

  st->A            = NULL;
  st->nmat         = 0;
  st->sigma        = 0.0;
  st->defsigma     = 0.0;
  st->matmode      = ST_MATMODE_COPY;
  st->str          = UNKNOWN_NONZERO_PATTERN;
  st->transform    = PETSC_FALSE;
  st->structured   = PETSC_FALSE;
  st->D            = NULL;
  st->Pmat         = NULL;
  st->Pmat_set     = PETSC_FALSE;
  st->Psplit       = NULL;
  st->nsplit       = 0;
  st->strp         = UNKNOWN_NONZERO_PATTERN;

  st->ksp          = NULL;
  st->usesksp      = PETSC_FALSE;
  st->nwork        = 0;
  st->work         = NULL;
  st->wb           = NULL;
  st->wht          = NULL;
  st->state        = ST_STATE_INITIAL;
  st->Astate       = NULL;
  st->T            = NULL;
  st->Op           = NULL;
  st->opseized     = PETSC_FALSE;
  st->opready      = PETSC_FALSE;
  st->P            = NULL;
  st->M            = NULL;
  st->sigma_set    = PETSC_FALSE;
  st->asymm        = PETSC_FALSE;
  st->aherm        = PETSC_FALSE;
  st->data         = NULL;

  *newst = st;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Checks whether the ST matrices are all symmetric or hermitian.
*/
static inline PetscErrorCode STMatIsSymmetricKnown(ST st,PetscBool *symm,PetscBool *herm)
{
  PetscInt       i;
  PetscBool      sbaij=PETSC_FALSE,set,flg=PETSC_FALSE;

  PetscFunctionBegin;
  /* check if problem matrices are all sbaij */
  for (i=0;i<st->nmat;i++) {
    PetscCall(PetscObjectTypeCompareAny((PetscObject)st->A[i],&sbaij,MATSEQSBAIJ,MATMPISBAIJ,""));
    if (!sbaij) break;
  }
  /* check if user has set the symmetric flag */
  *symm = PETSC_TRUE;
  for (i=0;i<st->nmat;i++) {
    PetscCall(MatIsSymmetricKnown(st->A[i],&set,&flg));
    if (!set || !flg) { *symm = PETSC_FALSE; break; }
  }
  if (sbaij) *symm = PETSC_TRUE;
#if defined(PETSC_USE_COMPLEX)
  /* check if user has set the hermitian flag */
  *herm = PETSC_TRUE;
  for (i=0;i<st->nmat;i++) {
    PetscCall(MatIsHermitianKnown(st->A[i],&set,&flg));
    if (!set || !flg) { *herm = PETSC_FALSE; break; }
  }
#else
  *herm = *symm;
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetMatrices - Sets the matrices associated with the eigenvalue problem.

   Collective

   Input Parameters:
+  st - the spectral transformation context
.  n  - number of matrices in array A
-  A  - the array of matrices associated with the eigensystem

   Notes:
   It must be called before STSetUp(). If it is called again after STSetUp() then
   the ST object is reset.

   Level: intermediate

.seealso: STGetMatrix(), STGetNumMatrices(), STSetUp(), STReset()
@*/
PetscErrorCode STSetMatrices(ST st,PetscInt n,Mat A[])
{
  PetscInt       i;
  PetscBool      same=PETSC_TRUE;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,n,2);
  PetscCheck(n>0,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"Must have one or more matrices, you have %" PetscInt_FMT,n);
  PetscAssertPointer(A,3);
  PetscCheckSameComm(st,1,*A,3);
  STCheckNotSeized(st,1);
  PetscCheck(!st->nsplit || st->nsplit==n,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"The number of matrices must be the same as in STSetSplitPreconditioner()");

  if (st->state) {
    if (n!=st->nmat) same = PETSC_FALSE;
    for (i=0;same&&i<n;i++) {
      if (A[i]!=st->A[i]) same = PETSC_FALSE;
    }
    if (!same) PetscCall(STReset(st));
  } else same = PETSC_FALSE;
  if (!same) {
    PetscCall(MatDestroyMatrices(PetscMax(2,st->nmat),&st->A));
    PetscCall(PetscCalloc1(PetscMax(2,n),&st->A));
    PetscCall(PetscFree(st->Astate));
    PetscCall(PetscMalloc1(PetscMax(2,n),&st->Astate));
  }
  for (i=0;i<n;i++) {
    PetscValidHeaderSpecific(A[i],MAT_CLASSID,3);
    PetscCall(PetscObjectReference((PetscObject)A[i]));
    PetscCall(MatDestroy(&st->A[i]));
    st->A[i] = A[i];
    st->Astate[i] = ((PetscObject)A[i])->state;
  }
  if (n==1) {
    st->A[1] = NULL;
    st->Astate[1] = 0;
  }
  st->nmat = n;
  if (same) st->state = ST_STATE_UPDATED;
  else st->state = ST_STATE_INITIAL;
  PetscCheck(!same || !st->Psplit,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"Support for changing the matrices while using a split preconditioner is not implemented yet");
  st->opready = PETSC_FALSE;
  if (!same) PetscCall(STMatIsSymmetricKnown(st,&st->asymm,&st->aherm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetMatrix - Gets the matrices associated with the original eigensystem.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  k  - the index of the requested matrix (starting in 0)

   Output Parameters:
.  A - the requested matrix

   Level: intermediate

.seealso: STSetMatrices(), STGetNumMatrices()
@*/
PetscErrorCode STGetMatrix(ST st,PetscInt k,Mat *A)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,k,2);
  PetscAssertPointer(A,3);
  STCheckMatrices(st,1);
  PetscCheck(k>=0 && k<st->nmat,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,st->nmat-1);
  PetscCheck(((PetscObject)st->A[k])->state==st->Astate[k],PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"Cannot retrieve original matrices (have been modified)");
  *A = st->A[k];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetMatrixTransformed - Gets the matrices associated with the transformed eigensystem.

   Not Collective

   Input Parameters:
+  st - the spectral transformation context
-  k  - the index of the requested matrix (starting in 0)

   Output Parameters:
.  T - the requested matrix

   Level: developer

.seealso: STGetMatrix(), STGetNumMatrices()
@*/
PetscErrorCode STGetMatrixTransformed(ST st,PetscInt k,Mat *T)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,k,2);
  PetscAssertPointer(T,3);
  STCheckMatrices(st,1);
  PetscCheck(k>=0 && k<st->nmat,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,st->nmat-1);
  PetscCheck(st->T,PetscObjectComm((PetscObject)st),PETSC_ERR_POINTER,"There are no transformed matrices");
  *T = st->T[k];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetNumMatrices - Returns the number of matrices stored in the ST.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameters:
.  n - the number of matrices passed in STSetMatrices()

   Level: intermediate

.seealso: STSetMatrices()
@*/
PetscErrorCode STGetNumMatrices(ST st,PetscInt *n)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(n,2);
  *n = st->nmat;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STResetMatrixState - Resets the stored state of the matrices in the ST.

   Logically Collective

   Input Parameter:
.  st - the spectral transformation context

   Note:
   This is useful in solvers where the user matrices are modified during
   the computation, as in nonlinear inverse iteration. The effect is that
   STGetMatrix() will retrieve the modified matrices as if they were
   the matrices originally provided by the user.

   Level: developer

.seealso: STGetMatrix(), EPSPowerSetNonlinear()
@*/
PetscErrorCode STResetMatrixState(ST st)
{
  PetscInt i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  for (i=0;i<st->nmat;i++) st->Astate[i] = ((PetscObject)st->A[i])->state;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetPreconditionerMat - Sets the matrix to be used to build the preconditioner.

   Collective

   Input Parameters:
+  st  - the spectral transformation context
-  mat - the matrix that will be used in constructing the preconditioner

   Notes:
   This matrix will be passed to the internal KSP object (via the last argument
   of KSPSetOperators()) as the matrix to be used when constructing the preconditioner.
   If no matrix is set or mat is set to NULL, A-sigma*B will be used
   to build the preconditioner, being sigma the value set by STSetShift().

   More precisely, this is relevant for spectral transformations that represent
   a rational matrix function, and use a KSP object for the denominator, called
   K in the description of STGetOperator(). It includes also the STPRECOND case.
   If the user has a good approximation to matrix K that can be used to build a
   cheap preconditioner, it can be passed with this function. Note that it affects
   only the Pmat argument of KSPSetOperators(), not the Amat argument.

   If a preconditioner matrix is set, the default is to use an iterative KSP
   rather than a direct method.

   An alternative to pass an approximation of A-sigma*B with this function is
   to provide approximations of A and B via STSetSplitPreconditioner(). The
   difference is that when sigma changes the preconditioner is recomputed.

   Use NULL to remove a previously set matrix.

   Level: advanced

.seealso: STGetPreconditionerMat(), STSetShift(), STGetOperator(), STSetSplitPreconditioner()
@*/
PetscErrorCode STSetPreconditionerMat(ST st,Mat mat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  if (mat) {
    PetscValidHeaderSpecific(mat,MAT_CLASSID,2);
    PetscCheckSameComm(st,1,mat,2);
  }
  STCheckNotSeized(st,1);
  PetscCheck(!mat || !st->Psplit,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"Cannot call both STSetPreconditionerMat and STSetSplitPreconditioner");
  if (mat) PetscCall(PetscObjectReference((PetscObject)mat));
  PetscCall(MatDestroy(&st->Pmat));
  st->Pmat     = mat;
  st->Pmat_set = mat? PETSC_TRUE: PETSC_FALSE;
  st->state    = ST_STATE_INITIAL;
  st->opready  = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetPreconditionerMat - Returns the matrix previously set by STSetPreconditionerMat().

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  mat - the matrix that will be used in constructing the preconditioner or
   NULL if no matrix was set by STSetPreconditionerMat().

   Level: advanced

.seealso: STSetPreconditionerMat()
@*/
PetscErrorCode STGetPreconditionerMat(ST st,Mat *mat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(mat,2);
  *mat = st->Pmat_set? st->Pmat: NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetSplitPreconditioner - Sets the matrices from which to build the preconditioner
   in split form.

   Collective

   Input Parameters:
+  st     - the spectral transformation context
.  n      - number of matrices
.  Psplit - array of matrices
-  strp   - structure flag for Psplit matrices

   Notes:
   The number of matrices passed here must be the same as in STSetMatrices().

   For linear eigenproblems, the preconditioner matrix is computed as
   Pmat(sigma) = A0-sigma*B0, where A0 and B0 are approximations of A and B
   (the eigenproblem matrices) provided via the Psplit array in this function.
   Compared to STSetPreconditionerMat(), this function allows setting a preconditioner
   in a way that is independent of the shift sigma. Whenever the value of sigma
   changes the preconditioner is recomputed.

   Similarly, for polynomial eigenproblems the matrix for the preconditioner
   is expressed as Pmat(sigma) = sum_i Psplit_i*phi_i(sigma), for i=1,...,n, where
   the phi_i's are the polynomial basis functions.

   The structure flag provides information about the relative nonzero pattern of the
   Psplit_i matrices, in the same way as in STSetMatStructure().

   Use n=0 to reset a previously set split preconditioner.

   Level: advanced

.seealso: STGetSplitPreconditionerTerm(), STGetSplitPreconditionerInfo(), STSetPreconditionerMat(), STSetMatrices(), STSetMatStructure()
@*/
PetscErrorCode STSetSplitPreconditioner(ST st,PetscInt n,Mat Psplit[],MatStructure strp)
{
  PetscInt       i,N=0,M,M0=0,mloc,nloc,mloc0=0;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,n,2);
  PetscCheck(n>=0,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"Negative value of n = %" PetscInt_FMT,n);
  PetscCheck(!n || !st->Pmat_set,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"Cannot call both STSetPreconditionerMat and STSetSplitPreconditioner");
  PetscCheck(!n || !st->nmat || st->nmat==n,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"The number of matrices must be the same as in STSetMatrices()");
  if (n) PetscAssertPointer(Psplit,3);
  PetscValidLogicalCollectiveEnum(st,strp,4);
  STCheckNotSeized(st,1);

  for (i=0;i<n;i++) {
    PetscValidHeaderSpecific(Psplit[i],MAT_CLASSID,3);
    PetscCheckSameComm(st,1,Psplit[i],3);
    PetscCall(MatGetSize(Psplit[i],&M,&N));
    PetscCall(MatGetLocalSize(Psplit[i],&mloc,&nloc));
    PetscCheck(M==N,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONG,"Psplit[%" PetscInt_FMT "] is a non-square matrix (%" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",i,M,N);
    PetscCheck(mloc==nloc,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONG,"Psplit[%" PetscInt_FMT "] does not have equal row and column local sizes (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,mloc,nloc);
    if (!i) { M0 = M; mloc0 = mloc; }
    PetscCheck(M==M0,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_INCOMP,"Dimensions of Psplit[%" PetscInt_FMT "] do not match with previous matrices (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,M,M0);
    PetscCheck(mloc==mloc0,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_INCOMP,"Local dimensions of Psplit[%" PetscInt_FMT "] do not match with previous matrices (%" PetscInt_FMT ", %" PetscInt_FMT ")",i,mloc,mloc0);
    PetscCall(PetscObjectReference((PetscObject)Psplit[i]));
  }

  if (st->Psplit) PetscCall(MatDestroyMatrices(st->nsplit,&st->Psplit));

  /* allocate space and copy matrices */
  if (n) {
    PetscCall(PetscMalloc1(n,&st->Psplit));
    for (i=0;i<n;i++) st->Psplit[i] = Psplit[i];
  }
  st->nsplit = n;
  st->strp   = strp;
  st->state  = ST_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetSplitPreconditionerTerm - Gets the matrices associated with
   the split preconditioner.

   Not Collective

   Input Parameters:
+  st - the spectral transformation context
-  k  - the index of the requested matrix (starting in 0)

   Output Parameter:
.  Psplit - the returned matrix

   Level: advanced

.seealso: STSetSplitPreconditioner(), STGetSplitPreconditionerInfo()
@*/
PetscErrorCode STGetSplitPreconditionerTerm(ST st,PetscInt k,Mat *Psplit)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,k,2);
  PetscAssertPointer(Psplit,3);
  PetscCheck(k>=0 && k<st->nsplit,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,st->nsplit-1);
  PetscCheck(st->Psplit,PetscObjectComm((PetscObject)st),PETSC_ERR_ORDER,"You have not called STSetSplitPreconditioner()");
  *Psplit = st->Psplit[k];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetSplitPreconditionerInfo - Returns the number of matrices of the split
   preconditioner, as well as the structure flag.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameters:
+  n    - the number of matrices passed in STSetSplitPreconditioner()
-  strp - the matrix structure flag passed in STSetSplitPreconditioner()

   Level: advanced

.seealso: STSetSplitPreconditioner(), STGetSplitPreconditionerTerm()
@*/
PetscErrorCode STGetSplitPreconditionerInfo(ST st,PetscInt *n,MatStructure *strp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  if (n)    *n    = st->nsplit;
  if (strp) *strp = st->strp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetShift - Sets the shift associated with the spectral transformation.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  shift - the value of the shift

   Notes:
   In some spectral transformations, changing the shift may have associated
   a lot of work, for example recomputing a factorization.

   This function is normally not directly called by users, since the shift is
   indirectly set by EPSSetTarget().

   Level: intermediate

.seealso: EPSSetTarget(), STGetShift(), STSetDefaultShift()
@*/
PetscErrorCode STSetShift(ST st,PetscScalar shift)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  PetscValidLogicalCollectiveScalar(st,shift,2);
  if (st->sigma != shift) {
    STCheckNotSeized(st,1);
    if (st->state==ST_STATE_SETUP) PetscTryTypeMethod(st,setshift,shift);
    st->sigma = shift;
  }
  st->sigma_set = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetShift - Gets the shift associated with the spectral transformation.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  shift - the value of the shift

   Level: intermediate

.seealso: STSetShift()
@*/
PetscErrorCode STGetShift(ST st,PetscScalar* shift)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(shift,2);
  *shift = st->sigma;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetDefaultShift - Sets the value of the shift that should be employed if
   the user did not specify one.

   Logically Collective

   Input Parameters:
+  st - the spectral transformation context
-  defaultshift - the default value of the shift

   Level: developer

.seealso: STSetShift()
@*/
PetscErrorCode STSetDefaultShift(ST st,PetscScalar defaultshift)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveScalar(st,defaultshift,2);
  if (st->defsigma != defaultshift) {
    st->defsigma = defaultshift;
    st->state    = ST_STATE_INITIAL;
    st->opready  = PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STScaleShift - Multiply the shift with a given factor.

   Logically Collective

   Input Parameters:
+  st     - the spectral transformation context
-  factor - the scaling factor

   Note:
   This function does not update the transformation matrices, as opposed to
   STSetShift().

   Level: developer

.seealso: STSetShift()
@*/
PetscErrorCode STScaleShift(ST st,PetscScalar factor)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveScalar(st,factor,2);
  st->sigma *= factor;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetBalanceMatrix - Sets the diagonal matrix to be used for balancing.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  D  - the diagonal matrix (represented as a vector)

   Notes:
   If this matrix is set, STApply will effectively apply D*OP*D^{-1}. Use NULL
   to reset a previously passed D.

   Balancing is usually set via EPSSetBalance, but the advanced user may use
   this function to bypass the usual balancing methods.

   Level: developer

.seealso: EPSSetBalance(), STApply(), STGetBalanceMatrix()
@*/
PetscErrorCode STSetBalanceMatrix(ST st,Vec D)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  if (st->D == D) PetscFunctionReturn(PETSC_SUCCESS);
  STCheckNotSeized(st,1);
  if (D) {
    PetscValidHeaderSpecific(D,VEC_CLASSID,2);
    PetscCheckSameComm(st,1,D,2);
    PetscCall(PetscObjectReference((PetscObject)D));
  }
  PetscCall(VecDestroy(&st->D));
  st->D = D;
  st->state   = ST_STATE_INITIAL;
  st->opready = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetBalanceMatrix - Gets the balance matrix used by the spectral transformation.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  D  - the diagonal matrix (represented as a vector)

   Note:
   If the matrix was not set, a null pointer will be returned.

   Level: developer

.seealso: STSetBalanceMatrix()
@*/
PetscErrorCode STGetBalanceMatrix(ST st,Vec *D)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(D,2);
  *D = st->D;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatCreateVecs - Get vector(s) compatible with the ST matrices.

   Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameters:
+  right - (optional) vector that the matrix can be multiplied against
-  left  - (optional) vector that the matrix vector product can be stored in

   Level: developer

.seealso: STMatCreateVecsEmpty()
@*/
PetscErrorCode STMatCreateVecs(ST st,Vec *right,Vec *left)
{
  PetscFunctionBegin;
  STCheckMatrices(st,1);
  PetscCall(MatCreateVecs(st->A[0],right,left));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatCreateVecsEmpty - Get vector(s) compatible with the ST matrices, i.e. with the same
   parallel layout, but without internal array.

   Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameters:
+  right - (optional) vector that the matrix can be multiplied against
-  left  - (optional) vector that the matrix vector product can be stored in

   Level: developer

.seealso: STMatCreateVecs(), MatCreateVecsEmpty()
@*/
PetscErrorCode STMatCreateVecsEmpty(ST st,Vec *right,Vec *left)
{
  PetscFunctionBegin;
  STCheckMatrices(st,1);
  PetscCall(MatCreateVecsEmpty(st->A[0],right,left));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatGetSize - Returns the number of rows and columns of the ST matrices.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameters:
+  m - the number of global rows
-  n - the number of global columns

   Level: developer

.seealso: STMatGetLocalSize()
@*/
PetscErrorCode STMatGetSize(ST st,PetscInt *m,PetscInt *n)
{
  PetscFunctionBegin;
  STCheckMatrices(st,1);
  PetscCall(MatGetSize(st->A[0],m,n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatGetLocalSize - Returns the number of local rows and columns of the ST matrices.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameters:
+  m - the number of local rows
-  n - the number of local columns

   Level: developer

.seealso: STMatGetSize()
@*/
PetscErrorCode STMatGetLocalSize(ST st,PetscInt *m,PetscInt *n)
{
  PetscFunctionBegin;
  STCheckMatrices(st,1);
  PetscCall(MatGetLocalSize(st->A[0],m,n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetOptionsPrefix - Sets the prefix used for searching for all
   ST options in the database.

   Logically Collective

   Input Parameters:
+  st     - the spectral transformation context
-  prefix - the prefix string to prepend to all ST option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: STAppendOptionsPrefix(), STGetOptionsPrefix()
@*/
PetscErrorCode STSetOptionsPrefix(ST st,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  if (!st->ksp) PetscCall(STGetKSP(st,&st->ksp));
  PetscCall(KSPSetOptionsPrefix(st->ksp,prefix));
  PetscCall(KSPAppendOptionsPrefix(st->ksp,"st_"));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)st,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STAppendOptionsPrefix - Appends to the prefix used for searching for all
   ST options in the database.

   Logically Collective

   Input Parameters:
+  st     - the spectral transformation context
-  prefix - the prefix string to prepend to all ST option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: STSetOptionsPrefix(), STGetOptionsPrefix()
@*/
PetscErrorCode STAppendOptionsPrefix(ST st,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)st,prefix));
  if (!st->ksp) PetscCall(STGetKSP(st,&st->ksp));
  PetscCall(KSPSetOptionsPrefix(st->ksp,((PetscObject)st)->prefix));
  PetscCall(KSPAppendOptionsPrefix(st->ksp,"st_"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetOptionsPrefix - Gets the prefix used for searching for all
   ST options in the database.

   Not Collective

   Input Parameters:
.  st - the spectral transformation context

   Output Parameters:
.  prefix - pointer to the prefix string used, is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: STSetOptionsPrefix(), STAppendOptionsPrefix()
@*/
PetscErrorCode STGetOptionsPrefix(ST st,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)st,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STView - Prints the ST data structure.

   Collective

   Input Parameters:
+  st - the ST context
-  viewer - optional visualization context

   Note:
   The available visualization contexts include
+     PETSC_VIEWER_STDOUT_SELF - standard output (default)
-     PETSC_VIEWER_STDOUT_WORLD - synchronized standard
         output where only the first processor opens
         the file.  All other processors send their
         data to the first processor to print.

   The user can open an alternative visualization contexts with
   PetscViewerASCIIOpen() (output to a specified file).

   Level: beginner

.seealso: EPSView()
@*/
PetscErrorCode STView(ST st,PetscViewer viewer)
{
  STType         cstr;
  char           str[50];
  PetscBool      isascii,isstring;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)st),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(st,1,viewer,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERSTRING,&isstring));
  if (isascii) {
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)st,viewer));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscTryTypeMethod(st,view,viewer);
    PetscCall(PetscViewerASCIIPopTab(viewer));
    PetscCall(SlepcSNPrintfScalar(str,sizeof(str),st->sigma,PETSC_FALSE));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  shift: %s\n",str));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of matrices: %" PetscInt_FMT "\n",st->nmat));
    switch (st->matmode) {
    case ST_MATMODE_COPY:
      break;
    case ST_MATMODE_INPLACE:
      PetscCall(PetscViewerASCIIPrintf(viewer,"  shifting the matrix and unshifting at exit\n"));
      break;
    case ST_MATMODE_SHELL:
      PetscCall(PetscViewerASCIIPrintf(viewer,"  using a shell matrix\n"));
      break;
    }
    if (st->nmat>1 && st->matmode != ST_MATMODE_SHELL) PetscCall(PetscViewerASCIIPrintf(viewer,"  nonzero pattern of the matrices: %s\n",MatStructures[st->str]));
    if (st->Psplit) PetscCall(PetscViewerASCIIPrintf(viewer,"  using split preconditioner matrices with %s\n",MatStructures[st->strp]));
    if (st->transform && st->nmat>2) PetscCall(PetscViewerASCIIPrintf(viewer,"  computing transformed matrices\n"));
    if (st->structured) PetscCall(PetscViewerASCIIPrintf(viewer,"  exploiting structure in the application of the operator\n"));
  } else if (isstring) {
    PetscCall(STGetType(st,&cstr));
    PetscCall(PetscViewerStringSPrintf(viewer," %-7.7s",cstr));
    PetscTryTypeMethod(st,view,viewer);
  }
  if (st->usesksp) {
    if (!st->ksp) PetscCall(STGetKSP(st,&st->ksp));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(KSPView(st->ksp,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STViewFromOptions - View from options

   Collective

   Input Parameters:
+  st   - the spectral transformation context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: STView(), STCreate()
@*/
PetscErrorCode STViewFromOptions(ST st,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)st,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   STRegister - Adds a method to the spectral transformation package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined transformation
-  function - routine to create method context

   Notes:
   STRegister() may be called multiple times to add several user-defined
   spectral transformations.

   Example Usage:
.vb
    STRegister("my_transform",MyTransformCreate);
.ve

   Then, your spectral transform can be chosen with the procedural interface via
$     STSetType(st,"my_transform")
   or at runtime via the option
$     -st_type my_transform

   Level: advanced

.seealso: STRegisterAll()
@*/
PetscErrorCode STRegister(const char *name,PetscErrorCode (*function)(ST))
{
  PetscFunctionBegin;
  PetscCall(STInitializePackage());
  PetscCall(PetscFunctionListAdd(&STList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}
