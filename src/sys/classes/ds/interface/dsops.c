/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   DS operations: DSSolve(), DSVectors(), etc
*/

#include <slepc/private/dsimpl.h>      /*I "slepcds.h" I*/

/*@
   DSGetLeadingDimension - Returns the leading dimension of the allocated
   matrices.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  ld - leading dimension (maximum allowed dimension for the matrices)

   Level: advanced

.seealso: [](sec:ds), `DSAllocate()`, `DSSetDimensions()`
@*/
PetscErrorCode DSGetLeadingDimension(DS ds,PetscInt *ld)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(ld,2);
  *ld = ds->ld;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetState - Change the state of the `DS` object.

   Logically Collective

   Input Parameters:
+  ds    - the direct solver context
-  state - the new state

   Notes:
   The state indicates that the dense system is in an initial state (raw),
   in an intermediate state (such as tridiagonal, Hessenberg or
   Hessenberg-triangular), in a condensed state (such as diagonal, Schur or
   generalized Schur), or in a truncated state.

   The state is automatically changed in functions such as `DSSolve()` or `DSTruncate()`.
   This function is normally used to return to the raw state when the
   condensed structure is destroyed, or to indicate that `DSSolve()` must
   start with a problem that already has an intermediate form.

   Level: advanced

.seealso: [](sec:ds), `DSGetState()`, `DSSolve()`, `DSTruncate()`
@*/
PetscErrorCode DSSetState(DS ds,DSStateType state)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(ds,state,2);
  switch (state) {
    case DS_STATE_RAW:
    case DS_STATE_INTERMEDIATE:
    case DS_STATE_CONDENSED:
    case DS_STATE_TRUNCATED:
      if (ds->state!=state) PetscCall(PetscInfo(ds,"State has changed from %s to %s\n",DSStateTypes[ds->state],DSStateTypes[state]));
      ds->state = state;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"Wrong state");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetState - Returns the current state.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  state - current state

   Level: advanced

.seealso: [](sec:ds), `DSSetState()`
@*/
PetscErrorCode DSGetState(DS ds,DSStateType *state)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(state,2);
  *state = ds->state;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetDimensions - Resize the matrices in the `DS` object.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
.  n  - the new size
.  l  - number of locked (inactive) leading columns
-  k  - intermediate dimension (e.g., position of arrow)

   Notes:
   The internal arrays are not reallocated.

   Some `DS` types have additional dimensions, e.g., the number of columns
   in `DSSVD`. For these, you should call a specific interface function.

   Use `PETSC_CURRENT` to leave any of the values unchanged. Use `PETSC_DETERMINE`
   to set `n` to the leading dimension, `l` to the minimum value (0), and `k` to `n/2`.

   Level: intermediate

.seealso: [](sec:ds), `DSGetDimensions()`, `DSAllocate()`, `DSSVDSetDimensions()`
@*/
PetscErrorCode DSSetDimensions(DS ds,PetscInt n,PetscInt l,PetscInt k)
{
  PetscInt       on,ol,ok;
  PetscBool      issvd;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  PetscValidLogicalCollectiveInt(ds,n,2);
  PetscValidLogicalCollectiveInt(ds,l,3);
  PetscValidLogicalCollectiveInt(ds,k,4);
  on = ds->n; ol = ds->l; ok = ds->k;
  if (n == PETSC_DETERMINE) {
    ds->n = ds->extrarow? ds->ld-1: ds->ld;
  } else if (n != PETSC_CURRENT) {
    PetscCheck(n>=0 && n<=ds->ld,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of n. Must be between 0 and ld");
    PetscCall(PetscObjectTypeCompareAny((PetscObject)ds,&issvd,DSSVD,DSGSVD,""));  /* SVD and GSVD have extra column instead of extra row */
    PetscCheck(!ds->extrarow || n<ds->ld || issvd,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"A value of n equal to ld leaves no room for extra row");
    ds->n = n;
  }
  ds->t = ds->n;   /* truncated length equal to the new dimension */
  if (l == PETSC_DETERMINE) {
    ds->l = 0;
  } else if (l != PETSC_CURRENT) {
    PetscCheck(l>=0 && l<=ds->n,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of l. Must be between 0 and n");
    ds->l = l;
  }
  if (k == PETSC_DETERMINE) {
    ds->k = ds->n/2;
  } else if (k != PETSC_CURRENT) {
    PetscCheck(k>=0 || k<=ds->n,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of k. Must be between 0 and n");
    ds->k = k;
  }
  if (on!=ds->n || ol!=ds->l || ok!=ds->k) PetscCall(PetscInfo(ds,"New dimensions are: n=%" PetscInt_FMT ", l=%" PetscInt_FMT ", k=%" PetscInt_FMT "\n",ds->n,ds->l,ds->k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetDimensions - Returns the current dimensions.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
+  n  - the current size
.  l  - number of locked (inactive) leading columns
.  k  - intermediate dimension (e.g., position of arrow)
-  t  - truncated length

   Notes:
   The `t` parameter makes sense only if `DSTruncate()` has been called.
   Otherwise its value equals `n`.

   Some `DS` types have additional dimensions, e.g., the number of columns
   in `DSSVD`. For these, you should call a specific interface function.

   Level: intermediate

.seealso: [](sec:ds), `DSSetDimensions()`, `DSTruncate()`, `DSSVDGetDimensions()`
@*/
PetscErrorCode DSGetDimensions(DS ds,PetscInt *n,PetscInt *l,PetscInt *k,PetscInt *t)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  if (n) *n = ds->n;
  if (l) *l = ds->l;
  if (k) *k = ds->k;
  if (t) *t = ds->t;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSTruncate - Truncates the system represented in the `DS` object.

   Logically Collective

   Input Parameters:
+  ds   - the direct solver context
.  n    - the new size
-  trim - a flag to indicate if the factorization must be trimmed

   Note:
   The new size is set to `n`. Note that in some cases the new size could
   be `n+1` or `n-1` to avoid breaking a 2x2 diagonal block (e.g., in real
   Schur form). In cases where the extra row is meaningful, the first
   `n` elements are kept as the extra row for the new system.

   If the flag `trim` is turned on, it resets the locked and intermediate
   dimensions to zero, see `DSSetDimensions()`, and sets the state to
   `DS_STATE_RAW`. It also cleans the extra row if being used.

   The typical usage of `trim=PETSC_TRUE` is to truncate the Schur decomposition
   at the end of a Krylov iteration. In this case, the state must be
   changed to `DS_STATE_RAW` so that `DSVectors()` computes eigenvectors from scratch.

   Level: advanced

.seealso: [](sec:ds), `DSSetDimensions()`, `DSSetExtraRow()`, `DSStateType`
@*/
PetscErrorCode DSTruncate(DS ds,PetscInt n,PetscBool trim)
{
  DSStateType    old;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckAlloc(ds,1);
  PetscValidLogicalCollectiveInt(ds,n,2);
  PetscValidLogicalCollectiveBool(ds,trim,3);
  PetscCheck(n>=ds->l && n<=ds->n,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of n (%" PetscInt_FMT "). Must be between l (%" PetscInt_FMT ") and n (%" PetscInt_FMT ")",n,ds->l,ds->n);
  PetscCall(PetscLogEventBegin(DS_Other,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,truncate,n,trim);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Other,ds,0,0,0));
  PetscCall(PetscInfo(ds,"Decomposition %s to size n=%" PetscInt_FMT "\n",trim?"trimmed":"truncated",ds->n));
  old = ds->state;
  ds->state = trim? DS_STATE_RAW: DS_STATE_TRUNCATED;
  if (old!=ds->state) PetscCall(PetscInfo(ds,"State has changed from %s to %s\n",DSStateTypes[old],DSStateTypes[ds->state]));
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSMatGetSize - Returns the numbers of rows and columns of one of the `DS` matrices.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
-  t  - the requested matrix

   Output Parameters:
+  n  - the number of rows
-  m  - the number of columns

   Note:
   This is equivalent to `MatGetSize()` on a matrix obtained with `DSGetMat()`.

   Level: developer

.seealso: [](sec:ds), `DSSetDimensions()`, `DSGetMat()`
@*/
PetscErrorCode DSMatGetSize(DS ds,DSMatType t,PetscInt *m,PetscInt *n)
{
  PetscInt       rows,cols;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckValidMat(ds,t,2);
  if (ds->ops->matgetsize) PetscUseTypeMethod(ds,matgetsize,t,&rows,&cols);
  else {
    if (ds->state==DS_STATE_TRUNCATED && t>=DS_MAT_Q) rows = ds->t;
    else rows = (t==DS_MAT_A && ds->extrarow)? ds->n+1: ds->n;
    if (t==DS_MAT_T) cols = PetscDefined(USE_COMPLEX)? 2: 3;
    else if (t==DS_MAT_D) cols = 1;
    else cols = ds->n;
  }
  if (m) *m = rows;
  if (n) *n = cols;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSMatIsHermitian - Checks if one of the `DS` matrices is known to be Hermitian.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
-  t  - the requested matrix

   Output Parameter:
.  flg - the Hermitian flag

   Note:
   Does not check the matrix values directly. The flag is set according to the
   problem structure. For instance, in `DSHEP` matrix `A` is Hermitian.

   Level: developer

.seealso: [](sec:ds), `DSGetMat()`
@*/
PetscErrorCode DSMatIsHermitian(DS ds,DSMatType t,PetscBool *flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckValidMat(ds,t,2);
  PetscAssertPointer(flg,3);
  *flg = PETSC_FALSE;
  PetscTryTypeMethod(ds,hermitian,t,flg);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DSGetTruncateSize_Default(DS ds,PetscInt l,PetscInt n,PetscInt *k)
{
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar val;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(MatGetValue(ds->omat[DS_MAT_A],l+(*k),l+(*k)-1,&val));
  if (val != 0.0) {
    if (l+(*k)<n-1) (*k)++;
    else (*k)--;
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetTruncateSize - Gets the correct size to be used in `DSTruncate()`
   to avoid breaking a 2x2 block.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
.  l  - the size of the locked part (set to 0 to use `ds->l`)
-  n  - the total matrix size (set to 0 to use `ds->n`)

   Output Parameter:
.  k  - the wanted truncation size (possibly modified)

   Notes:
   This should be called before `DSTruncate()` to make sure that the truncation
   does not break a 2x2 block corresponding to a complex conjugate eigenvalue.

   The total size is `n` (either user-provided or `ds->n` if 0 is passed). The
   size where the truncation is intended is equal to `l+k` (where `l` can be
   equal to the locked size `ds->l` if set to 0). Then if there is a 2x2 block
   at the `l+k` limit, the value of `k` is increased (or decreased) by 1.

   Level: advanced

.seealso: [](sec:ds), `DSTruncate()`, `DSSetDimensions()`
@*/
PetscErrorCode DSGetTruncateSize(DS ds,PetscInt l,PetscInt n,PetscInt *k)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  PetscValidLogicalCollectiveInt(ds,l,2);
  PetscValidLogicalCollectiveInt(ds,n,3);
  PetscAssertPointer(k,4);
  PetscUseTypeMethod(ds,gettruncatesize,l?l:ds->l,n?n:ds->n,k);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetMat - Returns a sequential dense `Mat` object containing the requested
   matrix.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
-  m  - the requested matrix

   Output Parameter:
.  A  - Mat object

   Notes:
   The returned `Mat` has sizes equal to the current `DS` dimensions (see `DSSetDimensions()`),
   and contains the values that would be obtained with `DSGetArray()`
   (not `DSGetArrayReal()`). If the `DS` was truncated, then the number of rows
   is equal to the dimension prior to truncation, see `DSTruncate()`.
   The communicator is always `PETSC_COMM_SELF`.

   It is implemented with `MatDenseGetSubMatrix()`, and when no longer needed
   the user must call `DSRestoreMat()` which will invoke `MatDenseRestoreSubMatrix()`.

   For matrices `DS_MAT_T` and `DS_MAT_D`, this function will return a `Mat` object
   that cannot be used directly for computations, since it uses compact storage
   (three and one diagonals for $T$ and $D$, respectively). In complex scalars, the
   internal array stores real values, so it is sufficient with two columns for $T$.

   Level: advanced

.seealso: [](sec:ds), `DSRestoreMat()`, `DSSetDimensions()`, `DSGetArray()`, `DSGetArrayReal()`, `DSTruncate()`, `DSGetMatAndColumn()`
@*/
PetscErrorCode DSGetMat(DS ds,DSMatType m,Mat *A)
{
  PetscInt  rows=0,cols=0;
  PetscBool flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  DSCheckValidMat(ds,m,2);
  PetscAssertPointer(A,3);

  PetscCall(DSMatGetSize(ds,m,&rows,&cols));
  PetscCheck(rows && cols,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"Must call DSSetDimensions() first");
  PetscCall(MatDenseGetSubMatrix(ds->omat[m],0,rows,0,cols,A));

  /* set Hermitian flag */
  PetscCall(DSMatIsHermitian(ds,m,&flg));
  PetscCall(MatSetOption(*A,MAT_SYMMETRIC,flg));
  PetscCall(MatSetOption(*A,MAT_HERMITIAN,flg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSRestoreMat - Restores the matrix after `DSGetMat()` was called.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
.  m  - the requested matrix
-  A  - the fetched `Mat` object

   Note:
   A call to this function must match a previous call of `DSGetMat()`.

   Level: advanced

.seealso: [](sec:ds), `DSGetMat()`, `DSRestoreArray()`, `DSRestoreArrayReal()`
@*/
PetscErrorCode DSRestoreMat(DS ds,DSMatType m,Mat *A)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  DSCheckValidMat(ds,m,2);
  PetscAssertPointer(A,3);

  PetscCall(MatDenseRestoreSubMatrix(ds->omat[m],A));
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetMatAndColumn - Returns a sequential dense `Mat` object containing the requested
   matrix and one of its columns as a `Vec`.

   Not Collective

   Input Parameters:
+  ds  - the direct solver context
.  m   - the requested matrix
-  col - the index of the requested column

   Output Parameters:
+  A   - `Mat` object
-  v   - `Vec` object (the column)

   Notes:
   This calls `DSGetMat()` and then it extracts the selected column.
   The user must call `DSRestoreMatAndColumn()` to recover the original state.
   For matrices `DS_MAT_T` and `DS_MAT_D`, in complex scalars this function implies
   copying from real values stored internally to scalar values in the Vec.

   Level: advanced

.seealso: [](sec:ds), `DSRestoreMatAndColumn()`, `DSGetMat()`
@*/
PetscErrorCode DSGetMatAndColumn(DS ds,DSMatType m,PetscInt col,Mat *A,Vec *v)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  DSCheckValidMat(ds,m,2);
  PetscAssertPointer(A,4);
  PetscAssertPointer(v,5);

  PetscCall(DSGetMat(ds,m,A));
  if (PetscDefined(USE_COMPLEX) && (m==DS_MAT_T || m==DS_MAT_D)) {
    const PetscScalar *as;
    PetscScalar       *vs;
    PetscReal         *ar;
    PetscInt          i,n,lda;
    PetscCall(MatCreateVecs(*A,NULL,v));
    PetscCall(VecGetSize(*v,&n));
    PetscCall(MatDenseGetLDA(*A,&lda));
    PetscCall(MatDenseGetArrayRead(*A,&as));
    PetscCall(VecGetArrayWrite(*v,&vs));
    ar = (PetscReal*)as;
    for (i=0;i<n;i++) vs[i] = ar[i+col*lda];
    PetscCall(VecRestoreArrayWrite(*v,&vs));
    PetscCall(MatDenseRestoreArrayRead(*A,&as));
  } else PetscCall(MatDenseGetColumnVec(*A,col,v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSRestoreMatAndColumn - Restores the matrix and vector after `DSGetMatAndColumn()`
   was called.

   Not Collective

   Input Parameters:
+  ds  - the direct solver context
.  m   - the requested matrix
.  col - the requested column
.  A   - the fetched `Mat` object
-  v   - the fetched `Vec` object

   Note:
   A call to this function must match a previous call of `DSGetMatAndColumn()`.

   Level: advanced

.seealso: [](sec:ds), `DSGetMatAndColumn()`, `DSRestoreMat()`
@*/
PetscErrorCode DSRestoreMatAndColumn(DS ds,DSMatType m,PetscInt col,Mat *A,Vec *v)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  DSCheckValidMat(ds,m,2);
  PetscAssertPointer(A,4);
  PetscAssertPointer(v,5);

  if (PetscDefined(USE_COMPLEX) && (m==DS_MAT_T || m==DS_MAT_D)) {
    const PetscScalar *vs;
    PetscScalar       *as;
    PetscReal         *ar;
    PetscInt          i,n,lda;
    PetscCall(VecGetSize(*v,&n));
    PetscCall(MatDenseGetLDA(*A,&lda));
    PetscCall(MatDenseGetArray(*A,&as));
    PetscCall(VecGetArrayRead(*v,&vs));
    ar = (PetscReal*)as;
    for (i=0;i<n;i++) ar[i+col*lda] = PetscRealPart(vs[i]);
    PetscCall(VecRestoreArrayRead(*v,&vs));
    PetscCall(MatDenseRestoreArray(*A,&as));
    PetscCall(VecDestroy(v));
  } else PetscCall(MatDenseRestoreColumnVec(*A,col,v));
  PetscCall(DSRestoreMat(ds,m,A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSGetArray - Returns a pointer to the internal array of one of the
   matrices. You MUST call `DSRestoreArray()` when you no longer
   need to access the array.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
-  m  - the requested matrix

   Output Parameter:
.  a  - pointer to the values

   Note:
   To get read-only access, use `DSGetMat()` followed by `MatDenseGetArrayRead()`.

   Level: advanced

.seealso: [](sec:ds), `DSRestoreArray()`, `DSGetArrayReal()`
@*/
PetscErrorCode DSGetArray(DS ds,DSMatType m,PetscScalar *a[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  DSCheckValidMat(ds,m,2);
  PetscAssertPointer(a,3);
  PetscCall(MatDenseGetArray(ds->omat[m],a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSRestoreArray - Restores the matrix after `DSGetArray()` was called.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
.  m  - the requested matrix
-  a  - pointer to the values

   Level: advanced

.seealso: [](sec:ds), `DSGetArray()`, `DSGetArrayReal()`
@*/
PetscErrorCode DSRestoreArray(DS ds,DSMatType m,PetscScalar *a[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  DSCheckValidMat(ds,m,2);
  PetscAssertPointer(a,3);
  PetscCall(MatDenseRestoreArray(ds->omat[m],a));
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSGetArrayReal - Returns a real pointer to the internal array of $T$ or $D$.
   You MUST call `DSRestoreArrayReal()` when you no longer need to access the array.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
-  m  - the requested matrix

   Output Parameter:
.  a  - pointer to the values

   Note:
   This function can be used only for `DS_MAT_T` and `DS_MAT_D`. These matrices always
   store real values, even in complex scalars, that is why the returned pointer is
   `PetscReal`.

   Level: advanced

.seealso: [](sec:ds), `DSRestoreArrayReal()`, `DSGetArray()`
@*/
PetscErrorCode DSGetArrayReal(DS ds,DSMatType m,PetscReal *a[])
{
#if defined(PETSC_USE_COMPLEX)
  PetscScalar *as;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  DSCheckValidMatReal(ds,m,2);
  PetscAssertPointer(a,3);
#if defined(PETSC_USE_COMPLEX)
  PetscCall(MatDenseGetArray(ds->omat[m],&as));
  *a = (PetscReal*)as;
#else
  PetscCall(MatDenseGetArray(ds->omat[m],a));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSRestoreArrayReal - Restores the matrix after `DSGetArrayReal()` was called.

   Not Collective

   Input Parameters:
+  ds - the direct solver context
.  m  - the requested matrix
-  a  - pointer to the values

   Level: advanced

.seealso: [](sec:ds), `DSGetArrayReal()`, `DSGetArray()`
@*/
PetscErrorCode DSRestoreArrayReal(DS ds,DSMatType m,PetscReal *a[])
{
#if defined(PETSC_USE_COMPLEX)
  PetscScalar *as;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  DSCheckValidMatReal(ds,m,2);
  PetscAssertPointer(a,3);
#if defined(PETSC_USE_COMPLEX)
  PetscCall(MatDenseRestoreArray(ds->omat[m],&as));
  *a = NULL;
#else
  PetscCall(MatDenseRestoreArray(ds->omat[m],a));
#endif
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSolve - Solve the problem.

   Logically Collective

   Input Parameters:
+  ds   - the direct solver context
.  eigr - array to store the computed eigenvalues (real part)
-  eigi - array to store the computed eigenvalues (imaginary part)

   Notes:
   This call brings the dense system to condensed form. No ordering
   of the eigenvalues is enforced (for this, call `DSSort()` afterwards).

   In some `DS` types, the arguments will hold singular values,
   instead of eigenvalues. Singular values are always real.

   Level: intermediate

.seealso: [](sec:ds), `DSSort()`, `DSStateType`
@*/
PetscErrorCode DSSolve(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckAlloc(ds,1);
  PetscAssertPointer(eigr,2);
  if (ds->state>=DS_STATE_CONDENSED) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCheck(ds->ops->solve[ds->method],PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The specified method number does not exist for this DS");
  PetscCall(PetscInfo(ds,"Starting solve with problem sizes: n=%" PetscInt_FMT ", l=%" PetscInt_FMT ", k=%" PetscInt_FMT "\n",ds->n,ds->l,ds->k));
  PetscCall(PetscLogEventBegin(DS_Solve,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,solve[ds->method],eigr,eigi);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Solve,ds,0,0,0));
  PetscCall(PetscInfo(ds,"State has changed from %s to CONDENSED\n",DSStateTypes[ds->state]));
  ds->state = DS_STATE_CONDENSED;
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSort - Sort the result of `DSSolve()` according to a given sorting
   criterion.

   Logically Collective

   Input Parameters:
+  ds   - the direct solver context
.  rr   - (optional) array containing auxiliary values (real part)
-  ri   - (optional) array containing auxiliary values (imaginary part)

   Input/Output Parameters:
+  eigr - array containing the computed eigenvalues (real part)
.  eigi - array containing the computed eigenvalues (imaginary part)
-  k    - (optional) number of elements in the leading group

   Notes:
   This routine sorts the arrays provided in `eigr` and `eigi`, and also
   sorts the dense system stored inside `ds` (assumed to be in condensed form).
   The sorting criterion is specified with `DSSetSlepcSC()`.

   If arrays `rr` and `ri` are provided, then a (partial) reordering based on these
   values rather than on the eigenvalues is performed. In symmetric problems
   a total order is obtained (parameter `k` is ignored), but otherwise the result
   is sorted only partially. In this latter case, it is only guaranteed that
   all the first `k` elements satisfy the comparison with any of the last `n-k`
   elements. The output value of parameter `k` is the final number of elements in
   the first set.

   Level: intermediate

.seealso: [](sec:ds), `DSSolve()`, `DSSetSlepcSC()`, `DSSortWithPermutation()`
@*/
PetscErrorCode DSSort(DS ds,PetscScalar eigr[],PetscScalar eigi[],PetscScalar rr[],PetscScalar ri[],PetscInt *k)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckSolved(ds,1);
  PetscAssertPointer(eigr,2);
  if (rr) PetscAssertPointer(rr,4);
  PetscCheck(ds->state<DS_STATE_TRUNCATED,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"Cannot sort a truncated DS");
  PetscCheck(ds->sc,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"Must provide a sorting criterion first");
  PetscCheck(!k || rr,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"Argument k can only be used together with rr");

  for (i=0;i<ds->n;i++) ds->perm[i] = i;   /* initialize to trivial permutation */
  PetscCall(PetscLogEventBegin(DS_Other,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,sort,eigr,eigi,rr,ri,k);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Other,ds,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscCall(PetscInfo(ds,"Finished sorting\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSortWithPermutation - Reorders the result of `DSSolve()` according to a given
   permutation.

   Logically Collective

   Input Parameters:
+  ds   - the direct solver context
-  perm - permutation that indicates the new ordering

   Input/Output Parameters:
+  eigr - array with the reordered eigenvalues (real part)
-  eigi - array with the reordered eigenvalues (imaginary part)

   Notes:
   This routine reorders the arrays provided in `eigr` and `eigi`, and also the dense
   system stored inside `ds` (assumed to be in condensed form). There is no sorting
   criterion, as opposed to `DSSort()`. Instead, the new ordering is given in argument
   `perm`.

   Level: advanced

.seealso: [](sec:ds), `DSSolve()`, `DSSort()`
@*/
PetscErrorCode DSSortWithPermutation(DS ds,PetscInt perm[],PetscScalar eigr[],PetscScalar eigi[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckSolved(ds,1);
  PetscAssertPointer(perm,2);
  PetscAssertPointer(eigr,3);
  PetscCheck(ds->state<DS_STATE_TRUNCATED,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"Cannot sort a truncated DS");

  PetscCall(PetscLogEventBegin(DS_Other,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,sortperm,perm,eigr,eigi);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Other,ds,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscCall(PetscInfo(ds,"Finished sorting\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSynchronize - Make sure that all processes have the same data, performing
   communication if necessary.

   Collective

   Input Parameter:
.  ds   - the direct solver context

   Input/Output Parameters:
+  eigr - (optional) array with the computed eigenvalues (real part)
-  eigi - (optional) array with the computed eigenvalues (imaginary part)

   Notes:
   When the `DS` has been created with a communicator with more than one process,
   the internal data, especially the computed matrices, may diverge in the
   different processes. This happens when using multithreaded BLAS and may
   cause numerical issues in some ill-conditioned problems. This function
   performs the necessary communication among the processes so that the
   internal data is exactly equal in all of them.

   Depending on the parallel mode as set with `DSSetParallel()`, this function
   will either do nothing or synchronize the matrices computed by `DSSolve()`
   and `DSSort()`. The arguments `eigr` and `eigi` are typically those used in the
   calls to `DSSolve()` and `DSSort()`.

   Level: developer

.seealso: [](sec:ds), `DSSetParallel()`, `DSSolve()`, `DSSort()`
@*/
PetscErrorCode DSSynchronize(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckAlloc(ds,1);
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ds),&size));
  if (size>1 && ds->pmode==DS_PARALLEL_SYNCHRONIZED) {
    PetscCall(PetscLogEventBegin(DS_Synchronize,ds,0,0,0));
    PetscUseTypeMethod(ds,synchronize,eigr,eigi);
    PetscCall(PetscLogEventEnd(DS_Synchronize,ds,0,0,0));
    PetscCall(PetscObjectStateIncrease((PetscObject)ds));
    PetscCall(PetscInfo(ds,"Synchronization completed (%s)\n",DSParallelTypes[ds->pmode]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSVectors - Compute vectors associated to the dense system such
   as eigenvectors.

   Logically Collective

   Input Parameters:
+  ds  - the direct solver context
-  mat - the matrix, used to indicate which vectors are required

   Input/Output Parameter:
.  j   - (optional) index of vector to be computed

   Output Parameter:
.  rnorm - (optional) computed residual norm

   Notes:
   Allowed values for mat are `DS_MAT_X`, `DS_MAT_Y`, `DS_MAT_U` and `DS_MAT_V`, to
   compute right or left eigenvectors, or left or right singular vectors,
   respectively.

   If `NULL` is passed in argument `j` then all vectors are computed,
   otherwise `j` indicates which vector must be computed. In real non-symmetric
   problems, on exit the index `j` will be incremented when a complex conjugate
   pair is found.

   This function can be invoked after the dense problem has been solved,
   to get the residual norm estimate of the associated Ritz pair. In that
   case, the relevant information is returned in `rnorm`.

   For computing eigenvectors, LAPACK's `_trevc` is used so the matrix must
   be in (quasi-)triangular form, or call `DSSolve()` first.

   Level: intermediate

.seealso: [](sec:ds), `DSSolve()`
@*/
PetscErrorCode DSVectors(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckAlloc(ds,1);
  PetscValidLogicalCollectiveEnum(ds,mat,2);
  PetscCheck(mat<DS_NUM_MAT,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"Invalid matrix");
  PetscCheck(!rnorm || j,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"Must give a value of j");
  if (!ds->omat[mat]) PetscCall(DSAllocateMat_Private(ds,mat));
  if (!j) PetscCall(PetscInfo(ds,"Computing all vectors on %s\n",DSMatName[mat]));
  PetscCall(PetscLogEventBegin(DS_Vectors,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,vectors,mat,j,rnorm);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Vectors,ds,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSUpdateExtraRow - Performs all necessary operations so that the extra
   row gets up-to-date after a call to `DSSolve()`.

   Logically Collective

   Input Parameter:
.  ds - the direct solver context

   Level: advanced

.seealso: [](sec:ds), `DSSolve()`, `DSSetExtraRow()`
@*/
PetscErrorCode DSUpdateExtraRow(DS ds)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckAlloc(ds,1);
  PetscCheck(ds->extrarow,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"Should have called DSSetExtraRow");
  PetscCall(PetscInfo(ds,"Updating extra row\n"));
  PetscCall(PetscLogEventBegin(DS_Other,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,update);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Other,ds,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSCond - Compute the condition number.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  cond - the computed condition number

   Notes:
   In standard eigenvalue problems, returns the $\infty$-norm condition number of the first
   matrix, computed as $\kappa_\infty(A) = \|A\|_\infty\|A^{-1}\|_\infty$.

   In GSVD problems, returns the maximum of $\kappa_2(A)$ and $\kappa_2(B)$, where $\kappa_2(\cdot)$ is
   computed as the ratio of the largest and smallest singular values.

   Does not take into account the extra row.

   Level: advanced

.seealso: [](sec:ds), `DSSolve()`, `DSSetExtraRow()`
@*/
PetscErrorCode DSCond(DS ds,PetscReal *cond)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckAlloc(ds,1);
  PetscAssertPointer(cond,2);
  PetscCall(PetscLogEventBegin(DS_Other,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,cond,cond);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Other,ds,0,0,0));
  PetscCall(PetscInfo(ds,"Computed condition number = %g\n",(double)*cond));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSTranslateHarmonic - Computes a translation of the dense system.

   Logically Collective

   Input Parameters:
+  ds      - the direct solver context
.  tau     - the translation amount
.  beta    - last component of vector $b$
-  recover - boolean flag to indicate whether to recover or not

   Output Parameters:
+  g       - the computed vector (optional)
-  gamma   - scale factor (optional)

   Notes:
   This function is intended for use in the context of Krylov methods only.
   It computes a translation of a Krylov decomposition in order to extract
   eigenpair approximations by harmonic Rayleigh-Ritz.
   The matrix is updated as $A + gb^*$ where $g = (A-\tau I)^{-*}b$ and
   vector $b$ is assumed to be $\beta e_n^*$.

   The $\gamma$ factor is defined as $\sqrt{1+g^*g}$ and can be interpreted as
   the factor by which the residual of the Krylov decomposition is scaled.

   If the `recover` flag is activated, the computed translation undoes the
   translation done previously. In that case, parameter `tau` is ignored.

   Level: developer

.seealso: [](sec:ds), `DSTranslateRKS()`
@*/
PetscErrorCode DSTranslateHarmonic(DS ds,PetscScalar tau,PetscReal beta,PetscBool recover,PetscScalar g[],PetscReal *gamma)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckAlloc(ds,1);
  if (recover) PetscCall(PetscInfo(ds,"Undoing the translation\n"));
  else PetscCall(PetscInfo(ds,"Computing the translation\n"));
  PetscCall(PetscLogEventBegin(DS_Other,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,transharm,tau,beta,recover,g,gamma);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Other,ds,0,0,0));
  ds->state = DS_STATE_RAW;
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSTranslateRKS - Computes a modification of the dense system corresponding
   to an update of the shift in a rational Krylov method.

   Logically Collective

   Input Parameters:
+  ds    - the direct solver context
-  alpha - the translation amount

   Notes:
   This function is intended for use in the context of Krylov methods only.
   It takes the leading $(k+1,k)$ submatrix of $A$, containing the truncated
   Rayleigh quotient of a Krylov-Schur relation computed from a shift
   $\sigma_1$ and transforms it to obtain a Krylov relation as if computed
   from a different shift $\sigma_2$. The new matrix is computed as
   $\alpha^{-1}(I-QR^{-1})$, where $QR=I-\alpha A$ and
   $\alpha = \sigma_1-\sigma_2$.

   Matrix $Q$ is placed in `DS_MAT_Q` so that it can be used to update the
   Krylov basis.

   Level: developer

.seealso: [](sec:ds), `DSTranslateHarmonic()`
@*/
PetscErrorCode DSTranslateRKS(DS ds,PetscScalar alpha)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidType(ds,1);
  DSCheckAlloc(ds,1);
  PetscCall(PetscInfo(ds,"Translating with alpha=%g\n",(double)PetscRealPart(alpha)));
  PetscCall(DSSetCompact(ds,PETSC_FALSE));
  PetscCall(PetscLogEventBegin(DS_Other,ds,0,0,0));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscUseTypeMethod(ds,transrks,alpha);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogEventEnd(DS_Other,ds,0,0,0));
  ds->state = DS_STATE_RAW;
  PetscCall(PetscObjectStateIncrease((PetscObject)ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}
