/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV operations, except those involving global communication
*/

#include <slepc/private/bvimpl.h>      /*I "slepcbv.h" I*/
#include <slepcds.h>

/*@
   BVMult - Computes $Y = \beta Y + \alpha X Q$.

   Logically Collective

   Input Parameters:
+  Y     - first basis vectors context (modified on output)
.  alpha - first scalar
.  beta  - second scalar
.  X     - second basis vectors context
-  Q     - (optional) sequential dense matrix

   Notes:
   `X` and `Y` must be different objects. The case `X=Y` can be addressed with
   `BVMultInPlace()`.

   If matrix `Q` is `NULL`, then an AXPY operation $Y = \beta Y + \alpha X$ is
   done (i.e., results as if $Q = I$). If provided,
   the matrix `Q` must be a sequential dense `Mat`, with all entries equal on
   all processes (otherwise each process will compute a different update).
   The dimensions of `Q` must be at least `m,n` where `m` is the number of active
   columns of `X` and `n` is the number of active columns of `Y`.

   The leading columns of `Y` are not modified. Also, if `X` has leading
   columns specified, then these columns do not participate in the computation.
   Hence, only rows (resp. columns) of `Q` starting from `lx` (resp. `ly`) are used,
   where `lx` (resp. `ly`) is the number of leading columns of `X` (resp. `Y`).

   Level: intermediate

.seealso: [](sec:bv), `BVMultVec()`, `BVMultColumn()`, `BVMultInPlace()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVMult(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Y,BV_CLASSID,1);
  PetscValidLogicalCollectiveScalar(Y,alpha,2);
  PetscValidLogicalCollectiveScalar(Y,beta,3);
  PetscValidHeaderSpecific(X,BV_CLASSID,4);
  if (Q) PetscValidHeaderSpecific(Q,MAT_CLASSID,5);
  PetscValidType(Y,1);
  BVCheckSizes(Y,1);
  BVCheckOp(Y,1,mult);
  PetscValidType(X,4);
  BVCheckSizes(X,4);
  if (Q) PetscValidType(Q,5);
  PetscCheckSameTypeAndComm(Y,1,X,4);
  PetscCheck(X!=Y,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_WRONG,"X and Y arguments must be different");
  if (Q) {
    SlepcMatCheckSeq(Q);
    PetscCall(MatGetSize(Q,&m,&n));
    PetscCheck(m>=X->k,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_SIZ,"Mat argument has %" PetscInt_FMT " rows, should have at least %" PetscInt_FMT,m,X->k);
    PetscCheck(n>=Y->k,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_SIZ,"Mat argument has %" PetscInt_FMT " columns, should have at least %" PetscInt_FMT,n,Y->k);
  }
  PetscCheck(X->n==Y->n,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %" PetscInt_FMT ", Y %" PetscInt_FMT,X->n,Y->n);

  PetscCall(PetscLogEventBegin(BV_Mult,X,Y,0,0));
  PetscUseTypeMethod(Y,mult,alpha,beta,X,Q);
  PetscCall(PetscLogEventEnd(BV_Mult,X,Y,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)Y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMultVec - Computes $y = \beta y + \alpha X q$.

   Logically Collective

   Input Parameters:
+  X     - a basis vectors object
.  alpha - first scalar
.  beta  - second scalar
.  y     - a vector (modified on output)
-  q     - an array of scalars

   Notes:
   This operation is the analogue of `BVMult()` but with a `BV` and a `Vec`,
   instead of two `BV`s. Note that arguments are listed in different order
   with respect to `BVMult()`.

   If `X` has leading columns specified, then these columns do not participate
   in the computation.

   The length of array `q` must be equal to the number of active columns of `X`
   minus the number of leading columns, i.e., the first entry of `q` multiplies
   the first non-leading column.

   Level: intermediate

.seealso: [](sec:bv), `BVMult()`, `BVMultColumn()`, `BVMultInPlace()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVMultVec(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar q[])
{
  PetscInt       n,N;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveScalar(X,alpha,2);
  PetscValidLogicalCollectiveScalar(X,beta,3);
  PetscValidHeaderSpecific(y,VEC_CLASSID,4);
  PetscAssertPointer(q,5);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  BVCheckOp(X,1,multvec);
  PetscValidType(y,4);
  PetscCheckSameComm(X,1,y,4);

  PetscCall(VecGetSize(y,&N));
  PetscCall(VecGetLocalSize(y,&n));
  PetscCheck(N==X->N && n==X->n,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Vec sizes (global %" PetscInt_FMT ", local %" PetscInt_FMT ") do not match BV sizes (global %" PetscInt_FMT ", local %" PetscInt_FMT ")",N,n,X->N,X->n);

  PetscCall(PetscLogEventBegin(BV_MultVec,X,y,0,0));
  PetscUseTypeMethod(X,multvec,alpha,beta,y,q);
  PetscCall(PetscLogEventEnd(BV_MultVec,X,y,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMultColumn - Computes $y = \beta y + \alpha X q$, where $y$ is the j-th column
   of $X$.

   Logically Collective

   Input Parameters:
+  X     - a basis vectors object
.  alpha - first scalar
.  beta  - second scalar
.  j     - the column index
-  q     - an array of scalars

   Notes:
   This operation is equivalent to `BVMultVec()` but it uses column `j` of `X`
   rather than taking a `Vec` as an argument. The number of active columns of
   `X` is set to `j` before the computation, and restored afterwards.
   If `X` has leading columns specified, then these columns do not participate
   in the computation. Therefore, the length of array `q` must be equal to `j`
   minus the number of leading columns.

   Developer Note:
   If `q` is `NULL`, then the coefficients are taken from position `nc+l` of the
   internal buffer vector, see `BVGetBufferVec()`.

   Level: advanced

.seealso: [](sec:bv), `BVMult()`, `BVMultVec()`, `BVMultInPlace()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVMultColumn(BV X,PetscScalar alpha,PetscScalar beta,PetscInt j,PetscScalar q[])
{
  PetscInt       ksave;
  Vec            y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveScalar(X,alpha,2);
  PetscValidLogicalCollectiveScalar(X,beta,3);
  PetscValidLogicalCollectiveInt(X,j,4);
  PetscValidType(X,1);
  BVCheckSizes(X,1);

  PetscCheck(j>=0,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j<X->m,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j,X->m);

  PetscCall(PetscLogEventBegin(BV_MultVec,X,0,0,0));
  ksave = X->k;
  X->k = j;
  if (!q && !X->buffer) PetscCall(BVGetBufferVec(X,&X->buffer));
  PetscCall(BVGetColumn(X,j,&y));
  PetscUseTypeMethod(X,multvec,alpha,beta,y,q);
  PetscCall(BVRestoreColumn(X,j,&y));
  X->k = ksave;
  PetscCall(PetscLogEventEnd(BV_MultVec,X,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMultInPlace - Update a set of vectors as $V_{:,s:e-1} = V Q_{:,s:e-1}$.

   Logically Collective

   Input Parameters:
+  Q - a sequential dense matrix
.  s - first column of `V` to be overwritten
-  e - first column of `V` not to be overwritten

   Input/Output Parameter:
.  V - basis vectors

   Notes:
   The matrix `Q` must be a sequential dense `Mat`, with all entries equal on
   all processes (otherwise each process will compute a different update).

   This function computes $V_{:,s:e-1} = V Q_{:,s:e-1}$, that is, given a set of
   vectors $V$, columns from $s$ to $e-1$ are overwritten with columns from $s$ to
   $e-1$ of the matrix-matrix product $V Q$. Only columns $s$ to $e-1$ of `Q` are
   referenced.

   Level: intermediate

.seealso: [](sec:bv), `BVMult()`, `BVMultVec()`, `BVMultInPlaceHermitianTranspose()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVMultInPlace(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidHeaderSpecific(Q,MAT_CLASSID,2);
  PetscValidLogicalCollectiveInt(V,s,3);
  PetscValidLogicalCollectiveInt(V,e,4);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidType(Q,2);
  SlepcMatCheckSeq(Q);

  PetscCheck(s>=V->l && s<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument s has wrong value %" PetscInt_FMT ", should be between %" PetscInt_FMT " and %" PetscInt_FMT,s,V->l,V->m);
  PetscCheck(e>=V->l && e<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument e has wrong value %" PetscInt_FMT ", should be between %" PetscInt_FMT " and %" PetscInt_FMT,e,V->l,V->m);
  PetscCall(MatGetSize(Q,&m,&n));
  PetscCheck(m>=V->k,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument has %" PetscInt_FMT " rows, should have at least %" PetscInt_FMT,m,V->k);
  PetscCheck(e<=n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument only has %" PetscInt_FMT " columns, the requested value of e is larger: %" PetscInt_FMT,n,e);

  PetscCall(PetscLogEventBegin(BV_MultInPlace,V,Q,0,0));
  PetscUseTypeMethod(V,multinplace,Q,s,e);
  PetscCall(PetscLogEventEnd(BV_MultInPlace,V,Q,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMultInPlaceHermitianTranspose - Update a set of vectors as $V_{:,s:e-1} = V Q_{:,s:e-1}^*$.

   Logically Collective

   Input Parameters:
+  Q - a sequential dense matrix
.  s - first column of `V` to be overwritten
-  e - first column of `V` not to be overwritten

   Input/Output Parameter:
.  V - basis vectors

   Notes:
   This is a variant of `BVMultInPlace()` where the conjugate transpose
   of `Q` is used.

   Level: intermediate

.seealso: [](sec:bv), `BVMultInPlace()`
@*/
PetscErrorCode BVMultInPlaceHermitianTranspose(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidHeaderSpecific(Q,MAT_CLASSID,2);
  PetscValidLogicalCollectiveInt(V,s,3);
  PetscValidLogicalCollectiveInt(V,e,4);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidType(Q,2);
  SlepcMatCheckSeq(Q);

  PetscCheck(s>=V->l && s<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument s has wrong value %" PetscInt_FMT ", should be between %" PetscInt_FMT " and %" PetscInt_FMT,s,V->l,V->m);
  PetscCheck(e>=V->l && e<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument e has wrong value %" PetscInt_FMT ", should be between %" PetscInt_FMT " and %" PetscInt_FMT,e,V->l,V->m);
  PetscCall(MatGetSize(Q,&m,&n));
  PetscCheck(n>=V->k,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument has %" PetscInt_FMT " columns, should have at least %" PetscInt_FMT,n,V->k);
  PetscCheck(e<=m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument only has %" PetscInt_FMT " rows, the requested value of e is larger: %" PetscInt_FMT,m,e);

  PetscCall(PetscLogEventBegin(BV_MultInPlace,V,Q,0,0));
  PetscUseTypeMethod(V,multinplacetrans,Q,s,e);
  PetscCall(PetscLogEventEnd(BV_MultInPlace,V,Q,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVScale - Multiply the `BV` entries by a scalar value.

   Logically Collective

   Input Parameters:
+  bv    - basis vectors
-  alpha - scaling factor

   Note:
   All active columns (except the leading ones) are scaled.

   Level: intermediate

.seealso: [](sec:bv), `BVScaleColumn()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVScale(BV bv,PetscScalar alpha)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveScalar(bv,alpha,2);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  if (alpha == (PetscScalar)1.0) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscLogEventBegin(BV_Scale,bv,0,0,0));
  PetscUseTypeMethod(bv,scale,-1,alpha);
  PetscCall(PetscLogEventEnd(BV_Scale,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVScaleColumn - Scale one column of a `BV`.

   Logically Collective

   Input Parameters:
+  bv    - basis vectors
.  j     - index of column to be scaled
-  alpha - scaling factor

   Level: intermediate

.seealso: [](sec:bv), `BVScale()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVScaleColumn(BV bv,PetscInt j,PetscScalar alpha)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidLogicalCollectiveScalar(bv,alpha,3);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCheck(j>=0 && j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %" PetscInt_FMT ", the number of columns is %" PetscInt_FMT,j,bv->m);
  if (alpha == (PetscScalar)1.0) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscLogEventBegin(BV_Scale,bv,0,0,0));
  PetscUseTypeMethod(bv,scale,j,alpha);
  PetscCall(PetscLogEventEnd(BV_Scale,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode BVSetRandomColumn_Private(BV bv,PetscInt k)
{
  PetscInt       i,low,high;
  PetscScalar    *px,t;
  Vec            x;

  PetscFunctionBegin;
  PetscCall(BVGetColumn(bv,k,&x));
  if (bv->rrandom) {  /* generate the same vector irrespective of number of processes */
    PetscCall(VecGetOwnershipRange(x,&low,&high));
    PetscCall(VecGetArray(x,&px));
    for (i=0;i<bv->N;i++) {
      PetscCall(PetscRandomGetValue(bv->rand,&t));
      if (i>=low && i<high) px[i-low] = t;
    }
    PetscCall(VecRestoreArray(x,&px));
  } else PetscCall(VecSetRandom(x,bv->rand));
  PetscCall(BVRestoreColumn(bv,k,&x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode BVSetRandomNormalColumn_Private(BV bv,PetscInt k,Vec w1,Vec w2)
{
  PetscInt       i,low,high;
  PetscScalar    *px,s,t;
  Vec            x;

  PetscFunctionBegin;
  PetscCall(BVGetColumn(bv,k,&x));
  if (bv->rrandom) {  /* generate the same vector irrespective of number of processes */
    PetscCall(VecGetOwnershipRange(x,&low,&high));
    PetscCall(VecGetArray(x,&px));
    for (i=0;i<bv->N;i++) {
      PetscCall(PetscRandomGetValue(bv->rand,&s));
      PetscCall(PetscRandomGetValue(bv->rand,&t));
      if (i>=low && i<high) {
#if defined(PETSC_USE_COMPLEX)
        px[i-low] = PetscCMPLX(PetscSqrtReal(-2.0*PetscLogReal(PetscRealPart(s)))*PetscCosReal(2.0*PETSC_PI*PetscRealPart(t)),PetscSqrtReal(-2.0*PetscLogReal(PetscImaginaryPart(s)))*PetscCosReal(2.0*PETSC_PI*PetscImaginaryPart(t)));
#else
        px[i-low] = PetscSqrtReal(-2.0*PetscLogReal(s))*PetscCosReal(2.0*PETSC_PI*t);
#endif
      }
    }
    PetscCall(VecRestoreArray(x,&px));
  } else PetscCall(VecSetRandomNormal(x,bv->rand,w1,w2));
  PetscCall(BVRestoreColumn(bv,k,&x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode BVSetRandomSignColumn_Private(BV bv,PetscInt k)
{
  PetscInt       i,low,high;
  PetscScalar    *px,t;
  Vec            x;

  PetscFunctionBegin;
  PetscCall(BVGetColumn(bv,k,&x));
  PetscCall(VecGetOwnershipRange(x,&low,&high));
  if (bv->rrandom) {  /* generate the same vector irrespective of number of processes */
    PetscCall(VecGetArray(x,&px));
    for (i=0;i<bv->N;i++) {
      PetscCall(PetscRandomGetValue(bv->rand,&t));
      if (i>=low && i<high) px[i-low] = (PetscRealPart(t)<0.5)? -1.0: 1.0;
    }
    PetscCall(VecRestoreArray(x,&px));
  } else {
    PetscCall(VecSetRandom(x,bv->rand));
    PetscCall(VecGetArray(x,&px));
    for (i=low;i<high;i++) {
      px[i-low] = (PetscRealPart(px[i-low])<0.5)? -1.0: 1.0;
    }
    PetscCall(VecRestoreArray(x,&px));
  }
  PetscCall(BVRestoreColumn(bv,k,&x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetRandom - Set the columns of a `BV` to random numbers.

   Logically Collective

   Input Parameter:
.  bv - basis vectors

   Note:
   All active columns (except the leading ones) are modified.

   Level: advanced

.seealso: [](sec:bv), `BVSetRandomContext()`, `BVSetRandomColumn()`, `BVSetRandomNormal()`, `BVSetRandomSign()`, `BVSetRandomCond()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVSetRandom(BV bv)
{
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCall(BVGetRandomContext(bv,&bv->rand));
  PetscCall(PetscLogEventBegin(BV_SetRandom,bv,0,0,0));
  for (k=bv->l;k<bv->k;k++) PetscCall(BVSetRandomColumn_Private(bv,k));
  PetscCall(PetscLogEventEnd(BV_SetRandom,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetRandomColumn - Set one column of a `BV` to random numbers.

   Logically Collective

   Input Parameters:
+  bv - basis vectors
-  j  - index of column to be set

   Level: advanced

.seealso: [](sec:bv), `BVSetRandomContext()`, `BVSetRandom()`, `BVSetRandomNormal()`, `BVSetRandomCond()`
@*/
PetscErrorCode BVSetRandomColumn(BV bv,PetscInt j)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscCheck(j>=0 && j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %" PetscInt_FMT ", the number of columns is %" PetscInt_FMT,j,bv->m);

  PetscCall(BVGetRandomContext(bv,&bv->rand));
  PetscCall(PetscLogEventBegin(BV_SetRandom,bv,0,0,0));
  PetscCall(BVSetRandomColumn_Private(bv,j));
  PetscCall(PetscLogEventEnd(BV_SetRandom,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetRandomNormal - Set the columns of a `BV` to random numbers with a normal
   distribution.

   Logically Collective

   Input Parameter:
.  bv - basis vectors

   Notes:
   All active columns (except the leading ones) are modified.

   Other functions such as `BVSetRandom()`, `BVSetRandomColumn()`, and `BVSetRandomCond()`
   produce random numbers with a uniform distribution. This function returns values
   that fit a normal distribution (Gaussian).

   Developer Notes:
   The current implementation obtains each of the columns by applying the Box-Muller
   transform on two random vectors with uniformly distributed entries.

   Level: advanced

.seealso: [](sec:bv), `BVSetRandomContext()`, `BVSetRandom()`, `BVSetRandomColumn()`, `BVSetRandomCond()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVSetRandomNormal(BV bv)
{
  PetscInt       k;
  Vec            w1=NULL,w2=NULL;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCall(BVGetRandomContext(bv,&bv->rand));
  if (!bv->rrandom) {
    PetscCall(BVCreateVec(bv,&w1));
    PetscCall(BVCreateVec(bv,&w2));
  }
  PetscCall(PetscLogEventBegin(BV_SetRandom,bv,0,0,0));
  for (k=bv->l;k<bv->k;k++) PetscCall(BVSetRandomNormalColumn_Private(bv,k,w1,w2));
  PetscCall(PetscLogEventEnd(BV_SetRandom,bv,0,0,0));
  if (!bv->rrandom) {
    PetscCall(VecDestroy(&w1));
    PetscCall(VecDestroy(&w2));
  }
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetRandomSign - Set the entries of a `BV` to values 1 or -1 with equal
   probability.

   Logically Collective

   Input Parameter:
.  bv - basis vectors

   Notes:
   All active columns (except the leading ones) are modified.

   This function is used, e.g., in contour integral methods when estimating
   the number of eigenvalues enclosed by the contour via an unbiased
   estimator of $\operatorname{tr}(f(A))$ {cite:p}`Bai96`.

   Developer Notes:
   The current implementation obtains random numbers and then replaces them
   with -1 or 1 depending on the value being less than 0.5 or not.

   Level: advanced

.seealso: [](sec:bv), `BVSetRandomContext()`, `BVSetRandom()`, `BVSetRandomColumn()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVSetRandomSign(BV bv)
{
  PetscScalar    low,high;
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCall(BVGetRandomContext(bv,&bv->rand));
  PetscCall(PetscRandomGetInterval(bv->rand,&low,&high));
  PetscCheck(PetscRealPart(low)==0.0 && PetscRealPart(high)==1.0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"The PetscRandom object in the BV must have interval [0,1]");
  PetscCall(PetscLogEventBegin(BV_SetRandom,bv,0,0,0));
  for (k=bv->l;k<bv->k;k++) PetscCall(BVSetRandomSignColumn_Private(bv,k));
  PetscCall(PetscLogEventEnd(BV_SetRandom,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetRandomCond - Set the columns of a `BV` to random numbers, in a way that
   the generated matrix has a given condition number.

   Logically Collective

   Input Parameters:
+  bv    - basis vectors
-  condn - condition number

   Note:
   All active columns (except the leading ones) are modified.

   Level: advanced

.seealso: [](sec:bv), `BVSetRandomContext()`, `BVSetRandom()`, `BVSetRandomColumn()`, `BVSetRandomNormal()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVSetRandomCond(BV bv,PetscReal condn)
{
  PetscInt       k,i;
  PetscScalar    *eig,*d;
  DS             ds;
  Mat            A,X,Xt,M,G;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCall(BVGetRandomContext(bv,&bv->rand));
  PetscCall(PetscLogEventBegin(BV_SetRandom,bv,0,0,0));
  /* B = rand(n,k) */
  for (k=bv->l;k<bv->k;k++) PetscCall(BVSetRandomColumn_Private(bv,k));
  PetscCall(DSCreate(PetscObjectComm((PetscObject)bv),&ds));
  PetscCall(DSSetType(ds,DSHEP));
  PetscCall(DSAllocate(ds,bv->m));
  PetscCall(DSSetDimensions(ds,bv->k,bv->l,bv->k));
  /* [V,S] = eig(B'*B) */
  PetscCall(DSGetMat(ds,DS_MAT_A,&A));
  PetscCall(BVDot(bv,bv,A));
  PetscCall(DSRestoreMat(ds,DS_MAT_A,&A));
  PetscCall(PetscMalloc1(bv->m,&eig));
  PetscCall(DSSolve(ds,eig,NULL));
  PetscCall(DSSynchronize(ds,eig,NULL));
  PetscCall(DSVectors(ds,DS_MAT_X,NULL,NULL));
  /* M = diag(linspace(1/condn,1,n)./sqrt(diag(S)))' */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,bv->k,bv->k,NULL,&M));
  PetscCall(MatZeroEntries(M));
  PetscCall(MatDenseGetArray(M,&d));
  for (i=0;i<bv->k;i++) d[i+i*bv->m] = (1.0/condn+(1.0-1.0/condn)/(bv->k-1)*i)/PetscSqrtScalar(eig[i]);
  PetscCall(MatDenseRestoreArray(M,&d));
  /* G = X*M*X' */
  PetscCall(DSGetMat(ds,DS_MAT_X,&X));
  PetscCall(MatTranspose(X,MAT_INITIAL_MATRIX,&Xt));
  PetscCall(MatProductCreate(Xt,M,NULL,&G));
  PetscCall(MatProductSetType(G,MATPRODUCT_PtAP));
  PetscCall(MatProductSetFromOptions(G));
  PetscCall(MatProductSymbolic(G));
  PetscCall(MatProductNumeric(G));
  PetscCall(MatProductClear(G));
  PetscCall(DSRestoreMat(ds,DS_MAT_X,&X));
  PetscCall(MatDestroy(&Xt));
  PetscCall(MatDestroy(&M));
  /* B = B*G */
  PetscCall(BVMultInPlace(bv,G,bv->l,bv->k));
  PetscCall(MatDestroy(&G));
  PetscCall(PetscFree(eig));
  PetscCall(DSDestroy(&ds));
  PetscCall(PetscLogEventEnd(BV_SetRandom,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMatMult - Computes the matrix-vector product for each column, $Y=AV$.

   Neighbor-wise Collective

   Input Parameters:
+  V - basis vectors context
-  A - the matrix

   Output Parameter:
.  Y - the result

   Notes:
   Both `V` and `Y` must be distributed in the same manner. Only active columns
   (excluding the leading ones) are processed.
   In the result `Y`, columns are overwritten starting from the leading ones.
   The number of active columns in `V` and `Y` should match, although they need
   not be the same columns.

   It is possible to choose whether the computation is done column by column
   or as a `Mat`-`Mat` product, see `BVSetMatMultMethod()`.

   Level: intermediate

.seealso: [](sec:bv), `BVCopy()`, `BVSetActiveColumns()`, `BVMatMultColumn()`, `BVMatMultTranspose()`, `BVMatMultHermitianTranspose()`, `BVSetMatMultMethod()`
@*/
PetscErrorCode BVMatMult(BV V,Mat A,BV Y)
{
  PetscInt       M,N,m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  BVCheckOp(V,1,matmult);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidType(A,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidType(Y,3);
  BVCheckSizes(Y,3);
  PetscCheckSameComm(V,1,A,2);
  PetscCheckSameTypeAndComm(V,1,Y,3);

  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));
  PetscCheck(M==Y->N,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_INCOMP,"Mismatching row dimension A %" PetscInt_FMT ", Y %" PetscInt_FMT,M,Y->N);
  PetscCheck(m==Y->n,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_INCOMP,"Mismatching local row dimension A %" PetscInt_FMT ", Y %" PetscInt_FMT,m,Y->n);
  PetscCheck(N==V->N,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching column dimension A %" PetscInt_FMT ", V %" PetscInt_FMT,N,V->N);
  PetscCheck(n==V->n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local column dimension A %" PetscInt_FMT ", V %" PetscInt_FMT,n,V->n);
  PetscCheck(V->k-V->l==Y->k-Y->l,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Y has %" PetscInt_FMT " active columns, should match %" PetscInt_FMT " active columns in V",Y->k-Y->l,V->k-V->l);

  PetscCall(PetscLogEventBegin(BV_MatMult,V,A,Y,0));
  PetscUseTypeMethod(V,matmult,A,Y);
  PetscCall(PetscLogEventEnd(BV_MatMult,V,A,Y,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)Y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMatMultTranspose - Computes the matrix-vector product with the transpose
   of a matrix for each column, $Y=A^T V$.

   Neighbor-wise Collective

   Input Parameters:
+  V - basis vectors context
-  A - the matrix

   Output Parameter:
.  Y - the result

   Notes:
   Both `V` and `Y` must be distributed in the same manner. Only active columns
   (excluding the leading ones) are processed.
   In the result `Y`, columns are overwritten starting from the leading ones.
   The number of active columns in `V` and `Y` should match, although they need
   not be the same columns.

   Currently implemented via `MatCreateTranspose()`.

   Level: intermediate

.seealso: [](sec:bv), `BVMatMult()`, `BVMatMultHermitianTranspose()`, `BVMatMultTransposeColumn()`
@*/
PetscErrorCode BVMatMultTranspose(BV V,Mat A,BV Y)
{
  PetscInt       M,N,m,n;
  Mat            AT;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidType(A,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidType(Y,3);
  BVCheckSizes(Y,3);
  PetscCheckSameComm(V,1,A,2);
  PetscCheckSameTypeAndComm(V,1,Y,3);

  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));
  PetscCheck(M==V->N,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching row dimension A %" PetscInt_FMT ", V %" PetscInt_FMT,M,V->N);
  PetscCheck(m==V->n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local row dimension A %" PetscInt_FMT ", V %" PetscInt_FMT,m,V->n);
  PetscCheck(N==Y->N,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_INCOMP,"Mismatching column dimension A %" PetscInt_FMT ", Y %" PetscInt_FMT,N,Y->N);
  PetscCheck(n==Y->n,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_INCOMP,"Mismatching local column dimension A %" PetscInt_FMT ", Y %" PetscInt_FMT,n,Y->n);
  PetscCheck(V->k-V->l==Y->k-Y->l,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Y has %" PetscInt_FMT " active columns, should match %" PetscInt_FMT " active columns in V",Y->k-Y->l,V->k-V->l);

  PetscCall(MatCreateTranspose(A,&AT));
  PetscCall(BVMatMult(V,AT,Y));
  PetscCall(MatDestroy(&AT));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMatMultHermitianTranspose - Computes the matrix-vector product with the
   conjugate transpose of a matrix for each column, $Y=A^* V$.

   Neighbor-wise Collective

   Input Parameters:
+  V - basis vectors context
-  A - the matrix

   Output Parameter:
.  Y - the result

   Note:
   Both `V` and `Y` must be distributed in the same manner. Only active columns
   (excluding the leading ones) are processed.
   In the result `Y`, columns are overwritten starting from the leading ones.
   The number of active columns in `V` and `Y` should match, although they need
   not be the same columns.

   Currently implemented via `MatCreateHermitianTranspose()`.

   Level: intermediate

.seealso: [](sec:bv), `BVMatMult()`, `BVMatMultTranspose()`, `BVMatMultHermitianTransposeColumn()`
@*/
PetscErrorCode BVMatMultHermitianTranspose(BV V,Mat A,BV Y)
{
  PetscInt       j,M,N,m,n;
  Vec            v,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidType(A,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidType(Y,3);
  BVCheckSizes(Y,3);
  PetscCheckSameComm(V,1,A,2);
  PetscCheckSameTypeAndComm(V,1,Y,3);

  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));
  PetscCheck(M==V->N,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching row dimension A %" PetscInt_FMT ", V %" PetscInt_FMT,M,V->N);
  PetscCheck(m==V->n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local row dimension A %" PetscInt_FMT ", V %" PetscInt_FMT,m,V->n);
  PetscCheck(N==Y->N,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_INCOMP,"Mismatching column dimension A %" PetscInt_FMT ", Y %" PetscInt_FMT,N,Y->N);
  PetscCheck(n==Y->n,PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_INCOMP,"Mismatching local column dimension A %" PetscInt_FMT ", Y %" PetscInt_FMT,n,Y->n);
  PetscCheck(V->k-V->l==Y->k-Y->l,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Y has %" PetscInt_FMT " active columns, should match %" PetscInt_FMT " active columns in V",Y->k-Y->l,V->k-V->l);

  /* TODO: recover this code if PETSc ever gets MATPRODUCT_AhB done
  PetscCall(MatCreateHermitianTranspose(A,&AH));
  PetscCall(BVMatMult(V,AH,Y));
  PetscCall(MatDestroy(&AH));
  */
  for (j=0;j<V->k-V->l;j++) {
    PetscCall(BVGetColumn(V,V->l+j,&v));
    PetscCall(BVGetColumn(Y,Y->l+j,&y));
    PetscCall(MatMultHermitianTranspose(A,v,y));
    PetscCall(BVRestoreColumn(V,V->l+j,&v));
    PetscCall(BVRestoreColumn(Y,Y->l+j,&y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMatMultColumn - Computes the matrix-vector product for a specified
   column, storing the result in the next column $v_{j+1}=Av_j$.

   Neighbor-wise Collective

   Input Parameters:
+  V - basis vectors context
.  A - the matrix
-  j - the index of the column

   Level: intermediate

.seealso: [](sec:bv), `BVMatMult()`, `BVMatMultTransposeColumn()`, `BVMatMultHermitianTransposeColumn()`
@*/
PetscErrorCode BVMatMultColumn(BV V,Mat A,PetscInt j)
{
  Vec            vj,vj1;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(V,1,A,2);
  PetscValidLogicalCollectiveInt(V,j,3);
  PetscCheck(j>=0,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j+1<V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Result should go in index j+1=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j+1,V->m);

  PetscCall(PetscLogEventBegin(BV_MatMultVec,V,A,0,0));
  PetscCall(BVGetColumn(V,j,&vj));
  PetscCall(BVGetColumn(V,j+1,&vj1));
  PetscCall(MatMult(A,vj,vj1));
  PetscCall(BVRestoreColumn(V,j,&vj));
  PetscCall(BVRestoreColumn(V,j+1,&vj1));
  PetscCall(PetscLogEventEnd(BV_MatMultVec,V,A,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMatMultTransposeColumn - Computes the transpose matrix-vector product for a
   specified column, storing the result in the next column $v_{j+1}=A^T v_j$.

   Neighbor-wise Collective

   Input Parameters:
+  V - basis vectors context
.  A - the matrix
-  j - the index of the column

   Level: intermediate

.seealso: [](sec:bv), `BVMatMult()`, `BVMatMultColumn()`, `BVMatMultHermitianTransposeColumn()`
@*/
PetscErrorCode BVMatMultTransposeColumn(BV V,Mat A,PetscInt j)
{
  Vec            vj,vj1;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(V,1,A,2);
  PetscValidLogicalCollectiveInt(V,j,3);
  PetscCheck(j>=0,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j+1<V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Result should go in index j+1=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j+1,V->m);

  PetscCall(PetscLogEventBegin(BV_MatMultVec,V,A,0,0));
  PetscCall(BVGetColumn(V,j,&vj));
  PetscCall(BVGetColumn(V,j+1,&vj1));
  PetscCall(MatMultTranspose(A,vj,vj1));
  PetscCall(BVRestoreColumn(V,j,&vj));
  PetscCall(BVRestoreColumn(V,j+1,&vj1));
  PetscCall(PetscLogEventEnd(BV_MatMultVec,V,A,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMatMultHermitianTransposeColumn - Computes the conjugate-transpose matrix-vector
   product for a specified column, storing the result in the next column $v_{j+1}=A^* v_j$.

   Neighbor-wise Collective

   Input Parameters:
+  V - basis vectors context
.  A - the matrix
-  j - the index of the column

   Level: intermediate

.seealso: [](sec:bv), `BVMatMult()`, `BVMatMultColumn()`, `BVMatMultTransposeColumn()`, `BVMatMultHermitianTranspose()`
@*/
PetscErrorCode BVMatMultHermitianTransposeColumn(BV V,Mat A,PetscInt j)
{
  Vec            vj,vj1;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(V,1,A,2);
  PetscValidLogicalCollectiveInt(V,j,3);
  PetscCheck(j>=0,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j+1<V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Result should go in index j+1=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j+1,V->m);

  PetscCall(PetscLogEventBegin(BV_MatMultVec,V,A,0,0));
  PetscCall(BVGetColumn(V,j,&vj));
  PetscCall(BVGetColumn(V,j+1,&vj1));
  PetscCall(MatMultHermitianTranspose(A,vj,vj1));
  PetscCall(BVRestoreColumn(V,j,&vj));
  PetscCall(BVRestoreColumn(V,j+1,&vj1));
  PetscCall(PetscLogEventEnd(BV_MatMultVec,V,A,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}
