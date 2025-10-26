/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV operations involving global communication
*/

#include <slepc/private/bvimpl.h>      /*I "slepcbv.h" I*/

/*
  BVDot for the particular case of non-standard inner product with
  matrix B, which is assumed to be symmetric (or complex Hermitian)
*/
static inline PetscErrorCode BVDot_Private(BV X,BV Y,Mat M)
{
  PetscObjectId  idx,idy;
  PetscInt       i,j,jend,m;
  PetscScalar    *marray;
  PetscBool      symm=PETSC_FALSE;
  Mat            B;
  Vec            z;

  PetscFunctionBegin;
  BVCheckOp(Y,1,dotvec);
  PetscCall(MatGetSize(M,&m,NULL));
  PetscCall(MatDenseGetArray(M,&marray));
  PetscCall(PetscObjectGetId((PetscObject)X,&idx));
  PetscCall(PetscObjectGetId((PetscObject)Y,&idy));
  B = Y->matrix;
  Y->matrix = NULL;
  if (idx==idy) symm=PETSC_TRUE;  /* M=X'BX is symmetric */
  jend = X->k;
  for (j=X->l;j<jend;j++) {
    if (symm) Y->k = j+1;
    PetscCall(BVGetColumn(X->cached,j,&z));
    PetscUseTypeMethod(Y,dotvec,z,marray+j*m+Y->l);
    PetscCall(BVRestoreColumn(X->cached,j,&z));
    if (symm) {
      for (i=X->l;i<j;i++)
        marray[j+i*m] = PetscConj(marray[i+j*m]);
    }
  }
  PetscCall(MatDenseRestoreArray(M,&marray));
  Y->matrix = B;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDot - Computes the 'block-dot' product of two basis vectors objects.

   Collective

   Input Parameters:
+  X - first basis vectors
-  Y - second basis vectors

   Output Parameter:
.  M - the resulting matrix

   Notes:
   This is the generalization of `VecDot()` for a collection of vectors, $M = Y^*X$.
   The result is a matrix $M$ whose entry $m_{ij}$ is equal to $y_i^* x_j$ (where $y_i^*$
   denotes the conjugate transpose of $y_i$).

   If a non-standard inner product has been specified with `BVSetMatrix()`,
   then the result is $M = Y^*BX$. In this case, both `X` and `Y` must have
   the same associated matrix.

   On entry, `M` must be a sequential dense `Mat` with dimensions `m,n` at least, where
   `m` is the number of active columns of `Y` and `n` is the number of active columns of `X`.
   Only rows (resp. columns) of `M` starting from `ly` (resp. `lx`) are computed,
   where `ly` (resp. `lx`) is the number of leading columns of `Y` (resp. `X`).

   `X` and `Y` need not be different objects.

   Level: intermediate

.seealso: [](sec:bv), `BVDotVec()`, `BVDotColumn()`, `BVSetActiveColumns()`, `BVSetMatrix()`
@*/
PetscErrorCode BVDot(BV X,BV Y,Mat M)
{
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidHeaderSpecific(Y,BV_CLASSID,2);
  PetscValidHeaderSpecific(M,MAT_CLASSID,3);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  PetscValidType(Y,2);
  BVCheckSizes(Y,2);
  PetscValidType(M,3);
  PetscCheckSameTypeAndComm(X,1,Y,2);
  SlepcMatCheckSeq(M);

  PetscCall(MatGetSize(M,&m,&n));
  PetscCheck(m>=Y->k,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Mat argument has %" PetscInt_FMT " rows, should have at least %" PetscInt_FMT,m,Y->k);
  PetscCheck(n>=X->k,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Mat argument has %" PetscInt_FMT " columns, should have at least %" PetscInt_FMT,n,X->k);
  PetscCheck(X->n==Y->n,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %" PetscInt_FMT ", Y %" PetscInt_FMT,X->n,Y->n);
  PetscCheck(X->matrix==Y->matrix,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"X and Y must have the same inner product matrix");
  if (X->l==X->k || Y->l==Y->k) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscLogEventBegin(BV_Dot,X,Y,0,0));
  if (X->matrix) { /* non-standard inner product */
    /* compute BX first */
    PetscCall(BV_IPMatMultBV(X));
    if (X->vmm==BV_MATMULT_VECS) {
      /* perform computation column by column */
      PetscCall(BVDot_Private(X,Y,M));
    } else PetscUseTypeMethod(X->cached,dot,Y,M);
  } else PetscUseTypeMethod(X,dot,Y,M);
  PetscCall(PetscLogEventEnd(BV_Dot,X,Y,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDotVec - Computes multiple dot products of a vector against all the
   column vectors of a `BV`.

   Collective

   Input Parameters:
+  X - basis vectors
-  y - a vector

   Output Parameter:
.  m - an array where the result must be placed

   Notes:
   This is analogue to `VecMDot()`, but using `BV` to represent a collection
   of vectors. The result is $m = X^*y$, so $m_i$ is equal to $x_i^* y$. Note
   that here `X` is transposed as opposed to `BVDot()`.

   If a non-standard inner product has been specified with `BVSetMatrix()`,
   then the result is $m = X^*By$.

   The length of array `m` must be equal to the number of active columns of `X`
   minus the number of leading columns, i.e., the first entry of `m` is the
   product of the first non-leading column with `y`.

   Level: intermediate

.seealso: [](sec:bv), `BVDot()`, `BVDotColumn()`, `BVSetActiveColumns()`, `BVSetMatrix()`
@*/
PetscErrorCode BVDotVec(BV X,Vec y,PetscScalar m[])
{
  PetscInt       n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidHeaderSpecific(y,VEC_CLASSID,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  BVCheckOp(X,1,dotvec);
  PetscValidType(y,2);
  PetscCheckSameTypeAndComm(X,1,y,2);

  PetscCall(VecGetLocalSize(y,&n));
  PetscCheck(X->n==n,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %" PetscInt_FMT ", y %" PetscInt_FMT,X->n,n);

  PetscCall(PetscLogEventBegin(BV_DotVec,X,y,0,0));
  PetscUseTypeMethod(X,dotvec,y,m);
  PetscCall(PetscLogEventEnd(BV_DotVec,X,y,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDotVecBegin - Starts a split phase dot product computation.

   Input Parameters:
+  X - basis vectors
.  y - a vector
-  m - an array where the result will go (can be `NULL`)

   Note:
   Each call to `BVDotVecBegin()` should be paired with a call to `BVDotVecEnd()`.

   Level: advanced

.seealso: [](sec:bv), `BVDotVecEnd()`, `BVDotVec()`
@*/
PetscErrorCode BVDotVecBegin(BV X,Vec y,PetscScalar *m)
{
  PetscInt            i,n,nv;
  PetscSplitReduction *sr;
  MPI_Comm            comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidHeaderSpecific(y,VEC_CLASSID,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  PetscValidType(y,2);
  PetscCheckSameTypeAndComm(X,1,y,2);

  PetscCall(VecGetLocalSize(y,&n));
  PetscCheck(X->n==n,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %" PetscInt_FMT ", y %" PetscInt_FMT,X->n,n);

  if (X->ops->dotvec_begin) PetscUseTypeMethod(X,dotvec_begin,y,m);
  else {
    BVCheckOp(X,1,dotvec_local);
    nv = X->k-X->l;
    PetscCall(PetscObjectGetComm((PetscObject)X,&comm));
    PetscCall(PetscSplitReductionGet(comm,&sr));
    PetscCheck(sr->state==STATE_BEGIN,PetscObjectComm((PetscObject)X),PETSC_ERR_ORDER,"Called before all VecxxxEnd() called");
    for (i=0;i<nv;i++) {
      if (sr->numopsbegin+i >= sr->maxops) PetscCall(PetscSplitReductionExtend(sr));
      sr->reducetype[sr->numopsbegin+i] = PETSC_SR_REDUCE_SUM;
      sr->invecs[sr->numopsbegin+i]     = (void*)X;
    }
    PetscCall(PetscLogEventBegin(BV_DotVec,X,y,0,0));
    PetscUseTypeMethod(X,dotvec_local,y,sr->lvalues+sr->numopsbegin);
    sr->numopsbegin += nv;
    PetscCall(PetscLogEventEnd(BV_DotVec,X,y,0,0));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDotVecEnd - Ends a split phase dot product computation.

   Input Parameters:
+  X - basis vectors
.  y - a vector
-  m - an array where the result will go

   Note:
   Each call to `BVDotVecBegin()` should be paired with a call to `BVDotVecEnd()`.

   Level: advanced

.seealso: [](sec:bv), `BVDotVecBegin()`, `BVDotVec()`
@*/
PetscErrorCode BVDotVecEnd(BV X,Vec y,PetscScalar *m)
{
  PetscInt            i,nv;
  PetscSplitReduction *sr;
  MPI_Comm            comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidType(X,1);
  BVCheckSizes(X,1);

  if (X->ops->dotvec_end) PetscUseTypeMethod(X,dotvec_end,y,m);
  else {
    nv = X->k-X->l;
    PetscCall(PetscObjectGetComm((PetscObject)X,&comm));
    PetscCall(PetscSplitReductionGet(comm,&sr));
    PetscCall(PetscSplitReductionEnd(sr));

    PetscCheck(sr->numopsend<sr->numopsbegin,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Called VecxxxEnd() more times than VecxxxBegin()");
    PetscCheck((void*)X==sr->invecs[sr->numopsend],PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Called BVxxxEnd() in a different order or with a different BV than BVxxxBegin()");
    PetscCheck(sr->reducetype[sr->numopsend]==PETSC_SR_REDUCE_SUM,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Wrong type of reduction");
    for (i=0;i<nv;i++) m[i] = sr->gvalues[sr->numopsend++];

    /* Finished getting all the results so reset to no outstanding requests */
    if (sr->numopsend == sr->numopsbegin) {
      sr->state       = STATE_BEGIN;
      sr->numopsend   = 0;
      sr->numopsbegin = 0;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDotColumn - Computes multiple dot products of a column against all the
   previous columns of a `BV`.

   Collective

   Input Parameters:
+  X - basis vectors
-  j - the column index

   Output Parameter:
.  q - an array where the result must be placed

   Notes:
   This operation is equivalent to `BVDotVec()` but it uses column `j` of `X`
   rather than taking a `Vec` as an argument. The number of active columns of
   `X` is set to `j` before the computation, and restored afterwards.
   If `X` has leading columns specified, then these columns do not participate
   in the computation. Therefore, the length of array `q` must be at least `j`
   minus the number of leading columns.

   Developer Notes:
   If `q` is `NULL`, then the result is written in position `nc+l` of the internal
   buffer vector, see `BVGetBufferVec()`.

   Level: intermediate

.seealso: [](sec:bv), `BVDot()`, `BVDotVec()`, `BVSetActiveColumns()`, `BVSetMatrix()`
@*/
PetscErrorCode BVDotColumn(BV X,PetscInt j,PetscScalar q[])
{
  PetscInt       ksave;
  Vec            y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(X,j,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  BVCheckOp(X,1,dotvec);

  PetscCheck(j>=0,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j<X->m,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j,X->m);

  PetscCall(PetscLogEventBegin(BV_DotVec,X,0,0,0));
  ksave = X->k;
  X->k = j;
  if (!q && !X->buffer) PetscCall(BVGetBufferVec(X,&X->buffer));
  PetscCall(BVGetColumn(X,j,&y));
  PetscUseTypeMethod(X,dotvec,y,q);
  PetscCall(BVRestoreColumn(X,j,&y));
  X->k = ksave;
  PetscCall(PetscLogEventEnd(BV_DotVec,X,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDotColumnBegin - Starts a split phase dot product computation.

   Input Parameters:
+  X - basis vectors
.  j - the column index
-  m - an array where the result will go (can be `NULL`)

   Note:
   Each call to `BVDotColumnBegin()` should be paired with a call to `BVDotColumnEnd()`.

   Level: advanced

.seealso: [](sec:bv), `BVDotColumnEnd()`, `BVDotColumn()`
@*/
PetscErrorCode BVDotColumnBegin(BV X,PetscInt j,PetscScalar *m)
{
  PetscInt            i,nv,ksave;
  PetscSplitReduction *sr;
  MPI_Comm            comm;
  Vec                 y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(X,j,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);

  PetscCheck(j>=0,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j<X->m,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j,X->m);
  ksave = X->k;
  X->k = j;
  PetscCall(BVGetColumn(X,j,&y));

  if (X->ops->dotvec_begin) PetscUseTypeMethod(X,dotvec_begin,y,m);
  else {
    BVCheckOp(X,1,dotvec_local);
    nv = X->k-X->l;
    PetscCall(PetscObjectGetComm((PetscObject)X,&comm));
    PetscCall(PetscSplitReductionGet(comm,&sr));
    PetscCheck(sr->state==STATE_BEGIN,PetscObjectComm((PetscObject)X),PETSC_ERR_ORDER,"Called before all VecxxxEnd() called");
    for (i=0;i<nv;i++) {
      if (sr->numopsbegin+i >= sr->maxops) PetscCall(PetscSplitReductionExtend(sr));
      sr->reducetype[sr->numopsbegin+i] = PETSC_SR_REDUCE_SUM;
      sr->invecs[sr->numopsbegin+i]     = (void*)X;
    }
    PetscCall(PetscLogEventBegin(BV_DotVec,X,0,0,0));
    PetscUseTypeMethod(X,dotvec_local,y,sr->lvalues+sr->numopsbegin);
    sr->numopsbegin += nv;
    PetscCall(PetscLogEventEnd(BV_DotVec,X,0,0,0));
  }
  PetscCall(BVRestoreColumn(X,j,&y));
  X->k = ksave;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDotColumnEnd - Ends a split phase dot product computation.

   Input Parameters:
+  X - basis vectors
.  j - the column index
-  m - an array where the result will go

   Notes:
   Each call to `BVDotColumnBegin()` should be paired with a call to `BVDotColumnEnd()`.

   Level: advanced

.seealso: [](sec:bv), `BVDotColumnBegin()`, `BVDotColumn()`
@*/
PetscErrorCode BVDotColumnEnd(BV X,PetscInt j,PetscScalar *m)
{
  PetscInt            i,nv,ksave;
  PetscSplitReduction *sr;
  MPI_Comm            comm;
  Vec                 y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(X,j,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);

  PetscCheck(j>=0,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j<X->m,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j,X->m);
  ksave = X->k;
  X->k = j;

  if (X->ops->dotvec_end) {
    PetscCall(BVGetColumn(X,j,&y));
    PetscUseTypeMethod(X,dotvec_end,y,m);
    PetscCall(BVRestoreColumn(X,j,&y));
  } else {
    nv = X->k-X->l;
    PetscCall(PetscObjectGetComm((PetscObject)X,&comm));
    PetscCall(PetscSplitReductionGet(comm,&sr));
    PetscCall(PetscSplitReductionEnd(sr));

    PetscCheck(sr->numopsend<sr->numopsbegin,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Called VecxxxEnd() more times than VecxxxBegin()");
    PetscCheck((void*)X==sr->invecs[sr->numopsend],PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Called BVxxxEnd() in a different order or with a different BV than BVxxxBegin()");
    PetscCheck(sr->reducetype[sr->numopsend]==PETSC_SR_REDUCE_SUM,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Wrong type of reduction");
    for (i=0;i<nv;i++) m[i] = sr->gvalues[sr->numopsend++];

    /* Finished getting all the results so reset to no outstanding requests */
    if (sr->numopsend == sr->numopsbegin) {
      sr->state       = STATE_BEGIN;
      sr->numopsend   = 0;
      sr->numopsbegin = 0;
    }
  }
  X->k = ksave;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode BVNorm_Private(BV bv,Vec z,NormType type,PetscReal *val)
{
  PetscScalar    p;

  PetscFunctionBegin;
  PetscCall(BV_IPMatMult(bv,z));
  PetscCall(VecDot(bv->Bx,z,&p));
  PetscCall(BV_SafeSqrt(bv,p,val));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode BVNorm_Begin_Private(BV bv,Vec z,NormType type,PetscReal *val)
{
  PetscScalar    p;

  PetscFunctionBegin;
  PetscCall(BV_IPMatMult(bv,z));
  PetscCall(VecDotBegin(bv->Bx,z,&p));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode BVNorm_End_Private(BV bv,Vec z,NormType type,PetscReal *val)
{
  PetscScalar    p;

  PetscFunctionBegin;
  PetscCall(VecDotEnd(bv->Bx,z,&p));
  PetscCall(BV_SafeSqrt(bv,p,val));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVNorm - Computes the matrix norm of the `BV`.

   Collective

   Input Parameters:
+  bv   - basis vectors
-  type - the norm type

   Output Parameter:
.  val  - the norm

   Notes:
   All active columns (except the leading ones) are considered as a matrix.
   The allowed norms are `NORM_1`, `NORM_FROBENIUS`, and `NORM_INFINITY`.

   This operation fails if a non-standard inner product has been
   specified with `BVSetMatrix()`.

   Level: intermediate

.seealso: [](sec:bv), `BVNormVec()`, `BVNormColumn()`, `BVNormalize()`, `BVSetActiveColumns()`, `BVSetMatrix()`
@*/
PetscErrorCode BVNorm(BV bv,NormType type,PetscReal *val)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveEnum(bv,type,2);
  PetscAssertPointer(val,3);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCheck(type!=NORM_2 && type!=NORM_1_AND_2,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");
  PetscCheck(!bv->matrix,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Matrix norm not available for non-standard inner product");

  PetscCall(PetscLogEventBegin(BV_Norm,bv,0,0,0));
  PetscUseTypeMethod(bv,norm,-1,type,val);
  PetscCall(PetscLogEventEnd(BV_Norm,bv,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVNormVec - Computes the norm of a given vector.

   Collective

   Input Parameters:
+  bv   - basis vectors
.  v    - the vector
-  type - the norm type

   Output Parameter:
.  val  - the norm

   Notes:
   This is the analogue of `BVNormColumn()` but for a vector that is not in the `BV`.
   If a non-standard inner product has been specified with `BVSetMatrix()`,
   then the returned value is $\sqrt{v^*Bv}$, where $B$ is the inner product
   matrix (argument `type` is ignored). Otherwise, `VecNorm()` is called.

   Level: developer

.seealso: [](sec:bv), `BVNorm()`, `BVNormColumn()`, `BVSetMatrix()`
@*/
PetscErrorCode BVNormVec(BV bv,Vec v,NormType type,PetscReal *val)
{
  PetscInt       n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(v,VEC_CLASSID,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscAssertPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscValidType(v,2);
  PetscCheckSameComm(bv,1,v,2);

  PetscCheck(type!=NORM_1_AND_2 || bv->matrix,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  PetscCall(PetscLogEventBegin(BV_NormVec,bv,0,0,0));
  if (bv->matrix) { /* non-standard inner product */
    PetscCall(VecGetLocalSize(v,&n));
    PetscCheck(bv->n==n,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension bv %" PetscInt_FMT ", v %" PetscInt_FMT,bv->n,n);
    PetscCall(BVNorm_Private(bv,v,type,val));
  } else PetscCall(VecNorm(v,type,val));
  PetscCall(PetscLogEventEnd(BV_NormVec,bv,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVNormVecBegin - Starts a split phase norm computation.

   Input Parameters:
+  bv   - basis vectors
.  v    - the vector
.  type - the norm type
-  val  - the norm

   Note:
   Each call to `BVNormVecBegin()` should be paired with a call to `BVNormVecEnd()`.

   Level: advanced

.seealso: [](sec:bv), `BVNormVecEnd()`, `BVNormVec()`
@*/
PetscErrorCode BVNormVecBegin(BV bv,Vec v,NormType type,PetscReal *val)
{
  PetscInt       n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(v,VEC_CLASSID,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscAssertPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscValidType(v,2);
  PetscCheckSameTypeAndComm(bv,1,v,2);

  PetscCheck(type!=NORM_1_AND_2 || bv->matrix,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  PetscCall(PetscLogEventBegin(BV_NormVec,bv,0,0,0));
  if (bv->matrix) { /* non-standard inner product */
    PetscCall(VecGetLocalSize(v,&n));
    PetscCheck(bv->n==n,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension bv %" PetscInt_FMT ", v %" PetscInt_FMT,bv->n,n);
    PetscCall(BVNorm_Begin_Private(bv,v,type,val));
  } else PetscCall(VecNormBegin(v,type,val));
  PetscCall(PetscLogEventEnd(BV_NormVec,bv,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVNormVecEnd - Ends a split phase norm computation.

   Input Parameters:
+  bv   - basis vectors
.  v    - the vector
.  type - the norm type
-  val  - the norm

   Note:
   Each call to `BVNormVecBegin()` should be paired with a call to `BVNormVecEnd()`.

   Level: advanced

.seealso: [](sec:bv), `BVNormVecBegin()`, `BVNormVec()`
@*/
PetscErrorCode BVNormVecEnd(BV bv,Vec v,NormType type,PetscReal *val)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscAssertPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCheck(type!=NORM_1_AND_2,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  if (bv->matrix) PetscCall(BVNorm_End_Private(bv,v,type,val));  /* non-standard inner product */
  else PetscCall(VecNormEnd(v,type,val));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVNormColumn - Computes the vector norm of a selected column.

   Collective

   Input Parameters:
+  bv   - basis vectors
.  j    - index of column to be used
-  type - the norm type

   Output Parameter:
.  val  - the norm

   Notes:
   The norm of $v_j$ is computed (`NORM_1`, `NORM_2`, or `NORM_INFINITY`).
   If a non-standard inner product has been specified with `BVSetMatrix()`,
   then the returned value is $\sqrt{v_j^*B v_j}$,
   where $B$ is the inner product matrix (argument `type` is ignored).

   Level: intermediate

.seealso: [](sec:bv), `BVNorm()`, `BVNormVec()`, `BVNormalize()`, `BVSetActiveColumns()`, `BVSetMatrix()`
@*/
PetscErrorCode BVNormColumn(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  Vec            z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscAssertPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCheck(j>=0 && j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %" PetscInt_FMT ", the number of columns is %" PetscInt_FMT,j,bv->m);
  PetscCheck(type!=NORM_1_AND_2 || bv->matrix,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  PetscCall(PetscLogEventBegin(BV_NormVec,bv,0,0,0));
  if (bv->matrix) { /* non-standard inner product */
    PetscCall(BVGetColumn(bv,j,&z));
    PetscCall(BVNorm_Private(bv,z,type,val));
    PetscCall(BVRestoreColumn(bv,j,&z));
  } else PetscUseTypeMethod(bv,norm,j,type,val);
  PetscCall(PetscLogEventEnd(BV_NormVec,bv,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVNormColumnBegin - Starts a split phase norm computation.

   Input Parameters:
+  bv   - basis vectors
.  j    - index of column to be used
.  type - the norm type
-  val  - the norm

   Note:
   Each call to `BVNormColumnBegin()` should be paired with a call to `BVNormColumnEnd()`.

   Level: advanced

.seealso: [](sec:bv), `BVNormColumnEnd()`, `BVNormColumn()`
@*/
PetscErrorCode BVNormColumnBegin(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscSplitReduction *sr;
  PetscReal           lresult;
  MPI_Comm            comm;
  Vec                 z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscAssertPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCheck(j>=0 && j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %" PetscInt_FMT ", the number of columns is %" PetscInt_FMT,j,bv->m);
  PetscCheck(type!=NORM_1_AND_2 || bv->matrix,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  PetscCall(PetscLogEventBegin(BV_NormVec,bv,0,0,0));
  PetscCall(BVGetColumn(bv,j,&z));
  if (bv->matrix) PetscCall(BVNorm_Begin_Private(bv,z,type,val)); /* non-standard inner product */
  else if (bv->ops->norm_begin) PetscUseTypeMethod(bv,norm_begin,j,type,val);
  else {
    BVCheckOp(bv,1,norm_local);
    PetscCall(PetscObjectGetComm((PetscObject)z,&comm));
    PetscCall(PetscSplitReductionGet(comm,&sr));
    PetscCheck(sr->state==STATE_BEGIN,PetscObjectComm((PetscObject)bv),PETSC_ERR_ORDER,"Called before all VecxxxEnd() called");
    if (sr->numopsbegin >= sr->maxops) PetscCall(PetscSplitReductionExtend(sr));
    sr->invecs[sr->numopsbegin] = (void*)bv;
    PetscUseTypeMethod(bv,norm_local,j,type,&lresult);
    if (type == NORM_2) lresult = lresult*lresult;
    if (type == NORM_MAX) sr->reducetype[sr->numopsbegin] = PETSC_SR_REDUCE_MAX;
    else sr->reducetype[sr->numopsbegin] = PETSC_SR_REDUCE_SUM;
    sr->lvalues[sr->numopsbegin++] = lresult;
  }
  PetscCall(BVRestoreColumn(bv,j,&z));
  PetscCall(PetscLogEventEnd(BV_NormVec,bv,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVNormColumnEnd - Ends a split phase norm computation.

   Input Parameters:
+  bv   - basis vectors
.  j    - index of column to be used
.  type - the norm type
-  val  - the norm

   Note:
   Each call to `BVNormColumnBegin()` should be paired with a call to `BVNormColumnEnd()`.

   Level: advanced

.seealso: [](sec:bv), `BVNormColumnBegin()`, `BVNormColumn()`
@*/
PetscErrorCode BVNormColumnEnd(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscSplitReduction *sr;
  MPI_Comm            comm;
  Vec                 z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscAssertPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCheck(type!=NORM_1_AND_2,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  PetscCall(BVGetColumn(bv,j,&z));
  if (bv->matrix) PetscCall(BVNorm_End_Private(bv,z,type,val)); /* non-standard inner product */
  else if (bv->ops->norm_end) PetscUseTypeMethod(bv,norm_end,j,type,val);
  else {
    PetscCall(PetscObjectGetComm((PetscObject)z,&comm));
    PetscCall(PetscSplitReductionGet(comm,&sr));
    PetscCall(PetscSplitReductionEnd(sr));

    PetscCheck(sr->numopsend<sr->numopsbegin,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Called VecxxxEnd() more times then VecxxxBegin()");
    PetscCheck((void*)bv==sr->invecs[sr->numopsend],PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Called VecxxxEnd() in a different order or with a different vector than VecxxxBegin()");
    PetscCheck(sr->reducetype[sr->numopsend]==PETSC_SR_REDUCE_MAX || type!=NORM_MAX,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Called BVNormEnd(,NORM_MAX,) on a reduction started with VecDotBegin() or NORM_1 or NORM_2");
    *val = PetscRealPart(sr->gvalues[sr->numopsend++]);
    if (type == NORM_2) *val = PetscSqrtReal(*val);
    if (sr->numopsend == sr->numopsbegin) {
      sr->state       = STATE_BEGIN;
      sr->numopsend   = 0;
      sr->numopsbegin = 0;
    }
  }
  PetscCall(BVRestoreColumn(bv,j,&z));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVNormalize - Normalize all columns (starting from the leading ones).

   Collective

   Input Parameters:
+  bv   - basis vectors
-  eigi - (optional) imaginary parts of eigenvalues

   Notes:
   On output, all columns will have unit norm. The normalization is done with
   respect to the 2-norm, or to the $B$-norm if a non-standard inner product has
   been specified with `BVSetMatrix()`, see `BVNormColumn()`.

   If the optional argument `eigi` is passed (taken into account only in real
   scalars) it is interpreted as the imaginary parts of the eigenvalues and
   the `BV` is supposed to contain the corresponding eigenvectors. Suppose the
   first three values are `eigi` = $\{ 0, \alpha, -\alpha \}$, then the first column
   is normalized as usual, but the second and third ones are normalized assuming
   that they contain the real and imaginary parts of a complex conjugate pair of
   eigenvectors.

   If `eigi` is passed, the inner-product matrix is ignored.

   If there are leading columns, they are not modified (are assumed to be already
   normalized).

   Level: intermediate

.seealso: [](sec:bv), `BVNormColumn()`, `BVSetMatrix()`
@*/
PetscErrorCode BVNormalize(BV bv,PetscScalar eigi[])
{
  PetscReal      norm;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  PetscCall(PetscLogEventBegin(BV_Normalize,bv,0,0,0));
  if (bv->matrix && !eigi) {
    for (i=bv->l;i<bv->k;i++) {
      PetscCall(BVNormColumn(bv,i,NORM_2,&norm));
      PetscCall(BVScaleColumn(bv,i,1.0/norm));
    }
  } else PetscTryTypeMethod(bv,normalize,eigi);
  PetscCall(PetscLogEventEnd(BV_Normalize,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute Y^H*A*X: right part column by column (with MatMult) and bottom
  part row by row (with MatMultHermitianTranspose); result placed in marray[*,ldm]
*/
static inline PetscErrorCode BVMatProject_Vec(BV X,Mat A,BV Y,PetscScalar *marray,PetscInt ldm,PetscBool symm)
{
  PetscInt       i,j,lx,ly,kx,ky,ulim;
  Vec            z,f;

  PetscFunctionBegin;
  lx = X->l; kx = X->k;
  ly = Y->l; ky = Y->k;
  PetscCall(BVCreateVec(X,&f));
  BVCheckOp(Y,3,dotvec);
  for (j=lx;j<kx;j++) {
    PetscCall(BVGetColumn(X,j,&z));
    PetscCall(MatMult(A,z,f));
    PetscCall(BVRestoreColumn(X,j,&z));
    ulim = PetscMin(ly+(j-lx)+1,ky);
    Y->l = 0; Y->k = ulim;
    PetscUseTypeMethod(Y,dotvec,f,marray+j*ldm);
    if (symm) {
      for (i=0;i<j;i++) marray[j+i*ldm] = PetscConj(marray[i+j*ldm]);
    }
  }
  if (!symm) {
    BVCheckOp(X,1,dotvec);
    PetscCall(BV_AllocateCoeffs(Y));
    for (j=ly;j<ky;j++) {
      PetscCall(BVGetColumn(Y,j,&z));
      PetscCall(MatMultHermitianTranspose(A,z,f));
      PetscCall(BVRestoreColumn(Y,j,&z));
      ulim = PetscMin(lx+(j-ly),kx);
      X->l = 0; X->k = ulim;
      PetscUseTypeMethod(X,dotvec,f,Y->h);
      for (i=0;i<ulim;i++) marray[j+i*ldm] = PetscConj(Y->h[i]);
    }
  }
  PetscCall(VecDestroy(&f));
  X->l = lx; X->k = kx;
  Y->l = ly; Y->k = ky;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute Y^H*A*X= [   --   | Y0'*W1 ]
                   [ Y1'*W0 | Y1'*W1 ]
  Allocates auxiliary BV to store the result of A*X, then one BVDot
  call for top-right part and another one for bottom part;
  result placed in marray[*,ldm]
*/
static inline PetscErrorCode BVMatProject_MatMult(BV X,Mat A,BV Y,PetscScalar *marray,PetscInt ldm)
{
  PetscInt          j,lx,ly,kx,ky;
  const PetscScalar *harray;
  Mat               H;
  BV                W;

  PetscFunctionBegin;
  lx = X->l; kx = X->k;
  ly = Y->l; ky = Y->k;
  PetscCall(BVDuplicate(X,&W));
  X->l = 0; X->k = kx;
  W->l = 0; W->k = kx;
  PetscCall(BVMatMult(X,A,W));

  /* top-right part, Y0'*AX1 */
  if (ly>0 && lx<kx) {
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ly,kx,NULL,&H));
    W->l = lx; W->k = kx;
    Y->l = 0;  Y->k = ly;
    PetscCall(BVDot(W,Y,H));
    PetscCall(MatDenseGetArrayRead(H,&harray));
    for (j=lx;j<kx;j++) PetscCall(PetscArraycpy(marray+j*ldm,harray+j*ly,ly));
    PetscCall(MatDenseRestoreArrayRead(H,&harray));
    PetscCall(MatDestroy(&H));
  }

  /* bottom part, Y1'*AX */
  if (kx>0 && ly<ky) {
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&H));
    W->l = 0;  W->k = kx;
    Y->l = ly; Y->k = ky;
    PetscCall(BVDot(W,Y,H));
    PetscCall(MatDenseGetArrayRead(H,&harray));
    for (j=0;j<kx;j++) PetscCall(PetscArraycpy(marray+j*ldm+ly,harray+j*ky+ly,ky-ly));
    PetscCall(MatDenseRestoreArrayRead(H,&harray));
    PetscCall(MatDestroy(&H));
  }
  PetscCall(BVDestroy(&W));
  X->l = lx; X->k = kx;
  Y->l = ly; Y->k = ky;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute Y^H*A*X= [   --   | Y0'*W1 ]
                   [ Y1'*W0 | Y1'*W1 ]
  First stage: allocate auxiliary BV to store A*X1, one BVDot for right part;
  Second stage: resize BV to accommodate A'*Y1, then call BVDot for transpose of
  bottom-left part; result placed in marray[*,ldm]
*/
static inline PetscErrorCode BVMatProject_MatMult_2(BV X,Mat A,BV Y,PetscScalar *marray,PetscInt ldm,PetscBool symm)
{
  PetscInt          i,j,lx,ly,kx,ky;
  const PetscScalar *harray;
  Mat               H;
  BV                W;

  PetscFunctionBegin;
  lx = X->l; kx = X->k;
  ly = Y->l; ky = Y->k;

  /* right part, Y'*AX1 */
  PetscCall(BVDuplicateResize(X,kx-lx,&W));
  if (ky>0 && lx<kx) {
    PetscCall(BVMatMult(X,A,W));
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ky,kx-lx,NULL,&H));
    Y->l = 0; Y->k = ky;
    PetscCall(BVDot(W,Y,H));
    PetscCall(MatDenseGetArrayRead(H,&harray));
    for (j=lx;j<kx;j++) PetscCall(PetscArraycpy(marray+j*ldm,harray+(j-lx)*ky,ky));
    PetscCall(MatDenseRestoreArrayRead(H,&harray));
    PetscCall(MatDestroy(&H));
  }

  /* bottom-left part, Y1'*AX0 */
  if (lx>0 && ly<ky) {
    if (symm) {
      /* do not compute, just copy symmetric elements */
      for (i=ly;i<ky;i++) {
        for (j=0;j<lx;j++) marray[i+j*ldm] = PetscConj(marray[j+i*ldm]);
      }
    } else {
      PetscCall(BVResize(W,ky-ly,PETSC_FALSE));
      Y->l = ly; Y->k = ky;
      PetscCall(BVMatMultHermitianTranspose(Y,A,W));
      PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,lx,ky-ly,NULL,&H));
      X->l = 0; X->k = lx;
      PetscCall(BVDot(W,X,H));
      PetscCall(MatDenseGetArrayRead(H,&harray));
      for (i=0;i<ky-ly;i++) {
        for (j=0;j<lx;j++) {
          marray[i+j*ldm+ly] = PetscConj(harray[j+i*(ky-ly)]);
        }
      }
      PetscCall(MatDenseRestoreArrayRead(H,&harray));
      PetscCall(MatDestroy(&H));
    }
  }
  PetscCall(BVDestroy(&W));
  X->l = lx; X->k = kx;
  Y->l = ly; Y->k = ky;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute Y^H*X = [   --   | Y0'*X1 ]     (X contains A*X):
                  [ Y1'*X0 | Y1'*X1 ]
  one BVDot call for top-right part and another one for bottom part;
  result placed in marray[*,ldm]
*/
static inline PetscErrorCode BVMatProject_Dot(BV X,BV Y,PetscScalar *marray,PetscInt ldm)
{
  PetscInt          j,lx,ly,kx,ky;
  const PetscScalar *harray;
  Mat               H;

  PetscFunctionBegin;
  lx = X->l; kx = X->k;
  ly = Y->l; ky = Y->k;

  /* top-right part, Y0'*X1 */
  if (ly>0 && lx<kx) {
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ly,kx,NULL,&H));
    X->l = lx; X->k = kx;
    Y->l = 0;  Y->k = ly;
    PetscCall(BVDot(X,Y,H));
    PetscCall(MatDenseGetArrayRead(H,&harray));
    for (j=lx;j<kx;j++) PetscCall(PetscArraycpy(marray+j*ldm,harray+j*ly,ly));
    PetscCall(MatDenseRestoreArrayRead(H,&harray));
    PetscCall(MatDestroy(&H));
  }

  /* bottom part, Y1'*X */
  if (kx>0 && ly<ky) {
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&H));
    X->l = 0;  X->k = kx;
    Y->l = ly; Y->k = ky;
    PetscCall(BVDot(X,Y,H));
    PetscCall(MatDenseGetArrayRead(H,&harray));
    for (j=0;j<kx;j++) PetscCall(PetscArraycpy(marray+j*ldm+ly,harray+j*ky+ly,ky-ly));
    PetscCall(MatDenseRestoreArrayRead(H,&harray));
    PetscCall(MatDestroy(&H));
  }
  X->l = lx; X->k = kx;
  Y->l = ly; Y->k = ky;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMatProject - Computes the projection of a matrix onto a subspace.

   Collective

   Input Parameters:
+  X - basis vectors
.  A - (optional) matrix to be projected
-  Y - left basis vectors, can be equal to `X`

   Output Parameter:
.  M - the resulting matrix

   Notes:
   If `A=NULL`, then it is assumed that `X` already contains $AX$.

   This operation is similar to `BVDot()`, with important differences.
   The goal is to compute the matrix resulting from the orthogonal projection
   of `A` onto the subspace spanned by the columns of `X`, $M = X^*AX$, or the
   oblique projection onto `X` along `Y`, $M = Y^*AX$.

   A difference with respect to `BVDot()` is that the standard inner product
   is always used, regardless of a non-standard inner product being specified
   with `BVSetMatrix()`.

   On entry, `M` must be a sequential dense `Mat` with dimensions `ky,kx` at least,
   where `ky` (resp. `kx`) is the number of active columns of `Y` (resp. `X`).
   Another difference with respect to `BVDot()` is that all entries of `M` are
   computed except the leading `ly,lx` part, where `ly` (resp. `lx`) is the
   number of leading columns of `Y` (resp. `X`). Hence, the leading columns of
   `X` and `Y` participate in the computation, as opposed to `BVDot()`.
   The leading part of `M` is assumed to be already available from previous
   computations.

   In the orthogonal projection case, `Y=X`, some computation can be saved if
   `A` is real symmetric (or complex Hermitian). In order to exploit this
   property, the symmetry flag of `A` must be set with `MatSetOption()`.

   Level: intermediate

.seealso: [](sec:bv), `BVDot()`, `BVSetActiveColumns()`, `BVSetMatrix()`, `MatSetOption()`
@*/
PetscErrorCode BVMatProject(BV X,Mat A,BV Y,Mat M)
{
  PetscBool      set,flg,symm=PETSC_FALSE;
  PetscInt       m,n,ldm;
  PetscScalar    *marray;
  Mat            Xmatrix,Ymatrix;
  PetscObjectId  idx,idy;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidHeaderSpecific(M,MAT_CLASSID,4);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  if (A) {
    PetscValidType(A,2);
    PetscCheckSameComm(X,1,A,2);
  }
  PetscValidType(Y,3);
  BVCheckSizes(Y,3);
  PetscValidType(M,4);
  PetscCheckSameTypeAndComm(X,1,Y,3);
  SlepcMatCheckSeq(M);

  PetscCall(MatGetSize(M,&m,&n));
  PetscCall(MatDenseGetLDA(M,&ldm));
  PetscCheck(m>=Y->k,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Matrix M has %" PetscInt_FMT " rows, should have at least %" PetscInt_FMT,m,Y->k);
  PetscCheck(n>=X->k,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Matrix M has %" PetscInt_FMT " columns, should have at least %" PetscInt_FMT,n,X->k);
  PetscCheck(X->n==Y->n,PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %" PetscInt_FMT ", Y %" PetscInt_FMT,X->n,Y->n);

  PetscCall(PetscLogEventBegin(BV_MatProject,X,A,Y,0));
  /* temporarily set standard inner product */
  Xmatrix = X->matrix;
  Ymatrix = Y->matrix;
  X->matrix = Y->matrix = NULL;

  PetscCall(PetscObjectGetId((PetscObject)X,&idx));
  PetscCall(PetscObjectGetId((PetscObject)Y,&idy));
  if (A && idx==idy) { /* check symmetry of M=X'AX */
    PetscCall(MatIsHermitianKnown(A,&set,&flg));
    symm = set? flg: PETSC_FALSE;
  }

  PetscCall(MatDenseGetArray(M,&marray));

  if (A) {
    if (X->vmm==BV_MATMULT_VECS) {
      /* perform computation column by column */
      PetscCall(BVMatProject_Vec(X,A,Y,marray,ldm,symm));
    } else {
      /* use BVMatMult, then BVDot */
      PetscCall(MatHasOperation(A,MATOP_MULT_TRANSPOSE,&flg));
      if (symm || (flg && X->l>=X->k/2 && Y->l>=Y->k/2)) PetscCall(BVMatProject_MatMult_2(X,A,Y,marray,ldm,symm));
      else PetscCall(BVMatProject_MatMult(X,A,Y,marray,ldm));
    }
  } else {
    /* use BVDot on subblocks */
    PetscCall(BVMatProject_Dot(X,Y,marray,ldm));
  }

  PetscCall(MatDenseRestoreArray(M,&marray));
  PetscCall(PetscLogEventEnd(BV_MatProject,X,A,Y,0));
  /* restore non-standard inner product */
  X->matrix = Xmatrix;
  Y->matrix = Ymatrix;
  PetscFunctionReturn(PETSC_SUCCESS);
}
