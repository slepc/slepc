/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic BV routines
*/

#include <slepc/private/bvimpl.h>      /*I "slepcbv.h" I*/

PetscBool         BVRegisterAllCalled = PETSC_FALSE;
PetscFunctionList BVList = NULL;

/*@
   BVSetType - Selects the type for the `BV` object.

   Logically Collective

   Input Parameters:
+  bv   - the basis vectors context
-  type - a known type

   Options Database Key:
.  -bv_type \<type\> - sets `BV` type

   Level: intermediate

.seealso: [](sec:bv), `BVGetType()`
@*/
PetscErrorCode BVSetType(BV bv,BVType type)
{
  PetscErrorCode (*r)(BV);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)bv,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscStrcmp(type,BVTENSOR,&match));
  PetscCheck(!match,PetscObjectComm((PetscObject)bv),PETSC_ERR_ORDER,"Use BVCreateTensor() to create a BV of type tensor");

  PetscCall(PetscFunctionListFind(BVList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested BV type %s",type);

  PetscTryTypeMethod(bv,destroy);
  PetscCall(PetscMemzero(bv->ops,sizeof(struct _BVOps)));

  PetscCall(PetscObjectChangeTypeName((PetscObject)bv,type));
  if (bv->n < 0 && bv->N < 0) {
    bv->ops->create = r;
  } else {
    PetscCall(PetscLogEventBegin(BV_Create,bv,0,0,0));
    PetscCall((*r)(bv));
    PetscCall(PetscLogEventEnd(BV_Create,bv,0,0,0));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetType - Gets the `BV` type name (as a string) from the `BV` context.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  type - name of the type of basis vectors

   Level: intermediate

.seealso: [](sec:bv), `BVSetType()`
@*/
PetscErrorCode BVGetType(BV bv,BVType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)bv)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetSizes - Sets the local and global sizes, and the number of columns.

   Collective

   Input Parameters:
+  bv - the basis vectors
.  n  - the local size (or `PETSC_DECIDE` to have it set)
.  N  - the global size (or `PETSC_DECIDE`)
-  m  - the number of columns

   Notes:
   `n` and `N` cannot be both `PETSC_DECIDE`.
   If one process calls this with `N` equal to `PETSC_DECIDE` then all processes must,
   otherwise the program will hang.

   Level: beginner

.seealso: [](sec:bv), `BVSetSizesFromVec()`, `BVGetSizes()`, `BVResize()`
@*/
PetscErrorCode BVSetSizes(BV bv,PetscInt n,PetscInt N,PetscInt m)
{
  PetscInt       ma;
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (N >= 0) PetscValidLogicalCollectiveInt(bv,N,3);
  PetscValidLogicalCollectiveInt(bv,m,4);
  PetscCheck(N<0 || n<=N,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Local size %" PetscInt_FMT " cannot be larger than global size %" PetscInt_FMT,n,N);
  PetscCheck(m>0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Number of columns %" PetscInt_FMT " must be positive",m);
  PetscCheck((bv->n<0 && bv->N<0) || (bv->n==n && bv->N==N),PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Cannot change/reset vector sizes to %" PetscInt_FMT " local %" PetscInt_FMT " global after previously setting them to %" PetscInt_FMT " local %" PetscInt_FMT " global",n,N,bv->n,bv->N);
  PetscCheck(bv->m<=0 || bv->m==m,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Cannot change the number of columns to %" PetscInt_FMT " after previously setting it to %" PetscInt_FMT "; use BVResize()",m,bv->m);
  PetscCheck(!bv->map,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Vector layout was already defined by a previous call to BVSetSizes/FromVec");
  bv->n = n;
  bv->N = N;
  bv->m = m;
  bv->k = m;
  /* create layout and get actual dimensions */
  PetscCall(PetscLayoutCreate(PetscObjectComm((PetscObject)bv),&bv->map));
  PetscCall(PetscLayoutSetSize(bv->map,bv->N));
  PetscCall(PetscLayoutSetLocalSize(bv->map,bv->n));
  PetscCall(PetscLayoutSetUp(bv->map));
  PetscCall(PetscLayoutGetSize(bv->map,&bv->N));
  PetscCall(PetscLayoutGetLocalSize(bv->map,&bv->n));
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)bv),&size));
  PetscCall(BVSetVecType(bv,(size==1)?VECSEQ:VECMPI));
  if (bv->matrix) {  /* check compatible dimensions of user-provided matrix */
    PetscCall(MatGetLocalSize(bv->matrix,&ma,NULL));
    PetscCheck(bv->n==ma,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Local dimension %" PetscInt_FMT " does not match that of matrix given at BVSetMatrix %" PetscInt_FMT,bv->n,ma);
  }
  if (bv->ops->create) {
    PetscCall(PetscLogEventBegin(BV_Create,bv,0,0,0));
    PetscUseTypeMethod(bv,create);
    PetscCall(PetscLogEventEnd(BV_Create,bv,0,0,0));
    bv->ops->create = NULL;
    bv->defersfo = PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetSizesFromVec - Sets the local and global sizes, and the number of columns.
   Local and global sizes are specified indirectly by passing a template vector.

   Collective

   Input Parameters:
+  bv - the basis vectors
.  t  - the template vector
-  m  - the number of columns

   Notes:
   Apart from the local and global sizes, the template vector `t` is also used
   to get the vector type, which will be set as the `BV` vector type with
   `BVSetVecType()`. This is relevant in cases where computation must be done
   on the GPU.

   If `t` is a nested vector (type `VECNEST`), then the vector type will be
   obtained from the first subvector. The rationale is that `BV` will operate
   internally with regular vectors, even though the user problem is defined via
   nested matrices and vectors.

   Level: beginner

.seealso: [](sec:bv), `BVSetSizes()`, `BVGetSizes()`, `BVResize()`, `BVSetVecType()`
@*/
PetscErrorCode BVSetSizesFromVec(BV bv,Vec t,PetscInt m)
{
  PetscInt       ma;
  PetscLayout    map;
  VecType        vtype;
  PetscBool      nest;
  Vec            v;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(t,VEC_CLASSID,2);
  PetscCheckSameComm(bv,1,t,2);
  PetscValidLogicalCollectiveInt(bv,m,3);
  PetscCheck(m>0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Number of columns %" PetscInt_FMT " must be positive",m);
  PetscCheck(!bv->map,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Vector layout was already defined by a previous call to BVSetSizes/FromVec");
  PetscCall(PetscObjectTypeCompare((PetscObject)t,VECNEST,&nest));
  if (nest) {
    PetscCall(VecNestGetSubVec(t,0,&v));
    PetscCall(VecGetType(v,&vtype));
  } else PetscCall(VecGetType(t,&vtype));
  PetscCall(BVSetVecType(bv,vtype));
  PetscCall(VecGetLayout(t,&map));
  PetscCall(PetscLayoutReference(map,&bv->map));
  PetscCall(VecGetSize(t,&bv->N));
  PetscCall(VecGetLocalSize(t,&bv->n));
  if (bv->matrix) {  /* check compatible dimensions of user-provided matrix */
    PetscCall(MatGetLocalSize(bv->matrix,&ma,NULL));
    PetscCheck(bv->n==ma,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Local dimension %" PetscInt_FMT " does not match that of matrix given at BVSetMatrix %" PetscInt_FMT,bv->n,ma);
  }
  bv->m = m;
  bv->k = m;
  if (bv->ops->create) {
    PetscUseTypeMethod(bv,create);
    bv->ops->create = NULL;
    bv->defersfo = PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetSizes - Returns the local and global sizes, and the number of columns.

   Not Collective

   Input Parameter:
.  bv - the basis vectors

   Output Parameters:
+  n  - the local size
.  N  - the global size
-  m  - the number of columns

   Note:
   Normal usage requires that `bv` has already been given its sizes, otherwise
   the call fails. However, this function can also be used to determine if
   a `BV` object has been initialized completely (sizes and type). For this,
   call with `n=NULL` and `N=NULL`, then a return value of `m=0` indicates that
   the `BV` object is not ready for use yet.

   Level: beginner

.seealso: [](sec:bv), `BVSetSizes()`, `BVSetSizesFromVec()`
@*/
PetscErrorCode BVGetSizes(BV bv,PetscInt *n,PetscInt *N,PetscInt *m)
{
  PetscFunctionBegin;
  if (!bv) {
    if (m && !n && !N) *m = 0;
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (n || N) BVCheckSizes(bv,1);
  if (n) *n = bv->n;
  if (N) *N = bv->N;
  if (m) *m = bv->m;
  if (m && !n && !N && !((PetscObject)bv)->type_name) *m = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetNumConstraints - Set the number of constraints.

   Logically Collective

   Input Parameters:
+  V  - basis vectors
-  nc - number of constraints

   Notes:
   This function sets the number of constraints to `nc` and marks all remaining
   columns as regular. Normal usage would be to call `BVInsertConstraints()` instead.

   If `nc` is smaller than the previously set value, then some of the constraints
   are discarded. In particular, using `nc=0` removes all constraints preserving
   the content of regular columns.

   Level: developer

.seealso: [](sec:bv), `BVInsertConstraints()`
@*/
PetscErrorCode BVSetNumConstraints(BV V,PetscInt nc)
{
  PetscInt       total,diff,i;
  Vec            x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(V,nc,2);
  PetscCheck(nc>=0,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Number of constraints (given %" PetscInt_FMT ") cannot be negative",nc);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscCheck(V->ci[0]==-V->nc-1 && V->ci[1]==-V->nc-1,PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Cannot call BVSetNumConstraints after BVGetColumn");

  diff = nc-V->nc;
  if (!diff) PetscFunctionReturn(PETSC_SUCCESS);
  total = V->nc+V->m;
  PetscCheck(total-nc>0,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Not enough columns for the given nc value");
  if (diff<0) {  /* lessen constraints, shift contents of BV */
    for (i=0;i<V->m;i++) {
      PetscCall(BVGetColumn(V,i,&x));
      PetscCall(BVGetColumn(V,i+diff,&y));
      PetscCall(VecCopy(x,y));
      PetscCall(BVRestoreColumn(V,i,&x));
      PetscCall(BVRestoreColumn(V,i+diff,&y));
    }
  }
  V->nc = nc;
  V->ci[0] = -V->nc-1;
  V->ci[1] = -V->nc-1;
  V->m = total-nc;
  V->l = PetscMin(V->l,V->m);
  V->k = PetscMin(V->k,V->m);
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetNumConstraints - Returns the number of constraints.

   Not Collective

   Input Parameter:
.  bv - the basis vectors

   Output Parameter:
.  nc - the number of constraints

   Level: advanced

.seealso: [](sec:bv), `BVGetSizes()`, `BVInsertConstraints()`
@*/
PetscErrorCode BVGetNumConstraints(BV bv,PetscInt *nc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(nc,2);
  *nc = bv->nc;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVResize - Change the number of columns.

   Collective

   Input Parameters:
+  bv   - the basis vectors
.  m    - the new number of columns
-  copy - a flag indicating whether current values should be kept

   Note:
   Internal storage is reallocated. If the `copy` flag is set to true, then
   the contents are copied to the leading part of the new space.

   Level: advanced

.seealso: [](sec:bv), `BVSetSizes()`, `BVSetSizesFromVec()`
@*/
PetscErrorCode BVResize(BV bv,PetscInt m,PetscBool copy)
{
  PetscScalar       *array;
  const PetscScalar *omega;
  Vec               v;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,m,2);
  PetscValidLogicalCollectiveBool(bv,copy,3);
  PetscValidType(bv,1);
  PetscCheck(m>0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Number of columns %" PetscInt_FMT " must be positive",m);
  PetscCheck(!bv->nc || bv->issplit,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Cannot resize a BV with constraints");
  if (bv->m == m) PetscFunctionReturn(PETSC_SUCCESS);
  BVCheckOp(bv,1,resize);

  PetscCall(PetscLogEventBegin(BV_Create,bv,0,0,0));
  PetscUseTypeMethod(bv,resize,m,copy);
  PetscCall(VecDestroy(&bv->buffer));
  PetscCall(BVDestroy(&bv->cached));
  PetscCall(PetscFree2(bv->h,bv->c));
  if (bv->omega) {
    if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
      PetscCall(VecCreateSeqCUDA(PETSC_COMM_SELF,m,&v));
#else
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_PLIB,"Something wrong happened");
#endif
    } else if (bv->hip) {
#if defined(PETSC_HAVE_HIP)
      PetscCall(VecCreateSeqHIP(PETSC_COMM_SELF,m,&v));
#else
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_PLIB,"Something wrong happened");
#endif
    } else PetscCall(VecCreateSeq(PETSC_COMM_SELF,m,&v));
    if (copy) {
      PetscCall(VecGetArray(v,&array));
      PetscCall(VecGetArrayRead(bv->omega,&omega));
      PetscCall(PetscArraycpy(array,omega,PetscMin(m,bv->m)));
      PetscCall(VecRestoreArrayRead(bv->omega,&omega));
      PetscCall(VecRestoreArray(v,&array));
    } else PetscCall(VecSet(v,1.0));
    PetscCall(VecDestroy(&bv->omega));
    bv->omega = v;
  }
  bv->m = m;
  bv->k = PetscMin(bv->k,m);
  bv->l = PetscMin(bv->l,m);
  PetscCall(PetscLogEventEnd(BV_Create,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetActiveColumns - Specify the columns that will be involved in operations.

   Logically Collective

   Input Parameters:
+  bv - the basis vectors context
.  l  - number of leading columns
-  k  - number of active columns

   Notes:
   In operations such as `BVMult()` or `BVDot()`, only the first `k` columns are
   considered. This is useful when the `BV` is filled from left to right, so
   the last `m-k` columns do not have relevant information.

   Also in operations such as `BVMult()` or `BVDot()`, the first `l` columns are
   normally not included in the computation. See the manual page of each
   operation.

   In orthogonalization operations, the first `l` columns are treated
   differently, they participate in the orthogonalization but the computed
   coefficients are not stored.

   Use `PETSC_CURRENT` to leave any of the values unchanged. Use `PETSC_DETERMINE`
   to set `l` to the minimum value (`0`) and `k` to the maximum (`m`).

   Level: intermediate

.seealso: [](sec:bv), `BVGetActiveColumns()`, `BVSetSizes()`
@*/
PetscErrorCode BVSetActiveColumns(BV bv,PetscInt l,PetscInt k)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,l,2);
  PetscValidLogicalCollectiveInt(bv,k,3);
  BVCheckSizes(bv,1);
  if (PetscUnlikely(k == PETSC_DETERMINE)) {
    bv->k = bv->m;
  } else if (k != PETSC_CURRENT) {
    PetscCheck(k>=0 && k<=bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of k (%" PetscInt_FMT "). Must be between 0 and m (%" PetscInt_FMT ")",k,bv->m);
    bv->k = k;
  }
  if (PetscUnlikely(l == PETSC_DETERMINE)) {
    bv->l = 0;
  } else if (l != PETSC_CURRENT) {
    PetscCheck(l>=0 && l<=bv->k,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of l (%" PetscInt_FMT "). Must be between 0 and k (%" PetscInt_FMT ")",l,bv->k);
    bv->l = l;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetActiveColumns - Returns the current active dimensions.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameters:
+  l  - number of leading columns
-  k  - number of active columns

   Level: intermediate

.seealso: [](sec:bv), `BVSetActiveColumns()`
@*/
PetscErrorCode BVGetActiveColumns(BV bv,PetscInt *l,PetscInt *k)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (l) *l = bv->l;
  if (k) *k = bv->k;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetMatrix - Specifies the inner product to be used in orthogonalization.

   Collective

   Input Parameters:
+  bv    - the basis vectors context
.  B     - a Hermitian matrix (may be `NULL`)
-  indef - a flag indicating if the matrix is indefinite

   Notes:
   This is used to specify a non-standard inner product, whose matrix
   representation is given by `B`. Then, all inner products required during
   orthogonalization are computed as $(x,y)_B=y^*Bx$ rather than the
   standard form $(x,y)=y^*x$.

   Matrix `B` must be real symmetric (or complex Hermitian). A genuine inner
   product requires that `B` is also positive (semi-)definite. However, we
   also allow for an indefinite `B` (setting `indef=PETSC_TRUE`), in which
   case the orthogonalization uses an indefinite inner product.

   This affects operations `BVDot()`, `BVNorm()`, `BVOrthogonalize()`, and variants.

   Setting `B=NULL` has the same effect as if the identity matrix was passed.

   Level: advanced

.seealso: [](sec:bv), `BVGetMatrix()`, `BVDot()`, `BVNorm()`, `BVOrthogonalize()`, `BVSetDefiniteTolerance()`
@*/
PetscErrorCode BVSetMatrix(BV bv,Mat B,PetscBool indef)
{
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveBool(bv,indef,3);
  if (B!=bv->matrix || (B && ((PetscObject)B)->id!=((PetscObject)bv->matrix)->id) || indef!=bv->indef) {
    if (B) {
      PetscValidHeaderSpecific(B,MAT_CLASSID,2);
      PetscCall(MatGetLocalSize(B,&m,&n));
      PetscCheck(m==n,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_SIZ,"Matrix must be square");
      PetscCheck(!bv->m || bv->n==n,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension BV %" PetscInt_FMT ", Mat %" PetscInt_FMT,bv->n,n);
    }
    if (B) PetscCall(PetscObjectReference((PetscObject)B));
    PetscCall(MatDestroy(&bv->matrix));
    bv->matrix = B;
    bv->indef  = indef;
    PetscCall(PetscObjectStateIncrease((PetscObject)bv));
    if (bv->Bx) PetscCall(PetscObjectStateIncrease((PetscObject)bv->Bx));
    if (bv->cached) PetscCall(PetscObjectStateIncrease((PetscObject)bv->cached));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetMatrix - Retrieves the matrix representation of the inner product.

   Not Collective

   Input Parameter:
.  bv    - the basis vectors context

   Output Parameters:
+  B     - the matrix of the inner product (may be `NULL`)
-  indef - the flag indicating if the matrix is indefinite

   Level: advanced

.seealso: [](sec:bv), `BVSetMatrix()`
@*/
PetscErrorCode BVGetMatrix(BV bv,Mat *B,PetscBool *indef)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (B)     *B     = bv->matrix;
  if (indef) *indef = bv->indef;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVApplyMatrix - Multiplies a vector by the matrix representation of the
   inner product.

   Neighbor-wise Collective

   Input Parameters:
+  bv - the basis vectors context
-  x  - the vector

   Output Parameter:
.  y  - the result

   Note:
   If no matrix was specified this function copies the vector.

   Level: advanced

.seealso: [](sec:bv), `BVSetMatrix()`, `BVApplyMatrixBV()`
@*/
PetscErrorCode BVApplyMatrix(BV bv,Vec x,Vec y)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  if (bv->matrix) {
    PetscCall(BV_IPMatMult(bv,x));
    PetscCall(VecCopy(bv->Bx,y));
  } else PetscCall(VecCopy(x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVApplyMatrixBV - Multiplies the `BV` vectors by the matrix representation
   of the inner product.

   Neighbor-wise Collective

   Input Parameter:
.  X - the basis vectors context

   Output Parameter:
.  Y - the basis vectors to store the result (optional)

   Note:
   This function computes $Y = B X$, where $B$ is the matrix given with
   `BVSetMatrix()`. This operation is computed as in `BVMatMult()`.
   If no matrix was specified, then it just copies $Y = X$.

   If no `Y` is given, the result is stored internally in the cached `BV`.

   Level: developer

.seealso: [](sec:bv), `BVSetMatrix()`, `BVApplyMatrix()`, `BVMatMult()`, `BVGetCachedBV()`
@*/
PetscErrorCode BVApplyMatrixBV(BV X,BV Y)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  if (Y) {
    PetscValidHeaderSpecific(Y,BV_CLASSID,2);
    if (X->matrix) PetscCall(BVMatMult(X,X->matrix,Y));
    else PetscCall(BVCopy(X,Y));
  } else PetscCall(BV_IPMatMultBV(X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetSignature - Sets the signature matrix to be used in orthogonalization.

   Logically Collective

   Input Parameters:
+  bv    - the basis vectors context
-  omega - a vector representing the diagonal of the signature matrix

   Note:
   The signature matrix $\Omega = V^*B V$ is relevant only for an indefinite $B$.

   Level: developer

.seealso: [](sec:bv), `BVSetMatrix()`, `BVGetSignature()`
@*/
PetscErrorCode BVSetSignature(BV bv,Vec omega)
{
  PetscInt          i,n;
  const PetscScalar *pomega;
  PetscScalar       *intern;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscValidHeaderSpecific(omega,VEC_CLASSID,2);
  PetscValidType(omega,2);

  PetscCall(VecGetSize(omega,&n));
  PetscCheck(n==bv->k,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_SIZ,"Vec argument has %" PetscInt_FMT " elements, should be %" PetscInt_FMT,n,bv->k);
  PetscCall(BV_AllocateSignature(bv));
  if (bv->indef) {
    PetscCall(VecGetArrayRead(omega,&pomega));
    PetscCall(VecGetArray(bv->omega,&intern));
    for (i=0;i<n;i++) intern[bv->nc+i] = pomega[i];
    PetscCall(VecRestoreArray(bv->omega,&intern));
    PetscCall(VecRestoreArrayRead(omega,&pomega));
  } else PetscCall(PetscInfo(bv,"Ignoring signature because BV is not indefinite\n"));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetSignature - Retrieves the signature matrix from last orthogonalization.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  omega - a vector representing the diagonal of the signature matrix

   Note:
   The signature matrix $\Omega = V^*B V$ is relevant only for an indefinite $B$.

   Level: developer

.seealso: [](sec:bv), `BVSetMatrix()`, `BVSetSignature()`
@*/
PetscErrorCode BVGetSignature(BV bv,Vec omega)
{
  PetscInt          i,n;
  PetscScalar       *pomega;
  const PetscScalar *intern;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscValidHeaderSpecific(omega,VEC_CLASSID,2);
  PetscValidType(omega,2);

  PetscCall(VecGetSize(omega,&n));
  PetscCheck(n==bv->k,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_SIZ,"Vec argument has %" PetscInt_FMT " elements, should be %" PetscInt_FMT,n,bv->k);
  if (bv->indef && bv->omega) {
    PetscCall(VecGetArray(omega,&pomega));
    PetscCall(VecGetArrayRead(bv->omega,&intern));
    for (i=0;i<n;i++) pomega[i] = intern[bv->nc+i];
    PetscCall(VecRestoreArrayRead(bv->omega,&intern));
    PetscCall(VecRestoreArray(omega,&pomega));
  } else PetscCall(VecSet(omega,1.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetBufferVec - Attach a vector object to be used as buffer space for
   several operations.

   Collective

   Input Parameters:
+  bv     - the basis vectors context)
-  buffer - the vector

   Notes:
   Use `BVGetBufferVec()` to retrieve the vector (for example, to free it
   at the end of the computations).

   The vector must be sequential of length $(n_c+m)m$, where $m$ is the number
   of columns of `bv` and $n_c$ is the number of constraints.

   Level: developer

.seealso: [](sec:bv), `BVGetBufferVec()`, `BVSetSizes()`, `BVGetNumConstraints()`
@*/
PetscErrorCode BVSetBufferVec(BV bv,Vec buffer)
{
  PetscInt       ld,n;
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(buffer,VEC_CLASSID,2);
  BVCheckSizes(bv,1);
  PetscCall(VecGetSize(buffer,&n));
  ld = bv->m+bv->nc;
  PetscCheck(n==ld*bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_SIZ,"Buffer size must be %" PetscInt_FMT,ld*bv->m);
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)buffer),&size));
  PetscCheck(size==1,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Buffer must be a sequential vector");

  PetscCall(PetscObjectReference((PetscObject)buffer));
  PetscCall(VecDestroy(&bv->buffer));
  bv->buffer = buffer;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetBufferVec - Obtain the buffer vector associated with the `BV` object.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  buffer - vector

   Notes:
   The vector is created if not available previously. It is a sequential vector
   of length $(n_c+m) m$, where $m$ is the number of columns of `bv` and $n_c$ is the number
   of constraints.

   Developer Notes:
   The buffer vector is viewed as a column-major matrix with leading dimension equal
   to $n_c+m$, and $m$ columns at most. In the most common usage, it has the structure
   $$
     \left[\begin{array}{c|c}
       s & \begin{array}{c} C \\\hline H \end{array}
     \end{array}\right]
   $$
   where $H$ is an upper Hessenberg matrix of order $m (m-1)$, $C$ contains coefficients
   related to orthogonalization against constraints (first $n_c$ rows), and $s$ is the
   first column that contains scratch values computed during Gram-Schmidt
   orthogonalization. In particular, `BVDotColumn()` and `BVMultColumn()` use $s$ to
   store the coefficients.

   Level: developer

.seealso: [](sec:bv), `BVSetBufferVec()`, `BVSetSizes()`, `BVGetNumConstraints()`, `BVDotColumn()`, `BVMultColumn()`
@*/
PetscErrorCode BVGetBufferVec(BV bv,Vec *buffer)
{
  PetscInt       ld;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(buffer,2);
  BVCheckSizes(bv,1);
  if (!bv->buffer) {
    ld = bv->m+bv->nc;
    PetscCall(VecCreate(PETSC_COMM_SELF,&bv->buffer));
    PetscCall(VecSetSizes(bv->buffer,PETSC_DECIDE,ld*bv->m));
    PetscCall(VecSetType(bv->buffer,bv->vtype));
  }
  *buffer = bv->buffer;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetRandomContext - Sets the `PetscRandom` object associated with the `BV`,
   to be used in operations that need random numbers.

   Collective

   Input Parameters:
+  bv   - the basis vectors context
-  rand - the random number generator context

   Level: advanced

.seealso: [](sec:bv), `BVGetRandomContext()`, `BVSetRandom()`, `BVSetRandomNormal()`, `BVSetRandomColumn()`, `BVSetRandomCond()`, `PetscRandomCreate()`
@*/
PetscErrorCode BVSetRandomContext(BV bv,PetscRandom rand)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(rand,PETSC_RANDOM_CLASSID,2);
  PetscCheckSameComm(bv,1,rand,2);
  PetscCall(PetscObjectReference((PetscObject)rand));
  PetscCall(PetscRandomDestroy(&bv->rand));
  bv->rand = rand;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetRandomContext - Gets the `PetscRandom` object associated with the `BV`.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  rand - the random number generator context

   Level: advanced

.seealso: [](sec:bv), `BVSetRandomContext()`, `BVSetRandom()`, `BVSetRandomNormal()`, `BVSetRandomColumn()`, `BVSetRandomCond()`
@*/
PetscErrorCode BVGetRandomContext(BV bv,PetscRandom* rand)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(rand,2);
  if (!bv->rand) {
    PetscCall(PetscRandomCreate(PetscObjectComm((PetscObject)bv),&bv->rand));
    if (bv->cuda) PetscCall(PetscRandomSetType(bv->rand,PETSCCURAND));
    if (bv->sfocalled) PetscCall(PetscRandomSetFromOptions(bv->rand));
    if (bv->rrandom) {
      PetscCall(PetscRandomSetSeed(bv->rand,0x12345678));
      PetscCall(PetscRandomSeed(bv->rand));
    }
  }
  *rand = bv->rand;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetFromOptions - Sets `BV` options from the options database.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Note:
   To see all options, run your program with the `-help` option.

   Level: beginner

.seealso: [](sec:bv), `BVSetOptionsPrefix()`
@*/
PetscErrorCode BVSetFromOptions(BV bv)
{
  char               type[256];
  PetscBool          flg1,flg2,flg3,flg4;
  PetscReal          r;
  BVOrthogType       otype;
  BVOrthogRefineType orefine;
  BVOrthogBlockType  oblock;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscCall(BVRegisterAll());
  PetscObjectOptionsBegin((PetscObject)bv);
    PetscCall(PetscOptionsFList("-bv_type","Basis Vectors type","BVSetType",BVList,(char*)(((PetscObject)bv)->type_name?((PetscObject)bv)->type_name:BVMAT),type,sizeof(type),&flg1));
    if (flg1) PetscCall(BVSetType(bv,type));
    else if (!((PetscObject)bv)->type_name) PetscCall(BVSetType(bv,BVMAT));

    otype = bv->orthog_type;
    PetscCall(PetscOptionsEnum("-bv_orthog_type","Orthogonalization method","BVSetOrthogonalization",BVOrthogTypes,(PetscEnum)otype,(PetscEnum*)&otype,&flg1));
    orefine = bv->orthog_ref;
    PetscCall(PetscOptionsEnum("-bv_orthog_refine","Iterative refinement mode during orthogonalization","BVSetOrthogonalization",BVOrthogRefineTypes,(PetscEnum)orefine,(PetscEnum*)&orefine,&flg2));
    r = bv->orthog_eta;
    PetscCall(PetscOptionsReal("-bv_orthog_eta","Parameter of iterative refinement during orthogonalization","BVSetOrthogonalization",r,&r,&flg3));
    oblock = bv->orthog_block;
    PetscCall(PetscOptionsEnum("-bv_orthog_block","Block orthogonalization method","BVSetOrthogonalization",BVOrthogBlockTypes,(PetscEnum)oblock,(PetscEnum*)&oblock,&flg4));
    if (flg1 || flg2 || flg3 || flg4) PetscCall(BVSetOrthogonalization(bv,otype,orefine,r,oblock));

    PetscCall(PetscOptionsEnum("-bv_matmult","Method for BVMatMult","BVSetMatMultMethod",BVMatMultTypes,(PetscEnum)bv->vmm,(PetscEnum*)&bv->vmm,NULL));

    PetscCall(PetscOptionsReal("-bv_definite_tol","Tolerance for checking a definite inner product","BVSetDefiniteTolerance",r,&r,&flg1));
    if (flg1) PetscCall(BVSetDefiniteTolerance(bv,r));

    /* undocumented option to generate random vectors that are independent of the number of processes */
    PetscCall(PetscOptionsGetBool(NULL,NULL,"-bv_reproducible_random",&bv->rrandom,NULL));

    if (bv->ops->create) bv->defersfo = PETSC_TRUE;   /* defer call to setfromoptions */
    else PetscTryTypeMethod(bv,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)bv,PetscOptionsObject));
  PetscOptionsEnd();
  bv->sfocalled = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetOrthogonalization - Specifies the method used for the orthogonalization of
   vectors (classical or modified Gram-Schmidt with or without refinement), and
   for the block-orthogonalization (simultaneous orthogonalization of a set of
   vectors).

   Logically Collective

   Input Parameters:
+  bv     - the basis vectors context
.  type   - the method of vector orthogonalization
.  refine - type of refinement
.  eta    - parameter for selective refinement
-  block  - the method of block orthogonalization

   Options Database Keys:
+  -bv_orthog_type \<type\>     - the vector orthogonalization, either `cgs` or `mgs`
.  -bv_orthog_refine \<refine\> - the refinement type, `never`, `ifneeded` (default) or `always`
.  -bv_orthog_eta \<eta\>       - the value of `eta`
-  -bv_orthog_block \<block\>   - the block-orthogonalization method

   Notes:
   The default settings work well for most problems.

   The parameter `eta` should be a real value between 0 and 1, that is used only when
   the refinement type is `ifneeded`. Use `PETSC_DETERMINE` to set a reasonable
   default value. Use `PETSC_CURRENT` to leave the current value unchanged.

   When using several processes, `BV_ORTHOG_MGS` is likely to result in bad scalability.

   If the method set for block orthogonalization is `BV_ORTHOG_BLOCK_GS`, then the computation
   is done column by column with the vector orthogonalization.

   Level: advanced

.seealso: [](sec:bv), `BVOrthogonalizeColumn()`, `BVGetOrthogonalization()`, `BVOrthogType`, `BVOrthogRefineType`, `BVOrthogBlockType`
@*/
PetscErrorCode BVSetOrthogonalization(BV bv,BVOrthogType type,BVOrthogRefineType refine,PetscReal eta,BVOrthogBlockType block)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveEnum(bv,type,2);
  PetscValidLogicalCollectiveEnum(bv,refine,3);
  PetscValidLogicalCollectiveReal(bv,eta,4);
  PetscValidLogicalCollectiveEnum(bv,block,5);
  switch (type) {
    case BV_ORTHOG_CGS:
    case BV_ORTHOG_MGS:
      bv->orthog_type = type;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Unknown orthogonalization type");
  }
  switch (refine) {
    case BV_ORTHOG_REFINE_NEVER:
    case BV_ORTHOG_REFINE_IFNEEDED:
    case BV_ORTHOG_REFINE_ALWAYS:
      bv->orthog_ref = refine;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Unknown refinement type");
  }
  if (eta == (PetscReal)PETSC_DETERMINE) {
    bv->orthog_eta = 0.7071;
  } else if (eta != (PetscReal)PETSC_CURRENT) {
    PetscCheck(eta>0.0 && eta<=1.0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Invalid eta value");
    bv->orthog_eta = eta;
  }
  switch (block) {
    case BV_ORTHOG_BLOCK_GS:
    case BV_ORTHOG_BLOCK_CHOL:
    case BV_ORTHOG_BLOCK_TSQR:
    case BV_ORTHOG_BLOCK_TSQRCHOL:
    case BV_ORTHOG_BLOCK_SVQB:
      bv->orthog_block = block;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Unknown block orthogonalization type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetOrthogonalization - Gets the orthogonalization settings from the `BV` object.

   Not Collective

   Input Parameter:
.  bv - basis vectors context

   Output Parameters:
+  type   - the method of vector orthogonalization
.  refine - type of refinement
.  eta    - parameter for selective refinement
-  block  - the method of block orthogonalization

   Level: advanced

.seealso: [](sec:bv), `BVOrthogonalizeColumn()`, `BVSetOrthogonalization()`, `BVOrthogType`, `BVOrthogRefineType`, `BVOrthogBlockType`
@*/
PetscErrorCode BVGetOrthogonalization(BV bv,BVOrthogType *type,BVOrthogRefineType *refine,PetscReal *eta,BVOrthogBlockType *block)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (type)   *type   = bv->orthog_type;
  if (refine) *refine = bv->orthog_ref;
  if (eta)    *eta    = bv->orthog_eta;
  if (block)  *block  = bv->orthog_block;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetMatMultMethod - Specifies the method used for the `BVMatMult()` operation.

   Logically Collective

   Input Parameters:
+  bv     - the basis vectors context
-  method - the method for the `BVMatMult()` operation

   Options Database Key:
.  -bv_matmult \<method\> - choose one of the methods: `vecs`, `mat`

   Note:
   The default is `BV_MATMULT_MAT` except in the case of `BVVECS`.

   Level: advanced

.seealso: [](sec:bv), `BVMatMult()`, `BVGetMatMultMethod()`, `BVMatMultType`
@*/
PetscErrorCode BVSetMatMultMethod(BV bv,BVMatMultType method)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveEnum(bv,method,2);
  switch (method) {
    case BV_MATMULT_VECS:
    case BV_MATMULT_MAT:
      bv->vmm = method;
      break;
    case BV_MATMULT_MAT_SAVE:
      PetscCall(PetscInfo(bv,"BV_MATMULT_MAT_SAVE is deprecated, using BV_MATMULT_MAT\n"));
      bv->vmm = BV_MATMULT_MAT;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Unknown matmult method");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetMatMultMethod - Gets the method used for the `BVMatMult()` operation.

   Not Collective

   Input Parameter:
.  bv - basis vectors context

   Output Parameter:
.  method - the method for the `BVMatMult()` operation

   Level: advanced

.seealso: [](sec:bv), `BVMatMult()`, `BVSetMatMultMethod()`, `BVMatMultType`
@*/
PetscErrorCode BVGetMatMultMethod(BV bv,BVMatMultType *method)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(method,2);
  *method = bv->vmm;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetColumn - Returns a `Vec` object that contains the entries of the
   requested column of the basis vectors object.

   Logically Collective

   Input Parameters:
+  bv - the basis vectors context
-  j  - the index of the requested column

   Output Parameter:
.  v  - vector containing the jth column

   Notes:
   The returned `Vec` must be seen as a reference (not a copy) of the `BV`
   column, that is, modifying the `Vec` will change the `BV` entries as well.

   The returned `Vec` must not be destroyed. `BVRestoreColumn()` must be
   called when it is no longer needed. At most, two columns can be fetched,
   that is, this function can only be called twice before the corresponding
   `BVRestoreColumn()` is invoked.

   A negative index `j` selects the `i`-th constraint, where `i=-j`. Constraints
   should not be modified.

   Level: beginner

.seealso: [](sec:bv), `BVRestoreColumn()`, `BVInsertConstraints()`
@*/
PetscErrorCode BVGetColumn(BV bv,PetscInt j,Vec *v)
{
  PetscInt       l;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  BVCheckOp(bv,1,getcolumn);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscCheck(j>=0 || -j<=bv->nc,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"You requested constraint %" PetscInt_FMT " but only %" PetscInt_FMT " are available",-j,bv->nc);
  PetscCheck(j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"You requested column %" PetscInt_FMT " but only %" PetscInt_FMT " are available",j,bv->m);
  PetscCheck(j!=bv->ci[0] && j!=bv->ci[1],PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Column %" PetscInt_FMT " already fetched in a previous call to BVGetColumn",j);
  l = BVAvailableVec;
  PetscCheck(l!=-1,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Too many requested columns; you must call BVRestoreColumn for one of the previously fetched columns");
  PetscUseTypeMethod(bv,getcolumn,j,v);
  bv->ci[l] = j;
  PetscCall(VecGetState(bv->cv[l],&bv->st[l]));
  PetscCall(PetscObjectGetId((PetscObject)bv->cv[l],&bv->id[l]));
  *v = bv->cv[l];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVRestoreColumn - Restore a column obtained with `BVGetColumn()`.

   Logically Collective

   Input Parameters:
+  bv - the basis vectors context
.  j  - the index of the column
-  v  - vector obtained with `BVGetColumn()`

   Note:
   The arguments must match the corresponding call to `BVGetColumn()`.

   Level: beginner

.seealso: [](sec:bv), `BVGetColumn()`
@*/
PetscErrorCode BVRestoreColumn(BV bv,PetscInt j,Vec *v)
{
  PetscObjectId    id;
  PetscObjectState st;
  PetscInt         l;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscAssertPointer(v,3);
  PetscValidHeaderSpecific(*v,VEC_CLASSID,3);
  PetscCheck(j>=0 || -j<=bv->nc,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"You requested constraint %" PetscInt_FMT " but only %" PetscInt_FMT " are available",-j,bv->nc);
  PetscCheck(j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"You requested column %" PetscInt_FMT " but only %" PetscInt_FMT " are available",j,bv->m);
  PetscCheck(j==bv->ci[0] || j==bv->ci[1],PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Column %" PetscInt_FMT " has not been fetched with a call to BVGetColumn",j);
  l = (j==bv->ci[0])? 0: 1;
  PetscCall(PetscObjectGetId((PetscObject)*v,&id));
  PetscCheck(id==bv->id[l],PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Argument 3 is not the same Vec that was obtained with BVGetColumn");
  PetscCall(VecGetState(*v,&st));
  if (st!=bv->st[l]) PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscUseTypeMethod(bv,restorecolumn,j,v);
  bv->ci[l] = -bv->nc-1;
  bv->st[l] = -1;
  bv->id[l] = 0;
  *v = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   BVGetArray - Returns a pointer to a contiguous array that contains this
   process' portion of the `BV` data.

   Logically Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  a  - location to put pointer to the array

   Notes:
   `BVRestoreArray()` must be called when access to the array is no longer needed.
   This operation may imply a data copy, for `BV` types that do not store
   data contiguously in memory.

   The pointer will normally point to the first entry of the first column,
   but if the `BV` has constraints then these go before the regular columns.

   Note that for manipulating the pointer to the `BV` array, one must take into
   account the leading dimension, which might be different from the local
   number of rows, see `BVGetLeadingDimension()`.

   Use `BVGetArrayRead()` for read-only access.

   Level: advanced

.seealso: [](sec:bv), `BVRestoreArray()`, `BVInsertConstraints()`, `BVGetLeadingDimension()`, `BVGetArrayRead()`, `BVType`
@*/
PetscErrorCode BVGetArray(BV bv,PetscScalar **a)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  BVCheckOp(bv,1,getarray);
  PetscUseTypeMethod(bv,getarray,a);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   BVRestoreArray - Restore the `BV` object after `BVGetArray()` has been called.

   Logically Collective

   Input Parameters:
+  bv - the basis vectors context
-  a  - location of pointer to array obtained from `BVGetArray()`

   Note:
   This operation may imply a data copy, for `BV` types that do not store
   data contiguously in memory.

   Level: advanced

.seealso: [](sec:bv), `BVGetArray()`
@*/
PetscErrorCode BVRestoreArray(BV bv,PetscScalar **a)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscTryTypeMethod(bv,restorearray,a);
  if (a) *a = NULL;
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   BVGetArrayRead - Returns a read-only pointer to a contiguous array that
   contains this processor's portion of the `BV` data.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  a  - location to put pointer to the array

   Notes:
   `BVRestoreArrayRead()` must be called when access to the array is no
   longer needed. This operation may imply a data copy, for `BV` types that
   do not store data contiguously in memory.

   The pointer will normally point to the first entry of the first column,
   but if the `BV` has constraints then these go before the regular columns.

   Level: advanced

.seealso: [](sec:bv), `BVRestoreArrayRead()`, `BVInsertConstraints()`, `BVGetLeadingDimension()`, `BVGetArray()`, `BVType`
@*/
PetscErrorCode BVGetArrayRead(BV bv,const PetscScalar **a)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  BVCheckOp(bv,1,getarrayread);
  PetscUseTypeMethod(bv,getarrayread,a);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   BVRestoreArrayRead - Restore the `BV` object after `BVGetArrayRead()` has
   been called.

   Not Collective

   Input Parameters:
+  bv - the basis vectors context
-  a  - location of pointer to array obtained from `BVGetArrayRead()`

   Level: advanced

.seealso: [](sec:bv), `BVGetArrayRead()`
@*/
PetscErrorCode BVRestoreArrayRead(BV bv,const PetscScalar **a)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscTryTypeMethod(bv,restorearrayread,a);
  if (a) *a = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCreateVec - Creates a new `Vec` object with the same type and dimensions
   as the columns of the basis vectors object.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  v  - the new vector

   Note:
   The user is responsible for destroying the returned vector.

   Level: beginner

.seealso: [](sec:bv), `BVCreateMat()`, `BVCreateVecEmpty()`
@*/
PetscErrorCode BVCreateVec(BV bv,Vec *v)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscAssertPointer(v,2);
  PetscCall(VecCreate(PetscObjectComm((PetscObject)bv),v));
  PetscCall(VecSetLayout(*v,bv->map));
  PetscCall(VecSetType(*v,bv->vtype));
  PetscCall(VecSetUp(*v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCreateVecEmpty - Creates a new `Vec` object with the same type and dimensions
   as the columns of the basis vectors object, but without internal array.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  v  - the new vector

   Note:
   This works as `BVCreateVec()`, but the new vector does not have the array allocated,
   so the intended usage is with `VecPlaceArray()`.

   Level: developer

.seealso: [](sec:bv), `BVCreateVec()`
@*/
PetscErrorCode BVCreateVecEmpty(BV bv,Vec *v)
{
  PetscBool  standard,cuda,hip,mpi;
  PetscInt   N,nloc,bs;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscAssertPointer(v,2);

  PetscCall(PetscStrcmpAny(bv->vtype,&standard,VECSEQ,VECMPI,""));
  PetscCall(PetscStrcmpAny(bv->vtype,&cuda,VECSEQCUDA,VECMPICUDA,""));
  PetscCall(PetscStrcmpAny(bv->vtype,&hip,VECSEQHIP,VECMPIHIP,""));
  if (standard || cuda || hip) {
    PetscCall(PetscStrcmpAny(bv->vtype,&mpi,VECMPI,VECMPICUDA,VECMPIHIP,""));
    PetscCall(PetscLayoutGetLocalSize(bv->map,&nloc));
    PetscCall(PetscLayoutGetSize(bv->map,&N));
    PetscCall(PetscLayoutGetBlockSize(bv->map,&bs));
    if (cuda) {
#if defined(PETSC_HAVE_CUDA)
      if (mpi) PetscCall(VecCreateMPICUDAWithArray(PetscObjectComm((PetscObject)bv),bs,nloc,N,NULL,v));
      else PetscCall(VecCreateSeqCUDAWithArray(PetscObjectComm((PetscObject)bv),bs,N,NULL,v));
#endif
    } else if (hip) {
#if defined(PETSC_HAVE_HIP)
      if (mpi) PetscCall(VecCreateMPIHIPWithArray(PetscObjectComm((PetscObject)bv),bs,nloc,N,NULL,v));
      else PetscCall(VecCreateSeqHIPWithArray(PetscObjectComm((PetscObject)bv),bs,N,NULL,v));
#endif
    } else {
      if (mpi) PetscCall(VecCreateMPIWithArray(PetscObjectComm((PetscObject)bv),bs,nloc,N,NULL,v));
      else PetscCall(VecCreateSeqWithArray(PetscObjectComm((PetscObject)bv),bs,N,NULL,v));
    }
  } else PetscCall(BVCreateVec(bv,v)); /* standard duplicate, with internal array */
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetVecType - Set the vector type to be used when creating vectors via `BVCreateVec()`.

   Collective

   Input Parameters:
+  bv - the basis vectors context
-  vtype - the vector type

   Level: advanced

   Note:
   This is not needed if the `BV` object is set up with `BVSetSizesFromVec()`, but may be
   required in the case of `BVSetSizes()` if one wants to work with non-standard vectors.

.seealso: [](sec:bv), `BVCreateVec()`, `BVGetVecType()`, `BVSetSizesFromVec()`, `BVSetSizes()`
@*/
PetscErrorCode BVSetVecType(BV bv,VecType vtype)
{
  PetscBool   std;
  PetscMPIInt size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscCall(PetscFree(bv->vtype));
  PetscCall(PetscStrcmp(vtype,VECSTANDARD,&std));
  if (std) {
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)bv),&size));
    PetscCall(PetscStrallocpy((size==1)?VECSEQ:VECMPI,(char**)&bv->vtype));
  } else PetscCall(PetscStrallocpy(vtype,(char**)&bv->vtype));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetVecType - Get the vector type to be used when creating vectors via `BVCreateVec()`.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  vtype - the vector type

   Level: advanced

.seealso: [](sec:bv), `BVSetVecType()`, `BVCreateVec()`
@*/
PetscErrorCode BVGetVecType(BV bv,VecType *vtype)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(vtype,2);
  *vtype = bv->vtype;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCreateMat - Creates a new `Mat` object of dense type and copies the contents
   of the `BV` object.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  A  - the new matrix

   Notes:
   The user is responsible for destroying the returned matrix.

   The matrix contains all columns of the `BV`, not just the active columns.

   Level: intermediate

.seealso: [](sec:bv), `BVCreateFromMat()`, `BVCreateVec()`, `BVGetMat()`
@*/
PetscErrorCode BVCreateMat(BV bv,Mat *A)
{
  PetscInt ksave,lsave;
  Mat      B;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscAssertPointer(A,2);

  PetscCall(MatCreateDenseFromVecType(PetscObjectComm((PetscObject)bv),bv->vtype,bv->n,PETSC_DECIDE,bv->N,bv->m,bv->ld,NULL,A));
  lsave = bv->l;
  ksave = bv->k;
  bv->l = 0;
  bv->k = bv->m;
  PetscCall(BVGetMat(bv,&B));
  PetscCall(MatCopy(B,*A,SAME_NONZERO_PATTERN));
  PetscCall(BVRestoreMat(bv,&B));
  bv->l = lsave;
  bv->k = ksave;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVGetMat_Default(BV bv,Mat *A)
{
  PetscScalar *vv,*aa;
  PetscBool   create=PETSC_FALSE;
  PetscInt    m,cols;

  PetscFunctionBegin;
  m = bv->k-bv->l;
  if (!bv->Aget) create=PETSC_TRUE;
  else {
    PetscCall(MatDenseGetArray(bv->Aget,&aa));
    PetscCheck(!aa,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"BVGetMat already called on this BV");
    PetscCall(MatGetSize(bv->Aget,NULL,&cols));
    if (cols!=m) {
      PetscCall(MatDestroy(&bv->Aget));
      create=PETSC_TRUE;
    }
  }
  PetscCall(BVGetArray(bv,&vv));
  if (create) {
    PetscCall(MatCreateDenseFromVecType(PetscObjectComm((PetscObject)bv),bv->vtype,bv->n,PETSC_DECIDE,bv->N,m,bv->ld,vv,&bv->Aget)); /* pass a pointer to avoid allocation of storage */
    PetscCall(MatDenseReplaceArray(bv->Aget,NULL));  /* replace with a null pointer, the value after BVRestoreMat */
  }
  PetscCall(MatDensePlaceArray(bv->Aget,vv+(bv->nc+bv->l)*bv->ld));  /* set the actual pointer */
  *A = bv->Aget;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetMat - Returns a `Mat` object of dense type that shares the memory of
   the `BV` object.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  A  - the matrix

   Notes:
   The returned matrix contains only the active columns. If the content of
   the `Mat` is modified, these changes are also done in the `BV` object. The
   user must call `BVRestoreMat()` when no longer needed.

   This operation implies a call to `BVGetArray()`, which may result in data
   copies.

   Level: advanced

.seealso: [](sec:bv), `BVRestoreMat()`, `BVCreateMat()`, `BVGetArray()`
@*/
PetscErrorCode BVGetMat(BV bv,Mat *A)
{
  char name[64];

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscAssertPointer(A,2);
  PetscUseTypeMethod(bv,getmat,A);
  if (((PetscObject)bv)->name) {   /* set A's name based on BV name */
    PetscCall(PetscStrncpy(name,"Mat_",sizeof(name)));
    PetscCall(PetscStrlcat(name,((PetscObject)bv)->name,sizeof(name)));
    PetscCall(PetscObjectSetName((PetscObject)*A,name));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVRestoreMat_Default(BV bv,Mat *A)
{
  PetscScalar *vv,*aa;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArray(bv->Aget,&aa));
  vv = aa-(bv->nc+bv->l)*bv->ld;
  PetscCall(MatDenseResetArray(bv->Aget));
  PetscCall(BVRestoreArray(bv,&vv));
  *A = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVRestoreMat - Restores the Mat obtained with `BVGetMat()`.

   Logically Collective

   Input Parameters:
+  bv - the basis vectors context
-  A  - the fetched matrix

   Note:
   A call to this function must match a previous call of `BVGetMat()`.
   The effect is that the contents of the `Mat` are copied back to the
   `BV` internal data structures.

   Level: advanced

.seealso: [](sec:bv), `BVGetMat()`, `BVRestoreArray()`
@*/
PetscErrorCode BVRestoreMat(BV bv,Mat *A)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscAssertPointer(A,2);
  PetscCheck(bv->Aget,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"BVRestoreMat must match a previous call to BVGetMat");
  PetscCheck(bv->Aget==*A,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Mat argument is not the same as the one obtained with BVGetMat");
  PetscUseTypeMethod(bv,restoremat,A);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Copy all user-provided attributes of V to another BV object W
 */
static inline PetscErrorCode BVDuplicate_Private(BV V,BV W)
{
  PetscFunctionBegin;
  PetscCall(PetscLayoutReference(V->map,&W->map));
  PetscCall(BVSetVecType(W,V->vtype));
  W->ld           = V->ld;
  PetscCall(BVSetType(W,((PetscObject)V)->type_name));
  W->orthog_type  = V->orthog_type;
  W->orthog_ref   = V->orthog_ref;
  W->orthog_eta   = V->orthog_eta;
  W->orthog_block = V->orthog_block;
  if (V->matrix) PetscCall(PetscObjectReference((PetscObject)V->matrix));
  W->matrix       = V->matrix;
  W->indef        = V->indef;
  W->vmm          = V->vmm;
  W->rrandom      = V->rrandom;
  W->deftol       = V->deftol;
  if (V->rand) PetscCall(PetscObjectReference((PetscObject)V->rand));
  W->rand         = V->rand;
  W->sfocalled    = V->sfocalled;
  PetscTryTypeMethod(V,duplicate,W);
  PetscCall(PetscObjectStateIncrease((PetscObject)W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDuplicate - Creates a new basis vectors object of the same type and
   dimensions as an existing one.

   Collective

   Input Parameter:
.  V - basis vectors context

   Output Parameter:
.  W - location to put the new `BV`

   Notes:
   The new `BV` has the same type and dimensions as `V`. Also, the inner
   product matrix and orthogonalization options are copied.

   `BVDuplicate()` DOES NOT COPY the entries, but rather allocates storage
   for the new basis vectors. Use `BVCopy()` to copy the content.

   Level: intermediate

.seealso: [](sec:bv), `BVDuplicateResize()`, `BVCreate()`, `BVSetSizesFromVec()`, `BVCopy()`
@*/
PetscErrorCode BVDuplicate(BV V,BV *W)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscAssertPointer(W,2);
  PetscCall(BVCreate(PetscObjectComm((PetscObject)V),W));
  (*W)->N = V->N;
  (*W)->n = V->n;
  (*W)->m = V->m;
  (*W)->k = V->m;
  PetscCall(BVDuplicate_Private(V,*W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDuplicateResize - Creates a new basis vectors object of the same type and
   dimensions as an existing one, but with possibly different number of columns.

   Collective

   Input Parameters:
+  V - basis vectors context
-  m - the new number of columns

   Output Parameter:
.  W - location to put the new BV

   Note:
   This is equivalent to a call to `BVDuplicate()` followed by `BVResize()`. The
   contents of `V` are not copied to `W`.

   Level: intermediate

.seealso: [](sec:bv), `BVDuplicate()`, `BVResize()`
@*/
PetscErrorCode BVDuplicateResize(BV V,PetscInt m,BV *W)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidLogicalCollectiveInt(V,m,2);
  PetscAssertPointer(W,3);
  PetscCall(BVCreate(PetscObjectComm((PetscObject)V),W));
  (*W)->N = V->N;
  (*W)->n = V->n;
  (*W)->m = m;
  (*W)->k = m;
  PetscCall(BVDuplicate_Private(V,*W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetCachedBV - Returns a `BV` object stored internally that holds the
   result of $BX$ after a call to `BVApplyMatrixBV()`.

   Collective

   Input Parameter:
.  bv    - the basis vectors context

   Output Parameter:
.  cached - the cached `BV`

   Note:
   The cached `BV` is created if not available previously.

   Level: developer

.seealso: [](sec:bv), `BVApplyMatrixBV()`
@*/
PetscErrorCode BVGetCachedBV(BV bv,BV *cached)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(cached,2);
  BVCheckSizes(bv,1);
  if (!bv->cached) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)bv),&bv->cached));
    bv->cached->N = bv->N;
    bv->cached->n = bv->n;
    bv->cached->m = bv->m;
    bv->cached->k = bv->m;
    PetscCall(BVDuplicate_Private(bv,bv->cached));
  }
  *cached = bv->cached;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCopy - Copies a basis vector object into another one, $W \leftarrow V$.

   Logically Collective

   Input Parameter:
.  V - basis vectors context

   Output Parameter:
.  W - the copy

   Note:
   Both `V` and `W` must be distributed in the same manner; local copies are
   done. Only active columns (excluding the leading ones) are copied.
   In the destination `W`, columns are overwritten starting from the leading ones.
   Constraints are not copied.

   Level: beginner

.seealso: [](sec:bv), `BVCopyVec()`, `BVCopyColumn()`, `BVDuplicate()`, `BVSetActiveColumns()`
@*/
PetscErrorCode BVCopy(BV V,BV W)
{
  PetscScalar       *womega;
  const PetscScalar *vomega;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  BVCheckOp(V,1,copy);
  PetscValidHeaderSpecific(W,BV_CLASSID,2);
  PetscValidType(W,2);
  BVCheckSizes(W,2);
  PetscCheckSameTypeAndComm(V,1,W,2);
  PetscCheck(V->n==W->n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension V %" PetscInt_FMT ", W %" PetscInt_FMT,V->n,W->n);
  PetscCheck(V->k-V->l<=W->m-W->l,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"W has %" PetscInt_FMT " non-leading columns, not enough to store %" PetscInt_FMT " columns",W->m-W->l,V->k-V->l);
  if (V==W || !V->n) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscLogEventBegin(BV_Copy,V,W,0,0));
  if (V->indef && V->matrix && V->indef==W->indef && V->matrix==W->matrix) {
    /* copy signature */
    PetscCall(BV_AllocateSignature(W));
    PetscCall(VecGetArrayRead(V->omega,&vomega));
    PetscCall(VecGetArray(W->omega,&womega));
    PetscCall(PetscArraycpy(womega+W->nc+W->l,vomega+V->nc+V->l,V->k-V->l));
    PetscCall(VecRestoreArray(W->omega,&womega));
    PetscCall(VecRestoreArrayRead(V->omega,&vomega));
  }
  PetscUseTypeMethod(V,copy,W);
  PetscCall(PetscLogEventEnd(BV_Copy,V,W,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCopyVec - Copies one of the columns of a basis vectors object into a `Vec`.

   Logically Collective

   Input Parameters:
+  V - basis vectors context
-  j - the index of the column to be copied

   Output Parameter:
.  w - the copied column

   Note:
   Both `V` and `w` must be distributed in the same manner; local copies are done.

   Level: beginner

.seealso: [](sec:bv), `BVCopy()`, `BVCopyColumn()`, `BVInsertVec()`
@*/
PetscErrorCode BVCopyVec(BV V,PetscInt j,Vec w)
{
  PetscInt       n,N;
  Vec            z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidLogicalCollectiveInt(V,j,2);
  PetscValidHeaderSpecific(w,VEC_CLASSID,3);
  PetscCheckSameComm(V,1,w,3);

  PetscCall(VecGetSize(w,&N));
  PetscCall(VecGetLocalSize(w,&n));
  PetscCheck(N==V->N && n==V->n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Vec sizes (global %" PetscInt_FMT ", local %" PetscInt_FMT ") do not match BV sizes (global %" PetscInt_FMT ", local %" PetscInt_FMT ")",N,n,V->N,V->n);

  PetscCall(PetscLogEventBegin(BV_Copy,V,w,0,0));
  PetscCall(BVGetColumn(V,j,&z));
  PetscCall(VecCopy(z,w));
  PetscCall(BVRestoreColumn(V,j,&z));
  PetscCall(PetscLogEventEnd(BV_Copy,V,w,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCopyColumn - Copies the values from one of the columns to another one.

   Logically Collective

   Input Parameters:
+  V - basis vectors context
.  j - the index of the source column
-  i - the index of the destination column

   Level: beginner

.seealso: [](sec:bv), `BVCopy()`, `BVCopyVec()`
@*/
PetscErrorCode BVCopyColumn(BV V,PetscInt j,PetscInt i)
{
  PetscScalar *omega;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidLogicalCollectiveInt(V,j,2);
  PetscValidLogicalCollectiveInt(V,i,3);
  if (j==i) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscLogEventBegin(BV_Copy,V,0,0,0));
  if (V->omega) {
    PetscCall(VecGetArray(V->omega,&omega));
    omega[i] = omega[j];
    PetscCall(VecRestoreArray(V->omega,&omega));
  }
  PetscUseTypeMethod(V,copycolumn,j,i);
  PetscCall(PetscLogEventEnd(BV_Copy,V,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetSplit_Private(BV bv,PetscBool left,BV *split)
{
  PetscInt  ncols;

  PetscFunctionBegin;
  ncols = left? bv->nc+bv->l: bv->m-bv->l;
  if (*split && (ncols!=(*split)->m || bv->N!=(*split)->N)) PetscCall(BVDestroy(split));
  if (!*split) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)bv),split));
    (*split)->issplit = left? 1: 2;
    (*split)->splitparent = bv;
    (*split)->N = bv->N;
    (*split)->n = bv->n;
    (*split)->m = bv->m;
    (*split)->k = bv->m;
    PetscCall(BVDuplicate_Private(bv,*split));
  }
  (*split)->l  = 0;
  (*split)->k  = left? bv->l: bv->k-bv->l;
  (*split)->nc = left? bv->nc: 0;
  (*split)->m  = ncols-(*split)->nc;
  if ((*split)->nc) {
    (*split)->ci[0] = -(*split)->nc-1;
    (*split)->ci[1] = -(*split)->nc-1;
  }
  if (left) PetscCall(PetscObjectStateGet((PetscObject)*split,&bv->lstate));
  else PetscCall(PetscObjectStateGet((PetscObject)*split,&bv->rstate));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetSplit - Splits the `BV` object into two `BV` objects that share the
   internal data, one of them containing the leading columns and the other
   one containing the remaining columns.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameters:
+  L - left `BV` containing leading columns (can be `NULL`)
-  R - right `BV` containing remaining columns (can be `NULL`)

   Notes:
   The columns are split in two sets. The leading columns (including the
   constraints) are assigned to the left `BV` and the remaining columns
   are assigned to the right `BV`. The number of leading columns, as
   specified with `BVSetActiveColumns()`, must be between `1` and `m-1` (to
   guarantee that both `L` and `R` have at least one column).

   The returned `BV`s must be seen as references (not copies) of the input
   `BV`, that is, modifying them will change the entries of `bv` as well.
   The returned `BV`s must not be destroyed. `BVRestoreSplit()` must be called
   when they are no longer needed.

   Pass `NULL` for any of the output `BV`s that is not needed.

   Level: advanced

.seealso: [](sec:bv), `BVRestoreSplit()`, `BVSetActiveColumns()`, `BVSetNumConstraints()`, `BVGetSplitRows()`
@*/
PetscErrorCode BVGetSplit(BV bv,BV *L,BV *R)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscCheck(bv->l,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Must indicate the number of leading columns with BVSetActiveColumns()");
  PetscCheck(bv->lsplit>=0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Cannot call BVGetSplit() after BVGetSplitRows()");
  PetscCheck(!bv->lsplit,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Cannot get the split BV's twice before restoring them with BVRestoreSplit()");
  bv->lsplit = bv->nc+bv->l;
  PetscCall(BVGetSplit_Private(bv,PETSC_TRUE,&bv->L));
  PetscCall(BVGetSplit_Private(bv,PETSC_FALSE,&bv->R));
  if (L) *L = bv->L;
  if (R) *R = bv->R;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVRestoreSplit - Restore the BV objects obtained with `BVGetSplit()`.

   Logically Collective

   Input Parameters:
+  bv - the basis vectors context
.  L  - left BV obtained with BVGetSplit()
-  R  - right BV obtained with BVGetSplit()

   Note:
   The arguments must match the corresponding call to `BVGetSplit()`.

   Level: advanced

.seealso: [](sec:bv), `BVGetSplit()`
@*/
PetscErrorCode BVRestoreSplit(BV bv,BV *L,BV *R)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscCheck(bv->lsplit>0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Must call BVGetSplit first");
  PetscCheck(!L || *L==bv->L,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Argument 2 is not the same BV that was obtained with BVGetSplit");
  PetscCheck(!R || *R==bv->R,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Argument 3 is not the same BV that was obtained with BVGetSplit");
  PetscCheck(!L || ((*L)->ci[0]<=(*L)->nc-1 && (*L)->ci[1]<=(*L)->nc-1),PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Argument 2 has unrestored columns, use BVRestoreColumn()");
  PetscCheck(!R || ((*R)->ci[0]<=(*R)->nc-1 && (*R)->ci[1]<=(*R)->nc-1),PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Argument 3 has unrestored columns, use BVRestoreColumn()");

  PetscTryTypeMethod(bv,restoresplit,L,R);
  bv->lsplit = 0;
  if (L) *L = NULL;
  if (R) *R = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Copy all user-provided attributes of V to another BV object W with different layout
 */
static inline PetscErrorCode BVDuplicateNewLayout_Private(BV V,BV W)
{
  PetscFunctionBegin;
  PetscCall(PetscLayoutCreateFromSizes(PetscObjectComm((PetscObject)V),W->n,W->N,1,&W->map));
  PetscCall(BVSetVecType(W,V->vtype));
  W->ld           = V->ld;
  PetscCall(BVSetType(W,((PetscObject)V)->type_name));
  W->orthog_type  = V->orthog_type;
  W->orthog_ref   = V->orthog_ref;
  W->orthog_eta   = V->orthog_eta;
  W->orthog_block = V->orthog_block;
  W->vmm          = V->vmm;
  W->rrandom      = V->rrandom;
  W->deftol       = V->deftol;
  if (V->rand) PetscCall(PetscObjectReference((PetscObject)V->rand));
  W->rand         = V->rand;
  W->sfocalled    = V->sfocalled;
  PetscTryTypeMethod(V,duplicate,W);
  PetscCall(PetscObjectStateIncrease((PetscObject)W));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVGetSplitRows_Private(BV bv,PetscBool top,IS is,BV *split)
{
  PetscInt  rstart,rend,lstart;
  PetscInt  N,n;
  PetscBool contig;

  PetscFunctionBegin;
  PetscCall(PetscLayoutGetRange(bv->map,&rstart,&rend));
  PetscCall(ISContiguousLocal(is,rstart,rend,&lstart,&contig));
  PetscCheck(contig,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"%s index set is not contiguous",(top==PETSC_TRUE)?"Upper":"Lower");
  if (top) PetscCheck(lstart==0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Upper index set should start at first local row");
  else bv->lsplit = -lstart;
  PetscCall(ISGetSize(is,&N));
  PetscCall(ISGetLocalSize(is,&n));
  if (*split && (bv->m!=(*split)->m || N!=(*split)->N)) PetscCall(BVDestroy(split));
  if (!*split) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)bv),split));
    (*split)->issplit = top? -1: -2;
    (*split)->splitparent = bv;
    (*split)->N = N;
    (*split)->n = n;
    (*split)->m = bv->m;
    PetscCall(BVDuplicateNewLayout_Private(bv,*split));
  }
  (*split)->k  = bv->k;
  (*split)->l  = bv->l;
  (*split)->nc = bv->nc;
  if ((*split)->nc) {
    (*split)->ci[0] = -(*split)->nc-1;
    (*split)->ci[1] = -(*split)->nc-1;
  }
  if (top) PetscCall(PetscObjectStateGet((PetscObject)*split,&bv->rstate));
  else PetscCall(PetscObjectStateGet((PetscObject)*split,&bv->lstate));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetSplitRows - Splits the `BV` object into two `BV` objects that share the
   internal data, using a disjoint horizontal splitting.

   Collective

   Input Parameters:
+  bv   - the basis vectors context
.  isup - the index set that defines the upper part of the horizontal splitting
-  islo - the index set that defines the lower part of the horizontal splitting

   Output Parameters:
+  U - the resulting `BV` containing the upper rows
-  L - the resulting `BV` containing the lower rows

   Notes:
   The index sets must be such that every MPI process can extract the selected
   rows from its local part of the input `BV`, and this part must be contiguous.
   With one process, `isup` will list contiguous indices starting from 0, and `islo`
   will contain the remaining indices, hence we refer to upper and lower part.
   However, with several processes the indices will be interleaved because
   `isup` will refer to the upper part of the local array.

   The intended use of this function is with matrices of `MATNEST` type, where
   `MatNestGetISs()` will return the appropriate index sets.

   The returned `BV`s must be seen as references (not copies) of the input
   `BV`, that is, modifying them will change the entries of `bv` as well.
   The returned `BV`s must not be destroyed. `BVRestoreSplitRows()` must be called
   when they are no longer needed.

   Pass `NULL` for any of the output `BV`s that is not needed.

   Level: advanced

.seealso: [](sec:bv), `BVRestoreSplitRows()`, `BVGetSplit()`
@*/
PetscErrorCode BVGetSplitRows(BV bv,IS isup,IS islo,BV *U,BV *L)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  if (U) PetscValidHeaderSpecific(isup,IS_CLASSID,2);
  if (L) PetscValidHeaderSpecific(islo,IS_CLASSID,3);
  PetscCheck(bv->lsplit<=0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Cannot call BVGetSplitRows() after BVGetSplit()");
  PetscCheck(!bv->lsplit,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Cannot get the split BV's twice before restoring them with BVRestoreSplitRows()");
  if (U) {
    PetscCall(BVGetSplitRows_Private(bv,PETSC_TRUE,isup,&bv->R));
    *U = bv->R;
  }
  if (L) {
    PetscCall(BVGetSplitRows_Private(bv,PETSC_FALSE,islo,&bv->L));
    *L = bv->L;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVRestoreSplitRows - Restore the BV objects obtained with `BVGetSplitRows()`.

   Logically Collective

   Input Parameters:
+  bv - the basis vectors context
.  isup - the index set that defines the upper part of the horizontal splitting
.  islo - the index set that defines the lower part of the horizontal splitting
.  U - the `BV` containing the upper rows
-  L - the `BV` containing the lower rows

   Note:
   The arguments must match the corresponding call to `BVGetSplitRows()`.

   Level: advanced

.seealso: [](sec:bv), `BVGetSplitRows()`
@*/
PetscErrorCode BVRestoreSplitRows(BV bv,IS isup,IS islo,BV *U,BV *L)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscCheck(bv->lsplit<0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Must call BVGetSplitRows first");
  PetscCheck(!U || ((*U)->ci[0]<=(*U)->nc-1 && (*U)->ci[1]<=(*U)->nc-1),PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"The upper BV has unrestored columns, use BVRestoreColumn()");
  PetscCheck(!L || ((*L)->ci[0]<=(*L)->nc-1 && (*L)->ci[1]<=(*L)->nc-1),PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"The lower BV has unrestored columns, use BVRestoreColumn()");
  PetscTryTypeMethod(bv,restoresplitrows,isup,islo,U,L);
  bv->lsplit = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetDefiniteTolerance - Set the tolerance to be used when checking a
   definite inner product.

   Logically Collective

   Input Parameters:
+  bv     - basis vectors
-  deftol - the tolerance

   Options Database Key:
.  -bv_definite_tol \<deftol\> - the tolerance

   Notes:
   When using a non-standard inner product, see `BVSetMatrix()`, the solver needs
   to compute $\sqrt{z^*B z}$ for various vectors $z$. If the inner product has not
   been declared indefinite, the value $z^*B z$ must be positive, but due to
   rounding error a tiny value may become negative. A tolerance is used to
   detect this situation. Likewise, in complex arithmetic $z^*B z$ should be
   real, and we use the same tolerance to check whether a nonzero imaginary part
   can be considered negligible.

   This function sets this tolerance, which defaults to `10*PETSC_MACHINE_EPSILON`.
   The default value should be good for most applications.

   Level: advanced

.seealso: [](sec:bv), `BVSetMatrix()`
@*/
PetscErrorCode BVSetDefiniteTolerance(BV bv,PetscReal deftol)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveReal(bv,deftol,2);
  if (deftol == (PetscReal)PETSC_DEFAULT || deftol == (PetscReal)PETSC_DECIDE) bv->deftol = 10*PETSC_MACHINE_EPSILON;
  else {
    PetscCheck(deftol>0.0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of deftol. Must be > 0");
    bv->deftol = deftol;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetDefiniteTolerance - Returns the tolerance for checking a definite
   inner product.

   Not Collective

   Input Parameter:
.  bv - the basis vectors

   Output Parameter:
.  deftol - the tolerance

   Level: advanced

.seealso: [](sec:bv), `BVSetDefiniteTolerance()`
@*/
PetscErrorCode BVGetDefiniteTolerance(BV bv,PetscReal *deftol)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(deftol,2);
  *deftol = bv->deftol;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetLeadingDimension - Set the leading dimension to be used for storing the `BV` data.

   Not Collective

   Input Parameters:
+  bv - basis vectors
-  ld - the leading dimension

   Notes:
   This parameter is relevant for `BVMAT`, though it might be employed in other types
   as well.

   When the internal data of the `BV` is stored as a dense matrix, the leading dimension
   has the same meaning as in `MatDenseSetLDA()`, i.e., the distance in number of
   elements from one entry of the matrix to the one in the next column at the same
   row. The leading dimension refers to the local array, and hence can be different
   in different processes.

   The user does not need to change this parameter. The default value is equal to the
   number of local rows, but this value may be increased a little to guarantee alignment
   (especially in the case of GPU storage).

   Level: advanced

.seealso: [](sec:bv), `BVGetLeadingDimension()`
@*/
PetscErrorCode BVSetLeadingDimension(BV bv,PetscInt ld)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,ld,2);
  PetscCheck((bv->n<0 && bv->N<0) || !bv->ops->create,PetscObjectComm((PetscObject)bv),PETSC_ERR_ORDER,"Must call BVSetLeadingDimension() before setting the BV type and sizes");
  if (ld == PETSC_DEFAULT || ld == PETSC_DECIDE) bv->ld = 0;
  else {
    PetscCheck(ld>0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ld. Must be > 0");
    bv->ld = ld;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetLeadingDimension - Returns the leading dimension of the `BV`.

   Not Collective

   Input Parameter:
.  bv - the basis vectors

   Output Parameter:
.  ld - the leading dimension

   Level: advanced

   Notes:
   The returned value may be different in different processes.

   The leading dimension must be used when accessing the internal array via
   `BVGetArray()` or `BVGetArrayRead()`.

.seealso: [](sec:bv), `BVSetLeadingDimension()`, `BVGetArray()`, `BVGetArrayRead()`
@*/
PetscErrorCode BVGetLeadingDimension(BV bv,PetscInt *ld)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(ld,2);
  *ld = bv->ld;
  PetscFunctionReturn(PETSC_SUCCESS);
}
