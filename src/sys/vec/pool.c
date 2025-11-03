/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Implementation of a pool of Vec using VecDuplicateVecs
*/

#include <slepc/private/vecimplslepc.h>       /*I "slepcvec.h" I*/

/*
   SlepcVecPoolCreate - Create a pool of Vec.

   Collective

   Input Parameters:
+  v - template vector.
-  init_size - first guess of maximum vectors.

   Output Parameter:
.  p - the pool context.

   Level: developer

.seealso: `SlepcVecPoolGetVecs()`, `SlepcVecPoolDestroy()`
*/
PetscErrorCode SlepcVecPoolCreate(Vec v,PetscInt init_size,VecPool *p)
{
  VecPool_       *pool;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  PetscValidLogicalCollectiveInt(v,init_size,2);
  PetscAssertPointer(p,3);
  PetscCheck(init_size>=0,PetscObjectComm((PetscObject)v),PETSC_ERR_ARG_WRONG,"init_size should be positive");
  PetscCall(PetscCalloc1(1,&pool));
  PetscCall(PetscObjectReference((PetscObject)v));
  pool->v     = v;
  pool->guess = init_size;
  *p = pool;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcVecPoolDestroy - Destroy the pool of Vec.

   Collective

   Input Parameter:
.  p - pool of Vec.

   Level: developer

.seealso: `SlepcVecPoolGetVecs()`, `SlepcVecPoolCreate()`
*/
PetscErrorCode SlepcVecPoolDestroy(VecPool *p)
{
  VecPool_       *pool = (VecPool_*)*p;

  PetscFunctionBegin;
  if (!*p) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(VecDestroy(&pool->v));
  PetscCall(VecDestroyVecs(pool->n,&pool->vecs));
  pool->n     = 0;
  pool->used  = 0;
  pool->guess = 0;
  PetscCall(SlepcVecPoolDestroy((VecPool*)&pool->next));
  PetscCall(PetscFree(pool));
  *p = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcVecPoolGetVecs - Get an array of Vec from the pool.

   Collective

   Input Parameters:
+  p - pool of Vec.
-  n - number of vectors.

   Output Parameter:
.  vecs - vectors

   Level: developer

.seealso: `SlepcVecPoolRestoreVecs()`
*/
PetscErrorCode SlepcVecPoolGetVecs(VecPool p,PetscInt n,Vec **vecs)
{
  VecPool_       *pool = (VecPool_*)p;

  PetscFunctionBegin;
  PetscAssertPointer(p,1);
  PetscAssertPointer(vecs,3);
  PetscCheck(n>=0,PetscObjectComm((PetscObject)pool->v),PETSC_ERR_ARG_OUTOFRANGE,"n should be positive");
  while (pool->next) pool = pool->next;
  if (pool->n-pool->used < n) {
    pool->guess = PetscMax(p->guess,pool->used+n);
    if (pool->vecs && pool->used == 0) PetscCall(VecDestroyVecs(pool->n,&pool->vecs));
    if (pool->vecs) {
      PetscCall(SlepcVecPoolCreate(p->v,pool->guess-pool->used,&pool->next));
      pool = pool->next;
    }
    pool->n = pool->guess;
    PetscCall(VecDuplicateVecs(p->v,pool->n,&pool->vecs));
  }
  *vecs = pool->vecs + pool->used;
  pool->used += n;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcVecPoolRestoreVecs - Get back an array of Vec previously returned by
   SlepcVecPoolGetVecs().

   Collective

   Input Parameters:
+  p    - pool of Vec.
.  n    - number of vectors.
-  vecs - vectors

   Level: developer

.seealso: `SlepcVecPoolGetVecs()`
*/
PetscErrorCode SlepcVecPoolRestoreVecs(VecPool p,PetscInt n,Vec **vecs)
{
  VecPool_       *pool = (VecPool_*)p, *pool0 = pool;

  PetscFunctionBegin;
  while (pool->next) pool = (pool0 = pool)->next;
  if (pool->used == 0 && pool0 != pool) {
    pool0->guess = pool0->used + pool->guess;
    PetscCall(SlepcVecPoolDestroy((VecPool*)&pool));
    pool = pool0;
    pool->next = NULL;
  }
  pool->used -= n;
  PetscCheck(pool->used>=0,PetscObjectComm((PetscObject)pool->v),PETSC_ERR_ARG_OUTOFRANGE,"Unmatched SlepcVecPoolRestoreVecs");
  PetscFunctionReturn(PETSC_SUCCESS);
}
