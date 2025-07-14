/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   MFN routines related to problem setup
*/

#include <slepc/private/mfnimpl.h>       /*I "slepcmfn.h" I*/

/*@
   MFNSetUp - Sets up all the internal data structures necessary for the
   execution of the matrix function solver.

   Collective

   Input Parameter:
.  mfn   - matrix function context

   Notes:
   This function need not be called explicitly in most cases, since MFNSolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: MFNCreate(), MFNSolve(), MFNDestroy()
@*/
PetscErrorCode MFNSetUp(MFN mfn)
{
  PetscInt       N;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);

  /* reset the convergence flag from the previous solves */
  mfn->reason = MFN_CONVERGED_ITERATING;

  if (mfn->setupcalled) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscLogEventBegin(MFN_SetUp,mfn,0,0,0));

  /* Set default solver type (MFNSetFromOptions was not called) */
  if (!((PetscObject)mfn)->type_name) PetscCall(MFNSetType(mfn,MFNKRYLOV));
  if (!mfn->fn) PetscCall(MFNGetFN(mfn,&mfn->fn));
  if (!((PetscObject)mfn->fn)->type_name) PetscCall(FNSetFromOptions(mfn->fn));

  /* Check problem dimensions */
  PetscCheck(mfn->A,PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_WRONGSTATE,"MFNSetOperator must be called first");
  PetscCall(MatGetSize(mfn->A,&N,NULL));
  if (mfn->ncv > N) mfn->ncv = N;

  /* call specific solver setup */
  PetscUseTypeMethod(mfn,setup);

  /* set tolerance if not yet set */
  if (mfn->tol==(PetscReal)PETSC_DETERMINE) mfn->tol = SLEPC_DEFAULT_TOL;

  PetscCall(PetscLogEventEnd(MFN_SetUp,mfn,0,0,0));
  mfn->setupcalled = 1;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNSetOperator - Sets the matrix for which the matrix function is to be computed.

   Collective

   Input Parameters:
+  mfn - the matrix function context
-  A   - the problem matrix

   Notes:
   It must be called before MFNSetUp(). If it is called again after MFNSetUp() then
   the MFN object is reset.

   Level: beginner

.seealso: MFNSolve(), MFNSetUp(), MFNReset()
@*/
PetscErrorCode MFNSetOperator(MFN mfn,Mat A)
{
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(mfn,1,A,2);

  PetscCall(MatGetSize(A,&m,&n));
  PetscCheck(m==n,PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_WRONG,"A is a non-square matrix");
  PetscCall(PetscObjectReference((PetscObject)A));
  if (mfn->setupcalled) PetscCall(MFNReset(mfn));
  else PetscCall(MatDestroy(&mfn->A));
  mfn->A = A;
  mfn->setupcalled = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNGetOperator - Gets the matrix associated with the MFN object.

   Collective

   Input Parameter:
.  mfn - the MFN context

   Output Parameters:
.  A  - the matrix for which the matrix function is to be computed

   Level: intermediate

.seealso: MFNSolve(), MFNSetOperator()
@*/
PetscErrorCode MFNGetOperator(MFN mfn,Mat *A)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscAssertPointer(A,2);
  *A = mfn->A;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MFNAllocateSolution - Allocate memory storage for common variables such
   as the basis vectors.

   Collective

   Input Parameters:
+  mfn   - matrix function context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developer Notes:
   This is SLEPC_EXTERN because it may be required by user plugin MFN
   implementations.

   Level: developer

.seealso: MFNSetUp()
@*/
PetscErrorCode MFNAllocateSolution(MFN mfn,PetscInt extra)
{
  PetscInt       oldsize,requested;
  Vec            t;

  PetscFunctionBegin;
  requested = mfn->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  PetscCall(BVGetSizes(mfn->V,NULL,NULL,&oldsize));

  /* allocate basis vectors */
  if (!mfn->V) PetscCall(MFNGetBV(mfn,&mfn->V));
  if (!oldsize) {
    if (!((PetscObject)mfn->V)->type_name) PetscCall(BVSetType(mfn->V,BVMAT));
    PetscCall(MatCreateVecsEmpty(mfn->A,&t,NULL));
    PetscCall(BVSetSizesFromVec(mfn->V,t,requested));
    PetscCall(VecDestroy(&t));
  } else PetscCall(BVResize(mfn->V,requested,PETSC_FALSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}
