/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   LME routines related to problem setup
*/

#include <slepc/private/lmeimpl.h>       /*I "slepclme.h" I*/

static inline PetscErrorCode LMESetUp_Lyapunov(LME lme)
{
  Mat            C1,C2,X1,X2;
  Vec            dc,dx;

  PetscFunctionBegin;
  PetscCall(MatLRCGetMats(lme->C,NULL,&C1,&dc,&C2));
  PetscCheck(C1==C2,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"Lyapunov matrix equation requires symmetric right-hand side C");
  PetscCheck(!dc,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"Lyapunov solvers currently require positive-definite right-hand side C");
  if (lme->X) {
    PetscCall(MatLRCGetMats(lme->X,NULL,&X1,&dx,&X2));
    PetscCheck(X1==X2,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"Lyapunov matrix equation requires symmetric solution X");
    PetscCheck(!dx,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"Lyapunov solvers currently assume a positive-definite solution X");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetUp - Sets up all the internal data structures necessary for the
   execution of the linear matrix equation solver.

   Collective

   Input Parameter:
.  lme   - linear matrix equation solver context

   Notes:
   This function need not be called explicitly in most cases, since LMESolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: LMECreate(), LMESolve(), LMEDestroy()
@*/
PetscErrorCode LMESetUp(LME lme)
{
  PetscInt       N;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);

  /* reset the convergence flag from the previous solves */
  lme->reason = LME_CONVERGED_ITERATING;

  if (lme->setupcalled) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscLogEventBegin(LME_SetUp,lme,0,0,0));

  /* Set default solver type (LMESetFromOptions was not called) */
  if (!((PetscObject)lme)->type_name) PetscCall(LMESetType(lme,LMEKRYLOV));

  /* Check problem dimensions */
  PetscCheck(lme->A,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"LMESetCoefficients must be called first");
  PetscCheck(lme->C,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"LMESetRHS must be called first");
  PetscCall(MatGetSize(lme->A,&N,NULL));
  if (lme->ncv > N) lme->ncv = N;

  /* setup options for the particular equation type */
  switch (lme->problem_type) {
    case LME_LYAPUNOV:
      PetscCall(LMESetUp_Lyapunov(lme));
      break;
    case LME_SYLVESTER:
      LMECheckCoeff(lme,lme->B,"B","Sylvester");
      break;
    case LME_GEN_LYAPUNOV:
      LMECheckCoeff(lme,lme->D,"D","Generalized Lyapunov");
      break;
    case LME_GEN_SYLVESTER:
      LMECheckCoeff(lme,lme->B,"B","Generalized Sylvester");
      LMECheckCoeff(lme,lme->D,"D","Generalized Sylvester");
      LMECheckCoeff(lme,lme->E,"E","Generalized Sylvester");
      break;
    case LME_DT_LYAPUNOV:
      break;
    case LME_STEIN:
      LMECheckCoeff(lme,lme->D,"D","Stein");
      break;
  }
  PetscCheck(lme->problem_type==LME_LYAPUNOV,PetscObjectComm((PetscObject)lme),PETSC_ERR_SUP,"There is no solver yet for this matrix equation type");

  /* call specific solver setup */
  PetscUseTypeMethod(lme,setup);

  /* set tolerance if not yet set */
  if (lme->tol==(PetscReal)PETSC_DETERMINE) lme->tol = SLEPC_DEFAULT_TOL;

  PetscCall(PetscLogEventEnd(LME_SetUp,lme,0,0,0));
  lme->setupcalled = 1;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode LMESetCoefficients_Private(LME lme,Mat A,Mat *lmeA)
{
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,&m,&n));
  PetscCheck(m==n,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONG,"Matrix is non-square");
  if (!lme->setupcalled) PetscCall(MatDestroy(lmeA));
  PetscCall(PetscObjectReference((PetscObject)A));
  *lmeA = A;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetCoefficients - Sets the coefficient matrices that define the linear matrix
   equation to be solved.

   Collective

   Input Parameters:
+  lme - the matrix function context
.  A   - first coefficient matrix
.  B   - second coefficient matrix
.  D   - third coefficient matrix
-  E   - fourth coefficient matrix

   Notes:
   The matrix equation takes the general form A*X*E+D*X*B=C, where matrix C is not
   provided here but with LMESetRHS(). Not all four matrices must be passed, some
   can be NULL instead, see LMESetProblemType() for details.

   It must be called before LMESetUp(). If it is called again after LMESetUp() then
   the LME object is reset.

   In order to delete a previously set matrix, pass a NULL in the corresponding
   argument.

   Level: beginner

.seealso: LMESolve(), LMESetUp(), LMESetRHS(), LMESetProblemType()
@*/
PetscErrorCode LMESetCoefficients(LME lme,Mat A,Mat B,Mat D,Mat E)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(lme,1,A,2);
  if (B) {
    PetscValidHeaderSpecific(B,MAT_CLASSID,3);
    PetscCheckSameComm(lme,1,B,3);
  }
  if (D) {
    PetscValidHeaderSpecific(D,MAT_CLASSID,4);
    PetscCheckSameComm(lme,1,D,4);
  }
  if (E) {
    PetscValidHeaderSpecific(E,MAT_CLASSID,5);
    PetscCheckSameComm(lme,1,E,5);
  }

  if (lme->setupcalled) PetscCall(LMEReset(lme));

  PetscCall(LMESetCoefficients_Private(lme,A,&lme->A));
  if (B) PetscCall(LMESetCoefficients_Private(lme,B,&lme->B));
  else if (!lme->setupcalled) PetscCall(MatDestroy(&lme->B));
  if (D) PetscCall(LMESetCoefficients_Private(lme,D,&lme->D));
  else if (!lme->setupcalled) PetscCall(MatDestroy(&lme->D));
  if (E) PetscCall(LMESetCoefficients_Private(lme,E,&lme->E));
  else if (!lme->setupcalled) PetscCall(MatDestroy(&lme->E));

  lme->setupcalled = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetCoefficients - Gets the coefficient matrices of the matrix equation.

   Collective

   Input Parameter:
.  lme - the LME context

   Output Parameters:
+  A   - first coefficient matrix
.  B   - second coefficient matrix
.  D   - third coefficient matrix
-  E   - fourth coefficient matrix

   Level: intermediate

.seealso: LMESolve(), LMESetCoefficients()
@*/
PetscErrorCode LMEGetCoefficients(LME lme,Mat *A,Mat *B,Mat *D,Mat *E)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (A) *A = lme->A;
  if (B) *B = lme->B;
  if (D) *D = lme->D;
  if (E) *E = lme->E;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetRHS - Sets the right-hand side of the matrix equation, as a low-rank
   matrix.

   Collective

   Input Parameters:
+  lme - the matrix function context
-  C   - the right-hand side matrix

   Notes:
   The matrix equation takes the general form A*X*E+D*X*B=C, where matrix C is
   given with this function. C must be a low-rank matrix of type MATLRC, that is,
   C = U*D*V' where D is diagonal of order k, and U, V are dense tall-skinny
   matrices with k columns. No sparse matrix must be provided when creating the
   MATLRC matrix.

   In equation types that require C to be symmetric, such as Lyapunov, C must be
   created with V=U (or V=NULL).

   Level: beginner

.seealso: LMESetSolution(), LMESetProblemType()
@*/
PetscErrorCode LMESetRHS(LME lme,Mat C)
{
  Mat            A;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidHeaderSpecific(C,MAT_CLASSID,2);
  PetscCheckSameComm(lme,1,C,2);
  PetscCheckTypeName(C,MATLRC);

  PetscCall(MatLRCGetMats(C,&A,NULL,NULL,NULL));
  PetscCheck(!A,PetscObjectComm((PetscObject)C),PETSC_ERR_SUP,"The MatLRC must not have a sparse matrix term");

  PetscCall(PetscObjectReference((PetscObject)C));
  PetscCall(MatDestroy(&lme->C));
  lme->C = C;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetRHS - Gets the right-hand side of the matrix equation.

   Collective

   Input Parameter:
.  lme - the LME context

   Output Parameters:
.  C   - the low-rank matrix

   Level: intermediate

.seealso: LMESolve(), LMESetRHS()
@*/
PetscErrorCode LMEGetRHS(LME lme,Mat *C)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(C,2);
  *C = lme->C;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMESetSolution - Sets the placeholder for the solution of the matrix
   equation, as a low-rank matrix.

   Collective

   Input Parameters:
+  lme - the matrix function context
-  X   - the solution matrix

   Notes:
   The matrix equation takes the general form A*X*E+D*X*B=C, where the solution
   matrix is of low rank and is written in factored form X = U*D*V'. This function
   provides a Mat object of type MATLRC that stores U, V and (optionally) D.
   These factors will be computed during LMESolve().

   In equation types whose solution X is symmetric, such as Lyapunov, X must be
   created with V=U (or V=NULL).

   If the user provides X with this function, then the solver will
   return a solution with rank at most the number of columns of U. Alternatively,
   it is possible to let the solver choose the rank of the solution, by
   setting X to NULL and then calling LMEGetSolution() after LMESolve().

   Level: intermediate

.seealso: LMEGetSolution(), LMESetRHS(), LMESetProblemType(), LMESolve()
@*/
PetscErrorCode LMESetSolution(LME lme,Mat X)
{
  Mat            A;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (X) {
    PetscValidHeaderSpecific(X,MAT_CLASSID,2);
    PetscCheckSameComm(lme,1,X,2);
    PetscCheckTypeName(X,MATLRC);
    PetscCall(MatLRCGetMats(X,&A,NULL,NULL,NULL));
    PetscCheck(!A,PetscObjectComm((PetscObject)X),PETSC_ERR_SUP,"The MatLRC must not have a sparse matrix term");
    PetscCall(PetscObjectReference((PetscObject)X));
  }
  PetscCall(MatDestroy(&lme->X));
  lme->X = X;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetSolution - Gets the solution of the matrix equation.

   Collective

   Input Parameter:
.  lme - the LME context

   Output Parameters:
.  X   - the low-rank matrix

   Level: intermediate

.seealso: LMESolve(), LMESetSolution()
@*/
PetscErrorCode LMEGetSolution(LME lme,Mat *X)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(X,2);
  *X = lme->X;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEAllocateSolution - Allocate memory storage for common variables such
   as the basis vectors.

   Collective

   Input Parameters:
+  lme   - linear matrix equation solver context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developer Notes:
   This is SLEPC_EXTERN because it may be required by user plugin LME
   implementations.

   Level: developer

.seealso: LMESetUp()
@*/
PetscErrorCode LMEAllocateSolution(LME lme,PetscInt extra)
{
  PetscInt       oldsize,requested;
  Vec            t;

  PetscFunctionBegin;
  requested = lme->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  PetscCall(BVGetSizes(lme->V,NULL,NULL,&oldsize));

  /* allocate basis vectors */
  if (!lme->V) PetscCall(LMEGetBV(lme,&lme->V));
  if (!oldsize) {
    if (!((PetscObject)lme->V)->type_name) PetscCall(BVSetType(lme->V,BVMAT));
    PetscCall(MatCreateVecsEmpty(lme->A,&t,NULL));
    PetscCall(BVSetSizesFromVec(lme->V,t,requested));
    PetscCall(VecDestroy(&t));
  } else PetscCall(BVResize(lme->V,requested,PETSC_FALSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}
