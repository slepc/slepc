/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   LME routines related to the solution process
*/

#include <slepc/private/lmeimpl.h>   /*I "slepclme.h" I*/
#include <slepcblaslapack.h>

/*@
   LMESolve - Solves the linear matrix equation.

   Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Options Database Keys:
+  -lme_view             - print information about the solver once the solve is complete
.  -lme_view_pre         - print information about the solver before the solve starts
.  -lme_view_mat_a       - view the A matrix
.  -lme_view_mat_b       - view the B matrix
.  -lme_view_mat_d       - view the D matrix
.  -lme_view_mat_e       - view the E matrix
.  -lme_view_rhs         - view the right hand side of the equation
.  -lme_view_solution    - view the computed solution
-  -lme_converged_reason - print reason for convergence/divergence, and number of iterations

   Notes:
   The matrix coefficients are specified with `LMESetCoefficients()`. The right-hand
   side is specified with `LMESetRHS()`. The placeholder for the solution is specified
   with `LMESetSolution()`.

   All the command-line options listed above admit an optional argument specifying
   the viewer type and options. For instance, use `-lme_view_mat_a binary:amatrix.bin`
   to save the $A$ matrix to a binary file.

   Level: beginner

.seealso: [](ch:lme), `LMECreate()`, `LMESetUp()`, `LMEDestroy()`, `LMESetTolerances()`, `LMESetCoefficients()`, `LMESetRHS()`, `LMESetSolution()`
@*/
PetscErrorCode LMESolve(LME lme)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);

  /* call setup */
  PetscCall(LMESetUp(lme));
  lme->its    = 0;
  lme->errest = 0.0;

  PetscCall(LMEViewFromOptions(lme,NULL,"-lme_view_pre"));

  /* call solver */
  PetscCheck(lme->ops->solve[lme->problem_type],PetscObjectComm((PetscObject)lme),PETSC_ERR_SUP,"The specified solver does not support equation type %s",LMEProblemTypes[lme->problem_type]);
  PetscCall(PetscLogEventBegin(LME_Solve,lme,0,0,0));
  PetscUseTypeMethod(lme,solve[lme->problem_type]);
  PetscCall(PetscLogEventEnd(LME_Solve,lme,0,0,0));

  PetscCheck(lme->reason,PetscObjectComm((PetscObject)lme),PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");

  PetscCheck(!lme->errorifnotconverged || lme->reason>=0,PetscObjectComm((PetscObject)lme),PETSC_ERR_NOT_CONVERGED,"LMESolve has not converged");

  /* various viewers */
  PetscCall(LMEViewFromOptions(lme,NULL,"-lme_view"));
  PetscCall(LMEConvergedReasonViewFromOptions(lme));
  PetscCall(MatViewFromOptions(lme->A,(PetscObject)lme,"-lme_view_mat_a"));
  if (lme->B) PetscCall(MatViewFromOptions(lme->B,(PetscObject)lme,"-lme_view_mat_b"));
  if (lme->D) PetscCall(MatViewFromOptions(lme->D,(PetscObject)lme,"-lme_view_mat_d"));
  if (lme->E) PetscCall(MatViewFromOptions(lme->E,(PetscObject)lme,"-lme_view_mat_e"));
  PetscCall(MatViewFromOptions(lme->C,(PetscObject)lme,"-lme_view_rhs"));
  PetscCall(MatViewFromOptions(lme->X,(PetscObject)lme,"-lme_view_solution"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetIterationNumber - Gets the current iteration number. If the
   call to `LMESolve()` is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  its - number of iterations

   Note:
   During the $i$-th iteration this call returns $i-1$. If `LMESolve()` is
   complete, then parameter `its` contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call `LMEGetConvergedReason()` to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: [](ch:lme), `LMEGetConvergedReason()`, `LMESetTolerances()`
@*/
PetscErrorCode LMEGetIterationNumber(LME lme,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(its,2);
  *its = lme->its;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetConvergedReason - Gets the reason why the `LMESolve()` iteration was
   stopped.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged, see
   `LMEConvergedReason` for the possible values

   Options Database Key:
.  -lme_converged_reason - print reason for convergence/divergence, and number of iterations

   Note:
   If this routine is called before or doing the `LMESolve()` the value of
   `LME_CONVERGED_ITERATING` is returned.

   Level: intermediate

.seealso: [](ch:lme), `LMESetTolerances()`, `LMESolve()`, `LMEConvergedReason`, `LMESetErrorIfNotConverged()`
@*/
PetscErrorCode LMEGetConvergedReason(LME lme,LMEConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(reason,2);
  *reason = lme->reason;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEGetErrorEstimate - Returns the error estimate obtained during the solve.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  errest - the error estimate

   Note:
   This is the error estimated internally by the solver. The actual
   error bound can be computed with `LMEComputeError()`. Note that some
   solvers may not be able to provide an error estimate.

   Level: advanced

.seealso: [](ch:lme), `LMEComputeError()`
@*/
PetscErrorCode LMEGetErrorEstimate(LME lme,PetscReal *errest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(errest,2);
  *errest = lme->errest;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   LMEComputeResidualNorm_Lyapunov - Computes the Frobenius norm of the residual matrix
   associated with the Lyapunov equation.
*/
static PetscErrorCode LMEComputeResidualNorm_Lyapunov(LME lme,PetscReal *norm)
{
  PetscInt          j,n,N,k,l,lda,ldb;
  PetscBLASInt      n_,N_,k_,l_,lda_,ldb_;
  PetscScalar       *Rarray,alpha=1.0,beta=0.0;
  const PetscScalar *A,*B;
  BV                W,AX,X1,C1;
  Mat               R,X1m,C1m;
  Vec               v,w;
  VecScatter        vscat;

  PetscFunctionBegin;
  PetscCall(MatLRCGetMats(lme->C,NULL,&C1m,NULL,NULL));
  PetscCall(MatLRCGetMats(lme->X,NULL,&X1m,NULL,NULL));
  PetscCall(BVCreateFromMat(C1m,&C1));
  PetscCall(BVSetFromOptions(C1));
  PetscCall(BVCreateFromMat(X1m,&X1));
  PetscCall(BVSetFromOptions(X1));
  PetscCall(BVGetSizes(X1,&n,&N,&k));
  PetscCall(BVGetSizes(C1,NULL,NULL,&l));
  PetscCall(PetscBLASIntCast(n,&n_));
  PetscCall(PetscBLASIntCast(N,&N_));
  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(l,&l_));

  /* create W to store a redundant copy of a BV in each process */
  PetscCall(BVCreate(PETSC_COMM_SELF,&W));
  PetscCall(BVSetSizes(W,N,N,k));
  PetscCall(BVSetFromOptions(W));
  PetscCall(BVGetColumn(X1,0,&v));
  PetscCall(VecScatterCreateToAll(v,&vscat,NULL));
  PetscCall(BVRestoreColumn(X1,0,&v));

  /* create AX to hold the product A*X1 */
  PetscCall(BVDuplicate(X1,&AX));
  PetscCall(BVMatMult(X1,lme->A,AX));

  /* create dense matrix to hold the residual R=C1*C1'+AX*X1'+X1*AX' */
  PetscCall(MatCreateDense(PetscObjectComm((PetscObject)lme),n,n,N,N,NULL,&R));

  /* R=C1*C1' */
  PetscCall(MatDenseGetArrayWrite(R,&Rarray));
  for (j=0;j<l;j++) {
    PetscCall(BVGetColumn(C1,j,&v));
    PetscCall(BVGetColumn(W,j,&w));
    PetscCall(VecScatterBegin(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(BVRestoreColumn(C1,j,&v));
    PetscCall(BVRestoreColumn(W,j,&w));
  }
  if (n) {
    PetscCall(BVGetArrayRead(C1,&A));
    PetscCall(BVGetLeadingDimension(C1,&lda));
    PetscCall(PetscBLASIntCast(lda,&lda_));
    PetscCall(BVGetArrayRead(W,&B));
    PetscCall(BVGetLeadingDimension(W,&ldb));
    PetscCall(PetscBLASIntCast(ldb,&ldb_));
    PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n_,&N_,&l_,&alpha,(PetscScalar*)A,&lda_,(PetscScalar*)B,&ldb_,&beta,Rarray,&n_));
    PetscCall(BVRestoreArrayRead(C1,&A));
    PetscCall(BVRestoreArrayRead(W,&B));
  }
  beta = 1.0;

  /* R+=AX*X1' */
  for (j=0;j<k;j++) {
    PetscCall(BVGetColumn(X1,j,&v));
    PetscCall(BVGetColumn(W,j,&w));
    PetscCall(VecScatterBegin(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(BVRestoreColumn(X1,j,&v));
    PetscCall(BVRestoreColumn(W,j,&w));
  }
  if (n) {
    PetscCall(BVGetArrayRead(AX,&A));
    PetscCall(BVGetLeadingDimension(AX,&lda));
    PetscCall(PetscBLASIntCast(lda,&lda_));
    PetscCall(BVGetArrayRead(W,&B));
    PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n_,&N_,&k_,&alpha,(PetscScalar*)A,&lda_,(PetscScalar*)B,&ldb_,&beta,Rarray,&n_));
    PetscCall(BVRestoreArrayRead(AX,&A));
    PetscCall(BVRestoreArrayRead(W,&B));
  }

  /* R+=X1*AX' */
  for (j=0;j<k;j++) {
    PetscCall(BVGetColumn(AX,j,&v));
    PetscCall(BVGetColumn(W,j,&w));
    PetscCall(VecScatterBegin(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(BVRestoreColumn(AX,j,&v));
    PetscCall(BVRestoreColumn(W,j,&w));
  }
  if (n) {
    PetscCall(BVGetArrayRead(X1,&A));
    PetscCall(BVGetLeadingDimension(AX,&lda));
    PetscCall(PetscBLASIntCast(lda,&lda_));
    PetscCall(BVGetArrayRead(W,&B));
    PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n_,&N_,&k_,&alpha,(PetscScalar*)A,&lda_,(PetscScalar*)B,&ldb_,&beta,Rarray,&n_));
    PetscCall(BVRestoreArrayRead(X1,&A));
    PetscCall(BVRestoreArrayRead(W,&B));
  }
  PetscCall(MatDenseRestoreArrayWrite(R,&Rarray));

  /* compute ||R||_F */
  PetscCall(MatNorm(R,NORM_FROBENIUS,norm));

  PetscCall(BVDestroy(&W));
  PetscCall(VecScatterDestroy(&vscat));
  PetscCall(BVDestroy(&AX));
  PetscCall(MatDestroy(&R));
  PetscCall(BVDestroy(&C1));
  PetscCall(BVDestroy(&X1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   LMEComputeError - Computes the error (based on the residual norm) associated
   with the last equation solved.

   Collective

   Input Parameter:
.  lme  - the linear matrix equation solver context

   Output Parameter:
.  error - the error

   Notes:
   This function is not scalable (in terms of memory or parallel communication),
   so it should not be called except in the case of small problem size. For
   large equations, use `LMEGetErrorEstimate()`.

   Level: advanced

.seealso: [](ch:lme), `LMESolve()`, `LMEGetErrorEstimate()`
@*/
PetscErrorCode LMEComputeError(LME lme,PetscReal *error)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscAssertPointer(error,2);

  PetscCall(PetscLogEventBegin(LME_ComputeError,lme,0,0,0));
  /* compute residual norm */
  switch (lme->problem_type) {
    case LME_LYAPUNOV:
      PetscCall(LMEComputeResidualNorm_Lyapunov(lme,error));
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_SUP,"Not implemented for equation type %s",LMEProblemTypes[lme->problem_type]);
  }

  /* compute error */
  /* currently we only support absolute error, so just return the norm */
  PetscCall(PetscLogEventEnd(LME_ComputeError,lme,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}
