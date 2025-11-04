/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to the solution process

   References:

       [1] C. Campos and J.E. Roman, "NEP: a module for the parallel solution
           of nonlinear eigenvalue problems in SLEPc", ACM Trans. Math. Soft.
           47(3), 23:1--23:29, 2021.
*/

#include <slepc/private/nepimpl.h>       /*I "slepcnep.h" I*/
#include <slepc/private/bvimpl.h>
#include <petscdraw.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-nep,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"{NEP}: a module for the parallel solution of nonlinear eigenvalue problems in {SLEPc}\",\n"
  "   journal = \"{ACM} Trans. Math. Software\",\n"
  "   volume = \"47\",\n"
  "   number = \"3\",\n"
  "   pages = \"23:1--23:29\",\n"
  "   year = \"2021\",\n"
  "   doi = \"10.1145/3447544\"\n"
  "}\n";

PetscErrorCode NEPComputeVectors(NEP nep)
{
  PetscFunctionBegin;
  NEPCheckSolved(nep,1);
  if (nep->state==NEP_STATE_SOLVED) PetscTryTypeMethod(nep,computevectors);
  nep->state = NEP_STATE_EIGENVECTORS;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPSolve - Solves the nonlinear eigenproblem.

   Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Options Database Keys:
+  -nep_view             - print information about the solver once the solve is complete
.  -nep_view_pre         - print information about the solver before the solve starts
.  -nep_view_matk        - view the split form matrix $A_k$ (replace `k` by an integer from 0 to `nt`-1)
.  -nep_view_fnk         - view the split form function $f_k$ (replace `k` by an integer from 0 to `nt`-1)
.  -nep_view_vectors     - view the computed eigenvectors
.  -nep_view_values      - view the computed eigenvalues
.  -nep_converged_reason - print reason for convergence/divergence, and number of iterations
.  -nep_error_absolute   - print absolute errors of each eigenpair
.  -nep_error_relative   - print relative errors of each eigenpair
-  -nep_error_backward   - print backward errors of each eigenpair

   Notes:
   `NEPSolve()` will return without generating an error regardless of whether
   all requested solutions were computed or not. Call `NEPGetConverged()` to get the
   actual number of computed solutions, and `NEPGetConvergedReason()` to determine if
   the solver converged or failed and why.

   All the command-line options listed above admit an optional argument specifying
   the viewer type and options. For instance, use `-nep_view_vectors binary:myvecs.bin`
   to save the eigenvectors to a binary file, `-nep_view_values draw` to draw the computed
   eigenvalues graphically, or `-nep_error_relative :myerr.m:ascii_matlab` to save
   the errors in a file that can be executed in Matlab.

   Level: beginner

.seealso: [](ch:nep), `NEPCreate()`, `NEPSetUp()`, `NEPDestroy()`, `NEPSetTolerances()`, `NEPGetConverged()`, `NEPGetConvergedReason()`
@*/
PetscErrorCode NEPSolve(NEP nep)
{
  PetscInt       i;
  char           str[16];

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (nep->state>=NEP_STATE_SOLVED) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscCitationsRegister(citation,&cited));
  PetscCall(PetscLogEventBegin(NEP_Solve,nep,0,0,0));

  /* call setup */
  PetscCall(NEPSetUp(nep));
  nep->nconv = 0;
  nep->its = 0;
  for (i=0;i<nep->ncv;i++) {
    nep->eigr[i]   = 0.0;
    nep->eigi[i]   = 0.0;
    nep->errest[i] = 0.0;
    nep->perm[i]   = i;
  }
  PetscCall(NEPViewFromOptions(nep,NULL,"-nep_view_pre"));
  PetscCall(RGViewFromOptions(nep->rg,NULL,"-rg_view"));

  /* call solver */
  PetscUseTypeMethod(nep,solve);
  PetscCheck(nep->reason,PetscObjectComm((PetscObject)nep),PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");
  nep->state = NEP_STATE_SOLVED;

  /* Only the first nconv columns contain useful information */
  PetscCall(BVSetActiveColumns(nep->V,0,nep->nconv));
  if (nep->twosided) PetscCall(BVSetActiveColumns(nep->W,0,nep->nconv));

  if (nep->refine==NEP_REFINE_SIMPLE && nep->rits>0 && nep->nconv>0) {
    PetscCall(NEPComputeVectors(nep));
    PetscCall(NEPNewtonRefinementSimple(nep,&nep->rits,nep->rtol,nep->nconv));
    nep->state = NEP_STATE_EIGENVECTORS;
  }

  /* sort eigenvalues according to nep->which parameter */
  PetscCall(SlepcSortEigenvalues(nep->sc,nep->nconv,nep->eigr,nep->eigi,nep->perm));
  PetscCall(PetscLogEventEnd(NEP_Solve,nep,0,0,0));

  /* various viewers */
  PetscCall(NEPViewFromOptions(nep,NULL,"-nep_view"));
  PetscCall(NEPConvergedReasonViewFromOptions(nep));
  PetscCall(NEPErrorViewFromOptions(nep));
  PetscCall(NEPValuesViewFromOptions(nep));
  PetscCall(NEPVectorsViewFromOptions(nep));
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    for (i=0;i<nep->nt;i++) {
      PetscCall(PetscSNPrintf(str,sizeof(str),"-nep_view_mat%" PetscInt_FMT,i));
      PetscCall(MatViewFromOptions(nep->A[i],(PetscObject)nep,str));
      PetscCall(PetscSNPrintf(str,sizeof(str),"-nep_view_fn%" PetscInt_FMT,i));
      PetscCall(FNViewFromOptions(nep->f[i],(PetscObject)nep,str));
    }
  }

  /* Remove the initial subspace */
  nep->nini = 0;

  /* Reset resolvent information */
  PetscCall(MatDestroy(&nep->resolvent));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPProjectOperator - Computes the projection of the nonlinear operator.

   Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
.  j0  - initial index
-  j1  - final index

   Notes:
   This is available for split operator only.

   The nonlinear operator $T(\lambda)$ is projected onto $\operatorname{span}(V)$, where $V$ is
   an orthonormal basis built internally by the solver. The projected
   operator is equal to $\sum_i V^* A_i V f_i(\lambda)$, so this function
   computes all matrices $E_i = V^* A_i V$, and stores them in the extra
   matrices inside `DS`, see `DSMatType`. Only rows/columns in the range [`j0`,`j1`-1] are computed,
   the previous ones are assumed to be available already.

   Level: developer

.seealso: [](ch:nep), `NEPSetSplitOperator()`
@*/
PetscErrorCode NEPProjectOperator(NEP nep,PetscInt j0,PetscInt j1)
{
  PetscInt       k;
  Mat            G;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,j0,2);
  PetscValidLogicalCollectiveInt(nep,j1,3);
  NEPCheckProblem(nep,1);
  NEPCheckSplit(nep,1);
  PetscCall(BVSetActiveColumns(nep->V,j0,j1));
  for (k=0;k<nep->nt;k++) {
    PetscCall(DSGetMat(nep->ds,DSMatExtra[k],&G));
    PetscCall(BVMatProject(nep->V,nep->A[k],nep->V,G));
    PetscCall(DSRestoreMat(nep->ds,DSMatExtra[k],&G));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPApplyFunction - Applies the nonlinear function $T(\lambda)$ to a given vector.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  lambda - scalar argument
.  x      - vector to be multiplied against
-  v      - workspace vector (used only in the case of split form)

   Output Parameters:
+  y   - result vector
.  A   - (optional) Function matrix, for callback interface only
-  B   - (unused) preconditioning matrix

   Note:
   If the nonlinear operator is represented in split form, the result
   $y = T(\lambda)x$ is computed without building $T(\lambda)$ explicitly. In
   that case, parameters `A` and `B` are not used. Otherwise, the matrix
   $T(\lambda)$ is built and the effect is the same as a call to
   `NEPComputeFunction()` followed by a `MatMult()`.

   Level: developer

.seealso: [](ch:nep), `NEPSetSplitOperator()`, `NEPComputeFunction()`, `NEPApplyAdjoint()`
@*/
PetscErrorCode NEPApplyFunction(NEP nep,PetscScalar lambda,Vec x,Vec v,Vec y,Mat A,Mat B)
{
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,lambda,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  if (v) PetscValidHeaderSpecific(v,VEC_CLASSID,4);
  PetscValidHeaderSpecific(y,VEC_CLASSID,5);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,6);
  if (B) PetscValidHeaderSpecific(B,MAT_CLASSID,7);

  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    PetscCall(VecSet(y,0.0));
    for (i=0;i<nep->nt;i++) {
      PetscCall(FNEvaluateFunction(nep->f[i],lambda,&alpha));
      PetscCall(MatMult(nep->A[i],x,v));
      PetscCall(VecAXPY(y,alpha,v));
    }
  } else {
    if (!A) A = nep->function;
    PetscCall(NEPComputeFunction(nep,lambda,A,A));
    PetscCall(MatMult(A,x,y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPApplyAdjoint - Applies the adjoint nonlinear function $T^*(\lambda)$ to a given vector.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  lambda - scalar argument
.  x      - vector to be multiplied against
-  v      - workspace vector (used only in the case of split form)

   Output Parameters:
+  y   - result vector
.  A   - (optional) Function matrix, for callback interface only
-  B   - (unused) preconditioning matrix

   Note:
   If the nonlinear operator is represented in split form, the result
   $y = T^*(\lambda)x$ is computed without building $T(\lambda)$ explicitly. In
   that case, parameters `A` and `B` are not used. Otherwise, the matrix
   $T(\lambda)$ is built and the effect is the same as a call to
   `NEPComputeFunction()` followed by a `MatMultHermitianTranspose()`.

   Level: developer

.seealso: [](ch:nep), `NEPSetSplitOperator()`, `NEPComputeFunction()`, `NEPApplyFunction()`
@*/
PetscErrorCode NEPApplyAdjoint(NEP nep,PetscScalar lambda,Vec x,Vec v,Vec y,Mat A,Mat B)
{
  PetscInt       i;
  PetscScalar    alpha;
  Vec            w;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,lambda,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  if (v) PetscValidHeaderSpecific(v,VEC_CLASSID,4);
  PetscValidHeaderSpecific(y,VEC_CLASSID,5);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,6);
  if (B) PetscValidHeaderSpecific(B,MAT_CLASSID,7);

  PetscCall(VecDuplicate(x,&w));
  PetscCall(VecCopy(x,w));
  PetscCall(VecConjugate(w));
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    PetscCall(VecSet(y,0.0));
    for (i=0;i<nep->nt;i++) {
      PetscCall(FNEvaluateFunction(nep->f[i],lambda,&alpha));
      PetscCall(MatMultHermitianTranspose(nep->A[i],w,v));
      PetscCall(VecAXPY(y,alpha,v));
    }
  } else {
    if (!A) A = nep->function;
    PetscCall(NEPComputeFunction(nep,lambda,A,A));
    PetscCall(MatMultHermitianTranspose(A,w,y));
  }
  PetscCall(VecDestroy(&w));
  PetscCall(VecConjugate(y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPApplyJacobian - Applies the nonlinear Jacobian $T'(\lambda)$ to a given vector.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  lambda - scalar argument
.  x      - vector to be multiplied against
-  v      - workspace vector (used only in the case of split form)

   Output Parameters:
+  y   - result vector
-  A   - (optional) Jacobian matrix, for callback interface only

   Note:
   If the nonlinear operator is represented in split form, the result
   $y = T'(\lambda)x$ is computed without building $T'(\lambda)$ explicitly. In
   that case, parameter `A` is not used. Otherwise, the matrix
   $T'(\lambda)$ is built and the effect is the same as a call to
   `NEPComputeJacobian()` followed by a `MatMult()`.

   Level: developer

.seealso: [](ch:nep), `NEPSetSplitOperator()`, `NEPComputeJacobian()`
@*/
PetscErrorCode NEPApplyJacobian(NEP nep,PetscScalar lambda,Vec x,Vec v,Vec y,Mat A)
{
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,lambda,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  if (v) PetscValidHeaderSpecific(v,VEC_CLASSID,4);
  PetscValidHeaderSpecific(y,VEC_CLASSID,5);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,6);

  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    PetscCall(VecSet(y,0.0));
    for (i=0;i<nep->nt;i++) {
      PetscCall(FNEvaluateDerivative(nep->f[i],lambda,&alpha));
      PetscCall(MatMult(nep->A[i],x,v));
      PetscCall(VecAXPY(y,alpha,v));
    }
  } else {
    if (!A) A = nep->jacobian;
    PetscCall(NEPComputeJacobian(nep,lambda,A));
    PetscCall(MatMult(A,x,y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetIterationNumber - Gets the current iteration number. If the
   call to `NEPSolve()` is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  its - number of iterations

   Note:
   During the $i$-th iteration this call returns $i-1$. If `NEPSolve()` is
   complete, then parameter `its` contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call `NEPGetConvergedReason()` to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: [](ch:nep), `NEPGetConvergedReason()`, `NEPSetTolerances()`
@*/
PetscErrorCode NEPGetIterationNumber(NEP nep,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(its,2);
  *its = nep->its;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetConverged - Gets the number of converged eigenpairs.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  nconv - number of converged eigenpairs

   Note:
   This function should be called after `NEPSolve()` has finished.

   The value `nconv` may be different from the number of requested solutions
   `nev`, but not larger than `ncv`, see `NEPSetDimensions()`.

   Level: beginner

.seealso: [](ch:nep), `NEPSetDimensions()`, `NEPSolve()`, `NEPGetEigenpair()`
@*/
PetscErrorCode NEPGetConverged(NEP nep,PetscInt *nconv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(nconv,2);
  NEPCheckSolved(nep,1);
  *nconv = nep->nconv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetConvergedReason - Gets the reason why the `NEPSolve()` iteration was
   stopped.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged, see
   `NEPConvergedReason` for the possible values

   Options Database Key:
.  -nep_converged_reason - print reason for convergence/divergence, and number of iterations

   Note:
   If this routine is called before or doing the `NEPSolve()` the value of
   `NEP_CONVERGED_ITERATING` is returned.

   Level: intermediate

.seealso: [](ch:nep), `NEPSetTolerances()`, `NEPSolve()`, `NEPConvergedReason`
@*/
PetscErrorCode NEPGetConvergedReason(NEP nep,NEPConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(reason,2);
  NEPCheckSolved(nep,1);
  *reason = nep->reason;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetEigenpair - Gets the `i`-th solution of the eigenproblem as computed by
   `NEPSolve()`. The solution consists in both the eigenvalue and the eigenvector.

   Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  i   - index of the solution

   Output Parameters:
+  eigr - real part of eigenvalue
.  eigi - imaginary part of eigenvalue
.  Vr   - real part of eigenvector
-  Vi   - imaginary part of eigenvector

   Notes:
   It is allowed to pass `NULL` for `Vr` and `Vi`, if the eigenvector is not
   required. Otherwise, the caller must provide valid `Vec` objects, i.e.,
   they must be created by the calling program with e.g. `MatCreateVecs()`.

   If the eigenvalue is real, then `eigi` and `Vi` are set to zero. If PETSc is
   configured with complex scalars the eigenvalue is stored
   directly in `eigr` (`eigi` is set to zero) and the eigenvector in `Vr` (`Vi`
   is set to zero). In any case, the user can pass `NULL` in any argument that
   is not required.

   The index `i` should be a value between 0 and `nconv`-1 (see `NEPGetConverged()`).
   Eigenpairs are indexed according to the ordering criterion established
   with `NEPSetWhichEigenpairs()`.

   The eigenvector is normalized to have unit norm.

   Level: beginner

.seealso: [](ch:nep), `NEPSolve()`, `NEPGetConverged()`, `NEPSetWhichEigenpairs()`, `NEPGetLeftEigenvector()`
@*/
PetscErrorCode NEPGetEigenpair(NEP nep,PetscInt i,PetscScalar *eigr,PetscScalar *eigi,Vec Vr,Vec Vi)
{
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,i,2);
  if (Vr) { PetscValidHeaderSpecific(Vr,VEC_CLASSID,5); PetscCheckSameComm(nep,1,Vr,5); }
  if (Vi) { PetscValidHeaderSpecific(Vi,VEC_CLASSID,6); PetscCheckSameComm(nep,1,Vi,6); }
  NEPCheckSolved(nep,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCheck(i<nep->nconv,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see NEPGetConverged()");

  PetscCall(NEPComputeVectors(nep));
  k = nep->perm[i];

  /* eigenvalue */
#if defined(PETSC_USE_COMPLEX)
  if (eigr) *eigr = nep->eigr[k];
  if (eigi) *eigi = 0;
#else
  if (eigr) *eigr = nep->eigr[k];
  if (eigi) *eigi = nep->eigi[k];
#endif

  /* eigenvector */
  PetscCall(BV_GetEigenvector(nep->V,k,nep->eigi[k],Vr,Vi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetLeftEigenvector - Gets the `i`-th left eigenvector as computed by `NEPSolve()`.

   Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  i   - index of the solution

   Output Parameters:
+  Wr   - real part of left eigenvector
-  Wi   - imaginary part of left eigenvector

   Notes:
   The caller must provide valid `Vec` objects, i.e., they must be created
   by the calling program with e.g. `MatCreateVecs()`.

   If the corresponding eigenvalue is real, then `Wi` is set to zero. If PETSc is
   configured with complex scalars the eigenvector is stored directly in `Wr`
   (`Wi` is set to zero). In any case, the user can pass `NULL` in `Wr` or `Wi` if
   one of them is not required.

   The index `i` should be a value between 0 and `nconv`-1 (see `NEPGetConverged()`).
   Eigensolutions are indexed according to the ordering criterion established
   with `NEPSetWhichEigenpairs()`.

   Left eigenvectors are available only if the twosided flag was set, see
   `NEPSetTwoSided()`.

   Level: intermediate

.seealso: [](ch:nep), `NEPGetEigenpair()`, `NEPGetConverged()`, `NEPSetWhichEigenpairs()`, `NEPSetTwoSided()`
@*/
PetscErrorCode NEPGetLeftEigenvector(NEP nep,PetscInt i,Vec Wr,Vec Wi)
{
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,i,2);
  if (Wr) { PetscValidHeaderSpecific(Wr,VEC_CLASSID,3); PetscCheckSameComm(nep,1,Wr,3); }
  if (Wi) { PetscValidHeaderSpecific(Wi,VEC_CLASSID,4); PetscCheckSameComm(nep,1,Wi,4); }
  NEPCheckSolved(nep,1);
  PetscCheck(nep->twosided,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONGSTATE,"Must request left vectors with NEPSetTwoSided");
  PetscCheck(i>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCheck(i<nep->nconv,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see NEPGetConverged()");
  PetscCall(NEPComputeVectors(nep));
  k = nep->perm[i];
  PetscCall(BV_GetEigenvector(nep->W,k,nep->eigi[k],Wr,Wi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPGetErrorEstimate - Returns the error estimate associated to the `i`-th
   computed eigenpair.

   Not Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  i   - index of eigenpair

   Output Parameter:
.  errest - the error estimate

   Note:
   This is the error estimate used internally by the eigensolver. The actual
   error bound can be computed with `NEPComputeError()`.

   Level: advanced

.seealso: [](ch:nep), `NEPComputeError()`
@*/
PetscErrorCode NEPGetErrorEstimate(NEP nep,PetscInt i,PetscReal *errest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscAssertPointer(errest,3);
  NEPCheckSolved(nep,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCheck(i<nep->nconv,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see NEPGetConverged()");
  *errest = nep->errest[nep->perm[i]];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   NEPComputeResidualNorm_Private - Computes the norm of the residual vector
   associated with an eigenpair.

   Input Parameters:
     adj    - whether the adjoint T^* must be used instead of T
     lambda - eigenvalue
     x      - eigenvector
     w      - array of work vectors (two vectors in split form, one vector otherwise)
*/
PetscErrorCode NEPComputeResidualNorm_Private(NEP nep,PetscBool adj,PetscScalar lambda,Vec x,Vec *w,PetscReal *norm)
{
  Vec            y,z=NULL;

  PetscFunctionBegin;
  y = w[0];
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) z = w[1];
  if (adj) PetscCall(NEPApplyAdjoint(nep,lambda,x,z,y,NULL,NULL));
  else PetscCall(NEPApplyFunction(nep,lambda,x,z,y,NULL,NULL));
  PetscCall(VecNorm(y,NORM_2,norm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPComputeError - Computes the error (based on the residual norm) associated
   with the `i`-th computed eigenpair.

   Collective

   Input Parameters:
+  nep  - the nonlinear eigensolver context
.  i    - the solution index
-  type - the type of error to compute, see `NEPErrorType`

   Output Parameter:
.  error - the error

   Notes:
   The error can be computed in various ways, all of them based on the residual
   norm computed as $\|T(\lambda)x\|_2$ where $(\lambda,x)$ is the approximate eigenpair.

   If the computation of left eigenvectors was enabled with `NEPSetTwoSided()`,
   then the error will be computed using the maximum of the value above and
   the left residual norm $\|y^*T(\lambda)\|_2$, where $y$ is the approximate left
   eigenvector.

   Level: beginner

.seealso: [](ch:nep), `NEPErrorType`, `NEPSolve()`, `NEPGetErrorEstimate()`, `NEPSetTwoSided()`
@*/
PetscErrorCode NEPComputeError(NEP nep,PetscInt i,NEPErrorType type,PetscReal *error)
{
  Vec            xr,xi=NULL;
  PetscInt       j,nwork,issplit=0;
  PetscScalar    kr,ki,s;
  PetscReal      er,z=0.0,errorl,nrm;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,i,2);
  PetscValidLogicalCollectiveEnum(nep,type,3);
  PetscAssertPointer(error,4);
  NEPCheckSolved(nep,1);

  /* allocate work vectors */
#if defined(PETSC_USE_COMPLEX)
  nwork = 2;
#else
  nwork = 3;
#endif
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    issplit = 1;
    nwork++;  /* need an extra work vector for NEPComputeResidualNorm_Private */
  }
  PetscCall(NEPSetWorkVecs(nep,nwork));
  xr = nep->work[issplit+1];
#if !defined(PETSC_USE_COMPLEX)
  xi = nep->work[issplit+2];
#endif

  /* compute residual norms */
  PetscCall(NEPGetEigenpair(nep,i,&kr,&ki,xr,xi));
#if !defined(PETSC_USE_COMPLEX)
  PetscCheck(ki==0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Not implemented for complex eigenvalues with real scalars");
#endif
  PetscCall(NEPComputeResidualNorm_Private(nep,PETSC_FALSE,kr,xr,nep->work,error));
  PetscCall(VecNorm(xr,NORM_2,&er));

  /* if two-sided, compute left residual norm and take the maximum */
  if (nep->twosided) {
    PetscCall(NEPGetLeftEigenvector(nep,i,xr,xi));
    PetscCall(NEPComputeResidualNorm_Private(nep,PETSC_TRUE,kr,xr,nep->work,&errorl));
    *error = PetscMax(*error,errorl);
  }

  /* compute error */
  switch (type) {
    case NEP_ERROR_ABSOLUTE:
      break;
    case NEP_ERROR_RELATIVE:
      *error /= PetscAbsScalar(kr)*er;
      break;
    case NEP_ERROR_BACKWARD:
      if (nep->fui!=NEP_USER_INTERFACE_SPLIT) {
        PetscCall(NEPComputeFunction(nep,kr,nep->function,nep->function));
        PetscCall(MatHasOperation(nep->function,MATOP_NORM,&flg));
        PetscCheck(flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
        PetscCall(MatNorm(nep->function,NORM_INFINITY,&nrm));
        *error /= nrm*er;
        break;
      }
      /* initialization of matrix norms */
      if (!nep->nrma[0]) {
        for (j=0;j<nep->nt;j++) {
          PetscCall(MatHasOperation(nep->A[j],MATOP_NORM,&flg));
          PetscCheck(flg,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
          PetscCall(MatNorm(nep->A[j],NORM_INFINITY,&nep->nrma[j]));
        }
      }
      for (j=0;j<nep->nt;j++) {
        PetscCall(FNEvaluateFunction(nep->f[j],kr,&s));
        z = z + nep->nrma[j]*PetscAbsScalar(s);
      }
      *error /= z*er;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid error type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPComputeFunction - Computes the function matrix $T(\lambda)$ that has been
   set with `NEPSetFunction()`.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
-  lambda - the scalar argument

   Output Parameters:
+  A   - Function matrix
-  B   - optional preconditioning matrix

   Note:
   `NEPComputeFunction()` is typically used within nonlinear eigensolvers
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.seealso: [](ch:nep), `NEPSetFunction()`, `NEPGetFunction()`, `NEPComputeJacobian()`
@*/
PetscErrorCode NEPComputeFunction(NEP nep,PetscScalar lambda,Mat A,Mat B)
{
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckProblem(nep,1);
  switch (nep->fui) {
  case NEP_USER_INTERFACE_CALLBACK:
    PetscCheck(nep->computefunction,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER,"Must call NEPSetFunction() first");
    PetscCall(PetscLogEventBegin(NEP_FunctionEval,nep,A,B,0));
    PetscCallBack("NEP user Function function",(*nep->computefunction)(nep,lambda,A,B,nep->functionctx));
    PetscCall(PetscLogEventEnd(NEP_FunctionEval,nep,A,B,0));
    break;
  case NEP_USER_INTERFACE_SPLIT:
    PetscCall(MatZeroEntries(A));
    if (A != B) PetscCall(MatZeroEntries(B));
    for (i=0;i<nep->nt;i++) {
      PetscCall(FNEvaluateFunction(nep->f[i],lambda,&alpha));
      PetscCall(MatAXPY(A,alpha,nep->A[i],nep->mstr));
      if (A != B) PetscCall(MatAXPY(B,alpha,nep->P[i],nep->mstrp));
    }
    break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPComputeJacobian - Computes the Jacobian matrix $T'(\lambda)$ that has been
   set with `NEPSetJacobian()`.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
-  lambda - the scalar argument

   Output Parameter:
.  A   - Jacobian matrix

   Notes:
   Most users should not need to explicitly call this routine, as it
   is used internally within the nonlinear eigensolvers.

   Level: developer

.seealso: [](ch:nep), `NEPSetJacobian()`, `NEPGetJacobian()`, `NEPComputeFunction()`
@*/
PetscErrorCode NEPComputeJacobian(NEP nep,PetscScalar lambda,Mat A)
{
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckProblem(nep,1);
  switch (nep->fui) {
  case NEP_USER_INTERFACE_CALLBACK:
    PetscCheck(nep->computejacobian,PetscObjectComm((PetscObject)nep),PETSC_ERR_USER,"Must call NEPSetJacobian() first");
    PetscCall(PetscLogEventBegin(NEP_JacobianEval,nep,A,0,0));
    PetscCallBack("NEP user Jacobian function",(*nep->computejacobian)(nep,lambda,A,nep->jacobianctx));
    PetscCall(PetscLogEventEnd(NEP_JacobianEval,nep,A,0,0));
    break;
  case NEP_USER_INTERFACE_SPLIT:
    PetscCall(MatZeroEntries(A));
    for (i=0;i<nep->nt;i++) {
      PetscCall(FNEvaluateDerivative(nep->f[i],lambda,&alpha));
      PetscCall(MatAXPY(A,alpha,nep->A[i],nep->mstr));
    }
    break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
