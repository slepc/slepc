/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   PEP routines related to the solution process

   References:

       [1] C. Campos and J.E. Roman, "Parallel iterative refinement in
           polynomial eigenvalue problems", Numer. Linear Algebra Appl.
           23(4):730-745, 2016.
*/

#include <slepc/private/pepimpl.h>       /*I "slepcpep.h" I*/
#include <slepc/private/bvimpl.h>
#include <petscdraw.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-pep-refine,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"Parallel iterative refinement in polynomial eigenvalue problems\",\n"
  "   journal = \"Numer. Linear Algebra Appl.\",\n"
  "   volume = \"23\",\n"
  "   number = \"4\",\n"
  "   pages = \"730--745\",\n"
  "   year = \"2016,\"\n"
  "   doi = \"https://doi.org/10.1002/nla.2052\"\n"
  "}\n";

PetscErrorCode PEPComputeVectors(PEP pep)
{
  PetscFunctionBegin;
  PEPCheckSolved(pep,1);
  if (pep->state==PEP_STATE_SOLVED) PetscTryTypeMethod(pep,computevectors);
  pep->state = PEP_STATE_EIGENVECTORS;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPExtractVectors(PEP pep)
{
  PetscFunctionBegin;
  PEPCheckSolved(pep,1);
  if (pep->state==PEP_STATE_SOLVED) PetscTryTypeMethod(pep,extractvectors);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPSolve - Solves the polynomial eigenproblem.

   Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Options Database Keys:
+  -pep_view             - print information about the solver once the solve is complete
.  -pep_view_pre         - print information about the solver before the solve starts
.  -pep_view_matk        - view the coefficient matrix $A_k$ (replace `k` by an integer from 0 to `nmat`-1)
.  -pep_view_vectors     - view the computed eigenvectors
.  -pep_view_values      - view the computed eigenvalues
.  -pep_converged_reason - print reason for convergence/divergence, and number of iterations
.  -pep_error_absolute   - print absolute errors of each eigenpair
.  -pep_error_relative   - print relative errors of each eigenpair
-  -pep_error_backward   - print backward errors of each eigenpair

   Notes:
   The problem matrices are specified with `PEPSetOperators()`.

   `PEPSolve()` will return without generating an error regardless of whether
   all requested solutions were computed or not. Call `PEPGetConverged()` to get the
   actual number of computed solutions, and `PEPGetConvergedReason()` to determine if
   the solver converged or failed and why.

   All the command-line options listed above admit an optional argument specifying
   the viewer type and options. For instance, use `-pep_view_mat0 binary:matrix0.bin`
   to save the $A_0$ matrix to a binary file, `-pep_view_values draw` to draw the computed
   eigenvalues graphically, or `-pep_error_relative :myerr.m:ascii_matlab` to save
   the errors in a file that can be executed in Matlab.

   Level: beginner

.seealso: [](ch:pep), `PEPCreate()`, `PEPSetUp()`, `PEPDestroy()`, `PEPSetTolerances()`, `PEPGetConverged()`, `PEPGetConvergedReason()`
@*/
PetscErrorCode PEPSolve(PEP pep)
{
  PetscInt       i,k;
  PetscBool      flg,islinear;
  char           str[16];

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (pep->state>=PEP_STATE_SOLVED) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscLogEventBegin(PEP_Solve,pep,0,0,0));

  /* call setup */
  PetscCall(PEPSetUp(pep));
  pep->nconv = 0;
  pep->its   = 0;
  k = pep->lineariz? pep->ncv: pep->ncv*(pep->nmat-1);
  for (i=0;i<k;i++) {
    pep->eigr[i]   = 0.0;
    pep->eigi[i]   = 0.0;
    pep->errest[i] = 0.0;
    pep->perm[i]   = i;
  }
  PetscCall(PEPViewFromOptions(pep,NULL,"-pep_view_pre"));
  PetscCall(RGViewFromOptions(pep->rg,NULL,"-rg_view"));

  /* Call solver */
  PetscUseTypeMethod(pep,solve);
  PetscCheck(pep->reason,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");
  pep->state = PEP_STATE_SOLVED;

  /* Only the first nconv columns contain useful information */
  PetscCall(BVSetActiveColumns(pep->V,0,pep->nconv));

  PetscCall(PetscObjectTypeCompare((PetscObject)pep,PEPLINEAR,&islinear));
  if (!islinear) {
    PetscCall(STPostSolve(pep->st));
    /* Map eigenvalues back to the original problem */
    PetscCall(STGetTransform(pep->st,&flg));
    if (flg) PetscTryTypeMethod(pep,backtransform);
  }

#if !defined(PETSC_USE_COMPLEX)
  /* reorder conjugate eigenvalues (positive imaginary first) */
  for (i=0;i<pep->nconv-1;i++) {
    if (pep->eigi[i] != 0) {
      if (pep->eigi[i] < 0) {
        pep->eigi[i] = -pep->eigi[i];
        pep->eigi[i+1] = -pep->eigi[i+1];
        /* the next correction only works with eigenvectors */
        PetscCall(PEPComputeVectors(pep));
        PetscCall(BVScaleColumn(pep->V,i+1,-1.0));
      }
      i++;
    }
  }
#endif

  if (pep->refine!=PEP_REFINE_NONE) PetscCall(PetscCitationsRegister(citation,&cited));

  if (pep->refine==PEP_REFINE_SIMPLE && pep->rits>0 && pep->nconv>0) {
    PetscCall(PEPComputeVectors(pep));
    PetscCall(PEPNewtonRefinementSimple(pep,&pep->rits,pep->rtol,pep->nconv));
  }

  /* sort eigenvalues according to pep->which parameter */
  PetscCall(SlepcSortEigenvalues(pep->sc,pep->nconv,pep->eigr,pep->eigi,pep->perm));
  PetscCall(PetscLogEventEnd(PEP_Solve,pep,0,0,0));

  /* various viewers */
  PetscCall(PEPViewFromOptions(pep,NULL,"-pep_view"));
  PetscCall(PEPConvergedReasonViewFromOptions(pep));
  PetscCall(PEPErrorViewFromOptions(pep));
  PetscCall(PEPValuesViewFromOptions(pep));
  PetscCall(PEPVectorsViewFromOptions(pep));
  for (i=0;i<pep->nmat;i++) {
    PetscCall(PetscSNPrintf(str,sizeof(str),"-pep_view_mat%" PetscInt_FMT,i));
    PetscCall(MatViewFromOptions(pep->A[i],(PetscObject)pep,str));
  }

  /* Remove the initial subspace */
  pep->nini = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetIterationNumber - Gets the current iteration number. If the
   call to `PEPSolve()` is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameter:
.  its - number of iterations

   Note:
   During the $i$-th iteration this call returns $i-1$. If `PEPSolve()` is
   complete, then parameter `its` contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call `PEPGetConvergedReason()` to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: [](ch:pep), `PEPGetConvergedReason()`, `PEPSetTolerances()`
@*/
PetscErrorCode PEPGetIterationNumber(PEP pep,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(its,2);
  *its = pep->its;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetConverged - Gets the number of converged eigenpairs.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameter:
.  nconv - number of converged eigenpairs

   Notes:
   This function should be called after `PEPSolve()` has finished.

   The value `nconv` may be different from the number of requested solutions
   `nev`, but not larger than `ncv`, see `PEPSetDimensions()`.

   Level: beginner

.seealso: [](ch:pep), `PEPSetDimensions()`, `PEPSolve()`, `PEPGetEigenpair()`
@*/
PetscErrorCode PEPGetConverged(PEP pep,PetscInt *nconv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(nconv,2);
  PEPCheckSolved(pep,1);
  *nconv = pep->nconv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetConvergedReason - Gets the reason why the `PEPSolve()` iteration was
   stopped.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged, see
   `PEPConvergedReason` for the possible values

   Options Database Key:
.  -pep_converged_reason - print reason for convergence/divergence, and number of iterations

   Note:
   If this routine is called before or doing the `PEPSolve()` the value of
   `PEP_CONVERGED_ITERATING` is returned.

   Level: intermediate

.seealso: [](ch:pep), `PEPSetTolerances()`, `PEPSolve()`, `PEPConvergedReason`
@*/
PetscErrorCode PEPGetConvergedReason(PEP pep,PEPConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(reason,2);
  PEPCheckSolved(pep,1);
  *reason = pep->reason;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetEigenpair - Gets the `i`-th solution of the eigenproblem as computed by
   `PEPSolve()`. The solution consists in both the eigenvalue and the eigenvector.

   Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
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
   directly in `eigr` (`eigi` is set to zero) and the eigenvector in `Vr` (`Vi` is
   set to zero). In any case, the user can pass `NULL` in `Vr` or `Vi` if one of
   them is not required.

   The index `i` should be a value between 0 and `nconv`-1 (see `PEPGetConverged()`).
   Eigenpairs are indexed according to the ordering criterion established
   with `PEPSetWhichEigenpairs()`.

   The eigenvector is normalized to have unit norm.

   Level: beginner

.seealso: [](ch:pep), `PEPSolve()`, `PEPGetConverged()`, `PEPSetWhichEigenpairs()`
@*/
PetscErrorCode PEPGetEigenpair(PEP pep,PetscInt i,PetscScalar *eigr,PetscScalar *eigi,Vec Vr,Vec Vi)
{
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,i,2);
  if (Vr) { PetscValidHeaderSpecific(Vr,VEC_CLASSID,5); PetscCheckSameComm(pep,1,Vr,5); }
  if (Vi) { PetscValidHeaderSpecific(Vi,VEC_CLASSID,6); PetscCheckSameComm(pep,1,Vi,6); }
  PEPCheckSolved(pep,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCheck(i<pep->nconv,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see PEPGetConverged()");

  PetscCall(PEPComputeVectors(pep));
  k = pep->perm[i];

  /* eigenvalue */
#if defined(PETSC_USE_COMPLEX)
  if (eigr) *eigr = pep->eigr[k];
  if (eigi) *eigi = 0;
#else
  if (eigr) *eigr = pep->eigr[k];
  if (eigi) *eigi = pep->eigi[k];
#endif

  /* eigenvector */
  PetscCall(BV_GetEigenvector(pep->V,k,pep->eigi[k],Vr,Vi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPGetErrorEstimate - Returns the error estimate associated to the `i`-th
   computed eigenpair.

   Not Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
-  i   - index of eigenpair

   Output Parameter:
.  errest - the error estimate

   Note:
   This is the error estimate used internally by the eigensolver. The actual
   error bound can be computed with `PEPComputeError()`.

   Level: advanced

.seealso: [](ch:pep), `PEPComputeError()`
@*/
PetscErrorCode PEPGetErrorEstimate(PEP pep,PetscInt i,PetscReal *errest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(errest,3);
  PEPCheckSolved(pep,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCheck(i<pep->nconv,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see PEPGetConverged()");
  *errest = pep->errest[pep->perm[i]];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPComputeResidualNorm_Private - Computes the norm of the residual vector
   associated with an eigenpair.

   Input Parameters:
     kr,ki - eigenvalue
     xr,xi - eigenvector
     z     - array of 4 work vectors (z[2],z[3] not referenced in complex scalars)
*/
PetscErrorCode PEPComputeResidualNorm_Private(PEP pep,PetscScalar kr,PetscScalar ki,Vec xr,Vec xi,Vec *z,PetscReal *norm)
{
  Mat            *A=pep->A;
  PetscInt       i,nmat=pep->nmat;
  PetscScalar    t[20],*vals=t,*ivals=NULL;
  Vec            u,w;
#if !defined(PETSC_USE_COMPLEX)
  Vec            ui,wi;
  PetscReal      ni;
  PetscBool      imag;
  PetscScalar    it[20];
#endif

  PetscFunctionBegin;
  u = z[0]; w = z[1];
  PetscCall(VecSet(u,0.0));
#if !defined(PETSC_USE_COMPLEX)
  ui = z[2]; wi = z[3];
  ivals = it;
#endif
  if (nmat>20) {
    PetscCall(PetscMalloc1(nmat,&vals));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(PetscMalloc1(nmat,&ivals));
#endif
  }
  PetscCall(PEPEvaluateBasis(pep,kr,ki,vals,ivals));
#if !defined(PETSC_USE_COMPLEX)
  if (ki == 0 || PetscAbsScalar(ki) < PetscAbsScalar(kr*PETSC_MACHINE_EPSILON))
    imag = PETSC_FALSE;
  else {
    imag = PETSC_TRUE;
    PetscCall(VecSet(ui,0.0));
  }
#endif
  for (i=0;i<nmat;i++) {
    if (vals[i]!=0.0) {
      PetscCall(MatMult(A[i],xr,w));
      PetscCall(VecAXPY(u,vals[i],w));
    }
#if !defined(PETSC_USE_COMPLEX)
    if (imag) {
      if (ivals[i]!=0 || vals[i]!=0) {
        PetscCall(MatMult(A[i],xi,wi));
        if (vals[i]==0) PetscCall(MatMult(A[i],xr,w));
      }
      if (ivals[i]!=0) {
        PetscCall(VecAXPY(u,-ivals[i],wi));
        PetscCall(VecAXPY(ui,ivals[i],w));
      }
      if (vals[i]!=0) PetscCall(VecAXPY(ui,vals[i],wi));
    }
#endif
  }
  PetscCall(VecNorm(u,NORM_2,norm));
#if !defined(PETSC_USE_COMPLEX)
  if (imag) {
    PetscCall(VecNorm(ui,NORM_2,&ni));
    *norm = SlepcAbsEigenvalue(*norm,ni);
  }
#endif
  if (nmat>20) {
    PetscCall(PetscFree(vals));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(PetscFree(ivals));
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPComputeError - Computes the error (based on the residual norm) associated
   with the `i`-th computed eigenpair.

   Collective

   Input Parameters:
+  pep  - the polynomial eigensolver context
.  i    - the solution index
-  type - the type of error to compute, see `PEPErrorType`

   Output Parameter:
.  error - the error

   Note:
   The error can be computed in various ways, all of them based on the residual
   norm $\|P(\lambda)x\|_2$ where $(\lambda,x)$ is the approximate eigenpair.

   Level: beginner

.seealso: [](ch:pep), `PEPErrorType`, `PEPSolve()`, `PEPGetErrorEstimate()`
@*/
PetscErrorCode PEPComputeError(PEP pep,PetscInt i,PEPErrorType type,PetscReal *error)
{
  Vec            xr,xi,w[4];
  PetscScalar    kr,ki;
  PetscReal      t,z=0.0;
  PetscInt       j;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,i,2);
  PetscValidLogicalCollectiveEnum(pep,type,3);
  PetscAssertPointer(error,4);
  PEPCheckSolved(pep,1);

  /* allocate work vectors */
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PEPSetWorkVecs(pep,3));
  xi   = NULL;
  w[2] = NULL;
  w[3] = NULL;
#else
  PetscCall(PEPSetWorkVecs(pep,6));
  xi   = pep->work[3];
  w[2] = pep->work[4];
  w[3] = pep->work[5];
#endif
  xr   = pep->work[0];
  w[0] = pep->work[1];
  w[1] = pep->work[2];

  /* compute residual norms */
  PetscCall(PEPGetEigenpair(pep,i,&kr,&ki,xr,xi));
  PetscCall(PEPComputeResidualNorm_Private(pep,kr,ki,xr,xi,w,error));

  /* compute error */
  switch (type) {
    case PEP_ERROR_ABSOLUTE:
      break;
    case PEP_ERROR_RELATIVE:
      *error /= SlepcAbsEigenvalue(kr,ki);
      break;
    case PEP_ERROR_BACKWARD:
      /* initialization of matrix norms */
      if (!pep->nrma[pep->nmat-1]) {
        for (j=0;j<pep->nmat;j++) {
          PetscCall(MatHasOperation(pep->A[j],MATOP_NORM,&flg));
          PetscCheck(flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
          PetscCall(MatNorm(pep->A[j],NORM_INFINITY,&pep->nrma[j]));
        }
      }
      t = SlepcAbsEigenvalue(kr,ki);
      for (j=pep->nmat-1;j>=0;j--) {
        z = z*t+pep->nrma[j];
      }
      *error /= z;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid error type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
