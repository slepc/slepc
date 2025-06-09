/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   EPS routines related to the solution process
*/

#include <slepc/private/epsimpl.h>   /*I "slepceps.h" I*/
#include <slepc/private/bvimpl.h>
#include <petscdraw.h>

PetscErrorCode EPSComputeVectors(EPS eps)
{
  PetscFunctionBegin;
  EPSCheckSolved(eps,1);
  if (eps->state==EPS_STATE_SOLVED) PetscTryTypeMethod(eps,computevectors);
  eps->state = EPS_STATE_EIGENVECTORS;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeValues(EPS eps)
{
  PetscBool      injective,iscomp,isfilter;
  PetscInt       i,n,aux,nconv0;
  Mat            A,B=NULL,G,Z;

  PetscFunctionBegin;
  switch (eps->categ) {
    case EPS_CATEGORY_KRYLOV:
    case EPS_CATEGORY_OTHER:
      PetscCall(STIsInjective(eps->st,&injective));
      if (injective) {
        /* one-to-one mapping: backtransform eigenvalues */
        PetscUseTypeMethod(eps,backtransform);
      } else {
        /* compute eigenvalues from Rayleigh quotient */
        PetscCall(DSGetDimensions(eps->ds,&n,NULL,NULL,NULL));
        if (!n) break;
        PetscCall(EPSGetOperators(eps,&A,&B));
        PetscCall(BVSetActiveColumns(eps->V,0,n));
        PetscCall(DSGetCompact(eps->ds,&iscomp));
        PetscCall(DSSetCompact(eps->ds,PETSC_FALSE));
        PetscCall(DSGetMat(eps->ds,DS_MAT_A,&G));
        PetscCall(BVMatProject(eps->V,A,eps->V,G));
        PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&G));
        if (B) {
          PetscCall(DSGetMat(eps->ds,DS_MAT_B,&G));
          PetscCall(BVMatProject(eps->V,B,eps->V,G));
          PetscCall(DSRestoreMat(eps->ds,DS_MAT_B,&G));
        }
        PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
        PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
        PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));
        PetscCall(DSSetCompact(eps->ds,iscomp));
        if (eps->ishermitian && (!eps->isgeneralized || eps->ispositive)) { /* V = V * Z */
          PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));
          PetscCall(DSGetMat(eps->ds,DS_MAT_X,&Z));
          PetscCall(BVMultInPlace(eps->V,Z,0,n));
          PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&Z));
        }
        /* in case of STFILTER discard computed eigenvalues that lie outside the wanted interval */
        PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilter));
        if (isfilter) {
          nconv0 = eps->nconv;
          for (i=0;i<eps->nconv;i++) {
            if (PetscRealPart(eps->eigr[eps->perm[i]])<eps->inta || PetscRealPart(eps->eigr[eps->perm[i]])>eps->intb) {
              eps->nconv--;
              if (i<eps->nconv) { SlepcSwap(eps->perm[i],eps->perm[eps->nconv],aux); i--; }
            }
          }
          if (nconv0>eps->nconv) PetscCall(PetscInfo(eps,"Discarded %" PetscInt_FMT " computed eigenvalues lying outside the interval\n",nconv0-eps->nconv));
        }
      }
      break;
    case EPS_CATEGORY_PRECOND:
    case EPS_CATEGORY_CONTOUR:
      /* eigenvalues already available as an output of the solver */
      break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSolve - Solves the eigensystem.

   Collective

   Input Parameter:
.  eps - eigensolver context obtained from EPSCreate()

   Options Database Keys:
+  -eps_view - print information about the solver used
.  -eps_view_mat0 - view the first matrix (A)
.  -eps_view_mat1 - view the second matrix (B)
.  -eps_view_vectors - view the computed eigenvectors
.  -eps_view_values - view the computed eigenvalues
.  -eps_converged_reason - print reason for convergence, and number of iterations
.  -eps_error_absolute - print absolute errors of each eigenpair
.  -eps_error_relative - print relative errors of each eigenpair
-  -eps_error_backward - print backward errors of each eigenpair

   Notes:
   All the command-line options listed above admit an optional argument specifying
   the viewer type and options. For instance, use '-eps_view_mat0 binary:amatrix.bin'
   to save the A matrix to a binary file, '-eps_view_values draw' to draw the computed
   eigenvalues graphically, or '-eps_error_relative :myerr.m:ascii_matlab' to save
   the errors in a file that can be executed in Matlab.

   Level: beginner

.seealso: EPSCreate(), EPSSetUp(), EPSDestroy(), EPSSetTolerances()
@*/
PetscErrorCode EPSSolve(EPS eps)
{
  PetscInt       i;
  PetscBool      hasname;
  STMatMode      matmode;
  Mat            A,B;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (eps->state>=EPS_STATE_SOLVED) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscLogEventBegin(EPS_Solve,eps,0,0,0));

  /* Call setup */
  PetscCall(EPSSetUp(eps));
  eps->nconv = 0;
  eps->its   = 0;
  for (i=0;i<eps->ncv;i++) {
    eps->eigr[i]   = 0.0;
    eps->eigi[i]   = 0.0;
    eps->errest[i] = 0.0;
    eps->perm[i]   = i;
  }
  PetscCall(EPSViewFromOptions(eps,NULL,"-eps_view_pre"));
  PetscCall(RGViewFromOptions(eps->rg,NULL,"-rg_view"));

  /* Call solver */
  PetscUseTypeMethod(eps,solve);
  PetscCheck(eps->reason,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");
  eps->state = EPS_STATE_SOLVED;

  /* Only the first nconv columns contain useful information (except in CISS) */
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
  if (eps->twosided) PetscCall(BVSetActiveColumns(eps->W,0,eps->nconv));

  /* If inplace, purify eigenvectors before reverting operator */
  PetscCall(STGetMatMode(eps->st,&matmode));
  if (matmode == ST_MATMODE_INPLACE && eps->ispositive) PetscCall(EPSComputeVectors(eps));
  PetscCall(STPostSolve(eps->st));

  /* Map eigenvalues back to the original problem if appropriate */
  PetscCall(EPSComputeValues(eps));

#if !defined(PETSC_USE_COMPLEX)
  /* Reorder conjugate eigenvalues (positive imaginary first) */
  for (i=0;i<eps->nconv-1;i++) {
    if (eps->eigi[i] != 0) {
      if (eps->eigi[i] < 0) {
        eps->eigi[i] = -eps->eigi[i];
        eps->eigi[i+1] = -eps->eigi[i+1];
        /* the next correction only works with eigenvectors */
        PetscCall(EPSComputeVectors(eps));
        PetscCall(BVScaleColumn(eps->V,i+1,-1.0));
        if (eps->W) PetscCall(BVScaleColumn(eps->W,i+1,-1.0));
      }
      i++;
    }
  }
#endif

  /* Sort eigenvalues according to eps->which parameter */
  PetscCall(SlepcSortEigenvalues(eps->sc,eps->nconv,eps->eigr,eps->eigi,eps->perm));
  PetscCall(PetscLogEventEnd(EPS_Solve,eps,0,0,0));

  /* Various viewers */
  PetscCall(EPSViewFromOptions(eps,NULL,"-eps_view"));
  PetscCall(EPSConvergedReasonViewFromOptions(eps));
  PetscCall(EPSErrorViewFromOptions(eps));
  PetscCall(EPSValuesViewFromOptions(eps));
  PetscCall(EPSVectorsViewFromOptions(eps));

  PetscCall(PetscOptionsHasName(((PetscObject)eps)->options,((PetscObject)eps)->prefix,"-eps_view_mat0",&hasname));
  if (hasname) {
    PetscCall(EPSGetOperators(eps,&A,NULL));
    PetscCall(MatViewFromOptions(A,(PetscObject)eps,"-eps_view_mat0"));
  }
  if (eps->isgeneralized) {
    PetscCall(PetscOptionsHasName(((PetscObject)eps)->options,((PetscObject)eps)->prefix,"-eps_view_mat1",&hasname));
    if (hasname) {
      PetscCall(EPSGetOperators(eps,NULL,&B));
      PetscCall(MatViewFromOptions(B,(PetscObject)eps,"-eps_view_mat1"));
    }
  }

  /* Remove deflation and initial subspaces */
  if (eps->nds) {
    PetscCall(BVSetNumConstraints(eps->V,0));
    eps->nds = 0;
  }
  eps->nini = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetIterationNumber - Gets the current iteration number. If the
   call to EPSSolve() is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  its - number of iterations

   Note:
   During the i-th iteration this call returns i-1. If EPSSolve() is
   complete, then parameter "its" contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call EPSGetConvergedReason() to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: EPSGetConvergedReason(), EPSSetTolerances()
@*/
PetscErrorCode EPSGetIterationNumber(EPS eps,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(its,2);
  *its = eps->its;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetConverged - Gets the number of converged eigenpairs.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  nconv - number of converged eigenpairs

   Note:
   This function should be called after EPSSolve() has finished.

   Level: beginner

.seealso: EPSSetDimensions(), EPSSolve(), EPSGetEigenpair()
@*/
PetscErrorCode EPSGetConverged(EPS eps,PetscInt *nconv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(nconv,2);
  EPSCheckSolved(eps,1);
  PetscCall(EPS_GetActualConverged(eps,nconv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetConvergedReason - Gets the reason why the EPSSolve() iteration was
   stopped.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged

   Options Database Key:
.  -eps_converged_reason - print the reason to a viewer

   Notes:
   Possible values for reason are
+  EPS_CONVERGED_TOL - converged up to tolerance
.  EPS_CONVERGED_USER - converged due to a user-defined condition
.  EPS_DIVERGED_ITS - required more than max_it iterations to reach convergence
.  EPS_DIVERGED_BREAKDOWN - generic breakdown in method
-  EPS_DIVERGED_SYMMETRY_LOST - pseudo-Lanczos was not able to keep symmetry

   Can only be called after the call to EPSSolve() is complete.

   Level: intermediate

.seealso: EPSSetTolerances(), EPSSolve(), EPSConvergedReason
@*/
PetscErrorCode EPSGetConvergedReason(EPS eps,EPSConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(reason,2);
  EPSCheckSolved(eps,1);
  *reason = eps->reason;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetInvariantSubspace - Gets an orthonormal basis of the computed invariant
   subspace.

   Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  v - an array of vectors

   Notes:
   This function should be called after EPSSolve() has finished.

   The user should provide in v an array of nconv vectors, where nconv is
   the value returned by EPSGetConverged().

   The first k vectors returned in v span an invariant subspace associated
   with the first k computed eigenvalues (note that this is not true if the
   k-th eigenvalue is complex and matrix A is real; in this case the first
   k+1 vectors should be used). An invariant subspace X of A satisfies Ax
   in X for all x in X (a similar definition applies for generalized
   eigenproblems).

   Level: intermediate

.seealso: EPSGetEigenpair(), EPSGetConverged(), EPSSolve()
@*/
PetscErrorCode EPSGetInvariantSubspace(EPS eps,Vec v[])
{
  PetscInt       i;
  BV             V=eps->V;
  Vec            w;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(v,2);
  PetscValidHeaderSpecific(*v,VEC_CLASSID,2);
  EPSCheckSolved(eps,1);
  PetscCheck(eps->ishermitian || eps->state!=EPS_STATE_EIGENVECTORS,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"EPSGetInvariantSubspace must be called before EPSGetEigenpair,EPSGetEigenvector or EPSComputeError");
  if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
    PetscCall(BVDuplicateResize(eps->V,eps->nconv,&V));
    PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
    PetscCall(BVCopy(eps->V,V));
    for (i=0;i<eps->nconv;i++) {
      PetscCall(BVGetColumn(V,i,&w));
      PetscCall(VecPointwiseDivide(w,w,eps->D));
      PetscCall(BVRestoreColumn(V,i,&w));
    }
    PetscCall(BVOrthogonalize(V,NULL));
  }
  for (i=0;i<eps->nconv;i++) PetscCall(BVCopyVec(V,i,v[i]));
  if (eps->balance!=EPS_BALANCE_NONE && eps->D) PetscCall(BVDestroy(&V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetEigenpair - Gets the i-th solution of the eigenproblem as computed by
   EPSSolve(). The solution consists in both the eigenvalue and the eigenvector.

   Collective

   Input Parameters:
+  eps - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  eigr - real part of eigenvalue
.  eigi - imaginary part of eigenvalue
.  Vr   - real part of eigenvector
-  Vi   - imaginary part of eigenvector

   Notes:
   It is allowed to pass NULL for Vr and Vi, if the eigenvector is not
   required. Otherwise, the caller must provide valid Vec objects, i.e.,
   they must be created by the calling program with e.g. MatCreateVecs().

   If the eigenvalue is real, then eigi and Vi are set to zero. If PETSc is
   configured with complex scalars the eigenvalue is stored
   directly in eigr (eigi is set to zero) and the eigenvector in Vr (Vi is
   set to zero). In both cases, the user can pass NULL in eigi and Vi.

   The index i should be a value between 0 and nconv-1 (see EPSGetConverged()).
   Eigenpairs are indexed according to the ordering criterion established
   with EPSSetWhichEigenpairs().

   The 2-norm of the eigenvector is one unless the problem is generalized
   Hermitian. In this case the eigenvector is normalized with respect to the
   norm defined by the B matrix.

   Level: beginner

.seealso: EPSGetEigenvalue(), EPSGetEigenvector(), EPSGetLeftEigenvector(), EPSSolve(),
          EPSGetConverged(), EPSSetWhichEigenpairs(), EPSGetInvariantSubspace()
@*/
PetscErrorCode EPSGetEigenpair(EPS eps,PetscInt i,PetscScalar *eigr,PetscScalar *eigi,Vec Vr,Vec Vi)
{
  PetscInt nconv;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);
  EPSCheckSolved(eps,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCall(EPS_GetActualConverged(eps,&nconv));
  PetscCheck(i<nconv,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see EPSGetConverged()");
  PetscCall(EPSGetEigenvalue(eps,i,eigr,eigi));
  if (Vr || Vi) PetscCall(EPSGetEigenvector(eps,i,Vr,Vi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetEigenvalue - Gets the i-th eigenvalue as computed by EPSSolve().

   Not Collective

   Input Parameters:
+  eps - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  eigr - real part of eigenvalue
-  eigi - imaginary part of eigenvalue

   Notes:
   If the eigenvalue is real, then eigi is set to zero. If PETSc is
   configured with complex scalars the eigenvalue is stored
   directly in eigr (eigi is set to zero).

   The index i should be a value between 0 and nconv-1 (see EPSGetConverged()).
   Eigenpairs are indexed according to the ordering criterion established
   with EPSSetWhichEigenpairs().

   Level: beginner

.seealso: EPSSolve(), EPSGetConverged(), EPSSetWhichEigenpairs(), EPSGetEigenpair()
@*/
PetscErrorCode EPSGetEigenvalue(EPS eps,PetscInt i,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscInt k,nconv;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  EPSCheckSolved(eps,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCall(EPS_GetActualConverged(eps,&nconv));
  PetscCheck(i<nconv,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see EPSGetConverged()");
  if (nconv==eps->nconv) {
    k = eps->perm[i];
#if defined(PETSC_USE_COMPLEX)
    if (eigr) *eigr = eps->eigr[k];
    if (eigi) *eigi = 0;
#else
    if (eigr) *eigr = eps->eigr[k];
    if (eigi) *eigi = eps->eigi[k];
#endif
  } else {
    PetscCheck(eps->problem_type==EPS_BSE,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Problem type should be BSE");
    /* BSE problem, even index is +lambda, odd index is -lambda */
    k = eps->perm[i/2];
#if defined(PETSC_USE_COMPLEX)
    if (eigr) *eigr = (i%2)? -eps->eigr[k]: eps->eigr[k];
    if (eigi) *eigi = 0;
#else
    if (eigr) *eigr = (i%2)? -eps->eigr[k]: eps->eigr[k];
    if (eigi) *eigi = eps->eigi[k];
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetEigenvector - Gets the i-th right eigenvector as computed by EPSSolve().

   Collective

   Input Parameters:
+  eps - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  Vr   - real part of eigenvector
-  Vi   - imaginary part of eigenvector

   Notes:
   The caller must provide valid Vec objects, i.e., they must be created
   by the calling program with e.g. MatCreateVecs().

   If the corresponding eigenvalue is real, then Vi is set to zero. If PETSc is
   configured with complex scalars the eigenvector is stored
   directly in Vr (Vi is set to zero). In any case, the user can pass NULL in Vr
   or Vi if one of them is not required.

   The index i should be a value between 0 and nconv-1 (see EPSGetConverged()).
   Eigenpairs are indexed according to the ordering criterion established
   with EPSSetWhichEigenpairs().

   The 2-norm of the eigenvector is one unless the problem is generalized
   Hermitian. In this case the eigenvector is normalized with respect to the
   norm defined by the B matrix.

   Level: beginner

.seealso: EPSSolve(), EPSGetConverged(), EPSSetWhichEigenpairs(), EPSGetEigenpair(), EPSGetLeftEigenvector()
@*/
PetscErrorCode EPSGetEigenvector(EPS eps,PetscInt i,Vec Vr,Vec Vi)
{
  PetscInt nconv;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);
  if (Vr) { PetscValidHeaderSpecific(Vr,VEC_CLASSID,3); PetscCheckSameComm(eps,1,Vr,3); }
  if (Vi) { PetscValidHeaderSpecific(Vi,VEC_CLASSID,4); PetscCheckSameComm(eps,1,Vi,4); }
  EPSCheckSolved(eps,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCall(EPS_GetActualConverged(eps,&nconv));
  PetscCheck(i<nconv,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see EPSGetConverged()");
  PetscCall(EPSComputeVectors(eps));
  PetscCall(EPS_GetEigenvector(eps,eps->V,i,Vr,Vi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetLeftEigenvector - Gets the i-th left eigenvector as computed by EPSSolve().

   Collective

   Input Parameters:
+  eps - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  Wr   - real part of left eigenvector
-  Wi   - imaginary part of left eigenvector

   Notes:
   The caller must provide valid Vec objects, i.e., they must be created
   by the calling program with e.g. MatCreateVecs().

   If the corresponding eigenvalue is real, then Wi is set to zero. If PETSc is
   configured with complex scalars the eigenvector is stored directly in Wr
   (Wi is set to zero). In any case, the user can pass NULL in Wr or Wi if
   one of them is not required.

   The index i should be a value between 0 and nconv-1 (see EPSGetConverged()).
   Eigensolutions are indexed according to the ordering criterion established
   with EPSSetWhichEigenpairs().

   Left eigenvectors are available only if the twosided flag was set, see
   EPSSetTwoSided().

   Level: intermediate

.seealso: EPSGetEigenvector(), EPSGetConverged(), EPSSetWhichEigenpairs(), EPSSetTwoSided()
@*/
PetscErrorCode EPSGetLeftEigenvector(EPS eps,PetscInt i,Vec Wr,Vec Wi)
{
  PetscInt    nconv;
  PetscBool   trivial;
  Mat         H;
  IS          is[2];
  Vec         v;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);
  if (Wr) { PetscValidHeaderSpecific(Wr,VEC_CLASSID,3); PetscCheckSameComm(eps,1,Wr,3); }
  if (Wi) { PetscValidHeaderSpecific(Wi,VEC_CLASSID,4); PetscCheckSameComm(eps,1,Wi,4); }
  EPSCheckSolved(eps,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCall(EPS_GetActualConverged(eps,&nconv));
  PetscCheck(i<nconv,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see EPSGetConverged()");

  trivial = (eps->problem_type==EPS_HEP || eps->problem_type==EPS_GHEP || eps->problem_type==EPS_BSE)? PETSC_TRUE: PETSC_FALSE;
  if (!trivial) PetscCheck(eps->twosided,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must request left vectors with EPSSetTwoSided");

  PetscCall(EPSComputeVectors(eps));
  if (trivial) {
    PetscCall(EPS_GetEigenvector(eps,eps->V,i,Wr,Wi));
    if (eps->problem_type==EPS_BSE) {   /* change sign of bottom part of the vector */
      PetscCall(STGetMatrix(eps->st,0,&H));
      PetscCall(MatNestGetISs(H,is,NULL));
      if (Wr) {
        PetscCall(VecGetSubVector(Wr,is[1],&v));
        PetscCall(VecScale(v,-1.0));
        PetscCall(VecRestoreSubVector(Wr,is[1],&v));
      }
#if !defined(PETSC_USE_COMPLEX)
      if (Wi) {
        PetscCall(VecGetSubVector(Wi,is[1],&v));
        PetscCall(VecScale(v,-1.0));
        PetscCall(VecRestoreSubVector(Wi,is[1],&v));
      }
#endif
    }
  } else {
    PetscCall(EPS_GetEigenvector(eps,eps->W,i,Wr,Wi));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSGetErrorEstimate - Returns the error estimate associated to the i-th
   computed eigenpair.

   Not Collective

   Input Parameters:
+  eps - eigensolver context
-  i   - index of eigenpair

   Output Parameter:
.  errest - the error estimate

   Notes:
   This is the error estimate used internally by the eigensolver. The actual
   error bound can be computed with EPSComputeError(). See also the users
   manual for details.

   Level: advanced

.seealso: EPSComputeError()
@*/
PetscErrorCode EPSGetErrorEstimate(EPS eps,PetscInt i,PetscReal *errest)
{
  PetscInt nconv;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(errest,3);
  EPSCheckSolved(eps,1);
  PetscCheck(i>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index cannot be negative");
  PetscCall(EPS_GetActualConverged(eps,&nconv));
  PetscCheck(i<nconv,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The index can be nconv-1 at most, see EPSGetConverged()");
  if (nconv==eps->nconv) {
    *errest = eps->errest[eps->perm[i]];
  } else {
    PetscCheck(eps->problem_type==EPS_BSE,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Problem type should be BSE");
    /* BSE problem, even index is +lambda, odd index is -lambda, assume both have same error */
    *errest = eps->errest[eps->perm[i/2]];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSComputeResidualNorm_Private - Computes the norm of the residual vector
   associated with an eigenpair.

   Input Parameters:
     trans - whether A' must be used instead of A
     kr,ki - eigenvalue
     xr,xi - eigenvector
     z     - three work vectors (the second one not referenced in complex scalars)
*/
PetscErrorCode EPSComputeResidualNorm_Private(EPS eps,PetscBool trans,PetscScalar kr,PetscScalar ki,Vec xr,Vec xi,Vec *z,PetscReal *norm)
{
  PetscInt       nmat;
  Mat            A,B;
  Vec            u,w;
  PetscScalar    alpha;
#if !defined(PETSC_USE_COMPLEX)
  Vec            v;
  PetscReal      ni,nr;
#endif
  PetscErrorCode (*matmult)(Mat,Vec,Vec) = trans? MatMultHermitianTranspose: MatMult;

  PetscFunctionBegin;
  u = z[0]; w = z[2];
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  if (nmat>1) PetscCall(STGetMatrix(eps->st,1,&B));

#if !defined(PETSC_USE_COMPLEX)
  v = z[1];
  if (ki == 0 || PetscAbsScalar(ki) < PetscAbsScalar(kr*PETSC_MACHINE_EPSILON)) {
#endif
    PetscCall((*matmult)(A,xr,u));                          /* u=A*x */
    if (PetscAbsScalar(kr) > PETSC_MACHINE_EPSILON) {
      if (nmat>1) PetscCall((*matmult)(B,xr,w));
      else PetscCall(VecCopy(xr,w));                        /* w=B*x */
      alpha = trans? -PetscConj(kr): -kr;
      PetscCall(VecAXPY(u,alpha,w));                        /* u=A*x-k*B*x */
    }
    PetscCall(VecNorm(u,NORM_2,norm));
#if !defined(PETSC_USE_COMPLEX)
  } else {
    PetscCall((*matmult)(A,xr,u));                          /* u=A*xr */
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      if (nmat>1) PetscCall((*matmult)(B,xr,v));
      else PetscCall(VecCopy(xr,v));                        /* v=B*xr */
      PetscCall(VecAXPY(u,-kr,v));                          /* u=A*xr-kr*B*xr */
      if (nmat>1) PetscCall((*matmult)(B,xi,w));
      else PetscCall(VecCopy(xi,w));                        /* w=B*xi */
      PetscCall(VecAXPY(u,trans?-ki:ki,w));                 /* u=A*xr-kr*B*xr+ki*B*xi */
    }
    PetscCall(VecNorm(u,NORM_2,&nr));
    PetscCall((*matmult)(A,xi,u));                          /* u=A*xi */
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      PetscCall(VecAXPY(u,-kr,w));                          /* u=A*xi-kr*B*xi */
      PetscCall(VecAXPY(u,trans?ki:-ki,v));                 /* u=A*xi-kr*B*xi-ki*B*xr */
    }
    PetscCall(VecNorm(u,NORM_2,&ni));
    *norm = SlepcAbsEigenvalue(nr,ni);
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSComputeError - Computes the error (based on the residual norm) associated
   with the i-th computed eigenpair.

   Collective

   Input Parameters:
+  eps  - the eigensolver context
.  i    - the solution index
-  type - the type of error to compute

   Output Parameter:
.  error - the error

   Notes:
   The error can be computed in various ways, all of them based on the residual
   norm ||Ax-kBx||_2 where k is the eigenvalue and x is the eigenvector.

   Level: beginner

.seealso: EPSErrorType, EPSSolve(), EPSGetErrorEstimate()
@*/
PetscErrorCode EPSComputeError(EPS eps,PetscInt i,EPSErrorType type,PetscReal *error)
{
  Mat            A,B;
  Vec            xr,xi,w[3];
  PetscReal      t,vecnorm=1.0,errorl;
  PetscScalar    kr,ki;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);
  PetscValidLogicalCollectiveEnum(eps,type,3);
  PetscAssertPointer(error,4);
  EPSCheckSolved(eps,1);

  /* allocate work vectors */
#if defined(PETSC_USE_COMPLEX)
  PetscCall(EPSSetWorkVecs(eps,3));
  xi   = NULL;
  w[1] = NULL;
#else
  PetscCall(EPSSetWorkVecs(eps,5));
  xi   = eps->work[3];
  w[1] = eps->work[4];
#endif
  xr   = eps->work[0];
  w[0] = eps->work[1];
  w[2] = eps->work[2];

  /* compute residual norm */
  PetscCall(EPSGetEigenpair(eps,i,&kr,&ki,xr,xi));
  PetscCall(EPSComputeResidualNorm_Private(eps,PETSC_FALSE,kr,ki,xr,xi,w,error));

  /* compute 2-norm of eigenvector */
  if (eps->problem_type==EPS_GHEP) PetscCall(VecNorm(xr,NORM_2,&vecnorm));

  /* if two-sided, compute left residual norm and take the maximum */
  if (eps->twosided) {
    PetscCall(EPSGetLeftEigenvector(eps,i,xr,xi));
    PetscCall(EPSComputeResidualNorm_Private(eps,PETSC_TRUE,kr,ki,xr,xi,w,&errorl));
    *error = PetscMax(*error,errorl);
  }

  /* compute error */
  switch (type) {
    case EPS_ERROR_ABSOLUTE:
      break;
    case EPS_ERROR_RELATIVE:
      *error /= SlepcAbsEigenvalue(kr,ki)*vecnorm;
      break;
    case EPS_ERROR_BACKWARD:
      /* initialization of matrix norms */
      if (!eps->nrma) {
        PetscCall(STGetMatrix(eps->st,0,&A));
        PetscCall(MatHasOperation(A,MATOP_NORM,&flg));
        PetscCheck(flg,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
        PetscCall(MatNorm(A,NORM_INFINITY,&eps->nrma));
      }
      if (eps->isgeneralized) {
        if (!eps->nrmb) {
          PetscCall(STGetMatrix(eps->st,1,&B));
          PetscCall(MatHasOperation(B,MATOP_NORM,&flg));
          PetscCheck(flg,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
          PetscCall(MatNorm(B,NORM_INFINITY,&eps->nrmb));
        }
      } else eps->nrmb = 1.0;
      t = SlepcAbsEigenvalue(kr,ki);
      *error /= (eps->nrma+t*eps->nrmb)*vecnorm;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid error type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSGetStartVector - Generate a suitable vector to be used as the starting vector
   for the recurrence that builds the right subspace.

   Collective

   Input Parameters:
+  eps - the eigensolver context
-  i   - iteration number

   Output Parameters:
.  breakdown - flag indicating that a breakdown has occurred

   Notes:
   The start vector is computed from another vector: for the first step (i=0),
   the first initial vector is used (see EPSSetInitialSpace()); otherwise a random
   vector is created. Then this vector is forced to be in the range of OP (only
   for generalized definite problems) and orthonormalized with respect to all
   V-vectors up to i-1. The resulting vector is placed in V[i].

   The flag breakdown is set to true if either i=0 and the vector belongs to the
   deflation space, or i>0 and the vector is linearly dependent with respect
   to the V-vectors.
*/
PetscErrorCode EPSGetStartVector(EPS eps,PetscInt i,PetscBool *breakdown)
{
  PetscReal      norm;
  PetscBool      lindep;
  Vec            w,z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);

  /* For the first step, use the first initial vector, otherwise a random one */
  if (i>0 || eps->nini==0) PetscCall(BVSetRandomColumn(eps->V,i));

  /* Force the vector to be in the range of OP for generalized problems with B-inner product */
  if (eps->ispositive || (eps->isgeneralized && eps->ishermitian)) {
    PetscCall(BVCreateVec(eps->V,&w));
    PetscCall(BVCopyVec(eps->V,i,w));
    PetscCall(BVGetColumn(eps->V,i,&z));
    PetscCall(STApply(eps->st,w,z));
    PetscCall(BVRestoreColumn(eps->V,i,&z));
    PetscCall(VecDestroy(&w));
  }

  /* Orthonormalize the vector with respect to previous vectors */
  PetscCall(BVOrthogonalizeColumn(eps->V,i,NULL,&norm,&lindep));
  if (breakdown) *breakdown = lindep;
  else if (lindep || norm == 0.0) {
    PetscCheck(i,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Initial vector is zero or belongs to the deflation space");
    PetscCheck(!i,PetscObjectComm((PetscObject)eps),PETSC_ERR_CONV_FAILED,"Unable to generate more start vectors");
  }
  PetscCall(BVScaleColumn(eps->V,i,1.0/norm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSGetLeftStartVector - Generate a suitable vector to be used as the left starting
   vector for the recurrence that builds the left subspace. See EPSGetStartVector().
*/
PetscErrorCode EPSGetLeftStartVector(EPS eps,PetscInt i,PetscBool *breakdown)
{
  PetscReal      norm;
  PetscBool      lindep;
  Vec            w,z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);

  /* For the first step, use the first initial vector, otherwise a random one */
  if (i>0 || eps->ninil==0) PetscCall(BVSetRandomColumn(eps->W,i));

  /* Force the vector to be in the range of OP' for generalized problems with B-inner product */
  if (eps->ispositive || (eps->isgeneralized && eps->ishermitian)) {
    PetscCall(BVCreateVec(eps->W,&w));
    PetscCall(BVCopyVec(eps->W,i,w));
    PetscCall(BVGetColumn(eps->W,i,&z));
    PetscCall(STApplyHermitianTranspose(eps->st,w,z));
    PetscCall(BVRestoreColumn(eps->W,i,&z));
    PetscCall(VecDestroy(&w));
  }

  /* Orthonormalize the vector with respect to previous vectors */
  PetscCall(BVOrthogonalizeColumn(eps->W,i,NULL,&norm,&lindep));
  if (breakdown) *breakdown = lindep;
  else if (lindep || norm == 0.0) {
    PetscCheck(i,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Left initial vector is zero");
    PetscCheck(!i,PetscObjectComm((PetscObject)eps),PETSC_ERR_CONV_FAILED,"Unable to generate more left start vectors");
  }
  PetscCall(BVScaleColumn(eps->W,i,1.0/norm));
  PetscFunctionReturn(PETSC_SUCCESS);
}
