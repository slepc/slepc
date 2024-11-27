/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file contains some simple default routines for common operations
*/

#include <slepc/private/epsimpl.h>   /*I "slepceps.h" I*/
#include <slepcvec.h>

PetscErrorCode EPSBackTransform_Default(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(STBackTransform(eps->st,eps->nconv,eps->eigr,eps->eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSComputeVectors_Hermitian - Copies the Lanczos vectors as eigenvectors
  using purification for generalized eigenproblems.
 */
PetscErrorCode EPSComputeVectors_Hermitian(EPS eps)
{
  PetscBool      iscayley,indef;
  Mat            B,C;

  PetscFunctionBegin;
  if (eps->purify) {
    PetscCall(EPS_Purify(eps,eps->nconv));
    PetscCall(BVNormalize(eps->V,NULL));
  } else {
    /* In the case of Cayley transform, eigenvectors need to be B-normalized */
    PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STCAYLEY,&iscayley));
    if (iscayley && eps->isgeneralized) {
      PetscCall(STGetMatrix(eps->st,1,&B));
      PetscCall(BVGetMatrix(eps->V,&C,&indef));
      PetscCheck(!indef,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"The inner product should not be indefinite");
      PetscCall(BVSetMatrix(eps->V,B,PETSC_FALSE));
      PetscCall(BVNormalize(eps->V,NULL));
      PetscCall(BVSetMatrix(eps->V,C,PETSC_FALSE));  /* restore original matrix */
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSComputeVectors_Indefinite - similar to the Schur version but
  for indefinite problems
 */
PetscErrorCode EPSComputeVectors_Indefinite(EPS eps)
{
  PetscInt       n;
  Mat            X;

  PetscFunctionBegin;
  PetscCall(DSGetDimensions(eps->ds,&n,NULL,NULL,NULL));
  PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));
  PetscCall(DSGetMat(eps->ds,DS_MAT_X,&X));
  PetscCall(BVMultInPlace(eps->V,X,0,n));
  PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&X));

  /* purification */
  if (eps->purify) PetscCall(EPS_Purify(eps,eps->nconv));

  /* normalization */
  PetscCall(BVNormalize(eps->V,eps->eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSComputeVectors_Twosided - Adjust left eigenvectors in generalized problems: y = B^-* y.
 */
PetscErrorCode EPSComputeVectors_Twosided(EPS eps)
{
  PetscInt       i;
  Vec            w,y;

  PetscFunctionBegin;
  if (!eps->twosided || !eps->isgeneralized) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(EPSSetWorkVecs(eps,1));
  w = eps->work[0];
  for (i=0;i<eps->nconv;i++) {
    PetscCall(BVCopyVec(eps->W,i,w));
    PetscCall(BVGetColumn(eps->W,i,&y));
    PetscCall(STMatSolveHermitianTranspose(eps->st,w,y));
    PetscCall(BVRestoreColumn(eps->W,i,&y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSComputeVectors_Schur - Compute eigenvectors from the vectors
  provided by the eigensolver. This version is intended for solvers
  that provide Schur vectors. Given the partial Schur decomposition
  OP*V=V*T, the following steps are performed:
      1) compute eigenvectors of T: T*Z=Z*D
      2) compute eigenvectors of OP: X=V*Z
 */
PetscErrorCode EPSComputeVectors_Schur(EPS eps)
{
  PetscInt       i;
  Mat            Z;
  Vec            z;

  PetscFunctionBegin;
  if (eps->ishermitian) {
    if (eps->isgeneralized && !eps->ispositive) PetscCall(EPSComputeVectors_Indefinite(eps));
    else PetscCall(EPSComputeVectors_Hermitian(eps));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* right eigenvectors */
  PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));

  /* V = V * Z */
  PetscCall(DSGetMat(eps->ds,DS_MAT_X,&Z));
  PetscCall(BVMultInPlace(eps->V,Z,0,eps->nconv));
  PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&Z));

  /* Purify eigenvectors */
  if (eps->purify) PetscCall(EPS_Purify(eps,eps->nconv));

  /* Fix eigenvectors if balancing was used */
  if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
    for (i=0;i<eps->nconv;i++) {
      PetscCall(BVGetColumn(eps->V,i,&z));
      PetscCall(VecPointwiseDivide(z,z,eps->D));
      PetscCall(BVRestoreColumn(eps->V,i,&z));
    }
  }

  /* normalize eigenvectors (when using purification or balancing) */
  if (eps->purify || (eps->balance!=EPS_BALANCE_NONE && eps->D)) PetscCall(BVNormalize(eps->V,eps->eigi));

  /* left eigenvectors */
  if (eps->twosided) {
    PetscCall(DSVectors(eps->ds,DS_MAT_Y,NULL,NULL));
    /* W = W * Z */
    PetscCall(DSGetMat(eps->ds,DS_MAT_Y,&Z));
    PetscCall(BVMultInPlace(eps->W,Z,0,eps->nconv));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Y,&Z));
    /* Fix left eigenvectors if balancing was used */
    if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
      for (i=0;i<eps->nconv;i++) {
        PetscCall(BVGetColumn(eps->W,i,&z));
        PetscCall(VecPointwiseMult(z,z,eps->D));
        PetscCall(BVRestoreColumn(eps->W,i,&z));
      }
    }
    PetscCall(EPSComputeVectors_Twosided(eps));
    /* normalize */
    PetscCall(BVNormalize(eps->W,eps->eigi));
#if !defined(PETSC_USE_COMPLEX)
    for (i=0;i<eps->nconv-1;i++) {
      if (eps->eigi[i] != 0.0) {
        if (eps->eigi[i] > 0.0) PetscCall(BVScaleColumn(eps->W,i+1,-1.0));
        i++;
      }
    }
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSSetWorkVecs - Sets a number of work vectors into an EPS object.

   Collective

   Input Parameters:
+  eps - eigensolver context
-  nw  - number of work vectors to allocate

   Developer Notes:
   This is SLEPC_EXTERN because it may be required by user plugin EPS
   implementations.

   Level: developer

.seealso: EPSSetUp()
@*/
PetscErrorCode EPSSetWorkVecs(EPS eps,PetscInt nw)
{
  Vec            t;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nw,2);
  PetscCheck(nw>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"nw must be > 0: nw = %" PetscInt_FMT,nw);
  if (eps->nwork < nw) {
    PetscCall(VecDestroyVecs(eps->nwork,&eps->work));
    eps->nwork = nw;
    PetscCall(BVGetColumn(eps->V,0,&t));
    PetscCall(VecDuplicateVecs(t,nw,&eps->work));
    PetscCall(BVRestoreColumn(eps->V,0,&t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSSetWhichEigenpairs_Default - Sets the default value for which,
  depending on the ST.
 */
PetscErrorCode EPSSetWhichEigenpairs_Default(EPS eps)
{
  PetscBool      target;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompareAny((PetscObject)eps->st,&target,STSINVERT,STCAYLEY,""));
  if (target) eps->which = EPS_TARGET_MAGNITUDE;
  else eps->which = EPS_LARGEST_MAGNITUDE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSConvergedRelative - Checks convergence relative to the eigenvalue.
*/
PetscErrorCode EPSConvergedRelative(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal w;

  PetscFunctionBegin;
  w = SlepcAbsEigenvalue(eigr,eigi);
  *errest = (w!=0.0)? res/w: PETSC_MAX_REAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSConvergedAbsolute - Checks convergence absolutely.
*/
PetscErrorCode EPSConvergedAbsolute(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = res;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSConvergedNorm - Checks convergence relative to the eigenvalue and
  the matrix norms.
*/
PetscErrorCode EPSConvergedNorm(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal w;

  PetscFunctionBegin;
  w = SlepcAbsEigenvalue(eigr,eigi);
  *errest = res / (eps->nrma + w*eps->nrmb);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSStoppingBasic - Default routine to determine whether the outer eigensolver
   iteration must be stopped.

   Collective

   Input Parameters:
+  eps    - eigensolver context obtained from EPSCreate()
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged eigenpairs
.  nev    - number of requested eigenpairs
-  ctx    - context (not used here)

   Output Parameter:
.  reason - result of the stopping test

   Notes:
   A positive value of reason indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, reason must be set to EPS_CONVERGED_ITERATING
   (zero).

   EPSStoppingBasic() will stop if all requested eigenvalues are converged, or if
   the maximum number of iterations has been reached.

   Use EPSSetStoppingTest() to provide your own test instead of using this one.

   Level: advanced

.seealso: EPSSetStoppingTest(), EPSStoppingThreshold(), EPSConvergedReason, EPSGetConvergedReason()
@*/
PetscErrorCode EPSStoppingBasic(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx)
{
  PetscFunctionBegin;
  *reason = EPS_CONVERGED_ITERATING;
  if (nconv >= nev) {
    PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: %" PetscInt_FMT " eigenpairs converged at iteration %" PetscInt_FMT "\n",nconv,its));
    *reason = EPS_CONVERGED_TOL;
  } else if (its >= max_it) {
    *reason = EPS_DIVERGED_ITS;
    PetscCall(PetscInfo(eps,"Linear eigensolver iteration reached maximum number of iterations (%" PetscInt_FMT ")\n",its));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSStoppingThreshold - Routine to determine whether the outer eigenvalue solver
   iteration must be stopped, according to some threshold for the computed values.

   Collective

   Input Parameters:
+  eps    - eigenvalue solver context obtained from EPSCreate()
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged eigenpairs (ignored here)
.  nev    - number of requested eigenpairs (ignored here)
-  ctx    - context containing additional data (EPSStoppingCtx)

   Output Parameter:
.  reason - result of the stopping test

   Notes:
   A positive value of reason indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, reason must be set to EPS_CONVERGED_ITERATING
   (zero).

   EPSStoppingThreshold() will stop when one of the computed eigenvalues is not
   above/below the threshold given at EPSSetThreshold(). If a number of wanted
   eigenvalues has been specified via EPSSetDimensions() then it is also taken into
   account, and the solver will stop when one of the two conditions (threshold or
   number of converged values) is met.

   Use EPSSetStoppingTest() to provide your own test instead of using this one.

   Level: advanced

.seealso: EPSSetStoppingTest(), EPSStoppingBasic(), EPSSetThreshold(), EPSSetDimensions(), EPSConvergedReason, EPSGetConvergedReason()
@*/
PetscErrorCode EPSStoppingThreshold(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx)
{
  PetscReal thres,firstev,lastev;
  PetscBool magnit,rel;
  EPSWhich  which;

  PetscFunctionBegin;
  *reason = EPS_CONVERGED_ITERATING;
  firstev = ((EPSStoppingCtx)ctx)->firstev;
  lastev  = ((EPSStoppingCtx)ctx)->lastev;
  thres   = ((EPSStoppingCtx)ctx)->thres;
  rel     = ((EPSStoppingCtx)ctx)->threlative;
  which   = ((EPSStoppingCtx)ctx)->which;
  magnit  = (which==EPS_SMALLEST_MAGNITUDE || which==EPS_LARGEST_MAGNITUDE || which==EPS_TARGET_MAGNITUDE)? PETSC_TRUE: PETSC_FALSE;
  if (nconv && magnit && which==EPS_TARGET_MAGNITUDE && ((rel && ((thres>1.0 && lastev>thres*firstev) || (thres<1.0 && lastev<thres*firstev))) || (!rel && lastev>thres))) {
    if (!rel) PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: the eigenvalue magnitude %g is above the threshold %g\n",(double)lastev,(double)thres));
    else if (thres>1.0) PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: the ratio %g/%g is above the threshold %g\n",(double)lastev,(double)firstev,(double)thres));
    else PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: the ratio %g/%g is below the threshold %g\n",(double)lastev,(double)firstev,(double)thres));
    *reason = EPS_CONVERGED_TOL;
  } else if (nconv && magnit && ((which==EPS_LARGEST_MAGNITUDE && ((rel && lastev<thres*firstev) || (!rel && lastev<thres))) || (which==EPS_SMALLEST_MAGNITUDE && lastev>thres))) {
    if (which==EPS_SMALLEST_MAGNITUDE) PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: the eigenvalue magnitude %g is above the threshold %g\n",(double)lastev,(double)thres));
    else if (!rel) PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: the eigenvalue magnitude %g is below the threshold %g\n",(double)lastev,(double)thres));
    else PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: the ratio %g/%g is below the threshold %g\n",(double)lastev,(double)firstev,(double)thres));
    *reason = EPS_CONVERGED_TOL;
  } else if (nconv && !magnit && ((which==EPS_LARGEST_REAL && lastev<thres) || (which==EPS_SMALLEST_REAL && lastev>thres))) {
    if (which==EPS_LARGEST_REAL) PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: eigenvalue %g is below the threshold %g\n",(double)lastev,(double)thres));
    else PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: eigenvalue %g is above the threshold %g\n",(double)lastev,(double)thres));
    *reason = EPS_CONVERGED_TOL;
  } else if (nev && nconv >= nev) {
    PetscCall(PetscInfo(eps,"Linear eigensolver finished successfully: %" PetscInt_FMT " eigenpairs converged at iteration %" PetscInt_FMT "\n",nconv,its));
    *reason = EPS_CONVERGED_TOL;
  } else if (its >= max_it) {
    *reason = EPS_DIVERGED_ITS;
    PetscCall(PetscInfo(eps,"Linear eigensolver iteration reached maximum number of iterations (%" PetscInt_FMT ")\n",its));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSComputeRitzVector - Computes the current Ritz vector.

  Simple case (complex scalars or real scalars with Zi=NULL):
    x = V*Zr  (V is a basis of nv vectors, Zr has length nv)

  Split case:
    x = V*Zr  y = V*Zi  (Zr and Zi have length nv)
*/
PetscErrorCode EPSComputeRitzVector(EPS eps,PetscScalar *Zr,PetscScalar *Zi,BV V,Vec x,Vec y)
{
  PetscInt       l,k;
  PetscReal      norm;
#if !defined(PETSC_USE_COMPLEX)
  Vec            z;
#endif

  PetscFunctionBegin;
  /* compute eigenvector */
  PetscCall(BVGetActiveColumns(V,&l,&k));
  PetscCall(BVSetActiveColumns(V,0,k));
  PetscCall(BVMultVec(V,1.0,0.0,x,Zr));

  /* purify eigenvector if necessary */
  if (eps->purify) {
    PetscCall(STApply(eps->st,x,y));
    if (eps->ishermitian) PetscCall(BVNormVec(eps->V,y,NORM_2,&norm));
    else PetscCall(VecNorm(y,NORM_2,&norm));
    PetscCall(VecScale(y,1.0/norm));
    PetscCall(VecCopy(y,x));
  }
  /* fix eigenvector if balancing is used */
  if (!eps->ishermitian && eps->balance!=EPS_BALANCE_NONE && eps->D) PetscCall(VecPointwiseDivide(x,x,eps->D));
#if !defined(PETSC_USE_COMPLEX)
  /* compute imaginary part of eigenvector */
  if (Zi) {
    PetscCall(BVMultVec(V,1.0,0.0,y,Zi));
    if (eps->ispositive) {
      PetscCall(BVCreateVec(V,&z));
      PetscCall(STApply(eps->st,y,z));
      PetscCall(VecNorm(z,NORM_2,&norm));
      PetscCall(VecScale(z,1.0/norm));
      PetscCall(VecCopy(z,y));
      PetscCall(VecDestroy(&z));
    }
    if (eps->balance!=EPS_BALANCE_NONE && eps->D) PetscCall(VecPointwiseDivide(y,y,eps->D));
  } else
#endif
    PetscCall(VecSet(y,0.0));

  /* normalize eigenvectors (when using balancing) */
  if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
#if !defined(PETSC_USE_COMPLEX)
    if (Zi) PetscCall(VecNormalizeComplex(x,y,PETSC_TRUE,NULL));
    else
#endif
    PetscCall(VecNormalize(x,NULL));
  }
  PetscCall(BVSetActiveColumns(V,l,k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSBuildBalance_Krylov - uses a Krylov subspace method to compute the
  diagonal matrix to be applied for balancing in non-Hermitian problems.
*/
PetscErrorCode EPSBuildBalance_Krylov(EPS eps)
{
  Vec               z,p,r;
  PetscInt          i,j;
  PetscReal         norma;
  PetscScalar       *pz,*pD;
  const PetscScalar *pr,*pp;
  PetscRandom       rand;

  PetscFunctionBegin;
  PetscCall(EPSSetWorkVecs(eps,3));
  PetscCall(BVGetRandomContext(eps->V,&rand));
  r = eps->work[0];
  p = eps->work[1];
  z = eps->work[2];
  PetscCall(VecSet(eps->D,1.0));

  for (j=0;j<eps->balance_its;j++) {

    /* Build a random vector of +-1's */
    PetscCall(VecSetRandom(z,rand));
    PetscCall(VecGetArray(z,&pz));
    for (i=0;i<eps->nloc;i++) {
      if (PetscRealPart(pz[i])<0.5) pz[i]=-1.0;
      else pz[i]=1.0;
    }
    PetscCall(VecRestoreArray(z,&pz));

    /* Compute p=DA(D\z) */
    PetscCall(VecPointwiseDivide(r,z,eps->D));
    PetscCall(STApply(eps->st,r,p));
    PetscCall(VecPointwiseMult(p,p,eps->D));
    if (eps->balance == EPS_BALANCE_TWOSIDE) {
      if (j==0) {
        /* Estimate the matrix inf-norm */
        PetscCall(VecAbs(p));
        PetscCall(VecMax(p,NULL,&norma));
      }
      /* Compute r=D\(A'Dz) */
      PetscCall(VecPointwiseMult(z,z,eps->D));
      PetscCall(STApplyHermitianTranspose(eps->st,z,r));
      PetscCall(VecPointwiseDivide(r,r,eps->D));
    }

    /* Adjust values of D */
    PetscCall(VecGetArrayRead(r,&pr));
    PetscCall(VecGetArrayRead(p,&pp));
    PetscCall(VecGetArray(eps->D,&pD));
    for (i=0;i<eps->nloc;i++) {
      if (eps->balance == EPS_BALANCE_TWOSIDE) {
        if (PetscAbsScalar(pp[i])>eps->balance_cutoff*norma && pr[i]!=0.0)
          pD[i] *= PetscSqrtReal(PetscAbsScalar(pr[i]/pp[i]));
      } else {
        if (pp[i]!=0.0) pD[i] /= PetscAbsScalar(pp[i]);
      }
    }
    PetscCall(VecRestoreArrayRead(r,&pr));
    PetscCall(VecRestoreArrayRead(p,&pp));
    PetscCall(VecRestoreArray(eps->D,&pD));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
