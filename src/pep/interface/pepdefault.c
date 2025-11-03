/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Simple default routines for common PEP operations
*/

#include <slepc/private/pepimpl.h>     /*I "slepcpep.h" I*/

/*@
   PEPSetWorkVecs - Sets a number of work vectors into a `PEP` object.

   Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
-  nw  - number of work vectors to allocate

   Developer Note:
   This is `SLEPC_EXTERN` because it may be required by user plugin `PEP`
   implementations.

   Level: developer

.seealso: [](ch:pep), `PEPSetUp()`
@*/
PetscErrorCode PEPSetWorkVecs(PEP pep,PetscInt nw)
{
  Vec            t;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,nw,2);
  PetscCheck(nw>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"nw must be > 0: nw = %" PetscInt_FMT,nw);
  if (pep->nwork < nw) {
    PetscCall(VecDestroyVecs(pep->nwork,&pep->work));
    pep->nwork = nw;
    PetscCall(BVGetColumn(pep->V,0,&t));
    PetscCall(VecDuplicateVecs(t,nw,&pep->work));
    PetscCall(BVRestoreColumn(pep->V,0,&t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  PEPConvergedRelative - Checks convergence relative to the eigenvalue.
*/
PetscErrorCode PEPConvergedRelative(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal w;

  PetscFunctionBegin;
  w = SlepcAbsEigenvalue(eigr,eigi);
  *errest = (w!=0.0)? res/w: PETSC_MAX_REAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  PEPConvergedNorm - Checks convergence relative to the matrix norms.
*/
PetscErrorCode PEPConvergedNorm(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal      w=0.0,t;
  PetscInt       j;
  PetscBool      flg;

  PetscFunctionBegin;
  /* initialization of matrix norms */
  if (!pep->nrma[pep->nmat-1]) {
    for (j=0;j<pep->nmat;j++) {
      PetscCall(MatHasOperation(pep->A[j],MATOP_NORM,&flg));
      PetscCheck(flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"The convergence test related to the matrix norms requires a matrix norm operation");
      PetscCall(MatNorm(pep->A[j],NORM_INFINITY,&pep->nrma[j]));
    }
  }
  t = SlepcAbsEigenvalue(eigr,eigi);
  for (j=pep->nmat-1;j>=0;j--) {
    w = w*t+pep->nrma[j];
  }
  *errest = res/w;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  PEPSetWhichEigenpairs_Default - Sets the default value for which,
  depending on the ST.
 */
PetscErrorCode PEPSetWhichEigenpairs_Default(PEP pep)
{
  PetscBool      target;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&target));
  if (target) pep->which = PEP_TARGET_MAGNITUDE;
  else pep->which = PEP_LARGEST_MAGNITUDE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  PEPConvergedAbsolute - Checks convergence absolutely.
*/
PetscErrorCode PEPConvergedAbsolute(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = res;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPStoppingBasic - Default routine to determine whether the outer eigensolver
   iteration must be stopped.

   Collective

   Input Parameters:
+  pep    - the polynomial eigensolver context
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged eigenpairs
.  nev    - number of requested eigenpairs
-  ctx    - context (not used here)

   Output Parameter:
.  reason - result of the stopping test

   Notes:
   `PEPStoppingBasic()` will stop if all requested eigenvalues are converged, or if
   the maximum number of iterations has been reached.

   This is the default stopping test.
   Use `PEPSetStoppingTest()` to provide your own test instead of using this one.

   Level: advanced

.seealso: [](ch:pep), `PEPSetStoppingTest()`, `PEPConvergedReason`, `PEPGetConvergedReason()`
@*/
PetscErrorCode PEPStoppingBasic(PEP pep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,PEPConvergedReason *reason,void *ctx)
{
  PetscFunctionBegin;
  *reason = PEP_CONVERGED_ITERATING;
  if (nconv >= nev) {
    PetscCall(PetscInfo(pep,"Polynomial eigensolver finished successfully: %" PetscInt_FMT " eigenpairs converged at iteration %" PetscInt_FMT "\n",nconv,its));
    *reason = PEP_CONVERGED_TOL;
  } else if (its >= max_it) {
    *reason = PEP_DIVERGED_ITS;
    PetscCall(PetscInfo(pep,"Polynomial eigensolver iteration reached maximum number of iterations (%" PetscInt_FMT ")\n",its));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPBackTransform_Default(PEP pep)
{
  PetscFunctionBegin;
  PetscCall(STBackTransform(pep->st,pep->nconv,pep->eigr,pep->eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPComputeVectors_Default(PEP pep)
{
  PetscInt       i;
  Vec            v;

  PetscFunctionBegin;
  PetscCall(PEPExtractVectors(pep));

  /* Fix eigenvectors if balancing was used */
  if ((pep->scale==PEP_SCALE_DIAGONAL || pep->scale==PEP_SCALE_BOTH) && pep->Dr && (pep->refine!=PEP_REFINE_MULTIPLE)) {
    for (i=0;i<pep->nconv;i++) {
      PetscCall(BVGetColumn(pep->V,i,&v));
      PetscCall(VecPointwiseMult(v,v,pep->Dr));
      PetscCall(BVRestoreColumn(pep->V,i,&v));
    }
  }

  /* normalization */
  PetscCall(BVNormalize(pep->V,pep->eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  PEPBuildDiagonalScaling - compute two diagonal matrices to be applied for balancing
  in polynomial eigenproblems.
*/
PetscErrorCode PEPBuildDiagonalScaling(PEP pep)
{
  PetscInt       it,i,j,k,nmat,nr,e,nz,lst,lend,nc=0,*cols,emax,emin,emaxl,eminl;
  const PetscInt *cidx,*ridx;
  Mat            M,*T,A;
  PetscMPIInt    n;
  PetscBool      cont=PETSC_TRUE,flg=PETSC_FALSE;
  PetscScalar    *array,*Dr,*Dl,t;
  PetscReal      l2,d,*rsum,*aux,*csum,w=1.0;
  MatStructure   str;
  MatInfo        info;

  PetscFunctionBegin;
  l2 = 2*PetscLogReal(2.0);
  nmat = pep->nmat;
  PetscCall(PetscMPIIntCast(pep->n,&n));
  PetscCall(STGetMatStructure(pep->st,&str));
  PetscCall(PetscMalloc1(nmat,&T));
  for (k=0;k<nmat;k++) PetscCall(STGetMatrixTransformed(pep->st,k,&T[k]));
  /* Form local auxiliary matrix M */
  PetscCall(PetscObjectBaseTypeCompareAny((PetscObject)T[0],&cont,MATMPIAIJ,MATSEQAIJ,""));
  PetscCheck(cont,PetscObjectComm((PetscObject)T[0]),PETSC_ERR_SUP,"Only for MPIAIJ or SEQAIJ matrix types");
  PetscCall(PetscObjectBaseTypeCompare((PetscObject)T[0],MATMPIAIJ,&cont));
  if (cont) {
    PetscCall(MatMPIAIJGetLocalMat(T[0],MAT_INITIAL_MATRIX,&M));
    flg = PETSC_TRUE;
  } else PetscCall(MatDuplicate(T[0],MAT_COPY_VALUES,&M));
  PetscCall(MatGetInfo(M,MAT_LOCAL,&info));
  nz = (PetscInt)info.nz_used;
  PetscCall(MatSeqAIJGetArray(M,&array));
  for (i=0;i<nz;i++) {
    t = PetscAbsScalar(array[i]);
    array[i] = t*t;
  }
  PetscCall(MatSeqAIJRestoreArray(M,&array));
  for (k=1;k<nmat;k++) {
    if (flg) PetscCall(MatMPIAIJGetLocalMat(T[k],MAT_INITIAL_MATRIX,&A));
    else {
      if (str==SAME_NONZERO_PATTERN) PetscCall(MatCopy(T[k],A,SAME_NONZERO_PATTERN));
      else PetscCall(MatDuplicate(T[k],MAT_COPY_VALUES,&A));
    }
    PetscCall(MatGetInfo(A,MAT_LOCAL,&info));
    nz = (PetscInt)info.nz_used;
    PetscCall(MatSeqAIJGetArray(A,&array));
    for (i=0;i<nz;i++) {
      t = PetscAbsScalar(array[i]);
      array[i] = t*t;
    }
    PetscCall(MatSeqAIJRestoreArray(A,&array));
    w *= pep->slambda*pep->slambda*pep->sfactor;
    PetscCall(MatAXPY(M,w,A,str));
    if (flg || str!=SAME_NONZERO_PATTERN || k==nmat-2) PetscCall(MatDestroy(&A));
  }
  PetscCall(MatGetRowIJ(M,0,PETSC_FALSE,PETSC_FALSE,&nr,&ridx,&cidx,&cont));
  PetscCheck(cont,PetscObjectComm((PetscObject)T[0]),PETSC_ERR_SUP,"It is not possible to compute scaling diagonals for these PEP matrices");
  PetscCall(MatGetInfo(M,MAT_LOCAL,&info));
  nz = (PetscInt)info.nz_used;
  PetscCall(VecGetOwnershipRange(pep->Dl,&lst,&lend));
  PetscCall(PetscMalloc4(nr,&rsum,pep->n,&csum,pep->n,&aux,PetscMin(pep->n-lend+lst,nz),&cols));
  PetscCall(VecSet(pep->Dr,1.0));
  PetscCall(VecSet(pep->Dl,1.0));
  PetscCall(VecGetArray(pep->Dl,&Dl));
  PetscCall(VecGetArray(pep->Dr,&Dr));
  PetscCall(MatSeqAIJGetArray(M,&array));
  PetscCall(PetscArrayzero(aux,pep->n));
  for (j=0;j<nz;j++) {
    /* Search non-zero columns outsize lst-lend */
    if (aux[cidx[j]]==0 && (cidx[j]<lst || lend<=cidx[j])) cols[nc++] = cidx[j];
    /* Local column sums */
    aux[cidx[j]] += PetscAbsScalar(array[j]);
  }
  for (it=0;it<pep->sits && cont;it++) {
    emaxl = 0; eminl = 0;
    /* Column sum  */
    if (it>0) { /* it=0 has been already done*/
      PetscCall(MatSeqAIJGetArray(M,&array));
      PetscCall(PetscArrayzero(aux,pep->n));
      for (j=0;j<nz;j++) aux[cidx[j]] += PetscAbsScalar(array[j]);
      PetscCall(MatSeqAIJRestoreArray(M,&array));
    }
    PetscCallMPI(MPIU_Allreduce(aux,csum,n,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)pep->Dr)));
    /* Update Dr */
    for (j=lst;j<lend;j++) {
      d = PetscLogReal(csum[j])/l2;
      e = -(PetscInt)((d < 0)?(d-0.5):(d+0.5));
      d = PetscPowReal(2.0,e);
      Dr[j-lst] *= d;
      aux[j] = d*d;
      emaxl = PetscMax(emaxl,e);
      eminl = PetscMin(eminl,e);
    }
    for (j=0;j<nc;j++) {
      d = PetscLogReal(csum[cols[j]])/l2;
      e = -(PetscInt)((d < 0)?(d-0.5):(d+0.5));
      d = PetscPowReal(2.0,e);
      aux[cols[j]] = d*d;
      emaxl = PetscMax(emaxl,e);
      eminl = PetscMin(eminl,e);
    }
    /* Scale M */
    PetscCall(MatSeqAIJGetArray(M,&array));
    for (j=0;j<nz;j++) {
      array[j] *= aux[cidx[j]];
    }
    PetscCall(MatSeqAIJRestoreArray(M,&array));
    /* Row sum */
    PetscCall(PetscArrayzero(rsum,nr));
    PetscCall(MatSeqAIJGetArray(M,&array));
    for (i=0;i<nr;i++) {
      for (j=ridx[i];j<ridx[i+1];j++) rsum[i] += PetscAbsScalar(array[j]);
      /* Update Dl */
      d = PetscLogReal(rsum[i])/l2;
      e = -(PetscInt)((d < 0)?(d-0.5):(d+0.5));
      d = PetscPowReal(2.0,e);
      Dl[i] *= d;
      /* Scale M */
      for (j=ridx[i];j<ridx[i+1];j++) array[j] *= d*d;
      emaxl = PetscMax(emaxl,e);
      eminl = PetscMin(eminl,e);
    }
    PetscCall(MatSeqAIJRestoreArray(M,&array));
    /* Compute global max and min */
    PetscCallMPI(MPIU_Allreduce(&emaxl,&emax,1,MPIU_INT,MPI_MAX,PetscObjectComm((PetscObject)pep->Dl)));
    PetscCallMPI(MPIU_Allreduce(&eminl,&emin,1,MPIU_INT,MPI_MIN,PetscObjectComm((PetscObject)pep->Dl)));
    if (emax<=emin+2) cont = PETSC_FALSE;
  }
  PetscCall(VecRestoreArray(pep->Dr,&Dr));
  PetscCall(VecRestoreArray(pep->Dl,&Dl));
  /* Free memory*/
  PetscCall(MatDestroy(&M));
  PetscCall(PetscFree4(rsum,csum,aux,cols));
  PetscCall(PetscFree(T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPComputeScaleFactor - compute sfactor as described in [Betcke 2008].
*/
PetscErrorCode PEPComputeScaleFactor(PEP pep)
{
  PetscBool      has0,has1,flg;
  PetscReal      norm0,norm1;
  Mat            T[2];
  PEPBasis       basis;
  PetscInt       i;

  PetscFunctionBegin;
  if (pep->scale==PEP_SCALE_NONE || pep->scale==PEP_SCALE_DIAGONAL) {  /* no scalar scaling */
    pep->sfactor = 1.0;
    pep->dsfactor = 1.0;
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  if (pep->sfactor_set) PetscFunctionReturn(PETSC_SUCCESS);  /* user provided value */
  pep->sfactor = 1.0;
  pep->dsfactor = 1.0;
  PetscCall(PEPGetBasis(pep,&basis));
  if (basis==PEP_BASIS_MONOMIAL) {
    PetscCall(STGetTransform(pep->st,&flg));
    if (flg) {
      PetscCall(STGetMatrixTransformed(pep->st,0,&T[0]));
      PetscCall(STGetMatrixTransformed(pep->st,pep->nmat-1,&T[1]));
    } else {
      T[0] = pep->A[0];
      T[1] = pep->A[pep->nmat-1];
    }
    if (pep->nmat>2) {
      PetscCall(MatHasOperation(T[0],MATOP_NORM,&has0));
      PetscCall(MatHasOperation(T[1],MATOP_NORM,&has1));
      if (has0 && has1) {
        PetscCall(MatNorm(T[0],NORM_INFINITY,&norm0));
        PetscCall(MatNorm(T[1],NORM_INFINITY,&norm1));
        pep->sfactor = PetscPowReal(norm0/norm1,1.0/(pep->nmat-1));
        pep->dsfactor = norm1;
        for (i=pep->nmat-2;i>0;i--) {
          PetscCall(STGetMatrixTransformed(pep->st,i,&T[1]));
          PetscCall(MatHasOperation(T[1],MATOP_NORM,&has1));
          if (has1) {
            PetscCall(MatNorm(T[1],NORM_INFINITY,&norm1));
            pep->dsfactor = pep->dsfactor*pep->sfactor+norm1;
          } else break;
        }
        if (has1) {
          pep->dsfactor = pep->dsfactor*pep->sfactor+norm0;
          pep->dsfactor = pep->nmat/pep->dsfactor;
        } else pep->dsfactor = 1.0;
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPBasisCoefficients - compute polynomial basis coefficients
*/
PetscErrorCode PEPBasisCoefficients(PEP pep,PetscReal *pbc)
{
  PetscReal *ca,*cb,*cg;
  PetscInt  k,nmat=pep->nmat;

  PetscFunctionBegin;
  ca = pbc;
  cb = pbc+nmat;
  cg = pbc+2*nmat;
  switch (pep->basis) {
  case PEP_BASIS_MONOMIAL:
    for (k=0;k<nmat;k++) {
      ca[k] = 1.0; cb[k] = 0.0; cg[k] = 0.0;
    }
    break;
  case PEP_BASIS_CHEBYSHEV1:
    ca[0] = 1.0; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = .5; cb[k] = 0.0; cg[k] = .5;
    }
    break;
  case PEP_BASIS_CHEBYSHEV2:
    ca[0] = .5; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = .5; cb[k] = 0.0; cg[k] = .5;
    }
    break;
  case PEP_BASIS_LEGENDRE:
    ca[0] = 1.0; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = k+1; cb[k] = -2*k; cg[k] = k;
    }
    break;
  case PEP_BASIS_LAGUERRE:
    ca[0] = -1.0; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = -(k+1); cb[k] = 2*k+1; cg[k] = -k;
    }
    break;
  case PEP_BASIS_HERMITE:
    ca[0] = .5; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = .5; cb[k] = 0.0; cg[k] = -k;
    }
    break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
