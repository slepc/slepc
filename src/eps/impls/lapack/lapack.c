/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the LAPACK eigenvalue subroutines.
   Generalized problems are transformed to standard ones only if necessary.
*/

#include <slepc/private/epsimpl.h>

static PetscErrorCode EPSSetUp_LAPACK(EPS eps)
{
  int            ierra,ierrb;
  PetscBool      isshift,flg,denseok=PETSC_FALSE;
  Mat            A,B,OP,shell,Ar,Br,Adense=NULL,Bdense=NULL,Ads,Bds;
  PetscScalar    shift;
  PetscInt       nmat;
  KSP            ksp;
  PC             pc;

  PetscFunctionBegin;
  if (eps->nev==0) eps->nev = 1;
  eps->ncv = eps->n;
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = 1;
  if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  PetscCheck(eps->which!=EPS_ALL || eps->inta==eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support interval computation");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_STOPPING);
  PetscCall(EPSAllocateSolution(eps,0));

  /* attempt to get dense representations of A and B separately */
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isshift));
  if (isshift) {
    PetscCall(STGetNumMatrices(eps->st,&nmat));
    PetscCall(STGetMatrix(eps->st,0,&A));
    PetscCall(MatHasOperation(A,MATOP_CREATE_SUBMATRICES,&flg));
    if (flg) {
      PetscCall(PetscPushErrorHandler(PetscReturnErrorHandler,NULL));
      ierra  = MatCreateRedundantMatrix(A,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Ar);
      if (!ierra) ierra |= MatConvert(Ar,MATSEQDENSE,MAT_INITIAL_MATRIX,&Adense);
      ierra |= MatDestroy(&Ar);
      PetscCall(PetscPopErrorHandler());
    } else ierra = 1;
    if (nmat>1) {
      PetscCall(STGetMatrix(eps->st,1,&B));
      PetscCall(MatHasOperation(B,MATOP_CREATE_SUBMATRICES,&flg));
      if (flg) {
        PetscCall(PetscPushErrorHandler(PetscReturnErrorHandler,NULL));
        ierrb  = MatCreateRedundantMatrix(B,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Br);
        if (!ierrb) ierrb |= MatConvert(Br,MATSEQDENSE,MAT_INITIAL_MATRIX,&Bdense);
        ierrb |= MatDestroy(&Br);
        PetscCall(PetscPopErrorHandler());
      } else ierrb = 1;
    } else ierrb = 0;
    denseok = PetscNot(ierra || ierrb);
  }

  /* setup DS */
  if (denseok) {
    if (eps->isgeneralized) {
      if (eps->ishermitian) {
        if (eps->ispositive) PetscCall(DSSetType(eps->ds,DSGHEP));
        else PetscCall(DSSetType(eps->ds,DSGNHEP)); /* TODO: should be DSGHIEP */
      } else PetscCall(DSSetType(eps->ds,DSGNHEP));
    } else {
      if (eps->ishermitian) PetscCall(DSSetType(eps->ds,DSHEP));
      else PetscCall(DSSetType(eps->ds,DSNHEP));
    }
  } else PetscCall(DSSetType(eps->ds,DSNHEP));
  PetscCall(DSAllocate(eps->ds,eps->ncv));
  PetscCall(DSSetDimensions(eps->ds,eps->ncv,0,0));

  if (denseok) {
    PetscCall(STGetShift(eps->st,&shift));
    if (shift != 0.0) {
      if (nmat>1) PetscCall(MatAXPY(Adense,-shift,Bdense,SAME_NONZERO_PATTERN));
      else PetscCall(MatShift(Adense,-shift));
    }
    /* use dummy pc and ksp to avoid problems when B is not positive definite */
    PetscCall(STGetKSP(eps->st,&ksp));
    PetscCall(KSPSetType(ksp,KSPPREONLY));
    PetscCall(KSPGetPC(ksp,&pc));
    PetscCall(PCSetType(pc,PCNONE));
  } else {
    PetscCall(PetscInfo(eps,"Using slow explicit operator\n"));
    PetscCall(STGetOperator(eps->st,&shell));
    PetscCall(MatComputeOperator(shell,MATDENSE,&OP));
    PetscCall(STRestoreOperator(eps->st,&shell));
    PetscCall(MatDestroy(&Adense));
    PetscCall(MatCreateRedundantMatrix(OP,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Adense));
    PetscCall(MatDestroy(&OP));
  }

  /* fill DS matrices */
  PetscCall(DSGetMat(eps->ds,DS_MAT_A,&Ads));
  PetscCall(MatCopy(Adense,Ads,SAME_NONZERO_PATTERN));
  PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&Ads));
  if (denseok && eps->isgeneralized) {
    PetscCall(DSGetMat(eps->ds,DS_MAT_B,&Bds));
    PetscCall(MatCopy(Bdense,Bds,SAME_NONZERO_PATTERN));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_B,&Bds));
  }
  PetscCall(DSSetState(eps->ds,DS_STATE_RAW));
  PetscCall(MatDestroy(&Adense));
  PetscCall(MatDestroy(&Bdense));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_LAPACK(EPS eps)
{
  PetscInt       n=eps->n,i,low,high;
  PetscScalar    *array,*pX,*pY;
  Vec            v,w;

  PetscFunctionBegin;
  PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
  PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
  PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

  /* right eigenvectors */
  PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));
  PetscCall(DSGetArray(eps->ds,DS_MAT_X,&pX));
  for (i=0;i<eps->ncv;i++) {
    PetscCall(BVGetColumn(eps->V,i,&v));
    PetscCall(VecGetOwnershipRange(v,&low,&high));
    PetscCall(VecGetArray(v,&array));
    PetscCall(PetscArraycpy(array,pX+i*n+low,high-low));
    PetscCall(VecRestoreArray(v,&array));
    PetscCall(BVRestoreColumn(eps->V,i,&v));
  }
  PetscCall(DSRestoreArray(eps->ds,DS_MAT_X,&pX));

  /* left eigenvectors */
  if (eps->twosided) {
    PetscCall(DSVectors(eps->ds,DS_MAT_Y,NULL,NULL));
    PetscCall(DSGetArray(eps->ds,DS_MAT_Y,&pY));
    for (i=0;i<eps->ncv;i++) {
      PetscCall(BVGetColumn(eps->W,i,&w));
      PetscCall(VecGetOwnershipRange(w,&low,&high));
      PetscCall(VecGetArray(w,&array));
      PetscCall(PetscArraycpy(array,pY+i*n+low,high-low));
      PetscCall(VecRestoreArray(w,&array));
      PetscCall(BVRestoreColumn(eps->W,i,&w));
    }
    PetscCall(DSRestoreArray(eps->ds,DS_MAT_Y,&pY));
  }

  eps->nconv  = eps->ncv;
  eps->its    = 1;
  eps->reason = EPS_CONVERGED_TOL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_LAPACK(EPS eps)
{
  PetscFunctionBegin;
  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_OTHER;

  eps->ops->solve          = EPSSolve_LAPACK;
  eps->ops->setup          = EPSSetUp_LAPACK;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->backtransform  = EPSBackTransform_Default;
  PetscFunctionReturn(PETSC_SUCCESS);
}
