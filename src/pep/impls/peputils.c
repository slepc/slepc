/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Common subroutines for several PEP solvers
*/

#include <slepc/private/pepimpl.h>    /*I "slepcpep.h" I*/
#include <slepcblaslapack.h>

/*
  Computes T_j=phy_idx(T) for a matrix T.
    Tp (if j>0) and Tpp (if j>1) are the evaluation
    of phy_(j-1) and phy(j-2)respectively.
*/
PetscErrorCode PEPEvaluateBasisMat(PEP pep,PetscInt k,PetscScalar *T,PetscInt ldt,PetscInt idx,PetscScalar *Tpp,PetscInt ldtpp,PetscScalar *Tp,PetscInt ldtp,PetscScalar *Tj,PetscInt ldtj)
{
  PetscInt       i;
  PetscReal      *ca,*cb,*cg;
  PetscScalar    g,a;
  PetscBLASInt   k_,ldt_,ldtj_,ldtp_;

  PetscFunctionBegin;
  if (idx==0) {
    for (i=0;i<k;i++) {
      PetscCall(PetscArrayzero(Tj+i*ldtj,k));
      Tj[i+i*ldtj] = 1.0;
    }
  } else {
    PetscCall(PetscBLASIntCast(ldt,&ldt_));
    PetscCall(PetscBLASIntCast(ldtj,&ldtj_));
    PetscCall(PetscBLASIntCast(ldtp,&ldtp_));
    PetscCall(PetscBLASIntCast(k,&k_));
    ca = pep->pbc; cb = pep->pbc+pep->nmat; cg = pep->pbc+2*pep->nmat;
    for (i=0;i<k;i++) T[i*ldt+i] -= cb[idx-1];
    if (idx>1) {
      for (i=0;i<k;i++) PetscCall(PetscArraycpy(Tj+i*ldtj,Tpp+i*ldtpp,k));
    }
    a = 1/ca[idx-1];
    g = (idx==1)?0.0:-cg[idx-1]/ca[idx-1];
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&a,T,&ldt_,Tp,&ldtp_,&g,Tj,&ldtj_));
    for (i=0;i<k;i++) T[i*ldt+i] += cb[idx-1];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPEvaluateBasis - evaluate the polynomial basis on a given parameter sigma
*/
PetscErrorCode PEPEvaluateBasis(PEP pep,PetscScalar sigma,PetscScalar isigma,PetscScalar *vals,PetscScalar *ivals)
{
  PetscInt   nmat=pep->nmat,k;
  PetscReal  *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat;

  PetscFunctionBegin;
  if (ivals) for (k=0;k<nmat;k++) ivals[k] = 0.0;
  vals[0] = 1.0;
  vals[1] = (sigma-b[0])/a[0];
#if !PetscDefined(USE_COMPLEX)
  if (ivals) ivals[1] = isigma/a[0];
#endif
  for (k=2;k<nmat;k++) {
    vals[k] = ((sigma-b[k-1])*vals[k-1]-g[k-1]*vals[k-2])/a[k-1];
    if (ivals) vals[k] -= isigma*ivals[k-1]/a[k-1];
#if !PetscDefined(USE_COMPLEX)
    if (ivals) ivals[k] = ((sigma-b[k-1])*ivals[k-1]+isigma*vals[k-1]-g[k-1]*ivals[k-2])/a[k-1];
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPEvaluateBasisDerivative - evaluate the derivative of the polynomial basis on a given parameter sigma
*/
PetscErrorCode PEPEvaluateBasisDerivative(PEP pep,PetscScalar sigma,PetscScalar isigma,PetscScalar *vals,PetscScalar *ivals)
{
  PetscInt       nmat=pep->nmat,k;
  PetscReal      *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat;

  PetscFunctionBegin;
  PetscCall(PEPEvaluateBasis(pep,sigma,isigma,vals,ivals));
  for (k=nmat-1;k>0;k--) {
    vals[k] = vals[k-1];
    if (ivals) ivals[k] = ivals[k-1];
  }
  vals[0] = 0.0;
  vals[1] = vals[1]/a[0];
#if !PetscDefined(USE_COMPLEX)
  if (ivals) ivals[1] = ivals[1]/a[0];
#endif
  for (k=2;k<nmat;k++) {
    vals[k] += (sigma-b[k-1])*vals[k-1]-g[k-1]*vals[k-2];
    if (ivals) vals[k] -= isigma*ivals[k-1];
    vals[k] /= a[k-1];
#if !PetscDefined(USE_COMPLEX)
    if (ivals) {
      ivals[k] += (sigma-b[k-1])*ivals[k-1]+isigma*vals[k-1]-g[k-1]*ivals[k-2];
      ivals[k] /= a[k-1];
    }
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
