/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Common subroutines for all Krylov-type PEP solvers
*/

#include <slepc/private/pepimpl.h>    /*I "slepcpep.h" I*/
#include <slepcblaslapack.h>
#include "pepkrylov.h"

PetscErrorCode PEPExtractVectors_TOAR(PEP pep)
{
  PetscInt          i,j,nq,deg=pep->nmat-1,lds,idxcpy=0,ldds,k,ld;
  PetscScalar       *X,*er,*ei,*SS,*vals,*ivals,sone=1.0,szero=0.0,*yi,*yr,*tr,*ti,alpha,*pS0;
  const PetscScalar *S;
  PetscBLASInt      k_,nq_,lds_,one=1,ldds_,cols,zero=0;
  PetscBool         flg;
  PetscReal         norm,max,t,factor=1.0,done=1.0;
  Vec               xr,xi,w[4];
  PEP_TOAR          *ctx = (PEP_TOAR*)pep->data;
  Mat               S0,MS;

  PetscFunctionBegin;
  PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
  PetscCall(MatDenseGetArrayRead(MS,&S));
  PetscCall(BVGetSizes(pep->V,NULL,NULL,&ld));
  PetscCall(BVGetActiveColumns(pep->V,NULL,&nq));
  k = pep->nconv;
  if (k==0) PetscFunctionReturn(PETSC_SUCCESS);
  lds = deg*ld;
  PetscCall(DSGetLeadingDimension(pep->ds,&ldds));
  PetscCall(PetscCalloc5(k,&er,k,&ei,nq*k,&SS,pep->nmat,&vals,pep->nmat,&ivals));
  PetscCall(STGetTransform(pep->st,&flg));
  if (flg) factor = pep->sfactor;
  for (i=0;i<k;i++) {
    er[i] = factor*pep->eigr[i];
    ei[i] = factor*pep->eigi[i];
  }
  PetscCall(STBackTransform(pep->st,k,er,ei));

  PetscCall(DSVectors(pep->ds,DS_MAT_X,NULL,NULL));
  PetscCall(DSGetArray(pep->ds,DS_MAT_X,&X));

  PetscCall(PetscBLASIntCast(k,&k_));
  PetscCall(PetscBLASIntCast(nq,&nq_));
  PetscCall(PetscBLASIntCast(lds,&lds_));
  PetscCall(PetscBLASIntCast(ldds,&ldds_));

  if (pep->extract==PEP_EXTRACT_NONE || pep->refine==PEP_REFINE_MULTIPLE) {
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&nq_,&k_,&k_,&sone,S,&lds_,X,&ldds_,&szero,SS,&nq_));
  } else {
    switch (pep->extract) {
    case PEP_EXTRACT_NONE:
      break;
    case PEP_EXTRACT_NORM:
      for (i=0;i<k;i++) {
        PetscCall(PEPEvaluateBasis(pep,er[i],ei[i],vals,ivals));
        max = 1.0;
        for (j=1;j<deg;j++) {
          norm = SlepcAbsEigenvalue(vals[j],ivals[j]);
          if (max < norm) { max = norm; idxcpy = j; }
        }
        PetscCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+idxcpy*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
#if !defined(PETSC_USE_COMPLEX)
        if (PetscRealPart(ei[i])!=0.0) {
          i++;
          PetscCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+idxcpy*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
        }
#endif
      }
      break;
    case PEP_EXTRACT_RESIDUAL:
      PetscCall(VecDuplicate(pep->work[0],&xr));
      PetscCall(VecDuplicate(pep->work[0],&w[0]));
      PetscCall(VecDuplicate(pep->work[0],&w[1]));
#if !defined(PETSC_USE_COMPLEX)
      PetscCall(VecDuplicate(pep->work[0],&w[2]));
      PetscCall(VecDuplicate(pep->work[0],&w[3]));
      PetscCall(VecDuplicate(pep->work[0],&xi));
#else
      xi = NULL;
#endif
      for (i=0;i<k;i++) {
        max = 0.0;
        for (j=0;j<deg;j++) {
          PetscCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+j*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
          PetscCall(BVMultVec(pep->V,1.0,0.0,xr,SS+i*nq));
#if !defined(PETSC_USE_COMPLEX)
          PetscCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+j*ld,&lds_,X+(i+1)*ldds,&one,&szero,SS+i*nq,&one));
          PetscCall(BVMultVec(pep->V,1.0,0.0,xi,SS+i*nq));
#endif
          PetscCall(PEPComputeResidualNorm_Private(pep,er[i],ei[i],xr,xi,w,&norm));
          if (norm>max) { max = norm; idxcpy=j; }
        }
        PetscCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+idxcpy*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
#if !defined(PETSC_USE_COMPLEX)
        if (PetscRealPart(ei[i])!=0.0) {
          i++;
          PetscCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+idxcpy*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
        }
#endif
      }
      PetscCall(VecDestroy(&xr));
      PetscCall(VecDestroy(&w[0]));
      PetscCall(VecDestroy(&w[1]));
#if !defined(PETSC_USE_COMPLEX)
      PetscCall(VecDestroy(&w[2]));
      PetscCall(VecDestroy(&w[3]));
      PetscCall(VecDestroy(&xi));
#endif
      break;
    case PEP_EXTRACT_STRUCTURED:
      PetscCall(PetscMalloc2(k,&tr,k,&ti));
      for (i=0;i<k;i++) {
        t = 0.0;
        PetscCall(PEPEvaluateBasis(pep,er[i],ei[i],vals,ivals));
        yr = X+i*ldds; yi = NULL;
#if !defined(PETSC_USE_COMPLEX)
        if (ei[i]!=0.0) { yr = tr; yi = ti; }
#endif
        for (j=0;j<deg;j++) {
          alpha = PetscConj(vals[j]);
#if !defined(PETSC_USE_COMPLEX)
          if (ei[i]!=0.0) {
            PetscCall(PetscArrayzero(tr,k));
            PetscCallBLAS("BLASaxpy",BLASaxpy_(&k_,&vals[j],X+i*ldds,&one,tr,&one));
            PetscCallBLAS("BLASaxpy",BLASaxpy_(&k_,&ivals[j],X+(i+1)*ldds,&one,tr,&one));
            PetscCall(PetscArrayzero(ti,k));
            PetscCallBLAS("BLASaxpy",BLASaxpy_(&k_,&vals[j],X+(i+1)*ldds,&one,ti,&one));
            alpha = -ivals[j];
            PetscCallBLAS("BLASaxpy",BLASaxpy_(&k_,&alpha,X+i*ldds,&one,ti,&one));
            alpha = 1.0;
          }
#endif
          PetscCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&alpha,S+j*ld,&lds_,yr,&one,&sone,SS+i*nq,&one));
          t += SlepcAbsEigenvalue(vals[j],ivals[j])*SlepcAbsEigenvalue(vals[j],ivals[j]);
          if (yi) {
            PetscCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&alpha,S+j*ld,&lds_,yi,&one,&sone,SS+(i+1)*nq,&one));
          }
        }
        cols = yi? 2: 1;
        PetscCallLAPACKInfo("LAPACKlascl",LAPACKlascl_("G",&zero,&zero,&t,&done,&nq_,&cols,SS+i*nq,&nq_,&info));
        if (yi) i++;
      }
      PetscCall(PetscFree2(tr,ti));
      break;
    }
  }

  /* update vectors V = V*S */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,nq,k,NULL,&S0));
  PetscCall(MatDenseGetArrayWrite(S0,&pS0));
  for (i=0;i<k;i++) PetscCall(PetscArraycpy(pS0+i*nq,SS+i*nq,nq));
  PetscCall(MatDenseRestoreArrayWrite(S0,&pS0));
  PetscCall(BVMultInPlace(pep->V,S0,0,k));
  PetscCall(MatDestroy(&S0));
  PetscCall(PetscFree5(er,ei,SS,vals,ivals));
  PetscCall(MatDenseRestoreArrayRead(MS,&S));
  PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   PEPKrylovConvergence - This is the analogue to EPSKrylovConvergence, but
   for polynomial Krylov methods.

   Differences:
   - Always non-symmetric
   - Does not check for STSHIFT
   - No correction factor
   - No support for true residual
*/
PetscErrorCode PEPKrylovConvergence(PEP pep,PetscBool getall,PetscInt kini,PetscInt nits,PetscReal beta,PetscInt *kout)
{
  PetscInt       k,newk,marker,inside;
  PetscScalar    re,im;
  PetscReal      resnorm;
  PetscBool      istrivial;

  PetscFunctionBegin;
  PetscCall(RGIsTrivial(pep->rg,&istrivial));
  marker = -1;
  if (pep->trackall) getall = PETSC_TRUE;
  for (k=kini;k<kini+nits;k++) {
    /* eigenvalue */
    re = pep->eigr[k];
    im = pep->eigi[k];
    if (PetscUnlikely(!istrivial)) {
      PetscCall(STBackTransform(pep->st,1,&re,&im));
      PetscCall(RGCheckInside(pep->rg,1,&re,&im,&inside));
      if (marker==-1 && inside<0) marker = k;
      re = pep->eigr[k];
      im = pep->eigi[k];
    }
    newk = k;
    PetscCall(DSVectors(pep->ds,DS_MAT_X,&newk,&resnorm));
    resnorm *= beta;
    /* error estimate */
    PetscCall((*pep->converged)(pep,re,im,resnorm,&pep->errest[k],pep->convergedctx));
    if (marker==-1 && pep->errest[k] >= pep->tol) marker = k;
    if (PetscUnlikely(newk==k+1)) {
      pep->errest[k+1] = pep->errest[k];
      k++;
    }
    if (marker!=-1 && !getall) break;
  }
  if (marker!=-1) k = marker;
  *kout = k;
  PetscFunctionReturn(PETSC_SUCCESS);
}
