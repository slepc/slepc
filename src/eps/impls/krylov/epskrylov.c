/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Common subroutines for all Krylov-type solvers
*/

#include <slepc/private/epsimpl.h>
#include <slepc/private/slepcimpl.h>
#include <slepcblaslapack.h>

/*
   EPSDelayedArnoldi - This function is equivalent to BVMatArnoldi but
   performs the computation in a different way. The main idea is that
   reorthogonalization is delayed to the next Arnoldi step. This version is
   more scalable but in some cases convergence may stagnate.
*/
PetscErrorCode EPSDelayedArnoldi(EPS eps,PetscScalar *H,PetscInt ldh,PetscInt k,PetscInt *M,PetscReal *beta,PetscBool *breakdown)
{
  PetscInt       i,j,m=*M;
  Vec            u,t;
  PetscScalar    shh[100],*lhh,dot,dot2;
  PetscReal      norm1=0.0,norm2=1.0;
  Vec            vj,vj1,vj2=NULL;

  PetscFunctionBegin;
  if (m<=100) lhh = shh;
  else PetscCall(PetscMalloc1(m,&lhh));
  PetscCall(BVCreateVec(eps->V,&u));
  PetscCall(BVCreateVec(eps->V,&t));

  PetscCall(BVSetActiveColumns(eps->V,0,m));
  for (j=k;j<m;j++) {
    PetscCall(BVGetColumn(eps->V,j,&vj));
    PetscCall(BVGetColumn(eps->V,j+1,&vj1));
    PetscCall(STApply(eps->st,vj,vj1));
    PetscCall(BVRestoreColumn(eps->V,j,&vj));
    PetscCall(BVRestoreColumn(eps->V,j+1,&vj1));

    PetscCall(BVDotColumnBegin(eps->V,j+1,H+ldh*j));
    if (j>k) {
      PetscCall(BVDotColumnBegin(eps->V,j,lhh));
      PetscCall(BVGetColumn(eps->V,j,&vj));
      PetscCall(VecDotBegin(vj,vj,&dot));
      if (j>k+1) {
        PetscCall(BVNormVecBegin(eps->V,u,NORM_2,&norm2));
        PetscCall(BVGetColumn(eps->V,j-2,&vj2));
        PetscCall(VecDotBegin(u,vj2,&dot2));
      }
      PetscCall(BVDotColumnEnd(eps->V,j+1,H+ldh*j));
      PetscCall(BVDotColumnEnd(eps->V,j,lhh));
      PetscCall(VecDotEnd(vj,vj,&dot));
      PetscCall(BVRestoreColumn(eps->V,j,&vj));
      if (j>k+1) {
        PetscCall(BVNormVecEnd(eps->V,u,NORM_2,&norm2));
        PetscCall(VecDotEnd(u,vj2,&dot2));
        PetscCall(BVRestoreColumn(eps->V,j-2,&vj2));
      }
      norm1 = PetscSqrtReal(PetscRealPart(dot));
      for (i=0;i<j;i++) H[ldh*j+i] = H[ldh*j+i]/norm1;
      H[ldh*j+j] = H[ldh*j+j]/dot;
      PetscCall(BVCopyVec(eps->V,j,t));
      PetscCall(BVScaleColumn(eps->V,j,1.0/norm1));
      PetscCall(BVScaleColumn(eps->V,j+1,1.0/norm1));
    } else PetscCall(BVDotColumnEnd(eps->V,j+1,H+ldh*j)); /* j==k */

    PetscCall(BVMultColumn(eps->V,-1.0,1.0,j+1,H+ldh*j));

    if (j>k) {
      PetscCall(BVSetActiveColumns(eps->V,0,j));
      PetscCall(BVMultVec(eps->V,-1.0,1.0,t,lhh));
      PetscCall(BVSetActiveColumns(eps->V,0,m));
      for (i=0;i<j;i++) H[ldh*(j-1)+i] += lhh[i];
    }

    if (j>k+1) {
      PetscCall(BVGetColumn(eps->V,j-1,&vj1));
      PetscCall(VecCopy(u,vj1));
      PetscCall(BVRestoreColumn(eps->V,j-1,&vj1));
      PetscCall(BVScaleColumn(eps->V,j-1,1.0/norm2));
      H[ldh*(j-2)+j-1] = norm2;
    }

    if (j<m-1) PetscCall(VecCopy(t,u));
  }

  PetscCall(BVNormVec(eps->V,t,NORM_2,&norm2));
  PetscCall(VecScale(t,1.0/norm2));
  PetscCall(BVGetColumn(eps->V,m-1,&vj1));
  PetscCall(VecCopy(t,vj1));
  PetscCall(BVRestoreColumn(eps->V,m-1,&vj1));
  H[ldh*(m-2)+m-1] = norm2;

  PetscCall(BVDotColumn(eps->V,m,lhh));

  PetscCall(BVMultColumn(eps->V,-1.0,1.0,m,lhh));
  for (i=0;i<m;i++)
    H[ldh*(m-1)+i] += lhh[i];

  PetscCall(BVNormColumn(eps->V,m,NORM_2,beta));
  PetscCall(BVScaleColumn(eps->V,m,1.0 / *beta));
  *breakdown = PETSC_FALSE;

  if (m>100) PetscCall(PetscFree(lhh));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&t));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSDelayedArnoldi1 - This function is similar to EPSDelayedArnoldi,
   but without reorthogonalization (only delayed normalization).
*/
PetscErrorCode EPSDelayedArnoldi1(EPS eps,PetscScalar *H,PetscInt ldh,PetscInt k,PetscInt *M,PetscReal *beta,PetscBool *breakdown)
{
  PetscInt       i,j,m=*M;
  PetscScalar    dot;
  PetscReal      norm=0.0;
  Vec            vj,vj1;

  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(eps->V,0,m));
  for (j=k;j<m;j++) {
    PetscCall(BVGetColumn(eps->V,j,&vj));
    PetscCall(BVGetColumn(eps->V,j+1,&vj1));
    PetscCall(STApply(eps->st,vj,vj1));
    PetscCall(BVRestoreColumn(eps->V,j+1,&vj1));
    if (j>k) {
      PetscCall(BVDotColumnBegin(eps->V,j+1,H+ldh*j));
      PetscCall(VecDotBegin(vj,vj,&dot));
      PetscCall(BVDotColumnEnd(eps->V,j+1,H+ldh*j));
      PetscCall(VecDotEnd(vj,vj,&dot));
      norm = PetscSqrtReal(PetscRealPart(dot));
      PetscCall(BVScaleColumn(eps->V,j,1.0/norm));
      H[ldh*(j-1)+j] = norm;
      for (i=0;i<j;i++) H[ldh*j+i] = H[ldh*j+i]/norm;
      H[ldh*j+j] = H[ldh*j+j]/dot;
      PetscCall(BVScaleColumn(eps->V,j+1,1.0/norm));
      *beta = norm;
    } else {  /* j==k */
      PetscCall(BVDotColumn(eps->V,j+1,H+ldh*j));
    }
    PetscCall(BVRestoreColumn(eps->V,j,&vj));
    PetscCall(BVMultColumn(eps->V,-1.0,1.0,j+1,H+ldh*j));
  }

  *breakdown = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSKrylovConvergence_Filter - Specialized version for STFILTER.
*/
static PetscErrorCode EPSKrylovConvergence_Filter(EPS eps,PetscBool getall,PetscInt kini,PetscInt nits,PetscReal beta,PetscReal gamma,PetscInt *kout)
{
  PetscInt       k,ninside,nconv;
  PetscScalar    re,im;
  PetscReal      resnorm;

  PetscFunctionBegin;
  ninside = 0;   /* count how many eigenvalues are located in the interval */
  for (k=kini;k<kini+nits;k++) {
    if (PetscRealPart(eps->eigr[k]) < gamma) break;
    ninside++;
  }
  if (eps->trackall) getall = PETSC_TRUE;
  eps->nev = ninside+kini;  /* adjust eigenvalue count */
  nconv = 0;   /* count how many eigenvalues satisfy the convergence criterion */
  for (k=kini;k<kini+(getall?nits:ninside);k++) {
    /* eigenvalue */
    re = eps->eigr[k];
    im = eps->eigi[k];
    PetscCall(DSVectors(eps->ds,DS_MAT_X,&k,&resnorm));
    resnorm *= beta;
    /* error estimate */
    PetscCall((*eps->converged)(eps,re,im,resnorm,&eps->errest[k],eps->convergedctx));
    if (eps->errest[k] < eps->tol) nconv++;
    else if (!getall) break;
  }
  *kout = kini+nconv;
  PetscCall(PetscInfo(eps,"Found %" PetscInt_FMT " eigenvalue approximations inside the interval (gamma=%g), k=%" PetscInt_FMT " nconv=%" PetscInt_FMT "\n",ninside,(double)gamma,k,nconv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSKrylovConvergence - Implements the loop that checks for convergence
   in Krylov methods.

   Input Parameters:
     eps   - the eigensolver; some error estimates are updated in eps->errest
     getall - whether all residuals must be computed
     kini  - initial value of k (the loop variable)
     nits  - number of iterations of the loop
     V     - set of basis vectors (used only if trueresidual is activated)
     nv    - number of vectors to process (dimension of Q, columns of V)
     beta  - norm of f (the residual vector of the Arnoldi/Lanczos factorization)
     corrf - correction factor for residual estimates (only in harmonic KS)

   Output Parameters:
     kout  - the first index where the convergence test failed
*/
PetscErrorCode EPSKrylovConvergence(EPS eps,PetscBool getall,PetscInt kini,PetscInt nits,PetscReal beta,PetscReal betat,PetscReal corrf,PetscInt *kout)
{
  PetscInt       k,newk,newk2,marker,ld,inside;
  PetscScalar    re,im,*Zr,*Zi,*X;
  PetscReal      resnorm,gamma,lerrest;
  PetscBool      isshift,isfilter,refined,istrivial;
  Vec            x=NULL,y=NULL,w[3];

  PetscFunctionBegin;
  if (PetscUnlikely(eps->which == EPS_ALL)) {
    PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilter));
    if (isfilter) {
      PetscCall(STFilterGetThreshold(eps->st,&gamma));
      PetscCall(EPSKrylovConvergence_Filter(eps,getall,kini,nits,beta,gamma,kout));
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }
  PetscCall(RGIsTrivial(eps->rg,&istrivial));
  if (PetscUnlikely(eps->trueres)) {
    PetscCall(BVCreateVec(eps->V,&x));
    PetscCall(BVCreateVec(eps->V,&y));
    PetscCall(BVCreateVec(eps->V,&w[0]));
    PetscCall(BVCreateVec(eps->V,&w[2]));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(BVCreateVec(eps->V,&w[1]));
#else
    w[1] = NULL;
#endif
  }
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(DSGetRefined(eps->ds,&refined));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isshift));
  marker = -1;
  if (eps->trackall) getall = PETSC_TRUE;
  for (k=kini;k<kini+nits;k++) {
    /* eigenvalue */
    re = eps->eigr[k];
    im = eps->eigi[k];
    if (!istrivial || eps->trueres || isshift || eps->conv==EPS_CONV_NORM) PetscCall(STBackTransform(eps->st,1,&re,&im));
    if (PetscUnlikely(!istrivial)) {
      PetscCall(RGCheckInside(eps->rg,1,&re,&im,&inside));
      if (marker==-1 && inside<0) marker = k;
      if (!(eps->trueres || isshift || eps->conv==EPS_CONV_NORM)) {  /* make sure eps->converged below uses the right value */
        re = eps->eigr[k];
        im = eps->eigi[k];
      }
    }
    newk = k;
    PetscCall(DSVectors(eps->ds,DS_MAT_X,&newk,&resnorm));
    if (PetscUnlikely(eps->trueres)) {
      PetscCall(DSGetArray(eps->ds,DS_MAT_X,&X));
      Zr = X+k*ld;
      if (newk==k+1) Zi = X+newk*ld;
      else Zi = NULL;
      PetscCall(EPSComputeRitzVector(eps,Zr,Zi,eps->V,x,y));
      PetscCall(DSRestoreArray(eps->ds,DS_MAT_X,&X));
      PetscCall(EPSComputeResidualNorm_Private(eps,PETSC_FALSE,re,im,x,y,w,&resnorm));
    }
    else if (!refined) resnorm *= beta*corrf;
    /* error estimate */
    PetscCall((*eps->converged)(eps,re,im,resnorm,&eps->errest[k],eps->convergedctx));
    if (marker==-1 && eps->errest[k] >= eps->tol) marker = k;
    if (PetscUnlikely(eps->twosided)) {
      newk2 = k;
      PetscCall(DSVectors(eps->ds,DS_MAT_Y,&newk2,&resnorm));
      resnorm *= betat;
      PetscCall((*eps->converged)(eps,re,im,resnorm,&lerrest,eps->convergedctx));
      eps->errest[k] = PetscMax(eps->errest[k],lerrest);
      if (marker==-1 && lerrest >= eps->tol) marker = k;
    }
    if (PetscUnlikely(newk==k+1)) {
      eps->errest[k+1] = eps->errest[k];
      k++;
    }
    if (marker!=-1 && !getall) break;
  }
  if (marker!=-1) k = marker;
  *kout = k;
  if (PetscUnlikely(eps->trueres)) {
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&y));
    PetscCall(VecDestroy(&w[0]));
    PetscCall(VecDestroy(&w[2]));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(VecDestroy(&w[1]));
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSPseudoLanczos(EPS eps,PetscReal *alpha,PetscReal *beta,PetscReal *omega,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscBool *symmlost,PetscReal *cos,Vec w)
{
  PetscInt       j,m = *M,i,ld,l;
  Vec            vj,vj1;
  PetscScalar    *hwork,lhwork[100];
  PetscReal      norm,norm1,norm2,t,sym=0.0,fro=0.0;
  PetscBLASInt   j_,one=1;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(DSGetDimensions(eps->ds,NULL,&l,NULL,NULL));
  if (cos) *cos = 1.0;
  if (m > 100) PetscCall(PetscMalloc1(m,&hwork));
  else hwork = lhwork;

  PetscCall(BVSetActiveColumns(eps->V,0,m));
  for (j=k;j<m;j++) {
    PetscCall(BVGetColumn(eps->V,j,&vj));
    PetscCall(BVGetColumn(eps->V,j+1,&vj1));
    PetscCall(STApply(eps->st,vj,vj1));
    PetscCall(BVRestoreColumn(eps->V,j,&vj));
    PetscCall(BVRestoreColumn(eps->V,j+1,&vj1));
    PetscCall(BVOrthogonalizeColumn(eps->V,j+1,hwork,&norm,breakdown));
    alpha[j] = PetscRealPart(hwork[j]);
    beta[j] = PetscAbsReal(norm);
    if (j==k) {
      PetscReal *f;

      PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&f));
      for (i=0;i<l;i++) hwork[i]  = 0.0;
      for (;i<j-1;i++)  hwork[i] -= f[2*ld+i];
      PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&f));
    }
    if (j>0) {
      hwork[j-1] -= beta[j-1];
      PetscCall(PetscBLASIntCast(j,&j_));
      sym = SlepcAbs(BLASnrm2_(&j_,hwork,&one),sym);
    }
    fro = SlepcAbs(fro,SlepcAbs(alpha[j],beta[j]));
    if (j>0) fro = SlepcAbs(fro,beta[j-1]);
    if (sym/fro>PetscMax(PETSC_SQRT_MACHINE_EPSILON,10*eps->tol)) { *symmlost = PETSC_TRUE; *M=j; break; }
    omega[j+1] = (norm<0.0)? -1.0: 1.0;
    PetscCall(BVScaleColumn(eps->V,j+1,1.0/norm));
    /* */
    if (cos) {
      PetscCall(BVGetColumn(eps->V,j+1,&vj1));
      PetscCall(VecNorm(vj1,NORM_2,&norm1));
      PetscCall(BVApplyMatrix(eps->V,vj1,w));
      PetscCall(BVRestoreColumn(eps->V,j+1,&vj1));
      PetscCall(VecNorm(w,NORM_2,&norm2));
      t = 1.0/(norm1*norm2);
      if (*cos>t) *cos = t;
    }
  }
  if (m > 100) PetscCall(PetscFree(hwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}
