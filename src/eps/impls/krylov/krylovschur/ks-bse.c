/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "krylovschur"

   Method: thick-restarted Lanczos for Bethe-Salpeter pseudo-Hermitan matrices

   References:

       [1] M. Shao et al, "A structure preserving Lanczos algorithm for computing
           the optical absorption spectrum", SIAM J. Matrix Anal. App. 39(2), 2018.

*/
#include <slepc/private/epsimpl.h>
#include "krylovschur.h"

static PetscErrorCode Orthog_Shao(Vec x,BV U,BV V,PetscInt j,PetscScalar *h,PetscScalar *c,PetscBool *breakdown)
{
  PetscInt i;

  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(U,0,j));
  PetscCall(BVSetActiveColumns(V,0,j));
  /* c = real(V^* x) ; c2 = imag(U^* x)*1i */
#if defined(PETSC_USE_COMPLEX)
  PetscCall(BVDotVecBegin(V,x,c));
  PetscCall(BVDotVecBegin(U,x,c+j));
  PetscCall(BVDotVecEnd(V,x,c));
  PetscCall(BVDotVecEnd(U,x,c+j));
#else
  PetscCall(BVDotVec(V,x,c));
#endif
#if defined(PETSC_USE_COMPLEX)
  for (i=0; i<j; i++) {
    c[i] = PetscRealPart(c[i]);
    c[j+i] = PetscCMPLX(0.0,PetscImaginaryPart(c[j+i]));
  }
#endif
  /* x = x-U*c-V*c2 */
  PetscCall(BVMultVec(U,-1.0,1.0,x,c));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(BVMultVec(V,-1.0,1.0,x,c+j));
#endif
  /* accumulate orthog coeffs into h */
  for (i=0; i<2*j; i++) h[i] += c[i];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Orthogonalize vector x against first j vectors in U and V
v is column j-1 of V */
static PetscErrorCode OrthogonalizeVector_Shao(Vec x,BV U,BV V,PetscInt j,Vec v,PetscReal *beta,PetscInt k,PetscScalar *h,PetscBool *breakdown)
{
  PetscReal alpha;
  PetscInt  i,l;

  PetscFunctionBegin;
  PetscCall(PetscArrayzero(h,2*j));

  /* Local orthogonalization */
  l = j==k+1?0:j-2;  /* 1st column to orthogonalize against */
  PetscCall(VecDotRealPart(x,v,&alpha));
  for (i=l; i<j-1; i++) h[i] = beta[i];
  h[j-1] = alpha;
  /* x = x-U(:,l:j-1)*h(l:j-1) */
  PetscCall(BVSetActiveColumns(U,l,j));
  PetscCall(BVMultVec(U,-1.0,1.0,x,h+l));

  /* Full orthogonalization */
  PetscCall(Orthog_Shao(x,U,V,j,h,h+2*j,breakdown));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSBSELanczos_Shao(EPS eps,BV U,BV V,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *M,PetscBool *breakdown)
{
  PetscInt       j,m = *M;
  Vec            v,x,y,w,f,g,vecs[2];
  Mat            H;
  IS             is[2];
  PetscReal      nrm;
  PetscScalar    *hwork,lhwork[100],gamma;

  PetscFunctionBegin;
  if (4*m > 100) PetscCall(PetscMalloc1(4*m,&hwork));
  else hwork = lhwork;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));

  /* create work vectors */
  PetscCall(BVGetColumn(V,0,&v));
  PetscCall(VecDuplicate(v,&w));
  vecs[0] = v;
  vecs[1] = w;
  PetscCall(VecCreateNest(PetscObjectComm((PetscObject)eps),2,is,vecs,&f));
  PetscCall(VecCreateNest(PetscObjectComm((PetscObject)eps),2,is,vecs,&g));
  PetscCall(BVRestoreColumn(V,0,&v));

  /* Normalize initial vector */
  if (k==0) {
    if (eps->nini==0) PetscCall(BVSetRandomColumn(eps->V,0));
    PetscCall(BVGetColumn(U,0,&x));
    PetscCall(BVGetColumn(V,0,&y));
    PetscCall(VecCopy(x,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecNestSetSubVec(f,0,x));
    PetscCall(VecNestSetSubVec(g,0,y));
    PetscCall(STApply(eps->st,f,g));
    PetscCall(VecDot(y,x,&gamma));
    nrm = PetscSqrtReal(PetscRealPart(gamma));
    PetscCall(VecScale(x,1.0/nrm));
    PetscCall(VecScale(y,1.0/nrm));
    PetscCall(BVRestoreColumn(U,0,&x));
    PetscCall(BVRestoreColumn(V,0,&y));
  }

  for (j=k;j<m;j++) {
    /* j+1 columns (indexes 0 to j) have been computed */
    PetscCall(BVGetColumn(V,j,&v));
    PetscCall(BVGetColumn(U,j+1,&x));
    PetscCall(BVGetColumn(V,j+1,&y));
    PetscCall(VecCopy(v,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecScale(w,-1.0));
    PetscCall(VecNestSetSubVec(f,0,v));
    PetscCall(VecNestSetSubVec(g,0,x));
    PetscCall(STApply(eps->st,f,g));
    PetscCall(OrthogonalizeVector_Shao(x,U,V,j+1,v,beta,k,hwork,breakdown));
    alpha[j] = PetscRealPart(hwork[j]);
    PetscCall(VecCopy(x,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecNestSetSubVec(f,0,x));
    PetscCall(VecNestSetSubVec(g,0,y));
    PetscCall(STApply(eps->st,f,g));
    PetscCall(VecDot(x,y,&gamma));
    beta[j] = PetscSqrtReal(PetscRealPart(gamma));
    PetscCall(VecScale(x,1.0/beta[j]));
    PetscCall(VecScale(y,1.0/beta[j]));
    PetscCall(BVRestoreColumn(V,j,&v));
    PetscCall(BVRestoreColumn(U,j+1,&x));
    PetscCall(BVRestoreColumn(V,j+1,&y));
  }
  if (4*m > 100) PetscCall(PetscFree(hwork));
  PetscCall(VecDestroy(&w));
  PetscCall(VecDestroy(&f));
  PetscCall(VecDestroy(&g));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_BSE_Shao(EPS eps)
{
  Mat         H;
  Vec         u1,v1;
  BV          U,V;
  IS          is[2];
  PetscInt    k;
  PetscScalar lambda;

  PetscFunctionBegin;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));
  for (k=0; k<eps->nconv; k++) {
    PetscCall(BVGetColumn(U,k,&u1));
    PetscCall(BVGetColumn(V,k,&v1));
    /* approx eigenvector is [    (eigr[k]*u1+v1)]
                             [conj(eigr[k]*u1-v1)]  */
    lambda = eps->eigr[k];
    PetscCall(STBackTransform(eps->st,1,&lambda,&eps->eigi[k]));
    PetscCall(VecAYPX(u1,lambda,v1));
    PetscCall(VecAYPX(v1,-2.0,u1));
    PetscCall(VecConjugate(v1));
    PetscCall(BVRestoreColumn(U,k,&u1));
    PetscCall(BVRestoreColumn(V,k,&v1));
  }
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
  /* Normalize eigenvectors */
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
  PetscCall(BVNormalize(eps->V,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode Orthog_Gruning(Vec x,BV U,BV V,BV HU,BV HV,PetscInt j,PetscScalar *h,PetscScalar *c,PetscBool s,PetscBool *breakdown)
{
  PetscInt i;

  PetscFunctionBegin;
  PetscCall(BVSetActiveColumns(U,0,j));
  PetscCall(BVSetActiveColumns(HU,0,j));
  if (s) {
    PetscCall(BVSetActiveColumns(V,0,j));
    PetscCall(BVSetActiveColumns(HV,0,j));
  } else {
    PetscCall(BVSetActiveColumns(V,0,j-1));
    PetscCall(BVSetActiveColumns(HV,0,j-1));
  }
#if defined(PETSC_USE_COMPLEX)
  PetscCall(BVDotVecBegin(HU,x,c));
  if (s || j>1) PetscCall(BVDotVecBegin(HV,x,c+j));
  PetscCall(BVDotVecEnd(HU,x,c));
  if (s || j>1) PetscCall(BVDotVecEnd(HV,x,c+j));
#else
  if (s) PetscCall(BVDotVec(HU,x,c));
  else PetscCall(BVDotVec(HV,x,c+j));
#endif
  for (i=0; i<j; i++) {
    if (s) {   /* c1 = 2*real(HU^* x) ; c2 = 2*imag(HV^* x)*1i */
#if defined(PETSC_USE_COMPLEX)
      c[i] = PetscRealPart(c[i]);
      c[j+i] = PetscCMPLX(0.0,PetscImaginaryPart(c[j+i]));
#else
      c[j+i] = 0.0;
#endif
    } else {   /* c1 = 2*imag(HU^* x)*1i ; c2 = 2*real(HV^* x) */
#if defined(PETSC_USE_COMPLEX)
      c[i] = PetscCMPLX(0.0,PetscImaginaryPart(c[i]));
      c[j+i] = PetscRealPart(c[j+i]);
#else
      c[i] = 0.0;
#endif
    }
  }
  /* x = x-U*c1-V*c2 */
#if defined(PETSC_USE_COMPLEX)
  PetscCall(BVMultVec(U,-2.0,1.0,x,c));
  PetscCall(BVMultVec(V,-2.0,1.0,x,c+j));
#else
  if (s) PetscCall(BVMultVec(U,-2.0,1.0,x,c));
  else PetscCall(BVMultVec(V,-2.0,1.0,x,c+j));
#endif
  /* accumulate orthog coeffs into h */
  for (i=0; i<2*j; i++) h[i] += 2*c[i];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Orthogonalize vector x against first j vectors in U and V */
static PetscErrorCode OrthogonalizeVector_Gruning(Vec x,BV U,BV V,BV HU,BV HV,PetscInt j,PetscReal *beta,PetscInt k,PetscScalar *h,PetscBool s,PetscBool *breakdown)
{
  PetscInt l,i;
  Vec      u;

  PetscFunctionBegin;
  PetscCall(PetscArrayzero(h,4*j));

  /* Local orthogonalization */
  if (s) {
    PetscCall(BVGetColumn(U,j-1,&u));
    PetscCall(VecAXPY(x,-*beta,u));
    PetscCall(BVRestoreColumn(U,j-1,&u));
    h[j-1] = *beta;
  } else {
    l = j==k+1?0:j-2;  /* 1st column to orthogonalize against */
    for (i=l; i<j-1; i++) h[j+i] = beta[i];
    /* x = x-V(:,l:j-2)*h(l:j-2) */
    PetscCall(BVSetActiveColumns(V,l,j-1));
    PetscCall(BVMultVec(V,-1.0,1.0,x,h+j+l));
  }

  /* Full orthogonalization */
  PetscCall(Orthog_Gruning(x,U,V,HU,HV,j,h,h+2*j,s,breakdown));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSBSELanczos_Gruning(EPS eps,BV U,BV V,BV HU,BV HV,PetscReal *beta1,PetscReal *beta2,PetscInt k,PetscInt *M,PetscBool *breakdown)
{
  PetscInt       j,m = *M;
  Vec            v,x,y,w,f,g,vecs[2];
  Mat            H;
  IS             is[2];
  PetscReal      nrm;
  PetscScalar    *hwork,lhwork[100],dot;

  PetscFunctionBegin;
  if (4*m > 100) PetscCall(PetscMalloc1(4*m,&hwork));
  else hwork = lhwork;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));

  /* create work vectors */
  PetscCall(BVGetColumn(V,0,&v));
  PetscCall(VecDuplicate(v,&w));
  vecs[0] = v;
  vecs[1] = w;
  PetscCall(VecCreateNest(PetscObjectComm((PetscObject)eps),2,is,vecs,&f));
  PetscCall(VecCreateNest(PetscObjectComm((PetscObject)eps),2,is,vecs,&g));
  PetscCall(BVRestoreColumn(V,0,&v));

  /* Normalize initial vector */
  if (k==0) {
    if (eps->nini==0) PetscCall(BVSetRandomColumn(eps->V,0));
    /* y = Hmult(v1,1) */
    PetscCall(BVGetColumn(U,k,&x));
    PetscCall(BVGetColumn(HU,k,&y));
    PetscCall(VecCopy(x,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecNestSetSubVec(f,0,x));
    PetscCall(VecNestSetSubVec(g,0,y));
    PetscCall(STApply(eps->st,f,g));
    /* nrm = sqrt(2*real(u1'*y)); */
    PetscCall(VecDot(x,y,&dot));
    nrm = PetscSqrtReal(PetscRealPart(2*dot));
    /* U(:,j) = u1/nrm; */
    /* HU(:,j) = y/nrm; */
    PetscCall(VecScale(x,1.0/nrm));
    PetscCall(VecScale(y,1.0/nrm));
    PetscCall(BVRestoreColumn(U,k,&x));
    PetscCall(BVRestoreColumn(HU,k,&y));
  }

  for (j=k;j<m;j++) {
    /* j+1 columns (indexes 0 to j) have been computed */
    PetscCall(BVGetColumn(HU,j,&x));
    PetscCall(BVGetColumn(V,j,&v));
    PetscCall(BVGetColumn(HV,j,&y));
    PetscCall(VecCopy(x,v));
    PetscCall(BVRestoreColumn(HU,j,&x));
    /* v = Orthogonalize HU(:,j) */
    PetscCall(OrthogonalizeVector_Gruning(v,U,V,HU,HV,j+1,beta2,k,hwork,PETSC_FALSE,breakdown));
    /* y = Hmult(v,-1) */
    PetscCall(VecCopy(v,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecScale(w,-1.0));
    PetscCall(VecNestSetSubVec(f,0,v));
    PetscCall(VecNestSetSubVec(g,0,y));
    PetscCall(STApply(eps->st,f,g));
    /* beta = sqrt(2*real(v'*y)); */
    PetscCall(VecDot(v,y,&dot));
    beta1[j] = PetscSqrtReal(PetscRealPart(2*dot)); //FIXME Check beta != 0?
    /* V(:,j) = v/beta1; */
    /* HV(:,j) = y/beta1; */
    PetscCall(VecScale(v,1.0/beta1[j]));
    PetscCall(VecScale(y,1.0/beta1[j]));
    PetscCall(BVRestoreColumn(V,j,&v));
    PetscCall(BVRestoreColumn(HV,j,&y));

    PetscCall(BVGetColumn(HV,j,&x));
    PetscCall(BVGetColumn(U,j+1,&v));
    PetscCall(BVGetColumn(HU,j+1,&y));
    PetscCall(VecCopy(x,v));
    PetscCall(BVRestoreColumn(HV,j,&x));
    /* v = Orthogonalize HV(:,j) */
    PetscCall(OrthogonalizeVector_Gruning(v,U,V,HU,HV,j+1,&beta1[j],k,hwork,PETSC_TRUE,breakdown));
    /* y = Hmult(v,1) */
    PetscCall(VecCopy(v,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecNestSetSubVec(f,0,v));
    PetscCall(VecNestSetSubVec(g,0,y));
    PetscCall(STApply(eps->st,f,g));
    /* beta = sqrt(2*real(v'*y)); */
    PetscCall(VecDot(v,y,&dot));
    beta2[j] = PetscSqrtReal(PetscRealPart(2*dot)); //FIXME Check beta != 0?
    /* U(:,j) = v/beta2; */
    /* HU(:,j) = y/beta2; */
    PetscCall(VecScale(v,1.0/beta2[j]));
    PetscCall(VecScale(y,1.0/beta2[j]));
    PetscCall(BVRestoreColumn(U,j+1,&v));
    PetscCall(BVRestoreColumn(HU,j+1,&y));
  }
  if (4*m > 100) PetscCall(PetscFree(hwork));
  PetscCall(VecDestroy(&w));
  PetscCall(VecDestroy(&f));
  PetscCall(VecDestroy(&g));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_BSE_Gruning(EPS eps)
{
  Mat         H;
  Vec         u1,v1;
  BV          U,V;
  IS          is[2];
  PetscInt    k;

  PetscFunctionBegin;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));
  /* approx eigenvector [x1] is [     u1+v1       ]
                        [x2]    [conj(u1)-conj(v1)]  */
  for (k=0; k<eps->nconv; k++) {
    PetscCall(BVGetColumn(U,k,&u1));
    PetscCall(BVGetColumn(V,k,&v1));
    /* x1 = u1 + v1 */
    PetscCall(VecAXPY(u1,1.0,v1));
    /* x2 = conj(u1) - conj(v1) = conj(u1 - v1) = conj((u1 + v1) - 2*v1) */
    PetscCall(VecAYPX(v1,-2.0,u1));
    PetscCall(VecConjugate(v1));
    PetscCall(BVRestoreColumn(U,k,&u1));
    PetscCall(BVRestoreColumn(V,k,&v1));
  }
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
  /* Normalize eigenvectors */
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
  PetscCall(BVNormalize(eps->V,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode Orthog_ProjectedBSE(Vec hx,Vec hy,BV X,BV Y,PetscInt j,PetscScalar *h,PetscScalar *c,PetscBool *breakdown)
{
  PetscInt          i;
  Mat               MX,MY,MXl,MYl;
  Vec               c1,c2,hxl,hyl,hz;
  PetscScalar       *cx1,*cx2;
  PetscMPIInt       len;

  PetscFunctionBegin;
  PetscCall(PetscMPIIntCast(j,&len));
  PetscCall(BVSetActiveColumns(X,0,j));
  PetscCall(BVSetActiveColumns(Y,0,j));
  /* BVTDotVec does not exist yet, implemented via MatMult operations */
  PetscCall(BVGetMat(X,&MX));
  PetscCall(BVGetMat(Y,&MY));
  PetscCall(MatDenseGetLocalMatrix(MX,&MXl));
  PetscCall(MatDenseGetLocalMatrix(MY,&MYl));
  PetscCall(MatCreateVecs(MXl,&c1,&hyl));
  PetscCall(MatCreateVecs(MXl,&c2,&hxl));

  /* c1 =  X^* hx - Y^* hy
   * c2 = -Y^T hx + X^T hy */

  PetscCall(VecGetLocalVector(hx,hxl));
  PetscCall(VecGetLocalVector(hy,hyl));
  PetscCall(VecDuplicate(hx,&hz));
  /* c1 = -(Y^* hy) */
  PetscCall(MatMultHermitianTranspose(MYl,hyl,c1));
  PetscCall(VecScale(c1,-1));
  /* c1 = c1 + X^* hx */
  PetscCall(MatMultHermitianTransposeAdd(MXl,hxl,c1,c1));
  PetscCall(VecGetArray(c1,&cx1));
  PetscCallMPI(MPIU_Allreduce(MPI_IN_PLACE,cx1,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)hx)));
  PetscCall(VecRestoreArray(c1,&cx1));
  /* c2 = -(Y^T hx) */
  PetscCall(MatMultTranspose(MYl,hxl,c2));
  PetscCall(VecScale(c2,-1));
  /* c2 = c2 + X^T hy */
  PetscCall(MatMultTransposeAdd(MXl,hyl,c2,c2));
  PetscCall(VecGetArray(c2,&cx2));
  PetscCallMPI(MPIU_Allreduce(MPI_IN_PLACE,cx2,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)hx)));
  PetscCall(VecRestoreArray(c2,&cx2));
  PetscCall(VecRestoreLocalVector(hx,hxl));
  PetscCall(VecRestoreLocalVector(hy,hyl));

  /* accumulate orthog coeffs into h */
  PetscCall(VecGetArrayRead(c1,(const PetscScalar**)&cx1));
  PetscCall(VecGetArrayRead(c2,(const PetscScalar**)&cx2));
  for (i=0; i<j; i++) h[i] += cx1[i];
  for (i=0; i<j; i++) h[i+j] += cx2[i];
  PetscCall(VecRestoreArrayRead(c1,(const PetscScalar**)&cx1));
  PetscCall(VecRestoreArrayRead(c2,(const PetscScalar**)&cx2));

  /* u = hx - X c1 - conj(Y) c2 */

  /* conj(Y) c2 */
  PetscCall(VecConjugate(c2));
  PetscCall(VecGetLocalVector(hz,hxl));
  PetscCall(MatMult(MYl,c2,hxl));
  PetscCall(VecConjugate(hxl));
  /* X c1 */
  PetscCall(MatMultAdd(MXl,c1,hxl,hxl));
  PetscCall(VecRestoreLocalVector(hz,hxl));
  PetscCall(VecAXPY(hx,-1,hz));

  PetscCall(BVRestoreMat(X,&MX));
  PetscCall(BVRestoreMat(Y,&MY));
  PetscCall(VecDestroy(&c1));
  PetscCall(VecDestroy(&c2));
  PetscCall(VecDestroy(&hxl));
  PetscCall(VecDestroy(&hyl));
  PetscCall(VecDestroy(&hz));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Orthogonalize vector x against first j vectors in X and Y */
static PetscErrorCode OrthogonalizeVector_ProjectedBSE(Vec hx,Vec hy,BV X,BV Y,PetscInt j,PetscReal *beta,PetscInt k,PetscScalar *h,PetscBool *breakdown)
{
  PetscInt    l,i;
  PetscScalar alpha,alpha1,alpha2;
  Vec         x,y;

  PetscFunctionBegin;
  PetscCall(PetscArrayzero(h,2*j));

  /* Local orthogonalization */
  l = j==k+1?0:j-2;  /* 1st column to orthogonalize against */
  /* alpha = X(:,j-1)'*hx-Y(:,j-1)'*hy */
  PetscCall(BVGetColumn(X,j-1,&x));
  PetscCall(BVGetColumn(Y,j-1,&y));
  PetscCall(VecDotBegin(hx,x,&alpha1));
  PetscCall(VecDotBegin(hy,y,&alpha2));
  PetscCall(VecDotEnd(hx,x,&alpha1));
  PetscCall(VecDotEnd(hy,y,&alpha2));
  PetscCall(BVRestoreColumn(X,j-1,&x));
  PetscCall(BVRestoreColumn(Y,j-1,&y));
  alpha = alpha1-alpha2;
  /* Store coeffs into h */
  for (i=l; i<j-1; i++) h[i] = h[j+i] = beta[i];
  h[j-1] = alpha;
  h[2*j-1] = alpha-1.0;
  /* Orthogonalize: hx = hx - X(:,l:j-1)*h1 - conj(Y(:,l:j-1))*h2 */
  PetscCall(BVSetActiveColumns(X,l,j));
  PetscCall(BVSetActiveColumns(Y,l,j));
  PetscCall(BVMultVec(X,-1.0,1.0,hx,h+l));
  PetscCall(VecConjugate(hx));
  for (i=j+l; i<2*j; i++) h[j+i] = PetscConj(h[i]);
  PetscCall(BVMultVec(Y,-1.0,1.0,hx,h+2*j+l));
  PetscCall(VecConjugate(hx));

  /* Full orthogonalization */
  PetscCall(VecCopy(hx,hy));
  PetscCall(VecConjugate(hy));
  PetscCall(Orthog_ProjectedBSE(hx,hy,X,Y,j,h,h+2*j,breakdown));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSBSELanczos_ProjectedBSE(EPS eps,BV X,BV Y,Vec v,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *M,PetscBool *breakdown)
{
  PetscInt       j,m = *M;
  Vec            u,x,y,w,f,g,vecs[2];
  Mat            H;
  IS             is[2];
  PetscReal      nrm;
  PetscScalar    *hwork,lhwork[100],gamma;

  PetscFunctionBegin;
  if (4*m > 100) PetscCall(PetscMalloc1(4*m,&hwork));
  else hwork = lhwork;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));

  /* create work vectors */
  PetscCall(BVGetColumn(Y,0,&u));
  PetscCall(VecDuplicate(u,&w));
  vecs[0] = u;
  vecs[1] = w;
  PetscCall(VecCreateNest(PetscObjectComm((PetscObject)eps),2,is,vecs,&f));
  PetscCall(VecCreateNest(PetscObjectComm((PetscObject)eps),2,is,vecs,&g));
  PetscCall(BVRestoreColumn(Y,0,&u));

  /* Normalize initial vector */
  if (k==0) {
    if (eps->nini==0) PetscCall(BVSetRandomColumn(eps->V,0));
    PetscCall(BVGetColumn(X,0,&x));
    /* v = Hmult(u,1) */
    PetscCall(BVGetColumn(Y,0,&y));
    PetscCall(VecCopy(x,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecNestSetSubVec(f,0,x));
    PetscCall(VecNestSetSubVec(g,0,y));
    PetscCall(STApply(eps->st,f,g));
    /* nrm = sqrt(real(u'v)) */
    PetscCall(VecDot(y,x,&gamma));
    nrm = PetscSqrtReal(PetscRealPart(gamma));
    /* u = u /(nrm*2) */
    PetscCall(VecScale(x,1.0/(2.0*nrm)));
    /* v = v /(nrm*2) */
    PetscCall(VecScale(y,1.0/(2.0*nrm)));
    PetscCall(VecCopy(y,v));
    /* X(:,1) = (u+v) */
    PetscCall(VecAXPY(x,1,y));
    /* Y(:,1) = conj(u-v) */
    PetscCall(VecAYPX(y,-2,x));
    PetscCall(VecConjugate(y));
    PetscCall(BVRestoreColumn(X,0,&x));
    PetscCall(BVRestoreColumn(Y,0,&y));
  }

  for (j=k;j<m;j++) {
    /* j+1 columns (indexes 0 to j) have been computed */
    PetscCall(BVGetColumn(X,j+1,&x));
    PetscCall(BVGetColumn(Y,j+1,&y));
    /* u = Hmult(v,-1)*/
    PetscCall(VecCopy(v,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecScale(w,-1.0));
    PetscCall(VecNestSetSubVec(f,0,v));
    PetscCall(VecNestSetSubVec(g,0,x));
    PetscCall(STApply(eps->st,f,g));
    /* hx = (u+v) */
    PetscCall(VecCopy(x,y));
    PetscCall(VecAXPY(x,1,v));
    /* hy = conj(u-v) */
    PetscCall(VecAXPY(y,-1,v));
    PetscCall(VecConjugate(y));
    /* [u,cd] = orthog(hx,hy,X(:,1:j),Y(:,1:j),opt)*/
    PetscCall(OrthogonalizeVector_ProjectedBSE(x,y,X,Y,j+1,beta,k,hwork,breakdown));
    /* alpha(j) = real(cd(j))-1/2 */
    alpha[j] = 2*(PetscRealPart(hwork[j]) - 0.5);
    /* v = Hmult(u,1) */
    PetscCall(VecCopy(x,w));
    PetscCall(VecConjugate(w));
    PetscCall(VecNestSetSubVec(f,0,x));
    PetscCall(VecNestSetSubVec(g,0,y));
    PetscCall(STApply(eps->st,f,g));
    /* nrm = sqrt(real(u'*v)) */
    /* beta(j) = nrm */
    PetscCall(VecDot(x,y,&gamma));
    beta[j] = 2.0*PetscSqrtReal(PetscRealPart(gamma));
    /* u = u/(nrm*2) */
    PetscCall(VecScale(x,1.0/beta[j]));
    /* v = v/(nrm*2) */
    PetscCall(VecScale(y,1.0/beta[j]));
    PetscCall(VecCopy(y,v));
    /* X(:,j+1) = (u+v) */
    PetscCall(VecAXPY(x,1,y));
    /* Y(:,j+1) = conj(u-v) */
    PetscCall(VecAYPX(y,-2,x));
    PetscCall(VecConjugate(y));
    /* Restore */
    PetscCall(BVRestoreColumn(X,j+1,&x));
    PetscCall(BVRestoreColumn(Y,j+1,&y));
  }
  if (4*m > 100) PetscCall(PetscFree(hwork));
  PetscCall(VecDestroy(&w));
  PetscCall(VecDestroy(&f));
  PetscCall(VecDestroy(&g));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_BSE_ProjectedBSE(EPS eps)
{
  Mat         H;
  Vec         x1,y1,cx1,cy1;
  BV          X,Y;
  IS          is[2];
  PetscInt    k;
  PetscScalar delta1,delta2,lambda;

  PetscFunctionBegin;
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&X,&Y));
  PetscCall(BVCreateVec(X,&cx1));
  PetscCall(BVCreateVec(Y,&cy1));
  for (k=0; k<eps->nconv; k++) {
    /* approx eigenvector is [ delta1*x1 + delta2*conj(y1) ]
                             [ delta1*y1 + delta2*conj(x1) ] */
    lambda = eps->eigr[k];
    PetscCall(STBackTransform(eps->st,1,&lambda,&eps->eigi[k]));
    delta1 = lambda+1.0;
    delta2 = lambda-1.0;
    PetscCall(BVGetColumn(X,k,&x1));
    PetscCall(BVGetColumn(Y,k,&y1));
    PetscCall(VecCopy(x1,cx1));
    PetscCall(VecCopy(y1,cy1));
    PetscCall(VecConjugate(cx1));
    PetscCall(VecConjugate(cy1));
    PetscCall(VecScale(x1,delta1));
    PetscCall(VecScale(y1,delta1));
    PetscCall(VecAXPY(x1,delta2,cy1));
    PetscCall(VecAXPY(y1,delta2,cx1));
    PetscCall(BVRestoreColumn(X,k,&x1));
    PetscCall(BVRestoreColumn(Y,k,&y1));
  }
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&X,&Y));
  PetscCall(BVSetActiveColumns(eps->V,0,eps->nconv));
  /* Normalize eigenvector */
  PetscCall(BVNormalize(eps->V,NULL));
  PetscCall(VecDestroy(&cx1));
  PetscCall(VecDestroy(&cy1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSetUp_KrylovSchur_BSE(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscBool       flg,sinvert;
  PetscInt        nev;

  PetscFunctionBegin;
  PetscCheck((eps->problem_type==EPS_BSE),PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Problem type should be BSE");
  EPSCheckUnsupportedCondition(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_EXTRACTION | EPS_FEATURE_BALANCE,PETSC_TRUE," with BSE structure");
  if (eps->nev==0 && eps->stop!=EPS_STOP_THRESHOLD) eps->nev = 1;
  nev = (eps->nev+1)/2;
  PetscCall(EPSSetDimensions_Default(eps,&nev,&eps->ncv,&eps->mpd));
  PetscCheck(eps->ncv<=nev+eps->mpd,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv)*((eps->stop==EPS_STOP_THRESHOLD)?10:1);

  PetscCall(PetscObjectTypeCompareAny((PetscObject)eps->st,&flg,STSINVERT,STSHIFT,""));
  PetscCheck(flg,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Krylov-Schur BSE only supports shift and shift-and-invert ST");
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&sinvert));
  if (!eps->which) {
    if (sinvert) eps->which = EPS_TARGET_MAGNITUDE;
    else eps->which = EPS_SMALLEST_MAGNITUDE;
  }

  if (!ctx->keep) ctx->keep = 0.5;
  PetscCall(STSetStructured(eps->st,PETSC_TRUE));

  PetscCall(EPSAllocateSolution(eps,1));
  switch (ctx->bse) {
    case EPS_KRYLOVSCHUR_BSE_SHAO:
      eps->ops->solve = EPSSolve_KrylovSchur_BSE_Shao;
      eps->ops->computevectors = EPSComputeVectors_BSE_Shao;
      PetscCall(DSSetType(eps->ds,DSHEP));
      PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      break;
    case EPS_KRYLOVSCHUR_BSE_GRUNING:
      eps->ops->solve = EPSSolve_KrylovSchur_BSE_Gruning;
      eps->ops->computevectors = EPSComputeVectors_BSE_Gruning;
      PetscCall(DSSetType(eps->ds,DSSVD));
      PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      break;
    case EPS_KRYLOVSCHUR_BSE_PROJECTEDBSE:
      eps->ops->solve = EPSSolve_KrylovSchur_BSE_ProjectedBSE;
      eps->ops->computevectors = EPSComputeVectors_BSE_ProjectedBSE;
      PetscCall(DSSetType(eps->ds,DSHEP));
      PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
      PetscCall(DSSetExtraRow(eps->ds,PETSC_TRUE));
      PetscCall(DSAllocate(eps->ds,eps->ncv+1));
      break;
    default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Unexpected error");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_BSE_Shao(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,k,l,ld,nv,nconv=0,nevsave;
  Mat             H,Q;
  BV              U,V;
  IS              is[2];
  PetscReal       *a,*b,beta;
  PetscBool       breakdown=PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));

  /* Extract matrix blocks */
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));

  /* Get the split bases */
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));

  nevsave  = eps->nev;
  eps->nev = (eps->nev+1)/2;
  l = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&a));
    b = a + ld;
    PetscCall(EPSBSELanczos_Shao(eps,U,V,a,b,eps->nconv+l,&nv,&breakdown));
    beta = b[nv-1];
    PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&a));
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSSetState(eps->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(eps->ds));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

    /* Check convergence */
    for (i=0;i<nv;i++) eps->eigr[i] = PetscSqrtReal(PetscRealPart(eps->eigr[i]));
    PetscCall(EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta,0.0,1.0,&k));
    EPSSetCtxThreshold(eps,eps->eigr,eps->eigi,k);
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(eps,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      PetscCheck(!breakdown,PetscObjectComm((PetscObject)eps),PETSC_ERR_CONV_FAILED,"Breakdown in BSE Krylov-Schur (beta=%g)",(double)beta);
      /* Prepare the Rayleigh quotient for restart */
      PetscCall(DSTruncate(eps->ds,k+l,PETSC_FALSE));
    }
    /* Update the corresponding vectors
       U(:,idx) = U*Q(:,idx),  V(:,idx) = V*Q(:,idx) */
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&Q));
    PetscCall(BVMultInPlace(U,Q,eps->nconv,k+l));
    PetscCall(BVMultInPlace(V,Q,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&Q));

    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      PetscCall(BVCopyColumn(eps->V,nv,k+l));
      if (eps->stop==EPS_STOP_THRESHOLD && nv-k<5) {  /* reallocate */
        eps->ncv = eps->mpd+k;
        PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
        PetscCall(EPSReallocateSolution(eps,eps->ncv+1));
        PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));
        for (i=nv;i<eps->ncv;i++) eps->perm[i] = i;
        PetscCall(DSReallocate(eps->ds,eps->ncv+1));
        PetscCall(DSGetLeadingDimension(eps->ds,&ld));
      }
    }
    eps->nconv = k;
    PetscCall(EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv));
  }

  eps->nev = nevsave;

  PetscCall(DSTruncate(eps->ds,eps->nconv,PETSC_TRUE));
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSConvergence_Gruning - convergence check based on SVDKrylovConvergence().
*/
static PetscErrorCode EPSConvergence_Gruning(EPS eps,PetscBool getall,PetscInt kini,PetscInt nits,PetscInt *kout)
{
  PetscInt       k,marker,ld;
  PetscReal      *alpha,*beta,resnorm;
  PetscBool      extra;

  PetscFunctionBegin;
  *kout = 0;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(DSGetExtraRow(eps->ds,&extra));
  PetscCheck(extra,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Only implemented for DS with extra row");
  marker = -1;
  if (eps->trackall) getall = PETSC_TRUE;
  PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&alpha));
  beta = alpha + ld;
  for (k=kini;k<kini+nits;k++) {
    resnorm = PetscAbsReal(beta[k]);
    PetscCall((*eps->converged)(eps,eps->eigr[k],eps->eigi[k],resnorm,&eps->errest[k],eps->convergedctx));
    if (marker==-1 && eps->errest[k] >= eps->tol) marker = k;
    if (marker!=-1 && !getall) break;
  }
  PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&alpha));
  if (marker!=-1) k = marker;
  *kout = k;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_BSE_Gruning(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,k,l,ld,nv,nconv=0,nevsave;
  Mat             H,Q,Z;
  BV              U,V,HU,HV;
  IS              is[2];
  PetscReal       *d1,*d2,beta;
  PetscBool       breakdown=PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));

  /* Extract matrix blocks */
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));

  /* Get the split bases */
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));

  /* Create HU and HV */
  PetscCall(BVDuplicate(U,&HU));
  PetscCall(BVDuplicate(V,&HV));

  nevsave  = eps->nev;
  eps->nev = (eps->nev+1)/2;
  l = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&d1));
    d2 = d1 + ld;
    PetscCall(EPSBSELanczos_Gruning(eps,U,V,HU,HV,d1,d2,eps->nconv+l,&nv,&breakdown));
    beta = d1[nv-1];
    PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&d1));

    /* Compute SVD */
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSSVDSetDimensions(eps->ds,nv));
    PetscCall(DSSetState(eps->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));

    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(eps->ds));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

    /* Check convergence */
    PetscCall(EPSConvergence_Gruning(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,&k));
    EPSSetCtxThreshold(eps,eps->eigr,eps->eigi,k);
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(eps,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      PetscCheck(!breakdown,PetscObjectComm((PetscObject)eps),PETSC_ERR_CONV_FAILED,"Breakdown in BSE Krylov-Schur (beta=%g)",(double)beta);
      /* Prepare the Rayleigh quotient for restart */
      PetscCall(DSTruncate(eps->ds,k+l,PETSC_FALSE));
    }
    /* Update the corresponding vectors
       U(:,idx) = U*Q(:,idx),  V(:,idx) = V*Z(:,idx) */
    PetscCall(DSGetMat(eps->ds,DS_MAT_U,&Q));
    PetscCall(DSGetMat(eps->ds,DS_MAT_V,&Z));
    PetscCall(BVMultInPlace(U,Z,eps->nconv,k+l));
    PetscCall(BVMultInPlace(V,Q,eps->nconv,k+l));
    PetscCall(BVMultInPlace(HU,Z,eps->nconv,k+l));
    PetscCall(BVMultInPlace(HV,Q,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_U,&Q));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_V,&Z));

    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      PetscCall(BVCopyColumn(U,nv,k+l));
      PetscCall(BVCopyColumn(HU,nv,k+l));
      if (eps->stop==EPS_STOP_THRESHOLD && nv-k<5) {  /* reallocate */
        eps->ncv = eps->mpd+k;
        PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
        PetscCall(EPSReallocateSolution(eps,eps->ncv+1));
        PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));
        PetscCall(BVResize(HU,eps->ncv+1,PETSC_TRUE));
        PetscCall(BVResize(HV,eps->ncv+1,PETSC_TRUE));
        for (i=nv;i<eps->ncv;i++) { eps->perm[i] = i; eps->eigi[i] = 0.0; }
        PetscCall(DSReallocate(eps->ds,eps->ncv+1));
        PetscCall(DSGetLeadingDimension(eps->ds,&ld));
      }
    }
    eps->nconv = k;
    PetscCall(EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv));
  }

  eps->nev = nevsave;

  PetscCall(DSTruncate(eps->ds,eps->nconv,PETSC_TRUE));
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));

  PetscCall(BVDestroy(&HU));
  PetscCall(BVDestroy(&HV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_BSE_ProjectedBSE(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,k,l,ld,nv,nconv=0,nevsave;
  Mat             H,Q;
  Vec             v;
  BV              U,V;
  IS              is[2];
  PetscReal       *a,*b,beta;
  PetscBool       breakdown=PETSC_FALSE;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));

  /* Extract matrix blocks */
  PetscCall(STGetMatrix(eps->st,0,&H));
  PetscCall(MatNestGetISs(H,is,NULL));

  /* Get the split bases */
  PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));

  /* Vector v */
  PetscCall(BVCreateVec(V,&v));

  nevsave  = eps->nev;
  eps->nev = (eps->nev+1)/2;
  l = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&a));
    b = a + ld;
    PetscCall(EPSBSELanczos_ProjectedBSE(eps,U,V,v,a,b,eps->nconv+l,&nv,&breakdown));
    beta = b[nv-1];
    PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&a));
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSSetState(eps->ds,l?DS_STATE_RAW:DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,eps->eigr,eps->eigi));
    PetscCall(DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(eps->ds));
    PetscCall(DSSynchronize(eps->ds,eps->eigr,eps->eigi));

    /* Check convergence */
    for (i=0;i<nv;i++) eps->eigr[i] = PetscSqrtReal(PetscRealPart(eps->eigr[i]));
    PetscCall(EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta,0.0,1.0,&k));
    EPSSetCtxThreshold(eps,eps->eigr,eps->eigi,k);
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx));
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(eps,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      PetscCheck(!breakdown,PetscObjectComm((PetscObject)eps),PETSC_ERR_CONV_FAILED,"Breakdown in BSE Krylov-Schur (beta=%g)",(double)beta);
      /* Prepare the Rayleigh quotient for restart */
      PetscCall(DSTruncate(eps->ds,k+l,PETSC_FALSE));
    }
    /* Update the corresponding vectors
       U(:,idx) = U*Q(:,idx),  V(:,idx) = V*Q(:,idx) */
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&Q));
    PetscCall(BVMultInPlace(U,Q,eps->nconv,k+l));
    PetscCall(BVMultInPlace(V,Q,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&Q));

    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      PetscCall(BVCopyColumn(eps->V,nv,k+l));
      if (eps->stop==EPS_STOP_THRESHOLD && nv-k<5) {  /* reallocate */
        eps->ncv = eps->mpd+k;
        PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));
        PetscCall(EPSReallocateSolution(eps,eps->ncv+1));
        PetscCall(BVGetSplitRows(eps->V,is[0],is[1],&U,&V));
        for (i=nv;i<eps->ncv;i++) eps->perm[i] = i;
        PetscCall(DSReallocate(eps->ds,eps->ncv+1));
        PetscCall(DSGetLeadingDimension(eps->ds,&ld));
      }
    }
    eps->nconv = k;
    PetscCall(EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv));
  }

  eps->nev = nevsave;

  PetscCall(DSTruncate(eps->ds,eps->nconv,PETSC_TRUE));
  PetscCall(BVRestoreSplitRows(eps->V,is[0],is[1],&U,&V));

  PetscCall(VecDestroy(&v));
  PetscFunctionReturn(PETSC_SUCCESS);
}
