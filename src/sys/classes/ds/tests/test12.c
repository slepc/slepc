/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSNEP.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  DS             ds;
  FN             f1,f2,f3,funs[3],qfun;
  SlepcSC        sc;
  PetscScalar    *Id,*A,*B,*wr,*wi,*X,*W,coeffs[2],auxr,alpha;
  PetscReal      tol,tau=0.001,radius=10,h,a=20,xi,re,im,nrm,aux;
  PetscInt       i,j,ii,jj,k,n=10,ld,nev,nfun;
  PetscViewer    viewer;
  PetscBool      verbose;
  RG             rg;
  DSMatType      mat[3]={DS_MAT_E0,DS_MAT_E1,DS_MAT_E2};
#if !PetscDefined(USE_COMPLEX)
  PetscScalar    auxi;
#endif

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-tau",&tau,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type NEP - dimension %" PetscInt_FMT ", tau=%g.\n",n,(double)tau));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-radius",&radius,NULL));

  /* Create DS object */
  PetscCall(DSCreate(PETSC_COMM_WORLD,&ds));
  PetscCall(DSSetType(ds,DSNEP));
  tol  = 1000*n*PETSC_MACHINE_EPSILON;
  PetscCall(DSNEPSetRefine(ds,tol,PETSC_DECIDE));
  PetscCall(DSSetFromOptions(ds));

  /* Set functions (prior to DSAllocate) */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f1));
  PetscCall(FNSetType(f1,FNRATIONAL));
  coeffs[0] = -1.0; coeffs[1] = 0.0;
  PetscCall(FNRationalSetNumerator(f1,2,coeffs));

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f2));
  PetscCall(FNSetType(f2,FNRATIONAL));
  coeffs[0] = 1.0;
  PetscCall(FNRationalSetNumerator(f2,1,coeffs));

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f3));
  PetscCall(FNSetType(f3,FNEXP));
  PetscCall(FNSetScale(f3,-tau,1.0));

  funs[0] = f1;
  funs[1] = f2;
  funs[2] = f3;
  PetscCall(DSNEPSetFN(ds,3,funs));

  /* Set dimensions */
  ld = n+2;  /* test leading dimension larger than n */
  PetscCall(DSAllocate(ds,ld));
  PetscCall(DSSetDimensions(ds,n,0,0));

  /* Set region (used only in method=1) */
  PetscCall(RGCreate(PETSC_COMM_WORLD,&rg));
  PetscCall(RGSetType(rg,RGELLIPSE));
  PetscCall(RGEllipseSetParameters(rg,0.0,radius,1.0));
  PetscCall(DSNEPSetRG(ds,rg));
  PetscCall(RGDestroy(&rg));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(DSView(ds,viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Show info about functions */
  PetscCall(DSNEPGetNumFN(ds,&nfun));
  for (i=0;i<nfun;i++) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Function %" PetscInt_FMT ":\n",i));
    PetscCall(DSNEPGetFN(ds,i,&qfun));
    PetscCall(FNView(qfun,NULL));
  }

  /* Fill matrices */
  PetscCall(DSGetArray(ds,DS_MAT_E0,&Id));
  for (i=0;i<n;i++) Id[i+i*ld]=1.0;
  PetscCall(DSRestoreArray(ds,DS_MAT_E0,&Id));
  h = PETSC_PI/(PetscReal)(n+1);
  PetscCall(DSGetArray(ds,DS_MAT_E1,&A));
  for (i=0;i<n;i++) A[i+i*ld]=-2.0/(h*h)+a;
  for (i=1;i<n;i++) {
    A[i+(i-1)*ld]=1.0/(h*h);
    A[(i-1)+i*ld]=1.0/(h*h);
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_E1,&A));
  PetscCall(DSGetArray(ds,DS_MAT_E2,&B));
  for (i=0;i<n;i++) {
    xi = (i+1)*h;
    B[i+i*ld] = -4.1+xi*(1.0-PetscExpReal(xi-PETSC_PI));
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_E2,&B));

  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  /* Solve */
  PetscCall(PetscCalloc2(n,&wr,n,&wi));
  PetscCall(DSGetSlepcSC(ds,&sc));
  sc->comparison    = SlepcCompareLargestMagnitude;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSSolve(ds,wr,wi));
  PetscCall(DSSort(ds,wr,wi,NULL,NULL,NULL));

  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After solve - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }
  PetscCall(DSGetDimensions(ds,NULL,NULL,NULL,&nev));

  /* Print computed eigenvalues */
  PetscCall(PetscMalloc1(ld*ld,&W));
  PetscCall(DSVectors(ds,DS_MAT_X,NULL,NULL));
  PetscCall(DSGetArray(ds,DS_MAT_X,&X));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed eigenvalues =\n"));
  for (i=0;i<nev;i++) {
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(wr[i]);
    im = PetscImaginaryPart(wr[i]);
#else
    re = wr[i];
    im = wi[i];
#endif
    /* Residual */
    PetscCall(PetscArrayzero(W,ld*ld));
    for (k=0;k<nfun;k++) {
      PetscCall(FNEvaluateFunction(funs[k],wr[i],&alpha));
      PetscCall(DSGetArray(ds,mat[k],&A));
      for (jj=0;jj<n;jj++) for (ii=0;ii<n;ii++) W[jj*ld+ii] += alpha*A[jj*ld+ii];
      PetscCall(DSRestoreArray(ds,mat[k],&A));
    }
    nrm = 0.0;
    for (k=0;k<n;k++) {
      auxr = 0.0;
#if !PetscDefined(USE_COMPLEX)
      auxi = 0.0;
#endif
      for (j=0;j<n;j++) {
        auxr += W[k+j*ld]*X[i*ld+j];
#if !PetscDefined(USE_COMPLEX)
        if (PetscAbs(wi[j])!=0.0) auxi += W[k+j*ld]*X[(i+1)*ld+j];
#endif
      }
      aux = SlepcAbsEigenvalue(auxr,auxi);
      nrm += aux*aux;
    }
    nrm = PetscSqrtReal(nrm);
    if (nrm/SlepcAbsEigenvalue(wr[i],wi[i])>tol) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the residual norm of the %" PetscInt_FMT "-th computed eigenpair %g\n",i,(double)nrm));
    if (PetscAbs(im)<1e-10) PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f\n",(double)re));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f%+.5fi\n",(double)re,(double)im));
  }
  PetscCall(PetscFree(W));
  PetscCall(DSRestoreArray(ds,DS_MAT_X,&X));
  PetscCall(DSRestoreArray(ds,DS_MAT_W,&W));
  PetscCall(PetscFree2(wr,wi));
  PetscCall(FNDestroy(&f1));
  PetscCall(FNDestroy(&f2));
  PetscCall(FNDestroy(&f3));
  PetscCall(DSDestroy(&ds));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      test:
         filter: grep -v "solving the problem"
         suffix: 1
      test:
         suffix: 2
         args: -ds_method 1 -radius 10 -ds_nep_refine_its 1
         filter: grep -v "solving the problem" | sed -e "s/[+-]0\.0*i//g" | sed -e "s/37411/37410/" | sed -e "s/tolerance [0-9]\.[0-9]*e[+-]\([0-9]*\)/tolerance removed/" | sed -e "s/tolerance [0-9]\.\([0-9]*\)/tolerance removed/"
         requires: complex

TEST*/
