/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSGNHEP.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  DS             ds;
  SlepcSC        sc;
  PetscScalar    *A,*B,*X,*wr,*wi;
  PetscReal      re,im,rnorm,aux;
  PetscInt       i,j,n=10,ld;
  PetscViewer    viewer;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type GNHEP - dimension %" PetscInt_FMT ".\n",n));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));

  /* Create DS object */
  PetscCall(DSCreate(PETSC_COMM_WORLD,&ds));
  PetscCall(DSSetType(ds,DSGNHEP));
  PetscCall(DSSetFromOptions(ds));
  ld = n+2;  /* test leading dimension larger than n */
  PetscCall(DSAllocate(ds,ld));
  PetscCall(DSSetDimensions(ds,n,0,0));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(DSView(ds,viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill A with Grcar matrix */
  PetscCall(DSGetArray(ds,DS_MAT_A,&A));
  PetscCall(PetscArrayzero(A,ld*n));
  for (i=1;i<n;i++) A[i+(i-1)*ld]=-1.0;
  for (j=0;j<4;j++) {
    for (i=0;i<n-j;i++) A[i+(i+j)*ld]=1.0;
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_A,&A));
  /* Fill B with an upper triangular matrix */
  PetscCall(DSGetArray(ds,DS_MAT_B,&B));
  PetscCall(PetscArrayzero(B,ld*n));
  B[0+0*ld]=-1.0;
  B[0+1*ld]=2.0;
  for (i=1;i<n;i++) B[i+i*ld]=1.0;
  PetscCall(DSRestoreArray(ds,DS_MAT_B,&B));

  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  /* Solve */
  PetscCall(PetscMalloc2(n,&wr,n,&wi));
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

  /* Print eigenvalues */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed eigenvalues =\n"));
  for (i=0;i<n;i++) {
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(wr[i]);
    im = PetscImaginaryPart(wr[i]);
#else
    re = wr[i];
    im = wi[i];
#endif
    if (PetscAbs(im)<1e-10) PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f\n",(double)re));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f%+.5fi\n",(double)re,(double)im));
  }

  /* Eigenvectors */
  j = 1;
  PetscCall(DSVectors(ds,DS_MAT_X,&j,&rnorm));  /* second eigenvector */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Value of rnorm for 2nd vector = %.3f\n",(double)rnorm));
  PetscCall(DSVectors(ds,DS_MAT_X,NULL,NULL));  /* all eigenvectors */
  j = 0;
  rnorm = 0.0;
  PetscCall(DSGetArray(ds,DS_MAT_X,&X));
  for (i=0;i<n;i++) {
#if PetscDefined(USE_COMPLEX)
    aux = PetscAbsScalar(X[i+j*ld]);
#else
    if (PetscAbs(wi[j])==0.0) aux = PetscAbsScalar(X[i+j*ld]);
    else aux = SlepcAbsEigenvalue(X[i+j*ld],X[i+(j+1)*ld]);
#endif
    rnorm += aux*aux;
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_X,&X));
  rnorm = PetscSqrtReal(rnorm);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Norm of 1st vector = %.3f\n",(double)rnorm));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After vectors - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  PetscCall(PetscFree2(wr,wi));
  PetscCall(DSDestroy(&ds));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      filter: sed -e "s/[+-]\([0-9]\.[0-9]*i\)/+-\\1/" | sed -e "s/+-0\.0*i//" | sed -e "s/1.29552/1.29551/"

TEST*/
