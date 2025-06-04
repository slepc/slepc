/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSSVD.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  DS             ds;
  SlepcSC        sc;
  PetscReal      sigma,rnorm,aux;
  PetscScalar    *A,*U,*w,d;
  PetscInt       i,j,k,n=15,m=10,m1,ld;
  PetscViewer    viewer;
  PetscBool      verbose,extrarow;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  k = PetscMin(n,m);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type SVD - dimension %" PetscInt_FMT "x%" PetscInt_FMT ".\n",n,m));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-extrarow",&extrarow));

  /* Create DS object */
  PetscCall(DSCreate(PETSC_COMM_WORLD,&ds));
  PetscCall(DSSetType(ds,DSSVD));
  PetscCall(DSSetFromOptions(ds));
  ld = PetscMax(n,m)+2;  /* test leading dimension larger than n */
  PetscCall(DSAllocate(ds,ld));
  PetscCall(DSSetDimensions(ds,n,0,0));
  PetscCall(DSSVDSetDimensions(ds,m));
  PetscCall(DSSVDGetDimensions(ds,&m1));
  PetscCheck(m1==m,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Inconsistent dimension value");
  PetscCall(DSSetExtraRow(ds,extrarow));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(DSView(ds,viewer));
  PetscCall(PetscViewerPopFormat(viewer));

  /* Fill with a rectangular Toeplitz matrix */
  PetscCall(DSGetArray(ds,DS_MAT_A,&A));
  for (i=0;i<k;i++) A[i+i*ld]=1.0;
  for (j=1;j<3;j++) {
    for (i=0;i<m-j;i++) {
      if ((i+j)<m) A[i+(i+j)*ld]=(PetscScalar)(j+1);
    }
  }
  for (j=1;j<n/2;j++) {
    for (i=0;i<n-j;i++) {
      if ((i+j)<n && i<m) A[(i+j)+i*ld]=-1.0;
    }
  }
  if (extrarow) { A[n-2+m*ld]=1.0; A[n-1+m*ld]=1.0; }  /* really an extra column */
  PetscCall(DSRestoreArray(ds,DS_MAT_A,&A));
  PetscCall(DSSetState(ds,DS_STATE_RAW));
  if (verbose) {
    PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n"));
  }
  PetscCall(DSView(ds,viewer));

  /* Solve */
  PetscCall(PetscMalloc1(k,&w));
  PetscCall(DSGetSlepcSC(ds,&sc));
  sc->comparison    = SlepcCompareLargestReal;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSSolve(ds,w,NULL));
  PetscCall(DSSort(ds,w,NULL,NULL,NULL,NULL));
  if (extrarow) PetscCall(DSUpdateExtraRow(ds));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After solve - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  /* Print singular values */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed singular values =\n"));
  for (i=0;i<k;i++) {
    sigma = PetscRealPart(w[i]);
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f\n",(double)sigma));
  }

  if (extrarow) {
    /* Check that extra column is correct */
    PetscCall(DSGetArray(ds,DS_MAT_A,&A));
    PetscCall(DSGetArray(ds,DS_MAT_U,&U));
    d = 0.0;
    for (i=0;i<n;i++) d += A[i+m*ld]-U[n-2+i*ld]-U[n-1+i*ld];
    if (PetscAbsScalar(d)>10*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: there is a mismatch in the extra row of %g\n",(double)PetscAbsScalar(d)));
    PetscCall(DSRestoreArray(ds,DS_MAT_A,&A));
    PetscCall(DSRestoreArray(ds,DS_MAT_U,&U));
  }

  /* Singular vectors */
  PetscCall(DSVectors(ds,DS_MAT_U,NULL,NULL));  /* all singular vectors */
  j = 0;
  rnorm = 0.0;
  PetscCall(DSGetArray(ds,DS_MAT_U,&U));
  for (i=0;i<n;i++) {
    aux = PetscAbsScalar(U[i+j*ld]);
    rnorm += aux*aux;
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_U,&U));
  rnorm = PetscSqrtReal(rnorm);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Norm of 1st U vector = %.3f\n",(double)rnorm));

  PetscCall(PetscFree(w));
  PetscCall(DSDestroy(&ds));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      args: -ds_method {{0 1}}
      suffix: 1
      filter: grep -v "solving the problem"
      requires: !single

   test:
      suffix: 2
      args: -extrarow -ds_method {{0 1}}
      filter: grep -v "solving the problem"
      requires: !single

TEST*/
