/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSGNHEP with upper quasi-triangular matrix pair.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  DS             ds;
  PetscScalar    *A,*B,*X;
  PetscReal      rnorm,aux;
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

  /* Fill A,B with upper quasi-triangular matrices */
  PetscCall(DSGetArray(ds,DS_MAT_A,&A));
  PetscCall(DSGetArray(ds,DS_MAT_B,&B));
  PetscCall(PetscArrayzero(A,ld*n));
  for (i=0;i<n;i++) A[i+i*ld]=2.0;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) A[i+(i+j)*ld]=0.001;
  }
  PetscCall(PetscArrayzero(B,ld*n));
  for (i=0;i<n;i++) B[i+i*ld]=1.0;
  B[1+0*ld]=B[0+1*ld]=PETSC_MACHINE_EPSILON;
  for (i=1;i<n;i+=3) {
    A[i+(i-1)*ld]=-A[(i-1)+i*ld];
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_A,&A));
  PetscCall(DSRestoreArray(ds,DS_MAT_B,&B));
  PetscCall(DSSetState(ds,DS_STATE_INTERMEDIATE));

  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  /* Eigenvectors */
  j = 0;
  PetscCall(DSVectors(ds,DS_MAT_X,&j,&rnorm));  /* first eigenvector */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Value of rnorm for 2nd vector = %.3f\n",(double)rnorm));
  PetscCall(DSVectors(ds,DS_MAT_X,NULL,NULL));  /* all eigenvectors */
  j = 0;
  rnorm = 0.0;
  PetscCall(DSGetArray(ds,DS_MAT_X,&X));
  for (i=0;i<n;i++) {
#if PetscDefined(USE_COMPLEX)
    aux = PetscAbsScalar(X[i+j*ld]);
#else
    aux = SlepcAbsEigenvalue(X[i+j*ld],X[i+(j+1)*ld]);
#endif
    rnorm += aux*aux;
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_X,&X));
  rnorm = PetscSqrtReal(rnorm);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Norm of 1st columns = %.3f\n",(double)rnorm));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After vectors - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  PetscCall(DSDestroy(&ds));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      filter: sed -e "s/-0\./0./"

TEST*/
