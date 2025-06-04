/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSGSVD.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  DS             ds;
  SlepcSC        sc;
  Mat            X;
  Vec            x0;
  PetscReal      sigma,rnorm,cond;
  PetscScalar    *A,*B,*w;
  PetscInt       i,j,k,n=15,m=10,p=10,m1,p1,ld;
  PetscViewer    viewer;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-p",&p,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type GSVD - dimension (%" PetscInt_FMT "+%" PetscInt_FMT ")x%" PetscInt_FMT ".\n",n,p,m));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));

  /* Create DS object */
  PetscCall(DSCreate(PETSC_COMM_WORLD,&ds));
  PetscCall(DSSetType(ds,DSGSVD));
  PetscCall(DSSetFromOptions(ds));
  ld   = PetscMax(PetscMax(p,m),n)+2;  /* test leading dimension larger than n */
  PetscCall(DSAllocate(ds,ld));
  PetscCall(DSSetDimensions(ds,n,0,0));
  PetscCall(DSGSVDSetDimensions(ds,m,p));
  PetscCall(DSGSVDGetDimensions(ds,&m1,&p1));
  PetscCheck(m1==m && p1==p,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Inconsistent dimension values");

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(DSView(ds,viewer));
  PetscCall(PetscViewerPopFormat(viewer));

  k = PetscMin(n,m);
  /* Fill A with a rectangular Toeplitz matrix */
  PetscCall(DSGetArray(ds,DS_MAT_A,&A));
  for (i=0;i<k;i++) A[i+i*ld]=1.0;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) {
      if ((i+j)<m) A[i+(i+j)*ld]=(PetscScalar)(j+1);
    }
  }
  for (j=1;j<n/2;j++) {
    for (i=0;i<n-j;i++) {
      if ((i+j)<n && i<m) A[(i+j)+i*ld]=-1.0;
    }
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_A,&A));

  k = PetscMin(p,m);
  /* Fill B with a shifted bidiagonal matrix */
  PetscCall(DSGetArray(ds,DS_MAT_B,&B));
  for (i=m-k;i<m;i++) {
    B[i-m+k+i*ld]=2.0-1.0/(PetscScalar)(i+1);
    if (i) B[i-1-m+k+i*ld]=1.0;
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_B,&B));

  PetscCall(DSSetState(ds,DS_STATE_RAW));
  if (verbose) {
    PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n"));
  }
  PetscCall(DSView(ds,viewer));

  /* Condition number */
  PetscCall(DSCond(ds,&cond));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Condition number = %.3f\n",(double)cond));

  /* Solve */
  PetscCall(PetscMalloc1(m,&w));
  PetscCall(DSGetSlepcSC(ds,&sc));
  sc->comparison    = SlepcCompareLargestReal;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSSolve(ds,w,NULL));
  PetscCall(DSSort(ds,w,NULL,NULL,NULL,NULL));
  PetscCall(DSSynchronize(ds,w,NULL));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After solve - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }
  /* Print singular values */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed singular values =\n"));
  PetscCall(DSGetDimensions(ds,NULL,NULL,NULL,&k));
  for (i=0;i<k;i++) {
    sigma = PetscRealPart(w[i]);
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %g\n",(double)sigma));
  }

  /* Singular vectors */
  PetscCall(DSVectors(ds,DS_MAT_X,NULL,NULL));  /* all singular vectors */
  PetscCall(DSGetMat(ds,DS_MAT_X,&X));
  PetscCall(MatCreateVecs(X,NULL,&x0));
  PetscCall(MatGetColumnVector(X,x0,0));
  PetscCall(VecNorm(x0,NORM_2,&rnorm));
  PetscCall(DSRestoreMat(ds,DS_MAT_X,&X));
  PetscCall(VecDestroy(&x0));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Norm of 1st X vector = %.3f\n",(double)rnorm));

  PetscCall(DSGetMat(ds,DS_MAT_U,&X));
  PetscCall(MatCreateVecs(X,NULL,&x0));
  PetscCall(MatGetColumnVector(X,x0,0));
  PetscCall(VecNorm(x0,NORM_2,&rnorm));
  PetscCall(DSRestoreMat(ds,DS_MAT_U,&X));
  PetscCall(VecDestroy(&x0));
  if (PetscAbs(rnorm-1.0)>10*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the 1st U vector has norm %g\n",(double)rnorm));

  PetscCall(DSGetMat(ds,DS_MAT_V,&X));
  PetscCall(MatCreateVecs(X,NULL,&x0));
  PetscCall(MatGetColumnVector(X,x0,0));
  PetscCall(VecNorm(x0,NORM_2,&rnorm));
  PetscCall(DSRestoreMat(ds,DS_MAT_V,&X));
  PetscCall(VecDestroy(&x0));
  if (PetscAbs(rnorm-1.0)>10*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the 1st V vector has norm %g\n",(double)rnorm));

  PetscCall(PetscFree(w));
  PetscCall(DSDestroy(&ds));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test21_1.out
      requires: !single
      nsize: {{1 2 3}}
      filter: grep -v "parallel operation mode" | grep -v " MPI process"
      test:
         suffix: 1
         args: -ds_parallel redundant
      test:
         suffix: 2
         args: -ds_parallel synchronized

TEST*/
