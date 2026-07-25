/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSPEP.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  DS             ds;
  SlepcSC        sc;
  Mat            X;
  Vec            x0;
  PetscScalar    *K,*C,*M,*wr,*wi,z=1.0;
  PetscReal      re,im,nrm,*pbc;
  PetscInt       i,j,n=10,d=2,ld;
  PetscViewer    viewer;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type PEP - n=%" PetscInt_FMT ".\n",n));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));

  /* Create DS object */
  PetscCall(DSCreate(PETSC_COMM_WORLD,&ds));
  PetscCall(DSSetType(ds,DSPEP));
  PetscCall(DSSetFromOptions(ds));
  PetscCall(DSPEPSetDegree(ds,d));

  /* Set dimensions */
  ld = n+2;  /* test leading dimension larger than n */
  PetscCall(DSAllocate(ds,ld));
  PetscCall(DSSetDimensions(ds,n,0,0));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(DSView(ds,viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill matrices */
  PetscCall(DSGetArray(ds,DS_MAT_E0,&K));
  for (i=0;i<n-1;i++) K[i+i*ld] = 2.0*n;
  K[n-1+(n-1)*ld] = 1.0*n;
  for (i=1;i<n;i++) {
    K[i+(i-1)*ld] = -1.0*n;
    K[(i-1)+i*ld] = -1.0*n;
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_E0,&K));
  PetscCall(DSGetArray(ds,DS_MAT_E1,&C));
  C[n-1+(n-1)*ld] = 2.0*PETSC_PI/z;
  PetscCall(DSRestoreArray(ds,DS_MAT_E1,&C));
  PetscCall(DSGetArray(ds,DS_MAT_E2,&M));
  for (i=0;i<n-1;i++) M[i+i*ld] = -4.0*PETSC_PI*PETSC_PI/n;
  M[i+i*ld] = -2.0*PETSC_PI*PETSC_PI/n;
  PetscCall(DSRestoreArray(ds,DS_MAT_E2,&M));

  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  /* Solve */
  PetscCall(PetscMalloc2(d*n,&wr,d*n,&wi));
  PetscCall(DSGetSlepcSC(ds,&sc));
  sc->comparison    = SlepcCompareLargestReal;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSSolve(ds,wr,wi));
  PetscCall(DSSort(ds,wr,wi,NULL,NULL,NULL));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After solve - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  /* Print polynomial coefficients */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Polynomial coefficients (alpha,beta,gamma) =\n"));
  PetscCall(DSPEPGetCoefficients(ds,&pbc));
  for (j=0;j<3;j++) {
    for (i=0;i<d+1;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f",(double)pbc[j+3*i]));
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  }
  PetscCall(PetscFree(pbc));

  /* Print eigenvalues */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed eigenvalues =\n"));
  for (i=0;i<d*n;i++) {
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
  PetscCall(DSVectors(ds,DS_MAT_X,NULL,NULL));  /* all eigenvectors */
  PetscCall(DSGetMat(ds,DS_MAT_X,&X));
  PetscCall(MatCreateVecs(X,NULL,&x0));
  PetscCall(MatGetColumnVector(X,x0,0));
  PetscCall(VecNorm(x0,NORM_2,&nrm));
  PetscCall(DSRestoreMat(ds,DS_MAT_X,&X));
  PetscCall(VecDestroy(&x0));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Norm of 1st vector = %.3f\n",(double)nrm));
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
      requires: !single

TEST*/
