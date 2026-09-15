/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV operations using internal buffer instead of array arguments.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec            t,v,z;
  BV             X;
  PetscInt       i,j,n=10,k=5,l=3;
  PetscReal      nrm;
  PetscViewer    view;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV with %" PetscInt_FMT " columns of dimension %" PetscInt_FMT ".\n",k,n));

  /* Create template vector */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&t));
  PetscCall(VecSetSizes(t,PETSC_DECIDE,n));
  PetscCall(VecSetFromOptions(t));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,k));
  PetscCall(BVSetFromOptions(X));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view));
  if (verbose) PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill X entries */
  for (j=0;j<k;j++) {
    PetscCall(BVGetColumn(X,j,&v));
    PetscCall(VecSet(v,0.0));
    for (i=0;i<4;i++) {
      if (i+j<n) PetscCall(VecSetValue(v,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(v));
    PetscCall(VecAssemblyEnd(v));
    PetscCall(BVRestoreColumn(X,j,&v));
  }
  if (verbose) PetscCall(BVView(X,view));

  /* Test BVDotColumn */
  PetscCall(BVDotColumn(X,2,NULL));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVDotColumn - - - - - - -\n"));
    PetscCall(BVGetBufferVec(X,&z));
    PetscCall(VecView(z,view));
  }
  /* Test BVMultColumn */
  PetscCall(BVMultColumn(X,-1.0,1.0,2,NULL));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVMultColumn - - - - - - - - -\n"));
    PetscCall(BVView(X,view));
  }

  PetscCall(BVNorm(X,NORM_FROBENIUS,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Frobenius Norm or X = %g\n",(double)nrm));

  PetscCall(BVDestroy(&X));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test13_1.out
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 1_cuda
         args: -bv_type {{svec mat}} -vec_type cuda
         requires: cuda
      test:
         suffix: 1_hip
         args: -bv_type {{svec mat}} -vec_type hip
         requires: hip
      test:
         suffix: 1_kokkos
         args: -bv_type mat -vec_type kokkos
         requires: hip kokkos_kernels

TEST*/
