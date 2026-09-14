/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV operations with non-standard inner product.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec            t,v;
  Mat            B,M;
  BV             X;
  PetscInt       i,j,n=10,k=5,Istart,Iend;
  PetscScalar    alpha;
  PetscReal      nrm;
  PetscViewer    view;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV with non-standard inner product (n=%" PetscInt_FMT ", k=%" PetscInt_FMT ").\n",n,k));

  /* Create inner product matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(B));
  PetscCall(PetscObjectSetName((PetscObject)B,"B"));

  PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(B,i,i-1,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(B,i,i+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(B,i,i,2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatCreateVecs(B,&t,NULL));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,k));
  PetscCall(BVSetFromOptions(X));
  PetscCall(BVSetMatrix(X,B,PETSC_FALSE));

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
  if (verbose) {
    PetscCall(MatView(B,view));
    PetscCall(BVView(X,view));
  }

  /* Test BVNormColumn */
  PetscCall(BVNormColumn(X,0,NORM_2,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"B-Norm of X[0] = %g\n",(double)nrm));

  /* Test BVOrthogonalizeColumn */
  for (j=0;j<k;j++) {
    PetscCall(BVOrthogonalizeColumn(X,j,NULL,&nrm,NULL));
    alpha = 1.0/nrm;
    PetscCall(BVScaleColumn(X,j,alpha));
  }
  if (verbose) PetscCall(BVView(X,view));

  /* Check orthogonality */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M));
  PetscCall(BVDot(X,X,M));
  PetscCall(MatShift(M,-1.0));
  PetscCall(MatNorm(M,NORM_1,&nrm));
  if (nrm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)nrm));

  /* Test BVNormVecBegin/End */
  PetscCall(BVGetColumn(X,0,&v));
  PetscCall(BVNormVecBegin(X,v,NORM_1,&nrm));
  PetscCall(BVNormVecEnd(X,v,NORM_1,&nrm));
  PetscCall(BVRestoreColumn(X,0,&v));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"B-Norm of X[0] = %g\n",(double)nrm));

  PetscCall(BVDestroy(&X));
  PetscCall(MatDestroy(&M));
  PetscCall(MatDestroy(&B));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test3_1.out
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 1_svec_vecs
         args: -bv_type svec -bv_matmult vecs
      test:
         suffix: 1_cuda
         args: -bv_type {{svec mat}} -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 1_hip
         args: -bv_type {{svec mat}} -mat_type aijhipsparse
         requires: hip
      test:
         suffix: 1_kokkos
         args: -bv_type mat -mat_type aijkokkos
         requires: hip kokkos_kernels
      test:
         suffix: 2
         nsize: 2
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 3
         nsize: 2
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_type mgs

TEST*/
