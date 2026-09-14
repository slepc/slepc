/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV operations with indefinite inner product.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec            t,v,w,omega;
  Mat            B,M;
  BV             X,Y;
  PetscInt       i,j,n=10,k=5,l,Istart,Iend;
  PetscScalar    alpha;
  PetscReal      nrm;
  PetscViewer    view;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV with indefinite inner product (n=%" PetscInt_FMT ", k=%" PetscInt_FMT ").\n",n,k));

  /* Create inner product matrix (standard involutionary permutation) */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(B));
  PetscCall(PetscObjectSetName((PetscObject)B,"B"));

  PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(B,i,n-i-1,1.0,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatCreateVecs(B,&t,NULL));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,k));
  PetscCall(BVSetFromOptions(X));
  PetscCall(BVSetMatrix(X,B,PETSC_TRUE));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view));
  if (verbose) PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill X entries */
  l = -3;
  for (j=0;j<k;j++) {
    PetscCall(BVGetColumn(X,j,&v));
    PetscCall(VecSet(v,-1.0));
    for (i=0;i<n/2;i++) {
      if (i+j<n) {
        l = (l + 3*i+j-2) % n;
        PetscCall(VecSetValue(v,i+j,(PetscScalar)l,INSERT_VALUES));
      }
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

  /* Create a copy on Y */
  PetscCall(BVDuplicate(X,&Y));
  PetscCall(PetscObjectSetName((PetscObject)Y,"Y"));
  PetscCall(BVCopy(X,Y));

  /* Check orthogonality */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M));
  PetscCall(BVDot(Y,Y,M));
  PetscCall(VecCreateSeq(PETSC_COMM_SELF,k,&omega));
  PetscCall(BVGetSignature(Y,omega));
  PetscCall(VecScale(omega,-1.0));
  PetscCall(MatDiagonalSet(M,omega,ADD_VALUES));
  PetscCall(MatNorm(M,NORM_1,&nrm));
  if (nrm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)nrm));

  /* Test BVSetSignature */
  PetscCall(VecScale(omega,-1.0));
  PetscCall(BVSetSignature(Y,omega));
  PetscCall(VecDestroy(&omega));

  /* Test BVApplyMatrix */
  PetscCall(VecDuplicate(t,&w));
  PetscCall(BVGetColumn(X,0,&v));
  PetscCall(BVApplyMatrix(X,v,w));
  PetscCall(BVApplyMatrix(X,w,t));
  PetscCall(VecAXPY(t,-1.0,v));
  PetscCall(BVRestoreColumn(X,0,&v));
  PetscCall(VecNorm(t,NORM_2,&nrm));
  PetscCheck(PetscAbsReal(nrm)<10*PETSC_MACHINE_EPSILON,PETSC_COMM_WORLD,PETSC_ERR_FP,"Wrong value, nrm = %g",(double)nrm);

  PetscCall(BVApplyMatrixBV(X,Y));
  PetscCall(BVGetColumn(Y,0,&v));
  PetscCall(VecAXPY(w,-1.0,v));
  PetscCall(BVRestoreColumn(Y,0,&v));
  PetscCall(VecNorm(w,NORM_2,&nrm));
  PetscCheck(PetscAbsReal(nrm)<10*PETSC_MACHINE_EPSILON,PETSC_COMM_WORLD,PETSC_ERR_FP,"Wrong value, nrm = %g",(double)nrm);

  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(MatDestroy(&M));
  PetscCall(MatDestroy(&B));
  PetscCall(VecDestroy(&w));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test5_1.out
      args: -bv_orthog_refine always
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}shared output}
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
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_type mgs
      test:
         suffix: 2_cuda
         args: -bv_type {{svec mat}} -mat_type aijcusparse -bv_orthog_type mgs
         requires: cuda
      test:
         suffix: 2_hip
         args: -bv_type {{svec mat}} -mat_type aijhipsparse -bv_orthog_type mgs
         requires: hip
      test:
         suffix: 2_kokkos
         args: -bv_type mat -mat_type aijkokkos -bv_orthog_type mgs
         requires: hip kokkos_kernels

TEST*/
