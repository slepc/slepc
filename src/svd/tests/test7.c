/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "SVD via the cyclic matrix with a user-provided EPS.\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = matrix rows.\n"
  "  -n <n>, where <n> = matrix columns (defaults to m+2).\n\n";

#include <slepcsvd.h>

/*
   This example computes the singular values of a rectangular bidiagonal matrix

              |  1  2                     |
              |     1  2                  |
              |        1  2               |
          A = |          .  .             |
              |             .  .          |
              |                1  2       |
              |                   1  2    |
 */

int main(int argc,char **argv)
{
  Mat                  A;
  SVD                  svd;
  EPS                  eps;
  ST                   st;
  KSP                  ksp;
  PC                   pc;
  PetscInt             m=20,n,Istart,Iend,i,col[2];
  PetscScalar          value[] = { 1, 2 };
  PetscBool            flg,expmat;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,&flg));
  if (!flg) n=m+2;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nRectangular bidiagonal matrix, m=%" PetscInt_FMT " n=%" PetscInt_FMT "\n\n",m,n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    col[0]=i; col[1]=i+1;
    if (i<n-1) PetscCall(MatSetValues(A,1,&i,2,col,value,INSERT_VALUES));
    else if (i==n-1) PetscCall(MatSetValue(A,i,col[0],value[0],INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Create a standalone EPS with appropriate settings
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));
  PetscCall(EPSSetTarget(eps,1.0));
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STSetType(st,STSINVERT));
  PetscCall(STSetShift(st,1.01));
  PetscCall(STGetKSP(st,&ksp));
  PetscCall(KSPSetType(ksp,KSPPREONLY));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCLU));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Compute singular values
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));
  PetscCall(SVDSetOperators(svd,A,NULL));
  PetscCall(SVDSetType(svd,SVDCYCLIC));
  PetscCall(SVDCyclicSetEPS(svd,eps));
  PetscCall(SVDCyclicSetExplicitMatrix(svd,PETSC_TRUE));
  PetscCall(SVDSetWhichSingularTriplets(svd,SVD_SMALLEST));
  PetscCall(SVDSetFromOptions(svd));
  PetscCall(PetscObjectTypeCompare((PetscObject)svd,SVDCYCLIC,&flg));
  if (flg) {
    PetscCall(SVDCyclicGetExplicitMatrix(svd,&expmat));
    if (expmat) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using explicit matrix with cyclic solver\n"));
  }
  PetscCall(SVDSolve(svd));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(SVDErrorView(svd,SVD_ERROR_ABSOLUTE,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(SVDDestroy(&svd));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -log_exclude svd

   testset:
      args: -log_exclude svd
      output_file: output/test7_1.out
      test:
         suffix: 2_cuda
         args: -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 2_hip
         args: -mat_type aijhipsparse
         requires: hip
      test:
         suffix: 2_kokkos
         args: -mat_type aijkokkos
         requires: hip kokkos_kernels

TEST*/
