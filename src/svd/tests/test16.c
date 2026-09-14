/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests multiple calls to SVDSolve with equal matrix size (GSVD).\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = number of rows of A.\n"
  "  -n <n>, where <n> = number of columns of A.\n"
  "  -p <p>, where <p> = number of rows of B.\n\n";

#include <slepcsvd.h>

/*
   This example solves two GSVD problems for the bidiagonal matrices

              |  1  2                     |       |  1                        |
              |     1  2                  |       |  2  1                     |
              |        1  2               |       |     2  1                  |
         A1 = |          .  .             |  A2 = |       .  .                |
              |             .  .          |       |          .  .             |
              |                1  2       |       |             2  1          |
              |                   1  2    |       |                2  1       |

   with B = tril(ones(p,n))
 */

int main(int argc,char **argv)
{
  Mat            A1,A2,B;
  SVD            svd;
  PetscInt       m=15,n=20,p=21,Istart,Iend,i,j,d,col[2];
  PetscScalar    valsa[] = { 1, 2 }, valsb[] = { 2, 1 };

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-p",&p,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized singular value decomposition, (%" PetscInt_FMT "+%" PetscInt_FMT ")x%" PetscInt_FMT "\n\n",m,p,n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                     Generate the matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A1));
  PetscCall(MatSetSizes(A1,PETSC_DECIDE,PETSC_DECIDE,m,n));
  PetscCall(MatSetFromOptions(A1));
  PetscCall(MatGetOwnershipRange(A1,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    col[0]=i; col[1]=i+1;
    if (i<n-1) PetscCall(MatSetValues(A1,1,&i,2,col,valsa,INSERT_VALUES));
    else if (i==n-1) PetscCall(MatSetValue(A1,i,col[0],valsa[0],INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A1,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A1,MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A2));
  PetscCall(MatSetSizes(A2,PETSC_DECIDE,PETSC_DECIDE,m,n));
  PetscCall(MatSetFromOptions(A2));
  PetscCall(MatGetOwnershipRange(A2,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    col[0]=i-1; col[1]=i;
    if (i==0) PetscCall(MatSetValue(A2,i,col[1],valsb[1],INSERT_VALUES));
    else if (i<n) PetscCall(MatSetValues(A2,1,&i,2,col,valsb,INSERT_VALUES));
    else if (i==n) PetscCall(MatSetValue(A2,i,col[0],valsb[0],INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A2,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A2,MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,p,n));
  PetscCall(MatSetFromOptions(B));
  PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  d = PetscMax(0,n-p);
  for (i=Istart;i<Iend;i++) {
    for (j=PetscMax(0,i-5);j<=PetscMin(i,n-1);j++) PetscCall(MatSetValue(B,i,j+d,1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Create the singular value solver, set options and solve
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));
  PetscCall(SVDSetOperators(svd,A1,B));
  PetscCall(SVDSetFromOptions(svd));
  PetscCall(SVDSolve(svd));
  PetscCall(SVDErrorView(svd,SVD_ERROR_NORM,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Solve second problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDSetOperators(svd,A2,B));
  PetscCall(SVDSolve(svd));
  PetscCall(SVDErrorView(svd,SVD_ERROR_NORM,NULL));

  /* Free work space */
  PetscCall(SVDDestroy(&svd));
  PetscCall(MatDestroy(&A1));
  PetscCall(MatDestroy(&A2));
  PetscCall(MatDestroy(&B));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -svd_nsv 3
      requires: !single
      output_file: output/test16_1.out
      test:
         suffix: 1_lapack
         args: -svd_type lapack
      test:
         suffix: 1_cross
         args: -svd_type cross -svd_cross_explicitmatrix {{0 1}}
      test:
         suffix: 1_cyclic
         args: -svd_type cyclic -svd_cyclic_explicitmatrix {{0 1}}
      test:
         suffix: 1_trlanczos
         args: -svd_type trlanczos -svd_trlanczos_gbidiag {{single lower}} -svd_trlanczos_ksp_rtol 1e-10
         requires: double
      test:
         suffix: 1_trlanczos_par
         nsize: 2
         args: -svd_type trlanczos -ds_parallel {{redundant synchronized}}

   testset:
      args: -svd_nsv 3 -mat_type aijcusparse
      requires: cuda !single
      output_file: output/test16_1.out
      test:
         suffix: 2_cross
         args: -svd_type cross -svd_cross_explicitmatrix {{0 1}}
      test:
         suffix: 2_cyclic
         args: -svd_type cyclic -svd_cyclic_explicitmatrix {{0 1}}
      test:
         suffix: 2_trlanczos
         args: -svd_type trlanczos -svd_trlanczos_gbidiag {{single lower}} -svd_trlanczos_ksp_rtol 1e-10
         requires: double

   testset:
      args: -svd_nsv 3 -mat_type aijhipsparse
      requires: hip !single
      output_file: output/test16_1.out
      test:
         suffix: 2_hip_cross
         args: -svd_type cross -svd_cross_explicitmatrix {{0 1}}
      test:
         suffix: 2_hip_cyclic
         args: -svd_type cyclic -svd_cyclic_explicitmatrix {{0 1}}
      test:
         suffix: 2_hip_trlanczos
         args: -svd_type trlanczos -svd_trlanczos_gbidiag {{single lower}} -svd_trlanczos_ksp_rtol 1e-10
         requires: double

   testset:
      args: -svd_nsv 3 -mat_type aijkokkos
      requires: hip kokkos_kernels !single
      output_file: output/test16_1.out
      test:
         suffix: 2_kokkos_cross
         args: -svd_type cross -svd_cross_explicitmatrix {{0 1}}
      test:
         suffix: 2_kokkos_cyclic
         args: -svd_type cyclic -svd_cyclic_explicitmatrix {{0 1}}
      test:
         suffix: 2_kokkos_trlanczos
         args: -svd_type trlanczos -svd_trlanczos_gbidiag {{single lower}} -svd_trlanczos_ksp_rtol 1e-10
         requires: double

TEST*/
