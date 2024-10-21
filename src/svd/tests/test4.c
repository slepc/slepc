/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test an SVD problem with more columns than rows.\n\n"
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
  Mat                  A,B;
  SVD                  svd;
  SVDConv              conv;
  SVDStop              stop;
  SVDWhich             which;
  SVDProblemType       ptype;
  SVDConvergedReason   reason;
  PetscInt             m=20,n,Istart,Iend,i,col[2],its;
  PetscScalar          value[] = { 1, 2 };
  PetscBool            flg,tmode;
  PetscViewerAndFormat *vf;
  const char           *ctest[] = { "absolute", "relative to the singular value", "user-defined" };
  const char           *stest[] = { "basic", "user-defined" };

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
        Compute singular values
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));
  PetscCall(SVDSetOperators(svd,A,NULL));

  /* test some interface functions */
  PetscCall(SVDGetOperators(svd,&B,NULL));
  PetscCall(MatView(B,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(SVDSetConvergenceTest(svd,SVD_CONV_ABS));
  PetscCall(SVDSetStoppingTest(svd,SVD_STOP_BASIC));
  /* test monitors */
  PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf));
  PetscCall(SVDMonitorSet(svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorFirst,vf,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy));
  /* PetscCall(SVDMonitorCancel(svd)); */
  PetscCall(SVDSetFromOptions(svd));

  /* query properties and print them */
  PetscCall(SVDGetProblemType(svd,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Problem type = %d",(int)ptype));
  PetscCall(SVDIsGeneralized(svd,&flg));
  if (flg) PetscCall(PetscPrintf(PETSC_COMM_WORLD," generalized"));
  PetscCall(SVDGetImplicitTranspose(svd,&tmode));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Transpose mode is %s\n",tmode?"implicit":"explicit"));
  PetscCall(SVDGetConvergenceTest(svd,&conv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Convergence test is %s\n",ctest[conv]));
  PetscCall(SVDGetStoppingTest(svd,&stop));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping test is %s\n",stest[stop]));
  PetscCall(SVDGetWhichSingularTriplets(svd,&which));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Which = %s\n",which?"smallest":"largest"));

  /* call the solver */
  PetscCall(SVDSolve(svd));
  PetscCall(SVDGetConvergedReason(svd,&reason));
  PetscCall(SVDGetIterationNumber(svd,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d\n",(int)reason));
  /* PetscCall(PetscPrintf(PETSC_COMM_WORLD," its = %" PetscInt_FMT "\n",its)); */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(SVDErrorView(svd,SVD_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(SVDDestroy(&svd));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -svd_monitor_cancel
      filter: grep -v "Transpose mode"
      output_file: output/test4_1.out
      test:
         suffix: 1_lanczos
         args: -svd_type lanczos
      test:
         suffix: 1_randomized
         args: -svd_type randomized
      test:
         suffix: 1_trlanczos
         args: -svd_type trlanczos -svd_ncv 12 -svd_trlanczos_restart 0.6
      test:
         suffix: 1_cross
         args: -svd_type cross
      test:
         suffix: 1_cross_exp
         args: -svd_type cross -svd_cross_explicitmatrix
      test:
         suffix: 1_cross_exp_imp
         args: -svd_type cross -svd_cross_explicitmatrix -svd_implicittranspose
         requires: !complex
      test:
         suffix: 1_cyclic
         args: -svd_type cyclic
      test:
         suffix: 1_cyclic_imp
         args: -svd_type cyclic -svd_implicittranspose
      test:
         suffix: 1_cyclic_exp
         args: -svd_type cyclic -svd_cyclic_explicitmatrix
      test:
         suffix: 1_lapack
         args: -svd_type lapack
      test:
         suffix: 1_scalapack
         args: -svd_type scalapack
         requires: scalapack

   testset:
      args: -svd_monitor_cancel -mat_type aijcusparse
      requires: cuda !single
      filter: grep -v "Transpose mode" | sed -e "s/seqaijcusparse/seqaij/"
      output_file: output/test4_1.out
      test:
         suffix: 2_cuda_lanczos
         args: -svd_type lanczos
      test:
         suffix: 2_cuda_trlanczos
         args: -svd_type trlanczos -svd_ncv 12
      test:
         suffix: 2_cuda_cross
         args: -svd_type cross

   test:
      suffix: 3
      nsize: 2
      args: -svd_type trlanczos -svd_ncv 14 -svd_monitor_cancel -ds_parallel synchronized

   testset:
      args: -svd_monitor_cancel -mat_type aijhipsparse
      requires: hip !single
      filter: grep -v "Transpose mode" | sed -e "s/seqaijhipsparse/seqaij/"
      output_file: output/test4_1.out
      test:
         suffix: 4_hip_lanczos
         args: -svd_type lanczos
      test:
         suffix: 4_hip_trlanczos
         args: -svd_type trlanczos -svd_ncv 12
      test:
         suffix: 4_hip_cross
         args: -svd_type cross

TEST*/
