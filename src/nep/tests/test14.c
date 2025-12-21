/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests a user-defined convergence test in NEP.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n";

/*
   Solve T(lambda)x=0 with T(lambda) = -D+sqrt(lambda)*I
      where D is the Laplacian operator in 1 dimension
*/

#include <slepcnep.h>

/*
  MyConvergedRel - Convergence test relative to the norm of D (given in ctx).
*/
PetscErrorCode MyConvergedRel(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal norm = *(PetscReal*)ctx;

  PetscFunctionBegin;
  *errest = res/norm;
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  NEP            nep;             /* nonlinear eigensolver context */
  Mat            A[2];
  PetscInt       n=100,Istart,Iend,i;
  PetscBool      terse;
  PetscReal      norm;
  FN             f[2];
  PetscScalar    coeffs;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSquare root eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear eigensolver, define problem in split form
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));

  /* Create matrices */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A[0]));
  PetscCall(MatSetSizes(A[0],PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A[0]));
  PetscCall(MatGetOwnershipRange(A[0],&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A[0],i,i-1,1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A[0],i,i+1,1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A[0],i,i,-2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A[0],MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A[0],MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,1.0,&A[1]));

  /* Define functions */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[0]));
  PetscCall(FNSetType(f[0],FNRATIONAL));
  coeffs = 1.0;
  PetscCall(FNRationalSetNumerator(f[0],1,&coeffs));
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[1]));
  PetscCall(FNSetType(f[1],FNSQRT));
  PetscCall(NEPSetSplitOperator(nep,2,A,f,SUBSET_NONZERO_PATTERN));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   Set some options and solve
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPSetTarget(nep,1.1));

  /* setup convergence test relative to the norm of D */
  PetscCall(MatNorm(A[0],NORM_1,&norm));
  PetscCall(NEPSetConvergenceTestFunction(nep,MyConvergedRel,&norm,NULL));

  PetscCall(NEPSetFromOptions(nep));
  PetscCall(NEPSolve(nep));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(NEPConvergedReasonView(nep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(NEPDestroy(&nep));
  PetscCall(MatDestroy(&A[0]));
  PetscCall(MatDestroy(&A[1]));
  PetscCall(FNDestroy(&f[0]));
  PetscCall(FNDestroy(&f[1]));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -nep_type slp -nep_nev 2 -terse
      output_file: output/test14_1.out
      test:
         suffix: 1
      test:
         suffix: 1_norm
         args: -nep_conv_norm

TEST*/
