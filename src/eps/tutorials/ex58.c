/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Linear Response eigenvalue problem.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = dimension of the blocks.\n"
  "  -reduced <0/1>, to use the reduced form of the LREP.\n\n";

#include <slepceps.h>

/*
   This example is similar to ex55.c, but with a real matrix

        H = [  A   B
              -B  -A ],

   where A,B are real symmetric with the following Toeplitz structure:

        A = pentadiag{a,b,c,b,a}
        B = tridiag{b,d,b}

   where a,b,c,d are real values.
*/

int main(int argc,char **argv)
{
  Mat            H,A,B,K,M;  /* problem matrices */
  EPS            eps;        /* eigenproblem solver context */
  PetscScalar    a,b,c,d;
  PetscReal      lev;
  PetscInt       n=24,Istart,Iend,i,nconv;
  PetscBool      terse,checkorthog,nest=PETSC_FALSE,reduced=PETSC_FALSE;
  Vec            t,*x,*y;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCheck(!PetscDefined(USE_COMPLEX),PETSC_COMM_SELF,PETSC_ERR_SUP,"This example cannot be used with complex scalars");

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-reduced",&reduced,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nLinear Response Eigenvalue Problem, n=%" PetscInt_FMT "%s\n\n",n,reduced?" (reduced)":""));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               Compute the problem matrices A and B
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  a = -0.1;
  b = 1.0;
  c = 4.5;
  d = 2.0;

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(B));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>1) PetscCall(MatSetValue(A,i,i-2,a,INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A,i,i-1,b,INSERT_VALUES));
    PetscCall(MatSetValue(A,i,i,c,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,i,i+1,b,INSERT_VALUES));
    if (i<n-2) PetscCall(MatSetValue(A,i,i+2,a,INSERT_VALUES));
  }

  PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(B,i,i-1,b,INSERT_VALUES));
    PetscCall(MatSetValue(B,i,i,d,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(B,i,i+1,b,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));

  if (reduced) {

    PetscInt ma,na,Ma,Na;
    const PetscScalar scal[] = { 1.0, -1.0 };

    PetscCall(MatGetSize(A,&Ma,&Na));
    PetscCall(MatGetLocalSize(A,&ma,&na));

    /* K = A-B */
    PetscCall(MatCreate(PetscObjectComm((PetscObject)A),&K));
    PetscCall(MatSetSizes(K,ma,na,Ma,Na));
    PetscCall(MatSetType(K,MATCOMPOSITE));
    PetscCall(MatCompositeAddMat(K,A));
    PetscCall(MatCompositeAddMat(K,B));
    PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatCompositeSetScalings(K,scal));

    /* M = A+B */
    PetscCall(MatCreate(PetscObjectComm((PetscObject)A),&M));
    PetscCall(MatSetSizes(M,ma,na,Ma,Na));
    PetscCall(MatSetType(M,MATCOMPOSITE));
    PetscCall(MatCompositeAddMat(M,A));
    PetscCall(MatCompositeAddMat(M,B));
    PetscCall(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY));

    PetscCall(MatCreateLREP(K,M,reduced,&H));
    PetscCall(MatDestroy(&K));
    PetscCall(MatDestroy(&M));

  } else PetscCall(MatCreateLREP(A,B,reduced,&H));

  /* if you prefer, set the vector type so that MatCreateVecs() returns nested vectors */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-nest",&nest,NULL));
  if (nest) PetscCall(MatNestSetVecType(H,VECNEST));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,H,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_LREP));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Solve the eigensystem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }

  /* check bi-orthogonality */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-checkorthog",&checkorthog));
  PetscCall(EPSGetConverged(eps,&nconv));
  if (checkorthog && nconv>0) {
    PetscCall(MatCreateVecs(H,&t,NULL));
    PetscCall(VecDuplicateVecs(t,nconv,&x));
    PetscCall(VecDuplicateVecs(t,nconv,&y));
    for (i=0;i<nconv;i++) {
      PetscCall(EPSGetEigenvector(eps,i,x[i],NULL));
      PetscCall(EPSGetLeftEigenvector(eps,i,y[i],NULL));
    }
    PetscCall(VecCheckOrthogonality(x,nconv,y,nconv,NULL,NULL,&lev));
    if (lev<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Level of bi-orthogonality of eigenvectors < 100*eps\n\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD," Level of bi-orthogonality of eigenvectors: %g\n\n",(double)lev));
    PetscCall(VecDestroy(&t));
    PetscCall(VecDestroyVecs(nconv,&x));
    PetscCall(VecDestroyVecs(nconv,&y));
  }

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(MatDestroy(&H));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   build:
      requires: !complex

   testset:
      args: -eps_nev 4 -eps_ncv 16 -terse -checkorthog -reduced {{0 1}} -nest {{0 1}}
      nsize: {{1 2}}
      filter: sed -e "s/ (reduced)//"
      output_file: output/ex58_1.out
      test:
         suffix: 1
      test:
         suffix: 1_dense
         args: -mat_type dense
      test:
         suffix: 1_cuda
         args: -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 1_hip
         args: -mat_type aijhipsparse
         requires: hip

   test:
      args: -n 90 -eps_threshold_absolute 2.4 -eps_ncv 10 -terse -checkorthog -reduced {{0 1}}
      filter: sed -e "s/ (reduced)//"
      suffix: 2

TEST*/
