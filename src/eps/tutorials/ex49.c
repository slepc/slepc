/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "User-defined split preconditioner when solving a generalized eigenproblem.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A,B,A0,B0,mats[2]; /* problem matrices and sparser approximations */
  EPS            eps;               /* eigenproblem solver context */
  ST             st;
  PetscInt       N,n=24,m,Istart,Iend,II,i,j;
  PetscBool      flag,terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  if (!flag) m=n;
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGHEP with split preconditioner, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,m));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Compute the problem matrices A and B
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(B));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,-0.2,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A,II,II+n,-0.2,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,-3.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,-3.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,7.0,INSERT_VALUES));
    PetscCall(MatSetValue(B,II,II,2.0,INSERT_VALUES));
  }
  if (Istart==0) {
    PetscCall(MatSetValue(B,0,0,6.0,INSERT_VALUES));
    PetscCall(MatSetValue(B,0,1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(B,1,0,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(B,1,1,1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Compute sparser approximations A0 and B0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A0));
  PetscCall(MatSetSizes(A0,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A0));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&B0));
  PetscCall(MatSetSizes(B0,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(B0));

  PetscCall(MatGetOwnershipRange(A0,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (j>0) PetscCall(MatSetValue(A0,II,II-1,-3.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A0,II,II+1,-3.0,INSERT_VALUES));
    PetscCall(MatSetValue(A0,II,II,7.0,INSERT_VALUES));
    PetscCall(MatSetValue(B0,II,II,2.0,INSERT_VALUES));
  }
  if (Istart==0) {
    PetscCall(MatSetValue(B0,0,0,6.0,INSERT_VALUES));
    PetscCall(MatSetValue(B0,1,1,1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A0,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A0,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(B0,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B0,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,B));
  PetscCall(EPSSetProblemType(eps,EPS_GHEP));
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STSetType(st,STSINVERT));
  mats[0] = A0; mats[1] = B0;
  PetscCall(STSetSplitPreconditioner(st,2,mats,SUBSET_NONZERO_PATTERN));
  PetscCall(EPSSetTarget(eps,0.0));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));
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
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(MatDestroy(&A0));
  PetscCall(MatDestroy(&B0));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_nev 4 -terse
      output_file: output/ex49_1.out
      requires: !single
      test:
         suffix: 1
      test:
         suffix: 1_jd
         args: -eps_type jd -st_type precond
      test:
         suffix: 1_lobpcg
         args: -eps_type lobpcg -st_type precond -eps_smallest_real -st_shift 0.2

   testset:
      args: -eps_type ciss -eps_all -rg_type ellipse -rg_ellipse_center 0 -rg_ellipse_radius 0.34 -rg_ellipse_vscale .2 -st_ksp_type gmres -terse
      output_file: output/ex49_2.out
      test:
         suffix: 2
      test:
         suffix: 2_split
         args: -eps_ciss_strategy split
         requires: !single
      test:
         TODO: bad residuals, only two eigenvalues
         suffix: 2_multishift
         args: -eps_ciss_strategy multishift
      test:
         suffix: 2_par
         nsize: 2
         args: -eps_ciss_partitions 2
         requires: !single

TEST*/
