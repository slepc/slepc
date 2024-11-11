/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 2 dimensions.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  EPSStop        stop;
  PetscReal      thres;
  PetscInt       N,n=10,m,Istart,Iend,II,nev,i,j;
  PetscBool      flag,terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  if (!flag) m=n;
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n2-D Laplacian Eigenproblem, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,m));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,-1.0,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A,II,II+n,-1.0,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,-1.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,4.0,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetStoppingTest(eps,&stop));
  if (stop!=EPS_STOP_THRESHOLD) {
    PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  } else {
    PetscCall(EPSGetThreshold(eps,&thres,NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using threshold: %.4g\n",(double)thres));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -n 72 -eps_nev 4 -eps_ncv 20 -terse
      output_file: output/ex2_1.out
      requires: !single
      test:
         suffix: 1
      test:
         suffix: 2
         requires: defined(PETSC_HAVE_DYNAMIC_LIBRARIES) defined(PETSC_USE_SHARED_LIBRARIES)
         args: -library_preload

   testset:
      args: -n 30 -eps_type ciss -eps_ciss_realmats -terse
      requires: !single
      output_file: output/ex2_ciss.out
      filter: grep -v method
      test:
         suffix: ciss_1
         nsize: 1
         args: -rg_type interval -rg_interval_endpoints 1.1,1.25,-.1,.1
         requires: complex
      test:
         suffix: ciss_1_hpddm
         nsize: 1
         args: -rg_type interval -rg_interval_endpoints 1.1,1.25 -st_ksp_type hpddm
         requires: hpddm
      test:
         suffix: ciss_2
         nsize: 2
         args: -rg_type ellipse -rg_ellipse_center 1.175 -rg_ellipse_radius 0.075 -eps_ciss_partitions 2
      test:
         suffix: ciss_2_block
         args: -rg_type ellipse -rg_ellipse_center 1.175 -rg_ellipse_radius 0.075 -eps_ciss_blocksize 3 -eps_ciss_moments 2
         requires: complex !__float128
      test:
         suffix: ciss_2_hpddm
         nsize: 2
         args: -rg_type ellipse -rg_ellipse_center 1.175 -rg_ellipse_radius 0.075 -eps_ciss_partitions 2 -eps_ciss_ksp_type hpddm
         requires: hpddm
      test:
         suffix: feast
         args: -eps_type feast -eps_interval 1.1,1.25 -eps_ncv 64 -options_left 0
         requires: feast

   testset:
      args: -n 30 -m 30 -eps_interval 3.9,4.15 -terse
      output_file: output/ex2_3.out
      filter: grep -v Solution
      requires: !single
      test:
         suffix: 3
         args: -st_type sinvert -st_pc_type cholesky
      test:
         suffix: 3_evsl
         args: -eps_type evsl -eps_evsl_slices 6
         requires: evsl

   testset:
      args: -n 45 -m 46 -eps_interval 4.54,4.57 -eps_ncv 24 -terse
      output_file: output/ex2_4.out
      filter: grep -v Solution
      requires: !single
      timeoutfactor: 2
      test:
         suffix: 4
         args: -st_type sinvert -st_pc_type cholesky
      test:
         suffix: 4_filter
         args: -eps_type {{krylovschur subspace}} -st_type filter -st_filter_degree 200
         requires: !__float128
      test:
         suffix: 4_filter_cuda
         args: -eps_type {{krylovschur subspace}} -st_type filter -st_filter_degree 200 -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 4_filter_hip
         args: -eps_type {{krylovschur subspace}} -st_type filter -st_filter_degree 200 -mat_type aijhipsparse
         requires: hip
      test:
         suffix: 4_evsl
         args: -eps_type evsl
         requires: evsl

   test:
      args: -n 25 -m 24 -eps_threshold_absolute .25 -eps_smallest_magnitude -eps_ncv 10 -terse
      suffix: 5
      requires: !single

   testset:
      args: -n 25 -m 24 -st_type sinvert -terse
      requires: double
      test:
         suffix: 6
         args: -eps_threshold_absolute .15 -eps_target 0.01
      test:
         suffix: 6_rel_large
         args: -eps_threshold_relative .98 -eps_target 8
      test:
         suffix: 6_rel_small
         args: -eps_threshold_relative 3

TEST*/
