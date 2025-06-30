/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests multiple calls to EPSSolve with different matrix.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A1,A2;       /* problem matrices */
  EPS            eps;         /* eigenproblem solver context */
  PetscReal      tol=PETSC_SMALL,v;
  Vec            d;
  PetscInt       n=30,i,Istart,Iend;
  PetscRandom    myrand;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nTridiagonal with random diagonal, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Create matrix tridiag([-1 0 -1])
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A1));
  PetscCall(MatSetSizes(A1,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A1));

  PetscCall(MatGetOwnershipRange(A1,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A1,i,i-1,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A1,i,i+1,-1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A1,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A1,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create two matrices by filling the diagonal with rand values
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatDuplicate(A1,MAT_COPY_VALUES,&A2));
  PetscCall(MatCreateVecs(A1,NULL,&d));
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&myrand));
  PetscCall(PetscRandomSetFromOptions(myrand));
  PetscCall(PetscRandomSetInterval(myrand,0.0,1.0));
  for (i=Istart;i<Iend;i++) {
    PetscCall(PetscRandomGetValueReal(myrand,&v));
    PetscCall(VecSetValue(d,i,v,INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(d));
  PetscCall(VecAssemblyEnd(d));
  PetscCall(MatDiagonalSet(A1,d,INSERT_VALUES));
  for (i=Istart;i<Iend;i++) {
    PetscCall(PetscRandomGetValueReal(myrand,&v));
    PetscCall(VecSetValue(d,i,v,INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(d));
  PetscCall(VecAssemblyEnd(d));
  PetscCall(MatDiagonalSet(A2,d,INSERT_VALUES));
  PetscCall(VecDestroy(&d));
  PetscCall(PetscRandomDestroy(&myrand));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Create the eigensolver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));
  PetscCall(EPSSetOperators(eps,A1,NULL));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve first eigenproblem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSSolve(eps));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," - - - First matrix - - -\n"));
  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve second eigenproblem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSSetOperators(eps,A2,NULL));
  PetscCall(EPSSolve(eps));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," - - - Second matrix - - -\n"));
  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A1));
  PetscCall(MatDestroy(&A2));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_nev 4
      requires: !single
      output_file: output/test3_1.out
      test:
         suffix: 1
         args: -eps_type {{krylovschur subspace arnoldi lapack}}
      test:
         suffix: 1_lanczos
         args: -eps_type lanczos -eps_lanczos_reorthog local
      test:
         suffix: 1_power
         args: -eps_type power -eps_max_it 20000
      test:
         suffix: 1_jd
         args: -eps_type jd -eps_jd_initial_size 7
      test:
         suffix: 1_gd
         args: -eps_type gd -eps_gd_initial_size 7
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion
      test:
         suffix: 1_arpack
         args: -eps_type arpack
         requires: arpack
      test:
         suffix: 1_primme
         args: -eps_type primme -eps_conv_abs -eps_primme_blocksize 4
         requires: primme
      test:
         suffix: 1_scalapack
         args: -eps_type scalapack
         requires: scalapack
      test:
         suffix: 1_elpa
         args: -eps_type elpa
         requires: elpa
      test:
         suffix: 1_elemental
         args: -eps_type elemental
         requires: elemental

   testset:
      args: -eps_nev 4 -eps_smallest_real -eps_max_it 500
      output_file: output/test3_2.out
      test:
         suffix: 2_rqcg
         args: -eps_type rqcg -eps_rqcg_reset 5 -eps_ncv 32
      test:
         suffix: 2_lobpcg
         args: -eps_type lobpcg -eps_lobpcg_blocksize 5 -st_pc_type none
      test:
         suffix: 2_lanczos
         args: -eps_type lanczos -eps_lanczos_reorthog local
         requires: !single
      test:
         suffix: 2_lanczos_delayed
         args: -eps_type lanczos -eps_lanczos_reorthog delayed -eps_tol 1e-8
         requires: !single
      test:
         suffix: 2_blopex
         args: -eps_type blopex -eps_conv_abs -st_shift -2
         requires: blopex

TEST*/
