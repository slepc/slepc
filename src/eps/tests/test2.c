/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests multiple calls to EPSSolve with the same matrix.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  ST             st;
  PetscReal      tol=PetscMax(1000*PETSC_MACHINE_EPSILON,1e-9);
  PetscInt       n=30,i,Istart,Iend;
  PetscBool      flg;
  EPSLanczosReorthogType reorth;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A,i,i-1,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,i,i+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,i,i,2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Create the eigensolver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));
  PetscCall(EPSSetFromOptions(eps));

  /* illustrate how to extract parameters from specific solver types */
  PetscCall(PetscObjectTypeCompare((PetscObject)eps,EPSLANCZOS,&flg));
  if (flg) {
    PetscCall(EPSLanczosGetReorthog(eps,&reorth));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Reorthogonalization type used in Lanczos: %s\n",EPSLanczosReorthogTypes[reorth]));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Solve for largest eigenvalues
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
  PetscCall(EPSSolve(eps));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," - - - Largest eigenvalues - - -\n"));
  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Solve for smallest eigenvalues
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));
  PetscCall(EPSSolve(eps));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," - - - Smallest eigenvalues - - -\n"));
  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Solve for interior eigenvalues (target=2.1)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));
  PetscCall(EPSSetTarget(eps,2.1));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps,EPSLANCZOS,&flg));
  if (flg) {
    PetscCall(EPSGetST(eps,&st));
    PetscCall(STSetType(st,STSINVERT));
  } else {
    PetscCall(PetscObjectTypeCompare((PetscObject)eps,EPSKRYLOVSCHUR,&flg));
    if (!flg) PetscCall(PetscObjectTypeCompare((PetscObject)eps,EPSARNOLDI,&flg));
    if (flg) PetscCall(EPSSetExtraction(eps,EPS_HARMONIC));
  }
  PetscCall(EPSSolve(eps));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," - - - Interior eigenvalues - - -\n"));
  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_nev 4
      requires: !single
      output_file: output/test2_1.out
      test:
         suffix: 1
         args: -eps_type {{arnoldi gd jd lapack}}
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion
         timeoutfactor: 2
      test:
         suffix: 1_krylovschur
         args: -eps_type krylovschur -eps_krylovschur_locking {{0 1}}
      test:
         suffix: 1_scalapack
         requires: scalapack !__float128
         args: -eps_type scalapack
      test:
         suffix: 1_elpa
         requires: elpa
         args: -eps_type elpa
      test:
         suffix: 1_elemental
         requires: elemental
         args: -eps_type elemental

   testset:
      args: -eps_type lanczos -eps_nev 4
      requires: !single
      filter: grep -v "Lanczos"
      output_file: output/test2_1.out
      test:
         suffix: 2
         args: -eps_lanczos_reorthog {{local full periodic partial}}
      test:
         suffix: 2_selective
         args: -eps_lanczos_reorthog selective
         requires: !defined(PETSCTEST_VALGRIND)

   testset:
      args: -n 32 -eps_nev 4
      requires: !single
      output_file: output/test2_3.out
      test:
         nsize: 2
         suffix: 3
         args: -eps_type {{krylovschur lapack}}
      test:
         nsize: 2
         suffix: 3_gd
         args: -eps_type gd -eps_gd_krylov_start
         timeoutfactor: 2
      test:
         suffix: 3_jd
         args: -eps_type jd -eps_jd_krylov_start -eps_ncv 18

   testset:
      args: -eps_nev 4 -eps_type {{krylovschur arnoldi gd jd}}
      requires: !single
      output_file: output/test2_1.out
      test:
         suffix: 4_cuda
         requires: cuda
         args: -mat_type aijcusparse
      test:
         suffix: 4_hip
         requires: hip
         args: -mat_type aijhipsparse

TEST*/
