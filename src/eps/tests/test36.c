/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests a HEP problem with Hermitian matrix.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;          /* matrix */
  EPS            eps;        /* eigenproblem solver context */
  PetscInt       N,n=20,m,Istart,Iend,II,i,j;
  PetscBool      flag;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  if (!flag) m=n;
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nHermitian Eigenproblem, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,m));
#if !defined(PETSC_USE_COMPLEX)
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This example requires complex scalars!");
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,-1.0-0.1*PETSC_i,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A,II,II+n,-1.0+0.1*PETSC_i,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,-1.0-0.1*PETSC_i,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,-1.0+0.1*PETSC_i,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,4.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSSolve(eps));
  PetscCall(EPSErrorView(eps,EPS_ERROR_BACKWARD,NULL));

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   build:
      requires: complex

   testset:
      args: -m 18 -n 19 -eps_nev 4 -eps_max_it 1000
      requires: !single complex
      output_file: output/test36_1.out
      test:
         suffix: 1
         args: -eps_type {{krylovschur subspace arnoldi gd lapack}}
      test:
         suffix: 1_jd
         args: -eps_type jd
         requires: !__float128
      test:
         suffix: 1_elemental
         args: -eps_type elemental
         requires: elemental

   test:
      suffix: 2
      args: -eps_nev 4 -eps_smallest_real -eps_type {{lobpcg rqcg lapack}}
      requires: !single complex

TEST*/
