/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Eigenvalue problem with Hamiltonian structure.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = dimension of the blocks.\n\n";

#include <slepceps.h>

/*
   This example computes eigenvalues of a matrix

        H = [ A    B
              C  -A^* ],

   where B and C are Hermitian or real symmetric. In particular, A, B and C have
   the following Toeplitz structure:

        A = pentadiag{a,b,c,d,e}
        B = tridiag{a,c,conj(a)}
        C = tridiag{b,e,conj(b)}

   where a,b,d are complex scalars, and c is real.
*/

int main(int argc,char **argv)
{
  Mat            H,A,B,C;    /* problem matrices */
  EPS            eps;        /* eigenproblem solver context */
  PetscScalar    a,b,c,d,e;
  PetscInt       n=24,Istart,Iend,i;
  PetscBool      terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nHamiltonian eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               Compute the problem matrices R and C
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#if defined(PETSC_USE_COMPLEX)
  a = PetscCMPLX(-0.1,0.2);
  b = PetscCMPLX(1.0,0.5);
  d = PetscCMPLX(2.0,0.2);
#else
  a = -0.1;
  b = 1.0;
  d = 2.0;
#endif
  c = 4.5;
  e = -2.5;

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(B));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&C));
  PetscCall(MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(C));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>1) PetscCall(MatSetValue(A,i,i-2,a,INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A,i,i-1,b,INSERT_VALUES));
    PetscCall(MatSetValue(A,i,i,c,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,i,i+1,d,INSERT_VALUES));
    if (i<n-2) PetscCall(MatSetValue(A,i,i+2,e,INSERT_VALUES));
  }

  PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(B,i,i-1,a,INSERT_VALUES));
    PetscCall(MatSetValue(B,i,i,c,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(B,i,i+1,PetscConj(a),INSERT_VALUES));
  }

  PetscCall(MatGetOwnershipRange(C,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(C,i,i-1,b,INSERT_VALUES));
    PetscCall(MatSetValue(C,i,i,e,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(C,i,i+1,PetscConj(b),INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreateHamiltonian(A,B,C,&H));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,H,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HAMILT));
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
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&H));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_nev 8 -eps_ncv 28 -terse
      nsize: {{1 2}}
      requires: double !complex
      output_file: output/ex56_1.out
      test:
         suffix: 1
      test:
         args: -eps_non_hermitian
         suffix: 1_nhep

   testset:
      args: -eps_nev 4 -eps_ncv 16 -terse
      nsize: {{1 2}}
      requires: double complex
      output_file: output/ex56_1_complex.out
      test:
         TODO: no support for complex scalars yet
         suffix: 1_complex
      test:
         args: -eps_non_hermitian
         suffix: 1_complex_nhep

TEST*/
