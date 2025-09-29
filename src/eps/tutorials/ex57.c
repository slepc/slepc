/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Another eigenvalue problem with Hamiltonian structure.\n\n"
  "Position and velocity control of a string of high-speed vehicles.\n"
  "The command line options are:\n"
  "  -m <m>, where <n> = number of vehicles.\n\n";

#include <slepceps.h>

/*
   This example computes eigenvalues of a matrix

        H = [ A    N
              K  -A^* ],

   where A, N and K are n-by-n matrices, with n=2*m-1 and

        N = diag( 1, 0, 1, 0,..., 0, 1)
        K = diag( 0,10, 0,10,...,10, 0)
        A = tridiag(b,d,c)
        d = [-1, 0,-1, 0,...,-1, 0,-1]
        b = [ 1, 0, 1, 0,..., 1, 0]
        c = [ 0,-1, 0,-1,..., 0,-1]

   References:

   [1] W.R. Ferng, Wen-Wei Lin, Chern-Shuh Wang, The shift-inverted J-Lanczos algorithm
       for the numerical solutions of large sparse algebraic Riccati equations,
       Computers & Mathematics with Applications, Volume 33, Issue 10, 1997,
       https://doi.org/10.1016/S0898-1221(97)00074-6.

   */

static PetscErrorCode eigenCompare (PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx)
{
  PetscReal abs1, abs2, tol=1e-12, r_ar, r_ai, r_br, r_bi;

  PetscFunctionBegin;
  #if defined(PETSC_USE_COMPLEX)
    r_ar = PetscRealPart(ar);
    r_ai = PetscImaginaryPart(ar);
    r_br = PetscRealPart(br);
    r_bi = PetscImaginaryPart(br);
  #else
    r_ar = ar;
    r_ai = ai;
    r_br = br;
    r_bi = bi;
  #endif
  if (PetscAbs(PetscAbs(r_ar)-PetscAbs(r_br)) < tol && PetscAbs(PetscAbs(r_ai)-PetscAbs(r_bi)) < tol) {
    /* Negative real part first*/
    if (r_ar<0 && r_br>0)
      *res = -1;
    else if (r_ar>0 && r_br<0)
      *res = 1;
    /* Positive imaginary part first*/
    else if (r_ai>0 && r_bi<0)
      *res = -1;
    else if (r_ai<0 && r_bi>0)
      *res = 1;
  }
  else {
    abs1 = SlepcAbs(r_ar,r_ai);
    abs2 = SlepcAbs(r_br,r_bi);
    if (abs1 > abs2)
      *res = -1;
    else if (abs1 < abs2)
      *res = 1;
    else
      *res = 0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Mat            H,A,N,K;    /* problem matrices */
  EPS            eps;        /* eigenproblem solver context */
  PetscInt       m=12,n,Istart,Iend,i;
  PetscBool      terse, sort_hamilt;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  n = 2*m-1;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nHamiltonian eigenproblem, n=%" PetscInt_FMT
    " (m=%" PetscInt_FMT " vehicles)\n\n",n,m));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               Compute the problem matrices A, N and K
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&N));
  PetscCall(MatSetSizes(N,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(N));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&K));
  PetscCall(MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(K));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i%2==0)
      PetscCall(MatSetValue(A,i,i,-1.0,INSERT_VALUES));
    else {
      if (i>0) PetscCall(MatSetValue(A,i,i-1,1.0,INSERT_VALUES));
      if (i<n-1) PetscCall(MatSetValue(A,i,i+1,-1.0,INSERT_VALUES));
    }
  }

  PetscCall(MatGetOwnershipRange(N,&Istart,&Iend));
  for (i=Istart+Istart%2;i<Iend;i+=2) {
    PetscCall(MatSetValue(N,i,i,1.0,INSERT_VALUES));
  }

  PetscCall(MatGetOwnershipRange(K,&Istart,&Iend));
  for (i=Istart+(Istart+1)%2;i<Iend;i+=2) {
    PetscCall(MatSetValue(K,i,i,10.0,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(N,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(N,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreateHamiltonian(A,N,K,&H));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,H,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HAMILT));

  PetscCall(PetscOptionsHasName(NULL,NULL,"-sort_hamilt",&sort_hamilt));
  if (sort_hamilt) {
    /* Adjust ordering of non-hermitian solver so that it is the same as in as in EPS_HAMILT solver */
    PetscCall(EPSSetEigenvalueComparison(eps,eigenCompare,NULL));
    PetscCall(EPSSetWhichEigenpairs(eps,EPS_WHICH_USER));
  }
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
  PetscCall(MatDestroy(&N));
  PetscCall(MatDestroy(&K));
  PetscCall(MatDestroy(&H));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_nev 8 -eps_ncv 28 -terse
      nsize: {{1 2}}
      output_file: output/ex57_1.out
      test:
         requires: double !complex
         suffix: 1
      test:
         requires: double !complex
         args: -eps_non_hermitian -sort_hamilt
         suffix: 1_nhep
      test:
         requires: double complex
         TODO: no support for complex scalars yet
         suffix: 1_complex
      test:
         requires: double complex
         args: -eps_non_hermitian -sort_hamilt
         suffix: 1_complex_nhep

   testset:
      args: -eps_nev 4 -eps_smallest_magnitude -eps_ncv 32 -terse
      output_file: output/ex57_2.out
      test:
         requires: double !complex
         suffix: 2

TEST*/
