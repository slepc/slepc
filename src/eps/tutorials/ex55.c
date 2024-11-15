/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Eigenvalue problem with Bethe-Salpeter structure.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = dimension of the blocks.\n\n";

#include <slepceps.h>

/*
   This example computes eigenvalues of a matrix

        H = [  R    C
              -C^H -R^T ],

   where R is Hermitian and C is complex symmetric. In particular, R and C have the
   following Toeplitz structure:

        R = pentadiag{a,b,c,conj(b),conj(a)}
        C = tridiag{b,d,b}

   where a,b,d are complex scalars, and c is real.
*/

int main(int argc,char **argv)
{
  Mat            H,R,C;      /* problem matrices */
  EPS            eps;        /* eigenproblem solver context */
  PetscScalar    a,b,c,d;
  PetscReal      lev;
  PetscInt       n=24,Istart,Iend,i,nconv;
  PetscBool      terse,checkorthog;
  Vec            t,*x,*y;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nBethe-Salpeter eigenproblem, n=%" PetscInt_FMT "\n\n",n));

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

  PetscCall(MatCreate(PETSC_COMM_WORLD,&R));
  PetscCall(MatSetSizes(R,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(R));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&C));
  PetscCall(MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(C));

  PetscCall(MatGetOwnershipRange(R,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>1) PetscCall(MatSetValue(R,i,i-2,a,INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(R,i,i-1,b,INSERT_VALUES));
    PetscCall(MatSetValue(R,i,i,c,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(R,i,i+1,PetscConj(b),INSERT_VALUES));
    if (i<n-2) PetscCall(MatSetValue(R,i,i+2,PetscConj(a),INSERT_VALUES));
  }

  PetscCall(MatGetOwnershipRange(C,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(C,i,i-1,b,INSERT_VALUES));
    PetscCall(MatSetValue(C,i,i,d,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(C,i,i+1,b,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreateBSE(R,C,&H));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,H,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_BSE));
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
  PetscCall(MatDestroy(&R));
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&H));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_nev 4 -eps_ncv 16 -eps_krylovschur_bse_type {{shao gruning projectedbse}} -terse -checkorthog
      filter: sed -e "s/17496/17495/g" | sed -e "s/38566/38567/g"
      nsize: {{1 2}}
      test:
         suffix: 1
         requires: complex
      test:
         suffix: 1_real
         requires: !complex

   testset:
      args: -eps_nev 4 -eps_ncv 16 -eps_krylovschur_bse_type {{shao gruning projectedbse}} -st_type sinvert -terse
      filter: sed -e "s/17496/17495/g" | sed -e "s/38566/38567/g"
      test:
         suffix: 1_sinvert
         requires: complex
      test:
         nsize: 4
         args: -mat_type scalapack
         suffix: 1_sinvert_scalapack
         requires: complex scalapack
         output_file: output/ex55_1_sinvert.out
      test:
         suffix: 1_real_sinvert
         requires: !complex
      test:
         nsize: 4
         args: -mat_type scalapack
         suffix: 1_real_sinvert_scalapack
         requires: !complex scalapack
         output_file: output/ex55_1_real_sinvert.out

   testset:
      args: -n 90 -eps_threshold_absolute 2.5 -eps_ncv {{10 24}} -eps_krylovschur_bse_type {{shao gruning projectedbse}} -eps_tol 1e-14 -terse -checkorthog
      test:
         suffix: 2
         requires: double complex
      test:
         suffix: 2_real
         requires: double !complex

TEST*/
