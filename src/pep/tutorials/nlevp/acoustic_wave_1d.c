/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example implements one of the problems found at
       NLEVP: A Collection of Nonlinear Eigenvalue Problems,
       The University of Manchester.
   The details of the collection can be found at:
       [1] T. Betcke et al., "NLEVP: A Collection of Nonlinear Eigenvalue
           Problems", ACM Trans. Math. Software 39(2), Article 7, 2013.

   The acoustic_wave_1d problem is a QEP from an acoustics application.
   Here we solve it with the eigenvalue scaled by the imaginary unit, to be
   able to use real arithmetic, so the computed eigenvalues should be scaled
   back.
*/

static char help[] = "Quadratic eigenproblem from an acoustics application (1-D).\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = dimension of the matrices.\n"
  "  -z <z>, where <z> = impedance (default 1.0).\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3];      /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       n=10,Istart,Iend,i;
  PetscScalar    z=1.0;
  char           str[50];
  PetscBool      terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-z",&z,NULL));
  PetscCall(SlepcSNPrintfScalar(str,sizeof(str),z,PETSC_FALSE));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nAcoustic wave 1-D, n=%" PetscInt_FMT " z=%s\n\n",n,str));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a tridiagonal */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&K));
  PetscCall(MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(K));

  PetscCall(MatGetOwnershipRange(K,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(K,i,i-1,-1.0*n,INSERT_VALUES));
    if (i<n-1) {
      PetscCall(MatSetValue(K,i,i,2.0*n,INSERT_VALUES));
      PetscCall(MatSetValue(K,i,i+1,-1.0*n,INSERT_VALUES));
    } else PetscCall(MatSetValue(K,i,i,1.0*n,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

  /* C is the zero matrix but one element*/
  PetscCall(MatCreate(PETSC_COMM_WORLD,&C));
  PetscCall(MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(C));

  PetscCall(MatGetOwnershipRange(C,&Istart,&Iend));
  if (n-1>=Istart && n-1<Iend) PetscCall(MatSetValue(C,n-1,n-1,-2*PETSC_PI/z,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY));

  /* M is a diagonal matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&M));
  PetscCall(MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(M));

  PetscCall(MatGetOwnershipRange(M,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i<n-1) PetscCall(MatSetValue(M,i,i,4*PETSC_PI*PETSC_PI/n,INSERT_VALUES));
    else PetscCall(MatSetValue(M,i,i,2*PETSC_PI*PETSC_PI/n,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  A[0] = K; A[1] = C; A[2] = M;
  PetscCall(PEPSetOperators(pep,3,A));
  PetscCall(PEPSetFromOptions(pep));
  PetscCall(PEPSolve(pep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(PEPConvergedReasonView(pep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(PEPDestroy(&pep));
  PetscCall(MatDestroy(&M));
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&K));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -pep_nev 4 -pep_tol 1e-7 -n 24 -terse
      output_file: output/acoustic_wave_1d_1.out
      requires: !single
      filter: sed -e "s/47426/47425/g"
      test:
         suffix: 1
         args: -st_type sinvert -st_transform -pep_type {{toar qarnoldi linear}}
      test:
         suffix: 1_stoar
         args: -st_type sinvert -st_transform -pep_type stoar -pep_hermitian -pep_stoar_locking 0 -pep_stoar_nev 11 -pep_ncv 10
      test:
         suffix: 2
         args: -st_type sinvert -st_transform -pep_type toar -pep_extract {{none norm residual}}
      test:
         suffix: 3
         args: -st_type sinvert -pep_type linear -pep_extract {{none norm residual}}
      test:
         suffix: 4
         args: -pep_type jd

TEST*/
