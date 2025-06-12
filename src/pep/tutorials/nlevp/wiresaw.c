/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example implements two of the problems found at
       NLEVP: A Collection of Nonlinear Eigenvalue Problems,
       The University of Manchester.
   The details of the collection can be found at:
       [1] T. Betcke et al., "NLEVP: A Collection of Nonlinear Eigenvalue
           Problems", ACM Trans. Math. Software 39(2), Article 7, 2013.

   WIRESAW1 is a gyroscopic QEP from vibration analysis of a wiresaw,
   where the parameter V represents the speed of the wire. When the
   parameter eta is nonzero, then it turns into the WIRESAW2 problem
   (with added viscous damping, e.g. eta=0.8).

       [2] S. Wei and I. Kao, "Vibration analysis of wire and frequency
           response in the modern wiresaw manufacturing process", J. Sound
           Vib. 213(5):1383-1395, 2000.
*/

static char help[] = "Vibration analysis of a wiresaw.\n\n"
  "The command line options are:\n"
  "  -n <n> ... dimension of the matrices (default 10).\n"
  "  -v <value> ... velocity of the wire (default 0.01).\n"
  "  -eta <value> ... viscous damping (default 0.0).\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            M,D,K,A[3];      /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       n=10,Istart,Iend,j,k;
  PetscReal      v=0.01,eta=0.0;
  PetscBool      terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-v",&v,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-eta",&eta,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nVibration analysis of a wiresaw, n=%" PetscInt_FMT " v=%g eta=%g\n\n",n,(double)v,(double)eta));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*D+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a diagonal matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&K));
  PetscCall(MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(K));

  PetscCall(MatGetOwnershipRange(K,&Istart,&Iend));
  for (j=Istart;j<Iend;j++) PetscCall(MatSetValue(K,j,j,(j+1)*(j+1)*PETSC_PI*PETSC_PI*(1.0-v*v),INSERT_VALUES));

  PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatScale(K,0.5));

  /* D is a tridiagonal */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&D));
  PetscCall(MatSetSizes(D,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(D));

  PetscCall(MatGetOwnershipRange(D,&Istart,&Iend));
  for (j=Istart;j<Iend;j++) {
    for (k=0;k<n;k++) {
      if ((j+k)%2) PetscCall(MatSetValue(D,j,k,8.0*(j+1)*(k+1)*v/((j+1)*(j+1)-(k+1)*(k+1)),INSERT_VALUES));
    }
  }

  PetscCall(MatAssemblyBegin(D,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(D,MAT_FINAL_ASSEMBLY));
  PetscCall(MatScale(D,0.5));

  /* M is a diagonal matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&M));
  PetscCall(MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(M));
  PetscCall(MatGetOwnershipRange(M,&Istart,&Iend));
  for (j=Istart;j<Iend;j++) PetscCall(MatSetValue(M,j,j,1.0,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY));
  PetscCall(MatScale(M,0.5));

  /* add damping */
  if (eta>0.0) {
    PetscCall(MatAXPY(K,eta,D,DIFFERENT_NONZERO_PATTERN)); /* K = K + eta*D */
    PetscCall(MatShift(D,eta)); /* D = D + eta*eye(n) */
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  A[0] = K; A[1] = D; A[2] = M;
  PetscCall(PEPSetOperators(pep,3,A));
  if (eta==0.0) PetscCall(PEPSetProblemType(pep,PEP_GYROSCOPIC));
  else PetscCall(PEPSetProblemType(pep,PEP_GENERAL));
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
  PetscCall(MatDestroy(&D));
  PetscCall(MatDestroy(&K));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -pep_nev 4 -terse
      requires: double
      output_file: output/wiresaw_1.out
      test:
         suffix: 1
         args: -pep_type {{toar qarnoldi}}
      test:
         suffix: 1_linear_h1
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization 1,0 -pep_linear_st_ksp_type bcgs -pep_linear_st_pc_type kaczmarz
      test:
         suffix: 1_linear_h2
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization 0,1
      test:
         suffix: 1_linear_other
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization .3,.7

   testset:
      args: -pep_nev 4 -eta 0.1 -terse
      requires: double
      output_file: output/wiresaw_2.out
      test:
         suffix: 2
         args: -pep_type {{toar qarnoldi}}
      test:
         suffix: 2_linear
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization {{1,0 0,1}} -pep_linear_st_ksp_type bcgs -pep_linear_st_pc_type kaczmarz
      test:
         suffix: 2_linear_other
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization .3,.7

TEST*/
