/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Block-diagonal orthogonal eigenproblem.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n"
  "  -seed <s>, where <s> = seed for random number generation.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  EPS            eps;
  Mat            A;
  PetscRandom    rand;
  PetscScalar    val,c,s;
  PetscInt       n=30,i,seed=0x12345678;
  PetscMPIInt    rank;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nOrthogonal eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&rand));
  PetscCall(PetscRandomSetFromOptions(rand));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-seed",&seed,NULL));
  PetscCall(PetscRandomSetSeed(rand,seed));
  PetscCall(PetscRandomSeed(rand));
  PetscCall(PetscRandomSetInterval(rand,0,2*PETSC_PI));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  if (!rank) {
    for (i=0;i<n/2;i++) {
      PetscCall(PetscRandomGetValue(rand,&val));
      c = PetscCosReal(PetscRealPart(val));
      s = PetscSinReal(PetscRealPart(val));
      PetscCall(MatSetValue(A,2*i,2*i,c,INSERT_VALUES));
      PetscCall(MatSetValue(A,2*i+1,2*i+1,c,INSERT_VALUES));
      PetscCall(MatSetValue(A,2*i,2*i+1,s,INSERT_VALUES));
      PetscCall(MatSetValue(A,2*i+1,2*i,-s,INSERT_VALUES));
    }
    if (n%2) PetscCall(MatSetValue(A,n-1,n-1,-1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSSolve(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscRandomDestroy(&rand));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      requires: complex double
      args: -eps_type ciss -eps_all -rg_type ring -rg_ring_center 0 -rg_ring_radius 1 -rg_ring_width 0.05 -rg_ring_startangle .93 -rg_ring_endangle .07
      filter: sed -e "s/[+-]\([0-9]\.[0-9]*i\)/+-\\1/g"
      test:
         suffix: 1_ring
         args: -eps_ciss_extraction {{ritz hankel}}

TEST*/
