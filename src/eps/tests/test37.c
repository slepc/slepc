/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests solving an eigenproblem defined with MatNest.\n\n"
  "Based on ex9.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = block dimension of the 2x2 block matrix.\n"
  "  -L <L>, where <L> = bifurcation parameter.\n"
  "  -alpha <alpha>, -beta <beta>, -delta1 <delta1>,  -delta2 <delta2>,\n"
  "       where <alpha> <beta> <delta1> <delta2> = model parameters.\n\n";

#include <slepceps.h>

/*
   This example computes the eigenvalues with largest real part of the
   following matrix

        A = [ tau1*T+(beta-1)*I     alpha^2*I
                  -beta*I        tau2*T-alpha^2*I ],

   where

        T = tridiag{1,-2,1}
        h = 1/(n+1)
        tau1 = delta1/(h*L)^2
        tau2 = delta2/(h*L)^2
 */

int main(int argc,char **argv)
{
  EPS            eps;
  Mat            A,T1,T2,D1,D2,mats[4];
  PetscScalar    alpha,beta,tau1,tau2,delta1,delta2,L,h;
  PetscInt       N=30,i,Istart,Iend;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model, n=%" PetscInt_FMT "\n\n",N));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  alpha  = 2.0;
  beta   = 5.45;
  delta1 = 0.008;
  delta2 = 0.004;
  L      = 0.51302;

  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-L",&L,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-alpha",&alpha,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-beta",&beta,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-delta1",&delta1,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-delta2",&delta2,NULL));

  h    = 1.0 / (PetscReal)(N+1);
  tau1 = delta1 / ((h*L)*(h*L));
  tau2 = delta2 / ((h*L)*(h*L));

  /* Create matrices T1, T2 */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&T1));
  PetscCall(MatSetSizes(T1,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(T1));
  PetscCall(MatGetOwnershipRange(T1,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(T1,i,i-1,1.0,INSERT_VALUES));
    if (i<N-1) PetscCall(MatSetValue(T1,i,i+1,1.0,INSERT_VALUES));
    PetscCall(MatSetValue(T1,i,i,-2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(T1,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(T1,MAT_FINAL_ASSEMBLY));

  PetscCall(MatDuplicate(T1,MAT_COPY_VALUES,&T2));
  PetscCall(MatScale(T1,tau1));
  PetscCall(MatShift(T1,beta-1.0));
  PetscCall(MatScale(T2,tau2));
  PetscCall(MatShift(T2,-alpha*alpha));

  /* Create matrices D1, D2 */
  PetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N,alpha*alpha,&D1));
  PetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N,-beta,&D2));

  /* Create the nest matrix */
  mats[0] = T1;
  mats[1] = D1;
  mats[2] = D2;
  mats[3] = T2;
  PetscCall(MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,mats,&A));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
  PetscCall(EPSSetTrueResidual(eps,PETSC_FALSE));
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSSolve(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&T1));
  PetscCall(MatDestroy(&T2));
  PetscCall(MatDestroy(&D1));
  PetscCall(MatDestroy(&D2));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -eps_nev 4
      requires: !single
      filter: sed -e "s/0.00019-2.13938i, 0.00019+2.13938i/0.00019+2.13938i, 0.00019-2.13938i/" | sed -e "s/-0.67192-2.52712i, -0.67192+2.52712i/-0.67192+2.52712i, -0.67192-2.52712i/"

TEST*/
