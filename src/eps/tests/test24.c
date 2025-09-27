/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Eigenproblem for the 1-D Laplacian with constraints.\n\n"
  "Based on ex1.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;
  EPS            eps;
  EPSType        type;
  Vec            *vi=NULL,*vc=NULL,t;
  PetscInt       n=30,nev=4,i,j,Istart,Iend,nini=0,ncon=0,bs;
  PetscReal      alpha,beta,restart;
  PetscBool      flg,lock;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nini",&nini,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ncon",&ncon,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%" PetscInt_FMT " nini=%" PetscInt_FMT " ncon=%" PetscInt_FMT "\n\n",n,nini,ncon));

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
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetType(eps,EPSLOBPCG));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));
  PetscCall(EPSSetConvergenceTest(eps,EPS_CONV_ABS));
  PetscCall(EPSSetDimensions(eps,nev,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(EPSLOBPCGSetBlockSize(eps,nev));
  PetscCall(EPSLOBPCGSetRestart(eps,0.7));
  PetscCall(EPSSetTolerances(eps,1e-8,1200));
  PetscCall(EPSSetFromOptions(eps));

  PetscCall(MatCreateVecs(A,&t,NULL));
  if (nini) {
    PetscCall(VecDuplicateVecs(t,nini,&vi));
    for (i=0;i<nini;i++) PetscCall(VecSetRandom(vi[i],NULL));
    PetscCall(EPSSetInitialSpace(eps,nini,vi));
  }
  if (ncon) {   /* constraints are exact eigenvectors of lowest eigenvalues */
    alpha = PETSC_PI/(n+1);
    beta  = PetscSqrtReal(2.0/(n+1));
    PetscCall(VecDuplicateVecs(t,ncon,&vc));
    for (i=0;i<ncon;i++) {
      for (j=0;j<n;j++) PetscCall(VecSetValue(vc[i],j,PetscSinReal(alpha*(j+1)*(i+1))*beta,INSERT_VALUES));
      PetscCall(VecAssemblyBegin(vc[i]));
      PetscCall(VecAssemblyEnd(vc[i]));
    }
    PetscCall(EPSSetDeflationSpace(eps,ncon,vc));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps,EPSLOBPCG,&flg));
  if (flg) {
    PetscCall(EPSLOBPCGGetLocking(eps,&lock));
    if (lock) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using soft locking\n"));
    PetscCall(EPSLOBPCGGetRestart(eps,&restart));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," LOBPCG Restart parameter=%.4g\n",(double)restart));
    PetscCall(EPSLOBPCGGetBlockSize(eps,&bs));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," LOBPCG Block size=%" PetscInt_FMT "\n",bs));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroyVecs(nini,&vi));
  PetscCall(VecDestroyVecs(ncon,&vc));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -ncon 2
      output_file: output/test24_1.out
      test:
         suffix: 1
         requires: !single
      test:
         suffix: 1_cuda
         args: -mat_type aijcusparse
         requires: cuda !single
      test:
         suffix: 1_hip
         args: -mat_type aijhipsparse
         requires: hip !single

TEST*/
