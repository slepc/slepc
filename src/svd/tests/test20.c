/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test setting and getting the signature in HSVD.\n\n"
  "Based on ex15.c. The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n"
  "  -mu <mu>, where <mu> = subdiagonal value.\n\n";

#include <slepcsvd.h>

int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  Vec            omega,v;         /* signature */
  SVD            svd;             /* singular value problem solver context */
  PetscReal      mu=PETSC_SQRT_MACHINE_EPSILON;
  PetscInt       n=50,i,j,Istart,Iend;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-mu",&mu,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nLauchli hyperbolic singular value decomposition, (%" PetscInt_FMT " x %" PetscInt_FMT ") mu=%g\n\n",n+1,n,(double)mu));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   Build Lauchli matrix and signature
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n+1,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i == 0) {
      for (j=0;j<n;j++) PetscCall(MatSetValue(A,0,j,1.0,INSERT_VALUES));
    } else PetscCall(MatSetValue(A,i,i-1,mu,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* signature omega = [ -1 1 1 ... 1 ] */
  PetscCall(MatCreateVecs(A,NULL,&omega));
  PetscCall(VecSet(omega,1.0));
  PetscCall(VecGetOwnershipRange(omega,&Istart,NULL));
  if (Istart==0) PetscCall(VecSetValue(omega,0,-1.0,INSERT_VALUES));
  PetscCall(VecAssemblyBegin(omega));
  PetscCall(VecAssemblyEnd(omega));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Create the singular value solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));
  PetscCall(SVDSetSignature(svd,omega));
  PetscCall(SVDSetOperators(svd,A,NULL));
  PetscCall(SVDSetType(svd,SVDTRLANCZOS));
  PetscCall(SVDSetFromOptions(svd));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the problem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDSolve(svd));
  PetscCall(SVDErrorView(svd,SVD_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));

  /* Test getting the signature */
  PetscCall(MatCreateVecs(A,NULL,&v));
  PetscCall(SVDGetSignature(svd,v));
  PetscCall(VecView(v,NULL));

  /* Free work space */
  PetscCall(SVDDestroy(&svd));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&omega));
  PetscCall(VecDestroy(&v));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      requires: double
      suffix: 1

TEST*/
