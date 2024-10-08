/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test changing MFN type.\n\n";

#include <slepcmfn.h>

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  MFN            mfn;
  FN             f;
  PetscReal      norm;
  PetscScalar    t=0.3;
  PetscInt       N,n=25,m,Istart,Iend,II,i,j;
  PetscBool      flag;
  Vec            v,y;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  if (!flag) m=n;
  N = n*m;
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-t",&t,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nMatrix exponential y=exp(t*A)*e, of the 2-D Laplacian, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,m));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                         Build the 2-D Laplacian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,-1.0,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A,II,II+n,-1.0,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,-1.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,4.0,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* set v = ones(n,1) */
  PetscCall(MatCreateVecs(A,&v,&y));
  PetscCall(VecSet(v,1.0));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f));
  PetscCall(FNSetType(f,FNEXP));

  PetscCall(MFNCreate(PETSC_COMM_WORLD,&mfn));
  PetscCall(MFNSetOperator(mfn,A));
  PetscCall(MFNSetType(mfn,MFNEXPOKIT));
  PetscCall(MFNSetDimensions(mfn,24));
  PetscCall(MFNSetTolerances(mfn,1e-5,1000));
  PetscCall(MFNSetFN(mfn,f));
  PetscCall(MFNSetErrorIfNotConverged(mfn,PETSC_TRUE));
  PetscCall(MFNSetFromOptions(mfn));
  PetscCall(MFNView(mfn,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Change MFN type and solve the problem, y=exp(t*A)*v
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MFNSetType(mfn,MFNKRYLOV));
  PetscCall(FNSetScale(f,t,1.0));
  PetscCall(MFNSolve(mfn,v,y));
  PetscCall(VecNorm(y,NORM_2,&norm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Computed vector at time t=%.4g has norm %g\n\n",(double)PetscRealPart(t),(double)norm));

  /*
     Free work space
  */
  PetscCall(MFNDestroy(&mfn));
  PetscCall(FNDestroy(&f));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&y));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:

TEST*/
