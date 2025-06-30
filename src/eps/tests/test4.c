/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test the solution of a HEP without calling EPSSetFromOptions (based on ex1.c).\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n"
  "  -type <eps_type> = eps type to test.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  PetscReal      tol=1000*PETSC_MACHINE_EPSILON;
  PetscInt       n=30,i,Istart,Iend,nev;
  PetscBool      flg,gd2;
  char           epstype[30] = "krylovschur";

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetString(NULL,NULL,"-type",epstype,sizeof(epstype),NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%" PetscInt_FMT "\n\n",n));

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
  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetDimensions(eps,4,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));

  /*
     Set solver parameters at runtime
  */
  PetscCall(PetscStrcmp(epstype,"gd2",&flg));
  if (flg) {
    PetscCall(EPSSetType(eps,EPSGD));
    PetscCall(EPSGDSetDoubleExpansion(eps,PETSC_TRUE));
    PetscCall(EPSGDGetDoubleExpansion(eps,&gd2));  /* not used */
  } else PetscCall(EPSSetType(eps,epstype));
  PetscCall(PetscStrcmp(epstype,"jd",&flg));
  if (flg) {
    PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));
    PetscCall(EPSSetTarget(eps,4.0));
  }
  PetscCall(PetscStrcmp(epstype,"lanczos",&flg));
  if (flg) PetscCall(EPSLanczosSetReorthog(eps,EPS_LANCZOS_REORTHOG_LOCAL));
  PetscCall(PetscObjectTypeCompareAny((PetscObject)eps,&flg,EPSRQCG,EPSLOBPCG,""));
  if (flg) PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      filter: sed -e "s/3.95905/3.95906/"
      output_file: output/test4_1.out
      test:
         suffix: 1
         args: -type {{krylovschur subspace arnoldi lanczos gd jd gd2 lapack}}
      test:
         suffix: 1_arpack
         args: -type arpack
         requires: arpack
      test:
         suffix: 1_primme
         args: -type primme
         requires: primme

   test:
      suffix: 2
      args: -type {{rqcg lobpcg}}

TEST*/
