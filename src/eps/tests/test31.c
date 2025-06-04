/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test STFILTER interface functions.\n\n"
  "Based on ex2.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;
  EPS            eps;
  ST             st;
  PetscInt       N,n=10,m,Istart,Iend,II,i,j,degree;
  PetscBool      flag,modify=PETSC_FALSE,terse;
  PetscReal      inta,intb,rleft,rright;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  if (!flag) m=n;
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n2-D Laplacian Eigenproblem, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,m));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-modify",&modify,&flag));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Create the 2-D Laplacian
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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetType(eps,EPSKRYLOVSCHUR));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_ALL));
  PetscCall(EPSSetInterval(eps,0.5,1.3));
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STSetType(st,STFILTER));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the problem and display the solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /* print filter information */
  PetscCall(PetscObjectTypeCompare((PetscObject)st,STFILTER,&flag));
  if (flag) {
    PetscCall(STFilterGetDegree(st,&degree));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Filter degree: %" PetscInt_FMT "\n",degree));
    PetscCall(STFilterGetInterval(st,&inta,&intb));
    PetscCall(STFilterGetRange(st,&rleft,&rright));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Requested interval: [%g,%g],  range: [%g,%g]\n\n",(double)inta,(double)intb,(double)rleft,(double)rright));
  }

  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              Solve the problem again after changing the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (modify) {
    PetscCall(MatSetValue(A,0,0,0.3,INSERT_VALUES));
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
    PetscCall(EPSSetOperators(eps,A,NULL));
    PetscCall(EPSSolve(eps));
    PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
    if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
    else {
      PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
      PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
      PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
      PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
    }
  }

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -terse
      filter: sed -e "s/0.161982,7.83797/0.162007,7.83897/"
      requires: !single

   test:
      suffix: 2
      args: -modify -st_filter_range -0.5,8 -terse
      requires: !single

   test:
      suffix: 3
      args: -modify -terse
      filter: sed -e "s/0.161982,7.83797/0.162007,7.83897/"
      requires: !single

   test:
      suffix: 4
      args: -terse -st_filter_type chebyshev -st_filter_damping {{none jackson}}
      filter: sed -e "s/0.161982,7.83797/0.162007,7.83897/"
      output_file: output/test31_1.out
      requires: !single

   test:
      suffix: 5
      args: -modify -terse -st_filter_range -0.5,8 -st_filter_type chebyshev -st_filter_damping {{lanczos fejer}}
      output_file: output/test31_2.out
      requires: !single

TEST*/
